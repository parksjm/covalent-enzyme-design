"""Build new constrained Gaussian inputs from a prepared PA6 XYZ manifest."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import logging
from pathlib import Path
import re
from typing import Any

from rdkit import Chem

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.ti1_reference_graph import validate_pa6_zwitterion

LOGGER = logging.getLogger(__name__)
TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9()+_.-]*$")


@dataclass(frozen=True)
class GaussianSettings:
    """Explicit geometry optimization settings for a restart attempt."""

    method: str
    basis: str
    solvent: str
    max_cycles: int
    cores: int
    memory_gb: int
    calc_fc: bool

    def __post_init__(self) -> None:
        """Reject incomplete or unsafe route settings."""
        for field in ("method", "basis"):
            if not TOKEN.fullmatch(getattr(self, field)):
                raise ValueError(f"Invalid Gaussian {field}: {getattr(self, field)!r}")
        if not self.solvent.startswith("SMD:") or not TOKEN.fullmatch(
            self.solvent.removeprefix("SMD:")
        ):
            raise ValueError("Solvent must be specified as SMD:<Gaussian solvent>")
        if self.max_cycles < 1 or self.cores < 1 or self.memory_gb < 1:
            raise ValueError("Max cycles, cores, and memory must all be positive")

    @property
    def route(self) -> str:
        """Return the validated Gaussian route section."""
        solvent = self.solvent.removeprefix("SMD:")
        calc_fc = ",CalcFC" if self.calc_fc else ""
        return (
            f"#p {self.method}/{self.basis} SCRF=(SMD,Solvent={solvent}) Int=UltraFine "
            f"Opt=(Tight{calc_fc},ModRedundant,MaxCycles={self.max_cycles}) "
            "SCF=(Tight,XQC) NoSymm"
        )


def prepare_pa6_resp_restart(
    source_manifest: Path, output_dir: Path, settings: GaussianSettings
) -> dict[str, Any]:
    """Create non-overwriting inputs from every prepared conformer in a manifest.

    Args:
        source_manifest: Existing manifest with canonical atom order, XYZ files,
            and frozen torsion declarations.
        output_dir: New directory for this independent restart attempt.
        settings: Fully explicit method, solvent, resource, and optimization data.

    Returns:
        New provenance manifest with input/source hashes and no job launch.

    Raises:
        FileExistsError: If the destination already exists.
        ValueError: If source provenance, XYZ atom order, or torsions are invalid.
    """
    if output_dir.exists():
        raise FileExistsError(f"Output already exists: {output_dir}")
    payload = json.loads(source_manifest.read_text())
    source_mol2 = Path(payload["source_mol2"])
    molecule = load_pa6(source_mol2).mol
    validate_pa6_zwitterion(molecule)
    atom_names = [atom.GetProp("name") for atom in molecule.GetAtoms()]
    if atom_names != payload.get("atom_names_in_canonical_order"):
        raise ValueError("Source manifest atom order differs from canonical MOL2")
    if (
        not isinstance(payload.get("frozen_torsions"), list)
        or len(payload["frozen_torsions"]) != 25
    ):
        raise ValueError("Source manifest must declare the same 25 frozen torsions")
    source_hash = _sha256(source_manifest)
    prepared: list[tuple[str, Path, list[str], list[list[float]]]] = []
    for conformer in payload.get("conformers", []):
        label = conformer["label"]
        xyz_path = source_manifest.parent / conformer["coordinates_xyz"]
        if _sha256(xyz_path) != conformer.get("coordinates_xyz_sha256"):
            raise ValueError(f"Prepared XYZ hash mismatch for {label}")
        elements, rows = _read_xyz(xyz_path, molecule)
        prepared.append((label, xyz_path, elements, rows))
    if not prepared:
        raise ValueError("Source manifest contains no conformers")

    output_dir.mkdir(parents=True, exist_ok=False)
    records = []
    for label, xyz_path, elements, rows in prepared:
        target_dir = output_dir / label
        target_dir.mkdir()
        copied_xyz = target_dir / "coordinates.xyz"
        target = target_dir / "opt_smd.gjf"
        copied_xyz.write_bytes(xyz_path.read_bytes())
        target.write_text(
            _gaussian_input(
                label, elements, rows, atom_names, payload["frozen_torsions"], settings
            )
        )
        records.append(
            {
                "label": label,
                "source_xyz": str(xyz_path.resolve()),
                "source_xyz_sha256": _sha256(xyz_path),
                "coordinates_xyz": str(copied_xyz.relative_to(output_dir)),
                "coordinates_xyz_sha256": _sha256(copied_xyz),
                "gaussian_input": str(target.relative_to(output_dir)),
                "gaussian_input_sha256": _sha256(target),
                "atom_count": len(atom_names),
            }
        )
    manifest = {
        "schema_version": 1,
        "purpose": "PA6 L4 constrained geometry restart inputs",
        "source_manifest": str(source_manifest.resolve()),
        "source_manifest_sha256": source_hash,
        "source_mol2": str(source_mol2.resolve()),
        "source_mol2_sha256": _sha256(source_mol2),
        "atom_names_in_canonical_order": atom_names,
        "atom_count": len(atom_names),
        "frozen_torsions": payload["frozen_torsions"],
        "gaussian": {
            **asdict(settings),
            "solvation_model": "SMD",
            "solvent": settings.solvent.removeprefix("SMD:"),
            "charge": 0,
            "multiplicity": 1,
            "route": settings.route,
            "jobs_launched": False,
        },
        "conformers": records,
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared %d PA6 restart inputs at %s", len(records), output_dir)
    return manifest


def _read_xyz(path: Path, molecule: Chem.Mol) -> tuple[list[str], list[list[float]]]:
    lines = path.read_text().splitlines()
    if len(lines) < 2 or int(lines[0]) != len(lines) - 2:
        raise ValueError(f"Malformed XYZ: {path}")
    elements, rows = [], []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"Malformed XYZ coordinate in {path}")
        elements.append(fields[0])
        rows.append([float(value) for value in fields[1:]])
    expected = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    if elements != expected:
        raise ValueError(f"XYZ element/order differs from canonical molecule: {path}")
    return elements, rows


def _gaussian_input(
    label: str,
    elements: list[str],
    rows: list[list[float]],
    names: list[str],
    torsions: list[dict[str, Any]],
    settings: GaussianSettings,
) -> str:
    lines = [
        f"%nprocshared={settings.cores}",
        f"%mem={settings.memory_gb}GB",
        f"%chk={label}_opt_smd_restart.chk",
        settings.route,
        "",
        f"PA6 L4 zwitterion {label}; independent constrained restart",
        "",
        "0 1",
    ]
    lines.extend(
        f"{element} {x:.10f} {y:.10f} {z:.10f}"
        for element, (x, y, z) in zip(elements, rows, strict=True)
    )
    lines.append("")
    indices = {name: index + 1 for index, name in enumerate(names)}
    lines.extend(
        "D " + " ".join(str(indices[name]) for name in item["atoms"]) + " F"
        for item in torsions
    )
    return "\n".join([*lines, "", ""])


def _sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()
