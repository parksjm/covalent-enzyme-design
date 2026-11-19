"""Prepare constrained PA6 L4 conformers for RESP electrostatic sampling."""

from __future__ import annotations

from dataclasses import dataclass
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
CORES = 8
MEMORY_GB = 16
ROUTE = (
    "#p wB97XD/def2TZVP SCRF=(SMD,Solvent=Water) Int=UltraFine "
    "Opt=(Tight,CalcFC,ModRedundant,MaxCycles=300) SCF=(Tight,XQC) NoSymm"
)
LABEL_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


@dataclass(frozen=True)
class Torsion:
    """A named four-atom heavy-atom torsion frozen during optimization.

    Attributes:
        name: Human-readable, unique torsion identifier.
        atoms: Four canonical PA6 atom names in bonded order.
    """

    name: str
    atoms: tuple[str, str, str, str]


def prepare_pa6_resp_conformers(
    source_mol2: Path,
    conformers: dict[str, Path],
    torsions_json: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Write canonical XYZ and constrained Gaussian inputs for PA6 conformers.

    Args:
        source_mol2: Canonical 79-atom PA6 L4 zwitterion MOL2.
        conformers: Explicit label-to-PDB mapping; each PDB must contain one
            complete residue named LIG with coordinates keyed by atom name.
        torsions_json: JSON list of objects containing ``name`` and four
            ``atoms``. Every atom must be heavy and consecutively bonded.
        output_dir: New directory for prepared, nonexecuted inputs.

    Returns:
        Serializable manifest written to ``manifest.json``.

    Raises:
        FileExistsError: If ``output_dir`` already exists.
        ValueError: If chemistry, mappings, labels, or torsions are invalid.
    """
    if output_dir.exists():
        raise FileExistsError(f"Output already exists: {output_dir}")
    if not conformers:
        raise ValueError("At least one explicit conformer is required")
    _validate_labels(conformers)
    substrate = load_pa6(source_mol2)
    validate_pa6_zwitterion(substrate.mol)
    atom_names = [_atom_name(atom) for atom in substrate.mol.GetAtoms()]
    if len(set(atom_names)) != len(atom_names):
        raise ValueError("Canonical MOL2 atom names are not unique")
    torsions = _load_torsions(torsions_json, substrate.mol, atom_names)

    mapped_conformers = [
        (label, pdb_path, _read_lig_coordinates(pdb_path, set(atom_names)))
        for label, pdb_path in sorted(conformers.items())
    ]
    output_dir.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, Any]] = []
    for label, pdb_path, coordinates in mapped_conformers:
        rows = [
            (
                atom.GetSymbol(),
                coordinates[atom_names[atom.GetIdx()]][0],
                coordinates[atom_names[atom.GetIdx()]][1],
                coordinates[atom_names[atom.GetIdx()]][2],
            )
            for atom in substrate.mol.GetAtoms()
        ]
        conformer_dir = output_dir / label
        conformer_dir.mkdir()
        xyz_path = conformer_dir / "coordinates.xyz"
        gaussian_path = conformer_dir / "opt_smd.gjf"
        xyz_path.write_text(_xyz_text(rows, label))
        gaussian_path.write_text(_gaussian_text(rows, torsions, atom_names, label))
        records.append(
            {
                "label": label,
                "source_pdb": str(pdb_path.resolve()),
                "source_pdb_sha256": _sha256(pdb_path),
                "coordinates_xyz": str(xyz_path.relative_to(output_dir)),
                "coordinates_xyz_sha256": _sha256(xyz_path),
                "gaussian_input": str(gaussian_path.relative_to(output_dir)),
                "gaussian_input_sha256": _sha256(gaussian_path),
            }
        )
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "chemistry": "PA6 L4 zwitterion, explicit NH3+/COO-, net charge 0",
        "atom_count": len(atom_names),
        "atom_names_in_canonical_order": atom_names,
        "source_mol2": str(source_mol2.resolve()),
        "source_mol2_sha256": _sha256(source_mol2),
        "torsions_json": str(torsions_json.resolve()),
        "torsions_json_sha256": _sha256(torsions_json),
        "frozen_torsions": [
            {"name": torsion.name, "atoms": list(torsion.atoms)} for torsion in torsions
        ],
        "gaussian": {
            "route": ROUTE,
            "cores": CORES,
            "memory_gb": MEMORY_GB,
            "charge": 0,
            "multiplicity": 1,
            "jobs_launched": False,
        },
        "conformers": records,
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared %d PA6 RESP conformers at %s", len(records), output_dir)
    return manifest


def _validate_labels(conformers: dict[str, Path]) -> None:
    for label in conformers:
        if not LABEL_PATTERN.fullmatch(label):
            raise ValueError(f"Invalid conformer label: {label!r}")


def _atom_name(atom: Chem.Atom) -> str:
    if not atom.HasProp("name"):
        raise ValueError(f"Canonical atom {atom.GetIdx() + 1} lacks a name")
    return atom.GetProp("name")


def _load_torsions(
    path: Path, molecule: Chem.Mol, atom_names: list[str]
) -> list[Torsion]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, list) or not payload:
        raise ValueError("Torsion JSON must be a nonempty list")
    by_name = {name: index for index, name in enumerate(atom_names)}
    torsions: list[Torsion] = []
    seen_names: set[str] = set()
    for record in payload:
        if not isinstance(record, dict) or set(record) != {"name", "atoms"}:
            raise ValueError("Each torsion must contain only 'name' and 'atoms'")
        name = record["name"]
        atoms = record["atoms"]
        if not isinstance(name, str) or not name or name in seen_names:
            raise ValueError(f"Torsion name must be nonempty and unique: {name!r}")
        if not isinstance(atoms, list) or len(atoms) != 4:
            raise ValueError(f"Torsion {name!r} must contain four atom names")
        if any(not isinstance(atom, str) or atom not in by_name for atom in atoms):
            raise ValueError(f"Torsion {name!r} references an unknown atom")
        indices = [by_name[atom] for atom in atoms]
        if len(set(indices)) != 4:
            raise ValueError(f"Torsion {name!r} repeats an atom")
        if any(molecule.GetAtomWithIdx(index).GetAtomicNum() == 1 for index in indices):
            raise ValueError(f"Torsion {name!r} must contain only heavy atoms")
        if any(
            molecule.GetBondBetweenAtoms(first, second) is None
            for first, second in zip(indices, indices[1:])
        ):
            raise ValueError(f"Torsion {name!r} atoms are not consecutively bonded")
        seen_names.add(name)
        torsions.append(Torsion(name, tuple(atoms)))
    return torsions


def _read_lig_coordinates(
    path: Path, expected_names: set[str]
) -> dict[str, tuple[float, float, float]]:
    coordinates: dict[str, tuple[float, float, float]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[17:20].strip() != "LIG":
            continue
        altloc = line[16:17]
        if altloc.strip():
            raise ValueError(f"Alternate location in {path}: {altloc!r}")
        name = line[12:16].strip()
        if name in coordinates:
            raise ValueError(f"Duplicate LIG atom name {name!r} in {path}")
        try:
            coordinates[name] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError as error:
            raise ValueError(f"Invalid coordinates for LIG atom {name!r}") from error
    found_names = set(coordinates)
    if found_names != expected_names:
        missing = sorted(expected_names - found_names)
        extra = sorted(found_names - expected_names)
        raise ValueError(
            f"LIG atom-name mismatch in {path}; missing={missing}, extra={extra}"
        )
    return coordinates


def _xyz_text(rows: list[tuple[str, float, float, float]], label: str) -> str:
    body = [f"{element} {x:.8f} {y:.8f} {z:.8f}" for element, x, y, z in rows]
    return "\n".join(
        [str(len(rows)), f"PA6 L4 zwitterion conformer {label}", *body, ""]
    )


def _gaussian_text(
    rows: list[tuple[str, float, float, float]],
    torsions: list[Torsion],
    atom_names: list[str],
    label: str,
) -> str:
    indices = {name: index + 1 for index, name in enumerate(atom_names)}
    lines = [
        f"%nprocshared={CORES}",
        f"%mem={MEMORY_GB}GB",
        f"%chk={label}_opt_smd.chk",
        ROUTE,
        "",
        f"PA6 L4 zwitterion {label}; constrained RESP geometry sampling",
        "",
        "0 1",
    ]
    lines.extend(f"{element} {x:.10f} {y:.10f} {z:.10f}" for element, x, y, z in rows)
    lines.append("")
    lines.extend(
        "D " + " ".join(str(indices[atom]) for atom in torsion.atoms) + " F"
        for torsion in torsions
    )
    return "\n".join([*lines, "", ""])


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
