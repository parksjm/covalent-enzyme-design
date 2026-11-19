"""Prepare matched embedded-QM calculations along an enzyme reaction path."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import logging
from pathlib import Path
import re
import shutil

LOGGER = logging.getLogger(__name__)
ANGSTROM_TO_BOHR = 1.8897261254578281
WINDOWS = ("pt_01", "pt_05", "pt_10", "pt_14", "pt_19")
QMMASK = (
    "@1621-1655,2727-2734,13905-13931,14555-14560,14574-14579,"
    "19914-19917,19925-19931,19937-19941"
)
ACTIVE_MASK = "(:113,186,929,930,968,970,1313-1316)<:12.0"
CHARGE_MULTIPLICITY = re.compile(r"^\s*(-?\d+)\s+(\d+)\s*$")


@dataclass(frozen=True)
class QmAtom:
    """One QM or link atom from an Amber-generated Gaussian input."""

    symbol: str
    x_angstrom: float
    y_angstrom: float
    z_angstrom: float


@dataclass(frozen=True)
class PointCharge:
    """One adjusted MM point charge from an Amber-generated Gaussian input."""

    x_angstrom: float
    y_angstrom: float
    z_angstrom: float
    charge_e: float


@dataclass(frozen=True)
class GaussianEmbedding:
    """QM/link geometry and external charges parsed from Gaussian input."""

    charge: int
    multiplicity: int
    atoms: tuple[QmAtom, ...]
    point_charges: tuple[PointCharge, ...]


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file.

    Args:
        path: Existing file.

    Returns:
        Lowercase hexadecimal digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_gaussian_embedding(path: Path) -> GaussianEmbedding:
    """Parse Amber's Gaussian QM geometry and point-charge sections.

    Args:
        path: Gaussian input retained by Sander's external interface.

    Returns:
        Parsed QM/link atoms and external MM charges.

    Raises:
        ValueError: If a required section is absent or malformed.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    charge_index = next(
        (index for index, line in enumerate(lines) if CHARGE_MULTIPLICITY.match(line)),
        None,
    )
    if charge_index is None:
        raise ValueError(f"Missing Gaussian charge/multiplicity line: {path}")
    match = CHARGE_MULTIPLICITY.match(lines[charge_index])
    assert match is not None
    atoms: list[QmAtom] = []
    index = charge_index + 1
    while index < len(lines) and lines[index].strip():
        fields = lines[index].split()
        if len(fields) != 4:
            raise ValueError(f"Malformed Gaussian QM coordinate: {lines[index]}")
        atoms.append(QmAtom(fields[0], *(float(value) for value in fields[1:])))
        index += 1
    index += 1
    charges: list[PointCharge] = []
    while index < len(lines) and lines[index].strip():
        fields = lines[index].split()
        if len(fields) != 4:
            raise ValueError(f"Malformed Gaussian point charge: {lines[index]}")
        x, y, z, charge = (float(value) for value in fields)
        charges.append(PointCharge(x, y, z, charge))
        index += 1
    if len(atoms) != 106:
        raise ValueError(f"Expected 106 QM/link atoms, found {len(atoms)}")
    if not charges:
        raise ValueError("Gaussian input contains no external point charges")
    return GaussianEmbedding(
        charge=int(match.group(1)),
        multiplicity=int(match.group(2)),
        atoms=tuple(atoms),
        point_charges=tuple(charges),
    )


def write_xtb_embedding(gaussian_input: Path, output: Path) -> dict[str, object]:
    """Convert one Amber/Gaussian embedding into standalone xTB inputs.

    Args:
        gaussian_input: Retained Gaussian input from the external QM/MM interface.
        output: Existing directory receiving xTB inputs.

    Returns:
        Input metadata and SHA-256 digests.
    """
    embedding = parse_gaussian_embedding(gaussian_input)
    geometry = output / "qm_link.xyz"
    point_charges = output / "point_charges.pc"
    control = output / "xtb.inp"
    atom_lines = [
        f"{atom.symbol:<2s} {atom.x_angstrom: .15f} {atom.y_angstrom: .15f} "
        f"{atom.z_angstrom: .15f}"
        for atom in embedding.atoms
    ]
    geometry.write_text(
        f"{len(atom_lines)}\nAmber QM atoms plus link hydrogens\n"
        + "\n".join(atom_lines)
        + "\n",
        encoding="utf-8",
    )
    charge_lines = [
        f"{charge.charge_e: .16f} "
        f"{charge.x_angstrom * ANGSTROM_TO_BOHR: .16f} "
        f"{charge.y_angstrom * ANGSTROM_TO_BOHR: .16f} "
        f"{charge.z_angstrom * ANGSTROM_TO_BOHR: .16f} 999.0"
        for charge in embedding.point_charges
    ]
    point_charges.write_text(
        f"{len(charge_lines)}\n" + "\n".join(charge_lines) + "\n",
        encoding="utf-8",
    )
    control.write_text(
        "$embedding\n  interface=legacy\n  input=point_charges.pc\n$end\n",
        encoding="utf-8",
    )
    metadata = {
        "charge": embedding.charge,
        "multiplicity": embedding.multiplicity,
        "qm_link_atom_count": len(embedding.atoms),
        "point_charge_count": len(embedding.point_charges),
        "point_charge_hardness": 999.0,
        "source_gaussian_input_sha256": sha256(gaussian_input),
        "sha256": {
            path.name: sha256(path) for path in (geometry, point_charges, control)
        },
    }
    (output / "xtb_input_manifest.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return metadata


def _mdin(qm_theory: str) -> str:
    gaussian = ""
    if qm_theory == "EXTERN":
        gaussian = """
&gau
 method='wB97XD', basis='def2SVP', mem='32GB',
 executable='/usr/local/g16/g16', guess='read',
 charge_analysis='MULLIKEN', scf_conv=8, ntpr=1,
 num_threads=16, debug=1, dipole=0, use_template=1, keep_files=1,
/
"""
    return f"""Matched embedded-QM single energy/gradient evaluation
&cntrl
 imin=1, maxcyc=1, ntmin=5, ncyc=1, dx0=0.02, drms=0.05,
 ntb=1, igb=0, ntc=1, ntf=1, cut=10.0, ifqnt=1, jfastw=4,
 ntr=0, ntpr=1, ntwx=1, ntwr=1, ioutfm=1, ntxo=2,
/
&qmmm
 qmmask='{QMMASK}', qmcharge=-2, spin=1, qm_theory='{qm_theory}',
 qmmm_int=1, adjust_q=2, qm_ewald=0, qm_pme=0, qmcut=20.0,
 scfconv=1.0d-10, lnk_dis=1.09, lnk_atomic_no=1, lnk_method=1,
/
{gaussian}&dlfind
 optalg='LBFGS', crdrep='CART', active='{ACTIVE_MASK}',
 tol=1.d-4, tole=1.d-5, maxcycle=1, maxstep=0.25, verbosity=1,
/
"""


def prepare_benchmark(coordinate_drive: Path, topology: Path, output: Path) -> None:
    """Create immutable five-window inputs for the matched method comparison.

    Args:
        coordinate_drive: Direct proton-transfer coordinate-drive directory.
        topology: Matching solvated Amber topology.
        output: New benchmark input directory.
    """
    coordinate_drive = coordinate_drive.resolve(strict=True)
    topology = topology.resolve(strict=True)
    output.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, object]] = []
    for window in WINDOWS:
        source = coordinate_drive / window / "optimized.ncrst"
        acceptance = coordinate_drive / window / "acceptance.json"
        if not source.is_file() or not acceptance.is_file():
            raise ValueError(f"Incomplete source window: {window}")
        window_dir = output / "windows" / window
        for method in ("dftb3_mm", "wb97xd_mm", "gfn2_xtb"):
            (window_dir / method).mkdir(parents=True)
        copied = window_dir / "input.ncrst"
        shutil.copy2(source, copied)
        shutil.copy2(acceptance, window_dir / "acceptance.json")
        (window_dir / "dftb3_mm" / "minimize.mdin").write_text(
            _mdin("DFTB3"), encoding="utf-8"
        )
        (window_dir / "wb97xd_mm" / "minimize.mdin").write_text(
            _mdin("EXTERN"), encoding="utf-8"
        )
        (window_dir / "wb97xd_mm" / "gau_job.tpl").write_text(
            "#P wB97XD/def2SVP Int=UltraFine SCF=(XQC,Tight,MaxCycle=512)\n\n",
            encoding="utf-8",
        )
        records.append(
            {
                "window": window,
                "source_restart": str(source),
                "source_restart_sha256": sha256(source),
                "copied_restart_sha256": sha256(copied),
                "acceptance": json.loads(acceptance.read_text()),
            }
        )
    manifest = {
        "benchmark": "direct_proton_transfer_embedded_method_comparison",
        "windows": list(WINDOWS),
        "topology": str(topology),
        "topology_sha256": sha256(topology),
        "qm_mask": QMMASK,
        "qm_atom_count_excluding_links": 98,
        "expected_link_atom_count": 8,
        "qm_charge": -2,
        "multiplicity": 1,
        "embedding": "adjusted and imaged MM point charges within 20 angstrom",
        "qm_ewald": False,
        "qm_pme": False,
        "methods": ["DFTB3/3OB", "GFN2-xTB 6.7.1", "omegaB97X-D/def2-SVP"],
        "records": records,
    }
    (output / "benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info("Prepared %d matched enzyme-path windows in %s", len(records), output)
