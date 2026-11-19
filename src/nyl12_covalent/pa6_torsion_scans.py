"""Prepare bidirectional relaxed DFT scans for PA6 bonded benchmarks."""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from nyl12_covalent.gaussian_checks import final_orientation

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class ScanSpec:
    """Define one relaxed proper or improper coordinate scan."""

    fragment: str
    label: str
    indices_zero_based: tuple[int, int, int, int]
    coordinate_kind: str
    steps: int
    increment_degrees: float
    parameter_target: str


SCANS = (
    ScanSpec(
        "internal_neutral",
        "amide_omega",
        (0, 1, 3, 4),
        "proper",
        24,
        15.0,
        "c3-c/ns amide proper terms",
    ),
    ScanSpec(
        "internal_neutral",
        "amine_side",
        (1, 3, 4, 5),
        "proper",
        24,
        15.0,
        "c3-ns-c3-c3 and ns-H coupled terms",
    ),
    ScanSpec(
        "internal_neutral",
        "central_methylene",
        (4, 5, 6, 7),
        "proper",
        24,
        15.0,
        "c3-c3-c3-c3 transferability control",
    ),
    ScanSpec(
        "internal_neutral",
        "carbonyl_side",
        (6, 7, 8, 9),
        "proper",
        24,
        15.0,
        "c3-c3-c3-c carbonyl-side torsion",
    ),
    ScanSpec(
        "internal_neutral",
        "amide_carbonyl_planarity",
        (0, 3, 1, 2),
        "improper",
        12,
        5.0,
        "c3-ns-c-o improper; parmchk2 penalty 6",
    ),
    ScanSpec(
        "internal_neutral",
        "amide_n_planarity",
        (1, 4, 3, 16),
        "improper",
        12,
        5.0,
        "c-c3-ns-hn improper; parmchk2 default",
    ),
    ScanSpec(
        "n_terminus_cation",
        "ammonium_rotation",
        (2, 1, 0, 10),
        "proper",
        24,
        15.0,
        "c3-c3-nz-hn terminal ammonium term",
    ),
    ScanSpec(
        "c_terminus_anion",
        "carboxylate_planarity",
        (8, 10, 9, 11),
        "improper",
        12,
        5.0,
        "c3-o-c-o improper; parmchk2 penalty 3",
    ),
)


def prepare_relaxed_scans(geometry_dir: Path, output_dir: Path) -> dict[str, Any]:
    """Create forward/reverse Gaussian scans from accepted optimizations.

    Args:
        geometry_dir: Fragment directory containing manifest, SDF seeds, and
            normally terminated Gaussian optimization logs.
        output_dir: New output directory for all scan inputs.

    Returns:
        Serializable scan manifest.

    Raises:
        FileExistsError: If ``output_dir`` exists.
        RuntimeError: If a geometry optimization is incomplete.
        ValueError: If atom identity/order or a scan declaration is invalid.
    """
    if output_dir.exists():
        raise FileExistsError(f"Scan output already exists: {output_dir}")
    source_manifest = json.loads((geometry_dir / "manifest.json").read_text())
    records_by_label = {
        record["label"]: record for record in source_manifest["fragments"]
    }
    cache: dict[str, tuple[list[int], np.ndarray, list[str], int]] = {}
    records = []
    output_dir.mkdir(parents=True)
    for scan in SCANS:
        if scan.fragment not in cache:
            cache[scan.fragment] = _accepted_geometry(
                geometry_dir, records_by_label[scan.fragment]
            )
        elements, xyz, names, charge = cache[scan.fragment]
        _validate_scan(scan, len(elements))
        directions = {}
        for direction, sign in (("forward", 1.0), ("reverse", -1.0)):
            target = output_dir / scan.fragment / scan.label / direction
            target.mkdir(parents=True)
            input_path = target / "scan.gjf"
            input_path.write_text(
                _scan_input(elements, xyz, charge, scan, direction, sign)
            )
            directions[direction] = {
                "input": str(input_path),
                "sha256": _sha256(input_path),
                "increment_degrees": sign * scan.increment_degrees,
            }
        records.append(
            {
                "fragment": scan.fragment,
                "label": scan.label,
                "coordinate_kind": scan.coordinate_kind,
                "atom_indices_zero_based": list(scan.indices_zero_based),
                "atom_indices_gaussian_one_based": [
                    index + 1 for index in scan.indices_zero_based
                ],
                "atom_names": [names[index] for index in scan.indices_zero_based],
                "steps_per_direction": scan.steps,
                "parameter_target": scan.parameter_target,
                "directions": directions,
            }
        )
    manifest = {
        "schema_version": 1,
        "scope": "Bidirectional PA6 relaxed DFT torsion/improper scans",
        "method": "wB97XD/def2SVP/SMD(water)",
        "scan_count": len(records),
        "gaussian_job_count": 2 * len(records),
        "source_manifest_sha256": _sha256(geometry_dir / "manifest.json"),
        "scans": records,
        "acceptance_requirements": [
            "all windows normally terminate and optimize",
            "forward/reverse profiles agree within declared hysteresis tolerance",
            "non-target coordinates have no discontinuous jumps",
            "def2-TZVP single points are evaluated on accepted scan geometries",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared %d bidirectional PA6 scans", len(records))
    return manifest


def _accepted_geometry(
    geometry_dir: Path, record: dict[str, Any]
) -> tuple[list[int], np.ndarray, list[str], int]:
    label = record["label"]
    fragment_dir = geometry_dir / label
    log_path = fragment_dir / "optimization.log"
    if not log_path.is_file():
        raise RuntimeError(f"Fragment optimization incomplete: {label}")
    text = log_path.read_text()
    if (
        text.count("Normal termination of Gaussian 16") != 1
        or "Optimization completed." not in text
        or "Error termination" in text
    ):
        raise RuntimeError(f"Fragment optimization incomplete: {label}")
    elements, xyz = final_orientation(text)
    molecules = list(
        Chem.SDMolSupplier(str(fragment_dir / "seed_mmff94s.sdf"), removeHs=False)
    )
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError(f"Invalid fragment SDF: {label}")
    molecule = molecules[0]
    expected = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    names = list(record["atom_names_in_order"])
    if len(names) != len(expected):
        raise ValueError(f"Manifest atom names differ from seed: {label}")
    if elements != expected or xyz.shape != (len(expected), 3):
        raise ValueError(f"Optimized atom order differs from seed: {label}")
    _validate_expected_bonds(molecule, xyz, label)
    return elements, xyz, names, int(record["charge"])


def _validate_expected_bonds(molecule: Chem.Mol, xyz: np.ndarray, label: str) -> None:
    """Reject optimized fragments with broken or collapsed declared bonds."""
    periodic = Chem.GetPeriodicTable()
    for bond in molecule.GetBonds():
        first = bond.GetBeginAtomIdx()
        second = bond.GetEndAtomIdx()
        expected = periodic.GetRcovalent(
            molecule.GetAtomWithIdx(first).GetAtomicNum()
        ) + periodic.GetRcovalent(molecule.GetAtomWithIdx(second).GetAtomicNum())
        distance = float(np.linalg.norm(xyz[first] - xyz[second]))
        if not expected - 0.25 <= distance <= expected + 0.25:
            raise ValueError(
                f"Optimized bond geometry invalid for {label}: "
                f"{first}-{second} is {distance:.4f} A"
            )


def _validate_scan(scan: ScanSpec, atom_count: int) -> None:
    if scan.coordinate_kind not in {"proper", "improper"}:
        raise ValueError(f"Unknown coordinate kind: {scan.coordinate_kind}")
    if len(set(scan.indices_zero_based)) != 4 or any(
        index < 0 or index >= atom_count for index in scan.indices_zero_based
    ):
        raise ValueError(f"Invalid scan atom indices: {scan.label}")
    if scan.steps < 1 or not np.isfinite(scan.increment_degrees):
        raise ValueError(f"Invalid scan schedule: {scan.label}")


def _scan_input(
    elements: list[int],
    xyz: np.ndarray,
    charge: int,
    scan: ScanSpec,
    direction: str,
    sign: float,
) -> str:
    periodic = Chem.GetPeriodicTable()
    checkpoint = f"{scan.fragment}_{scan.label}_{direction}.chk"
    lines = [
        "%nprocshared=4",
        "%mem=8GB",
        f"%chk={checkpoint}",
        "#p wB97XD/def2SVP SCRF=(SMD,Solvent=Water) "
        "Opt=(ModRedundant,Tight,MaxCycles=200) SCF=(Tight,XQC) NoSymm",
        "",
        f"PA6 benchmark {scan.fragment} {scan.label} {direction}",
        "",
        f"{charge} 1",
    ]
    lines.extend(
        f"{periodic.GetElementSymbol(number):2s} {x: .10f} {y: .10f} {z: .10f}"
        for number, (x, y, z) in zip(elements, xyz, strict=True)
    )
    indices = " ".join(str(index + 1) for index in scan.indices_zero_based)
    lines.extend(
        [
            "",
            f"D {indices} S {scan.steps} {sign * scan.increment_degrees:.6f}",
            "",
        ]
    )
    return "\n".join(lines)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
