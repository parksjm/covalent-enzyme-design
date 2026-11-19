"""Prepare non-overwriting recovery inputs for the PA6 amide omega scan."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from nyl12_covalent.gaussian_checks import final_orientation

LOGGER = logging.getLogger(__name__)

EXPECTED_INDICES = (0, 1, 3, 4)
EXPECTED_STEPS = 24
EXPECTED_INCREMENT = 15.0
FORWARD_STEPS = 45
FORWARD_INCREMENT = 5.0
REVERSE_STEPS = 72
REVERSE_INCREMENT = -5.0
EXPECTED_OPTIMIZATIONS = 10
EXPECTED_ERROR_TERMINATIONS = 1
METHOD = "wB97XD/def2SVP/SMD(water)"
ROUTE = (
    "#p wB97XD/def2SVP SCRF=(SMD,Solvent=Water) "
    "Opt=(ModRedundant,Tight,CalcFC,MaxCycles=300) "
    "SCF=(Tight,XQC) NoSymm"
)


def prepare_amide_omega_recovery(
    source_scan_dir: Path, source_manifest_path: Path, output_dir: Path
) -> dict[str, Any]:
    """Prepare passed forward-continuation and fresh reverse scan inputs.

    Args:
        source_scan_dir: Original ``amide_omega`` scan directory.
        source_manifest_path: Parent relaxed-scan manifest.
        output_dir: New directory for recovery inputs and manifest.

    Returns:
        Serializable recovery manifest with source identities and schedules.

    Raises:
        FileExistsError: If ``output_dir`` already exists.
        FileNotFoundError: If any required source artifact is missing.
        RuntimeError: If the original scan log fails strict recovery checks.
        ValueError: If source metadata, route, indices, or geometry is invalid.
    """
    if output_dir.exists():
        raise FileExistsError(f"Recovery output already exists: {output_dir}")

    source_scan_dir = source_scan_dir.resolve()
    source_manifest_path = source_manifest_path.resolve()
    forward_gjf = source_scan_dir / "forward" / "scan.gjf"
    forward_log = source_scan_dir / "forward" / "scan.log"
    for source in (source_manifest_path, forward_gjf, forward_log):
        if not source.is_file():
            raise FileNotFoundError(f"Required recovery source is missing: {source}")

    parent_manifest = json.loads(source_manifest_path.read_text())
    parent_record = _validated_parent_record(
        parent_manifest, source_manifest_path, source_scan_dir, forward_gjf
    )
    seed_text = forward_gjf.read_text()
    seed_elements, seed_xyz, charge, multiplicity = _parse_gaussian_input(seed_text)
    _validate_source_route(seed_text)

    log_text = forward_log.read_text()
    last_optimization = _validate_forward_log(log_text)
    accepted_elements, accepted_xyz = final_orientation(
        log_text[: last_optimization.end()]
    )
    if accepted_elements != seed_elements or accepted_xyz.shape != seed_xyz.shape:
        raise ValueError(
            "Last accepted orientation differs from original seed atom order"
        )

    forward_text = _render_scan(
        accepted_elements,
        accepted_xyz,
        charge,
        multiplicity,
        "forward",
        FORWARD_STEPS,
        FORWARD_INCREMENT,
    )
    reverse_text = _render_scan(
        seed_elements,
        seed_xyz,
        charge,
        multiplicity,
        "reverse",
        REVERSE_STEPS,
        REVERSE_INCREMENT,
    )
    manifest = _recovery_manifest(
        source_manifest_path,
        forward_gjf,
        forward_log,
        parent_record,
        forward_text,
        reverse_text,
    )

    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "forward").mkdir()
    (output_dir / "reverse").mkdir()
    (output_dir / "forward" / "scan.gjf").write_text(forward_text)
    (output_dir / "reverse" / "scan.gjf").write_text(reverse_text)
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared passed amide omega recovery inputs in %s", output_dir)
    return manifest


def _validated_parent_record(
    manifest: dict[str, Any],
    manifest_path: Path,
    source_scan_dir: Path,
    forward_gjf: Path,
) -> dict[str, Any]:
    if manifest.get("method") != METHOD:
        raise ValueError("Parent manifest method is not the required wB97XD/SMD water")
    matches = [
        record
        for record in manifest.get("scans", [])
        if record.get("fragment") == "internal_neutral"
        and record.get("label") == "amide_omega"
    ]
    if len(matches) != 1:
        raise ValueError("Parent manifest must contain exactly one amide_omega record")
    record = matches[0]
    if record.get("coordinate_kind") != "proper":
        raise ValueError("Parent amide_omega coordinate must be a proper torsion")
    if tuple(record.get("atom_indices_zero_based", ())) != EXPECTED_INDICES:
        raise ValueError("Parent amide_omega atom indices must be 1, 2, 4, 5")
    if record.get("atom_indices_gaussian_one_based") != [1, 2, 4, 5]:
        raise ValueError("Parent Gaussian atom indices must be 1, 2, 4, 5")
    if record.get("steps_per_direction") != EXPECTED_STEPS:
        raise ValueError("Parent amide_omega schedule must have 24 steps")
    directions = record.get("directions", {})
    if directions.get("forward", {}).get("increment_degrees") != EXPECTED_INCREMENT:
        raise ValueError("Parent forward schedule must increment by +15 degrees")
    if directions.get("reverse", {}).get("increment_degrees") != -EXPECTED_INCREMENT:
        raise ValueError("Parent reverse schedule must increment by -15 degrees")
    declared_input = Path(directions.get("forward", {}).get("input", ""))
    if not declared_input.is_absolute():
        declared_input = manifest_path.parent / declared_input
    if declared_input.resolve() != forward_gjf:
        raise ValueError(
            "Parent manifest forward input does not identify the source GJF"
        )
    if directions.get("forward", {}).get("sha256") != _sha256(forward_gjf):
        raise ValueError("Parent manifest forward GJF hash does not match source")
    if source_scan_dir.name != "amide_omega":
        raise ValueError("Source scan directory must be named amide_omega")
    return record


def _validate_forward_log(log_text: str) -> re.Match[str]:
    optimizations = list(re.finditer(r"Optimization completed\.", log_text))
    error_terminations = list(
        re.finditer(r"^\s*Error termination via\b[^\n]*", log_text, re.MULTILINE)
    )
    normal_terminations = re.findall(
        r"Normal termination of Gaussian(?: 16)?", log_text
    )
    if len(optimizations) != EXPECTED_OPTIMIZATIONS:
        raise RuntimeError(
            f"Expected exactly {EXPECTED_OPTIMIZATIONS} completed optimizations; "
            f"found {len(optimizations)}"
        )
    if len(error_terminations) != EXPECTED_ERROR_TERMINATIONS:
        raise RuntimeError(
            "Expected exactly one Gaussian error-termination event; "
            f"found {len(error_terminations)}"
        )
    if normal_terminations:
        raise RuntimeError("Original forward scan must have no normal terminations")
    if error_terminations[0].start() < optimizations[-1].end():
        raise RuntimeError("Last optimization completed after the error termination")
    return optimizations[-1]


def _parse_gaussian_input(text: str) -> tuple[list[int], np.ndarray, int, int]:
    lines = text.splitlines()
    charge_line = next(
        (
            index
            for index, line in enumerate(lines)
            if re.fullmatch(r"\s*-?\d+\s+\d+\s*", line)
        ),
        None,
    )
    if charge_line is None:
        raise ValueError("Gaussian input has no charge and multiplicity line")
    charge, multiplicity = (int(value) for value in lines[charge_line].split())
    elements: list[int] = []
    coordinates: list[list[float]] = []
    periodic = Chem.GetPeriodicTable()
    for line in lines[charge_line + 1 :]:
        if not line.strip():
            break
        fields = line.split()
        if len(fields) != 4:
            raise ValueError("Malformed Cartesian geometry in Gaussian input")
        try:
            atomic_number = periodic.GetAtomicNumber(fields[0])
            point = [float(value) for value in fields[1:]]
        except (RuntimeError, ValueError) as exc:
            raise ValueError("Invalid atom or coordinate in Gaussian input") from exc
        if atomic_number <= 0 or not all(math.isfinite(value) for value in point):
            raise ValueError("Gaussian input geometry must contain finite known atoms")
        elements.append(atomic_number)
        coordinates.append(point)
    xyz = np.asarray(coordinates, dtype=float)
    if not elements or xyz.shape != (len(elements), 3):
        raise ValueError("Gaussian input contains no complete Cartesian geometry")
    return elements, xyz, charge, multiplicity


def _validate_source_route(text: str) -> None:
    route_lines = [
        line.strip() for line in text.splitlines() if line.lstrip().startswith("#")
    ]
    route = " ".join(route_lines)
    required = (
        "wB97XD/def2SVP",
        "SCRF=(SMD,Solvent=Water)",
        "Opt=(ModRedundant,Tight,MaxCycles=200)",
        "SCF=(Tight,XQC)",
        "NoSymm",
    )
    if not route or any(item not in route for item in required):
        raise ValueError("Source forward GJF route differs from the parent scan method")
    if not re.search(
        r"^D\s+1\s+2\s+4\s+5\s+S\s+24\s+15(?:\.0+)?\s*$", text, re.MULTILINE
    ):
        raise ValueError("Source forward GJF must schedule D 1 2 4 5, 24 x +15 degrees")


def _render_scan(
    elements: list[int],
    xyz: np.ndarray,
    charge: int,
    multiplicity: int,
    direction: str,
    steps: int,
    increment: float,
) -> str:
    periodic = Chem.GetPeriodicTable()
    lines = [
        "%nprocshared=4",
        "%mem=8GB",
        f"%chk=amide_omega_recovery_{direction}.chk",
        ROUTE,
        "",
        f"PA6 amide_omega recovery {direction}",
        "",
        f"{charge} {multiplicity}",
    ]
    lines.extend(
        f"{periodic.GetElementSymbol(number):2s} {x: .10f} {y: .10f} {z: .10f}"
        for number, (x, y, z) in zip(elements, xyz, strict=True)
    )
    lines.extend(["", f"D 1 2 4 5 S {steps} {increment:.6f}", ""])
    return "\n".join(lines)


def _recovery_manifest(
    source_manifest: Path,
    forward_gjf: Path,
    forward_log: Path,
    parent_record: dict[str, Any],
    forward_text: str,
    reverse_text: str,
) -> dict[str, Any]:
    source_hashes = {
        "parent_manifest_sha256": _sha256(source_manifest),
        "forward_gjf_sha256": _sha256(forward_gjf),
        "forward_log_sha256": _sha256(forward_log),
    }
    # Gaussian's ``S n step`` schedule optimizes the starting geometry and
    # then performs ``n`` increments, yielding ``n + 1`` completed points.
    forward_point_count = FORWARD_STEPS + 1
    reverse_point_count = REVERSE_STEPS + 1
    return {
        "schema_version": 1,
        "preparation_status": "passed",
        "scope": "PA6 internal_neutral amide_omega recovery scans",
        "method": METHOD,
        "route": ROUTE,
        "source_paths": {
            "parent_manifest": str(source_manifest),
            "forward_gjf": str(forward_gjf),
            "forward_log": str(forward_log),
        },
        "source_hashes": source_hashes,
        "parent": {
            "fragment": "internal_neutral",
            "label": "amide_omega",
            "atom_indices_zero_based": parent_record["atom_indices_zero_based"],
            "atom_indices_gaussian_one_based": [1, 2, 4, 5],
            "steps_per_direction": EXPECTED_STEPS,
            "forward_increment_degrees": EXPECTED_INCREMENT,
            "reverse_increment_degrees": -EXPECTED_INCREMENT,
            "optimization_completed_count": EXPECTED_OPTIMIZATIONS,
            "error_termination_event_count": EXPECTED_ERROR_TERMINATIONS,
            "normal_termination_count": 0,
            "schedule_total_degrees": EXPECTED_STEPS * EXPECTED_INCREMENT,
        },
        "recovery_schedules": {
            "forward": {
                "starts_from": "final orientation before last completed optimization",
                "steps": FORWARD_STEPS,
                "increment_degrees": FORWARD_INCREMENT,
                "swept_degrees": FORWARD_STEPS * FORWARD_INCREMENT,
                "expected_points": forward_point_count,
                "input_sha256": _sha256_text(forward_text),
            },
            "reverse": {
                "starts_from": "original forward GJF seed geometry",
                "steps": REVERSE_STEPS,
                "increment_degrees": REVERSE_INCREMENT,
                "swept_degrees": REVERSE_STEPS * abs(REVERSE_INCREMENT),
                "expected_points": reverse_point_count,
                "input_sha256": _sha256_text(reverse_text),
            },
            "total_expected_points": forward_point_count + reverse_point_count,
        },
        "expected_outputs": [
            "forward/scan.gjf",
            "reverse/scan.gjf",
            "manifest.json",
        ],
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()
