"""Extract canonical protein CA coordinates and calculate pairwise diversity."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from itertools import combinations
from typing import Any

import numpy as np

from src.nyl12_covalent.catalytic_alignment_geometry import fit_coordinates

FOCUSED_ALPHA = (
    50,
    51,
    52,
    73,
    76,
    101,
    111,
    112,
    114,
    115,
    116,
    151,
    152,
    153,
    185,
    187,
    188,
)
ACTIVE_ALPHA = frozenset((*FOCUSED_ALPHA, 113, 156, 186))
ACTIVE_BETA = frozenset((1, 39, 40, 41, 42))
CHAIN_LENGTHS = {**dict.fromkeys("ABCD", 232), **dict.fromkeys("EFGH", 96)}
THREE_TO_ONE = dict(
    zip(
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
        "A R N D C Q E G H I L K M F P S T W Y V".split(),
        strict=True,
    )
)
_active_indices: list[int] = []
_offset = 0
for _chain, _length in CHAIN_LENGTHS.items():
    _positions = ACTIVE_ALPHA if _chain in "ABCD" else ACTIVE_BETA
    _active_indices.extend(_offset + position - 1 for position in _positions)
    _offset += _length
ACTIVE_INDICES = np.asarray(sorted(_active_indices), dtype=int)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of an existing input file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract_ca(path: Path, alpha: str, beta: str) -> tuple[np.ndarray, dict[str, str]]:
    """Read exactly one CA per canonical A-H position and verify sequence."""
    residues: dict[str, dict[int, tuple[str, np.ndarray]]] = {
        chain: {} for chain in CHAIN_LENGTHS
    }
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if line[:6] != "ATOM  " or line[12:16].strip() != "CA":
            continue
        chain = line[21]
        if chain not in CHAIN_LENGTHS:
            continue
        if line[16].strip() or line[26].strip():
            raise ValueError(f"Alternate/insertion CA at {path}:{line_number}")
        try:
            position = int(line[22:26])
            name = THREE_TO_ONE[line[17:20].strip()]
            xyz = np.asarray([float(line[i : i + 8]) for i in (30, 38, 46)])
        except (ValueError, KeyError) as exc:
            raise ValueError(f"Invalid canonical CA at {path}:{line_number}") from exc
        if position in residues[chain]:
            raise ValueError(f"Duplicate CA {chain}:{position} in {path}")
        if not np.isfinite(xyz).all():
            raise ValueError(f"Nonfinite CA coordinate at {path}:{line_number}")
        residues[chain][position] = name, xyz

    sequences: dict[str, str] = {}
    coordinates: list[np.ndarray] = []
    for chain, length in CHAIN_LENGTHS.items():
        observed = residues[chain]
        if set(observed) != set(range(1, length + 1)):
            raise ValueError(
                f"Incomplete canonical CA positions for chain {chain} in {path}"
            )
        sequence = "".join(observed[pos][0] for pos in range(1, length + 1))
        if sequence != (alpha if chain in "ABCD" else beta):
            raise ValueError(f"Chain {chain} sequence differs from manifest in {path}")
        sequences[chain] = sequence
        coordinates.extend(observed[pos][1] for pos in range(1, length + 1))
    result = np.asarray(coordinates, dtype=float)
    if result.shape != (1312, 3) or len(ACTIVE_INDICES) != 100:
        raise ValueError(f"Invalid CA inventory in {path}: {result.shape}")
    if not np.isfinite(result).all():
        raise ValueError(f"Nonfinite CA coordinates in {path}")
    return result, sequences


def verify_fit_helper() -> dict[str, float]:
    """Verify known rigid-transform recovery and proper-rotation handedness."""
    mobile = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    angle = 0.731
    expected = np.asarray(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    target = mobile @ expected + [4.0, -2.0, 7.0]
    rotation, _, _, rmsd = fit_coordinates(mobile, target)
    determinant = float(np.linalg.det(rotation))
    if rmsd > 1.0e-10 or determinant < 0.999999999:
        raise ValueError("Shared fit_coordinates failed known proper-rigid recovery")
    return {"known_transform_rmsd_angstrom": rmsd, "rotation_determinant": determinant}


def load_models(boltz_root: Path) -> tuple[list[dict[str, Any]], dict[str, str]]:
    """Validate 45 model receipts and load/hashes raw and relaxed coordinates."""
    candidate_path = boltz_root / "prepared/candidate_manifest.json"
    data = json.loads(candidate_path.read_text())
    candidates, seeds = {item["id"]: item for item in data["candidates"]}, data["seeds"]
    if len(candidates) != 9 or len(seeds) != 5:
        raise ValueError("Expected WT plus eight variants and five seeds")
    status_path = boltz_root / "status.json"
    status = json.loads(status_path.read_text())
    if status.get("status") != "completed" or any(
        status.get(key) != 45
        for key in ("planned_models", "normalized_models", "relaxed_models")
    ):
        raise ValueError(f"Campaign status is not completed at 45 models: {status}")

    records: list[dict[str, Any]] = []
    hashes: dict[str, str] = {}
    for variant, candidate in candidates.items():
        alpha, beta = candidate["alpha_sequence"], candidate["beta_sequence"]
        if (len(alpha), len(beta)) != (232, 96):
            raise ValueError(f"Unexpected family sequence length for {variant}")
        for index, seed in enumerate(seeds, start=1):
            label = f"seed_{index:02d}"
            model_root = boltz_root / "models" / variant / label
            raw = model_root / "raw/normalized.pdb"
            request_path = model_root / "request.json"
            acceptance_path = model_root / "relaxation/acceptance.json"
            relaxed = model_root / "relaxation/relaxed.pdb"
            for path in (raw, request_path, acceptance_path, relaxed):
                if not path.is_file():
                    raise FileNotFoundError(path)
                hashes[str(path.resolve())] = sha256(path)
            request = json.loads(request_path.read_text())
            if (request.get("variant"), request.get("seed"), request.get("sha256")) != (
                variant,
                seed,
                hashes[str(raw.resolve())],
            ):
                raise ValueError(
                    f"Raw model lineage/hash mismatch for {variant}/{label}"
                )
            acceptance = json.loads(acceptance_path.read_text())
            if (
                acceptance.get("status") != "minimization_complete"
                or acceptance.get("process_exit_code") != 0
            ):
                raise ValueError(f"Relaxation process incomplete for {variant}/{label}")
            converged, gradient = (
                acceptance.get("force_converged"),
                acceptance.get("final_rms_gradient"),
            )
            if not isinstance(converged, bool) or not isinstance(
                gradient, (int, float)
            ):
                raise ValueError(
                    f"Missing force-convergence data for {variant}/{label}"
                )
            raw_xyz, raw_seq = extract_ca(raw, alpha, beta)
            relaxed_xyz, relaxed_seq = extract_ca(relaxed, alpha, beta)
            if raw_seq != relaxed_seq:
                raise ValueError(
                    f"Raw/relaxed CA identity differs for {variant}/{label}"
                )
            records.append(
                {
                    "variant": variant,
                    "seed": seed,
                    "force_converged": converged,
                    "final_rms_gradient": float(gradient),
                    "raw_path": raw,
                    "relaxed_path": relaxed,
                    "raw_xyz": raw_xyz,
                    "relaxed_xyz": relaxed_xyz,
                    "raw_sha256": hashes[str(raw.resolve())],
                    "relaxed_sha256": hashes[str(relaxed.resolve())],
                }
            )
    if len(records) != 45:
        raise ValueError(f"Expected 45 models, found {len(records)}")
    return records, hashes


def pairwise_rmsds(models: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Calculate symmetric global-fit RMSDs for ten pairs per variant/stage."""
    rows: list[dict[str, Any]] = []
    for variant in sorted({model["variant"] for model in models}):
        group = sorted(
            (m for m in models if m["variant"] == variant), key=lambda m: m["seed"]
        )
        if len(group) != 5:
            raise ValueError(f"Variant {variant} does not have five models")
        for stage in ("raw", "relaxed"):
            for first, second in combinations(group, 2):
                left, right = first[f"{stage}_xyz"], second[f"{stage}_xyz"]
                rotation, center_left, center_right, global_rmsd = fit_coordinates(
                    left, right
                )
                reverse_rmsd = fit_coordinates(right, left)[3]
                if abs(global_rmsd - reverse_rmsd) > 1.0e-8:
                    raise ValueError(
                        f"Kabsch pair symmetry failed for {variant}/{stage}"
                    )
                aligned = (left - center_left) @ rotation + center_right
                delta = aligned[ACTIVE_INDICES] - right[ACTIVE_INDICES]
                active_rmsd = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
                if not np.isfinite((global_rmsd, active_rmsd)).all():
                    raise ValueError(f"Nonfinite RMSD for {variant}/{stage}")
                rows.append(
                    {
                        "variant": variant,
                        "stage": stage,
                        "seed_left": first["seed"],
                        "seed_right": second["seed"],
                        "force_converged_left": first["force_converged"]
                        if stage == "relaxed"
                        else "",
                        "force_converged_right": second["force_converged"]
                        if stage == "relaxed"
                        else "",
                        "both_force_converged": (
                            first["force_converged"] and second["force_converged"]
                            if stage == "relaxed"
                            else ""
                        ),
                        "global_ca_rmsd_angstrom": float(global_rmsd),
                        "active_ca_rmsd_angstrom": active_rmsd,
                    }
                )
    return rows
