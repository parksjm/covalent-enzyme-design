"""Retain accepted relay coordinates and resample their atom-labeled QM path."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7

from .native_sander_cli import source_sha256
from .relay_case import RelayCase, load_case


@dataclass
class RelayTrace:
    """Validated dense ES-to-acyl active-coordinate trace and provenance."""

    model: RelayCase
    paths: list[Path]
    active_coordinates: np.ndarray
    audits: list[dict[str, Any]]
    sha256: dict[str, str]


def proton_audit(model: RelayCase, xyz: np.ndarray) -> dict[str, Any]:
    """Reject compressed/missing real QM H while preserving every atom label.

    Args:
        model: Explicit source case, QM real-H identities and mapping.
        xyz: Full-system coordinates in angstrom.

    Returns:
        All 54 nearest H parents, proton distance extrema and acid inventory.

    Raises:
        ValueError: If H-heavy distances leave [0.65,1.45] A, H-H falls below
            0.65 A, or the retained ASH proton changes its declared parent.
    """
    heavy = np.array(
        [i for i in model.qm if model.topology.atoms[i].atomic_number != 1]
    )
    delta = xyz[model.hydrogens, None] - xyz[heavy][None]
    distances = np.linalg.norm(delta, axis=2)
    closest = np.argmin(distances, axis=1)
    nearest = distances[np.arange(len(closest)), closest]
    hh = np.linalg.norm(xyz[model.hydrogens, None] - xyz[model.hydrogens][None], axis=2)
    np.fill_diagonal(hh, np.inf)
    m = model.atom_map
    acid = float(np.linalg.norm(xyz[m["ash_proton"] - 1] - xyz[m["ash_oxygen"] - 1]))
    if nearest.min() < 0.65 or nearest.max() > 1.45 or hh.min() < 0.65:
        raise ValueError(
            f"Compressed/unassigned proton: H-heavy range {nearest.min():.4f}..{nearest.max():.4f}, H-H min {hh.min():.4f} A"
        )
    if not 0.8 <= acid <= 1.3:
        raise ValueError("Retained ASH274 proton parent geometry changed")
    return dict(
        real_H_count=len(model.hydrogens),
        nearest_H_heavy_min_A=float(nearest.min()),
        nearest_H_heavy_max_A=float(nearest.max()),
        min_H_H_A=float(hh.min()),
        ash_HD2_OD2_A=acid,
        hydrogen_parents_1based={
            str(int(h + 1)): int(heavy[p] + 1)
            for h, p in zip(model.hydrogens, closest, strict=True)
        },
    )


def stage_trace_paths(campaign: Path) -> list[Path]:
    """Return38 converged endpoints plus the explicit window025 transition bridge.

    Args:
        campaign: Completed relay endpoint/reverse-path campaign directory.

    Returns:
        Paths preserving stage/checkpoint time ordering and original H labels.

    Raises:
        ValueError: If any stage lacks process completion or force convergence.
    """
    stages = [campaign / "02_acyl/02_unbiased_release"]
    stages += [campaign / f"03_reverse_path/window_{i:03d}" for i in range(1, 37)]
    stages += [campaign / "03_reverse_path/unbiased_ES_release"]
    paths = []
    for stage in stages:
        process = json.loads((stage / "process.json").read_text())
        acceptance = json.loads((stage / "scientific_acceptance.json").read_text())
        if (
            process["completed"] is not True
            or acceptance["force_converged"] is not True
        ):
            raise ValueError(f"Source stage not process/force accepted: {stage}")
        checkpoints = sorted(stage.glob("accepted_[0-9][0-9][0-9][0-9][0-9].ncrst"))
        if not checkpoints:
            raise ValueError(f"Missing accepted checkpoint trace: {stage}")
        if stage.name == "window_025":
            paths.extend(checkpoints)
        paths.append(stage / "optimized.ncrst")
    if not any(
        p.parent.name == "window_025" and p.name == "accepted_00025.ncrst"
        for p in paths
    ):
        raise ValueError("Explicit window025 transition trace is missing")
    return paths[::-1]


def load_trace(case_path: Path) -> RelayTrace:
    """Load fixed-exterior accepted source coordinates without renumbering atoms.

    Args:
        case_path: Validated ASP272/ASH274 campaign case JSON.

    Returns:
        Dense active-coordinate trace retaining all source identities and hashes.

    Raises:
        ValueError: If model, box, frozen exterior, proton or stage gates fail.
    """
    model = load_case(case_path)
    qm_set = set(model.qm.tolist())
    links = sum(
        (b.atom1.idx in qm_set) != (b.atom2.idx in qm_set) for b in model.topology.bonds
    )
    if links != 8 or int(model.active.sum()) != 4662:
        raise ValueError("Expected unchanged eight-link / 4662-active model")
    water = model.topology.atoms[model.atom_map["water_o"] - 1]
    if water.residue.idx + 1 != 1576 or water.residue.name != "WAT":
        raise ValueError("Expected original explicit QM water1576")
    paths = stage_trace_paths(case_path.parent)
    expected_end = case_path.parent / "02_acyl/02_unbiased_release/optimized.ncrst"
    if paths[-1] != expected_end:
        raise ValueError("Final source must be exact optimized unbiased acyl endpoint")
    coordinates, audits, hashes = [], [], {}
    hashes[str(case_path)] = source_sha256(case_path)
    hashes[model.config["topology"]] = source_sha256(Path(model.config["topology"]))
    for path in paths:
        rst = Rst7.open(str(path))
        xyz = np.asarray(rst.coordinates, dtype=float).reshape(-1, 3)
        if xyz.shape != model.coordinates.shape or not np.isfinite(xyz).all():
            raise ValueError(f"Source atom count/finiteness mismatch: {path}")
        if not np.array_equal(rst.box, model.restart.box):
            raise ValueError(f"Source box changed: {path}")
        if not np.array_equal(xyz[~model.active], model.coordinates[~model.active]):
            raise ValueError(f"Frozen exterior changed: {path}")
        if (path.parent / "one_cycle.mdin").read_text() != model.config["mdin"]:
            raise ValueError("Source Hamiltonian/native input differs from fixed case")
        audits.append({"path": str(path), **proton_audit(model, xyz)})
        coordinates.append(xyz[model.active].copy())
        hashes[str(path)] = source_sha256(path)
    for audit, parent in ((audits[0], "thr_og"), (audits[-1], "leaving_n")):
        if audit["hydrogen_parents_1based"]["20318"] != model.atom_map[parent]:
            raise ValueError("Returned E-S/acyl source proton label20318 changed")
        if (
            audit["hydrogen_parents_1based"][str(model.atom_map["proton"])]
            != model.atom_map["water_o"]
        ):
            raise ValueError(
                "Original HG1 must remain on source water at both endpoints"
            )
    return RelayTrace(model, paths, np.stack(coordinates), audits, hashes)


def resample_trace(
    frames: np.ndarray,
    metric_indices: np.ndarray,
    count: int,
) -> tuple[np.ndarray, list[dict[str, Any]], dict[str, Any]]:
    """Sample a continuous labeled polyline by QM Cartesian arc length.

    Args:
        frames: Finite (frames, active_atoms, 3) coordinates in angstrom.
        metric_indices: Zero-based positions of real QM atoms within active order.
        count: Required number of samples including both endpoints.

    Returns:
        Samples, exact bracketing-frame/fraction metadata and both arc lengths.

    Raises:
        ValueError: If coordinates, metric identities, count or total arc fail.
    """
    frames, ids = np.asarray(frames, dtype=float), np.asarray(metric_indices)
    if (
        frames.ndim != 3
        or frames.shape[2] != 3
        or not np.isfinite(frames).all()
        or len(frames) < 2
        or type(count) is not int
        or count < 2
    ):
        raise ValueError("Expected finite polyline and at least two samples")
    if (
        ids.ndim != 1
        or not len(ids)
        or ids.dtype.kind not in "iu"
        or ids.min() < 0
        or ids.max() >= frames.shape[1]
        or len(set(ids.tolist())) != len(ids)
    ):
        raise ValueError("Invalid distinct QM metric indices")
    steps = np.linalg.norm(
        np.diff(frames[:, ids], axis=0).reshape(len(frames) - 1, -1), axis=1
    )
    active_steps = np.linalg.norm(
        np.diff(frames, axis=0).reshape(len(frames) - 1, -1), axis=1
    )
    arc = np.r_[0.0, np.cumsum(steps)]
    if arc[-1] <= 0:
        raise ValueError("Zero QM-coordinate path arc")
    selected, metadata = [], []
    for target in np.linspace(0, arc[-1], count):
        left = min(int(np.searchsorted(arc, target, side="right")) - 1, len(frames) - 2)
        left = max(left, 0)
        while left < len(frames) - 2 and arc[left + 1] == arc[left]:
            left += 1
        fraction = (
            float((target - arc[left]) / (arc[left + 1] - arc[left]))
            if arc[left + 1] > arc[left]
            else 0.0
        )
        selected.append((1 - fraction) * frames[left] + fraction * frames[left + 1])
        metadata.append(
            dict(
                left_frame=left,
                right_frame=left + 1,
                fraction=fraction,
                qm_arc_A=float(target),
            )
        )
    selected[0], selected[-1] = frames[0].copy(), frames[-1].copy()
    return (
        np.stack(selected),
        metadata,
        dict(
            qm_arc_A=float(arc[-1]),
            active_arc_A=float(active_steps.sum()),
            source_adjacent_qm_arcs_A=steps.tolist(),
            source_adjacent_active_arcs_A=active_steps.tolist(),
        ),
    )
