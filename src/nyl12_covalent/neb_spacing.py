"""Immutable native-band reconstruction and bounded NEB spacing predictions."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7

from .native_sander_cli import source_sha256
from .neb_xyz_snapshot import reconstruct_xyz_band
from .relay_case import RelayCase, load_case


def load_spacing_source(
    config_path: Path,
    snapshot: Path,
) -> tuple[RelayCase, Rst7, np.ndarray, dict[str, Any]]:
    """Reconstruct24 images solely from a retained verified native XYZ update.

    Args:
        config_path: Original running NEB configuration, read-only.
        snapshot: Retained accepted monitor update containing gzip native XYZ.

    Returns:
        Explicit model, reference restart, active24-image band and provenance.

    Raises:
        ValueError: If immutable payload/fixed-input hashes or mapping gates fail.
    """
    cfg = json.loads(config_path.read_text())
    meta = json.loads((snapshot / "band_snapshot.json").read_text())
    for key in ("reference", "topology"):
        expected = cfg["reactant"] if key == "reference" else cfg["topology"]
        if (
            Path(meta[key]["path"]) != Path(expected)
            or source_sha256(Path(expected)) != meta[key]["sha256"]
        ):
            raise ValueError(f"Snapshot fixed-input identity mismatch: {key}")
    if meta["restraints"] or meta["active_mask"] != cfg["active_mask"]:
        raise ValueError("Snapshot restraint/active-mask contract mismatch")
    payload = gzip.decompress((snapshot / "source_nebpath.xyz.gz").read_bytes())
    if hashlib.sha256(payload).hexdigest() != meta["source_xyz_sha256"]:
        raise ValueError("Retained native XYZ payload hash mismatch")
    model = load_case(Path(cfg["case_path"]))
    selected = set(model.qm.tolist())
    links = sum(
        (b.atom1.idx in selected) != (b.atom2.idx in selected)
        for b in model.topology.bonds
    )
    if (
        links != 8
        or len(model.qm) != 102
        or len(model.hydrogens) != 54
        or model.config["qm_charge"] != -1
    ):
        raise ValueError("Spacing source requires102QM54H8links and charge-1")
    reference = Rst7.open(cfg["reactant"])
    full, active = reconstruct_xyz_band(
        payload, model.topology, reference, cfg["active_mask"], 24, 4662
    )
    if (
        active != meta["active_indices_0based"]
        or active != np.flatnonzero(model.active).tolist()
    ):
        raise ValueError("Snapshot active identities differ from validated case")
    coordinates = np.stack([xyz[active] for xyz in full])
    evidence = dict(
        snapshot=str(snapshot.resolve()),
        cycle=meta["report"]["cycle"],
        report=meta["report"],
        source_xyz_sha256=meta["source_xyz_sha256"],
        snapshot_metadata_sha256=source_sha256(snapshot / "band_snapshot.json"),
        topology_sha256=meta["topology"]["sha256"],
        reference_sha256=meta["reference"]["sha256"],
        config_sha256=source_sha256(config_path),
        active_atom_ids=(np.array(active) + 1).tolist(),
        numbered_native_restarts_used=False,
        native_endpoints_preserved=True,
    )
    return model, reference, coordinates, evidence


def spacing_controls(
    band: np.ndarray,
    qm_positions: np.ndarray,
    *,
    focus_weight: float = 1.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Sample13 labeled controls from24 source images by weighted QM arc.

    Args:
        band: Immutable24-image active-coordinate source band in angstrom.
        qm_positions: Real QM atom positions in active ordering.
        focus_weight: Positive sampling density multiplier on source image9→10.
            One gives uniform QM arc; larger values emphasize this single segment.

    Returns:
        Exactly13 controls preserving both endpoint arrays bitwise, plus brackets.

    Raises:
        ValueError: If coordinates, metric selection or focus density are invalid.
    """
    xyz, qm = np.asarray(band, dtype=float), np.asarray(qm_positions)
    if (
        xyz.ndim != 3
        or xyz.shape[0] != 24
        or xyz.shape[2] != 3
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Expected24 finite active-coordinate source images")
    if (
        qm.ndim != 1
        or not len(qm)
        or qm.dtype.kind not in "iu"
        or qm.min() < 0
        or qm.max() >= xyz.shape[1]
    ):
        raise ValueError("Invalid QM metric positions")
    if not np.isfinite(focus_weight) or focus_weight <= 0:
        raise ValueError("Focus density must be positive and finite")
    steps = np.linalg.norm(np.diff(xyz[:, qm], axis=0).reshape(23, -1), axis=1)
    active_steps = np.linalg.norm(np.diff(xyz, axis=0).reshape(23, -1), axis=1)
    weighted = steps.copy()
    weighted[8] *= focus_weight
    cumulative = np.r_[0.0, np.cumsum(weighted)]
    if cumulative[-1] <= 0:
        raise ValueError("Zero QM path arc")
    controls, brackets = [], []
    for target in np.linspace(0, cumulative[-1], 13):
        left = max(
            0, min(int(np.searchsorted(cumulative, target, side="right")) - 1, 22)
        )
        while left < 22 and cumulative[left + 1] == cumulative[left]:
            left += 1
        fraction = (
            float(
                (target - cumulative[left]) / (cumulative[left + 1] - cumulative[left])
            )
            if cumulative[left + 1] > cumulative[left]
            else 0.0
        )
        controls.append((1 - fraction) * xyz[left] + fraction * xyz[left + 1])
        brackets.append(
            dict(source_images_1based=[left + 1, left + 2], fraction=fraction)
        )
    controls[0], controls[-1] = xyz[0].copy(), xyz[-1].copy()
    return np.stack(controls), dict(
        focus_segment_1based=[9, 10],
        focus_weight=focus_weight,
        qm_arc_A=float(steps.sum()),
        active_arc_A=float(active_steps.sum()),
        source_qm_steps_A=steps.tolist(),
        control_brackets=brackets,
        endpoints_bitwise_preserved=True,
    )


def predict_spacing_band(
    controls: np.ndarray,
    qm_positions: np.ndarray,
    regular_images: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Predict unchanged native all-active alignment/rounded variable-size band.

    Args:
        controls: Exactly13 active-coordinate input controls.
        qm_positions: QM metric positions within active ordering.
        regular_images: Explicit regular-band count,24 or48; CI is an extra slot.

    Returns:
        Predicted regular images and exact slot/alignment/endpoint drift evidence.

    Raises:
        ValueError: If inputs or native mapping cannot retain endpoints.
    """
    xyz, qm = np.asarray(controls, dtype=float), np.asarray(qm_positions)
    if (
        xyz.ndim != 3
        or xyz.shape[0] != 13
        or xyz.shape[2] != 3
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Expected13 finite controls")
    if (
        regular_images not in (24, 48)
        or qm.ndim != 1
        or not len(qm)
        or qm.dtype.kind not in "iu"
        or qm.min() < 0
        or qm.max() >= xyz.shape[1]
    ):
        raise ValueError("Expected24/48 regular images and valid QM positions")
    aligned = xyz.copy()
    reference = xyz[0] - xyz[0].mean(axis=0)
    for i in range(1, 13):
        centered = xyz[i] - xyz[i].mean(axis=0)
        u, _, vt = np.linalg.svd(centered.T @ reference)
        if np.linalg.det(u @ vt) < 0:
            u[:, -1] *= -1
        aligned[i] = centered @ (u @ vt) + xyz[0].mean(axis=0)
    arc = np.r_[
        0.0,
        np.cumsum(
            np.linalg.norm(np.diff(aligned[:, qm], axis=0).reshape(12, -1), axis=1)
        ),
    ]
    if arc[-1] < 0.000529177:
        raise ValueError("Zero aligned QM arc")
    slots = np.zeros(regular_images + 1, dtype=int)
    for i in range(12, 0, -1):
        slots[i] = int(np.floor(1 + arc[i] / arc[-1] * (regular_images - 1) + 0.5)) - 1
    band = np.zeros((regular_images + 1, *xyz.shape[1:]))
    band[:13] = aligned
    used = np.zeros(regular_images + 1, dtype=bool)
    used[:13] = True
    dropped = []
    for i in range(12, 0, -1):
        others = [j for j in range(12, 0, -1) if j != i]
        unique = not any(slots[j] == slots[i] for j in others)
        if not unique:
            unique = not any(slots[j] == i for j in others)
            if unique:
                slots[i] = i
        if not unique:
            slots[i], used[i] = -1, False
            dropped.append(i)
            continue
        target = int(slots[i])
        if target != i:
            band[target] = band[i]
            used[target], used[i] = True, False
    if slots[12] != regular_images - 1 or not used[0] or not used[regular_images - 1]:
        raise ValueError("Native rounded mapping lost exact final endpoint control")
    anchors = np.flatnonzero(used[:regular_images])
    for left, right in zip(anchors[:-1], anchors[1:], strict=True):
        for i in range(int(left) + 1, int(right)):
            fraction = (i - left) / (right - left)
            band[i] = (1 - fraction) * band[left] + fraction * band[right]
    if not np.array_equal(band[regular_images - 1], aligned[-1]):
        raise ValueError("Native final slot differs from aligned exact endpoint")
    shift = np.linalg.norm(aligned - xyz, axis=2)
    drift = float(np.linalg.norm(band[regular_images - 1] - xyz[-1], axis=1).max())
    return band[:regular_images], dict(
        regular_images=regular_images,
        neb_nimage=regular_images + 1,
        slots_1based=[1]
        + [int(slots[i] + 1) if slots[i] >= 0 else None for i in range(1, 13)],
        dropped_controls_zero_based=dropped,
        max_active_alignment_shift_A=float(shift.max()),
        max_qm_alignment_shift_A=float(shift[:, qm].max()),
        endpoint_max_drift_A=drift,
        endpoint_within_1e_minus6_A=drift <= 1e-6,
        actual_native_startup_verification_required=True,
    )
