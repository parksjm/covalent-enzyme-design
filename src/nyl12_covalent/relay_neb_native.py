"""Predict installed DL-FIND Cartesian NEB seeding and render supported input."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any

import numpy as np

from .native_sander_cli import ANGSTROM_PER_BOHR


def native_initial_band(
    controls: np.ndarray,
    qm_positions: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Mirror all-active alignment and Fortran rounded 24-slot initialization.

    Args:
        controls: Exactly 13 full-active Cartesian coordinate control structures.
        qm_positions: QM atom positions within the original active atom order.

    Returns:
        Predicted 24 regular Cartesian images and alignment/slot-map evidence.
        Actual native path.xyz remains authoritative at startup.

    Raises:
        ValueError: If controls are invalid, degenerate or native drops endpoints.
    """
    xyz = np.asarray(controls, dtype=float)
    qm = np.asarray(qm_positions)
    if (
        xyz.ndim != 3
        or xyz.shape[0] != 13
        or xyz.shape[2] != 3
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Expected exactly13 finite active-coordinate controls")
    if (
        qm.ndim != 1
        or not len(qm)
        or qm.dtype.kind not in "iu"
        or qm.min() < 0
        or qm.max() >= xyz.shape[1]
    ):
        raise ValueError("Invalid QM metric positions")
    aligned = xyz.copy()
    reference = xyz[0] - xyz[0].mean(axis=0)
    shifts = []
    for i in range(1, 13):
        centered = xyz[i] - xyz[i].mean(axis=0)
        left, _, right = np.linalg.svd(centered.T @ reference)
        rotation = left @ right
        if np.linalg.det(rotation) < 0:
            left[:, -1] *= -1
            rotation = left @ right
        aligned[i] = centered @ rotation + xyz[0].mean(axis=0)
        shift = np.linalg.norm(aligned[i] - xyz[i], axis=1)
        shifts.append(
            dict(
                control=i,
                active_max_shift_A=float(shift.max()),
                qm_max_shift_A=float(shift[qm].max()),
            )
        )
    arc = np.r_[
        0.0,
        np.cumsum(
            np.linalg.norm(np.diff(aligned[:, qm], axis=0).reshape(12, -1), axis=1)
        ),
    ]
    if arc[-1] < 1e-3 * ANGSTROM_PER_BOHR:
        raise ValueError("Zero aligned QM path arc")
    # Fortran NINT is half away from zero; all these slot values are positive.
    mapped = np.zeros(25, dtype=int)
    for i in range(12, 0, -1):
        mapped[i] = int(np.floor(1 + arc[i] / arc[-1] * 23 + 0.5)) - 1
    band = np.zeros((25, *xyz.shape[1:]))
    band[:13] = aligned
    used = np.zeros(25, dtype=bool)
    used[:13] = True
    dropped = []
    for i in range(12, 0, -1):
        others = [j for j in range(12, 0, -1) if j != i]
        unique = not any(mapped[j] == mapped[i] for j in others)
        if not unique:
            unique = not any(mapped[j] == i for j in others)
            if unique:
                mapped[i] = i
        if not unique:
            mapped[i], used[i] = -1, False
            dropped.append(i)
            continue
        target = int(mapped[i])
        if target != i:
            band[target] = band[i]
            used[target], used[i] = True, False
    if dropped:
        raise ValueError(
            f"Native initialization would drop source control(s): {dropped}"
        )
    if mapped[12] != 23:
        raise ValueError("Native initialization lost exact acyl endpoint control")
    if not used[0] or not used[23]:
        raise ValueError("Native seed slot mapping lost a regular endpoint")
    anchors = np.flatnonzero(used[:24])
    for left, right in zip(anchors[:-1], anchors[1:], strict=True):
        for i in range(int(left) + 1, int(right)):
            fraction = (i - left) / (right - left)
            band[i] = (1 - fraction) * band[left] + fraction * band[right]
    if not np.array_equal(band[23], aligned[-1]):
        raise ValueError(
            "Native slot24 differs from exact aligned final source control"
        )
    return band[:24], dict(
        algorithm="dlf_cartesian_align all-active unweighted proper Kabsch; dlf_neb_definepath QM weights/NINT slot map; Cartesian linear gaps",
        source_controls_to_regular_slots_1based=[1]
        + [int(mapped[i] + 1) if mapped[i] >= 0 else None for i in range(1, 13)],
        dropped_controls_zero_based=dropped,
        alignment_shifts=shifts,
        max_active_alignment_shift_A=max(row["active_max_shift_A"] for row in shifts),
        max_qm_alignment_shift_A=max(row["qm_max_shift_A"] for row in shifts),
        aligned_control_qm_arc_A=float(arc[-1]),
        native_slot24_exact_aligned_last_control=True,
        actual_native_path_xyz_required=True,
        numbered_restart_files_authoritative=False,
    )


def native_neb_input(base: str, images: list[Path]) -> str:
    """Preserve the current Hamiltonian and supply13 controls to native NEB.

    Args:
        base: Current case full one-cycle native input text.
        images: Reference plus exactly12 coords2 full-system restarts.

    Returns:
        Supported frozen-endpoint24regular-image native input; no RC restraints.

    Raises:
        ValueError: If source constraints, masks or control counts are invalid.
    """
    if len(images) != 13:
        raise ValueError(
            "Installed Sander requires13 controls: reference plus12coords2"
        )
    prefix, separator, _ = base.partition("&dlfind")
    if not separator or re.search(r"DISANG|&rst|restraintmask", base, re.I):
        raise ValueError("Expected unrestrained current-case native base input")
    for key in ("nmropt", "ntr"):
        matches = re.findall(rf"\b{key}\s*=\s*(\d+)", prefix, re.I)
        if matches != ["0"]:
            raise ValueError(f"Expected explicit {key}=0")
    masks = {
        key: re.findall(rf"\b{key}\s*=\s*'([^']+)'", base, re.I)
        for key in ("qmmask", "active")
    }
    if any(len(values) != 1 for values in masks.values()):
        raise ValueError("Expected one current QM and active mask")
    prefix, count = re.subn(r"\bmaxcyc\s*=\s*\d+", "maxcyc=100000", prefix, flags=re.I)
    if count != 1:
        raise ValueError("Expected one Amber cycle budget")
    controls = (
        "&dlfind\n neb='FROZEN', neb_nimage=25, nebk=0.01,\n"
        " neb_climb_test=0.0, neb_freeze_test=0.0,\n"
        f" refcrd='{images[0]}', wtmask='{masks['qmmask'][0]}',\n"
        + "".join(f" coords2({i})='{path}',\n" for i, path in enumerate(images[1:], 1))
        + " optalg='LBFGS', crdrep='CART', dimer=.false.,\n"
        + f" active='{masks['active'][0]}', tol=1.d-4, tole=1.d-5,\n"
        " maxcycle=16000, maxstep=0.05, verbosity=4,\n/\n"
    )
    return prefix + controls
