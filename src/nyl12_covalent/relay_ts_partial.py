"""Bounded partial-QM curvature screen on the unchanged relay QM/MM potential."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from pathlib import Path
import re

import numpy as np
from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask

from .stationary_hessian import worker_column
from .stationary_native import NativeState, digest, write_json

LOGGER = logging.getLogger(__name__)


def lift_mode(
    columns: np.ndarray, vectors: np.ndarray, inverse_masses: np.ndarray
) -> np.ndarray:
    """Lift a partial mass-weighted eigenvector into full active Cartesian space.

    Args:
        columns: Full active-space Cartesian indices used by the partial block.
        vectors: One normalized partial mass-weighted mode.
        inverse_masses: Full-space inverse square-root atomic masses, repeated xyz.

    Returns:
        Unit Cartesian displacement across every original active coordinate.
    """
    columns = np.asarray(columns)
    if (
        columns.ndim != 1
        or vectors.shape != columns.shape
        or len(np.unique(columns)) != len(columns)
        or not len(columns)
        or np.min(columns) < 0
        or np.max(columns) >= len(inverse_masses)
        or not np.isfinite(vectors).all()
        or not np.isfinite(inverse_masses).all()
        or np.any(inverse_masses <= 0)
    ):
        raise ValueError("Invalid partial mode inventory")
    result = np.zeros(len(inverse_masses))
    result[columns] = vectors * inverse_masses[columns]
    length = np.linalg.norm(result)
    if length <= 1e-14:
        raise ValueError("Degenerate partial mode")
    return result / length


def partial_screen(
    model: Path, restart: Path, output: Path, workers: int = 8, step: float = 0.000125
) -> dict[str, object]:
    """Compute real-QM block curvature and probe a lifted negative direction.

    All force calls retain the original Hamiltonian, full topology, charge,
    proton inventory and original fixed exterior. Partial curvature is a
    chemical diagnostic and mode seed; it does not establish full-space index.

    Args:
        model: Immutable original model with calibrated native precision library.
        restart: Full-system candidate coordinates.
        output: Fresh persistent derivative/evidence directory.
        workers: Maximum single-thread native worker processes.
        step: Central Cartesian finite-difference step in Angstrom.

    Returns:
        Partial curvature and step-checked full-gradient directional evidence.
    """
    if not 1 <= workers <= 8 or not 0 < step <= 0.001:
        raise ValueError("Invalid derivative resource or step settings")
    output.mkdir(parents=True, exist_ok=False)
    state = NativeState(model, restart)
    match = re.findall(r"\bqmmask\s*=\s*'([^']+)'", state.mdin_source, re.I)
    if len(match) != 1:
        raise ValueError("Require exactly one real-QM mask")
    parm = AmberParm(str(state.topology))
    qm = np.asarray(AmberMask(parm, match[0]).Selection(), dtype=bool)
    if qm.sum() != 98 or np.any(qm & ~state.active):
        raise ValueError(
            "Require unchanged 98 real QM atoms within original active shell"
        )
    atoms = np.flatnonzero(qm[state.active])
    columns = (3 * atoms[:, None] + np.arange(3)).reshape(-1)
    np.save(output / "partial_columns.npy", columns, allow_pickle=False)
    write_json(
        output / "specification.json",
        {
            "restart": str(restart),
            "restart_sha256": digest(restart),
            "topology_sha256": digest(state.topology),
            "qm_atoms": 98,
            "active_atoms": 4165,
            "dimension": len(columns),
            "step_A": step,
            "workers": workers,
            "full_active_index_verified": False,
            "scope": "Partial real-QM block; full-system embedding and analytic gradients",
        },
    )
    energy, gradient = state.gradient(state.xyz, output / "center")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(worker_column, (model, output, restart, int(c), step))
            for c in columns
        ]
        for count, future in enumerate(as_completed(futures), 1):
            future.result()
            if count == 1 or count % 20 == 0 or count == len(columns):
                write_json(
                    output / "progress.json",
                    {"completed_columns": count, "total_columns": len(columns)},
                )
                LOGGER.info("Partial QM Hessian: %s/%s", count, len(columns))
    raw = np.column_stack(
        [
            np.load(output / "columns" / f"column_{c:05d}" / "column.npy")[columns]
            for c in columns
        ]
    )
    symmetric = (raw + raw.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    np.save(output / "hessian_raw.npy", raw, allow_pickle=False)
    np.save(output / "eigenvalues.npy", values, allow_pickle=False)
    np.save(output / "eigenvectors.npy", vectors, allow_pickle=False)
    response = {
        "partial_lowest_eigenvalues": values[:10].tolist(),
        "partial_negative_count": int(np.count_nonzero(values < 0)),
        "antisymmetry_max_abs": float(np.max(np.abs(raw - raw.T))),
        "energy_hartree": energy,
        "maximum_gradient_Eh_per_bohr": float(np.max(np.abs(gradient)) * 0.529177),
        "full_active_index_verified": False,
        "partial_negative_direction_verified": False,
        "interpretation": "Partial QM block only; no full-shell index or basin claim.",
    }
    write_json(output / "scientific_acceptance.json", response)
    if values[0] >= 0:
        raise RuntimeError("Partial QM block has no negative mode; preserve diagnostic")
    direction = lift_mode(columns, vectors[:, 0], state.inv_sqrt)
    np.save(output / "lifted_cartesian_direction.npy", direction, allow_pickle=False)
    curvatures = []
    slopes = []
    for index, delta in enumerate((step, step / 2)):
        derivatives = []
        for sign, label in ((1, "plus"), (-1, "minus")):
            xyz = state.xyz.copy()
            xyz[state.active] += sign * delta * direction.reshape((-1, 3))
            _, grad = state.gradient(xyz, output / f"direction_{index}" / label)
            derivatives.append(float(np.dot(grad.reshape(-1), direction)))
        curvatures.append((derivatives[0] - derivatives[1]) / (2 * delta))
        slopes.append(derivatives)
    relative = abs(curvatures[0] - curvatures[1]) / max(abs(curvatures[1]), 1e-14)
    passed = bool(max(curvatures) < 0 and relative <= 1e-4)
    response.update(
        {
            "directional_curvatures_Eh_per_A2": curvatures,
            "directional_slopes_Eh_per_A": slopes,
            "relative_half_step_difference": relative,
            "partial_negative_direction_verified": passed,
            "mode_source": str(output / "lifted_cartesian_direction.npy"),
        }
    )
    write_json(output / "scientific_acceptance.json", response)
    if not passed:
        raise RuntimeError("Lifted negative-direction half-step verification failed")
    return response
