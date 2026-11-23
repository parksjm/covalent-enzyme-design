"""Restartable complete active-region finite-difference Hessians and inertia."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import logging
import os
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh

from nyl12_covalent.stationary_native import (
    NativeState,
    compute_column,
    digest,
    write_json,
)

LOGGER = logging.getLogger(__name__)


def worker_column(args: tuple[Path, Path, Path, int, float]) -> str:
    """Compute one column using one lower-priority CPU worker."""
    if os.getpriority(os.PRIO_PROCESS, 0) < 5:
        os.nice(5 - os.getpriority(os.PRIO_PROCESS, 0))
    return compute_column(*args)


def validate_inertia(
    values: np.ndarray, residuals: np.ndarray, expected_negative: int, error: float
) -> dict[str, object]:
    """Accept inertia only when all signs clear the stated uncertainty gate.

    Args:
        values: All mass-weighted eigenvalues, sorted ascending.
        residuals: Full-space residual norms for saved lowest eigenvectors.
        expected_negative: Required negative eigenvalue count.
        error: Absolute numerical sign-uncertainty threshold.

    Returns:
        Signed counts and a conservative numerical acceptance decision.
    """
    values = np.asarray(values)
    if values.ndim != 1 or not np.isfinite(values).all() or error <= 0:
        raise ValueError("Invalid inertia input")
    negative = int(np.sum(values < -error))
    unresolved = int(np.sum(np.abs(values) <= error))
    return {
        "passed": bool(
            negative == expected_negative
            and unresolved == 0
            and np.max(residuals) <= 1e-8
        ),
        "negative": negative,
        "unresolved": unresolved,
        "positive": int(np.sum(values > error)),
        "expected_negative": expected_negative,
        "sign_threshold_hartree_per_angstrom2_amu": error,
        "interpretation": "Complete declared active-region numerical Hessian inertia; finite-difference truncation is qualified by recorded step probes, not an analytic interval proof.",
    }


def full_hessian(
    root: Path,
    restart: Path,
    stage: Path,
    workers: int,
    step: float,
    expected_negative: int,
) -> dict[str, object]:
    """Compute all active Cartesian columns and classify complete inertia.

    Args:
        root: Immutable campaign directory with inputs and isolated library.
        restart: Stationary-point coordinates being tested.
        stage: Persistent derivative and matrix checkpoint directory.
        workers: Maximum native single-thread workers.
        step: Cartesian central-difference displacement in Angstrom.
        expected_negative: One for saddle; zero for endpoint minimum.

    Returns:
        Full-space numerical inertia acceptance and diagnostics.
    """
    stage.mkdir(exist_ok=True)
    state = NativeState(root, restart)
    dimension = len(state.indices) * 3
    specification = {
        "restart_sha256": digest(restart),
        "step": step,
        "dimension": dimension,
        "workers": workers,
        "library_sha256": digest(root / "library_build/libdlfind.so"),
    }
    manifest = stage / "specification.json"
    if manifest.exists() and json.loads(manifest.read_text()) != specification:
        raise ValueError("Hessian continuation specification changed")
    write_json(manifest, specification)
    started = time.monotonic()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(worker_column, (root, stage, restart, col, step))
            for col in range(dimension)
        ]
        try:
            for future in as_completed(futures):
                future.result()
                done += 1
                if done == 1 or done % 20 == 0 or done == dimension:
                    elapsed = time.monotonic() - started
                    write_json(
                        stage / "progress.json",
                        {
                            "completed_columns": done,
                            "total_columns": dimension,
                            "elapsed_seconds": elapsed,
                            "columns_per_hour": done * 3600 / elapsed,
                            "remaining_hours_at_current_rate": (dimension - done)
                            * elapsed
                            / done
                            / 3600,
                        },
                    )
                    LOGGER.info(
                        "%s: completed %s/%s columns", stage.name, done, dimension
                    )
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    raw = np.lib.format.open_memmap(
        stage / "hessian_raw.npy",
        mode="w+",
        dtype="float64",
        shape=(dimension, dimension),
    )
    for col in range(dimension):
        raw[:, col] = np.load(stage / "columns" / f"column_{col:05d}" / "column.npy")
    raw.flush()
    symmetric = np.asarray((raw + raw.T) / 2)
    antisymmetric = float(np.max(np.sum(np.abs(raw - raw.T), axis=1)) / 2)
    np.save(stage / "hessian_symmetric.npy", symmetric, allow_pickle=False)
    values = eigh(symmetric, eigvals_only=True, check_finite=True, driver="evr")
    low_values, modes = eigh(symmetric, subset_by_index=[0, 5], driver="evr")
    residuals = np.linalg.norm(symmetric @ modes - modes * low_values[None, :], axis=0)
    np.save(stage / "eigenvalues.npy", values, allow_pickle=False)
    np.save(stage / "lowest_eigenvectors.npy", modes, allow_pickle=False)
    acceptance = validate_inertia(
        values,
        residuals,
        expected_negative,
        max(
            np.finfo(float).eps * np.linalg.norm(symmetric, ord=np.inf) * dimension,
            2 * antisymmetric,
        ),
    )
    acceptance.update(
        {
            "dimension": dimension,
            "antisymmetric_infinity_norm": antisymmetric,
            "lowest_eigenvalues": low_values.tolist(),
            "lowest_residuals": residuals.tolist(),
            "raw_hessian_sha256": digest(stage / "hessian_raw.npy"),
            "symmetric_hessian_sha256": digest(stage / "hessian_symmetric.npy"),
            "method": "Complete Cartesian analytic-gradient central differences, mass-weighted after differentiation; fixed exterior and full 4165-atom active region retained.",
        }
    )
    write_json(stage / "scientific_acceptance.json", acceptance)
    return acceptance
