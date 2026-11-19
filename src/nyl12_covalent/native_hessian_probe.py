"""Bounded finite-difference mass-weighted Hessian spectrum probe."""

from __future__ import annotations

import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh

from nyl12_covalent import native_sander_cli
from nyl12_covalent.block_hessian_solver import (
    cache_hvp_direction,
    estimate_block_low_modes,
    project_cached_hvp_subspace,
)
from nyl12_covalent.native_sander_cli import (
    FORCE_SOURCE_FILES,
    HARTREE_TO_KCAL_PER_MOL,
    NativeSanderEvaluator,
    load_calibration,
    prepare_native_input,
)
from scripts.prepare_dimer_continuation import read_xyz_frames

LOGGER = logging.getLogger(__name__)
FloatArray = np.ndarray


class HvpBudgetExhausted(RuntimeError):
    """Raised when a bounded Hessian probe reaches its HVP budget."""


def finite_difference_hvp(
    vector: FloatArray,
    coordinates: FloatArray,
    active_indices: np.ndarray,
    inverse_sqrt_masses: FloatArray,
    gradient: Callable[[FloatArray], FloatArray],
    step_angstrom_sqrt_amu: float,
) -> FloatArray:
    """Apply the mass-weighted Hessian by a central gradient difference.

    Args:
        vector: Flattened vector in active mass-weighted Cartesian coordinates.
        coordinates: Full-system reference coordinates in Angstrom.
        active_indices: Zero-based active atom indices in topology order.
        inverse_sqrt_masses: One inverse square-root mass per active atom.
        gradient: Callable returning full-system energy gradients in Hartree/A.
        step_angstrom_sqrt_amu: Displacement length in mass-weighted coordinates.

    Returns:
        Hessian-vector product in Hartree/(A^2 amu).
    """
    shape = (len(active_indices), 3)
    direction = np.asarray(vector, dtype=float).reshape(shape)
    norm = float(np.linalg.norm(direction))
    if not np.isfinite(direction).all() or norm <= 0.0:
        raise ValueError("Hessian direction must be finite and nonzero")
    if not np.isfinite(step_angstrom_sqrt_amu) or step_angstrom_sqrt_amu <= 0:
        raise ValueError("Finite-difference step must be finite and positive")
    displacement = (
        step_angstrom_sqrt_amu * direction / norm * inverse_sqrt_masses[:, None]
    )
    plus = np.array(coordinates, dtype=float, copy=True)
    minus = np.array(coordinates, dtype=float, copy=True)
    plus[active_indices] += displacement
    minus[active_indices] -= displacement
    grad_plus = gradient(plus)[active_indices]
    grad_minus = gradient(minus)[active_indices]
    result = (grad_plus - grad_minus) * inverse_sqrt_masses[:, None]
    result *= norm / (2.0 * step_angstrom_sqrt_amu)
    return result.reshape(-1)


def estimate_low_modes(
    matvec: Callable[[FloatArray], FloatArray],
    dimension: int,
    *,
    count: int = 2,
    tolerance: float = 1.0e-3,
    max_iterations: int = 80,
    subspace_size: int = 12,
    seed: int = 20261113,
    initial_vector: FloatArray | None = None,
) -> dict[str, Any]:
    """Estimate lowest algebraic eigenpairs with deterministic Lanczos.

    Args:
        matvec: Symmetric Hessian-vector product.
        dimension: Active Cartesian dimension.
        count: Number of lowest algebraic eigenpairs to estimate.
        tolerance: Relative ARPACK convergence tolerance.
        max_iterations: Maximum Lanczos iterations.
        subspace_size: ARPACK Krylov subspace dimension.
        seed: Seed for the deterministic starting vector.
        initial_vector: Optional physically informed start vector.

    Returns:
        Eigenvalues, eigenvectors, residual norms and operation count.
    """
    if count < 1 or dimension <= count:
        raise ValueError("Require 1 <= count < dimension")
    if max_iterations < 1 or not 0 < tolerance < 1:
        raise ValueError("Invalid eigensolver convergence controls")
    ncv = min(dimension, max(subspace_size, count + 2))
    if initial_vector is None:
        rng = np.random.default_rng(seed).normal(size=dimension)
    else:
        rng = np.asarray(initial_vector, dtype=float).reshape(-1)
        if rng.shape != (dimension,) or not np.isfinite(rng).all():
            raise ValueError(
                "Initial Lanczos vector must be finite and match dimension"
            )
        if np.linalg.norm(rng) <= 1.0e-12:
            raise ValueError("Initial Lanczos vector must be nonzero")
    rng = rng / np.linalg.norm(rng)
    operator = LinearOperator((dimension, dimension), matvec=matvec, dtype=float)
    status = "converged"
    try:
        values, vectors = eigsh(
            operator,
            k=count,
            which="SA",
            tol=tolerance,
            maxiter=max_iterations,
            ncv=ncv,
            v0=rng,
        )
    except ArpackNoConvergence as error:
        values = error.eigenvalues
        vectors = error.eigenvectors
        status = "partial_convergence"
        if values is None or vectors is None:
            values = np.array([], dtype=float)
            vectors = np.empty((dimension, 0), dtype=float)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    residuals = (
        [
            float(np.linalg.norm(matvec(vectors[:, i]) - values[i] * vectors[:, i]))
            for i in range(len(values))
        ]
        if status == "converged"
        else []
    )
    return {
        "status": status,
        "eigenvalues_hartree_per_angstrom2_amu": values.tolist(),
        "eigenvectors": vectors,
        "residual_norms_hartree_per_angstrom2_amu": residuals,
        "ncv": ncv,
    }


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_probe(
    *,
    topology_path: Path,
    restart_path: Path,
    mdin_path: Path,
    mode_path: Path,
    calibration_path: Path,
    sander_path: Path,
    output_dir: Path,
    step_sizes: tuple[float, ...] = (0.01, 0.005),
    max_iterations: int = 80,
    tolerance: float = 1.0e-3,
    max_hvp_per_step: int = 200,
    subspace_size: int = 12,
    solver: str = "arpack",
    initial_spectrum_path: Path | None = None,
    count: int = 2,
) -> Path:
    """Run the bounded low-mode probe using native analytic gradients.

    Args:
        topology_path: Exact Amber topology for the stationary-point candidate.
        restart_path: Exact 0 K optimized candidate coordinates.
        mdin_path: Original dimer input supplying the 4165-atom active mask.
        mode_path: Saved converged dimer direction in active atom order.
        calibration_path: Passing native-gradient finite-difference acceptance.
        sander_path: Native Sander executable used by the validated adapter.
        output_dir: New directory for immutable probe inputs and results.
        step_sizes: Mass-weighted central-difference displacement sizes.
        max_iterations: Maximum ARPACK iterations per step size.
        tolerance: Relative ARPACK convergence tolerance.
        max_hvp_per_step: Hard upper bound on HVPs per displacement size.
        subspace_size: ARPACK Krylov subspace dimension.
        solver: `arpack` or dimer-seeded `block_lobpcg`.
        initial_spectrum_path: Prior spectrum JSON supplying Ritz warm starts.
        count: Number of low modes requested from the selected solver.

    Returns:
        The created output directory containing manifest and spectrum JSON.
    """
    if output_dir.exists():
        raise FileExistsError(output_dir)
    if not step_sizes or any(not np.isfinite(x) or x <= 0 for x in step_sizes):
        raise ValueError("All finite-difference steps must be finite and positive")
    if max_hvp_per_step < 1:
        raise ValueError("max_hvp_per_step must be positive")
    if count < 1:
        raise ValueError("Mode count must be positive")
    if subspace_size < 4:
        raise ValueError("subspace_size must be at least four")
    if solver not in {"arpack", "block_lobpcg"}:
        raise ValueError("Unknown low-mode eigensolver")
    if initial_spectrum_path is not None and solver != "block_lobpcg":
        raise ValueError("Ritz warm starts require block_lobpcg")
    scale, calibration = load_calibration(calibration_path)
    topology = AmberParm(str(topology_path))
    restart = Rst7.open(str(restart_path))
    coordinates = np.asarray(restart.coordinates, dtype=float).reshape((-1, 3))
    if len(coordinates) != len(topology.atoms):
        raise ValueError("Restart and topology atom counts differ")
    mdin_source = mdin_path.read_text(encoding="utf-8")
    match = re.search(r"\bactive\s*=\s*'([^']+)'", mdin_source, re.IGNORECASE)
    if match is None:
        raise ValueError("Source input lacks quoted active mask")
    active = np.asarray(AmberMask(topology, match.group(1)).Selection(), dtype=bool)
    active_indices = np.flatnonzero(active)
    if len(active_indices) != 4165:
        raise ValueError(f"Expected 4165 active atoms, found {len(active_indices)}")
    masses = np.asarray([atom.mass for atom in topology.atoms], dtype=float)[active]
    if not np.isfinite(masses).all() or np.any(masses <= 0):
        raise ValueError("Active atom masses must be finite and positive")
    inverse_sqrt_masses = 1.0 / np.sqrt(masses)
    mode = read_xyz_frames(mode_path, len(active_indices))[-1]
    if not np.isfinite(mode).all() or np.linalg.norm(mode) <= 1.0e-12:
        raise ValueError("Saved dimer mode must be finite and nonzero")
    initial_vector = (np.sqrt(masses[:, None]) * mode).reshape(-1)
    initial_vector /= np.linalg.norm(initial_vector)
    initial_vectors = None
    if initial_spectrum_path is not None:
        prior = json.loads(initial_spectrum_path.read_text(encoding="utf-8"))
        prior_results = prior.get("results_by_step", [])
        if not prior_results or len(prior_results[-1].get("modes", [])) != 2:
            raise ValueError("Prior spectrum must contain two saved Ritz modes")
        initial_vectors = np.asarray(prior_results[-1]["modes"], dtype=float)
        if initial_vectors.shape != (2, len(active_indices), 3):
            raise ValueError("Prior Ritz modes do not match the active atom order")
        initial_vectors = initial_vectors.reshape((2, -1)).T
        if not np.isfinite(initial_vectors).all():
            raise ValueError("Prior Ritz modes contain nonfinite values")
    mdin = prepare_native_input(mdin_source)
    output_dir.mkdir(parents=True)
    evaluator = NativeSanderEvaluator(
        topology_path.resolve(),
        restart,
        active,
        mdin,
        sander_path.resolve(),
        output_dir / "evaluations",
        scale,
    )
    (output_dir / "evaluations").mkdir()
    source_paths = [
        topology_path,
        restart_path,
        mdin_path,
        mode_path,
        calibration_path,
        *([initial_spectrum_path] if initial_spectrum_path is not None else []),
        sander_path,
        Path(__file__),
        Path(__file__).resolve().parents[2] / "scripts/run_native_hessian_probe.py",
        Path(native_sander_cli.__file__),
        Path(__file__).with_name("block_hessian_solver.py"),
        *FORCE_SOURCE_FILES,
    ]
    manifest = {
        "purpose": "lowest mass-weighted Hessian eigenvalues by FD HVP",
        "active_atom_count": len(active_indices),
        "dimension": 3 * len(active_indices),
        "fixed_exterior": True,
        "temperature_kelvin": 0,
        "scipy_version": scipy.__version__,
        "force_calibration": calibration,
        "finite_difference_steps_angstrom_sqrt_amu": list(step_sizes),
        "max_hvp_per_step": max_hvp_per_step,
        "maximum_native_gradient_evaluations": (2 * len(step_sizes) * max_hvp_per_step),
        "eigensolver": {
            "method": solver,
            "which": "lowest algebraic",
            "count": count,
            "tolerance": tolerance,
            "max_iterations": max_iterations,
            "seed": 20261113,
            "initial_vector": "final saved dimer mode multiplied by sqrt(active masses)",
            "warm_start_spectrum": (
                str(initial_spectrum_path.resolve())
                if initial_spectrum_path is not None
                else None
            ),
        },
        "source_sha256": {str(p.resolve()): sha256(p) for p in source_paths},
        "prepared_mdin_sha256": hashlib.sha256(mdin.encode()).hexdigest(),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    results = []
    for step in step_sizes:

        def gradient(coords: FloatArray) -> FloatArray:
            _, grad_kcal = evaluator(coords)
            return grad_kcal / HARTREE_TO_KCAL_PER_MOL

        hvp_count = 0
        cached_basis: list[FloatArray] = []
        cached_products: list[FloatArray] = []
        cached_ritz_history: list[dict[str, Any]] = []

        def matvec(vector: FloatArray) -> FloatArray:
            nonlocal hvp_count
            if hvp_count >= max_hvp_per_step:
                raise HvpBudgetExhausted(
                    f"HVP budget exhausted at {max_hvp_per_step} for step {step}"
                )
            hvp_count += 1
            product = finite_difference_hvp(
                vector,
                coordinates,
                active_indices,
                inverse_sqrt_masses,
                gradient,
                step,
            )
            cache_hvp_direction(cached_basis, cached_products, vector, product)
            checkpoint_interval = max(1, 5 * count)
            if (
                solver == "block_lobpcg"
                and hvp_count % checkpoint_interval == 0
                and len(cached_basis) >= count
            ):
                checkpoint = project_cached_hvp_subspace(
                    cached_basis,
                    cached_products,
                    count=count,
                    tolerance=tolerance,
                )
                checkpoint.pop("eigenvectors")
                checkpoint["hvp_count"] = hvp_count
                cached_ritz_history.append(checkpoint)
            return product

        try:
            if solver == "block_lobpcg":
                spectrum = estimate_block_low_modes(
                    matvec,
                    3 * len(active_indices),
                    initial_vector,
                    count=count,
                    initial_vectors=initial_vectors,
                    tolerance=tolerance,
                    max_iterations=max_iterations,
                )
            else:
                spectrum = estimate_low_modes(
                    matvec,
                    3 * len(active_indices),
                    count=count,
                    tolerance=tolerance,
                    max_iterations=max_iterations,
                    seed=20261113,
                    initial_vector=initial_vector,
                    subspace_size=subspace_size,
                )
        except HvpBudgetExhausted as error:
            partial = (
                project_cached_hvp_subspace(
                    cached_basis,
                    cached_products,
                    count=count,
                    tolerance=tolerance,
                )
                if solver == "block_lobpcg" and len(cached_basis) >= count
                else {"eigenpairs": []}
            )
            partial["status"] = "hvp_budget_exhausted"
            partial["step_angstrom_sqrt_amu"] = step
            partial["hvp_count_including_residual_checks"] = hvp_count
            partial["gradient_evaluations"] = len(evaluator.evidence)
            partial["diagnostic"] = str(error)
            partial["cached_ritz_residual_history"] = cached_ritz_history
            if "eigenvectors" in partial:
                vectors = partial.pop("eigenvectors")
                partial["modes"] = [
                    vectors[:, i].reshape((-1, 3)).tolist()
                    for i in range(vectors.shape[1])
                ]
                partial["eigenvalues_kcal_per_mol_angstrom2_amu"] = (
                    np.asarray(partial["eigenvalues_hartree_per_angstrom2_amu"])
                    * HARTREE_TO_KCAL_PER_MOL
                ).tolist()
            results.append(partial)
            break
        spectrum["step_angstrom_sqrt_amu"] = step
        vectors = spectrum.pop("eigenvectors")
        values = np.asarray(spectrum["eigenvalues_hartree_per_angstrom2_amu"])
        spectrum["modes"] = [
            vectors[:, i].reshape((-1, 3)).tolist() for i in range(vectors.shape[1])
        ]
        spectrum["eigenvalues_kcal_per_mol_angstrom2_amu"] = (
            values * HARTREE_TO_KCAL_PER_MOL
        ).tolist()
        spectrum["gradient_evaluations"] = len(evaluator.evidence)
        spectrum["hvp_count_including_residual_checks"] = hvp_count
        results.append(spectrum)
        LOGGER.info("Completed HVP spectrum at step %.6g", step)
    payload = {
        "results_by_step": results,
        "status": (
            "converged"
            if results and all(item.get("status") == "converged" for item in results)
            else "incomplete_or_partial"
        ),
        "step_sensitivity_hartree_per_angstrom2_amu": (
            np.asarray(results[0]["eigenvalues_hartree_per_angstrom2_amu"])
            - np.asarray(results[-1]["eigenvalues_hartree_per_angstrom2_amu"])
        ).tolist()
        if len(results) > 1
        and all("eigenvalues_hartree_per_angstrom2_amu" in item for item in results)
        else None,
        "gradient_evaluations": evaluator.evidence,
        "interpretation": "Finite-difference low-mode estimate; assess residuals and step sensitivity before TS acceptance.",
    }
    (output_dir / "spectrum.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return output_dir
