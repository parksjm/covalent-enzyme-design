"""Analytic occupancy/distance biases and active-coordinate optimizer adapter."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from .proton_coordination import FloatArray, coordination_transfer

HARTREE_TO_KCAL_MOL = 627.5094740631
PhysicalEvaluator = Callable[[FloatArray], tuple[float, FloatArray]]


@dataclass(frozen=True)
class OccupancyBiasConfig:
    """Zero-based atom mapping and harmonic targets (Å and kcal/mol).

    The occupancy coordinate is CN(acceptor,H) minus CN(donor,H).
    Springs have units kcal/mol per squared coordinate unit.
    """

    donor: int
    acceptor: int
    hydrogens: Sequence[int]
    attack_atoms: tuple[int, int]
    q_target: float
    r_target: float
    k_q: float
    k_r: float
    radius: float = 1.25


@dataclass(frozen=True)
class BiasAssessment:
    """Occupancy, distance, individual energies and full Cartesian derivative."""

    cn_donor: float
    cn_acceptor: float
    q: float
    r: float
    energy_q: float
    energy_r: float
    energy: float
    gradient: FloatArray


@dataclass(frozen=True)
class BiasedEvaluation:
    """Physical and biased evaluation with all energies in kcal/mol."""

    physical_energy: float
    energy: float
    gradient: FloatArray
    bias: BiasAssessment


def evaluate_bias(
    coordinates: FloatArray, config: OccupancyBiasConfig
) -> BiasAssessment:
    """Compute two harmonic biases and exact analytic Cartesian derivatives.

    Args:
        coordinates: Finite full-system Cartesian array (N,3), in Å.
        config: Explicit atom mapping, targets and nonnegative springs.

    Returns:
        Smooth occupancies, transfer coordinate, distance and bias derivatives.

    Raises:
        ValueError: If mapping, targets, springs or coordinates are invalid.
    """
    values = [config.q_target, config.r_target, config.k_q, config.k_r]
    if not np.isfinite(values).all() or min(config.k_q, config.k_r) < 0:
        raise ValueError("Targets and springs must be finite; springs nonnegative")
    if config.r_target <= 0:
        raise ValueError("Distance target must be positive")
    xyz = np.asarray(coordinates, dtype=float)
    transfer = coordination_transfer(
        xyz, config.donor, config.acceptor, config.hydrogens, config.radius
    )
    raw_pair = np.asarray(config.attack_atoms)
    if raw_pair.shape != (2,) or not np.issubdtype(raw_pair.dtype, np.integer):
        raise ValueError("Attack atoms must be two integer indices")
    first, second = raw_pair
    if first == second or min(first, second) < 0 or max(first, second) >= len(xyz):
        raise ValueError("Attack atoms must be distinct valid indices")
    delta = xyz[first] - xyz[second]
    distance = float(np.linalg.norm(delta))
    if distance == 0:
        raise ValueError("Attack distance must be nonzero for its derivative")
    dq, dr = transfer.value - config.q_target, distance - config.r_target
    energy_q, energy_r = 0.5 * config.k_q * dq**2, 0.5 * config.k_r * dr**2
    gradient = config.k_q * dq * transfer.gradient
    distance_gradient = config.k_r * dr * delta / distance
    gradient[first] += distance_gradient
    gradient[second] -= distance_gradient
    return BiasAssessment(
        transfer.donor,
        transfer.acceptor,
        transfer.value,
        distance,
        energy_q,
        energy_r,
        energy_q + energy_r,
        gradient,
    )


def evaluate_biased(
    coordinates: FloatArray,
    physical_evaluator: PhysicalEvaluator,
    config: OccupancyBiasConfig,
) -> BiasedEvaluation:
    """Combine Hartree physical energy and kcal/mol/Å physical gradient.

    Args:
        coordinates: Full-system Cartesian array (N,3), in Å.
        physical_evaluator: Callable returning Hartree energy and full gradient
            already in kcal/mol/Å, as supplied by the native Sander evaluator.
        config: Explicit harmonic bias definition.

    Returns:
        Consistent total kcal/mol energy and kcal/mol/Å gradient.

    Raises:
        ValueError: If physical energy or gradient is nonfinite or malformed.
    """
    xyz = np.array(coordinates, dtype=float, copy=True)
    bias = evaluate_bias(xyz, config)
    energy_hartree, physical_gradient = physical_evaluator(xyz)
    physical_gradient = np.asarray(physical_gradient, dtype=float)
    if not np.isscalar(energy_hartree) or not np.isfinite(energy_hartree):
        raise ValueError("Physical energy must be a finite scalar in Hartree")
    if physical_gradient.shape != xyz.shape or not np.isfinite(physical_gradient).all():
        raise ValueError("Physical gradient must be a finite full-system (N,3) array")
    physical_energy = float(energy_hartree) * HARTREE_TO_KCAL_MOL
    return BiasedEvaluation(
        physical_energy,
        physical_energy + bias.energy,
        physical_gradient + bias.gradient,
        bias,
    )


class ActiveBiasCallback:
    """Optimizer callback that freezes outer atoms by active-coordinate selection."""

    def __init__(
        self,
        physical_evaluator: PhysicalEvaluator,
        reference_coordinates: FloatArray,
        active_indices: Sequence[int],
        config: OccupancyBiasConfig,
    ) -> None:
        """Store immutable reference positions and explicit active atom indices.

        Args:
            physical_evaluator: Hartree-energy, kcal/mol/Å-gradient callable.
            reference_coordinates: Initial full-system Cartesian (N,3) array.
            active_indices: Distinct zero-based movable atom indices.
            config: Explicit bias mapping, targets and springs.

        Raises:
            ValueError: If reference positions or active selection are invalid.
        """
        self.reference_coordinates = np.array(
            reference_coordinates, dtype=float, copy=True
        )
        evaluate_bias(self.reference_coordinates, config)
        active = np.asarray(active_indices)
        if (
            active.ndim != 1
            or not len(active)
            or not np.issubdtype(active.dtype, np.integer)
        ):
            raise ValueError("Active indices must be a nonempty integer vector")
        if (
            len(set(active)) != len(active)
            or active.min() < 0
            or active.max() >= len(self.reference_coordinates)
        ):
            raise ValueError("Active indices must be distinct valid atom indices")
        self.active_indices = active.astype(int, copy=True)
        self.reference_coordinates.setflags(write=False)
        self.active_indices.setflags(write=False)
        self.physical_evaluator = physical_evaluator
        self.config = config
        self.last_evaluation: BiasedEvaluation | None = None

    def coordinates(self, active_coordinates: FloatArray) -> FloatArray:
        """Reconstruct full positions while retaining exact frozen coordinates.

        Args:
            active_coordinates: Flat array of 3 times the active atom count.

        Returns:
            Full-system Cartesian coordinates in Å.

        Raises:
            ValueError: If active positions have wrong shape or nonfinite values.
        """
        active = np.asarray(active_coordinates, dtype=float)
        if (
            active.shape != (3 * len(self.active_indices),)
            or not np.isfinite(active).all()
        ):
            raise ValueError("Active positions must be a finite flat 3*n_active vector")
        xyz = self.reference_coordinates.copy()
        xyz[self.active_indices] = active.reshape(-1, 3)
        return xyz

    def __call__(self, active_coordinates: FloatArray) -> tuple[float, FloatArray]:
        """Return total kcal/mol energy and selected kcal/mol/Å gradient.

        Args:
            active_coordinates: Flat active Cartesian positions in Å.

        Returns:
            Energy and matching flat gradient for an optimizer using jac=True.
        """
        self.last_evaluation = None
        evaluation = evaluate_biased(
            self.coordinates(active_coordinates), self.physical_evaluator, self.config
        )
        self.last_evaluation = evaluation
        return evaluation.energy, evaluation.gradient[
            self.active_indices
        ].ravel().copy()
