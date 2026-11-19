"""Composable analytic reaction-coordinate biases; no proton identity restraints."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, replace

import numpy as np

from .occupancy_bias import (
    HARTREE_TO_KCAL_MOL,
    ActiveBiasCallback,
    PhysicalEvaluator,
)
from .proton_coordination import (
    FloatArray,
    coordination_transfer,
    smooth_coordination,
)


@dataclass(frozen=True)
class DistanceCV:
    """Distance in Å between distinct zero-based atom indices."""

    first: int
    second: int


@dataclass(frozen=True)
class CoordinationCV:
    """Smooth hydrogen coordination at one center, with switching radius in Å."""

    center: int
    hydrogens: Sequence[int]
    radius: float = 1.25


@dataclass(frozen=True)
class CoordinationDifferenceCV:
    """CN(acceptor,H) minus CN(donor,H), using all declared real QM hydrogens."""

    donor: int
    acceptor: int
    hydrogens: Sequence[int]
    radius: float = 1.25


PrimitiveCV = DistanceCV | CoordinationCV | CoordinationDifferenceCV


@dataclass(frozen=True)
class WeightedCV:
    """Primitive CV and finite multiplier converting its units to output units."""

    cv: PrimitiveCV
    weight: float

    def __post_init__(self) -> None:
        if not isinstance(
            self.cv, (DistanceCV, CoordinationCV, CoordinationDifferenceCV)
        ):
            raise ValueError(
                "Weighted terms require primitive CVs; nesting is forbidden"
            )
        if not np.isscalar(self.weight) or not np.isfinite(self.weight):
            raise ValueError("CV weights must be finite scalars")
        if isinstance(self.cv, (CoordinationCV, CoordinationDifferenceCV)):
            object.__setattr__(
                self, "cv", replace(self.cv, hydrogens=tuple(self.cv.hydrogens))
            )


@dataclass(frozen=True)
class LinearCombinationCV:
    """Weighted sum with explicit output units and no nested combinations.

    Each weight converts its primitive coordinate unit into ``output_unit``.
    For an Å-valued sum, distance weights are dimensionless and coordination
    weights have units Å. No unit conversions are inferred or applied silently.
    """

    terms: tuple[WeightedCV, ...]
    output_unit: str

    def __post_init__(self) -> None:
        terms = tuple(self.terms)
        if not terms or not all(isinstance(term, WeightedCV) for term in terms):
            raise ValueError("Linear combination requires nonempty WeightedCV terms")
        if not isinstance(self.output_unit, str) or not self.output_unit.strip():
            raise ValueError("Linear combination requires explicit output units")
        object.__setattr__(self, "terms", terms)


CollectiveVariable = PrimitiveCV | LinearCombinationCV


@dataclass(frozen=True)
class HarmonicTarget:
    """One 0.5*spring*(CV-target)^2 term, with energy in kcal/mol.

    Spring units are kcal/mol per squared CV unit. Coordination targets are
    dimensionless. A distance target must be positive; springs may be zero.
    """

    cv: CollectiveVariable
    target: float
    spring: float

    def to_metadata(self) -> dict:
        """Return JSON-compatible metadata including the explicit CV type.

        Returns:
            Atom mapping, switching radius, target and spring in plain JSON types.
        """

        def metadata(cv: CollectiveVariable) -> dict:
            if isinstance(cv, LinearCombinationCV):
                return {
                    "type": type(cv).__name__,
                    "output_unit": cv.output_unit,
                    "terms": [
                        {"cv": metadata(term.cv), "weight": float(term.weight)}
                        for term in cv.terms
                    ],
                }
            values = asdict(cv)
            for key, value in values.items():
                if isinstance(value, (tuple, list, np.ndarray)):
                    values[key] = [int(index) for index in value]
                elif isinstance(value, np.generic):
                    values[key] = value.item()
            return {"type": type(cv).__name__, **values}

        return {
            "cv": metadata(self.cv),
            "target": float(self.target),
            "spring": float(self.spring),
        }


@dataclass(frozen=True)
class ReactionBiasAssessment:
    """Ordered CV values, targets, term energies and full analytic derivative."""

    values: tuple[float, ...]
    targets: tuple[float, ...]
    energies: tuple[float, ...]
    energy: float
    gradient: FloatArray


@dataclass(frozen=True)
class ReactionBiasedEvaluation:
    """Separate physical/total kcal/mol energies and full kcal/mol/Å gradient."""

    physical_energy: float
    energy: float
    gradient: FloatArray
    bias: ReactionBiasAssessment


def _coordinate(xyz: FloatArray, cv: CollectiveVariable) -> tuple[float, FloatArray]:
    if isinstance(cv, LinearCombinationCV):
        value, gradient = 0.0, np.zeros_like(xyz)
        for term in cv.terms:
            component, derivative = _coordinate(xyz, term.cv)
            value += term.weight * component
            gradient += term.weight * derivative
        return float(value), gradient
    if isinstance(cv, CoordinationCV):
        return smooth_coordination(xyz, cv.center, cv.hydrogens, cv.radius)
    if isinstance(cv, CoordinationDifferenceCV):
        result = coordination_transfer(
            xyz, cv.donor, cv.acceptor, cv.hydrogens, cv.radius
        )
        return result.value, result.gradient
    if not isinstance(cv, DistanceCV):
        raise ValueError("Unsupported collective variable type")
    pair = np.asarray((cv.first, cv.second))
    if (
        not np.issubdtype(pair.dtype, np.integer)
        or pair.min() < 0
        or pair.max() >= len(xyz)
        or cv.first == cv.second
    ):
        raise ValueError("Distance atoms must be distinct valid integer indices")
    delta = xyz[cv.first] - xyz[cv.second]
    distance = float(np.linalg.norm(delta))
    if distance == 0:
        raise ValueError("Distance must be nonzero for its derivative")
    gradient = np.zeros_like(xyz)
    gradient[cv.first] = delta / distance
    gradient[cv.second] = -delta / distance
    return distance, gradient


def evaluate_bias(
    coordinates: FloatArray, targets: Sequence[HarmonicTarget]
) -> ReactionBiasAssessment:
    """Evaluate independent harmonic targets and sum exact Cartesian gradients.

    Args:
        coordinates: Finite full-system Cartesian (N,3) array in Å.
        targets: Nonempty ordered terms; each carries its own CV and target.

    Returns:
        Ordered coordinate/target/energy tuples and summed energy/gradient.

    Raises:
        ValueError: If positions, atom mappings, targets or springs are invalid.
    """
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError("Coordinates must be a finite (N,3) array")
    if not len(targets):
        raise ValueError("At least one harmonic target is required")
    values, target_values, energies = [], [], []
    gradient = np.zeros_like(xyz)
    for term in targets:
        if not isinstance(term, HarmonicTarget):
            raise ValueError("Each term must be a HarmonicTarget")
        if (
            not np.isscalar(term.target)
            or not np.isscalar(term.spring)
            or not np.isfinite([term.target, term.spring]).all()
            or term.spring < 0
        ):
            raise ValueError("Targets/springs must be finite; springs nonnegative")
        if isinstance(term.cv, DistanceCV) and term.target <= 0:
            raise ValueError("Distance targets must be positive")
        value, derivative = _coordinate(xyz, term.cv)
        delta = value - term.target
        values.append(value)
        target_values.append(float(term.target))
        energies.append(float(0.5 * term.spring * delta**2))
        gradient += term.spring * delta * derivative
    return ReactionBiasAssessment(
        tuple(values), tuple(target_values), tuple(energies), sum(energies), gradient
    )


def evaluate_biased(
    coordinates: FloatArray,
    physical_evaluator: PhysicalEvaluator,
    targets: Sequence[HarmonicTarget],
) -> ReactionBiasedEvaluation:
    """Combine native physical energy (Hartree) and gradients (kcal/mol/Å).

    Args:
        coordinates: Full-system Cartesian (N,3) array in Å.
        physical_evaluator: Callable returning Hartree energy and full gradient.
        targets: Explicit harmonic terms with independent targets.

    Returns:
        Physical and total kcal/mol energies, total gradient and bias assessment.

    Raises:
        ValueError: If bias inputs or physical outputs are invalid/nonfinite.
    """
    xyz = np.array(coordinates, dtype=float, copy=True)
    bias = evaluate_bias(xyz, targets)
    energy, gradient = physical_evaluator(xyz)
    gradient = np.asarray(gradient, dtype=float)
    if not np.isscalar(energy) or not np.isfinite(energy):
        raise ValueError("Physical energy must be a finite scalar in Hartree")
    if gradient.shape != xyz.shape or not np.isfinite(gradient).all():
        raise ValueError("Physical gradient must be a finite full-system (N,3) array")
    physical = float(energy) * HARTREE_TO_KCAL_MOL
    return ReactionBiasedEvaluation(
        physical, physical + bias.energy, gradient + bias.gradient, bias
    )


class ReactionBiasCallback:
    """Active-coordinate optimizer adapter with exact frozen reference positions."""

    # Reuse the validated reconstruction/shape contract without occupancy config.
    coordinates = ActiveBiasCallback.coordinates

    def __init__(
        self,
        physical_evaluator: PhysicalEvaluator,
        reference_coordinates: FloatArray,
        active_indices: Sequence[int],
        targets: Sequence[HarmonicTarget],
    ) -> None:
        """Validate and retain explicit terms and immutable atom selections.

        Args:
            physical_evaluator: Hartree-energy, kcal/mol/Å-gradient callable.
            reference_coordinates: Full-system Cartesian (N,3) reference in Å.
            active_indices: Nonempty distinct movable zero-based atom indices.
            targets: Nonempty independent harmonic terms.

        Raises:
            ValueError: If reference coordinates, terms or selection are invalid.
        """
        self.reference_coordinates = np.array(reference_coordinates, dtype=float)
        # Snapshot mutable hydrogen lists as well as the sequence of terms.
        frozen = []
        for term in targets:
            if isinstance(term, HarmonicTarget) and isinstance(
                term.cv, (CoordinationCV, CoordinationDifferenceCV)
            ):
                term = replace(
                    term, cv=replace(term.cv, hydrogens=tuple(term.cv.hydrogens))
                )
            frozen.append(term)
        self.targets = tuple(frozen)
        evaluate_bias(self.reference_coordinates, self.targets)
        active = np.asarray(active_indices)
        if (
            active.ndim != 1
            or not len(active)
            or not np.issubdtype(active.dtype, np.integer)
            or len(set(active)) != len(active)
            or active.min() < 0
            or active.max() >= len(self.reference_coordinates)
        ):
            raise ValueError("Active indices must be nonempty distinct valid integers")
        self.active_indices = active.astype(int, copy=True)
        self.active_indices.setflags(write=False)
        self.reference_coordinates.setflags(write=False)
        self.physical_evaluator = physical_evaluator
        self.last_evaluation: ReactionBiasedEvaluation | None = None

    def __call__(self, active_coordinates: FloatArray) -> tuple[float, FloatArray]:
        """Return total kcal/mol energy and corresponding flat active gradient.

        Args:
            active_coordinates: Finite flat active Cartesian positions in Å.

        Returns:
            Total energy and kcal/mol/Å gradient for a jac=True optimizer.
        """
        self.last_evaluation = None
        result = evaluate_biased(
            self.coordinates(active_coordinates), self.physical_evaluator, self.targets
        )
        self.last_evaluation = result
        return result.energy, result.gradient[self.active_indices].ravel().copy()
