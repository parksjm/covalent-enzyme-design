"""Independent D2 geometry checks for ordered mature repeats AE/BF/CG/DH.

Coordinates are row vectors. An operator acts as ``(x-center) @ R.T+center``.
Fitting estimates a frame, never qualifies a source or alters an input artifact.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
AtomKey = tuple[str, int, str]
COUPLED_REPEATS = (("A", "E"), ("B", "F"), ("C", "G"), ("D", "H"))


def _array(value: object, shape: tuple[int, ...], name: str) -> FloatArray:
    result = np.asarray(value, dtype=float)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError(f"{name} must have finite values and shape {shape}")
    return result


def _coordinates(value: object) -> FloatArray:
    result = np.asarray(value, dtype=float)
    if result.ndim != 3 or result.shape[0] != 4 or result.shape[2] != 3:
        raise ValueError("Coordinates must have shape (4, n_atoms, 3)")
    if result.shape[1] < 4 or not np.isfinite(result).all():
        raise ValueError("At least four finite corresponding atoms are required")
    return result


def validate_correspondence(
    atom_keys: Sequence[Sequence[AtomKey]], n_atoms: int
) -> None:
    """Require identical, unique alpha/beta atom keys in all four repeats.

    Args:
        atom_keys: Four ordered lists of (alpha/beta, canonical residue, atom).
        n_atoms: Number of atoms per coupled repeat.

    Raises:
        ValueError: Mapping is incomplete, ambiguous, or inconsistent.
    """
    keys = [tuple(repeat) for repeat in atom_keys]
    if len(keys) != 4 or any(len(repeat) != n_atoms for repeat in keys):
        raise ValueError("Four complete repeat mappings are required")
    for repeat in keys:
        if any(
            len(key) != 3
            or key[0] not in {"alpha", "beta"}
            or not isinstance(key[1], int)
            or key[1] < 1
            or not isinstance(key[2], str)
            or not key[2].strip()
            for key in repeat
        ):
            raise ValueError("Invalid segment, residue or atom identity")
        if len(set(repeat)) != n_atoms:
            raise ValueError("Duplicate atom correspondence")
        if {key[0] for key in repeat} != {"alpha", "beta"}:
            raise ValueError("Each repeat must include both alpha and beta")
    if any(repeat != keys[0] for repeat in keys[1:]):
        raise ValueError("Correspondence order/identities differ between repeats")


@dataclass(frozen=True)
class D2Operators:
    """Exact D2 rotations in AE/BF/CG/DH order about a shared center."""

    rotations: FloatArray
    center: FloatArray

    def __post_init__(self) -> None:
        rotations = _array(self.rotations, (4, 3, 3), "rotations").copy()
        center = _array(self.center, (3,), "center").copy()
        if not np.allclose(rotations[0], np.eye(3), atol=1e-10, rtol=0):
            raise ValueError("First operator must be identity")
        for index, rotation in enumerate(rotations):
            if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10, rtol=0):
                raise ValueError("Operators must be orthogonal")
            if not np.isclose(np.linalg.det(rotation), 1, atol=1e-10, rtol=0):
                raise ValueError("Operators must be proper rotations")
            if index and not np.isclose(np.trace(rotation), -1, atol=1e-10, rtol=0):
                raise ValueError("D2 nonidentity operators must be half turns")
        for first in range(4):
            for second in range(4):
                if not np.allclose(
                    rotations[first] @ rotations[second],
                    rotations[first ^ second],
                    atol=1e-10,
                    rtol=0,
                ):
                    raise ValueError("D2 operator multiplication/closure failed")
        rotations.setflags(write=False)
        center.setflags(write=False)
        object.__setattr__(self, "rotations", rotations)
        object.__setattr__(self, "center", center)


def fit_d2_operators(
    coordinates: FloatArray, atom_keys: Sequence[Sequence[AtomKey]]
) -> D2Operators:
    """Estimate exact D2 operators from joint alpha/beta pair-fit axes.

    Args:
        coordinates: Matching atoms, shape (4, n_atoms, 3), in AE/BF/CG/DH order.
        atom_keys: Explicit ordered atom identity correspondence for each repeat.

    Returns:
        Exact operators from polar-orthogonalized axes and least-squares center.

    Raises:
        ValueError: Inputs are nonfinite, mapping is invalid or geometry degenerate.
    """
    xyz = _coordinates(coordinates)
    validate_correspondence(atom_keys, xyz.shape[1])
    centered = xyz - xyz.mean(axis=1, keepdims=True)
    if any(np.linalg.matrix_rank(repeat) < 3 for repeat in centered):
        raise ValueError("Repeat coordinates must span three dimensions")
    axes = []
    for target in centered[1:]:
        left, _, right = np.linalg.svd(centered[0].T @ target)
        correction = np.diag([1.0, 1.0, np.linalg.det(left @ right)])
        rotation = (left @ correction @ right).T
        # The eigenvector of the symmetric part with largest eigenvalue is
        # the rotation axis, including rotations slightly away from pi.
        _, vectors = np.linalg.eigh((rotation + rotation.T) / 2)
        axes.append(vectors[:, -1])
    left, singular, right = np.linalg.svd(np.stack(axes, axis=1))
    if singular[-1] < 1e-6:
        raise ValueError("Fitted axes do not define a D2 frame")
    frame = left @ right
    rotations = np.stack(
        [np.eye(3)] + [2 * np.outer(axis, axis) - np.eye(3) for axis in frame.T]
    )
    matrices = np.concatenate([np.eye(3) - rot for rot in rotations[1:]])
    centroids = xyz.mean(axis=1)
    offsets = np.concatenate(
        [centroids[i] - rot @ centroids[0] for i, rot in enumerate(rotations[1:], 1)]
    )
    center, _, _, _ = np.linalg.lstsq(matrices, offsets, rcond=None)
    return D2Operators(rotations, center)


def reference_projection(coordinates: FloatArray, operators: D2Operators) -> FloatArray:
    """Generate exact copies from repeat AE without moving that reference repeat.

    Args:
        coordinates: Source coordinates of shape (4, n_atoms, 3).
        operators: Qualified frame operators in the source coordinate frame.

    Returns:
        New array; no source coordinates are mutated.
    """
    xyz = _coordinates(coordinates)
    return np.stack(
        [
            (xyz[0] - operators.center) @ rot.T + operators.center
            for rot in operators.rotations
        ]
    )


def operator_closure_metrics(operators: D2Operators) -> dict[str, float]:
    """Return maximum dimensionless rotation and Angstrom translation errors."""
    rotations, center = operators.rotations, operators.center
    translations = center - np.einsum("kij,j->ki", rotations, center)
    rotation_errors, translation_errors = [], []
    for first in range(4):
        for second in range(4):
            target = first ^ second
            rotation_errors.append(
                np.linalg.norm(rotations[first] @ rotations[second] - rotations[target])
            )
            translation_errors.append(
                np.linalg.norm(
                    rotations[first] @ translations[second]
                    + translations[first]
                    - translations[target]
                )
            )
    return {
        "rotation_max_error": float(max(rotation_errors)),
        "translation_max_error_A": float(max(translation_errors)),
    }


def displacement_metrics(source: FloatArray, output: FloatArray) -> dict[str, float]:
    """Measure absolute atom displacement, with no alignment or atom remapping.

    Args:
        source: Nonempty finite atom coordinates with last dimension three.
        output: Coordinates with identical shape and atom ordering.

    Returns:
        RMS and maximum Euclidean displacement in Angstrom.
    """
    before = np.asarray(source, dtype=float)
    if before.ndim < 2 or before.shape[-1] != 3 or before.size == 0:
        raise ValueError("Source must contain atom coordinates")
    before = _array(before, before.shape, "source")
    after = _array(output, before.shape, "output")
    distances = np.linalg.norm(after - before, axis=-1)
    return {
        "rmsd_A": float(np.sqrt(np.mean(distances**2))),
        "max_deviation_A": float(distances.max()),
    }


def measure_d2_residual(
    coordinates: FloatArray,
    operators: D2Operators,
    *,
    rmsd_limit_A: float,
    max_deviation_limit_A: float,
) -> dict[str, object]:
    """Audit a fixed frame against shared consensus; never fit chains separately.

    Args:
        coordinates: Matching backbone atoms in AE/BF/CG/DH order.
        operators: Predeclared source-frame operators, never refit to output.
        rmsd_limit_A: Positive predeclared RMS displacement tolerance.
        max_deviation_limit_A: Positive predeclared maximum displacement tolerance.

    Returns:
        Residuals to the least-squares shared repeat and explicit tolerance result.
        This result alone is not a full backbone symmetry acceptance receipt.
    """
    xyz = _coordinates(coordinates)
    limits = np.array([rmsd_limit_A, max_deviation_limit_A], dtype=float)
    if not np.isfinite(limits).all() or np.any(limits <= 0):
        raise ValueError("Residual tolerances must be finite and positive")
    local = np.stack(
        [
            (repeat - operators.center) @ rot
            for repeat, rot in zip(xyz, operators.rotations, strict=True)
        ]
    )
    consensus = local.mean(axis=0)
    projected = np.stack(
        [consensus @ rot.T + operators.center for rot in operators.rotations]
    )
    metrics = displacement_metrics(xyz, projected)
    per_repeat = [
        displacement_metrics(a, b) for a, b in zip(xyz, projected, strict=True)
    ]
    return {
        "backbone_rmsd_A": metrics["rmsd_A"],
        "backbone_max_deviation_A": metrics["max_deviation_A"],
        "backbone_rmsd_A_limit": float(rmsd_limit_A),
        "backbone_max_deviation_A_limit": float(max_deviation_limit_A),
        "within_tolerances": bool(
            metrics["rmsd_A"] <= rmsd_limit_A
            and metrics["max_deviation_A"] <= max_deviation_limit_A
        ),
        "per_repeat": per_repeat,
    }


def measure_tied_displacement_residual(
    source: FloatArray,
    output: FloatArray,
    operators: D2Operators,
    *,
    rmsd_limit_A: float,
    max_deviation_limit_A: float,
) -> dict[str, object]:
    """Audit rotation-tied output-minus-source motion without source projection.

    Args:
        source: Immutable matching atoms in AE/BF/CG/DH order, shape (4,n,3).
        output: Final matching coordinates in the identical frame and ordering.
        operators: Predeclared rotations; translation cancels for displacement.
        rmsd_limit_A: Positive predeclared RMS displacement residual tolerance.
        max_deviation_limit_A: Positive predeclared maximum residual tolerance.

    Returns:
        Tied-motion residuals after backrotation and common averaging. Absolute
        coordinate D2 symmetry is a separate diagnostic, not established here.
    """
    before = _coordinates(source)
    after = _array(output, before.shape, "output")
    result = measure_d2_residual(
        after - before,
        D2Operators(operators.rotations, np.zeros(3)),
        rmsd_limit_A=rmsd_limit_A,
        max_deviation_limit_A=max_deviation_limit_A,
    )
    return {
        key.replace("backbone_", "tied_displacement_"): value
        for key, value in result.items()
    }
