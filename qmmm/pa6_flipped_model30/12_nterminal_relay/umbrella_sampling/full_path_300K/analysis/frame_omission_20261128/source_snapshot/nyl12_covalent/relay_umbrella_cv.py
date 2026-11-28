"""Distance-difference CVs and native Amber bias helpers for relay umbrellas.

The functions here are adapters for the existing QM/MM umbrella workflow.
They retain the historical ``q_acyl``/``q_h`` default while allowing the
N-terminal relay charts to choose any two named features from the same set.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from collections.abc import Iterable, Iterator, Mapping, Sequence

import numpy as np

from nyl12_covalent.neb_umbrella_observables import _distances


@dataclass(frozen=True)
class DistanceDifference:
    """A difference between two interatomic distances, using 1-based IDs."""

    name: str
    atom_ids: tuple[int, int, int, int]

    def __post_init__(self) -> None:
        if not self.name or len(self.atom_ids) != 4:
            raise ValueError("A named feature needs exactly four atom IDs")
        if any(type(atom_id) is not int or atom_id < 1 for atom_id in self.atom_ids):
            raise ValueError("Feature atom IDs must be positive one-based integers")

    def value(self, xyz: np.ndarray, box: Sequence[float] | None = None) -> float:
        """Calculate d(i,j)-d(k,l) in Angstrom from one full coordinate frame."""
        coordinates = np.asarray(xyz, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("Coordinates must have shape (atoms, 3)")
        if not np.isfinite(coordinates).all():
            raise ValueError("Coordinates must be finite")
        if max(self.atom_ids) > len(coordinates):
            raise ValueError("Feature atom ID exceeds coordinate array")
        i, j, k, fourth = (atom_id - 1 for atom_id in self.atom_ids)
        if box is not None:
            cell = np.asarray(box, dtype=float)
            if cell.shape != (6,) or not np.isfinite(cell).all():
                raise ValueError("Periodic cell must contain six finite Amber values")
            first = float(_distances(coordinates, [i + 1], [j + 1], cell)[0, 0])
            second = float(_distances(coordinates, [k + 1], [fourth + 1], cell)[0, 0])
            return first - second
        return float(
            np.linalg.norm(coordinates[i] - coordinates[j])
            - np.linalg.norm(coordinates[k] - coordinates[fourth])
        )

    def amber_restraint(
        self, center_A: float, force_constant: float, outer_extent_A: float = 10.0
    ) -> str:
        """Render one weighted-distance Amber ``&rst`` record.

        Amber's central energy is ``rk*(xi-center)**2`` in kcal/mol, without
        a factor of one half. Amber applies tangent-linear tails outside the
        ``center +/- outer_extent`` bounds.
        """
        values = (center_A, force_constant, outer_extent_A)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Restraint parameters must be finite")
        if force_constant <= 0 or outer_extent_A <= 0:
            raise ValueError("Force constant and outer extent must be positive")
        first, second, third, fourth = self.atom_ids
        atom_list = f"{first},{second},{third},{fourth},0,0,0,0"
        low, high = center_A - outer_extent_A, center_A + outer_extent_A
        return (
            f"&rst iat={atom_list}, rstwt=1,-1,0,0, ialtd=0, "
            f"r1={low:.12g}, r2={center_A:.12g}, r3={center_A:.12g}, "
            f"r4={high:.12g}, rk2={force_constant:.12g}, "
            f"rk3={force_constant:.12g}, /\n"
        )


FEATURES: Mapping[str, DistanceDifference] = {
    "q_acyl": DistanceDifference("q_acyl", (19916, 19929, 19916, 13916)),
    "q_act": DistanceDifference("q_act", (13916, 13917, 13905, 13917)),
    "q_del": DistanceDifference("q_del", (13905, 13917, 19929, 13917)),
    # Preserve the closed campaign's q_H definition and default ordering.
    "q_h": DistanceDifference("q_h", (13916, 13917, 19929, 13917)),
}
RAW_FEATURE_NAMES = ("q_acyl", "q_act", "q_del")
DEFAULT_BIAS_FEATURES = ("q_acyl", "q_h")


@dataclass(frozen=True)
class RelayBiasWindow:
    """An immutable pair of named Amber distance-difference restraints."""

    feature_names: tuple[str, str] = DEFAULT_BIAS_FEATURES
    centers_A: tuple[float, float] = (0.0, 0.0)
    force_constants: tuple[float, float] = (100.0, 100.0)
    outer_extent_A: float = 10.0

    def __post_init__(self) -> None:
        if len(self.feature_names) != 2 or self.feature_names[0] == self.feature_names[1]:
            raise ValueError("Each window needs two distinct named CVs")
        if any(name not in FEATURES for name in self.feature_names):
            raise ValueError("Window refers to an unsupported distance-difference CV")
        numbers = (*self.centers_A, *self.force_constants, self.outer_extent_A)
        if not all(math.isfinite(value) for value in numbers):
            raise ValueError("Window parameters must be finite")
        if min(*self.force_constants, self.outer_extent_A) <= 0:
            raise ValueError("Force constants and outer extent must be positive")

    def restraint_text(self) -> str:
        """Return both native Amber restraint records in declared CV order."""
        return "".join(
            FEATURES[name].amber_restraint(center, force, self.outer_extent_A)
            for name, center, force in zip(
                self.feature_names,
                self.centers_A,
                self.force_constants,
                strict=True,
            )
        )


def raw_features(
    xyz: np.ndarray, box: Sequence[float] | None = None
) -> dict[str, float]:
    """Evaluate all three relay features and legacy q_H on one frame.

    Supply the Amber cell to use the same minimum-image helper as the existing
    umbrella observables. Omitting it retains the explicit nonperiodic fixture
    API.
    """
    return {
        name: FEATURES[name].value(xyz, box)
        for name in (*RAW_FEATURE_NAMES, "q_h")
    }


def iter_raw_features(
    coordinate_frames: Iterable[np.ndarray],
    boxes: Iterable[Sequence[float]] | None = None,
) -> Iterator[dict[str, float]]:
    """Yield all raw features, optionally minimum-image corrected per frame."""
    if boxes is None:
        for xyz in coordinate_frames:
            yield raw_features(xyz)
        return
    for xyz, box in zip(coordinate_frames, boxes, strict=True):
        yield raw_features(xyz, box)


def amber_bias_energy(value_A: float, center_A: float, force: float, extent_A: float) -> float:
    """Return native Amber's harmonic center and tangent-linear tail energy."""
    if not all(math.isfinite(value) for value in (value_A, center_A, force, extent_A)):
        raise ValueError("Bias inputs must be finite")
    if force <= 0 or extent_A <= 0:
        raise ValueError("Force and extent must be positive")
    delta = abs(value_A - center_A)
    if delta <= extent_A:
        return force * delta**2
    return force * extent_A * (2.0 * delta - extent_A)


def cross_bias_energies(
    features_by_frame: Mapping[str, Sequence[float]],
    windows: Sequence[RelayBiasWindow],
) -> np.ndarray:
    """Evaluate every window bias on every frame, in kcal/mol.

    All features named by any window must be provided as same-length finite
    arrays. The returned matrix is indexed ``[window, frame]``; multiply by
    beta to obtain the bias-only reduced-potential matrix for MBAR.
    """
    if not windows:
        raise ValueError("At least one bias window is required")
    names = {name for window in windows for name in window.feature_names}
    missing = names - features_by_frame.keys()
    if missing:
        raise ValueError(f"Missing frame features for cross-bias: {sorted(missing)}")
    arrays = {name: np.asarray(features_by_frame[name], dtype=float) for name in names}
    lengths = {array.size for array in arrays.values()}
    if len(lengths) != 1 or not lengths or any(array.ndim != 1 for array in arrays.values()):
        raise ValueError("Each cross-bias feature must be a same-length 1D array")
    if any(not np.isfinite(array).all() for array in arrays.values()):
        raise ValueError("Cross-bias features must be finite")
    frame_count = lengths.pop()
    energy = np.zeros((len(windows), frame_count), dtype=float)
    for row, window in enumerate(windows):
        for name, center, force in zip(
            window.feature_names, window.centers_A, window.force_constants, strict=True
        ):
            energy[row] += np.fromiter(
                (
                    amber_bias_energy(value, center, force, window.outer_extent_A)
                    for value in arrays[name]
                ),
                dtype=float,
                count=frame_count,
            )
    return energy
