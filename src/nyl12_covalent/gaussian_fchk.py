"""Strict reader for the Gaussian formatted-checkpoint fields used here."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

LOGGER = logging.getLogger(__name__)
FloatArray: TypeAlias = NDArray[np.float64]
IntArray: TypeAlias = NDArray[np.int64]

_HEADER = re.compile(r"^\s*([IRCL])\s+(?:N=\s*(\d+)|(.*?))\s*$")
_INT = re.compile(r"^[+-]?\d+$")
_FLOAT = re.compile(r"^[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[DEde][+-]?\d+)?$")
_REQUIRED = {
    "numberofatoms",
    "charge",
    "multiplicity",
    "atomicnumbers",
    "currentcartesiancoordinates",
    "cartesianforceconstants",
}


@dataclass(frozen=True)
class GaussianFchk:
    """Selected molecular data from an fchk file, in Gaussian-native units.

    Attributes:
        atomic_numbers: Atomic numbers, one per atom.
        coordinates_bohr: Cartesian coordinates with shape ``(N, 3)`` in bohr.
        hessian_eh_per_bohr2: Cartesian Hessian with shape ``(3N, 3N)`` in
            Hartree per bohr squared. Gaussian's packed lower triangle is
            expanded symmetrically without mass weighting.
        charge: Molecular charge in elementary charge units.
        multiplicity: Spin multiplicity.
        gradient_eh_per_bohr: Optional Cartesian gradient with shape ``(N, 3)``
            in Hartree per bohr; ``None`` when the record is absent.
    """

    atomic_numbers: IntArray
    coordinates_bohr: FloatArray
    hessian_eh_per_bohr2: FloatArray
    charge: int
    multiplicity: int
    gradient_eh_per_bohr: FloatArray | None = None


def _label_key(label: str) -> str:
    """Normalize a fixed-record label for matching and duplicate detection."""
    return "".join(label.split()).lower()


def _record_header(line: str) -> tuple[str, str, int | str] | None:
    """Parse a fixed-column fchk header, returning label, type, and payload."""
    if len(line) < 44:
        return None
    label = line[:40].strip()
    match = _HEADER.fullmatch(line[40:])
    if not label or match is None:
        return None
    kind, count, scalar = match.groups()
    return label, kind, int(count) if count is not None else (scalar or "")


def _as_integer(token: str, label: str) -> int:
    if not _INT.fullmatch(token):
        raise ValueError(f"{label} contains a malformed integer: {token!r}")
    return int(token)


def _as_real(token: str, label: str) -> float:
    if not _FLOAT.fullmatch(token):
        raise ValueError(f"{label} contains a malformed real: {token!r}")
    value = float(token.replace("D", "E").replace("d", "e"))
    if not np.isfinite(value):
        raise ValueError(f"{label} contains a non-finite real")
    return value


def _parse_records(text: str) -> dict[str, tuple[str, list[str]]]:
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError("fchk must contain title, route, and data records")
    records: dict[str, tuple[str, list[str]]] = {}
    index = 2
    while index < len(lines):
        parsed = _record_header(lines[index])
        if parsed is None:
            raise ValueError(f"invalid fchk record header on line {index + 1}")
        label, kind, payload = parsed
        key = _label_key(label)
        required = key in _REQUIRED or key == "cartesiangradient"
        if required and key in records:
            raise ValueError(f"duplicate fchk field {label!r}")
        index += 1
        if isinstance(payload, int):
            if payload < 0:
                raise ValueError(f"negative array length for {label!r}")
            values: list[str] = []
            value_count = 0
            while value_count < payload:
                if index >= len(lines):
                    raise ValueError(f"truncated array for {label!r}")
                if _record_header(lines[index]) is not None:
                    raise ValueError(f"truncated array for {label!r}")
                if kind == "C":
                    line_values = [
                        lines[index][start : start + 12].strip()
                        for start in range(0, len(lines[index]), 12)
                    ]
                else:
                    line_values = lines[index].split()
                if not line_values:
                    index += 1
                    continue
                remaining = payload - value_count
                if len(line_values) > remaining:
                    extra = line_values[remaining:]
                    if kind != "C" or any(extra):
                        raise ValueError(f"too many values for array {label!r}")
                    line_values = line_values[:remaining]
                value_count += len(line_values)
                if required:
                    values.extend(line_values)
                index += 1
        else:
            scalar_values = [payload] if kind == "C" else payload.split()
            if kind != "C" and len(scalar_values) != 1:
                raise ValueError(f"invalid scalar value for {label!r}")
            values = scalar_values
        if required:
            records[key] = (kind, values)
    missing = _REQUIRED - records.keys()
    if missing:
        raise ValueError(f"missing required fchk fields: {', '.join(sorted(missing))}")
    return records


def _optional_gradient(
    records: dict[str, tuple[str, list[str]]], dimension: int
) -> FloatArray | None:
    """Parse an optional Cartesian gradient in Hartree per bohr."""
    key = "cartesiangradient"
    if key not in records:
        return None
    kind, values = records[key]
    if kind != "R":
        raise ValueError("Cartesian Gradient must use real-valued records")
    if len(values) != dimension:
        raise ValueError(
            f"Cartesian Gradient has {len(values)} values; expected {dimension}"
        )
    gradient = np.asarray(
        [_as_real(token, "Cartesian Gradient") for token in values],
        dtype=np.float64,
    )
    return gradient.reshape((dimension // 3, 3))


def parse_gaussian_fchk(
    text: str,
    *,
    expected_atomic_numbers: tuple[int, ...] | None = None,
    expected_charge: int | None = None,
    expected_multiplicity: int | None = None,
) -> GaussianFchk:
    """Parse required geometry and Hessian fields from fchk text.

    Args:
        text: Complete Gaussian formatted-checkpoint text.
        expected_atomic_numbers: Optional exact ordered atomic-number check.
        expected_charge: Optional required molecular charge.
        expected_multiplicity: Optional required spin multiplicity.

    Returns:
        Parsed atom numbers, coordinates, Cartesian Hessian, and state.

    Raises:
        ValueError: If records are absent, malformed, inconsistent, or fail
            any supplied molecular identity checks.
    """
    records = _parse_records(text)

    def scalar_int(key: str) -> int:
        kind, values = records[key]
        if kind != "I" or len(values) != 1:
            raise ValueError(f"{key} must be an integer scalar")
        return _as_integer(values[0], key)

    natoms = scalar_int("numberofatoms")
    charge = scalar_int("charge")
    multiplicity = scalar_int("multiplicity")
    if natoms < 1 or multiplicity < 1:
        raise ValueError("number of atoms and multiplicity must be positive")
    atom_kind, atom_values = records["atomicnumbers"]
    coord_kind, coord_values = records["currentcartesiancoordinates"]
    hess_kind, hess_values = records["cartesianforceconstants"]
    if atom_kind != "I" or coord_kind != "R" or hess_kind != "R":
        raise ValueError("atomic numbers must be I; coordinates and Hessian must be R")
    atomic_numbers = np.asarray(
        [_as_integer(value, "Atomic numbers") for value in atom_values], dtype=np.int64
    )
    coordinates = np.asarray(
        [_as_real(value, "Current cartesian coordinates") for value in coord_values],
        dtype=np.float64,
    )
    packed = np.asarray(
        [_as_real(value, "Cartesian Force Constants") for value in hess_values],
        dtype=np.float64,
    )
    dimension = 3 * natoms
    expected_packed = dimension * (dimension + 1) // 2
    if atomic_numbers.size != natoms:
        raise ValueError("Atomic numbers length does not match Number of atoms")
    if coordinates.size != dimension:
        raise ValueError("coordinate length does not match 3N")
    if packed.size != expected_packed:
        raise ValueError("Cartesian Force Constants length is not the 3N triangle")
    if np.any(atomic_numbers <= 0):
        raise ValueError("atomic numbers must be positive")
    if expected_atomic_numbers is not None and tuple(atomic_numbers) != tuple(
        expected_atomic_numbers
    ):
        raise ValueError("atomic numbers do not match the expected ordered identity")
    if expected_charge is not None and charge != expected_charge:
        raise ValueError(f"charge {charge} does not match expected {expected_charge}")
    if expected_multiplicity is not None and multiplicity != expected_multiplicity:
        raise ValueError(
            f"multiplicity {multiplicity} does not match expected {expected_multiplicity}"
        )

    hessian = np.zeros((dimension, dimension), dtype=np.float64)
    lower = np.tril_indices(dimension)
    hessian[lower] = packed
    hessian[(lower[1], lower[0])] = packed
    LOGGER.debug("Parsed Gaussian fchk for %d atoms", natoms)
    return GaussianFchk(
        atomic_numbers=atomic_numbers,
        coordinates_bohr=coordinates.reshape((natoms, 3)),
        hessian_eh_per_bohr2=hessian,
        charge=charge,
        multiplicity=multiplicity,
        gradient_eh_per_bohr=_optional_gradient(records, dimension),
    )


def read_gaussian_fchk(
    path: str | Path,
    *,
    expected_atomic_numbers: tuple[int, ...] | None = None,
    expected_charge: int | None = None,
    expected_multiplicity: int | None = None,
) -> GaussianFchk:
    """Read and validate a Gaussian formatted-checkpoint file.

    Args:
        path: Path to a UTF-8 fchk file.
        expected_atomic_numbers: Optional exact ordered atomic-number check.
        expected_charge: Optional required molecular charge.
        expected_multiplicity: Optional required spin multiplicity.

    Returns:
        Parsed data in bohr and Hartree per bohr squared.

    Raises:
        OSError: If the file cannot be read.
        ValueError: If its contents fail strict fchk or identity validation.
    """
    return parse_gaussian_fchk(
        Path(path).read_text(encoding="utf-8"),
        expected_atomic_numbers=expected_atomic_numbers,
        expected_charge=expected_charge,
        expected_multiplicity=expected_multiplicity,
    )
