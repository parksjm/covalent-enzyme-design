"""Bounded display selection and independent cpptraj DCD checks."""

from __future__ import annotations
import hashlib
import struct
from pathlib import Path
from typing import Any, BinaryIO
import numpy as np
from parmed.residue import AminoAcidResidue


def sha(path: Path) -> str:
    """Return an artifact SHA256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(handle: BinaryIO) -> bytes:
    """Read one little-endian Fortran record from cpptraj's DCD output."""
    header = handle.read(4)
    if len(header) != 4:
        raise ValueError("Incomplete DCD record")
    size = struct.unpack("<i", header)[0]
    data = handle.read(size)
    if handle.read(4) != header or len(data) != size:
        raise ValueError("DCD record marker mismatch")
    return data


def validate_dcd(
    path: Path, states: list[tuple[np.ndarray, np.ndarray]]
) -> dict[str, Any]:
    """Independently check all cpptraj DCD frames against source coordinates/boxes."""
    max_error = 0.0
    with path.open("rb") as handle:
        header = record(handle)
        if len(header) != 84 or header[:4] != b"CORD":
            raise ValueError("Unsupported DCD header")
        control = struct.unpack("<20i", header[4:])
        if control[0] != len(states) or not control[10]:
            raise ValueError("DCD frame count or periodic box flag mismatch")
        record(handle)
        atoms = struct.unpack("<i", record(handle))[0]
        if atoms != len(states[0][0]):
            raise ValueError("DCD atom count mismatch")
        for expected, expected_box in states:
            cell = np.frombuffer(record(handle), dtype="<f8")
            angles = cell[[4, 3, 1]].copy()
            if np.all(np.abs(angles) <= 1):
                angles = np.degrees(np.arccos(angles))
            box = np.r_[cell[[0, 2, 5]], angles]
            xyz = np.stack(
                [np.frombuffer(record(handle), dtype="<f4") for _ in range(3)], axis=1
            )
            if not np.isfinite(xyz).all() or not np.array_equal(
                xyz, expected.astype(np.float32)
            ):
                raise ValueError(
                    "DCD does not match source float32 display coordinates"
                )
            if not np.allclose(box, expected_box, atol=5e-5, rtol=0):
                raise ValueError("DCD periodic box differs from source")
            max_error = max(
                max_error, float(np.max(np.abs(xyz.astype(float) - expected)))
            )
        if handle.read(1):
            raise ValueError("Unexpected extra DCD records")
    return {
        "frames": len(states),
        "atoms": atoms,
        "maximum_coordinate_rounding_error_A": max_error,
        "float32_coordinates_verified": True,
        "all_periodic_boxes_verified": True,
    }


def select_display_atoms(
    parm: Any, states: list, qm_indices: list[int]
) -> tuple[np.ndarray, list[int]]:
    """Keep protein/LIG and a fixed whole-water union within 4 A of QM heavy atoms."""
    waters = [res for res in parm.residues if res.name == "WAT"]
    oxygens = []
    for residue in waters:
        oxygen = [atom.idx for atom in residue.atoms if atom.atomic_number == 8]
        if (
            len(oxygen) != 1
            or len(residue.atoms) != 3
            or sum(atom.atomic_number == 1 for atom in residue.atoms) != 2
        ):
            raise ValueError("Expected complete three-atom WAT residues")
        oxygens.append(oxygen[0])
    heavy = [index for index in qm_indices if parm.atoms[index].atomic_number != 1]
    retained = np.zeros(len(waters), dtype=bool)
    for xyz, _ in states:
        water_xyz = xyz[oxygens]
        for index in heavy:
            retained |= np.sum((water_xyz - xyz[index]) ** 2, axis=1) <= 16.0
    water_ids = [res.idx + 1 for res, keep in zip(waters, retained) if keep]
    residues = set(water_ids)
    for residue in parm.residues:
        if AminoAcidResidue.has(residue.name) or residue.name == "LIG":
            residues.add(residue.idx + 1)
        elif residue.name not in {"WAT", "Na+", "Cl-", "NA", "CL"}:
            raise ValueError(f"Undeclared nonprotein component: {residue.name}")
    atoms = np.array(
        [atom.idx for atom in parm.atoms if atom.residue.idx + 1 in residues]
    )
    if not set(qm_indices) <= set(atoms):
        raise ValueError("Stripping removed a QM atom")
    return atoms, water_ids
