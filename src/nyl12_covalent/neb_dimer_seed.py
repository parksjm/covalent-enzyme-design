"""Helpers for preparing a local dimer seed from adjacent native NEB images."""

from __future__ import annotations

import re
from typing import Any

import numpy as np


def final_neb_rows(
    mdout: str, expected_cycle: int, image_count: int
) -> list[tuple[int, float, float, float]]:
    """Validate and return energies and forces from the final native NEB report.

    Args:
        mdout: Complete native optimizer output text.
        expected_cycle: Required final cycle number.
        image_count: Number of regular images in the band.

    Returns:
        Tuples of one-based image, Hartree energy, tangential force, and
        perpendicular force, in native report order.

    Raises:
        ValueError: If the final report or its ordered save loop is incomplete.
    """
    reports = list(re.finditer(r"NEB Report\n", mdout))
    if not reports:
        raise ValueError("No complete NEB report found")
    start = reports[-1].end()
    end = mdout.find("Total path length", start)
    if end < 0:
        raise ValueError("Final NEB report is incomplete")
    rows = [
        (int(image), float(energy), float(tangent), float(perpendicular))
        for image, energy, tangent, perpendicular in re.findall(
            r"^Img\s+(\d+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)",
            mdout[start:end],
            re.M,
        )
    ]
    if [row[0] for row in rows] != list(range(1, image_count + 1)):
        raise ValueError("Final NEB report has an invalid image mapping")
    prefix = mdout[: reports[-1].start()]
    cycles = re.findall(r"No climbing image used in cycle\s+(\d+)", prefix)
    if not cycles or int(cycles[-1]) != expected_cycle:
        raise ValueError("Final NEB report cycle does not match the expected cycle")
    tail = mdout[end:]
    closure = re.search(rf"Testing convergence\s+in cycle\s+{expected_cycle}\b", tail)
    if closure is None:
        raise ValueError("Final NEB save loop has no convergence closure")
    saved = re.findall(r"Saving img\s+(\d+) to restart:", tail[: closure.start()])
    if [int(image) for image in saved] != list(range(1, image_count + 1)):
        raise ValueError("Final NEB save loop has an invalid image mapping")
    return rows


def build_local_dimer_seed(
    frames: list[np.ndarray],
    active_indices: list[int],
    rotation_indices: np.ndarray,
    left_image: int = 7,
    right_image: int = 8,
    displacement_angstrom: float = 0.01,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Interpolate adjacent images and make a normalized local dimer mode.

    Image numbers are one-based. The seed is the midpoint; the mode points
    forward from the left image toward the right image.

    Args:
        frames: Full-system coordinates reconstructed from native XYZ.
        active_indices: Zero-based topology indices mapped by the active mask.
        rotation_indices: Zero-based topology indices used to normalize the mode.
        left_image: One-based starting image.
        right_image: One-based ending image.
        displacement_angstrom: Mode displacement magnitude on rotation indices.

    Returns:
        Midpoint coordinates, displaced mode coordinates, and tangent norm.

    Raises:
        ValueError: If the mapping, image pair, masks, or displacement is invalid.
    """
    if not frames or not active_indices:
        raise ValueError("Frames and active atom mapping are required")
    if right_image != left_image + 1 or left_image < 1 or right_image > len(frames):
        raise ValueError("Images must be an adjacent in-range one-based pair")
    if not np.isfinite(displacement_angstrom) or displacement_angstrom <= 0:
        raise ValueError("Mode displacement must be finite and positive")
    if len(set(active_indices)) != len(active_indices):
        raise ValueError("Active atom mapping contains duplicate topology indices")
    if len(set(rotation_indices.tolist())) != len(rotation_indices):
        raise ValueError("Rotation mask contains duplicate topology indices")
    if not len(rotation_indices) or not set(rotation_indices).issubset(active_indices):
        raise ValueError("Rotation mask must be nonempty and contained in active atoms")
    shape = frames[0].shape
    if (
        len(shape) != 2
        or shape[1] != 3
        or any(frame.shape != shape for frame in frames)
    ):
        raise ValueError("Frame coordinate shapes are inconsistent")
    if not np.isfinite(np.asarray(frames)).all():
        raise ValueError("Frame coordinates must be finite")
    tangent = frames[right_image - 1] - frames[left_image - 1]
    norm = float(np.linalg.norm(tangent[rotation_indices]))
    if not np.isfinite(norm) or norm <= 1e-10:
        raise ValueError("Local NEB tangent has zero or invalid norm")
    seed = 0.5 * (frames[left_image - 1] + frames[right_image - 1])
    mode = seed + displacement_angstrom * tangent / norm
    return seed, mode, norm


def validate_reactive_atom_mapping(
    topology: Any,
    active_indices: list[int],
    declared_mapping: dict[str, dict[str, int | str]],
    required_mapping: dict[int, dict[str, int | str]],
) -> dict[int, int]:
    """Validate declared atom identities against topology serials and mask.

    Args:
        topology: ParmEd topology with ordered atoms and residues.
        active_indices: Zero-based indices selected by the native active mask.
        declared_mapping: Saved preparation mapping keyed by one-based serial.
        required_mapping: Required one-based serial identities and elements.

    Returns:
        One-based serials mapped to zero-based topology indices.

    Raises:
        ValueError: If any declared or required atom identity is inconsistent.
    """
    combined = dict(declared_mapping)
    combined.update(
        {str(serial): identity for serial, identity in required_mapping.items()}
    )
    active = set(active_indices)
    validated: dict[int, int] = {}
    for serial_text, identity in combined.items():
        serial = int(serial_text)
        index = serial - 1
        if index < 0 or index >= len(topology.atoms) or index not in active:
            raise ValueError(f"Mapped reactive atom {serial} is absent or inactive")
        atom = topology.atoms[index]
        expected_residue = int(identity["residue"])
        if (
            atom.name != identity["name"]
            or atom.residue.idx + 1 != expected_residue
            or (
                "atomic_number" in identity
                and atom.atomic_number != int(identity["atomic_number"])
            )
        ):
            raise ValueError(f"Topology identity mismatch for atom {serial}")
        validated[serial] = index
    return validated


def local_heavy_atom_clashes(
    topology: Any,
    coordinates: np.ndarray,
    reactive_indices: set[int],
    *,
    neighborhood_angstrom: float = 6.0,
    contact_cutoff_angstrom: float = 2.0,
    severe_cutoff_angstrom: float = 1.8,
) -> dict[str, Any]:
    """Report short local nonbonded heavy-atom contacts around reactive atoms.

    Direct bonds, 1-3 pairs, and pairings within the mapped reactive core are
    omitted because their distances are interpreted as reaction coordinates.

    Args:
        topology: ParmEd topology in coordinate order.
        coordinates: Full-system coordinates in Angstrom.
        reactive_indices: Zero-based topology indices defining the reaction core.
        neighborhood_angstrom: Radius around any reactive atom for local review.
        contact_cutoff_angstrom: Maximum reported heavy-atom contact distance.
        severe_cutoff_angstrom: Threshold used to flag severe contacts.

    Returns:
        JSON-serializable thresholds, severe status, and sorted contact pairs.

    Raises:
        ValueError: If coordinate shapes, indices, or cutoffs are invalid.
    """
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.shape != (len(topology.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Full-system coordinates do not match topology")
    if (
        not reactive_indices
        or min(reactive_indices) < 0
        or max(reactive_indices) >= len(xyz)
    ):
        raise ValueError("Reactive atom indices are empty or out of range")
    if (
        not 0
        < severe_cutoff_angstrom
        <= contact_cutoff_angstrom
        < neighborhood_angstrom
    ):
        raise ValueError("Clash audit cutoffs are invalid")
    heavy = [i for i, atom in enumerate(topology.atoms) if atom.atomic_number > 1]
    centers = xyz[list(reactive_indices)]
    local = [
        index
        for index in heavy
        if float(np.linalg.norm(centers - xyz[index], axis=1).min())
        <= neighborhood_angstrom
    ]
    bonded = {i: set() for i in range(len(topology.atoms))}
    for bond in topology.bonds:
        first, second = bond.atom1.idx, bond.atom2.idx
        bonded[first].add(second)
        bonded[second].add(first)
    excluded = {index: set(neighbors) for index, neighbors in bonded.items()}
    for index, neighbors in bonded.items():
        for neighbor in neighbors:
            excluded[index].update(bonded[neighbor])
    reactive_pairs = {
        frozenset((first, second))
        for first in reactive_indices
        for second in reactive_indices
        if first != second
    }
    contacts = []
    for position, first in enumerate(local):
        for second in local[position + 1 :]:
            if (
                second in excluded[first]
                or frozenset((first, second)) in reactive_pairs
            ):
                continue
            distance = float(np.linalg.norm(xyz[first] - xyz[second]))
            if distance < contact_cutoff_angstrom:
                contacts.append(
                    {
                        "distance_angstrom": distance,
                        "first": _atom_identity(topology.atoms[first], first),
                        "second": _atom_identity(topology.atoms[second], second),
                    }
                )
    contacts.sort(key=lambda contact: contact["distance_angstrom"])
    return {
        "neighborhood_angstrom": neighborhood_angstrom,
        "contact_cutoff_angstrom": contact_cutoff_angstrom,
        "severe_cutoff_angstrom": severe_cutoff_angstrom,
        "severe_contact_count": sum(
            pair["distance_angstrom"] < severe_cutoff_angstrom for pair in contacts
        ),
        "contacts": contacts,
        "interpretation": "Reactive-core pairs and 1-2/1-3 topology pairs omitted; review reported contacts before launch.",
    }


def _atom_identity(atom: Any, index: int) -> dict[str, int | str]:
    """Return concise atom identity for a contact report."""
    return {
        "serial": index + 1,
        "name": atom.name,
        "residue": atom.residue.idx + 1,
        "residue_name": atom.residue.name,
    }
