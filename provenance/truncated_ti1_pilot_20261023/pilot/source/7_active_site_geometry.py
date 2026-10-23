"""Measure Nyl12 active-site geometry using the four-site reference math.

Mathematical provenance: ``/data/jpc/enzyme-substrate-modeling/scripts/lib/nyl12_four_site_geometry_math.py``.
The source file SHA256 is calculated for each report and included as
``geometry_source_sha256`` so the exact borrowed implementation can be audited.
The metric definitions and signed dihedral convention are ported from that
source; no geometric acceptance thresholds are applied here.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Mapping

import numpy as np

Vector = tuple[float, float, float]
Coordinates = Mapping[str, Vector]

GEOMETRY_SOURCE = Path(
    "/data/jpc/enzyme-substrate-modeling/scripts/lib/nyl12_four_site_geometry_math.py"
)
TARGET_ATOMS: dict[str, dict[str, str]] = {
    "pa6": {
        "carbonyl_c": "C54",
        "carbonyl_o": "O49",
        "amide_n": "N59",
        "adjacent_c": "C62",
    },
    "poly": {
        "carbonyl_c": "C033",
        "carbonyl_o": "O034",
        "amide_n": "N035",
        "adjacent_c": "C032",
    },
}
ACTIVE_SITE_ATOMS = {
    "thr_og1": "Thr929.OG1",
    "phe_n": "Phe113.N",
    "asn_nd2": "Asn186.ND2",
}
AMBER_RESIDUES = {
    "thr_og1": (929, "THR", "OG1"),
    "phe_n": (113, "PHE", "N"),
    "asn_nd2": (186, "ASN", "ND2"),
}


def _sub(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vector, b: Vector) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vector, b: Vector) -> Vector:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vector) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Vector, label: str) -> Vector:
    length = _norm(a)
    if not math.isfinite(length) or length <= 1e-12:
        raise ValueError(f"Degenerate {label}")
    return (a[0] / length, a[1] / length, a[2] / length)


def _coordinate(coordinates: Coordinates, key: str) -> Vector:
    if key not in coordinates:
        raise ValueError(f"Missing required atom coordinate: {key}")
    point = coordinates[key]
    if len(point) != 3 or not all(math.isfinite(float(v)) for v in point):
        raise ValueError(f"Atom coordinate must contain three finite values: {key}")
    return (float(point[0]), float(point[1]), float(point[2]))


def distance(a: Vector, b: Vector) -> float:
    """Return Euclidean distance in the coordinate unit (normally angstroms)."""
    value = _norm(_sub(a, b))
    if not math.isfinite(value):
        raise ValueError("Distance is not finite")
    return value


def angle_degrees(a: Vector, vertex: Vector, c: Vector) -> float:
    """Return the angle A-vertex-C in degrees, failing on zero-length arms."""
    left, right = _sub(a, vertex), _sub(c, vertex)
    cosine = max(
        -1.0, min(1.0, _dot(_unit(left, "angle arm"), _unit(right, "angle arm")))
    )
    return math.degrees(math.acos(cosine))


def projected_herschlag_phi(
    og1: Vector, oxygen: Vector, carbon: Vector, nitrogen: Vector, adjacent_c: Vector
) -> float:
    """Return signed projected Herschlag phi in degrees.

    The carbonyl carbon is projected onto the plane (O, substrate N, adjacent C),
    then the reference signed-dihedral atan2 convention is applied to
    OG1, O, projected carbon, adjacent C. Degenerate planes or dihedrals fail.
    """
    plane_a, plane_b = _sub(nitrogen, oxygen), _sub(adjacent_c, oxygen)
    normal = _unit(_cross(plane_a, plane_b), "Herschlag projection plane")
    carbon_from_o = _sub(carbon, oxygen)
    normal_component = _dot(carbon_from_o, normal)
    projected = _sub(
        carbon_from_o,
        tuple(component * normal_component for component in normal),
    )
    projected_carbon = (
        oxygen[0] + projected[0],
        oxygen[1] + projected[1],
        oxygen[2] + projected[2],
    )

    v1 = _sub(oxygen, og1)
    v2 = _sub(projected_carbon, oxygen)
    v3 = _sub(adjacent_c, projected_carbon)
    n1, n2 = _cross(v1, v2), _cross(v2, v3)
    if _norm(n1) <= 1e-12 or _norm(n2) <= 1e-12:
        raise ValueError("Degenerate projected Herschlag dihedral")
    x = _dot(n1, n2)
    y = _dot(_cross(n1, _unit(v2, "Herschlag dihedral central bond")), n2)
    if (
        not math.isfinite(x)
        or not math.isfinite(y)
        or (abs(x) <= 1e-12 and abs(y) <= 1e-12)
    ):
        raise ValueError("Degenerate projected Herschlag dihedral")
    return math.degrees(math.atan2(y, x))


def source_sha256(path: Path = GEOMETRY_SOURCE) -> str:
    """Return SHA256 for the exact reference math source file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coordinates_from_amber(topology: object, substrate: str) -> dict[str, Vector]:
    """Validate an Amber model and collect the unique active-site coordinates.

    Residue numbers are one-based Amber residue positions. The imported model
    layout has 1,314 residues, with the covalent ligand at residue 1,313.
    """
    if substrate not in TARGET_ATOMS:
        raise ValueError(f"Unsupported substrate: {substrate}")
    residues = list(topology.residues)
    atoms = list(topology.atoms)
    xyz = np.asarray(topology.coordinates, dtype=float)
    if len(residues) != 1314:
        raise ValueError(f"expected 1314 Amber residues, found {len(residues)}")
    if xyz.shape != (len(atoms), 3):
        raise ValueError(
            f"coordinates must have shape ({len(atoms)}, 3), got {xyz.shape}"
        )
    if not np.isfinite(xyz).all():
        raise ValueError("coordinates contain non-finite values")
    found: dict[str, Vector] = {}
    requirements = dict(AMBER_RESIDUES)
    requirements.update(
        {
            "substrate": (
                1313,
                "LIG",
                TARGET_ATOMS[substrate]["carbonyl_c"],
            )
        }
    )
    for label, (
        residue_number,
        expected_resname,
        required_atom,
    ) in requirements.items():
        residue = residues[residue_number - 1]
        residue_name = str(residue.name).strip().upper()
        if residue_name != expected_resname:
            raise ValueError(
                f"Amber residue {residue_number} must be {expected_resname}; "
                f"found {residue_name}"
            )
        residue_atoms = list(residue.atoms)
        names = [str(atom.name).strip() for atom in residue_atoms]
        if names.count(required_atom) != 1:
            raise ValueError(
                f"Amber residue {residue_number} must contain exactly one "
                f"{required_atom}; found {names.count(required_atom)}"
            )
        atom = residue_atoms[names.index(required_atom)]
        found[
            ACTIVE_SITE_ATOMS[label]
            if label != "substrate"
            else f"substrate.{required_atom}"
        ] = tuple(float(value) for value in xyz[atom.idx])

    ligand = residues[1312]
    ligand_atoms = list(ligand.atoms)
    for semantic, atom_name in TARGET_ATOMS[substrate].items():
        names = [str(atom.name).strip() for atom in ligand_atoms]
        if names.count(atom_name) != 1:
            raise ValueError(
                f"Amber residue 1313 must contain exactly one {atom_name}; "
                f"found {names.count(atom_name)}"
            )
        atom = ligand_atoms[names.index(atom_name)]
        found[f"substrate.{atom_name}"] = tuple(float(value) for value in xyz[atom.idx])
    return found


def measure_active_site(
    coordinates: Coordinates,
    substrate: str,
    state: str,
    source_path: Path = GEOMETRY_SOURCE,
) -> dict[str, object]:
    """Measure active-site distances, attack angle, and projected phi.

    Coordinates use keys ``Thr929.OG1``, ``Phe113.N``, ``Asn186.ND2`` and
    ``substrate.<atom-name>``. Supported substrate labels are ``pa6`` and
    ``poly``. State must explicitly be ``ti1`` or ``es``; the report records
    that this is descriptive geometry and carries no TI1/ES classification.
    """
    if state not in {"ti1", "es"}:
        raise ValueError("state must be explicitly set to 'ti1' or 'es'")
    if substrate not in TARGET_ATOMS:
        raise ValueError(f"Unsupported substrate: {substrate}")
    target = TARGET_ATOMS[substrate]
    og1 = _coordinate(coordinates, ACTIVE_SITE_ATOMS["thr_og1"])
    phe_n = _coordinate(coordinates, ACTIVE_SITE_ATOMS["phe_n"])
    asn_nd2 = _coordinate(coordinates, ACTIVE_SITE_ATOMS["asn_nd2"])
    carbon = _coordinate(coordinates, f"substrate.{target['carbonyl_c']}")
    oxygen = _coordinate(coordinates, f"substrate.{target['carbonyl_o']}")
    nitrogen = _coordinate(coordinates, f"substrate.{target['amide_n']}")
    adjacent_c = _coordinate(coordinates, f"substrate.{target['adjacent_c']}")
    return {
        "substrate": substrate,
        "state": state,
        "state_interpretation": f"Explicitly supplied state={state}; metrics do not classify TI1 versus ES.",
        "atom_names": {**ACTIVE_SITE_ATOMS, **target},
        "units": {"distances": "angstrom", "angles": "degree"},
        "metrics": {
            "d1_thr929_og1_scissile_c": distance(og1, carbon),
            "d2_phe113_n_scissile_o": distance(phe_n, oxygen),
            "d3_asn186_nd2_scissile_o": distance(asn_nd2, oxygen),
            "a_BD_og1_c_o": angle_degrees(og1, carbon, oxygen),
            "herschlag_phi_projected": projected_herschlag_phi(
                og1, oxygen, carbon, nitrogen, adjacent_c
            ),
        },
        "geometry_source_path": str(source_path),
        "geometry_source_sha256": source_sha256(source_path),
    }
