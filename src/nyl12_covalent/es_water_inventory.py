"""Read-only explicit-water geometry inventory for a mapped ES structure."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm

_ROLES = {
    "thr_n_terminus": ("N", 7),
    "nucleophile": ("OG1", 8),
    "asp272_od1": ("OD1", 8),
    "asp272_od2": ("OD2", 8),
    "asp274_od1": ("OD1", 8),
    "asp274_od2": ("OD2", 8),
}
_WATER_NAMES = {"WAT", "HOH", "TIP3", "TIP3P", "TP3"}


def _select_roles(parm: Any, roles: dict[str, dict[str, int | str]]) -> dict:
    if set(roles) != set(_ROLES):
        raise ValueError(f"Roles must be exactly {sorted(_ROLES)}")
    selected = {}
    for role, (name, element) in _ROLES.items():
        row = roles[role]
        if not isinstance(row, dict) or set(row) != {"residue", "resname", "atom"}:
            raise ValueError(f"Invalid explicit mapping for {role}")
        number = row["residue"]
        if type(number) is not int or not 1 <= number <= len(parm.residues):
            raise ValueError(f"Invalid topology residue for {role}")
        residue = parm.residues[number - 1]
        if not isinstance(row["resname"], str) or residue.name != row["resname"]:
            raise ValueError(f"Residue identity mismatch for {role}")
        if row["atom"] != name:
            raise ValueError(f"Expected atom {name} for {role}")
        atoms = [atom for atom in residue.atoms if atom.name == name]
        if len(atoms) != 1 or atoms[0].atomic_number != element:
            raise ValueError(f"Missing, ambiguous or wrong-element atom for {role}")
        selected[role] = atoms[0]
    if len({atom.idx for atom in selected.values()}) != len(selected):
        raise ValueError("Roles must identify distinct atoms")
    for left, right in [
        ("thr_n_terminus", "nucleophile"),
        ("asp272_od1", "asp272_od2"),
        ("asp274_od1", "asp274_od2"),
    ]:
        if selected[left].residue is not selected[right].residue:
            raise ValueError(f"Role pair {left}/{right} must share a residue")
    if len({atom.residue.idx for atom in selected.values()}) != 3:
        raise ValueError("Thr and the two Asp selections require three residues")
    return selected


def _hydrogens(atom: Any) -> list:
    return sorted(
        [partner for partner in atom.bond_partners if partner.atomic_number == 1],
        key=lambda partner: partner.idx,
    )


def _atom_record(atom: Any, xyz: np.ndarray) -> dict:
    return {
        "index": atom.idx + 1,
        "name": atom.name,
        "residue_index": atom.residue.idx + 1,
        "residue_name": atom.residue.name,
        "coordinates_A": xyz[atom.idx].tolist(),
    }


def _donor_geometry(donor: Any, acceptor: Any, xyz: np.ndarray) -> dict:
    observations = []
    for hydrogen in _hydrogens(donor):
        first = xyz[donor.idx] - xyz[hydrogen.idx]
        second = xyz[acceptor.idx] - xyz[hydrogen.idx]
        norm = float(np.linalg.norm(first) * np.linalg.norm(second))
        if norm <= 0:
            raise ValueError("Coincident atoms make a donor-H-acceptor angle undefined")
        angle = float(
            np.degrees(np.arccos(np.clip(np.dot(first, second) / norm, -1, 1)))
        )
        observations.append(
            {
                "hydrogen_index": hydrogen.idx + 1,
                "angle_degrees": angle,
                "hydrogen_acceptor_distance_A": float(np.linalg.norm(second)),
            }
        )
    return {
        "donor_index": donor.idx + 1,
        "acceptor_index": acceptor.idx + 1,
        "explicit_hydrogen_count": len(observations),
        "observations": observations,
        "best_angle_degrees": max(
            (item["angle_degrees"] for item in observations), default=None
        ),
        "angle_at_least_120_degrees_count": sum(
            item["angle_degrees"] >= 120 for item in observations
        ),
    }


def inventory_es_waters(
    topology: Path,
    coordinates: Path,
    roles: dict[str, dict[str, int | str]],
) -> dict:
    """Describe explicit local-water contacts without assigning a mechanism.

    Args:
        topology: Amber topology, read without modification.
        coordinates: Matching Amber restart coordinates in an imaged frame.
        roles: Exact one-based topology residue, residue-name and atom-name
            mappings for thr_n_terminus, nucleophile and each Asp OD1/OD2.

    Returns:
        JSON-compatible atom inventories, Cartesian distances and explicit
        donor-H-acceptor angles. Missing donor H yields a null best angle.

    Raises:
        FileNotFoundError: An input file does not exist.
        ValueError: Coordinates, mappings or a selected local water are invalid.
    """
    for path in (topology, coordinates):
        if not path.is_file():
            raise FileNotFoundError(path)
    parm = AmberParm(str(topology), xyz=str(coordinates))
    xyz = np.asarray(parm.coordinates, dtype=float)
    if xyz.shape != (len(parm.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Coordinates must be a finite array of shape (atoms, 3)")
    selected = _select_roles(parm, roles)
    nitrogen, og1 = selected["thr_n_terminus"], selected["nucleophile"]
    aspartates = {key: atom for key, atom in selected.items() if key.startswith("asp")}
    waters = []
    for residue in parm.residues:
        if residue.name not in _WATER_NAMES:
            continue
        oxygens = [atom for atom in residue.atoms if atom.atomic_number == 8]
        nearby = [
            oxygen
            for oxygen in oxygens
            if min(
                np.linalg.norm(xyz[oxygen.idx] - xyz[nitrogen.idx]),
                np.linalg.norm(xyz[oxygen.idx] - xyz[og1.idx]),
            )
            <= 4.0
        ]
        if not nearby:
            continue
        if len(oxygens) != 1:
            raise ValueError(f"Local water residue {residue.idx + 1} has ambiguous O")
        oxygen = nearby[0]
        hydrogens = _hydrogens(oxygen)
        if (
            len(hydrogens) != 2
            or len(residue.atoms) != 3
            or any(atom.residue is not residue for atom in oxygen.bond_partners)
            or set(oxygen.bond_partners) != set(hydrogens)
            or any(len(atom.bond_partners) != 1 for atom in hydrogens)
        ):
            raise ValueError(f"Local water residue {residue.idx + 1} lacks bonded H2O")
        distances = {
            key: float(np.linalg.norm(xyz[oxygen.idx] - xyz[atom.idx]))
            for key, atom in selected.items()
        }
        nearest = min(aspartates, key=lambda key: distances[key])
        waters.append(
            {
                "residue_index": residue.idx + 1,
                "oxygen": _atom_record(oxygen, xyz),
                "hydrogens": [_atom_record(atom, xyz) for atom in hydrogens],
                "explicit_hydrogen_count": len(hydrogens),
                "oxygen_distances_A": distances,
                "nearest_asp_oxygen": nearest,
                "nearest_asp_oxygen_distance_A": distances[nearest],
                "geometric_only_both_contacts_within_3_5_A": (
                    distances["thr_n_terminus"] <= 3.5
                    and distances["nucleophile"] <= 3.5
                ),
                "donor_h_acceptor_geometry": {
                    "thr_og1_h_water_o": _donor_geometry(og1, oxygen, xyz),
                    "water_o_h_thr_n": _donor_geometry(oxygen, nitrogen, xyz),
                    "thr_n_h_water_o": _donor_geometry(nitrogen, oxygen, xyz),
                    "water_o_h_thr_og1": _donor_geometry(oxygen, og1, xyz),
                },
            }
        )
    return {
        "index_convention": "one-based Amber topology atom/residue indices",
        "geometry_policy": "Cartesian imaged coordinates; no minimum-image correction",
        "interpretation": "Descriptive geometry only; no mechanism assigned",
        "selection_cutoff_A": 4.0,
        "selected_atoms": {
            key: _atom_record(atom, xyz) for key, atom in selected.items()
        },
        "local_water_count": len(waters),
        "geometric_only_candidate_count": sum(
            water["geometric_only_both_contacts_within_3_5_A"] for water in waters
        ),
        "waters": waters,
    }
