"""Read-only Amber enzyme–substrate connectivity and parameter inventory."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import parmed

_LABELS = {
    "nucleophile",
    "scissile_carbon",
    "scissile_oxygen",
    "leaving_nitrogen",
    "oxyanion_donor_backbone",
    "oxyanion_donor_sidechain",
    "thr_n_terminus",
}


def _finite(value: Any, description: str) -> float:
    if value is None:
        raise ValueError(f"Missing {description}")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite {description}")
    return number


def _identity(atom: Any) -> dict[str, Any]:
    return {
        "index": atom.idx + 1,
        "name": atom.name,
        "element": atom.element_name,
        "atomic_number": atom.atomic_number,
        "residue": atom.residue.idx + 1,
        "resname": atom.residue.name,
    }


def _inventory(atom: Any) -> dict[str, Any]:
    if not atom.name or not atom.type or atom.atomic_number <= 0:
        raise ValueError(f"Missing atom identity/type/element at atom {atom.idx + 1}")
    if not isinstance(atom.atom_type, parmed.topologyobjects.AtomType):
        raise ValueError(f"Missing LJ atom type at atom {atom.idx + 1}")
    rmin = _finite(atom.atom_type.rmin, "LJ rmin")
    epsilon = _finite(atom.atom_type.epsilon, "LJ epsilon")
    if rmin < 0 or epsilon < 0:
        raise ValueError(f"Negative LJ parameter at atom {atom.idx + 1}")
    neighbors = sorted(atom.bond_partners, key=lambda neighbor: neighbor.idx)
    return {
        **_identity(atom),
        "type": atom.type,
        "charge": _finite(atom.charge, "atomic charge"),
        "lj": {"rmin_angstrom": rmin, "epsilon_kcal_per_mol": epsilon},
        "bonded_neighbors": [_identity(neighbor) for neighbor in neighbors],
        "explicit_hydrogen_neighbors": [
            neighbor.idx + 1 for neighbor in neighbors if neighbor.atomic_number == 1
        ],
    }


def _residue(structure: Any, number: Any) -> Any:
    if type(number) is not int or not 1 <= number <= len(structure.residues):
        raise ValueError(f"Invalid 1-based topology residue index: {number!r}")
    return structure.residues[number - 1]


def audit_es_topology(
    topology_path: Path,
    atom_map: dict[str, dict[str, int | str]],
    ligand_residue: int,
) -> dict[str, Any]:
    """Inventory Amber ES parameters and observed bonds without modification.

    Args:
        topology_path: Amber topology file, read without coordinates.
        atom_map: Required role labels with exact residue (1-based topology
            index), resname, and atom keys for every row.
        ligand_residue: 1-based topology residue index of the intact substrate.

    Returns:
        Atomic and residue inventories, ligand graph, and observed H counts.
        Neither bond orders nor charge methods nor mechanistic states are inferred.

    Raises:
        ValueError: Invalid mapping, parameters, or ES connectivity.
        FileNotFoundError: Missing topology file.
    """
    if not topology_path.is_file():
        raise FileNotFoundError(topology_path)
    if not isinstance(atom_map, dict) or set(atom_map) != _LABELS:
        raise ValueError(f"atom_map must contain exactly {sorted(_LABELS)}")
    structure = parmed.amber.AmberParm(str(topology_path))
    ligand = _residue(structure, ligand_residue)
    selected: dict[str, Any] = {}
    for label, row in atom_map.items():
        if not isinstance(row, dict) or set(row) != {"residue", "resname", "atom"}:
            raise ValueError(f"Invalid atom-map row for {label}")
        residue = _residue(structure, row["residue"])
        if not isinstance(row["resname"], str) or residue.name != row["resname"]:
            raise ValueError(f"Residue name mismatch for {label}")
        if not isinstance(row["atom"], str) or not row["atom"]:
            raise ValueError(f"Invalid atom name for {label}")
        matches = [atom for atom in residue.atoms if atom.name == row["atom"]]
        if len(matches) != 1:
            raise ValueError(f"Missing or ambiguous atom for {label}: {row['atom']}")
        selected[label] = matches[0]
    carbon = selected["scissile_carbon"]
    for label in ("scissile_carbon", "scissile_oxygen", "leaving_nitrogen"):
        if selected[label].residue is not ligand:
            raise ValueError(f"{label} is outside the ligand residue")
    if selected["nucleophile"].residue is ligand:
        raise ValueError("Nucleophile belongs to the ligand residue")
    if selected["nucleophile"] in carbon.bond_partners:
        raise ValueError("Covalent nucleophile–scissile-carbon bond in ES topology")
    for label in ("scissile_oxygen", "leaving_nitrogen"):
        if selected[label] not in carbon.bond_partners:
            raise ValueError(f"Scissile carbon is not bonded to {label}")
    ligand_ids = {atom.idx for atom in ligand.atoms}
    remaining = set(ligand_ids)
    components: list[list[int]] = []
    while remaining:
        pending = [min(remaining)]
        component: set[int] = set()
        while pending:
            index = pending.pop()
            if index in component:
                continue
            component.add(index)
            pending.extend(
                atom.idx
                for atom in structure.atoms[index].bond_partners
                if atom.idx in ligand_ids and atom.idx not in component
            )
        remaining -= component
        components.append(sorted(index + 1 for index in component))
    if len(components) != 1:
        raise ValueError(f"Ligand has {len(components)} connected components")
    external = [
        (atom.idx + 1, partner.idx + 1)
        for atom in ligand.atoms
        for partner in atom.bond_partners
        if partner.idx not in ligand_ids
    ]
    if external:
        raise ValueError(f"Ligand has external covalent bonds: {external}")
    residues = []
    for number in sorted({atom.residue.idx + 1 for atom in selected.values()}):
        residue = _residue(structure, number)
        inventory = [_inventory(atom) for atom in residue.atoms]
        residues.append(
            {
                "residue": number,
                "resname": residue.name,
                "net_charge": math.fsum(atom["charge"] for atom in inventory),
                "atom_count": len(inventory),
                "explicit_hydrogen_count": sum(
                    atom["atomic_number"] == 1 for atom in inventory
                ),
            }
        )
    ligand_atoms = [_inventory(atom) for atom in ligand.atoms]
    return {
        "topology_path": str(topology_path),
        "selected_atoms": {label: _inventory(atom) for label, atom in selected.items()},
        "ligand": {
            "residue": ligand_residue,
            "resname": ligand.name,
            "atoms": ligand_atoms,
            "net_charge": math.fsum(atom["charge"] for atom in ligand_atoms),
            "connected_components": components,
            "bonds": sorted(
                [
                    [atom.idx + 1, partner.idx + 1]
                    for atom in ligand.atoms
                    for partner in atom.bond_partners
                    if atom.idx < partner.idx
                ]
            ),
        },
        "mapped_residues": residues,
        "protonation_observations": {
            label: {
                "explicit_hydrogen_count": sum(
                    partner.atomic_number == 1 for partner in atom.bond_partners
                )
            }
            for label, atom in selected.items()
        },
        "limitations": [
            "Hydrogen counts describe observed topology bonds only.",
            "Bond orders, charge method, and mechanistic state are not inferred.",
        ],
    }
