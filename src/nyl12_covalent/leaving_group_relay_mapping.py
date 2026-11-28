"""Explicit topology-owned atom and proton inventories for relay screening."""

from __future__ import annotations

from typing import Any


def resolve_atom(topology: Any, descriptor: dict) -> int:
    """Resolve an atom from verified one-based residue position and name.

    Args:
        topology: ParmEd topology.
        descriptor: residue_index1, resname, atom, optional chain.

    Returns:
        Zero-based atom index.
    """
    index = int(descriptor["residue_index1"]) - 1
    if index < 0 or index >= len(topology.residues):
        raise ValueError("Residue index outside topology")
    residue = topology.residues[index]
    if residue.name != descriptor["resname"]:
        raise ValueError(f"Residue identity mismatch at {index + 1}: {residue.name}")
    if "chain" in descriptor and residue.chain != descriptor["chain"]:
        raise ValueError("Chain identity mismatch")
    atoms = [a for a in residue.atoms if a.name == descriptor["atom"]]
    if len(atoms) != 1:
        raise ValueError(f"Atom mapping ambiguous or absent: {descriptor}")
    return atoms[0].idx


def resolve_roles(topology: Any, config: dict) -> tuple[list[dict], dict[str, int]]:
    """Validate explicit donor/acceptor states, bonded protons and water roles.

    Args:
        topology: Independent topology for this system.
        config: Case with roles and geometry_atoms mappings.

    Returns:
        Chemical role records and geometry atom indices.
    """
    roles = []
    names = set()
    for role in config["roles"]:
        label, state = role["label"], role["state"]
        if not state or label in names:
            raise ValueError(
                "Role requires unique label and explicit protonation state"
            )
        names.add(label)
        index = resolve_atom(topology, role)
        atom = topology.atoms[index]
        owned = {a.name: a.idx for a in atom.bond_partners if a.atomic_number == 1}
        expected = role["hydrogens"]
        if set(expected) != set(owned) or len(expected) != len(owned):
            raise ValueError(f"Hydrogen ownership/inventory mismatch: {label}")
        if role.get("start", False) and not owned:
            raise ValueError(f"Starting donor lacks explicit proton: {label}")
        if state == "NH3+" and (len(owned) != 3 or role["acceptor"]):
            raise ValueError("NH3+ must have three hydrogens and cannot accept")
        if state == "NH2" and len(owned) != 2:
            raise ValueError("NH2 requires two hydrogens")
        if state in {"COO-", "deprotonated"} and owned:
            raise ValueError("Deprotonated atom cannot own a donor hydrogen")
        roles.append(
            {
                **role,
                "index": index,
                "h_indices": list(owned.values()),
                "water": False,
                "conditional": role["conditional"],
            }
        )
    for residue in topology.residues:
        if residue.name not in config.get("water_resnames", ["WAT"]):
            continue
        oxygen = [a for a in residue.atoms if a.atomic_number == 8]
        if len(oxygen) != 1 or len(residue.atoms) != 3:
            raise ValueError("Expected explicit three-atom water")
        hydrogens = [a.idx for a in oxygen[0].bond_partners if a.atomic_number == 1]
        if len(hydrogens) != 2:
            raise ValueError("Water hydrogen ownership mismatch")
        roles.append(
            {
                "label": f"water:{residue.idx + 1}",
                "index": oxygen[0].idx,
                "h_indices": hydrogens,
                "water": True,
                "acceptor": True,
                "start": False,
                "internal": True,
                "state": "H2O",
                "conditional": False,
            }
        )
    return roles, {
        key: resolve_atom(topology, value)
        for key, value in config["geometry_atoms"].items()
    }
