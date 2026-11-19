"""Independent atom, proton and boundary validation for declared PA6 ES cores."""

from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _identity(atom: Any) -> dict[str, Any]:
    return {
        "index": atom.idx + 1,
        "residue": atom.residue.idx + 1,
        "resname": atom.residue.name,
        "name": atom.name,
        "element": int(atom.atomic_number),
        "mass": float(atom.mass),
        "type": atom.type,
        "charge": float(atom.charge),
    }


def _selected_atoms(parm: Any, groups: list[dict]) -> tuple[set[int], dict[int, str]]:
    selected, owners = set(), {}
    names = set()
    for group in groups:
        if (
            not isinstance(group["name"], str)
            or not group["name"]
            or group["name"] in names
        ):
            raise ValueError("Chemical group names must be unique nonempty strings")
        names.add(group["name"])
        if type(group["formal_charge"]) is not int or not group["atoms"]:
            raise ValueError(
                "Each nonempty chemical group needs an integer formal charge"
            )
        for row in group["atoms"]:
            number = row["index"]
            if (
                type(number) is not int
                or not 1 <= number <= len(parm.atoms)
                or number in selected
            ):
                raise ValueError("QM atom indices must be unique one-based integers")
            actual = _identity(parm.atoms[number - 1])
            for key in ["index", "residue", "resname", "name", "element", "type"]:
                if row[key] != actual[key]:
                    raise ValueError(f"QM atom identity mismatch: {number} {key}")
            for key in ["mass", "charge"]:
                if (
                    not np.isfinite(row[key])
                    or abs(float(row[key]) - actual[key]) > 1e-8
                ):
                    raise ValueError(
                        f"QM atom numeric identity mismatch: {number} {key}"
                    )
            selected.add(number)
            owners[number] = group["name"]
    return selected, owners


def _components(parm: Any, selected: set[int], owners: dict[int, str]) -> None:
    for number in selected:
        atom = parm.atoms[number - 1]
        for partner in atom.bond_partners:
            other = partner.idx + 1
            if (
                atom.atomic_number == 1 or partner.atomic_number == 1
            ) and other not in selected:
                raise ValueError(
                    "QM selection must include every bonded hydrogen and its parent"
                )
            if other in selected and owners[number] != owners[other]:
                raise ValueError(
                    "Declared ES chemical components must not be covalently joined"
                )
    for group in set(owners.values()):
        members = {number for number, owner in owners.items() if owner == group}
        visited, pending = set(), [min(members)]
        while pending:
            number = pending.pop()
            if number in visited:
                continue
            visited.add(number)
            pending.extend(
                a.idx + 1
                for a in parm.atoms[number - 1].bond_partners
                if a.idx + 1 in members and a.idx + 1 not in visited
            )
        if visited != members:
            raise ValueError(f"Disconnected declared chemical component: {group}")


def _roles(
    parm: Any,
    selected: set[int],
    roles: dict[str, int],
    role_identities: dict[str, dict[str, Any]] | None = None,
) -> dict:
    defaults = {
        "thr_n_terminus": ("N", 7, 2),
        "nucleophile": ("OG1", 8, 1),
        "oxyanion_donor_backbone": ("N", 7, 1),
        "oxyanion_donor_sidechain": ("ND2", 7, 2),
        "asp272_od1": ("OD1", 8, 0),
        "asp272_od2": ("OD2", 8, 0),
        "asp274_od1": ("OD1", 8, 0),
        "asp274_od2": ("OD2", 8, 0),
        "scissile_carbon": ("C54", 6, 0),
        "scissile_oxygen": ("O49", 8, 0),
        "leaving_nitrogen": ("N59", 7, 1),
        "donor_carbonyl_c": ("C", 6, 0),
        "donor_carbonyl_o": ("O", 8, 0),
    }
    if role_identities is None:
        required = defaults
    else:
        if set(role_identities) != set(defaults):
            raise ValueError("Explicit role identities must cover every ES role")
        required = {}
        for role, identity in role_identities.items():
            if set(identity) != {"name", "element", "bonded_hydrogens"}:
                raise ValueError(f"Invalid explicit role identity schema: {role}")
            name = identity["name"]
            element = identity["element"]
            hydrogens = identity["bonded_hydrogens"]
            if (
                not isinstance(name, str)
                or not name
                or type(element) is not int
                or type(hydrogens) is not int
                or hydrogens < 0
            ):
                raise ValueError(f"Invalid explicit role identity: {role}")
            required[role] = (name, element, hydrogens)
    if set(roles) != set(required) or len(set(roles.values())) != len(required):
        raise ValueError("Explicit distinct ES role indices are required")
    observations = {}
    for role, (name, element, count) in required.items():
        number = roles[role]
        if type(number) is not int or number not in selected:
            raise ValueError(f"Role outside explicit QM selection: {role}")
        atom = parm.atoms[number - 1]
        hydrogens = [a.idx + 1 for a in atom.bond_partners if a.atomic_number == 1]
        if (
            atom.name != name
            or atom.atomic_number != element
            or len(hydrogens) != count
        ):
            raise ValueError(f"ES role/proton inventory mismatch: {role}")
        observations[role] = {
            "atom_index": number,
            "hydrogen_indices": sorted(hydrogens),
        }
    for left, right in [
        ("scissile_carbon", "scissile_oxygen"),
        ("scissile_carbon", "leaving_nitrogen"),
        ("donor_carbonyl_c", "donor_carbonyl_o"),
        ("donor_carbonyl_c", "oxyanion_donor_backbone"),
    ]:
        if (
            parm.atoms[roles[right] - 1]
            not in parm.atoms[roles[left] - 1].bond_partners
        ):
            raise ValueError("Scissile and donor amide bonds must remain intact")
    if (
        parm.atoms[roles["scissile_carbon"] - 1]
        in parm.atoms[roles["nucleophile"] - 1].bond_partners
    ):
        raise ValueError("Declared ES must not contain a Thr-substrate covalent bond")
    return observations


def audit_es_qm_region(
    topology: Path, coordinates: Path, spec: dict[str, Any]
) -> dict[str, Any]:
    """Independently validate a source-pinned explicit ES atom and cap proposal.

    Args:
        topology: Original full-solvent Amber topology; never modified.
        coordinates: Exact source-pinned imaged Amber coordinates.
        spec: Explicit atom identities, six chemical groups, role indices,
            formal charge, multiplicity, expected boundaries and source hashes.

    Returns:
        Atom mask, proton and boundary inventories, formal/electronic accounting
        and original partial-charge diagnostics. This does not validate embedding.

    Raises:
        FileNotFoundError: An input is missing.
        ValueError: Source identity, chemistry or atom/boundary accounting fails.
    """
    for path, key in [
        (topology, "topology_sha256"),
        (coordinates, "coordinates_sha256"),
    ]:
        if not path.is_file():
            raise FileNotFoundError(path)
        if _sha256(path) != spec[key]:
            raise ValueError(f"Source hash mismatch: {path}")
    if (
        type(spec["qmcharge"]) is not int
        or spec["qmcharge"] != -2
        or spec["multiplicity"] != 1
    ):
        raise ValueError(
            "This explicit ES case requires formal charge-2 and singlet multiplicity"
        )
    if (
        type(spec["multiplicity"]) is not int
        or spec["allow_multiple_link_parents"] is not False
    ):
        raise ValueError("Singlet and unique link parents must be explicitly declared")
    groups = spec["groups"]
    if len(groups) != 6 or sorted(group["formal_charge"] for group in groups) != [
        -1,
        -1,
        0,
        0,
        0,
        0,
    ]:
        raise ValueError("Expected four neutral components and two carboxylates")
    if sum(group["formal_charge"] for group in groups) != spec["qmcharge"]:
        raise ValueError("Formal group charges must sum to qmcharge")
    parm = AmberParm(str(topology), xyz=str(coordinates))
    xyz = np.asarray(parm.coordinates, dtype=float)
    if xyz.shape != (len(parm.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError(
            "Source coordinates must be finite and match all topology atoms"
        )
    selected, owners = _selected_atoms(parm, groups)
    _components(parm, selected, owners)
    roles = _roles(
        parm,
        selected,
        spec["roles"],
        spec.get("role_identities"),
    )
    masses = {1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999}
    for number in selected:
        atom = parm.atoms[number - 1]
        if (
            atom.atomic_number not in masses
            or abs(atom.mass - masses[atom.atomic_number]) > 0.05
            or not np.isfinite(atom.charge)
        ):
            raise ValueError(
                "Unsupported QM element, inconsistent mass or nonfinite charge"
            )
    actual_boundaries = []
    for bond in parm.bonds:
        left, right = bond.atom1, bond.atom2
        if (left.idx + 1 in selected) == (right.idx + 1 in selected):
            continue
        qm, mm = (left, right) if left.idx + 1 in selected else (right, left)
        if (
            qm.atomic_number != 6
            or qm.type not in {"CT", "CX", "XC", "2C", "c3"}
            or mm.atomic_number not in {6, 7}
            or len(qm.bond_partners) != 4
        ):
            raise ValueError("Link cuts require declared saturated carbon QM endpoints")
        actual_boundaries.append(
            {
                "qm_index": qm.idx + 1,
                "mm_index": mm.idx + 1,
                "qm_atom": _identity(qm),
                "mm_atom": _identity(mm),
                "source_bond_length_A": float(
                    np.linalg.norm(xyz[qm.idx] - xyz[mm.idx])
                ),
                "mm_bonded_neighbors": [_identity(a) for a in mm.bond_partners],
            }
        )
    expected = [(row["qm_index"], row["mm_index"]) for row in spec["boundaries"]]
    actual = [(row["qm_index"], row["mm_index"]) for row in actual_boundaries]
    if len(set(expected)) != len(expected) or set(expected) != set(actual):
        raise ValueError("Declared boundary pairs do not equal actual crossing bonds")
    if any(count != 1 for count in Counter(q for q, _ in actual).values()):
        raise ValueError("Multiple links at one QM atom are forbidden for this pilot")
    if any(count != 1 for count in Counter(m for _, m in actual).values()):
        raise ValueError("One MM boundary atom cannot be reused for multiple links")
    electrons = (
        sum(parm.atoms[number - 1].atomic_number for number in selected)
        + len(actual)
        - spec["qmcharge"]
    )
    if electrons % 2:
        raise ValueError(
            "Declared singlet requires an even electron count including links"
        )
    source_charge = float(sum(atom.charge for atom in parm.atoms))
    if not np.isfinite(source_charge):
        raise ValueError("Nonfinite full-system source charge")
    return {
        "name": spec["name"],
        "state": "ES; not optimized TI/TS",
        "qmmask": "@" + ",".join(map(str, sorted(selected))),
        "qmcharge": spec["qmcharge"],
        "multiplicity": 1,
        "real_atom_count": len(selected),
        "link_H_count": len(actual),
        "total_QM_nuclei": len(selected) + len(actual),
        "electron_count_with_links": electrons,
        "qm_atoms": [_identity(parm.atoms[number - 1]) for number in sorted(selected)],
        "boundaries": actual_boundaries,
        "proton_inventory": roles,
        "formal_charge_accounting": {
            group["name"]: group["formal_charge"] for group in groups
        },
        "source_QM_partial_charge_sum": float(
            sum(parm.atoms[number - 1].charge for number in selected)
        ),
        "source_total_partial_charge": source_charge,
        "source_charge_integer_residual": source_charge - round(source_charge),
        "embedding_and_link_distance_validated": False,
    }
