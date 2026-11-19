"""Identity-based ASH272 remapping for QM/MM selections and formal charge."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

ELEMENT_NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15}
ASP_STATE_CHARGE = {"ASP": -1, "ASH": 0}
ASH_PROTON = "HD2"


@dataclass(frozen=True)
class QMRegionAudit:
    """Audited real-atom mask, link count, formal charge, and electron parity."""

    atom_indices: tuple[int, ...]
    mask: str
    real_atom_count: int
    link_h_count: int
    formal_charge: int
    nuclear_charge_with_links: int
    electron_count: int
    multiplicity: int


def compact_mask(indices_one_based: Sequence[int]) -> str:
    """Compress sorted one-based atom indices into an Amber atom mask."""
    values = sorted(set(indices_one_based))
    if not values or values[0] < 1:
        raise ValueError("An atom mask requires positive one-based indices")
    ranges: list[str] = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return "@" + ",".join(ranges)


def compact_residue_mask(topology: Any, indices_zero_based: Sequence[int]) -> str:
    """Encode a selection as residue ranges when every selected residue is whole.

    Args:
        topology: ParmEd-like topology with residues and indexed atoms.
        indices_zero_based: Selected zero-based atom indices.

    Returns:
        An Amber residue mask with exactly the same atom membership.

    Raises:
        ValueError: If the selection is empty or contains a partial residue.
    """
    selected = set(indices_zero_based)
    if not selected:
        raise ValueError("A residue mask requires selected atoms")
    residue_numbers: list[int] = []
    covered: set[int] = set()
    for residue in topology.residues:
        atom_indices = {atom.idx for atom in residue.atoms}
        overlap = atom_indices & selected
        if overlap and overlap != atom_indices:
            raise ValueError(
                f"Cannot encode partial residue selection: {residue.idx + 1}"
            )
        if overlap:
            residue_numbers.append(residue.idx + 1)
            covered.update(overlap)
    if covered != selected:
        raise ValueError("Selection contains indices absent from the topology")
    return ":" + compact_mask(residue_numbers).removeprefix("@")


def _atom_key(atom: Any) -> tuple[int, str]:
    """Return the stable residue-position and atom-name identity."""
    return int(atom.residue.idx), str(atom.name).strip()


def _atomic_number(atom: Any) -> int:
    """Return an atom's validated nuclear charge."""
    number = int(getattr(atom, "atomic_number", 0) or 0)
    if number:
        return number
    element = str(getattr(atom, "element_name", "")).strip().upper()
    if element not in ELEMENT_NUMBERS:
        raise ValueError(f"Unsupported element for QM atom {atom.name}: {element!r}")
    return ELEMENT_NUMBERS[element]


def _bonded(topology: Any, first: Any, second: Any) -> bool:
    """Check whether two topology atoms share a bond."""
    return any(
        {bond.atom1.idx, bond.atom2.idx} == {first.idx, second.idx}
        for bond in topology.bonds
    )


def remap_indices(
    old_topology: Any,
    new_topology: Any,
    old_indices_zero_based: Sequence[int],
    protonation_states: Mapping[int, str],
    *,
    include_ash_protons: bool,
    require_complete_asp_sidechain: bool = True,
) -> tuple[int, ...]:
    """Remap atom indices by residue position/name and include added ASH protons.

    Args:
        old_topology: Topology defining the source atom identities.
        new_topology: Topology with the declared protonation state.
        old_indices_zero_based: Selected indices in ``old_topology``.
        protonation_states: One-based Amber residue position to expected ASP/ASH.
        include_ash_protons: Include HD2 for a selected ASP-to-ASH state change.

    Returns:
        Sorted zero-based atom indices in ``new_topology``.

    Raises:
        ValueError: A source atom is missing/ambiguous or protonation mismatches.
    """
    lookup: dict[tuple[int, str], list[Any]] = {}
    for atom in new_topology.atoms:
        lookup.setdefault(_atom_key(atom), []).append(atom)
    selected: set[int] = set()
    old_selected_by_residue: dict[int, set[str]] = {}
    for index in old_indices_zero_based:
        if index < 0 or index >= len(old_topology.atoms):
            raise ValueError(f"Source atom index outside topology: {index}")
        atom = old_topology.atoms[index]
        residue_index = int(atom.residue.idx)
        names = old_selected_by_residue.setdefault(residue_index, set())
        names.add(str(atom.name).strip())
        matches = lookup.get(_atom_key(atom), [])
        if len(matches) != 1:
            raise ValueError(
                f"Atom identity is missing or ambiguous: {_atom_key(atom)}"
            )
        target = matches[0]
        if _atomic_number(target) != _atomic_number(atom):
            raise ValueError(f"Element changed for atom identity: {_atom_key(atom)}")
        selected.add(int(target.idx))

    for residue_position, expected_name in protonation_states.items():
        if expected_name not in ASP_STATE_CHARGE:
            raise ValueError(f"Unsupported declared Asp state: {expected_name}")
        if residue_position < 1 or residue_position > len(new_topology.residues):
            raise ValueError(
                f"Declared residue is outside topology: {residue_position}"
            )
        old_residue = old_topology.residues[residue_position - 1]
        new_residue = new_topology.residues[residue_position - 1]
        if str(new_residue.name).strip().upper() != expected_name:
            raise ValueError(
                f"Residue {residue_position} is {new_residue.name}, "
                f"expected declared {expected_name}"
            )
        residue_index = residue_position - 1
        if residue_index not in old_selected_by_residue:
            continue
        old_name = str(old_residue.name).strip().upper()
        if old_name == expected_name:
            continue
        if (old_name, expected_name) != ("ASP", "ASH"):
            raise ValueError(
                f"Unsupported residue-state transition at {residue_position}: "
                f"{old_name} to {expected_name}"
            )
        required = {"CB", "HB2", "HB3", "CG", "OD1", "OD2"}
        if (
            require_complete_asp_sidechain
            and not required <= old_selected_by_residue[residue_index]
        ):
            raise ValueError(
                f"QM mask does not contain the full Asp sidechain at {residue_position}"
            )
        proton = [atom for atom in new_residue.atoms if atom.name == ASH_PROTON]
        oxygen = [atom for atom in new_residue.atoms if atom.name == "OD2"]
        if (
            len(proton) != 1
            or len(oxygen) != 1
            or not _bonded(new_topology, proton[0], oxygen[0])
        ):
            raise ValueError(
                f"ASH proton inventory/bond is invalid at {residue_position}"
            )
        if include_ash_protons:
            selected.add(int(proton[0].idx))

    return tuple(sorted(selected))


def audit_qm_region(
    old_topology: Any,
    new_topology: Any,
    old_indices_zero_based: Sequence[int],
    protonation_states: Mapping[int, str],
    *,
    reference_charge: int,
    expected_link_h_count: int,
    multiplicity: int = 1,
) -> QMRegionAudit:
    """Audit charge and electron parity after an explicit protonation change.

    Formal residue-charge differences are applied only when the changed Asp
    sidechain lies inside the source QM region. MM partial charges are not used.
    """
    indices = remap_indices(
        old_topology,
        new_topology,
        old_indices_zero_based,
        protonation_states,
        include_ash_protons=True,
    )
    selected_residues = {
        int(old_topology.atoms[index].residue.idx) for index in old_indices_zero_based
    }
    charge = reference_charge
    for residue_position, expected_name in protonation_states.items():
        residue_index = residue_position - 1
        if residue_index not in selected_residues:
            continue
        old_name = str(old_topology.residues[residue_index].name).strip().upper()
        if old_name != expected_name:
            charge += ASP_STATE_CHARGE[expected_name] - ASP_STATE_CHARGE[old_name]

    selected_set = set(indices)
    crossing_bonds = sum(
        (bond.atom1.idx in selected_set) != (bond.atom2.idx in selected_set)
        for bond in new_topology.bonds
    )
    if crossing_bonds != expected_link_h_count:
        raise ValueError(
            f"Expected {expected_link_h_count} QM/MM cuts, found {crossing_bonds}"
        )
    nuclear_charge = sum(_atomic_number(new_topology.atoms[index]) for index in indices)
    nuclear_charge += crossing_bonds
    electron_count = nuclear_charge - charge
    if multiplicity < 1 or (electron_count + multiplicity) % 2 != 1:
        raise ValueError("QM electron count is incompatible with multiplicity")
    return QMRegionAudit(
        atom_indices=indices,
        mask=compact_mask([index + 1 for index in indices]),
        real_atom_count=len(indices),
        link_h_count=crossing_bonds,
        formal_charge=charge,
        nuclear_charge_with_links=nuclear_charge,
        electron_count=electron_count,
        multiplicity=multiplicity,
    )


def remap_active_indices(
    old_topology: Any,
    new_topology: Any,
    old_indices_zero_based: Sequence[int],
    protonation_states: Mapping[int, str],
) -> tuple[int, ...]:
    """Remap a movable shell and include added HD2 when OD2 is mobile."""
    selected = set(
        remap_indices(
            old_topology,
            new_topology,
            old_indices_zero_based,
            protonation_states,
            include_ash_protons=False,
            require_complete_asp_sidechain=False,
        )
    )
    old_selected_by_residue: dict[int, set[str]] = {}
    for index in old_indices_zero_based:
        atom = old_topology.atoms[index]
        old_selected_by_residue.setdefault(int(atom.residue.idx), set()).add(
            str(atom.name).strip()
        )
    for residue_position, expected_name in protonation_states.items():
        residue_index = residue_position - 1
        names = old_selected_by_residue.get(residue_index, set())
        if expected_name != "ASH" or "OD2" not in names:
            continue
        residue = new_topology.residues[residue_index]
        proton = [atom for atom in residue.atoms if atom.name == ASH_PROTON]
        if len(proton) != 1:
            raise ValueError(f"Expected one HD2 at residue {residue_position}")
        selected.add(int(proton[0].idx))
    return tuple(sorted(selected))
