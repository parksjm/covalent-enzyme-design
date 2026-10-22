"""Build an explicit atom crosswalk between the TI1 model and native graph."""

from __future__ import annotations

from typing import Any

from rdkit import Chem


_SAME_ORIGIN_IDS = (
    "THR:N",
    "THR:CA",
    "THR:CB",
    "THR:OG1",
    "THR:CG2",
    "THR:C",
    "THR:O",
    "NME:N",
    "NME:CH3",
)
_REACTIVE_ORIGINS = {
    "TIL:C*": "PA6:C54",
    "TIL:O": "PA6:O49",
    "TIL:N": "PA6:N59",
    "TIL:acylCH3": "PA6:C62",
}
_REQUIRED_ORIGINS = (*_SAME_ORIGIN_IDS, *_REACTIVE_ORIGINS)
_REACTIVE_BONDS = (
    ("THR:OG1", "TIL:C*"),
    ("TIL:C*", "TIL:O"),
    ("TIL:C*", "TIL:N"),
    ("TIL:C*", "TIL:acylCH3"),
    ("TIL:N", "TIL:leavingCH3"),
)


def _origin_index(molecule: Chem.Mol, label: str) -> dict[str, int]:
    """Return unique atom-origin indices, requiring an origin on every atom."""
    result: dict[str, int] = {}
    for atom in molecule.GetAtoms():
        if not atom.HasProp("origin_id"):
            raise ValueError(f"Atom index {atom.GetIdx()} lacks origin_id")
        origin = atom.GetProp("origin_id")
        if origin in result:
            raise ValueError(f"Duplicate {label} origin_id: {origin}")
        result[origin] = atom.GetIdx()
    return result


def _formal_charge(molecule: Chem.Mol) -> int:
    return sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())


def _bond_is_single(molecule: Chem.Mol, first: int, second: int) -> bool:
    bond = molecule.GetBondBetweenAtoms(first, second)
    return bond is not None and bond.GetBondType() == Chem.BondType.SINGLE


def _heavy_origins(molecule: Chem.Mol) -> list[str]:
    return sorted(
        atom.GetProp("origin_id")
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1
    )


def _graph_summary(molecule: Chem.Mol) -> dict[str, int]:
    return {
        "atom_count": molecule.GetNumAtoms(),
        "heavy_atom_count": sum(
            atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms()
        ),
        "formal_charge": _formal_charge(molecule),
    }


def build_ti1_model_crosswalk(surrogate: Chem.Mol, native: Chem.Mol) -> dict[str, Any]:
    """Validate TI1 atom correspondence without transferring coordinates/charges.

    Args:
        surrogate: The 37-atom capped anionic TI1 model.
        native: The 103-atom native PA6L4/threonine/NME/water graph.

    Returns:
        A mapping report with atom correspondences, charge and graph summaries,
        and explicit limitations on charge and coordinate transfer.

    Raises:
        ValueError: If identifiers, required chemical states, connectivity,
            elements, or mapped stereochemistry do not satisfy the contract.
    """
    if surrogate.GetNumAtoms() != 37 or native.GetNumAtoms() != 103:
        raise ValueError("TI1 crosswalk requires 37-atom and 103-atom graphs")
    surrogate_ids = _origin_index(surrogate, "surrogate")
    native_ids = _origin_index(native, "native")

    missing_surrogate = (set(_REQUIRED_ORIGINS) | {"TIL:leavingCH3"}) - set(
        surrogate_ids
    )
    if missing_surrogate:
        raise ValueError(f"Missing surrogate origins: {sorted(missing_surrogate)}")
    missing_native = set(_SAME_ORIGIN_IDS) - set(native_ids)
    missing_native.update(
        target for target in _REACTIVE_ORIGINS.values() if target not in native_ids
    )
    if missing_native:
        raise ValueError(f"Missing native origins: {sorted(missing_native)}")

    native_leaving_candidates = []
    native_n_index = native_ids["PA6:N59"]
    for neighbor in native.GetAtomWithIdx(native_n_index).GetNeighbors():
        if neighbor.GetAtomicNum() > 1 and neighbor.GetProp("origin_id") != "PA6:C54":
            if neighbor.GetAtomicNum() != 6:
                raise ValueError("PA6:N59 non-C54 heavy neighbor must be carbon")
            native_leaving_candidates.append(neighbor.GetIdx())
    if len(native_leaving_candidates) != 1:
        raise ValueError("PA6:N59 must have exactly one heavy neighbor besides PA6:C54")
    leaving_idx = native_leaving_candidates[0]
    leaving_origin = native.GetAtomWithIdx(leaving_idx).GetProp("origin_id")
    if leaving_origin in native_ids and native_ids[leaving_origin] != leaving_idx:
        raise ValueError(f"Ambiguous native leaving-side origin: {leaving_origin}")

    mapping = {origin: origin for origin in _SAME_ORIGIN_IDS}
    mapping.update(_REACTIVE_ORIGINS)
    mapping["TIL:leavingCH3"] = leaving_origin
    rows: list[dict[str, Any]] = []
    for source, target in mapping.items():
        source_atom = surrogate.GetAtomWithIdx(surrogate_ids[source])
        target_atom = native.GetAtomWithIdx(native_ids[target])
        if source_atom.GetAtomicNum() != target_atom.GetAtomicNum():
            raise ValueError(f"Element mismatch for {source} -> {target}")
        rows.append(
            {
                "surrogate_origin_id": source,
                "native_origin_id": target,
                "surrogate_index_zero_based": source_atom.GetIdx(),
                "native_index_zero_based": target_atom.GetIdx(),
                "surrogate_index_one_based": source_atom.GetIdx() + 1,
                "native_index_one_based": target_atom.GetIdx() + 1,
                "element": source_atom.GetSymbol(),
                "surrogate_formal_charge": source_atom.GetFormalCharge(),
                "native_formal_charge": target_atom.GetFormalCharge(),
                "correspondence": "origin/connectivity correspondence only",
                "native_charge_transfer_eligible": False,
                "forcefield_validated": False,
            }
        )

    if _formal_charge(surrogate) != -1 or _formal_charge(native) != 0:
        raise ValueError("Unexpected total formal charge for TI1 model state")
    for molecule, ids, label, expected_charge, expected_h in (
        (surrogate, surrogate_ids, "surrogate", 0, 1),
        (native, native_ids, "native", 1, 3),
    ):
        thr_n = molecule.GetAtomWithIdx(ids["THR:N"])
        explicit_h = sum(n.GetAtomicNum() == 1 for n in thr_n.GetNeighbors())
        if thr_n.GetFormalCharge() != expected_charge or explicit_h != expected_h:
            raise ValueError(f"Unexpected {label} threonine N protonation state")
    for molecule, ids, label in (
        (surrogate, surrogate_ids, "surrogate"),
        (native, native_ids, "native"),
    ):
        oxygen_id = "TIL:O" if label == "surrogate" else "PA6:O49"
        if molecule.GetAtomWithIdx(ids[oxygen_id]).GetFormalCharge() != -1:
            raise ValueError(f"Unexpected {label} oxyanion charge")
        og1_id = "THR:OG1"
        og1 = molecule.GetAtomWithIdx(ids[og1_id])
        if any(neighbor.GetAtomicNum() == 1 for neighbor in og1.GetNeighbors()):
            raise ValueError(f"{label} THR:OG1 must have no explicit hydrogen")

    for first, second in _REACTIVE_BONDS:
        first_surrogate = surrogate_ids[first]
        second_surrogate = surrogate_ids[second]
        native_first = native_ids[mapping[first]]
        native_second = native_ids[mapping[second]]
        if not _bond_is_single(surrogate, first_surrogate, second_surrogate):
            raise ValueError(f"Surrogate bond {first}-{second} must be single")
        if not _bond_is_single(native, native_first, native_second):
            raise ValueError(f"Native mapped bond {first}-{second} must be single")

    surrogate_c = surrogate.GetAtomWithIdx(surrogate_ids["TIL:C*"])
    native_c = native.GetAtomWithIdx(native_ids["PA6:C54"])
    surrogate_cip = (
        surrogate_c.GetProp("_CIPCode") if surrogate_c.HasProp("_CIPCode") else None
    )
    native_cip = native_c.GetProp("_CIPCode") if native_c.HasProp("_CIPCode") else None
    if surrogate_cip not in {"R", "S"} or native_cip != surrogate_cip:
        raise ValueError("Mapped C* stereochemistry must have matching R/S CIP codes")

    for origin in ("TIL:acylCH3", "TIL:leavingCH3"):
        source_atom = surrogate.GetAtomWithIdx(surrogate_ids[origin])
        target_atom = native.GetAtomWithIdx(native_ids[mapping[origin]])
        source_h = sum(n.GetAtomicNum() == 1 for n in source_atom.GetNeighbors())
        target_h = sum(n.GetAtomicNum() == 1 for n in target_atom.GetNeighbors())
        if (source_h, target_h) != (3, 2):
            raise ValueError(f"Unexpected methyl-to-methylene boundary at {origin}")

    water_ids = {"WAT:O", "WAT:H1", "WAT:H2"}
    if not water_ids.issubset(native_ids) or water_ids.intersection(surrogate_ids):
        raise ValueError("Native-only WAT:O/H1/H2 must exist only in the native graph")
    return {
        "correspondences": rows,
        "unmapped_surrogate_heavy_atom_origins": sorted(
            set(_heavy_origins(surrogate)) - set(mapping)
        ),
        "unmapped_native_heavy_atom_origins": sorted(
            set(_heavy_origins(native)) - set(mapping.values())
        ),
        "omitted_hydrogen_policy": (
            "No hydrogen mapping or equivalence is asserted; protonation and "
            "methyl-to-methylene boundaries differ."
        ),
        "state_differences": [
            "THR:N protonation changes from neutral NH to cationic NH3+.",
            "Acyl and leaving methyl groups map to native methylene groups.",
            "Native-only water is WAT:O, WAT:H1, and WAT:H2.",
        ],
        "native_charge_transfer_eligible": False,
        "forcefield_validated": False,
        "native_coordinates_generated": False,
        "surrogate_graph": _graph_summary(surrogate),
        "native_graph": _graph_summary(native),
    }
