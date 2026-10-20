"""Audit a proposed Amber QM region without starting a QM/MM calculation."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Any, Sequence

from parmed.amber import AmberParm, AmberMask

LOGGER = logging.getLogger(__name__)
QM_MASK = ":929,1313,1314|(:930&!@C,O)"
QM_FORMAL_CHARGE = 0
_EXPECTED_RESIDUES = {929: "THR", 930: "THR", 1313: "LIG", 1314: "WAT"}
_LINK_H_ATOMIC_NUMBER = 1


def _sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atom_residue_number(atom: Any) -> int:
    """Return the one-based Amber residue number for an atom."""
    return int(atom.residue.idx) + 1


def _atom_record(atom: Any) -> dict[str, object]:
    """Serialize the relevant identity and charge fields for an atom."""
    return {
        "index": int(atom.idx) + 1,
        "residue_number": _atom_residue_number(atom),
        "residue_name": str(atom.residue.name),
        "name": str(atom.name),
        "element": str(atom.element_name),
        "atomic_number": int(atom.atomic_number),
        "partial_charge": float(atom.charge),
    }


def audit_qm_region(
    topology_path: Path,
    source_paths: Sequence[Path] = (),
    qmcharge: int = QM_FORMAL_CHARGE,
) -> dict[str, object]:
    """Audit the fixed QM selection and its sole link-H boundary.

    Args:
        topology_path: Amber prmtop whose residue and bond layout is audited.
        source_paths: Additional source files whose hashes belong in provenance.
        qmcharge: Requested integer formal charge for the QM region.

    Returns:
        A structural selection report with atom identities, charges, hashes,
        link boundary, and spin-parity check.

    Raises:
        ValueError: If residue layout, hydrogens, boundary, or electron parity
            does not match the requested singlet region.
        FileNotFoundError: If any topology or source file is missing.
    """
    topology_path = Path(topology_path).resolve(strict=True)
    sources = [Path(path).resolve(strict=True) for path in source_paths]
    if type(qmcharge) is not int or qmcharge != QM_FORMAL_CHARGE:
        raise ValueError(
            f"This audited TI1 state requires integer qmcharge={QM_FORMAL_CHARGE}"
        )
    topology = AmberParm(str(topology_path))
    atoms_by_residue: dict[int, list[Any]] = {}
    for atom in topology.atoms:
        atoms_by_residue.setdefault(_atom_residue_number(atom), []).append(atom)

    for number, expected_name in _EXPECTED_RESIDUES.items():
        residue_atoms = atoms_by_residue.get(number, [])
        if not residue_atoms:
            raise ValueError(f"Required residue {number} {expected_name} is missing")
        actual_names = {str(atom.residue.name) for atom in residue_atoms}
        if actual_names != {expected_name}:
            raise ValueError(
                f"Residue {number} must be {expected_name}; found {sorted(actual_names)}"
            )

    selected = {
        atom
        for atom in topology.atoms
        if _atom_residue_number(atom) in {929, 1313, 1314}
        or (_atom_residue_number(atom) == 930 and str(atom.name) not in {"C", "O"})
    }
    if not selected:
        raise ValueError("QM atom selection is empty")
    mask_selected = set(AmberMask(topology, QM_MASK).Selected())
    explicit_selected = {int(atom.idx) for atom in selected}
    if mask_selected != explicit_selected:
        missing = sorted(explicit_selected - mask_selected)
        extra = sorted(mask_selected - explicit_selected)
        raise ValueError(
            "Explicit atom selection disagrees with AmberMask "
            f"(missing from mask={missing}, extra in mask={extra})"
        )

    cut_bonds: list[tuple[Any, Any]] = []
    for atom in selected:
        for partner in atom.bond_partners:
            if partner not in selected:
                if atom.idx < partner.idx:
                    cut_bonds.append((atom, partner))
                elif partner.idx < atom.idx:
                    cut_bonds.append((partner, atom))
        if int(atom.atomic_number) == 1:
            unselected_partners = [p for p in atom.bond_partners if p not in selected]
            if unselected_partners:
                raise ValueError(
                    f"Selected hydrogen {atom.residue.name}{_atom_residue_number(atom)}:{atom.name} "
                    "is bonded to an unselected atom"
                )
    unique_cut_bonds = {
        (int(first.idx), int(second.idx)): (first, second)
        for first, second in cut_bonds
    }
    cut_bonds = list(unique_cut_bonds.values())
    if len(cut_bonds) != 1:
        raise ValueError(
            f"Expected exactly one QM/MM boundary bond; found {len(cut_bonds)}"
        )
    first, second = cut_bonds[0]
    first_record, second_record = _atom_record(first), _atom_record(second)
    if {first_record["element"], second_record["element"]} != {"C"}:
        raise ValueError("The sole QM/MM cut must be a carbon-carbon bond")
    first_selected = first in selected
    second_selected = second in selected
    if first_selected == second_selected:
        raise ValueError(
            "QM/MM boundary must have one selected and one unselected endpoint"
        )
    qm_link, mm_boundary = (first, second) if first_selected else (second, first)
    if _atom_residue_number(qm_link) != 930 or _atom_residue_number(mm_boundary) != 930:
        raise ValueError("The QM/MM cut must originate at residue 930")
    if str(qm_link.name) != "CA" or str(mm_boundary.name) != "C":
        raise ValueError("Expected the QM CA930-MM C930 carbon-carbon boundary")

    selected_heavy = {atom for atom in selected if int(atom.atomic_number) > 1}
    for atom in selected_heavy:
        for partner in atom.bond_partners:
            if int(partner.atomic_number) == 1 and partner not in selected:
                raise ValueError(
                    f"Hydrogen {partner.name} bonded to selected heavy atom "
                    f"{atom.residue.name}{_atom_residue_number(atom)}:{atom.name} is unselected"
                )

    records = sorted(
        (_atom_record(atom) for atom in selected), key=lambda row: int(row["index"])
    )
    partial_charge = sum(float(atom.charge) for atom in selected)
    nuclear_charge = sum(int(atom.atomic_number) for atom in selected)
    electron_count = nuclear_charge + _LINK_H_ATOMIC_NUMBER - qmcharge
    if electron_count % 2:
        raise ValueError(
            f"QM region has odd electron count {electron_count} after link-H addition; "
            "a singlet is not allowed"
        )

    report: dict[str, object] = {
        "audit_type": "structural_qm_region_selection_only",
        "topology": str(topology_path),
        "topology_sha256": _sha256(topology_path),
        "sources": [{"path": str(path), "sha256": _sha256(path)} for path in sources],
        "qmmask": QM_MASK,
        "qmcharge": qmcharge,
        "spin_assumption": "singlet",
        "selected_atom_count": len(selected),
        "selected_atoms": records,
        "selected_atom_indices_amber_1based": [int(row["index"]) for row in records],
        "selected_partial_charge": partial_charge,
        "qm_formal_charge_is_independent_of_partial_charge": True,
        "selected_nuclear_charge": nuclear_charge,
        "link_hydrogen_nuclear_charge": _LINK_H_ATOMIC_NUMBER,
        "electron_count_after_link_hydrogen": electron_count,
        "electron_count_parity": "even",
        "link_boundary": {
            "qm_atom": _atom_record(qm_link),
            "mm_atom": _atom_record(mm_boundary),
            "bond_elements": ["C", "C"],
            "link_atom": "H",
            "link_hydrogen": "H",
            "link_hydrogen_adds_one_proton_and_one_electron": True,
        },
        "limitations": [
            "This validates atom selection and boundary structure only.",
            "It does not validate DFTB, GB, or QM/MM electrostatic embedding behavior.",
            "It does not validate force-field charges or topology energetics.",
        ],
    }
    LOGGER.info(
        "Audited %d QM atoms; partial charge %.6f; formal charge %d; electrons %d",
        len(selected),
        partial_charge,
        qmcharge,
        electron_count,
    )
    return report
