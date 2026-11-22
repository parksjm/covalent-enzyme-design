"""Formal protein charge and salt counterion accounting for design models."""

from __future__ import annotations

from typing import Any

PROTEIN_RESIDUE_COUNT = 1312
N_TERMINI = frozenset((1, 233, 465, 697))
C_TERMINI = frozenset((232, 464, 696, 928, 1024, 1120, 1216, 1312))
CHARGE_TOLERANCE = 1.0e-4
STANDARD_RESIDUES = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)
ALLOWED_RESIDUES = STANDARD_RESIDUES | {"HID", "HIE"}
SIDECHAIN_CHARGES = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}


def protein_formal_charge(residues: list[Any]) -> int:
    """Validate and sum formal charges for the first 1,312 protein residues.

    Args:
        residues: ParmEd-like residues in topology order, including all protein
            chains before solvent and ions.

    Returns:
        Integer formal charge after side-chain and specified terminal states.

    Raises:
        ValueError: If the protein inventory, microstate, or residue charge is
            unsupported or inconsistent with its formal charge.
    """
    if len(residues) < PROTEIN_RESIDUE_COUNT:
        raise ValueError(
            f"Expected at least {PROTEIN_RESIDUE_COUNT} protein residues, "
            f"found {len(residues)}"
        )
    total = 0
    for position, residue in enumerate(residues[:PROTEIN_RESIDUE_COUNT], start=1):
        name = residue.name.strip().upper()
        if name not in ALLOWED_RESIDUES:
            raise ValueError(f"Unsupported protein microstate {name!r} at {position}")
        expected = SIDECHAIN_CHARGES.get(name, 0)
        expected += int(position in N_TERMINI) - int(position in C_TERMINI)
        observed = float(sum(atom.charge for atom in residue.atoms))
        if abs(observed - expected) > CHARGE_TOLERANCE:
            raise ValueError(
                f"Protein residue {position} {name} charge {observed:.8f}; "
                f"expected {expected:+d} within {CHARGE_TOLERANCE:g}"
            )
        total += expected
    return total


def ion_counts(salt_pairs: int, formal_charge: int) -> tuple[int, int]:
    """Return sodium/chloride counts for neutralization plus salt pairs.

    Args:
        salt_pairs: Number of extra sodium/chloride pairs for salt concentration.
        formal_charge: Integer formal charge of the solute before ions.

    Returns:
        Sodium and chloride residue counts, respectively.

    Raises:
        ValueError: If either argument is not an integer or salt_pairs is negative.
    """
    if (
        isinstance(salt_pairs, bool)
        or not isinstance(salt_pairs, int)
        or salt_pairs < 0
    ):
        raise ValueError("salt_pairs must be a nonnegative integer")
    if isinstance(formal_charge, bool) or not isinstance(formal_charge, int):
        raise ValueError("formal_charge must be an integer")
    return salt_pairs + max(-formal_charge, 0), salt_pairs + max(formal_charge, 0)
