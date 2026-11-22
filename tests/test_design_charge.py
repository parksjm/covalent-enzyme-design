"""Formal protein charge and counterion accounting tests."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import pytest

from nyl12_covalent.design_charge import ion_counts, protein_formal_charge


N_TERMINI = {1, 233, 465, 697}
C_TERMINI = {232, 464, 696, 928, 1024, 1120, 1216, 1312}
SIDECHAIN = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}


def _protein(mutations: dict[int, str] | None = None) -> list[Any]:
    """Build 1,312 one-atom residues with exact formal residue charges."""
    mutations = mutations or {}
    residues = []
    for position in range(1, 1313):
        name = mutations.get(position, "ALA")
        charge = SIDECHAIN.get(name, 0)
        charge += int(position in N_TERMINI) - int(position in C_TERMINI)
        residues.append(
            SimpleNamespace(
                name=name,
                atoms=[SimpleNamespace(charge=float(charge))],
            )
        )
    return residues


def test_wt_and_symmetric_lys_to_asp_formal_charges() -> None:
    """Count the WT baseline and four equivalent K59D charge shifts."""
    wt = _protein(
        {
            **dict.fromkeys(range(20, 28), "ASP"),
            **dict.fromkeys((100, 332, 564, 796), "LYS"),
        }
    )
    assert protein_formal_charge(wt) == -8

    lys_to_asp = {100: "ASP", 332: "ASP", 564: "ASP", 796: "ASP"}
    k59d = _protein(
        {
            **dict.fromkeys(range(20, 28), "ASP"),
            **lys_to_asp,
        }
    )
    assert protein_formal_charge(k59d) == -16


@pytest.mark.parametrize(
    ("changes", "expected"),
    [
        ({100: "ARG", 332: "ARG", 564: "ARG", 796: "ARG"}, -8),
        ({20: "ASN", 21: "ASN", 22: "ASN", 23: "ASN"}, -4),
        ({20: "ARG", 21: "ARG", 22: "ARG", 23: "ARG"}, 0),
    ],
)
def test_symmetric_positive_or_neutral_sidechain_changes(
    changes: dict[int, str], expected: int
) -> None:
    residues = _protein(
        {
            **dict.fromkeys(range(20, 28), "ASP"),
            **dict.fromkeys((100, 332, 564, 796), "LYS"),
            **changes,
        }
    )
    assert protein_formal_charge(residues) == expected


def test_rejects_wrong_residue_charge_unsupported_microstate_and_terminal_inventory() -> (
    None
):
    residues = _protein()
    residues[0].atoms[0].charge = 0.0
    with pytest.raises(ValueError, match="residue 1 ALA charge"):
        protein_formal_charge(residues)

    for unsupported in ("HIP", "ASH", "LYN", "CYM"):
        with pytest.raises(ValueError, match="Unsupported protein microstate"):
            protein_formal_charge(_protein({10: unsupported}))

    with pytest.raises(ValueError, match="at least 1312 protein residues"):
        protein_formal_charge(_protein()[:-1])


def test_allows_explicit_neutral_histidine_tautomers() -> None:
    assert protein_formal_charge(_protein({10: "HID", 11: "HIE"})) == -4


@pytest.mark.parametrize(
    ("pairs", "charge", "expected"),
    [
        (9, -8, (17, 9)),
        (9, -16, (25, 9)),
        (9, 5, (9, 14)),
        (0, 0, (0, 0)),
    ],
)
def test_ion_counts_neutralize_formal_charge_and_add_salt_pairs(
    pairs: int, charge: int, expected: tuple[int, int]
) -> None:
    assert ion_counts(pairs, charge) == expected


@pytest.mark.parametrize("pairs", [-1, 1.5, True])
def test_ion_counts_reject_invalid_pair_count(pairs: Any) -> None:
    with pytest.raises(ValueError, match="salt_pairs"):
        ion_counts(pairs, -8)


def test_ion_counts_reject_noninteger_formal_charge() -> None:
    with pytest.raises(ValueError, match="formal_charge"):
        ion_counts(0, True)
