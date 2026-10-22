"""Tests for the independent AF3 water/contact geometry audit."""

from __future__ import annotations

from dataclasses import replace
import math
from typing import Callable

import pytest

from nyl12_covalent.af3_postprocess import AtomSite
from nyl12_covalent.af3_water_audit import audit_sites


def _site(
    serial: int,
    element: str,
    name: str,
    residue: str,
    chain: str,
    sequence: int,
    xyz: tuple[float, float, float],
    confidence: float = 85.0,
) -> AtomSite:
    """Build one immutable heavy-atom record."""
    return AtomSite(serial, element, name, residue, chain, sequence, xyz, confidence)


def _valid_sites() -> list[AtomSite]:
    """Create the canonical A/E/I/J heavy-atom inventory with valid geometry."""
    sites: list[AtomSite] = []
    serial = 1
    for sequence in range(1, 233):
        residue = "PHE" if sequence == 113 else "ASN" if sequence == 186 else "ALA"
        sites.append(
            _site(serial, "C", "CA", residue, "A", sequence, (500 + sequence, 100, 0))
        )
        serial += 1
        if sequence == 113:
            sites.append(_site(serial, "N", "N", residue, "A", sequence, (10, 10, 0)))
            serial += 1
        if sequence == 186:
            sites.append(_site(serial, "N", "ND2", residue, "A", sequence, (20, 10, 0)))
            serial += 1
    for sequence in range(1, 97):
        residue = "THR" if sequence == 1 else "ALA"
        sites.append(
            _site(serial, "C", "CA", residue, "E", sequence, (700 + sequence, 0, 0))
        )
        serial += 1
        if sequence == 1:
            sites.extend(
                [
                    _site(serial, "N", "N", "THR", "E", 1, (2.8, 0, 0)),
                    _site(serial + 1, "O", "OG1", "THR", "E", 1, (0, 2.8, 0)),
                ]
            )
            serial += 2

    ligand_atoms = [
        ("C54", "C", (0, 4.4, 0)),
        ("O49", "O", (1.5, 4.4, 0)),
        ("N59", "N", (-1.5, 4.4, 0.7)),
        ("C62", "C", (0, 4.4, 1.5)),
    ]
    ligand_atoms.extend(
        (f"X{index:02d}", "C", (900 + index, 0, 0)) for index in range(29)
    )
    for name, element, xyz in ligand_atoms:
        sites.append(_site(serial, element, name, "PA6TI1", "I", 1, xyz))
        serial += 1
    sites.append(_site(serial, "O", "O", "HOH", "J", 1, (0, 0, 0)))
    return sites


def test_valid_model_reports_geometry_contacts_and_confidence() -> None:
    """A valid screened model reports counts and retains the QMMM rejection."""
    report = audit_sites(_valid_sites())

    assert report["inventory_passed"] is True
    assert report["chain_heavy_atom_counts"] == {"A": 234, "E": 98, "I": 33, "J": 1}
    assert report["bridge_distances"] == {"N_water_O_A": 2.8, "OG1_water_O_A": 2.8}
    assert report["bridge_distance_screen_passed"] is True
    assert report["water_heavy_contact_screen_passed"] is True
    assert report["reactive_bond_screen_passed"] is True
    assert report["heavy_atom_screen_passed"] is True
    assert report["accepted_for_qmmm"] is False
    assert report["confidence_from_rounded_cif"]["A"] == {
        "count": 234,
        "mean_plddt": 85.0,
        "minimum_plddt": 85.0,
        "maximum_plddt": 85.0,
    }
    assert len(report["reactive_angles"]) == 6
    assert "active_site" in report
    assert report["limitations"]


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda sites: sites + [sites[0]], "Duplicate atom identity"),
        (
            lambda sites: [
                replace(site, chain="K") if site.chain == "J" else site
                for site in sites
            ],
            "Expected exactly chains A, E, I, J",
        ),
        (
            lambda sites: [
                site for site in sites if not (site.chain == "A" and site.sequence == 1)
            ],
            "Unexpected residue inventory for chain A",
        ),
        (
            lambda sites: [
                site for site in sites if not (site.chain == "I" and site.name == "X28")
            ],
            "requires 33 ligand heavy atoms",
        ),
        (
            lambda sites: [
                replace(site, residue="WAT") if site.chain == "J" else site
                for site in sites
            ],
            "Expected one HOH oxygen",
        ),
        (
            lambda sites: [
                replace(site, element="H")
                if site.chain == "A" and site.sequence == 1
                else site
                for site in sites
            ],
            "Unsupported heavy element",
        ),
        (
            lambda sites: [
                replace(site, xyz=(math.nan, 0, 0)) if site.chain == "J" else site
                for site in sites
            ],
            "Coordinates must be finite",
        ),
        (
            lambda sites: [
                replace(site, confidence=101) if site.chain == "J" else site
                for site in sites
            ],
            "Invalid per-atom confidence",
        ),
        (
            lambda sites: [
                site for site in sites if not (site.chain == "E" and site.name == "OG1")
            ],
            "Missing or misidentified required atom",
        ),
    ],
)
def test_invalid_inventory_and_atoms_fail_explicitly(
    mutate: Callable[[list[AtomSite]], list[AtomSite]], message: str
) -> None:
    """Duplicate, missing, malformed, or misidentified sites are rejected."""
    with pytest.raises(ValueError, match=message):
        audit_sites(mutate(_valid_sites()))


def test_bridge_and_contact_threshold_failures_are_reported() -> None:
    """Bridge and elemental contact cutoffs affect screening, not input validity."""
    sites = _valid_sites()
    sites = [
        replace(site, xyz=(2.0, 0.0, 0.0)) if site.chain == "J" else site
        for site in sites
    ]

    report = audit_sites(sites)

    assert report["bridge_distance_screen_passed"] is False
    assert report["water_heavy_contact_screen_passed"] is False
    assert report["water_maximum_contact_shortfall_A"] == pytest.approx(1.7)
    assert report["heavy_atom_screen_passed"] is False
    nitrogen_contact = next(
        row
        for row in report["water_neighbors_within_5A"]
        if row["atom"].endswith("E:1:THR:N")
    )
    assert nitrogen_contact["minimum_A"] == 2.5


def test_elemental_contact_limits_and_tolerance() -> None:
    """N/O, C, and S use their specified minima and 0.10-A tolerance."""
    sites = [
        replace(site, element="S", xyz=(0.0, 0.0, 2.90))
        if site.chain == "I" and site.name == "X28"
        else site
        for site in _valid_sites()
    ]

    report = audit_sites(sites)
    contacts = {row["atom"]: row for row in report["water_neighbors_within_5A"]}

    assert contacts["E:1:THR:N"]["minimum_A"] == 2.5
    assert contacts["E:1:THR:OG1"]["minimum_A"] == 2.5
    assert contacts["I:1:PA6TI1:C54"]["minimum_A"] == 2.75
    assert contacts["I:1:PA6TI1:O49"]["minimum_A"] == 2.5
    assert contacts["I:1:PA6TI1:X28"]["minimum_A"] == 2.95
    assert report["water_maximum_contact_shortfall_A"] == pytest.approx(0.05)
    assert report["water_heavy_contact_screen_passed"] is True


def test_reactive_bond_outside_screen_is_reported() -> None:
    """The covalent-bond distance screen fails without changing the model."""
    sites = [
        replace(site, xyz=(3.0, 4.4, 0.0))
        if site.chain == "I" and site.name == "O49"
        else site
        for site in _valid_sites()
    ]

    report = audit_sites(sites)

    assert report["reactive_bond_screen_passed"] is False
    assert report["reactive_bonds"][1]["distance_A"] == pytest.approx(3.0)
    assert report["accepted_for_qmmm"] is False
