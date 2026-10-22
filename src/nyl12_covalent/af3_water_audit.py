"""Independent heavy-atom screening of the minimal Nyl12 AF3 water test."""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import math

import numpy as np

from .active_site_geometry import angle_degrees, distance, measure_active_site
from .af3_postprocess import AtomSite


def audit_sites(sites: list[AtomSite]) -> dict:
    """Audit water placement and TI1 geometry without altering coordinates.

    Args:
        sites: Heavy atoms from the single AF3 alpha-A/beta-E/PA6-I/water-J model.

    Returns:
        Independent screening results, descriptive metrics, and limitations.
        A heavy-atom pass never authorizes QM/MM without hydrogen/topology checks.

    Raises:
        ValueError: If atom inventory, identities, coordinates or confidence fail.
    """
    atoms = {(a.chain, a.sequence, a.name): a for a in sites}
    if len(atoms) != len(sites):
        raise ValueError("Duplicate atom identity")
    if {a.chain for a in sites} != {"A", "E", "I", "J"}:
        raise ValueError("Expected exactly chains A, E, I, J")
    for chain, count in (("A", 232), ("E", 96), ("I", 1), ("J", 1)):
        if {a.sequence for a in sites if a.chain == chain} != set(range(1, count + 1)):
            raise ValueError(f"Unexpected residue inventory for chain {chain}")
    if sum(a.chain == "I" for a in sites) != 33:
        raise ValueError("PA6 L4 TI1 requires 33 ligand heavy atoms")
    water_sites = [a for a in sites if a.chain == "J"]
    if len(water_sites) != 1 or (
        water_sites[0].name,
        water_sites[0].element,
        water_sites[0].residue,
    ) != ("O", "O", "HOH"):
        raise ValueError("Expected one HOH oxygen J1:O")
    for atom in sites:
        if atom.element not in {"C", "N", "O", "S"}:
            raise ValueError(f"Unsupported heavy element {atom.element}")
        if len(atom.xyz) != 3 or not all(math.isfinite(v) for v in atom.xyz):
            raise ValueError("Coordinates must be finite three-vectors")
        if not math.isfinite(atom.confidence) or not 0 <= atom.confidence <= 100:
            raise ValueError("Invalid per-atom confidence")
    required = {
        ("E", 1, "N"): "THR",
        ("E", 1, "OG1"): "THR",
        ("A", 113, "N"): "PHE",
        ("A", 186, "ND2"): "ASN",
    }
    required.update({("I", 1, name): "PA6TI1" for name in ("C54", "O49", "N59", "C62")})
    for identity, residue in required.items():
        if identity not in atoms or atoms[identity].residue != residue:
            raise ValueError(f"Missing or misidentified required atom {identity}")
    water = water_sites[0]
    nitrogen, og1 = atoms["E", 1, "N"], atoms["E", 1, "OG1"]
    oxygen, carbon = atoms["I", 1, "O49"], atoms["I", 1, "C54"]
    contacts = []
    for atom in sites:
        if atom.chain == "J":
            continue
        separation = distance(water.xyz, atom.xyz)
        threshold = {"C": 2.75, "N": 2.50, "O": 2.50, "S": 2.95}[atom.element]
        if separation <= 5.0:
            contacts.append(
                {
                    "atom": f"{atom.chain}:{atom.sequence}:{atom.residue}:{atom.name}",
                    "distance_A": separation,
                    "minimum_A": threshold,
                    "shortfall_A": max(0.0, threshold - separation),
                }
            )
    contacts.sort(key=lambda row: row["distance_A"])
    maximum = max((row["shortfall_A"] for row in contacts), default=0.0)
    bridge = {
        "N_water_O_A": distance(nitrogen.xyz, water.xyz),
        "OG1_water_O_A": distance(og1.xyz, water.xyz),
    }
    bridge_pass = all(2.6 <= value <= 3.2 for value in bridge.values())
    neighbors = [og1, oxygen, atoms["I", 1, "N59"], atoms["I", 1, "C62"]]
    bonds = [
        {
            "atoms": ["I:1:C54", f"{a.chain}:{a.sequence}:{a.name}"],
            "distance_A": distance(carbon.xyz, a.xyz),
        }
        for a in neighbors
    ]
    bonds_pass = all(1.2 <= row["distance_A"] <= 1.8 for row in bonds)
    angles = [
        {
            "atoms": [a.name, "C54", b.name],
            "angle_deg": angle_degrees(a.xyz, carbon.xyz, b.xyz),
        }
        for a, b in combinations(neighbors, 2)
    ]
    vectors = np.array([np.subtract(a.xyz, carbon.xyz) for a in neighbors[:3]])
    norms = np.linalg.norm(vectors, axis=1)
    if np.any(norms < 1e-10):
        raise ValueError("Coincident reactive atoms")
    volume = float(np.linalg.det(vectors / norms[:, None]))
    coordinates = {
        "Thr929.OG1": og1.xyz,
        "Phe113.N": atoms["A", 113, "N"].xyz,
        "Asn186.ND2": atoms["A", 186, "ND2"].xyz,
    }
    coordinates.update(
        {
            f"substrate.{name}": atoms["I", 1, name].xyz
            for name in ("C54", "O49", "N59", "C62")
        }
    )
    active = measure_active_site(coordinates, "pa6", "ti1")
    confidence = {}
    for chain in ("A", "E", "I", "J"):
        values = np.array([a.confidence for a in sites if a.chain == chain])
        confidence[chain] = {
            "count": len(values),
            "mean_plddt": float(values.mean()),
            "minimum_plddt": float(values.min()),
            "maximum_plddt": float(values.max()),
        }
    return {
        "inventory_passed": True,
        "chain_heavy_atom_counts": dict(Counter(a.chain for a in sites)),
        "bridge_distances": bridge,
        "bridge_restraint_range_A": [2.6, 3.2],
        "bridge_distance_screen_passed": bridge_pass,
        "N_water_O_OG1_angle_deg": angle_degrees(nitrogen.xyz, water.xyz, og1.xyz),
        "water_oxyanion_distance_A": distance(water.xyz, oxygen.xyz),
        "water_neighbors_within_5A": contacts,
        "water_maximum_contact_shortfall_A": maximum,
        "water_contact_tolerance_A": 0.10,
        "water_heavy_contact_screen_passed": maximum <= 0.10,
        "reactive_bonds": bonds,
        "reactive_bond_range_A": [1.2, 1.8],
        "reactive_bond_screen_passed": bonds_pass,
        "reactive_angles": angles,
        "signed_normalized_OG1_O49_N59_volume": volume,
        "active_site": active,
        "confidence_from_rounded_cif": confidence,
        "heavy_atom_screen_passed": bridge_pass and maximum <= 0.10 and bonds_pass,
        "accepted_for_qmmm": False,
        "limitations": [
            "Heavy atoms only: neither water H direction nor terminal protonation is established.",
            "Covalent connectivity is declared in AF3 input; distance screening does not validate bond order or parameters.",
            "Angles and signed volume are descriptive; no stereochemical or stationary-point acceptance is inferred.",
            "Restrained prediction and a single alpha/beta pair do not establish unrestrained binding or octamer compatibility.",
        ],
    }
