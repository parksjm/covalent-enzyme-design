"""Chemistry and geometric acceptance for the alpha/beta TI1 repair experiment."""

from __future__ import annotations

import numpy as np
from parmed.amber import AmberMask, AmberParm

from .active_site_geometry import angle_degrees
from .local_repair_model import QM_MASK
from .relay_water import _contact_limits


def audit_local_model(model: AmberParm) -> dict:
    """Check dimer chemistry and report water geometry independently of energy.

    Args:
        model: Hydrogenated 330-residue alpha/beta PA6 L4 TI1 topology/coordinates.

    Returns:
        Inventory acceptance, bridge/contact geometry, and reactive bond screens.

    Raises:
        ValueError: If the covalent or proton inventory is incompatible with TI1.
    """
    xyz = np.asarray(model.coordinates)
    if (
        len(model.residues) != 330
        or xyz.shape != (len(model.atoms), 3)
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Invalid dimer inventory or coordinates")
    thr, ligand, water = [model.residues[i] for i in (232, 328, 329)]
    t, lig, wat = [{a.name: a for a in r.atoms} for r in (thr, ligand, water)]
    if (thr.name, ligand.name, water.name) != ("THR", "LIG", "WAT") or len(
        ligand.atoms
    ) != 79:
        raise ValueError("Unexpected reacting residues or PA6 L4 inventory")
    if set(wat) != {"O", "H1", "H2"} or "HG1" in t:
        raise ValueError("Incorrect water/OG1 proton inventory")
    if {a.name for a in t["N"].bond_partners if a.atomic_number == 1} != {
        "H1",
        "H2",
        "H3",
    }:
        raise ValueError("Reacting terminal N must be NH3+")
    carbon = lig["C54"]
    if set(carbon.bond_partners) != {t["OG1"], lig["O49"], lig["N59"], lig["C62"]}:
        raise ValueError("TI1 tetrahedral connectivity changed")
    if lig["O49"].bond_partners != [carbon] or any(
        a.atomic_number == 1 for a in t["OG1"].bond_partners
    ):
        raise ValueError("Incorrect TI1 oxyanion or OG1 connectivity")
    if set(wat["O"].bond_partners) != {wat["H1"], wat["H2"]}:
        raise ValueError("Water covalent inventory changed")
    selected = set(AmberMask(model, QM_MASK).Selected())
    cuts = [
        b for b in model.bonds if (b.atom1.idx in selected) != (b.atom2.idx in selected)
    ]
    if (
        len(selected) != 109
        or len(cuts) != 1
        or {(a.residue.idx, a.name) for a in (cuts[0].atom1, cuts[0].atom2)}
        != {(233, "CA"), (233, "C")}
    ):
        raise ValueError("QM selection/boundary mismatch")
    bonds = []
    pairs = [(carbon, a, 1.2, 1.8) for a in carbon.bond_partners]
    pairs += [(t["N"], t[h], 0.8, 1.25) for h in ("H1", "H2", "H3")]
    pairs += [(wat["O"], wat[h], 0.8, 1.2) for h in ("H1", "H2")]
    for first, second, lower, upper in pairs:
        length = float(np.linalg.norm(xyz[first.idx] - xyz[second.idx]))
        bonds.append(
            {
                "atoms": [first.idx + 1, second.idx + 1],
                "names": [first.name, second.name],
                "distance_A": length,
                "accepted": lower <= length <= upper,
            }
        )
    environment = [a.idx for a in model.atoms if a.residue.idx != 329]
    o_limits, h_limits = _contact_limits(model, environment)
    distances = np.linalg.norm(
        xyz[[wat[n].idx for n in ("O", "H1", "H2")]][:, None, :]
        - xyz[environment][None, :, :],
        axis=2,
    )
    shortfall = np.maximum(np.stack([o_limits, h_limits, h_limits]) - distances, 0)
    contacts = []
    for wi, ei in np.argwhere(shortfall > 0.1):
        atom = model.atoms[environment[ei]]
        contacts.append(
            {
                "water_atom": ("O", "H1", "H2")[wi],
                "atom": f"{atom.residue.idx + 1}:{atom.name}",
                "distance_A": float(distances[wi, ei]),
                "shortfall_A": float(shortfall[wi, ei]),
            }
        )
    bridge = {
        name: float(np.linalg.norm(xyz[t[name].idx] - xyz[wat["O"].idx]))
        for name in ("N", "OG1")
    }
    donor_angles = [
        angle_degrees(xyz[t["N"].idx], xyz[t[h].idx], xyz[wat["O"].idx])
        for h in ("H1", "H2", "H3")
    ]
    acceptor_angles = [
        angle_degrees(xyz[wat["O"].idx], xyz[wat[h].idx], xyz[t["OG1"].idx])
        for h in ("H1", "H2")
    ]
    max_shortfall = float(shortfall.max())
    bond_pass = all(b["accepted"] for b in bonds)
    bridge_pass = all(2.5 <= d <= 3.5 for d in bridge.values())
    angle_pass = max(donor_angles) >= 130 and max(acceptor_angles) >= 130
    chirality = {}
    for center, neighbors in (
        (t["CA"], [t[n] for n in ("N", "C", "CB")]),
        (t["CB"], [t[n] for n in ("CA", "OG1", "CG2")]),
        (carbon, [t["OG1"], lig["O49"], lig["N59"]]),
    ):
        vectors = np.array([xyz[a.idx] - xyz[center.idx] for a in neighbors])
        chirality[f"{center.residue.idx + 1}:{center.name}"] = float(
            np.linalg.det(vectors)
        )
    return {
        "inventory_passed": True,
        "qm_atoms": len(selected),
        "bonds": bonds,
        "bond_screen_passed": bond_pass,
        "bridge_distances_A": bridge,
        "bridge_screen_passed": bridge_pass,
        "N_H_water_angles_deg": donor_angles,
        "water_H_OG1_angles_deg": acceptor_angles,
        "directional_hbond_screen_passed": angle_pass,
        "water_maximum_contact_shortfall_A": max_shortfall,
        "water_contact_tolerance_A": 0.1,
        "water_contacts_over_tolerance": contacts,
        "signed_volumes_A3": chirality,
        "geometry_screen_passed": bond_pass
        and bridge_pass
        and angle_pass
        and max_shortfall <= 0.1,
        "stationary_point_validated": False,
    }


def local_restraint_selection(model: AmberParm) -> dict:
    """Select nearby residues and outer heavy atoms for positional restraints.

    Args:
        model: Prepared dimer with fixed atom ordering.

    Returns:
        Explicit atom indices and masks; all hydrogens remain mobile.
    """
    xyz = np.asarray(model.coordinates)
    centers = [a.idx for a in model.atoms if a.residue.idx in (232, 328, 329)]
    distances = np.linalg.norm(xyz[:, None, :] - xyz[centers][None, :, :], axis=2).min(
        axis=1
    )
    residues = {
        a.residue.idx
        for a in model.atoms
        if a.atomic_number > 1 and distances[a.idx] <= 5.0
    }
    residues.update((232, 233, 328, 329))
    outer = [
        a.idx + 1
        for a in model.atoms
        if a.atomic_number > 1 and a.residue.idx not in residues
    ]
    all_heavy = [
        a.idx + 1 for a in model.atoms if a.atomic_number > 1 and a.residue.idx != 329
    ]
    if not outer:
        raise ValueError("Outer positional restraint selection is empty")
    return {
        "local_residues_1based": sorted(i + 1 for i in residues),
        "local_radius_A": 5.0,
        "outer_heavy_atom_indices_1based": outer,
        "initial_heavy_atom_indices_1based": all_heavy,
        "bridge_atoms_1based": [
            [
                next(a.idx + 1 for a in model.residues[232].atoms if a.name == name),
                model.residues[329].atoms[0].idx + 1,
            ]
            for name in ("N", "OG1")
        ],
        "local_sidechains_and_backbones_mobile": True,
        "all_hydrogens_mobile": True,
    }
