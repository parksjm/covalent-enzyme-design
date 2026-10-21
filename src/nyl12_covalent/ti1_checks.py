"""Read-only TI1 inventory, geometric screening, and Amber convergence gates."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm

from .relay_water import _contact_limits

TARGET_ATOMS = {
    "pa6_l2": ("C54", "O49", "N59", 79),
    "poly_xama_maxa": ("C033", "O034", "N035", 151),
}


def check_geometry(topology: Path, coordinates: Path, target: str) -> dict[str, Any]:
    """Screen TI1 inventory and coordinates without modifying either input.

    Args:
        topology: Previously product-gated Amber topology.
        coordinates: Amber coordinates in topology order.
        target: Measured-product register identifier.

    Returns:
        Inventory, bond-length, and relay-contact results. These are screening
        gates, not a validation of parameters or of a stationary point.

    Raises:
        ValueError: If topology inventory or coordinate shape is invalid.
    """
    carbon_name, oxygen_name, nitrogen_name, ligand_count = TARGET_ATOMS[target]
    model = AmberParm(str(topology), xyz=str(coordinates))
    xyz = np.asarray(model.coordinates)
    if xyz.shape != (len(model.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Coordinates must be finite and match topology atom count")
    if len(model.residues) != 1314:
        raise ValueError("Expected exactly 1314 TI1 residues")
    thr, ligand, water = (model.residues[i] for i in (928, 1312, 1313))
    if (thr.name, ligand.name, water.name) != ("THR", "LIG", "WAT"):
        raise ValueError("Unexpected catalytic residue, substrate, or water")
    if len(ligand.atoms) != ligand_count:
        raise ValueError("Ligand inventory does not match the product-gated target")
    t, lig, wat = ({a.name: a for a in r.atoms} for r in (thr, ligand, water))
    if len(lig) != len(ligand.atoms) or set(wat) != {"O", "H1", "H2"}:
        raise ValueError("Duplicate ligand names or non-H2O relay inventory")
    if "HG1" in t or not {"H1", "H2", "H3"}.issubset(t):
        raise ValueError("TI1 requires beta NH3+ and no OG1 proton")
    if {a.name for a in t["N"].bond_partners if a.atomic_number == 1} != {
        "H1",
        "H2",
        "H3",
    }:
        raise ValueError("TI1 beta N proton connectivity changed")
    for index in (1024, 1120, 1216):
        residue = model.residues[index]
        atoms = {a.name: a for a in residue.atoms}
        if residue.name != "THR" or "H3" in atoms:
            raise ValueError("Nonreacting beta terminus must remain neutral NH2")
        if {a.name for a in atoms["N"].bond_partners if a.atomic_number == 1} != {
            "H1",
            "H2",
        } or abs(sum(a.charge for a in residue.atoms)) > 0.002:
            raise ValueError("Nonreacting beta NH2 connectivity or charge changed")
    carbon, oxygen, nitrogen = (
        lig[n] for n in (carbon_name, oxygen_name, nitrogen_name)
    )
    if len(carbon.bond_partners) != 4 or not {t["OG1"], oxygen, nitrogen}.issubset(
        carbon.bond_partners
    ):
        raise ValueError("TI1 tetrahedral carbon lost a required bond")
    if oxygen.bond_partners != [carbon] or any(
        a.atomic_number == 1 for a in t["OG1"].bond_partners
    ):
        raise ValueError("TI1 oxyanion or OG1 protonation changed")
    if (
        set(wat["O"].bond_partners) != {wat["H1"], wat["H2"]}
        or abs(sum(a.charge for a in water.atoms)) > 0.001
    ):
        raise ValueError("Relay water must retain neutral H2O connectivity")

    bonds = []
    pairs = [(carbon, partner, 1.20, 1.80) for partner in carbon.bond_partners]
    pairs += [(t["N"], t[h], 0.80, 1.25) for h in ("H1", "H2", "H3")]
    pairs += [(wat["O"], wat[h], 0.80, 1.20) for h in ("H1", "H2")]
    for first, second, lower, upper in pairs:
        length = float(np.linalg.norm(xyz[first.idx] - xyz[second.idx]))
        bonds.append(
            {
                "atoms": [first.idx + 1, second.idx + 1],
                "names": [first.name, second.name],
                "distance_A": length,
                "range_A": [lower, upper],
                "accepted": lower <= length <= upper,
            }
        )
    environment = [a.idx for a in model.atoms if a.residue.idx != 1313]
    o_limits, h_limits = _contact_limits(model, environment)
    distances = np.linalg.norm(
        xyz[[wat[n].idx for n in ("O", "H1", "H2")]][:, None, :]
        - xyz[environment][None, :, :],
        axis=2,
    )
    shortfall = np.maximum(np.stack([o_limits, h_limits, h_limits]) - distances, 0)
    maximum = float(shortfall.max())
    relay_distances = {
        name: float(np.linalg.norm(xyz[t[name].idx] - xyz[wat["O"].idx]))
        for name in ("N", "OG1")
    }
    relay_ok = all(2.5 <= distance <= 3.5 for distance in relay_distances.values())
    return {
        "accepted": all(b["accepted"] for b in bonds) and maximum <= 0.10 and relay_ok,
        "inventory_passed": True,
        "bonds": bonds,
        "water_maximum_contact_shortfall_A": maximum,
        "water_allowed_contact_shortfall_A": 0.10,
        "relay_distances_A": relay_distances,
        "relay_distance_range_A": [2.5, 3.5],
        "limitations": "Geometric screen only; no energetic or mechanism validation. "
        "Product identity is inherited from the upstream product-gated topology.",
    }


def check_convergence(mdout: Path, drms: float = 0.05) -> dict[str, Any]:
    """Require a finite final Amber minimization RMS below the specified target.

    Args:
        mdout: Sander text output.
        drms: RMS gradient convergence threshold used in the input.

    Returns:
        Explicit convergence result; absent final results fail closed.
    """
    text = mdout.read_text()
    final = text.rsplit("FINAL RESULTS", 1)
    matches = (
        []
        if len(final) != 2
        else re.findall(
            r"NSTEP\s+ENERGY\s+RMS\s+GMAX[^\n]*\n\s*"
            r"(\d+)\s+(\S+)\s+(\S+)\s+(\S+)",
            final[1],
        )
    )
    values = None
    if matches:
        try:
            values = [float(v.replace("D", "E")) for v in matches[-1][1:]]
        except ValueError:
            values = None  # Asterisks/overflow are a recorded convergence failure.
    if values is not None and not np.isfinite(values).all():
        values = None
    failures = re.findall(
        r"(?im)^.*(?:SCF[^\n]*(?:fail|not converg)|QM[^\n]*error|"
        r"SANDER BOMB|abnormal termination).*$",
        text,
    )
    accepted = bool(
        values is not None
        and np.isfinite(values).all()
        and values[2] >= 0
        and 0 <= values[1] <= drms
        and not failures
    )
    return {
        "accepted": accepted,
        "final_energy_rms_gmax": values,
        "drms": drms,
        "engine_failure_lines": failures,
        "final_results_present": len(final) == 2,
    }
