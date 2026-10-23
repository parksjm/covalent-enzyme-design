"""Assess one native frozen-scaffold optimization without certifying a force field."""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import GetDihedralDeg

from nyl12_covalent.active_site_geometry import angle_degrees
from nyl12_covalent.gaussian_checks import final_orientation, validate_geometry

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Write separate scientific acceptance and a full-hydrogen final structure."""
    root = Path(__file__).resolve().parent
    text = (root / "optimization.log").read_text()
    molecule = Chem.SDMolSupplier(str(root / "native_seed.sdf"), removeHs=False)[0]
    if molecule is None:
        raise ValueError("Cannot read native graph")
    policy = json.loads((root / "pilot_policy.json").read_text())
    origins = json.loads((root / "atom_map.json").read_text())["atom_index_by_origin"]
    start = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (root / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    elements, xyz = final_orientation(text)
    geometry = validate_geometry(molecule, elements, xyz)
    frozen = policy["frozen_atom_indices_0based"]
    frozen_shift = float(np.max(np.linalg.norm(xyz[frozen] - start[frozen], axis=1)))

    def conformer(points: np.ndarray) -> Chem.Conformer:
        result = Chem.Conformer(103)
        result.Set3D(True)
        for idx, point in enumerate(points):
            result.SetAtomPosition(idx, tuple(float(v) for v in point))
        return result

    indices = [i - 1 for i in policy["chi1_constraint_indices_1based"]]
    delta = GetDihedralDeg(conformer(xyz), *indices) - GetDihedralDeg(
        conformer(start), *indices
    )
    chi_error = abs((delta + 180) % 360 - 180)
    checks = policy["prospective_checks"]
    distances = {}
    for label in checks["reactive_bond_windows_A"]:
        first, second = label.split("_")
        left = "THR:OG1" if first == "OG1" else f"PA6:{first}"
        distances[label] = float(
            np.linalg.norm(xyz[origins[left]] - xyz[origins[f"PA6:{second}"]])
        )
    hydrogen = [a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum() == 1]
    min_oh = float(
        np.min(np.linalg.norm(xyz[hydrogen] - xyz[origins["PA6:O49"]], axis=1))
    )
    water_o = xyz[origins["WAT:O"]]
    bridge = {
        name: float(np.linalg.norm(xyz[origins[f"THR:{name}"]] - water_o))
        for name in ("N", "OG1")
    }
    donor_angles = [
        angle_degrees(xyz[origins["THR:N"]], xyz[origins[h]], water_o)
        for h in ("THR:H1", "THR:H2", "THR:HG1")
    ]
    acceptor_angles = [
        angle_degrees(water_o, xyz[origins[h]], xyz[origins["THR:OG1"]])
        for h in ("WAT:H1", "WAT:H2")
    ]
    relay = {
        "distances_A": bridge,
        "N_H_water_angles_deg": donor_angles,
        "water_H_OG1_angles_deg": acceptor_angles,
        "accepted": all(2.5 <= d <= 3.5 for d in bridge.values())
        and max(donor_angles) >= 130
        and max(acceptor_angles) >= 130,
    }
    passed = {
        "normal_termination": text.rstrip()
        .splitlines()[-1]
        .strip()
        .startswith("Normal termination of Gaussian 16"),
        "optimized": "Optimization completed." in text,
        "electronic_state": bool(
            re.search(r"Charge\s*=\s*0\s+Multiplicity\s*=\s*1", text)
        ),
        "graph_and_stereo": geometry["accepted"],
        "relay_water_geometry": relay["accepted"],
        "fixed_coordinates": frozen_shift <= checks["frozen_atom_tolerance_A"],
        "chi1": chi_error <= checks["chi1_tolerance_deg"],
        "reaction_center": all(
            checks["reactive_bond_windows_A"][name][0]
            <= value
            <= checks["reactive_bond_windows_A"][name][1]
            for name, value in distances.items()
        ),
        "oxyanion_not_protonated_screen": min_oh >= checks["oxyanion_H_minimum_A"],
    }
    report = {
        "accepted": all(passed.values()),
        "scope": policy["scope"],
        "checks": passed,
        "geometry": geometry,
        "relay_water": relay,
        "reactive_distances_A": distances,
        "max_frozen_displacement_A": frozen_shift,
        "chi1_error_deg": chi_error,
        "minimum_oxyanion_H_distance_A": min_oh,
        "forcefield_accepted": False,
        "charge_fit_reference_accepted": False,
    }
    (root / "scientific_acceptance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    molecule.RemoveAllConformers()
    molecule.AddConformer(conformer(xyz))
    with Chem.SDWriter(str(root / "final_structure.sdf")) as writer:
        writer.write(molecule)
    lines = ["103", "Native TI1 frozen-scaffold pilot; see scientific_acceptance.json"]
    lines.extend(
        f"{a.GetSymbol():2s} " + " ".join(f"{v:.10f}" for v in point)
        for a, point in zip(molecule.GetAtoms(), xyz, strict=True)
    )
    (root / "final_structure.xyz").write_text("\n".join(lines) + "\n")
    LOGGER.info("Native pilot scientific acceptance: %s", report["accepted"])
    if not report["accepted"]:
        raise SystemExit(3)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
