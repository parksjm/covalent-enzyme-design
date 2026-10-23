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
    if (
        molecule.GetNumAtoms() != 63
        or not {"WAT:O", "WAT:H1", "WAT:H2"} <= origins.keys()
    ):
        raise ValueError("Expected the explicit-H 63-atom donor/water model")
    elements, xyz = final_orientation(text)
    geometry = validate_geometry(molecule, elements, xyz)
    frozen = policy["frozen_atom_indices_0based"]
    frozen_shift = float(np.max(np.linalg.norm(xyz[frozen] - start[frozen], axis=1)))

    def conformer(points: np.ndarray) -> Chem.Conformer:
        result = Chem.Conformer(molecule.GetNumAtoms())
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
    states = re.findall(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", text)
    ammonium = {
        name: {
            "N_H_A": float(np.linalg.norm(xyz[origins[name]] - xyz[origins["THR:N"]])),
            "OG1_H_A": float(
                np.linalg.norm(xyz[origins[name]] - xyz[origins["THR:OG1"]])
            ),
        }
        for name in ("THR:H1", "THR:H2", "THR:HG1")
    }
    limits = checks["ammonium_proton_limits"]
    water_o = xyz[origins["WAT:O"]]
    water_dist = {name: float(np.linalg.norm(xyz[origins[name]]-water_o)) for name in ("THR:N", "THR:OG1")}
    nh_angle = max(angle_degrees(xyz[origins["THR:N"]], xyz[origins[h]], water_o) for h in ("THR:H1", "THR:H2", "THR:HG1"))
    oh_angle = max(angle_degrees(water_o, xyz[origins[h]], xyz[origins["THR:OG1"]]) for h in ("WAT:H1", "WAT:H2"))
    passed = {
        "water_bridge": all(2.5 <= d <= 3.3 for d in water_dist.values()) and nh_angle >= 140 and oh_angle >= 140,
        "ammonium_retained": all(
            limits["N_H_minimum_A"] <= row["N_H_A"] <= limits["N_H_maximum_A"]
            and row["OG1_H_A"] >= limits["OG1_H_minimum_A"]
            for row in ammonium.values()
        ),
        "normal_termination": text.rstrip()
        .splitlines()[-1]
        .strip()
        .startswith("Normal termination of Gaussian 16"),
        "optimized": "Optimization completed." in text,
        "electronic_state": bool(states)
        and all(state == ("0", "1") for state in states),
        "graph_and_stereo": geometry["accepted"],
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
        "water_bridge_distances_A": water_dist,
        "N_H_Owater_angle_deg": nh_angle,
        "Owater_H_OG1_angle_deg": oh_angle,
        "reactive_distances_A": distances,
        "ammonium_proton_distances_A": ammonium,
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
    lines = [
        str(molecule.GetNumAtoms()),
        "Native TI1 frozen-scaffold pilot; see scientific_acceptance.json",
    ]
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
