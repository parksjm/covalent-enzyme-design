"""Classify final ASH272 QM/MM E-S chemistry and catalytic geometry."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from scripts.queue_modular_es_qmmm_after_mm import endpoint_geometry  # noqa: E402

BASE = ROOT / "qmmm/pa6_flipped_model30/06_ash272_comparison"
CASES = ("replica_1_tier_A", "replica_1_tier_B", "replica_2_tier_A", "replica_2_tier_B")


def main() -> None:
    """Write one JSON report with current Tier and chemical-retention checks."""
    topology = BASE / "00_system/ash272_solvated.prmtop"
    output: list[dict[str, object]] = []
    for name in CASES:
        case = BASE / "06_qmmm_optimizations" / name
        calculation = case / ("attempt_01/calculation" if "tier_B" in name else "")
        restart = calculation / "optimized.ncrst"
        parm = AmberParm(str(topology), xyz=str(restart))

        def point(residue: int, atom_name: str) -> np.ndarray:
            atoms = [a for a in parm.residues[residue - 1].atoms if a.name == atom_name]
            if len(atoms) != 1:
                raise ValueError(f"{name}: {residue} {atom_name} count {len(atoms)}")
            return np.asarray(parm.coordinates[atoms[0].idx], dtype=float)

        def distance(a: tuple[int, str], b: tuple[int, str]) -> float:
            return float(np.linalg.norm(point(*a) - point(*b)))

        bonds = {
            "thr_og1_hg1_A": distance((929, "OG1"), (929, "HG1")),
            "ash272_od2_hd2_A": distance((968, "OD2"), (968, "HD2")),
            "asp274_od2_ash272_hd2_A": distance((970, "OD2"), (968, "HD2")),
            "thr_n_h1_A": distance((929, "N"), (929, "H1")),
            "thr_n_h2_A": distance((929, "N"), (929, "H2")),
            "substrate_c_n_A": distance((1314, "C"), (1315, "N")),
        }
        chemistry = (
            bonds["thr_og1_hg1_A"] < 1.25
            and bonds["ash272_od2_hd2_A"] < 1.25
            and bonds["thr_n_h1_A"] < 1.25
            and bonds["thr_n_h2_A"] < 1.25
            and bonds["substrate_c_n_A"] < 1.6
        )
        geometry = endpoint_geometry(topology, restart)
        output.append(
            {
                "case": name,
                "tier": geometry["tier"],
                "geometry": geometry["metrics"],
                "bonds": bonds,
                "reactant_chemistry_retained": chemistry,
                "ash272_proton_transferred_to_asp274": (
                    bonds["ash272_od2_hd2_A"] > 1.3
                    and bonds["asp274_od2_ash272_hd2_A"] < 1.25
                ),
                "restart": str(restart),
            }
        )
    destination = BASE / "06_qmmm_optimizations/endpoint_chemistry_geometry.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
