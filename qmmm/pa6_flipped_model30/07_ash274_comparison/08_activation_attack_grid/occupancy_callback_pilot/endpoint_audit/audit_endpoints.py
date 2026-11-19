"""Audit retained occupancy-pilot endpoints without invoking an engine."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
import sys
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "endpoint_audit"
sys.path[:0] = [str(ROOT / "src"), str(ROOT), str(BASE)]
from runtime_support import SOURCE, configuration, load_model, state_record  # noqa: E402
from nyl12_covalent.occupancy_bias import evaluate_bias  # noqa: E402

EH = 627.5094740631
HASHES: dict[str, str] = {}


def track(path: Path) -> Path:
    """Record the exact source identity before reading evidence."""
    HASHES[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read_json(path: Path) -> Any:
    """Read hashed JSON evidence."""
    return json.loads(track(path).read_text())


def coordinates(path: Path) -> np.ndarray:
    """Read hashed Amber restart coordinates in Angstroms."""
    return np.asarray(Rst7.open(str(track(path))).coordinates).reshape(-1, 3)


def distance(x: np.ndarray, a: int, b: int) -> float:
    """Return a distance using one-based Amber atom identifiers."""
    return float(np.linalg.norm(x[a - 1] - x[b - 1]))


def angle(x: np.ndarray, a: int, b: int, c: int) -> float:
    """Return the a-b-c angle in degrees."""
    u, v = x[a - 1] - x[b - 1], x[c - 1] - x[b - 1]
    return float(
        np.degrees(
            np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))
        )
    )


def geometry(x: np.ndarray, neighbors: list[int]) -> dict[str, Any]:
    """Measure reactive geometry and carbonyl tetrahedralization explicitly."""
    center = x[19919]
    vectors = x[np.array(neighbors) - 1] - center
    normal = np.cross(vectors[1] - vectors[0], vectors[2] - vectors[0])
    height = float(np.dot(-vectors[0], normal) / np.linalg.norm(normal))
    pairs = {
        f"{a}-{b}": angle(x, a, 19920, b)
        for a, b in combinations(neighbors + [13916], 2)
    }
    distances = {
        f"{a}-{b}": distance(x, a, b)
        for a, b in [
            (13916, 19920),
            (19920, 19921),
            (19920, 19933),
            (19921, 1638),
            (19921, 2732),
            (13916, 20317),
            (13905, 20317),
            (13905, 14559),
            (13905, 14560),
            (14560, 14579),
            (14578, 14580),
            (14579, 14580),
        ]
    }
    for h in [13906, 13907, 13917, 20318, 20319]:
        for a in [13905, 13916, 20317, 14559, 14560]:
            distances[f"{h}-{a}"] = distance(x, h, a)
    return {
        "distances_A": distances,
        "burgi_dunitz_OG1_C_O_deg": angle(x, 13916, 19920, 19921),
        "CB_OG1_C_deg": angle(x, 13910, 13916, 19920),
        "carbonyl_original_neighbors_1based": neighbors,
        "carbonyl_signed_neighbor_plane_height_A": height,
        "carbonyl_signed_volume_A3": float(np.linalg.det(vectors) / 6),
        "carbonyl_original_neighbor_angle_sum_deg": sum(
            angle(x, a, 19920, b) for a, b in combinations(neighbors, 2)
        ),
        "carbonyl_four_neighbor_angles_deg": pairs,
    }


def main() -> None:
    """Write authoritative JSON and compact CSV for original ES and A-D."""
    model = load_model()
    for path in [
        Path(__file__),
        BASE / "runtime_support.py",
        ROOT / "src/nyl12_covalent/occupancy_bias.py",
        *model.source_paths,
    ]:
        track(path)
    original = coordinates(SOURCE / "optimized.ncrst")
    original_state = state_record(model, original)
    original_owners = {
        h["hydrogen_1based"]: h["nearest_qm_heavy_1based"]
        for h in original_state["proton_owners"]
    }
    neighbors = sorted(
        a.idx + 1
        for a in model.topology.atoms[19919].bond_partners
        if a.atomic_number != 1
    )
    assert neighbors == [19919, 19921, 19933]
    mdout = track(SOURCE / "minimized.mdout").read_text()
    es_energy = float(re.findall(r"Etot\s*=\s*([-\d.]+)", mdout)[-1])
    assert es_energy == -450493.0781
    original_eh = float(
        re.findall(r"Energy calculation finished, energy:\s*([-+\d.E]+)", mdout)[-1]
    )
    points: dict[str, Any] = {
        "original_ES": {
            "geometry": geometry(original, neighbors),
            "state": original_state,
            "physical_kcal_mol": es_energy,
            "native_Eh": original_eh,
        }
    }
    summaries = []
    first_a = None
    final_a = None
    for name, target in zip("ABCD", [3.125178400033067, 2.6, 2.1, 1.6], strict=True):
        directory = BASE / "pilot" / f"point_{name}"
        accepted = read_json(directory / "scientific_acceptance.json")
        process = read_json(directory / "process.json")
        paths = list(directory.glob("accepted_*_continuation_final.ncrst"))
        assert len(paths) == 1
        xyz = coordinates(paths[0])
        records = [
            json.loads(line)
            for line in track(directory / "accepted_iterations.jsonl")
            .read_text()
            .splitlines()
        ]
        record = [r for r in records if r["label"] == "continuation_final"][-1]
        native = [
            json.loads(line)
            for line in track(directory / "native_evidence.jsonl")
            .read_text()
            .splitlines()
        ]
        evidence = [
            r for r in native if r["evaluation"] == record["native_evaluations"]
        ][-1]
        number = evidence["evaluation"]
        for filename, sha in evidence["artifacts_sha256"].items():
            path = track(directory / filename)
            assert HASHES[str(path)] == sha
        assert np.array_equal(
            xyz, coordinates(directory / f"eval_{number:05d}.input.ncrst")
        )
        parsed_eh = float(
            re.findall(
                r"Energy calculation finished, energy:\s*([-+\d.E]+)",
                (directory / f"eval_{number:05d}.mdout").read_text(),
            )[-1]
        )
        assert parsed_eh == evidence["energy_hartree"]
        assert abs(record["physical_energy_kcal_mol"] - parsed_eh * EH) < 1e-7
        bias = evaluate_bias(xyz, configuration(model, 2.1, target))
        assert abs(bias.energy - record["bias_energy_kcal_mol"]) < 1e-9
        state = state_record(model, xyz)
        assert state == accepted["endpoint"]
        heavy_bonds = []
        for b in state["original_qm_bond_distances"]:
            ids = b["atoms_1based"]
            if all(model.topology.atoms[i - 1].atomic_number != 1 for i in ids):
                before = distance(original, *ids)
                heavy_bonds.append(
                    dict(
                        b,
                        original_A=before,
                        change_A=b["distance_A"] - before,
                        distance_ratio=b["distance_A"] / before,
                    )
                )
        owner_changes = [
            dict(h, original_owner_1based=original_owners[h["hydrogen_1based"]])
            for h in state["proton_owners"]
            if h["nearest_qm_heavy_1based"] != original_owners[h["hydrogen_1based"]]
        ]
        unassigned = [
            h for h in state["proton_owners"] if h["classification"] != "covalent"
        ]
        geom = geometry(xyz, neighbors)
        if name == "A":
            first_a = records[0]["physical_energy_kcal_mol"]
            final_a = record["physical_energy_kcal_mol"]
        energy = {k: v for k, v in record.items() if "energy" in k}
        energy.update(
            physical_relative_ES_kcal_mol=record["physical_energy_kcal_mol"]
            - es_energy,
            biased_relative_ES_kcal_mol=record["total_energy_kcal_mol"] - es_energy,
            physical_relative_first_A_kcal_mol=record["physical_energy_kcal_mol"]
            - first_a,
            physical_relative_converged_A_kcal_mol=record["physical_energy_kcal_mol"]
            - final_a,
        )
        gates = {
            "process_complete": process["completed"],
            "native_SCC_passed": process["native_process_and_SCC_passed"],
            "biased_force_converged": accepted["biased_force_converged"],
            "all_54_H_covalently_assigned": len(state["proton_owners"]) == 54
            and not unassigned,
            "original_heavy_bonds_no_distance_above_1p8A": all(
                b["distance_A"] <= 1.8 for b in heavy_bonds
            ),
            "TI1_identity_validated": False,
        }
        points[name] = {
            "endpoint_path": str(paths[0]),
            "process": process,
            "acceptance": accepted,
            "geometry": geom,
            "energy": energy,
            "proton_owner_changes_vs_ES": owner_changes,
            "unassigned_or_shared_H": unassigned,
            "original_qm_heavy_bonds": heavy_bonds,
            "state": state,
            "gates": gates,
        }
        summaries.append(
            dict(
                point=name,
                attack_A=geom["distances_A"]["13916-19920"],
                CO_A=geom["distances_A"]["19920-19921"],
                CN_A=geom["distances_A"]["19920-19933"],
                BD_deg=geom["burgi_dunitz_OG1_C_O_deg"],
                carbonyl_height_A=geom["carbonyl_signed_neighbor_plane_height_A"],
                original_angle_sum_deg=geom["carbonyl_original_neighbor_angle_sum_deg"],
                max_biased_gradient=accepted["max_active_gradient_kcal_mol_A"],
                owner_changes=len(owner_changes),
                unassigned_shared=len(unassigned),
                **energy,
            )
        )
    result = {
        "scope": "Retained data only: no QM/MM calls; all IDs are 1-based Amber topology indices",
        "model": {"real_QM": 102, "real_QM_H": 54, "charge": -1, "active": 4662},
        "definitions": {
            "pyramidalization": "C19920 signed distance from plane through original heavy neighbors; signed volume det(neighbor vectors)/6; planar angle sum 360 degrees",
            "heavy_bond_screen": "Distance <=1.8 A is a retention screen, not bond-order validation",
            "proton_ownership": "Nearest real QM heavy; covalent d1<=1.25A and d2>1.35A; shared d2<=1.35A",
            "energy": "Physical native Eh*627.5094740631; biased adds actual harmonic bias; ES reference uses printed final Etot, precision 0.0001 kcal/mol; neither is a barrier or PMF",
        },
        "mapping": {
            str(i): {
                "name": model.topology.atoms[i - 1].name,
                "residue": model.topology.atoms[i - 1].residue.name,
                "amber_residue": model.topology.atoms[i - 1].residue.idx + 1,
            }
            for i in [
                1638,
                2732,
                13905,
                13910,
                13912,
                13916,
                14559,
                14560,
                14578,
                14579,
                14580,
                19920,
                19921,
                19933,
                20317,
            ]
        },
        "points": points,
        "source_sha256": HASHES,
    }
    for filename, payload in [("audit.json", result), ("sourcehashes.json", HASHES)]:
        with (OUT / filename).open("x") as stream:
            json.dump(payload, stream, indent=2, allow_nan=False)
            stream.write("\n")
    with (OUT / "summary.csv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    logging.info(
        "Audited four converged endpoints and original ES; no engine execution"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
