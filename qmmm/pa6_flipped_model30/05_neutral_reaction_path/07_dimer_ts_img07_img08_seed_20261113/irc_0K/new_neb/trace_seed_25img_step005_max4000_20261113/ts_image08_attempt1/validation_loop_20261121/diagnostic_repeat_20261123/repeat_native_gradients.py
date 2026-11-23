"""Audit same-point gradient repeatability and independently recheck endpoints."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import subprocess
import sys

import numpy as np

from nyl12_covalent.dftb_scc_validation import audit_dftb_scc
from nyl12_covalent.stationary_endpoints import classify_distances
from nyl12_covalent.stationary_native import NativeState, digest, write_json

LOGGER = logging.getLogger(__name__)
ROOT = Path("/data/jpc/covalent-enzyme-design")


def journal(phase: str, action: str, result: str) -> None:
    """Append the scientific action before/after a native evaluation."""
    subprocess.run([
        sys.executable, str(ROOT / "scripts/record_agent_action.py"),
        "--log", str(ROOT / "output_agent_actions/20261021.jsonl"),
        "--agent", "ts_same_point_diagnostic", "--phase", phase,
        "--action", action, "--permission-mode", "escalated",
        "--reason", "PI-authorized TSverificationdiagnosis; recoveredSCC semanticsfixed",
        "--result", result, "--rollback", "Preserve exact inputs/native evidence; no optimization",
    ], check=True, stdout=subprocess.DEVNULL)


def main() -> None:
    """Run exactly three TS repeats and two fresh endpoint gradients."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args.output.mkdir(exist_ok=False)
    source = NativeState(args.source)
    if "qmcharge=-2" not in source.mdin.replace(" ", ""):
        raise ValueError("Expected frozen QM charge -2")
    specification = {
        "gradient_calls_maximum": 5, "active_atoms": len(source.indices),
        "qmcharge": -2, "electronic_temperature_K": 0,
        "geometry_optimization": False, "full_hessian_recomputed": False,
        "formal_IRC": False, "verified_ts": False,
        "source_restart_sha256": digest(source.restart_path),
        "topology_sha256": digest(source.topology),
        "purpose": "Same-coordinate determinism and endpoint fresh-force/chemical assessment only",
    }
    write_json(args.output / "specification.json", specification)
    records = {}
    gradients = []
    tasks = [(f"TS_repeat_{i:02d}", source.restart_path) for i in range(1, 4)]
    tasks += [(sign, args.previous / "candidate_descents" / sign / "optimized.ncrst")
              for sign in ("minus", "plus")]
    for label, restart in tasks:
        if not restart.is_file():
            raise FileNotFoundError(restart)
        state = NativeState(args.source, restart=restart)
        action = f"Native analytic gradient {label}"
        journal("intent", action, f"RestartSHA256={digest(restart)}; unique freshSander process")
        record = {"restart": str(restart), "restart_sha256": digest(restart)}
        if label in {"minus", "plus"}:
            history = audit_dftb_scc((restart.parent / "minimize.mdout").read_text())
            record["original_descent_electronic_history"] = history
        try:
            energy, gradient = state.gradient(state.xyz, args.output / label)
        except Exception as error:
            LOGGER.exception("Native point failed: %s", label)
            record.update(process_accepted=False, error=str(error))
            records[label] = record
            write_json(args.output / "progress.json", records)
            journal("outcome", action, f"Failed: {error}")
            continue
        maximum = float(np.max(np.abs(gradient)) * 0.529177)
        record.update(process_accepted=True, energy_hartree=energy,
                      max_gradient_hartree_per_bohr=maximum,
                      stationary_force_gate_passed=maximum <= 1e-5,
                      coordinates_sha256=digest(restart),
                      gradient_sha256=digest(args.output / label / "gradient.npy"))
        if label.startswith("TS_"):
            gradients.append(gradient)
        else:
            pairs = {"OG_C": (13915, 19915), "OG_H": (13915, 13916),
                     "N_H": (19928, 13916), "C_N": (19915, 19928),
                     "Nterm_HG1": (13904, 13916),
                     "Nterm_H1": (13904, 13905), "Nterm_H2": (13904, 13906),
                     "carbonyl_C_O": (19915, 19916)}
            distances = {name: float(np.linalg.norm(state.xyz[i] - state.xyz[j]))
                         for name, (i, j) in pairs.items()}
            record.update(distances_angstrom=distances,
                          geometry_class=classify_distances(distances))
        records[label] = record
        write_json(args.output / "progress.json", records)
        journal("outcome", action, json.dumps(record))
    repeats = [records[f"TS_repeat_{i:02d}"] for i in range(1, 4)]
    all_accepted = all(r["process_accepted"] for r in records.values())
    result = {"records": records, "all_five_native_points_accepted": all_accepted,
              "verified_ts": False, "formal_IRC": False,
              "global_hessian_acceptance_changed": False}
    if len(gradients) == 3:
        differences = [gradients[i] - gradients[0] for i in (1, 2)]
        result["TS_repeatability"] = {
            "max_gradient_difference_hartree_per_angstrom": max(float(np.max(np.abs(g))) for g in differences),
            "max_RMS_gradient_difference_hartree_per_angstrom": max(float(np.sqrt(np.mean(g * g))) for g in differences),
            "energy_span_hartree": max(r["energy_hartree"] for r in repeats) - min(r["energy_hartree"] for r in repeats),
            "identical_gradient_hashes": len({r["gradient_sha256"] for r in repeats}) == 1,
            "interpretation": "Determinism test only; identical repeats do not bound systematic Hessian error",
        }
        ts_energy = repeats[0]["energy_hartree"]
        for sign in ("minus", "plus"):
            if records[sign]["process_accepted"]:
                records[sign]["relative_to_TS_kcal_mol"] = (records[sign]["energy_hartree"] - ts_energy) * 627.509474
    write_json(args.output / "diagnosis.json", result)
    LOGGER.info("Native points accepted: %s/5; no TS/IRC claim", sum(r["process_accepted"] for r in records.values()))
    if not all_accepted:
        raise RuntimeError("Some independent native points failed; retained all evidence")


if __name__ == "__main__":
    main()
