"""Bounded same-state native probes of a completed candidate Hessian."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import logging
import os
from pathlib import Path

import numpy as np

from nyl12_covalent.stationary_candidate_branch import candidate_branch
from nyl12_covalent.stationary_native import NativeState, digest, write_json
from nyl12_covalent.stationary_probe_math import curvature_record

LOGGER = logging.getLogger(__name__)


def probe(args: tuple[Path, Path, str, int, float, int]) -> dict[str, object]:
    """Evaluate one signed displacement with frozen native-frame/SCC checks.

    Args:
        args: Source, output, probe type, coordinate/mode index, step, sign.

    Returns:
        Energy and exact retained-gradient evidence location.
    """
    source, output, kind, index, step, sign = args
    os.nice(max(0, 10 - os.getpriority(os.PRIO_PROCESS, 0)))
    state = NativeState(source)
    xyz = state.xyz.copy()
    if kind == "mode":
        vector = np.load(source / "ts_hessian/lowest_eigenvectors.npy")[:, index]
        xyz[state.active] += sign * step * (vector * state.inv_sqrt).reshape(-1, 3)
    else:
        xyz = state.displaced(index, sign * step)
    evidence = output / "probes" / f"{kind}_{index}_{step:g}_{sign:+d}"
    energy, _ = state.gradient(xyz, evidence)
    return {
        "kind": kind,
        "index": index,
        "step": step,
        "sign": sign,
        "energy": energy,
        "evidence": str(evidence),
    }


def run(source: Path, output: Path) -> None:
    """Run bounded diagnostics and conditionally queue candidate descents.

    Args:
        source: Frozen campaign with complete Hessian and native inputs.
        output: New independent diagnostic campaign directory.
    """
    state = NativeState(source)
    if "qmcharge=-2" not in state.mdin.replace(" ", ""):
        raise ValueError("Explicit frozen QM charge -2 required")
    modes = np.load(source / "ts_hessian/lowest_eigenvectors.npy")
    eigenvalues = np.load(source / "ts_hessian/eigenvalues.npy")
    write_json(
        output / "scientific_specification.json",
        {
            "mode": "candidate_diagnostic_only",
            "verified_ts": False,
            "source_campaign": str(source),
            "hamiltonian_changed": False,
            "active_atoms": len(state.indices),
            "qmcharge": -2,
            "protonation": "Unchanged frozen neutral Thr N / ASP272 / ASP274",
            "workers": 2,
            "threads_per_worker": 1,
            "maximum_probe_gradients": 27,
            "mode_steps_mass_weighted": [0.01, 0.025, 0.05],
            "pair_refined_steps_angstrom": [0.000125, 0.0000625],
            "acceptance": "Report independent curvatures; do not replace original full-Hessian acceptance.",
            "descents": "Only after negative curvatures and both 0.05 seeds >1e-7Eh downhill; candidate connectivity, not IRC or verified TS.",
        },
    )
    center_energy, _ = state.gradient(state.xyz, output / "center")
    tasks = [
        (source, output, "mode", mode, step, sign)
        for mode in range(3)
        for step in (0.01, 0.025, 0.05)
        for sign in (-1, 1)
    ]
    tasks += [
        (source, output, "column", col, step, sign)
        for col in (6586, 6587)
        for step in (0.000125, 0.0000625)
        for sign in (-1, 1)
    ]
    records = []
    failures = []
    with ProcessPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(probe, task): task for task in tasks}
        for future in as_completed(futures):
            try:
                records.append(future.result())
            except Exception as exc:
                failures.append(
                    {"task": [str(x) for x in futures[future]], "error": str(exc)}
                )
            write_json(
                output / "progress.json",
                {
                    "completed": len(records),
                    "failed": failures,
                    "total": len(tasks),
                    "center_energy": center_energy,
                },
            )
    write_json(
        output / "probe_receipts.json", {"records": records, "failures": failures}
    )
    if failures:
        raise RuntimeError("Native probes failed; preserved independent successes")
    table = {(x["kind"], x["index"], x["step"], x["sign"]): x for x in records}
    curvatures = []
    for mode in range(3):
        direction = (modes[:, mode] * state.inv_sqrt).reshape(-1, 3)
        for step in (0.01, 0.025, 0.05):
            plus = table[("mode", mode, step, 1)]
            minus = table[("mode", mode, step, -1)]
            record = curvature_record(
                np.load(Path(plus["evidence"]) / "gradient.npy"),
                np.load(Path(minus["evidence"]) / "gradient.npy"),
                direction,
                step,
                plus["energy"],
                minus["energy"],
                center_energy,
                float(eigenvalues[mode]),
            )
            curvatures.append(dict(record, mode=mode))
    pair_records = []
    for step in (0.00025, 0.000125, 0.0000625):
        columns = []
        for col in (6586, 6587):
            if step == 0.00025:
                folder = source / "ts_hessian/columns" / f"column_{col:05d}"
                receipt = json.loads((folder / "column.json").read_text())
                if (
                    receipt["sha256"] != digest(folder / "column.npy")
                    or receipt["step"] != step
                ):
                    raise ValueError("Original pair-column hash/step mismatch")
                column = np.load(folder / "column.npy")
            else:
                plus = Path(table[("column", col, step, 1)]["evidence"])
                minus = Path(table[("column", col, step, -1)]["evidence"])
                column = (
                    (
                        (
                            np.load(plus / "gradient.npy")
                            - np.load(minus / "gradient.npy")
                        )
                        / (2 * step)
                    ).reshape(-1)
                    * state.inv_sqrt
                    * state.inv_sqrt[col]
                )
            columns.append(column)
        pair_records.append(
            {
                "step_angstrom": step,
                "H_yz": float(columns[1][6586]),
                "H_zy": float(columns[0][6587]),
                "half_difference": float((columns[1][6586] - columns[0][6587]) / 2),
            }
        )
    negative = [r for r in curvatures if r["mode"] == 0]
    seed = next(r for r in negative if r["step_mass_weighted"] == 0.05)
    ready = all(r["gradient_curvature"] < 0 for r in negative) and all(
        seed[f"{sign}_energy_difference_hartree"] < -1e-7 for sign in ("minus", "plus")
    )
    write_json(
        output / "curvature_diagnostic.json",
        {
            "curvatures": curvatures,
            "dominant_asymmetry_pair": pair_records,
            "candidate_descents_ready": ready,
            "verified_ts": False,
            "global_inertia_acceptance_unchanged": True,
            "energy_precision": "Energy second differences may lose significance; gradient directional curvature is primary.",
        },
    )
    if ready:
        with ProcessPoolExecutor(max_workers=2) as executor:
            branches = {
                executor.submit(
                    candidate_branch, (source, output, sign, center_energy)
                ): sign
                for sign in ("minus", "plus")
            }
            results = {}
            for future in as_completed(branches):
                sign = branches[future]
                try:
                    results[sign] = future.result()
                except Exception as exc:
                    results[sign] = {"error": str(exc), "verified_ts": False}
                write_json(
                    output / "candidate_connectivity_diagnostic.json",
                    {"branches": results, "verified_ts": False, "formal_IRC": False},
                )


def main() -> None:
    """Parse explicit frozen-source and output paths for a diagnostic job."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.source, args.output)


if __name__ == "__main__":
    main()
