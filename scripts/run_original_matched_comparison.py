"""Coordinate a frozen six-case original-backbone comparison on Violet."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import time
from typing import Any

SEEDS = (20261126, 20261127)
QUEUES = (
    (("WT", SEEDS[0]), ("relay_prot101_06", SEEDS[0]), ("relay_prot101_03", SEEDS[1])),
    (("WT", SEEDS[1]), ("relay_prot101_03", SEEDS[0]), ("relay_prot101_06", SEEDS[1])),
)


def digest(path: Path) -> str:
    """Return a file's SHA256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    """Write a new immutable JSON artifact, refusing overwrites."""
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def verify_pins(pins: dict[str, str]) -> None:
    """Require unchanged absolute source paths before execution."""
    if not pins:
        raise ValueError("Explicit source pins required")
    for name, expected in pins.items():
        path = Path(name)
        if not path.is_absolute() or digest(path) != expected:
            raise ValueError(f"Frozen source changed: {name}")


def journal(config: dict[str, Any], phase: str, action: str) -> None:
    """Append required controller intent or outcome to the active journal."""
    root = Path(config["repository_root"])
    subprocess.run(
        [
            config["python"],
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "original_matched_controller",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Frozen matched comparison; independent failures retained",
            "--rollback",
            "Preserve every input and attempt; no process signaling or deletion",
        ],
        check=True,
    )


def wait_searches(config: dict[str, Any]) -> None:
    """Wait for both pre-existing search jobs, failing on nonzero exits."""
    root = Path(config["campaign_root"])
    while True:
        done = True
        for name in ("relay_prot101_06", "relay_prot101_03"):
            path = root / "jobs" / f"msa_{name}" / "status.json"
            state = json.loads(path.read_text())
            code = state["exit_status"]
            if code is not None and code != 0:
                raise RuntimeError(f"Search failed: {name}, exit {code}")
            done &= code is not None
        if done:
            return
        logging.info("Waiting for both existing precursor searches")
        time.sleep(30)


def invoke(config: dict[str, Any], argv: list[str], name: str) -> int:
    """Execute one logged stage and retain its command, output and exit status."""
    verify_pins(config["source_sha256"])
    root = Path(config["campaign_root"]) / "controller"
    journal(config, "intent", f"Execute {name}: {json.dumps(argv)}")
    env = {
        **os.environ,
        "NYL12_PROJECT_ROOT": config["repository_root"],
        "PYTHONPATH": str(Path(config["source_snapshot"]) / "src"),
        "OMP_NUM_THREADS": "2",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "2",
    }
    with (root / f"{name}.log").open("x") as log:
        result = subprocess.run(
            argv, stdout=log, stderr=subprocess.STDOUT, env=env, check=False
        )
    write_json(root / f"{name}.json", {"argv": argv, "exit_status": result.returncode})
    journal(config, "outcome", f"Stage {name} exit {result.returncode}")
    return result.returncode


def case_name(candidate: str, seed: int) -> str:
    """Return the stable native case identity."""
    return f"{candidate}_seed_{seed}"


def make_plans(config: dict[str, Any]) -> tuple[list[Path], list[Path]]:
    """Freeze two GPU-specific prediction and complete six-case Amber plans."""
    root = Path(config["campaign_root"])
    pred = json.loads(Path(config["prediction_plan"]).read_text())
    template = json.loads(Path(config["amber_template"]).read_text())
    candidates = {c["candidate"]: c for c in pred["candidates"]}
    if set(candidates) != {"WT", "relay_prot101_03", "relay_prot101_06"}:
        raise ValueError("Require exactly WT, relay_prot101_03 and relay_prot101_06")
    if pred["seeds"] != list(SEEDS):
        raise ValueError("Require fixed matched seeds")
    gpus = [candidates[n]["gpu_uuid"] for n in ("relay_prot101_06", "relay_prot101_03")]
    if len(set(gpus)) != 2:
        raise ValueError("Require two distinct reserved GPU UUIDs")
    template["cases"] = []
    for candidate in candidates:
        for seed in SEEDS:
            name = case_name(candidate, seed)
            work = root / "prediction/cases" / name
            case = {
                "candidate": candidate,
                "seed": seed,
                "case": name,
                "case_input": str(work / "case.json"),
                "native_receipt": str(work / "native_prediction.json"),
            }
            template["cases"].append(case)
            template["input_sha256"][case["case_input"]] = digest(
                Path(case["case_input"])
            )
    prediction_paths, amber_paths = [], []
    for worker, gpu in enumerate(gpus):
        prediction, amber = deepcopy(pred), deepcopy(template)
        for item in prediction["candidates"]:
            item["gpu_uuid"] = gpu
        amber["relaxation_config"]["gpu"] = gpu
        verify_pins(amber["input_sha256"])
        pp = root / "controller" / f"worker{worker}_prediction_plan.json"
        ap = root / "controller" / f"worker{worker}_amber_plan.json"
        write_json(pp, prediction)
        write_json(ap, amber)
        prediction_paths.append(pp)
        amber_paths.append(ap)
    return prediction_paths, amber_paths


def run(config: dict[str, Any]) -> dict[str, Any]:
    """Prepare and execute all native cases before matched Amber follow-through.

    Args:
        config: Absolute frozen paths and source pins for this new attempt.

    Returns:
        Six-case process and scientific outcome census.
    """
    for key in (
        "repository_root",
        "campaign_root",
        "python",
        "source_snapshot",
        "prediction_plan",
        "amber_template",
    ):
        if not Path(config[key]).is_absolute():
            raise ValueError(f"Require absolute {key}")
    verify_pins(config["source_sha256"])
    journal(config, "intent", "Initialize controller and wait for existing searches")
    root = Path(config["campaign_root"])
    (root / "controller").mkdir(exist_ok=False)
    write_json(root / "controller/config.json", config)
    wait_searches(config)
    script = Path(config["source_snapshot"]) / "scripts"
    prediction = str(root / "prediction")
    command = [config["python"], str(script / "run_matched_original_boltz.py")]
    for candidate in ("WT", "relay_prot101_06", "relay_prot101_03"):
        if invoke(
            config,
            command
            + [
                config["prediction_plan"],
                prediction,
                "--stage",
                "prepare",
                "--candidate",
                candidate,
            ],
            f"prepare_{candidate}",
        ):
            raise RuntimeError(f"Preparation failed: {candidate}")
    prediction_plans, amber_plans = make_plans(config)
    journal(config, "outcome", "Prepared six immutable cases and four GPU worker plans")

    def native_worker(worker: int) -> dict[str, int]:
        codes = {}
        for candidate, seed in QUEUES[worker]:
            name = case_name(candidate, seed)
            codes[name] = invoke(
                config,
                command
                + [
                    str(prediction_plans[worker]),
                    prediction,
                    "--stage",
                    "predict",
                    "--candidate",
                    candidate,
                    "--seed",
                    str(seed),
                ],
                f"predict_{name}",
            )
        return codes

    with ThreadPoolExecutor(max_workers=2) as pool:
        native_results = list(pool.map(native_worker, range(2)))
    native = {k: v for result in native_results for k, v in result.items()}

    def amber_worker(worker: int) -> dict[str, Any]:
        results = {}
        for candidate, seed in QUEUES[worker]:
            name = case_name(candidate, seed)
            if native[name]:
                results[name] = {
                    "status": "skipped_native_failure",
                    "native_exit_status": native[name],
                }
                write_json(
                    root / "controller" / f"skip_amber_{name}.json", results[name]
                )
                continue
            code = invoke(
                config,
                [
                    config["python"],
                    str(script / "run_matched_original_amber.py"),
                    "--plan",
                    str(amber_plans[worker]),
                    "--output",
                    str(root / "relaxation"),
                    "--case",
                    name,
                ],
                f"amber_{name}",
            )
            path = root / "relaxation" / f"failure_census_{name}.json"
            if path.is_file():
                results[name] = {
                    **json.loads(path.read_text())["cases"][name],
                    "amber_exit_status": code,
                }
            else:
                results[name] = {
                    "status": "failed_missing_amber_census",
                    "amber_exit_status": code,
                    "missing_path": str(path),
                }
        return results

    with ThreadPoolExecutor(max_workers=2) as pool:
        amber_results = list(pool.map(amber_worker, range(2)))
    cases = {k: v for result in amber_results for k, v in result.items()}
    report = {
        "cases": cases,
        "native_exit_status": native,
        "failure_count": sum(c["status"] != "completed" for c in cases.values()),
    }
    write_json(root / "controller/summary.json", report)
    fields = (
        "case",
        "status",
        "positive_ab_before",
        "positive_ab_after",
        "force_converged",
        "final_rms_gradient",
    )
    with (root / "controller/comparison.csv").open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name, result in cases.items():
            writer.writerow(
                {k: name if k == "case" else result.get(k, "") for k in fields}
            )
    journal(
        config,
        "outcome",
        f"Six-case census retained; failures={report['failure_count']}",
    )
    return report


def main() -> None:
    """Execute one explicitly configured persistent controller attempt."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    report = run(json.loads(args.config.read_text()))
    if report["failure_count"]:
        raise SystemExit(1)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
