"""Prepare or run two further phase I increments on the original path."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = Path("/data/jpc/covalent-enzyme-design")
OUTPUT = ROOT / "output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2"
RUNNER_PATH = ROOT / "provenance/asp272_relay_20261025/run_relay_warning_only.py"
INSPECTION_CONFIG = ROOT / "configs/asp272_relay_inspection_20261025.json"
PLOT_SCRIPT = ROOT / "scripts/plot_reaction_distances.py"
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.path_inspection import publish_inspection  # noqa: E402

spec = importlib.util.spec_from_file_location("relay_warning_only", RUNNER_PATH)
runner = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(runner)

START = OUTPUT / "stage1_16/optimized.ncrst"
ORIGINAL_COMPLETED = [f"stage1_{i:02d}" for i in range(17)]


def schedule(attack_start: float, q1_start: float) -> list[dict[str, float | str]]:
    """Continue the original phase I line by its next two interval steps."""
    return [
        {
            "label": f"stage1_{index:02d}",
            "attack_target_A": 1.55 + ((index - 16) / 16) * (1.55 - attack_start),
            "q1_target_A": 0.8 + ((index - 16) / 16) * (0.8 - q1_start),
        }
        for index in (17, 18)
    ]


def validate_status(output: Path) -> tuple[dict, dict]:
    """Require the exact completed phase I status and restart handoff."""
    state = json.loads((output / "status.json").read_text())
    if state.get("phase1_complete") is not True or state.get("state") != "stopped":
        raise ValueError("Phase I must be complete and stopped")
    if state.get("completed_windows") != ORIGINAL_COMPLETED:
        raise ValueError("Expected exactly stage1_00 through stage1_16")
    acceptance = json.loads((output / "stage1_16/acceptance.json").read_text())
    process = json.loads((output / "stage1_16/process.json").read_text())
    if not (
        acceptance["native_converged"]
        and acceptance["finite_coordinates_and_fixed_box"]
    ):
        raise ValueError("stage1_16 did not pass native acceptance")
    if process.get("exit_status") != 0:
        raise ValueError("stage1_16 process did not complete successfully")
    result_flag = process["command"].index("-r") + 1
    if Path(process["command"][result_flag]).resolve() != START.resolve():
        raise ValueError("stage1_16 process does not identify required restart")
    return state, acceptance


def prepare(output: Path) -> None:
    """Preserve status and write the immutable two-window continuation plan."""
    state, _acceptance = validate_status(output)
    if (output / "status_before_phase1_extension_20261025.json").exists():
        raise FileExistsError("Phase I extension status backup already exists")
    plan_path = output / "phase1_extension_plan.json"
    if plan_path.exists():
        raise FileExistsError("Phase I extension plan already exists")
    if not START.is_file():
        raise FileNotFoundError(START)
    inputs = json.loads((output / "inputs.json").read_text())
    attack_start = float(inputs["initial_distances"]["attack"])
    q1_start = float(
        inputs["initial_distances"]["OG_H"] - inputs["initial_distances"]["Asp_H"]
    )
    process = json.loads((output / "stage1_16/process.json").read_text())
    if (
        process["attack_target_A"] != 1.55
        or process["proton_difference_target_A"] != 0.8
    ):
        raise ValueError("Unexpected original phase-I final targets")
    for item in schedule(attack_start, q1_start):
        if (output / str(item["label"])).exists():
            raise FileExistsError("Never overwrite an existing extension window")
    for path, expected in inputs["sha256"].items():
        if runner.old.digest(Path(path)) != expected:
            raise ValueError(f"Original scientific input changed: {path}")
    sources = [
        RUNNER_PATH,
        ROOT / "provenance/pa6_ti1_driving_20261024/run_path.py",
        output / "source/run_relay_warning_only.py",
        output / "source/continuation.mdin",
        INSPECTION_CONFIG,
        PLOT_SCRIPT,
        ROOT / "src/nyl12_covalent/path_inspection.py",
        START,
        output / "stage1_16/acceptance.json",
        output / "stage1_16/process.json",
        output / "inputs.json",
        ROOT / "src/nyl12_covalent/path_distances.py",
        ROOT / "src/nyl12_covalent/path_publish.py",
    ]
    hashes = {str(path): runner.old.digest(path) for path in sources}
    frozen = output / "source/phase1_extension_20261025"
    frozen.mkdir(exist_ok=False)
    for index, path in enumerate(sources):
        target = frozen / f"{index:02d}_{path.name}"
        shutil.copy2(path, target)
        if runner.old.digest(target) != hashes[str(path)]:
            raise ValueError("Frozen extension source copy differs")
    shutil.copy2(__file__, output / "source/extend_phase1_20261025.py")
    runner.old.write(output / "status_before_phase1_extension_20261025.json", state)
    runner.old.write(
        plan_path,
        {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "start_restart": str(START),
            "start_restart_sha256": hashes[str(START)],
            "q1_definition": "d(OG1-HG1)-d(Asp272 OD1-HG1)",
            "attack_start_A": attack_start,
            "q1_start_A": q1_start,
            "original_interval_denominator": 16,
            "window_count": 2,
            "schedule": schedule(attack_start, q1_start),
            "continuation_policy": "same stage 1 two-target line, restraints, spring, tolerance, reference, and warning-only tracking",
            "source_sha256": hashes,
            "frozen_source_directory": str(frozen),
            "wrapper_sha256": runner.old.digest(Path(__file__)),
        },
    )


def run(output: Path) -> None:
    """Run both windows and record the final chemical gate separately."""
    plan_path = output / "phase1_extension_plan.json"
    if (
        not plan_path.is_file()
        or not (output / "status_before_phase1_extension_20261025.json").is_file()
    ):
        raise ValueError("Run prepare mode first")
    plan = json.loads(plan_path.read_text())
    if runner.old.digest(Path(__file__)) != plan["wrapper_sha256"]:
        raise ValueError("Extension controller changed after preparation")
    if runner.old.digest(START) != plan["start_restart_sha256"]:
        raise ValueError("Required stage1_16 restart changed")
    for path, digest in plan["source_sha256"].items():
        if runner.old.digest(Path(path)) != digest:
            raise ValueError(f"Frozen source changed: {path}")
    state = json.loads((output / "status.json").read_text())
    if state.get("completed_windows") != ORIGINAL_COMPLETED:
        raise ValueError("Extension must start from preserved phase I status")
    inputs = json.loads((output / "inputs.json").read_text())
    restart = START
    state["state"] = "running"
    state.pop("reason", None)
    state["phase1_extension_complete"] = False
    state["phase1_extension_started_utc"] = datetime.now(timezone.utc).isoformat()
    runner.old.write(output / "status.json", state)
    try:
        publish_inspection(output, INSPECTION_CONFIG, PLOT_SCRIPT)
        for item in plan["schedule"]:
            label = str(item["label"])
            restart = runner.window(
                output,
                label,
                float(item["attack_target_A"]),
                float(item["q1_target_A"]),
                1,
                restart,
                inputs,
            )
            state["completed_windows"].append(label)
            runner.old.write(output / "status.json", state)
            publish_inspection(output, INSPECTION_CONFIG, PLOT_SCRIPT)
        try:
            runner.gate(output, "stage1_18", 1)
            state["phase1_extension_chemical_gate"] = "passed"
        except ValueError as error:
            state["phase1_extension_chemical_gate"] = "failed"
            state["phase1_extension_chemical_gate_reason"] = str(error)
        state.update(
            state="stopped",
            phase1_extension_complete=True,
            reason="Two original-interval phase-I windows completed; review energy and geometry",
        )
    except Exception as error:
        state.update(
            state="stopped", phase1_extension_complete=False, reason=str(error)
        )
        raise
    finally:
        state["updated_utc"] = datetime.now(timezone.utc).isoformat()
        runner.old.write(output / "status.json", state)
        publish_inspection(output, INSPECTION_CONFIG, PLOT_SCRIPT)


def main() -> None:
    """Prepare provenance or explicitly run the phase I extension."""
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    output = args.output.resolve()
    prepare(output) if args.prepare else run(output)


if __name__ == "__main__":
    main()
