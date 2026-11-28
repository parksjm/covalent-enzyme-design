#!/usr/bin/env python3
"""Review and adopt completed heating ramps without rerunning native MD."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil
import sys
import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_umbrella_campaign import digest, write_json  # noqa:E402
from nyl12_covalent.neb_umbrella_segments import _checkpoint  # noqa:E402
from nyl12_covalent.relay_umbrella_artifacts import verify_prepared  # noqa:E402
from nyl12_covalent.relay_umbrella_runtime import native_acceptance  # noqa:E402


def prepare(
    prepared: Path, original: Path, output: Path, workers: int, initial: int
) -> dict:
    """Freeze reviewed heat adoptions, preserving old plans and failure records."""
    base = prepared.parent
    cfg = verify_prepared(base)
    old = json.loads(original.read_text())
    if digest(original) != original.with_suffix(".sha256").read_text().strip():
        raise ValueError("Original execution plan changed")
    if digest(prepared) != old["prepared_sha256"]:
        raise ValueError("Scientific preparation changed")
    if len(cfg["windows"]) != 57 or not 1 <= initial <= workers <= 57:
        raise ValueError("Reviewed57state/budget requirement changed")
    output.mkdir(exist_ok=False)
    (output / "heat_acceptance").mkdir()
    (output / "original_failed_status").mkdir()
    top = parmed.load_file(str(base / cfg["topology"]))
    mask = np.zeros(len(top.atoms), dtype=bool)
    selected = list(AmberMask(top, cfg["system"]["active_mask"]).Selected())
    if len(selected) != 4165:
        raise ValueError("Mobile atom inventory changed")
    mask[selected] = True
    adopted = {}
    status_updates = []
    temperatures = []
    for row in cfg["windows"]:
        status = base / row["label"] / "window_status.json"
        if not status.exists():
            continue
        failure = json.loads(status.read_text())
        if (
            failure.get("status") != "failed"
            or failure.get("error")
            != "ValueError: Final temperature outside explicit gate"
        ):
            raise ValueError(f"Unexpected existing window status: {status}")
        spec = row["stages"][0]
        stage = base / Path(spec["input"]).parent
        if spec["kind"] != "heat" or spec["steps"] != 2000:
            raise ValueError("Unexpected heating protocol")
        if any(
            (base / Path(s["input"]).parent / "process.json").exists()
            for s in row["stages"][1:]
        ):
            raise ValueError("Recovery cannot adopt already-attempted later MD")
        process = json.loads((stage / "process.json").read_text())
        if process.get("exit_status") != 0 or process.get("ended_utc") is None:
            raise ValueError("Heating process was not successfully closed")
        identity = json.loads((stage / "input_identity.json").read_text())
        expected = {
            str(p): digest(p)
            for p in (
                original,
                base / cfg["topology"],
                base / row["seed"],
                stage / "run.mdin",
                stage / "umbrella.RST",
            )
        }
        if identity["sha256"] != expected:
            raise ValueError("Original heat input/predecessor lineage mismatch")
        text = (stage / "run.mdin").read_text()
        if "type='TEMP0', istep1=0, istep2=2000, value1=25.0, value2=300.0" not in text:
            raise ValueError("Reviewed1ps25→300K heating ramp changed")
        native = (stage / "run.mdout").read_text()
        dof = re.findall(
            r"TOTAL # of degrees of freedom \(RNDF\)\s*=\s*([\d.]+)", native
        )
        if dof != ["12495."]:
            raise ValueError("Native temperature DOF differs from3x4165")
        receipt = native_acceptance(stage, 2000, [0, 400])
        seed = Rst7.open(str(base / row["seed"]))
        checkpoint, _ = _checkpoint(
            stage / "restart.ncrst",
            np.asarray(seed.coordinates).reshape(-1, 3),
            mask,
            seed.box,
        )
        receipt.update(checkpoint)
        failed_copy = output / "original_failed_status" / f"{row['label']}.json"
        shutil.copy2(status, failed_copy)
        receipt.update(
            window=row["label"],
            stage=str(stage),
            temperature_acceptance_role="completed_nonstationary_heating_ramp",
            thermal_equilibrium_asserted=False,
            native_mobile_degrees_of_freedom=12495,
            original_input_plan=str(original),
            original_input_plan_sha256=digest(original),
            original_failed_status=str(failed_copy),
            original_failed_status_sha256=digest(failed_copy),
            recovery_authorization="Root2026-11-25: heatprocess[0,400]K only; continueunchanged5ps300Kequil with225–375Kinstantfinalgate",
            native_source_sha256={
                str(p): digest(p)
                for p in (
                    stage / "process.json",
                    stage / "input_identity.json",
                    stage / "run.mdout",
                    stage / "run.mdin",
                    stage / "umbrella.RST",
                    stage / "trajectory.nc",
                    stage / "restart.ncrst",
                    base / row["seed"],
                    base / cfg["topology"],
                )
            },
        )
        path = output / "heat_acceptance" / f"{row['label']}.json"
        write_json(path, receipt)
        adopted[row["label"]] = dict(receipt=str(path), sha256=digest(path))
        temperatures.append(receipt["final_temperature_K"])
        status_updates.append(
            (
                status,
                dict(
                    status="reviewed_heat_recovery_ready",
                    prior_failure=str(failed_copy),
                    heat_acceptance=str(path),
                ),
            )
        )
    if len(adopted) != 24:
        raise ValueError(f"Expected24reviewedcompletedheats, found{len(adopted)}")
    snapshot = output / "source_snapshot"
    shutil.copytree(
        ROOT / "src/nyl12_covalent",
        snapshot / "src/nyl12_covalent",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    (snapshot / "scripts").mkdir()
    for name in (
        "run_relay_umbrella.py",
        "guard_relay_umbrella.py",
        "inspect_relay_umbrella.py",
        "sync_relay_umbrella.py",
        "prepare_relay_heat_recovery.py",
    ):
        shutil.copy2(ROOT / "scripts" / name, snapshot / "scripts" / name)
    plan = dict(old)
    plan.update(
        workers=workers,
        initial_workers=initial,
        session_prefix="relay_us_heat_recovery_20261125",
        heating_temperature_gate_K=[0, 400],
        temperature_gate_K=[225, 375],
        adopted_heat=adopted,
        original_execution_plan=str(original),
        original_execution_plan_sha256=digest(original),
        source_sha256={str(p): digest(p) for p in snapshot.rglob("*.py")},
        root_review_authorized=True,
        root_review="Rootexplicit2026-11-25:57states,28→40→57workersbenchmark; noheatMDredo; unchanged5ps300Kequil",
    )
    path = output / "plan.json"
    write_json(path, plan)
    path.with_suffix(".sha256").write_text(digest(path) + "\n")
    for status, state in status_updates:
        write_json(status, state)
    summary = dict(
        ready=True,
        windows=57,
        completed_heat_adoptions=24,
        reused_heat_steps=48000,
        unstarted_windows=33,
        original_failures_preserved=True,
        heating_final_instantaneous_temperature_range_K=[
            min(temperatures),
            max(temperatures),
        ],
        heat_role="Completedramp,notthermalized",
        equilibration_temperature_gate_K=[225, 375],
        workers_budget=workers,
        initial_workers=initial,
        source_plan=str(path),
        source_plan_sha256=digest(path),
    )
    write_json(output / "recovery_preflight.json", summary)
    return summary


def main() -> None:
    """Require every recovery input and resource budget explicitly."""
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("prepared", "original-plan", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--initial-workers", type=int, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(
                args.prepared.resolve(strict=True),
                args.original_plan.resolve(strict=True),
                args.output.absolute(),
                args.workers,
                args.initial_workers,
            )
        )
    )


if __name__ == "__main__":
    main()
