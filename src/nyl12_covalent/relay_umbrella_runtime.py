"""Run immutable relay umbrellas through explicit production block boundaries."""

from __future__ import annotations

import json
from datetime import datetime, timezone
import math
from pathlib import Path
import re
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from .amber_md_records import instantaneous_md_text
from .dftb_scc_validation import audit_dftb_scc
from .relay_umbrella_recovery import heat_adoption
from .neb_umbrella_campaign import digest, write_json
from .neb_umbrella_runtime import _run_stage
from .neb_umbrella_segments import _checkpoint


def verify_plan(path: Path) -> tuple[dict[str, Any], dict[str, Any], Path]:
    """Require unchanged plan, scientific files and host runtime."""
    if digest(path) != path.with_suffix(".sha256").read_text().strip():
        raise ValueError("Execution plan changed")
    plan = json.loads(path.read_text())
    prepared = Path(plan["prepared"])
    if digest(prepared) != plan["prepared_sha256"]:
        raise ValueError("Scientific preparation changed")
    cfg = json.loads(prepared.read_text())
    base = prepared.parent
    for relative, expected in cfg["frozen_files_sha256"].items():
        if digest(base / relative) != expected:
            raise ValueError(f"Frozen scientific input changed: {relative}")
    for field in ("runtime_sha256", "source_sha256"):
        for name, expected in plan[field].items():
            if digest(Path(name)) != expected:
                raise ValueError(f"Pinned file changed: {name}")
    return plan, cfg, base


def native_acceptance(stage: Path, steps: int, gate: list[float]) -> dict[str, Any]:
    """Require finite completed MD and resolved SCC for every integration step."""
    text = (stage / "run.mdout").read_text()
    electronic_history = audit_dftb_scc(text)
    if not electronic_history["verified_scc_history"]:
        raise ValueError("Unresolved or unverified electronic history")
    body = instantaneous_md_text(text)
    if re.search(
        r"\b(?:NaN|Infinity)\b|SHAKE cannot be accomplished|FATAL|ERROR:", body, re.I
    ):
        raise ValueError("Native numerical/engine failure")
    rows = re.findall(
        r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([-+\d.EeDd]+)\s+TEMP\(K\)\s*=\s*([-+\d.EeDd]+)",
        body,
    )
    energies = [
        float(v.replace("D", "E")) for v in re.findall(r"EPtot\s*=\s*(\S+)", body)
    ]
    temps = [float(r[2].replace("D", "E")) for r in rows]
    if not rows or int(rows[-1][0]) != steps:
        raise ValueError("Native MD did not complete expected step count")
    if not energies or not all(math.isfinite(v) for v in [*energies, *temps]):
        raise ValueError("Missing/nonfinite energy or temperature")
    if not gate[0] <= temps[-1] <= gate[1]:
        raise ValueError("Final temperature outside explicit gate")
    converged = {
        int(v) for v in re.findall(r"SCC-DFTB for step\s+(\d+) converged", body)
    }
    missing = set(range(1, steps + 1)) - converged
    if missing:
        raise ValueError(f"Unresolved SCC integration steps: {sorted(missing)[:12]}")
    return dict(
        accepted=True,
        completed_steps=steps,
        final_temperature_K=temps[-1],
        temperature_gate_K=gate,
        finite_energy_temperature=True,
        every_integration_step_SCC_converged=True,
        electronic_history=electronic_history,
    )


def ensure_diagnostics(
    cfg: dict[str, Any], row: dict[str, Any], stage: Path, topology: Any
) -> bool:
    """Retain native acceptance while recording or retrying derived-only analysis."""
    outcome = stage / "analysis_outcome.json"
    if outcome.exists() and json.loads(outcome.read_text())["status"] == "completed":
        return True
    if outcome.exists():
        archive = stage / (
            "failed_analysis_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        )
        archive.mkdir()
        for name in (
            "analysis_outcome.json",
            "analysis_request.json",
            "raw_features.csv",
            "samples_and_cross_bias.npz",
            "observables",
            "block_acceptance.json",
            "analysis_acceptance.json",
        ):
            path = stage / name
            if path.exists():
                path.rename(archive / name)
    hashes = {
        str(stage / name): digest(stage / name)
        for name in (
            "acceptance.json",
            "process.json",
            "restart.ncrst",
            "trajectory.nc",
        )
    }
    try:
        from .relay_umbrella_diagnostics import block_diagnostics

        block_diagnostics(cfg, row, stage, topology)
        for name in (
            "analysis_request.json",
            "raw_features.csv",
            "samples_and_cross_bias.npz",
        ):
            if not (stage / name).is_file():
                raise ValueError(f"Missing derived artifact: {name}")
        write_json(outcome, dict(status="completed", native_sha256=hashes))
        return True
    except Exception as error:
        write_json(
            outcome,
            dict(
                status="failed",
                error=f"{type(error).__name__}: {error}",
                native_sha256=hashes,
            ),
        )
        return False


def run_segment(plan_path: Path, label: str, limit: int) -> None:
    """Run a fresh window or extend an accepted prefix to a block limit."""
    plan, cfg, base = verify_plan(plan_path)
    if not 1 <= limit <= 6:
        raise ValueError("Production block limit must be 1..6")
    row = next(w for w in cfg["windows"] if w["label"] == label)
    folder = base / label
    topology_path = base / "inputs/system.prmtop"
    topology = parmed.load_file(str(topology_path))
    initial = Rst7.open(str(base / row["seed"]))
    baseline = np.asarray(initial.coordinates).reshape(-1, 3)
    active = np.asarray(
        list(AmberMask(topology, cfg["system"]["active_mask"]).Selected())
    )
    if len(active) != 4165 or len(baseline) != 104922:
        raise ValueError("Audited atom inventory changed")
    movable = np.zeros(len(baseline), dtype=bool)
    movable[active] = True
    restart = base / row["seed"]
    count = 0
    production_steps = 0
    status = folder / "window_status.json"
    if status.exists() and json.loads(status.read_text())["status"] == "failed":
        raise ValueError("Failed attempt requires explicit reviewed recovery")
    try:
        for spec in row["stages"]:
            if spec["kind"] == "production" and count >= limit:
                break
            stage = base / Path(spec["input"]).parent
            acceptance_path = stage / "acceptance.json"
            input_plan = plan_path
            if spec["kind"] == "heat":
                adoption = heat_adoption(plan, label, stage)
                if adoption is not None:
                    acceptance_path, input_plan = adoption
            gate = (
                plan["heating_temperature_gate_K"]
                if spec["kind"] == "heat"
                else plan["temperature_gate_K"]
            )
            if acceptance_path.exists():
                receipt = json.loads(acceptance_path.read_text())
                identity = json.loads((stage / "input_identity.json").read_text())
                expected = {
                    str(p): digest(p)
                    for p in (
                        input_plan,
                        topology_path,
                        restart,
                        stage / "run.mdin",
                        stage / "umbrella.RST",
                    )
                }
                if identity["sha256"] != expected:
                    raise ValueError("Accepted input/predecessor lineage changed")
                if json.loads((stage / "process.json").read_text())["exit_status"] != 0:
                    raise ValueError("Accepted stage lacks successful process")
                restart = stage / "restart.ncrst"
                check, _ = _checkpoint(restart, baseline, movable, initial.box)
                if check["checkpoint_sha256"] != receipt["checkpoint_sha256"]:
                    raise ValueError("Accepted checkpoint changed")
                native_acceptance(stage, spec["steps"], gate)

            else:
                if any(
                    (stage / n).exists()
                    for n in (
                        "process.json",
                        "input_identity.json",
                        "run.mdout",
                        "restart.ncrst",
                        "input.ncrst",
                    )
                ):
                    raise FileExistsError(f"Attempt already exists: {stage}")
                write_json(
                    status,
                    dict(
                        status="running", stage=spec["label"], production_blocks=count
                    ),
                )
                write_json(
                    stage / "input_identity.json",
                    dict(
                        sha256={
                            str(p): digest(p)
                            for p in (
                                plan_path,
                                topology_path,
                                restart,
                                stage / "run.mdin",
                                stage / "umbrella.RST",
                            )
                        },
                        environment=plan["environment"],
                        seed=spec["seed"],
                    ),
                )
                restart = _run_stage(
                    stage,
                    Path(plan["sander"]),
                    topology_path,
                    restart,
                    plan["environment"],
                )
                accepted = native_acceptance(stage, spec["steps"], gate)
                checkpoint, _ = _checkpoint(restart, baseline, movable, initial.box)
                accepted.update(checkpoint)
                accepted["temperature_acceptance_role"] = (
                    "completed_nonstationary_heating_ramp"
                    if spec["kind"] == "heat"
                    else "thermalized_stage_final_instantaneous_temperature"
                )
                write_json(acceptance_path, accepted)
            if spec["kind"] == "production":
                ensure_diagnostics(cfg, row, stage, topology)
            count += spec["kind"] == "production"
            if spec["kind"] == "production":
                production_steps += spec["steps"]
        if count != limit:
            raise ValueError("Prepared stages do not reach requested production limit")
        write_json(
            status,
            dict(
                status="completed" if count == 6 else "paused_at_block_boundary",
                production_blocks=count,
                production_ps=production_steps * 0.0005,
                sampling_converged=False,
            ),
        )
        write_json(
            folder / f"segment_receipt_limit_{limit:02d}.json",
            dict(
                status="accepted",
                production_blocks=count,
                checkpoint=str(restart),
                checkpoint_sha256=digest(restart),
                sampling_converged=False,
            ),
        )
    except Exception as error:
        write_json(
            status,
            dict(
                status="failed",
                production_blocks=count,
                error=f"{type(error).__name__}: {error}",
            ),
        )
        raise
