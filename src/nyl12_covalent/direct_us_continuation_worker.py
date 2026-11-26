"""Dynamic native block-boundary execution for direct-US continuation."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from . import direct_us_coordinates, direct_us_continuation as continuation
from . import neb_umbrella_campaign as campaign
from . import neb_umbrella_diagnostics as diagnostics
from . import neb_umbrella_runtime as runtime
from . import neb_umbrella_segments as segments

LOGGER = logging.getLogger(__name__)
ATTEMPT_FILES = (
    "process.json",
    "input_identity.json",
    "restart.ncrst",
    "run.mdout",
    "trajectory.nc",
    "engine.log",
    "input.ncrst",
    "topology.prmtop",
)


def _worker_event(log_path: Path, event: str, **details: Any) -> None:
    """Append one durable worker claim/status event to its private log."""
    record = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "event": event,
        **details,
    }
    with log_path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def _stage_by_block(plan: dict[str, Any], label: str, block: int) -> dict[str, Any]:
    return next(
        record
        for record in plan["stages"]
        if record["window"] == label and record["block"] == block
    )


def _window_baseline(
    cfg: dict[str, Any], row: dict[str, Any], topology: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load original seed and fixed-boundary references for one window."""
    seed = Rst7.open(str(row["seed"]))
    baseline = np.asarray(seed.coordinates).reshape(-1, 3).copy()
    if baseline.shape != (len(topology.atoms), 3) or not np.isfinite(baseline).all():
        raise ValueError(f"Invalid original seed coordinates: {row['label']}")
    active = np.asarray(list(AmberMask(topology, cfg["active_mask"]).Selected()))
    if len(active) != 4165 or len(np.unique(active)) != 4165:
        raise ValueError("Frozen-boundary active atom count changed from 4165")
    movable = np.zeros(len(baseline), dtype=bool)
    movable[active] = True
    return baseline, movable, np.asarray(seed.box).copy()


def _validate_accepted_prefix(
    plan: dict[str, Any],
    cfg: dict[str, Any],
    row: dict[str, Any],
    topology: Any,
    next_block: int,
) -> Path:
    """Validate pinned block06 and each accepted continuation block before use."""
    label = row["label"]
    origin = plan["prepared_input_identity"][label]
    for key, path_key in (
        ("process_sha256", "stage"),
        ("acceptance_sha256", "stage"),
        ("restart_sha256", "restart"),
    ):
        path = Path(origin[path_key])
        if key == "process_sha256":
            path /= "process.json"
        elif key == "acceptance_sha256":
            path /= "acceptance.json"
        if campaign.digest(path) != origin[key]:
            raise ValueError(f"Pinned 30 ps {key} changed for {label}")
    baseline, movable, box = _window_baseline(cfg, row, topology)
    old_restart = Path(origin["restart"])
    segments._checkpoint(old_restart, baseline, movable, box)
    restart = old_restart
    for block in range(continuation.FIRST_BLOCK, next_block):
        record = _stage_by_block(plan, label, block)
        stage = Path(record["stage"])
        process_path = stage / "process.json"
        acceptance_path = stage / "acceptance.json"
        process = json.loads(process_path.read_text())
        acceptance = json.loads(acceptance_path.read_text())
        marker = json.loads((stage / "continuation_status.json").read_text())
        restart = stage / "restart.ncrst"
        if process.get("exit_status") != 0 or marker.get("status") != "accepted":
            raise ValueError(f"Native stage process failed: {stage}")
        if acceptance.get("checkpoint_sha256") != campaign.digest(restart):
            raise ValueError(f"Accepted checkpoint hash mismatch: {stage}")
        checkpoint, xyz = segments._checkpoint(restart, baseline, movable, box)
        if acceptance.get("completed_steps") != 10_000 or not acceptance.get(
            "finite_energy_temperature"
        ):
            raise ValueError(f"Accepted production record is invalid: {stage}")
        if (
            len(xyz) != origin["atom_count"]
            or checkpoint["checkpoint_sha256"] != acceptance["checkpoint_sha256"]
        ):
            raise ValueError(f"Checkpoint identity or atom count changed: {stage}")
    return restart


def _run_one_stage(
    plan_path: Path,
    plan: dict[str, Any],
    cfg: dict[str, Any],
    row: dict[str, Any],
    topology: Any,
    block: int,
    worker_id: int,
) -> None:
    """Run one sequential 5 ps production block and write its acceptance evidence."""
    plan, cfg_checked = continuation.verify_plan(plan_path)
    if cfg_checked["output"] != cfg["output"]:
        raise ValueError("Prepared campaign changed during worker lifetime")
    label = row["label"]
    record = _stage_by_block(plan, label, block)
    stage = Path(record["stage"])
    if campaign.digest(stage / "run.mdin") != record["mdin_sha256"]:
        raise ValueError(f"Pinned continuation run.mdin changed: {stage}")
    if campaign.digest(stage / "umbrella.RST") != record["restraint_sha256"]:
        raise ValueError(f"Pinned continuation restraint changed: {stage}")
    restart = _validate_accepted_prefix(plan, cfg, row, topology, block)
    if (stage / "continuation_status.json").exists():
        raise FileExistsError(f"Continuation stage has already been attempted: {stage}")
    if any((stage / name).exists() for name in ATTEMPT_FILES):
        raise FileExistsError(
            f"Continuation stage has native attempt artifacts: {stage}"
        )

    environment = {
        "AMBERHOME": str(Path(cfg["sander"]).parent.parent),
        "LD_LIBRARY_PATH": str(Path(cfg["sander"]).parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    identity = {
        "continuation_plan_sha256": campaign.digest(plan_path),
        "prepared_sha256": plan["prepared_sha256"],
        "original_block06_restart_sha256": plan["prepared_input_identity"][label][
            "restart_sha256"
        ],
        "source_sha256": plan["source_sha256"],
        "previous_restart_sha256": campaign.digest(restart),
        "mdin_sha256": record["mdin_sha256"],
        "restraint_sha256": record["restraint_sha256"],
        "environment": environment,
        "worker_id": worker_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    _worker_event(
        plan_path.parent / "workers" / f"worker_{worker_id:02d}.jsonl",
        "stage_claimed",
        window=label,
        block=block,
        stage=str(stage),
        previous_restart_sha256=identity["previous_restart_sha256"],
    )
    campaign.write_json(
        stage / "continuation_status.json", {"status": "running", **identity}
    )
    try:
        campaign.write_json(stage / "input_identity.json", identity)
        new_restart = runtime._run_stage(
            stage,
            Path(cfg["sander"]),
            Path(cfg["topology"]),
            restart,
            environment,
        )
        acceptance = runtime.stage_acceptance(
            stage, "production", cfg["temperature_gate_K"]
        )
        baseline, movable, box = _window_baseline(cfg, row, topology)
        checkpoint, xyz = segments._checkpoint(new_restart, baseline, movable, box)
        acceptance.update(checkpoint)
        acceptance["observables"] = runtime.observables(
            xyz, topology, cfg["roles_1based"]
        )
        campaign.write_json(stage / "acceptance.json", acceptance)
        status = {
            "status": "accepted",
            "block": block,
            "production_ps_per_window": block * continuation.BLOCK_PS,
            "continuation_ps_per_window": (block - 6) * continuation.BLOCK_PS,
            "checkpoint_sha256": campaign.digest(new_restart),
            "sampling_convergence": "Not established by completion",
            "scc_step_verification": plan["scc_step_verification"],
            "accepted_utc": datetime.now(timezone.utc).isoformat(),
        }
        campaign.write_json(stage / "continuation_status.json", {**identity, **status})
        _worker_event(
            plan_path.parent / "workers" / f"worker_{worker_id:02d}.jsonl",
            "stage_accepted",
            window=label,
            block=block,
            checkpoint_sha256=status["checkpoint_sha256"],
        )
        try:
            diagnostics.block_diagnostics(cfg, row, stage, topology)
            aug_dir = Path(plan["augmentation_root"]) / label / stage.name
            augmentation = direct_us_coordinates.augment_stage(stage, aug_dir)
            campaign.write_json(
                stage / "postprocess.json",
                {
                    "status": "completed",
                    "native_stage_accepted": True,
                    "checkpoint_sha256": status["checkpoint_sha256"],
                    "cross_bias_sha256": campaign.digest(
                        stage / "samples_and_cross_bias.npz"
                    ),
                    "augmentation": str(aug_dir.resolve()),
                    "augmentation_csv_sha256": augmentation["csv_sha256"],
                    "completed_utc": datetime.now(timezone.utc).isoformat(),
                },
            )
            _worker_event(
                plan_path.parent / "workers" / f"worker_{worker_id:02d}.jsonl",
                "stage_postprocess_completed",
                window=label,
                block=block,
            )
        except Exception as error:
            campaign.write_json(
                stage / "postprocess.json",
                {
                    "status": "failed",
                    "native_stage_accepted": True,
                    "checkpoint_sha256": status["checkpoint_sha256"],
                    "error": f"{type(error).__name__}: {error}",
                    "failed_utc": datetime.now(timezone.utc).isoformat(),
                },
            )
            _worker_event(
                plan_path.parent / "workers" / f"worker_{worker_id:02d}.jsonl",
                "stage_postprocess_failed_native_accepted",
                window=label,
                block=block,
                error=f"{type(error).__name__}: {error}",
            )
    except Exception as error:
        failure = {
            "status": "failed",
            "block": block,
            "error": f"{type(error).__name__}: {error}",
            "failed_utc": datetime.now(timezone.utc).isoformat(),
        }
        campaign.write_json(stage / "continuation_status.json", {**identity, **failure})
        _worker_event(
            plan_path.parent / "workers" / f"worker_{worker_id:02d}.jsonl",
            "stage_failed",
            window=label,
            block=block,
            error=failure["error"],
        )
        raise
