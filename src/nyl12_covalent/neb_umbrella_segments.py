"""Resume native umbrella windows at accepted production block boundaries."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
import math
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from . import neb_umbrella_campaign as campaign
from . import neb_umbrella_diagnostics as diagnostics
from . import neb_umbrella_runtime as runtime

LOGGER = logging.getLogger(__name__)


def _checkpoint(
    path: Path, baseline: np.ndarray, movable: np.ndarray, box: np.ndarray
) -> tuple[dict[str, Any], np.ndarray]:
    """Validate one native checkpoint against the immutable seed boundary."""
    restart = Rst7.open(str(path))
    xyz = np.asarray(restart.coordinates).reshape(-1, 3)
    velocities = np.asarray(restart.velocities)
    if (
        xyz.shape != baseline.shape
        or not np.isfinite(xyz).all()
        or not restart.hasvels
        or not np.isfinite(velocities).all()
    ):
        raise ValueError(f"Invalid/nonfinite checkpoint coordinates/velocities: {path}")
    if not np.allclose(restart.box, box, rtol=0, atol=1e-5):
        raise ValueError(f"Fixed periodic box changed: {path}")
    drift = float(np.max(np.abs(xyz[~movable] - baseline[~movable])))
    if not math.isfinite(drift) or drift > 1e-6:
        raise ValueError(
            f"Frozen exterior moved beyond native output precision: {path}"
        )
    return {
        "checkpoint_sha256": campaign.digest(path),
        "frozen_outer_max_displacement_A": drift,
    }, xyz


def _adoption_matches(
    adoption_path: Path, stage: Path, process_path: Path
) -> dict[str, Any]:
    """Verify external scheduler evidence without changing its native process record."""
    original_bytes = process_path.read_bytes()
    process = json.loads(original_bytes)
    adoption_bytes = adoption_path.read_bytes()
    adoption = json.loads(adoption_bytes)
    if (
        adoption.get("stage") != str(stage.resolve())
        or adoption.get("native_exit_status") != 0
        or adoption.get("native_outputs_complete") is not True
        or adoption.get("process_json_sha256") != campaign.digest(process_path)
    ):
        raise ValueError("Scheduler adoption does not match attempted native stage")
    identity = adoption.get("native_process_identity")
    if not isinstance(identity, dict) or not identity:
        raise ValueError("Scheduler adoption lacks native process identity")
    if process.get("command") != identity.get("argv") or process.get(
        "cwd"
    ) != identity.get("cwd"):
        raise ValueError(
            "Scheduler process identity differs from native process record"
        )
    if "exit_status" in process and process["exit_status"] != 0:
        raise ValueError("Native process record reports nonzero exit status")
    observed = adoption.get("observed_utc")
    if not isinstance(observed, str) or not observed:
        raise ValueError("Scheduler adoption lacks observed_utc")
    raw_status = adoption.get("raw_native_wait_status")
    if raw_status not in (0, "0"):
        raise ValueError("Scheduler adoption lacks successful native wait status")
    archive = stage / "process_before_controller_handoff.json"
    with archive.open("xb") as stream:
        stream.write(original_bytes)
    process.update(
        exit_status=0,
        ended_utc=observed,
        completion_source="verified_linux_zombie_exit_status",
        adoption_receipt_sha256=campaign.digest(adoption_path),
    )
    campaign.write_json(process_path, process)
    return process


def _accepted_prefix(
    cfg: dict[str, Any],
    row: dict[str, Any],
    folder: Path,
    production_limit: int,
    baseline: np.ndarray,
    movable: np.ndarray,
    box: np.ndarray,
) -> tuple[list[str], Path, int, list[dict[str, Any]]]:
    """Validate prior acceptance records and return the first unaccepted stage."""
    completed: list[str] = []
    records: list[dict[str, Any]] = []
    restart = Path(row["seed"])
    production = 0
    gap = False
    attempted_gap = False
    for spec in row["stages"]:
        stage = folder / spec["label"]
        process_path = stage / "process.json"
        acceptance_path = stage / "acceptance.json"
        if acceptance_path.exists():
            if gap:
                raise ValueError("Accepted umbrella stages do not form a prefix")
            if not process_path.is_file():
                raise ValueError(f"Accepted stage lacks process record: {stage}")
            process = json.loads(process_path.read_text())
            if process.get("exit_status") != 0:
                raise ValueError(
                    f"Accepted stage lacks successful native process: {stage}"
                )
            acceptance = json.loads(acceptance_path.read_text())
            checkpoint_path = stage / "restart.ncrst"
            if not checkpoint_path.is_file():
                raise ValueError(f"Accepted stage lacks restart: {stage}")
            if acceptance.get("checkpoint_sha256") != campaign.digest(checkpoint_path):
                raise ValueError(f"Accepted checkpoint digest mismatch: {stage}")
            runtime.stage_acceptance(
                stage,
                spec["kind"],
                cfg["heating_temperature_gate_K"]
                if spec["kind"] == "heat"
                else cfg["temperature_gate_K"],
            )
            _checkpoint(checkpoint_path, baseline, movable, box)
            completed.append(spec["label"])
            restart = checkpoint_path
            production += spec["kind"] == "production"
            records.append({"stage": spec["label"], "state": "accepted_existing"})
            continue
        has_attempt = process_path.exists() or any(
            (stage / name).exists()
            for name in (
                "input_identity.json",
                "restart.ncrst",
                "run.mdout",
                "trajectory.nc",
                "engine.log",
                "input.ncrst",
                "topology.prmtop",
            )
        )
        if has_attempt:
            if spec["kind"] == "production" and production >= production_limit:
                raise ValueError(
                    "An attempted production stage exceeds requested block limit"
                )
            if gap or attempted_gap:
                raise ValueError(
                    "Attempted stage is not the next stage after accepted prefix"
                )
            attempted_gap = True
        gap = True
    if production > production_limit:
        raise ValueError("Existing accepted production blocks exceed requested limit")
    return completed, restart, production, records


def _validate_seed(
    cfg: dict[str, Any], row: dict[str, Any]
) -> tuple[Any, Any, np.ndarray, np.ndarray, np.ndarray]:
    """Apply the original topology, movable-mask, box, coordinate and velocity gates."""
    topology = parmed.load_file(cfg["topology"])
    source = Rst7.open(row["seed"])
    baseline = np.asarray(source.coordinates).reshape(-1, 3).copy()
    if baseline.shape != (len(topology.atoms), 3) or not np.isfinite(baseline).all():
        raise ValueError("Invalid/nonfinite original seed coordinates")
    active = np.array(list(AmberMask(topology, cfg["active_mask"]).Selected()))
    if len(active) != 4165 or len(set(active.tolist())) != 4165:
        raise ValueError("Frozen-boundary active atom count changed from 4165")
    movable = np.zeros(len(baseline), dtype=bool)
    movable[active] = True
    return topology, source, baseline, movable, np.asarray(source.box).copy()


def run_segment(
    prepared: Path,
    label: str,
    production_limit: int,
    adoption_path: Path | None = None,
) -> dict[str, Any]:
    """Run or adopt one window through an absolute production block limit.

    Args:
        prepared: Immutable prepared campaign JSON.
        label: Existing umbrella window label.
        production_limit: Absolute completed production block count, from one to six.
        adoption_path: Optional scheduler receipt for one completed native stage.

    Returns:
        Segment receipt with status, completed stages, production time and hashes.
    """
    if type(production_limit) is not int or not 1 <= production_limit <= 6:
        raise ValueError("production_limit must be an integer from 1 through 6")
    cfg = campaign.frozen(prepared)
    row = next((item for item in cfg["windows"] if item["label"] == label), None)
    if row is None:
        raise ValueError(f"Unknown umbrella window: {label}")
    output = prepared.parent
    folder = output / label
    status_path = folder / "window_status.json"
    had_status_on_entry = status_path.exists()
    backup_path = folder / "coverage_first_previous_window_status.json"
    receipt_path = folder / f"segment_receipt_limit_{production_limit:02d}.json"
    if receipt_path.exists():
        raise FileExistsError(f"Segment receipt already exists: {receipt_path}")
    if status_path.exists():
        old_status = json.loads(status_path.read_text())
        if old_status.get("status") == "failed":
            raise ValueError("Previously failed window status cannot be resumed")
    topology, source, baseline, movable, box = _validate_seed(cfg, row)
    completed, restart, production, records = _accepted_prefix(
        cfg, row, folder, production_limit, baseline, movable, box
    )
    if status_path.exists() and old_status.get("status") == "running":
        running_label = old_status.get("stage")
        running_stage = folder / str(running_label)
        if (
            running_label not in completed
            and not (running_stage / "process.json").is_file()
        ):
            raise FileExistsError(
                "Running status marks a stage attempted without native process evidence"
            )
    if production > production_limit:
        raise ValueError("Existing accepted production blocks exceed requested limit")
    environment = {
        "AMBERHOME": str(Path(cfg["sander"]).parent.parent),
        "LD_LIBRARY_PATH": str(Path(cfg["sander"]).parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    source_hashes = {
        str(path): campaign.digest(path) for path in campaign.source_paths()
    }
    source_hashes[str(Path(__file__).resolve())] = campaign.digest(Path(__file__))
    adoption_used = False
    for spec in row["stages"]:
        if spec["label"] in completed:
            continue
        if spec["kind"] == "production" and production >= production_limit:
            break
        stage = folder / spec["label"]
        stage.mkdir(parents=True, exist_ok=True)
        process_path = stage / "process.json"
        acceptance_path = stage / "acceptance.json"
        if acceptance_path.exists():
            raise ValueError(f"Accepted stage is outside the validated prefix: {stage}")
        if process_path.exists():
            if adoption_path is None or adoption_used:
                raise FileExistsError(
                    f"Attempted native stage cannot be rerun: {stage}"
                )
            _adoption_matches(adoption_path, stage, process_path)
            adoption_used = True
        elif any(
            (stage / name).exists()
            for name in (
                "input_identity.json",
                "run.mdout",
                "restart.ncrst",
                "trajectory.nc",
                "engine.log",
                "input.ncrst",
                "topology.prmtop",
            )
        ):
            raise FileExistsError(
                f"Attempted native outputs lack process evidence: {stage}"
            )
        if had_status_on_entry and status_path.exists() and not backup_path.exists():
            with backup_path.open("xb") as stream:
                stream.write(status_path.read_bytes())
        campaign.write_json(
            status_path,
            {"status": "running", "stage": spec["label"], "completed": completed},
        )
        was_adopted = process_path.exists()
        try:
            if not was_adopted:
                campaign.write_json(
                    stage / "input_identity.json",
                    {
                        "source_sha256": {
                            str(path): campaign.digest(path)
                            for path in [
                                prepared,
                                Path(cfg["topology"]),
                                restart,
                                stage / "run.mdin",
                                stage / "umbrella.RST",
                            ]
                        },
                        "environment": environment,
                        "started_utc": datetime.now(timezone.utc).isoformat(),
                    },
                )
                restart = runtime._run_stage(
                    stage,
                    Path(cfg["sander"]),
                    Path(cfg["topology"]),
                    restart,
                    environment,
                )
            else:
                restart = stage / "restart.ncrst"
            gate = (
                cfg["heating_temperature_gate_K"]
                if spec["kind"] == "heat"
                else cfg["temperature_gate_K"]
            )
            acceptance = runtime.stage_acceptance(stage, spec["kind"], gate)
            checkpoint_record, xyz = _checkpoint(restart, baseline, movable, box)
            acceptance.update(checkpoint_record)
            acceptance["observables"] = runtime.observables(
                xyz, topology, cfg["roles_1based"]
            )
            campaign.write_json(acceptance_path, acceptance)
            if spec["kind"] == "production":
                diagnostics.block_diagnostics(cfg, row, stage, topology)
                diagnostics.export_view(cfg, row, folder, topology)
        except Exception as error:
            campaign.write_json(
                status_path,
                {
                    "status": "failed",
                    "completed": completed,
                    "error": f"{type(error).__name__}: {error}",
                },
            )
            raise
        completed.append(spec["label"])
        records.append(
            {"stage": spec["label"], "state": "adopted" if was_adopted else "completed"}
        )
        if spec["kind"] == "production":
            production += 1
    if adoption_path is not None and not adoption_used:
        raise ValueError(
            "Scheduler adoption did not identify an unaccepted attempted stage"
        )
    status = "completed" if production == 6 else "paused_at_block_boundary"
    production_ps = production * 5
    status_record = {
        "status": status,
        "completed": completed,
        "production_blocks": production,
        "production_ps": production_ps,
        "sampling_convergence": "Not established by completion",
    }
    if had_status_on_entry and status_path.exists() and not backup_path.exists():
        with backup_path.open("xb") as stream:
            stream.write(status_path.read_bytes())
    if (
        status_path.exists()
        and status_path.read_text()
        != json.dumps(status_record, indent=2, allow_nan=False) + "\n"
    ):
        campaign.write_json(status_path, status_record)
    receipt = {
        **status_record,
        "label": label,
        "production_limit": production_limit,
        "completed_stage_records": records,
        "last_checkpoint": str(restart),
        "last_checkpoint_sha256": campaign.digest(restart),
        "prepared_sha256": campaign.digest(prepared),
        "source_sha256": source_hashes,
        "adoption_receipt_sha256": (
            campaign.digest(adoption_path) if adoption_path is not None else None
        ),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    campaign.write_json(receipt_path, receipt)
    LOGGER.info(
        "Umbrella segment %s stopped at %d production blocks", label, production
    )
    return receipt
