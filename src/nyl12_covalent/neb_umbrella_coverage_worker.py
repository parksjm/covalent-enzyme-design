"""Adopt a finishing native block, fill coverage, then resume balanced sampling."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import time

from .neb_umbrella_campaign import digest, write_json
from .neb_umbrella_coverage import context, event
from .neb_umbrella_handoff import native_completion, release_driver
from .neb_umbrella_segments import run_segment


def completed_blocks(root: Path, label: str) -> int:
    """Count only accepted production blocks in a window's original directory."""
    return len(list((root / label).glob("03_production_block_*_5ps/acceptance.json")))


def worker(schedule: Path, host: str, index: int, interval: int = 30) -> None:
    """Finish the active block, run the queued first block, and obey a global gate."""
    _, _, prepared = context(schedule, host)
    root = prepared.parent
    inventory = json.loads((schedule.parent / host / "inventory.json").read_text())
    row = inventory["workers"][index]
    folder = schedule.parent / host / f"worker{index + 1:02d}"
    status = folder / "status.json"
    try:
        write_json(
            status,
            dict(
                status="waiting_for_current_native_block",
                window=row["current"],
                native_pid=row["native"]["pid"],
            ),
        )
        while (code := native_completion(row["native"])) is None:
            time.sleep(interval)
        if code != 0:
            raise RuntimeError(f"Original native block failed: exit {code}")
        stage = Path(row["native"]["cwd"])
        adoption = folder / "native_block_adoption.json"
        write_json(
            adoption,
            dict(
                stage=str(stage),
                native_exit_status=code,
                raw_native_wait_status=0,
                native_process_identity=row["native"],
                process_json_sha256=digest(stage / "process.json"),
                native_outputs_complete=True,
                observed_utc=datetime.now(timezone.utc).isoformat(),
            ),
        )
        event(
            folder.parent,
            "native_block_completed",
            window=row["current"],
            adoption=str(adoption),
        )
        release_driver(row["driver"])
        current_limit = int(stage.name.split("_")[3])
        run_segment(prepared, row["current"], current_limit, adoption_path=adoption)
        for label in row["queued"]:
            write_json(status, dict(status="running_first_block", window=label))
            run_segment(prepared, label, 1)
        write_json(
            status,
            dict(
                status="waiting_for_all_63_first_blocks",
                windows=[row["current"], *row["queued"]],
            ),
        )
        gate = schedule.parent / host / "coverage_gate.json"
        while not gate.exists():
            time.sleep(interval)
        release = json.loads(gate.read_text())
        if (
            release["status"] != "all_windows_have_first_block"
            or release["windows"] != 63
        ):
            raise ValueError("Invalid coverage release")
        labels = [*row["queued"], row["current"]]
        while min(completed_blocks(root, label) for label in labels) < 6:
            label = min(labels, key=lambda label: completed_blocks(root, label))
            limit = completed_blocks(root, label) + 1
            write_json(
                status,
                dict(
                    status="extending_after_global_coverage",
                    window=label,
                    target_blocks=limit,
                ),
            )
            run_segment(prepared, label, limit)
        write_json(
            status,
            dict(
                status="completed",
                windows=labels,
                production_ps_per_window=30,
                sampling_convergence="Not implied by duration",
            ),
        )
    except BaseException as error:
        write_json(
            status, dict(status="failed", error=f"{type(error).__name__}: {error}")
        )
        raise
