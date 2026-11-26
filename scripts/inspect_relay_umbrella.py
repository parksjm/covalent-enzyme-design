#!/usr/bin/env python3
"""Write a compact process/progress observation without touching native inputs."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re


def inspect(campaign: Path, execution: Path) -> dict:
    """Observe exact campaign Sander processes and interval progress rates."""
    now = datetime.now(timezone.utc)
    target = execution / "runtime_summary.json"
    previous = json.loads(target.read_text()) if target.exists() else {}
    old = {row["pid"]: row for row in previous.get("sander_processes", [])}
    processes = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            argv = (entry / "cmdline").read_bytes().split(b"\0")
            cwd = (entry / "cwd").resolve()
            if (
                not argv
                or argv[0].decode(errors="replace")
                != "/data/jpc/code/ambertools26/venv/bin/sander"
                or not cwd.is_relative_to(campaign)
            ):
                continue
            record = dict(
                pid=int(entry.name), cwd=str(cwd), argv=[v.decode() for v in argv if v]
            )
            mdinfo = cwd / "mdinfo"
            rows = (
                re.findall(
                    r"NSTEP\s*=\s*(\d+).*?TEMP\(K\)\s*=\s*([-+\d.EeDd]+)",
                    mdinfo.read_text(),
                )
                if mdinfo.exists()
                else []
            )
            native_output = cwd / "run.mdout"
            actual_steps = []
            if native_output.exists():
                with native_output.open("rb") as stream:
                    stream.seek(max(0, native_output.stat().st_size - 65536))
                    tail = stream.read().decode(errors="replace")
                actual_steps = re.findall(r"SCC-DFTB for step\s+(\d+) converged", tail)
            if actual_steps:
                cycles = re.findall(r"converged in\s+(\d+) cycles", tail)
                record["last_SCC_cycles"] = int(cycles[-1]) if cycles else None
                record["recent_unresolved_warning_markers"] = len(
                    re.findall(r"DID NOT CONVERGE|FATAL|ERROR:", tail, re.I)
                )
            statm = (entry / "statm").read_text().split()
            record["rss_bytes"] = int(statm[1]) * 4096
            if rows or actual_steps:
                record.update(
                    step=int(actual_steps[-1]) if actual_steps else int(rows[-1][0]),
                    progress_source="native_SCC_completed_step"
                    if actual_steps
                    else "native_mdinfo",
                    temperature_K=float(rows[-1][1]) if rows else None,
                    observed_utc=now.isoformat(),
                )
                prior = old.get(int(entry.name))
                if prior and "step" in prior and record["step"] > prior["step"]:
                    elapsed = (
                        now - datetime.fromisoformat(prior["observed_utc"])
                    ).total_seconds()
                    record.update(
                        interval_steps=record["step"] - prior["step"],
                        interval_wall_s=elapsed,
                        interval_s_per_step=elapsed / (record["step"] - prior["step"]),
                    )
            processes.append(record)
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
    cfg = json.loads((campaign / "prepared.json").read_text())
    workers = [
        json.loads(p.read_text())
        for p in (execution / "workers").glob("*/job/status.json")
    ]
    windows = {
        row["label"]: json.loads(p.read_text())
        for row in cfg["windows"]
        if (p := campaign / row["label"] / "window_status.json").exists()
    }
    data = dict(
        observed_utc=now.isoformat(),
        window_count=len(cfg["windows"]),
        worker_count=len(workers),
        active_sander_count=len(processes),
        sander_processes=processes,
        window_status=windows,
        queue_claims=[
            json.loads(p.read_text()) for p in (execution / "claims").glob("*.json")
        ],
        worker_failures=[
            v for v in workers if v["ended_utc"] and v["exit_status"] != 0
        ],
        sampling_converged=False,
    )
    pending = target.with_suffix(".pending." + str(__import__("os").getpid()))
    pending.write_text(json.dumps(data, indent=2) + "\n")
    pending.replace(target)
    stable = campaign / "runtime_summary.json"
    temporary = stable.with_suffix(".pending." + str(__import__("os").getpid()))
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(stable)
    plan = execution / "plan.json"
    pointer = dict(
        observed_utc=now.isoformat(),
        execution_directory=str(execution),
        plan=str(plan),
        runtime_summary=str(target),
        active_sander_count=len(processes),
        role="Currentnativeumbrellaexecution; heatingacceptanceisramp-only, equilibriumrequires300Kequil/productiongates",
        sampling_converged=False,
    )
    pointer_path = campaign / "current_execution.json"
    temporary = pointer_path.with_suffix(".pending." + str(__import__("os").getpid()))
    temporary.write_text(json.dumps(pointer, indent=2) + "\n")
    temporary.replace(pointer_path)
    return data


def main() -> None:
    """Inspect one explicit campaign and controller path."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--execution", type=Path, required=True)
    args = parser.parse_args()
    data = inspect(
        args.campaign.resolve(strict=True), args.execution.resolve(strict=True)
    )
    print(
        json.dumps(
            {
                k: data[k]
                for k in (
                    "observed_utc",
                    "active_sander_count",
                    "worker_count",
                    "window_count",
                    "worker_failures",
                )
            }
        )
    )


if __name__ == "__main__":
    main()
