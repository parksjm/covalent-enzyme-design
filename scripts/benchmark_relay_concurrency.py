#!/usr/bin/env python3
"""Benchmark reviewed worker expansion without restarting any native MD."""

from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import time
from inspect_relay_umbrella import inspect

ROOT = Path("/data/jpc/covalent-enzyme-design")


def journal(phase: str, action: str) -> None:
    """Record benchmark and affinity operations before/after execution."""
    subprocess.run(
        [
            "python3",
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_resource_benchmark",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Rootauthorised28→40→57streams andphysical-vs-SMTcomparison",
            "--rollback",
            "NoMDrestart; originalnativeCPUaffinity restored aftercomparison",
        ],
        check=True,
    )


def physical_cpu_ids(text: str) -> list[int]:
    """Choose one logical CPU per physical socket/core from lscpu CSV."""
    first = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        cpu, core, socket = (int(part) for part in line.split(","))
        first.setdefault((socket, core), cpu)
    if not first:
        raise ValueError("No physicalCPUtopology records")
    return sorted(first.values())


def record(folder: Path, name: str, data: dict) -> None:
    """Retain one immutable benchmark observation."""
    with (folder / (name + ".json")).open("x") as stream:
        json.dump(data, stream, indent=2)
        stream.write("\n")


def sample(args: argparse.Namespace, count: int, name: str, deadline: float) -> dict:
    """Measure at least100nativeMDsteps per stableprocess, excludingstartup."""
    while True:
        start = inspect(args.campaign, args.plan.parent)
        rows = start["sander_processes"]
        if start["worker_failures"]:
            raise RuntimeError("Nativeworkerfailure duringbenchmark")
        if len(rows) == count and all("step" in r for r in rows):
            break
        if time.monotonic() > deadline:
            raise TimeoutError("Nativeprocessstartupdeadline")
        time.sleep(20)
    before = {r["pid"]: r for r in rows}
    begun = time.monotonic()
    while True:
        time.sleep(20)
        after = inspect(args.campaign, args.plan.parent)
        current = {r["pid"]: r for r in after["sander_processes"]}
        if set(current) != set(before):
            raise RuntimeError(
                "NativePIDcohortchanged; retainpartialbenchmark forreview"
            )
        deltas = {pid: current[pid]["step"] - before[pid]["step"] for pid in before}
        if min(deltas.values()) >= 100:
            break
        if time.monotonic() > deadline:
            raise TimeoutError("100stepbenchmarkdeadline")
    elapsed = time.monotonic() - begun
    data = dict(
        observed_utc=datetime.now(timezone.utc).isoformat(),
        workers=count,
        elapsed_s=elapsed,
        delta_steps=deltas,
        total_delta_steps=sum(deltas.values()),
        aggregate_steps_per_s=sum(deltas.values()) / elapsed,
        seconds_per_step_per_worker=elapsed * count / sum(deltas.values()),
        rss_GiB=sum(r["rss_bytes"] for r in current.values()) / 1024**3,
        load_average=os.getloadavg(),
        native_processes=list(current.values()),
        affinity={pid: sorted(os.sched_getaffinity(pid)) for pid in current},
        sampling_converged=False,
    )
    record(args.output, name, data)
    journal(
        "outcome",
        f"ActualMDresourcebenchmark{name}: "
        + json.dumps(
            {
                k: v
                for k, v in data.items()
                if k not in ("delta_steps", "native_processes", "affinity")
            }
        ),
    )
    return data


def main() -> None:
    """Execute explicit57stream expansion and restore resource-affinity changes."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "campaign", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--duration-minutes", type=float, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    plan = json.loads(args.plan.read_text())
    if plan["workers"] != 57 or plan["initial_workers"] != 28:
        raise ValueError("Rootreviewed57budget/28initialrequired")
    deadline = time.monotonic() + args.duration_minutes * 60
    controller = args.plan.parent / "source_snapshot/scripts/run_relay_umbrella.py"
    journal("intent", "StartboundedactualMDresourcebenchmark28→40→57")
    sample(args, 28, "workers_28", deadline)
    for first, count, total in ((28, 12, 40), (40, 17, 57)):
        journal(
            "intent",
            f"Expandnativeworkerbudget{first}→{total}; unchangedscience/RNG/checkpoints",
        )
        subprocess.run(
            [
                plan["python"],
                str(controller),
                "expand",
                "--plan",
                str(args.plan),
                "--first",
                str(first),
                "--count",
                str(count),
            ],
            check=True,
        )
        journal("outcome", f"Launchednativeworkers{first}..{total - 1}")
        sample(args, total, f"workers_{total}", deadline)
    current = inspect(args.campaign, args.plan.parent)
    original = {
        r["pid"]: os.sched_getaffinity(r["pid"]) for r in current["sander_processes"]
    }
    topology = subprocess.run(
        ["lscpu", "-p=CPU,CORE,SOCKET"], capture_output=True, text=True, check=True
    ).stdout
    physical = physical_cpu_ids(topology)
    record(
        args.output,
        "affinity_before",
        dict(
            original={pid: sorted(mask) for pid, mask in original.items()},
            physical_mask=physical,
        ),
    )
    journal(
        "intent",
        "Temporarilyrestrict57verifiednativeSanderprocesses tooneCPUperphysicalcore "
        + str(physical),
    )
    try:
        for pid in original:
            argv = (Path("/proc") / str(pid) / "cmdline").read_bytes().split(b"\0")
            if argv[0].decode() != "/data/jpc/code/ambertools26/venv/bin/sander":
                raise ValueError("NativePIDidentitychanged")
            os.sched_setaffinity(pid, physical)
        journal("outcome", "57nativeSanderCPUaffinitiesrestricted; noMDrestart")
        sample(args, 57, "workers_57_physical_only", deadline)
    finally:
        journal("intent", "RestoreallverifiednativeCPUaffinities frombeforecomparison")
        for pid, mask in original.items():
            argv = (Path("/proc") / str(pid) / "cmdline").read_bytes().split(b"\0")
            if argv[0].decode() != "/data/jpc/code/ambertools26/venv/bin/sander":
                raise ValueError("NativePIDidentitychangedbeforerestore")
            os.sched_setaffinity(pid, mask)
        journal("outcome", "OriginalnativeCPUaffinities restored; all57streamscontinue")
    record(
        args.output,
        "completed",
        dict(
            completed_utc=datetime.now(timezone.utc).isoformat(),
            native_streams=57,
            affinities_restored=True,
        ),
    )


if __name__ == "__main__":
    main()
