#!/usr/bin/env python3
"""Choose a material measured affinity improvement while retaining57MDstreams."""

from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace


def choose_affinity(before: float, physical: float, after: float) -> dict:
    """Require gain above3percent and twiceA-B-A baseline drift."""
    if min(before, physical, after) <= 0:
        raise ValueError("MeasuredaggregateMDthroughput mustbepositive")
    drift = abs(after - before) / ((before + after) / 2)
    threshold = max(0.03, 2 * drift)
    gain = physical / max(before, after) - 1
    return dict(
        winner="physical_28" if gain > threshold else "all_56_logical",
        physical_gain_vs_faster_default=gain,
        baseline_relative_drift=drift,
        material_gain_threshold=threshold,
        inconclusive=gain <= threshold,
        comparison="Boundedmixed-stageA-B-Aresourcebenchmark; notscientificsamplingconvergence",
    )


def main() -> None:
    """Wait forsafe restoredcomparison, repeatbaseline, then applyreviewedpolicy."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("operations", "plan", "campaign", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--duration-minutes", type=float, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    deadline = time.monotonic() + args.duration_minutes * 60
    while True:
        status = json.loads((args.operations / "benchmark_job/status.json").read_text())
        if status["ended_utc"] is not None:
            if status["exit_status"] != 0:
                raise RuntimeError(
                    "Initialresourcebenchmarkfailed; leaveMDresourcepolicyforreview"
                )
            break
        if time.monotonic() > deadline:
            raise TimeoutError("Initialbenchmarkcompletiondeadline")
        time.sleep(20)
    source = args.operations / "scripts/benchmark_relay_concurrency.py"
    __import__("sys").path.insert(0, str(source.parent))
    spec = importlib.util.spec_from_file_location("pinned_benchmark", source)
    if spec is None or spec.loader is None:
        raise ValueError("Pinnedbenchmarkmoduleunavailable")
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    results = args.operations / "benchmark_results"
    completion = json.loads((results / "completed.json").read_text())
    if completion["affinities_restored"] is not True:
        raise ValueError("Originalaffinitiesmustberestoredfirst")
    repeated = benchmark.sample(
        SimpleNamespace(campaign=args.campaign, plan=args.plan, output=results),
        57,
        "workers_57_all_logical_repeat",
        deadline,
    )
    first = json.loads((results / "workers_57.json").read_text())
    physical = json.loads((results / "workers_57_physical_only.json").read_text())
    decision = choose_affinity(
        first["aggregate_steps_per_s"],
        physical["aggregate_steps_per_s"],
        repeated["aggregate_steps_per_s"],
    )
    topology = json.loads((results / "affinity_before.json").read_text())
    chosen = (
        topology["physical_mask"]
        if decision["winner"] == "physical_28"
        else sorted(os.sched_getaffinity(0))
    )
    cfg = json.loads(args.plan.read_text())
    controller = str(args.plan.parent / "source_snapshot/scripts/run_relay_umbrella.py")
    originals = {}
    identities = {}
    for row in repeated["native_processes"]:
        pid = row["pid"]
        path = Path("/proc") / str(pid)
        argv = [p.decode() for p in (path / "cmdline").read_bytes().split(b"\0") if p]
        if argv[0] != cfg["sander"] or not (path / "cwd").resolve().is_relative_to(
            args.campaign
        ):
            raise ValueError("NativePIDidentitychangedbeforepolicyselection")
        parent = int((path / "stat").read_text().split(") ", 1)[1].split()[1])
        parent_argv = [
            p.decode()
            for p in (Path("/proc") / str(parent) / "cmdline").read_bytes().split(b"\0")
            if p
        ]
        if (
            len(parent_argv) < 2
            or parent_argv[1] != controller
            or str(args.plan) not in parent_argv
        ):
            raise ValueError("Nativeworkerparentidentitymismatch")
        for target, command in ((pid, argv), (parent, parent_argv)):
            originals[target] = os.sched_getaffinity(target)
            identities[target] = command
    decision.update(
        chosen_cpu_ids=chosen,
        native_streams=57,
        targets=identities,
        previous_affinities={p: sorted(mask) for p, mask in originals.items()},
    )
    (args.output / "decision_before_apply.json").write_text(
        json.dumps(decision, indent=2) + "\n"
    )
    benchmark.journal(
        "intent",
        "Applymeasuredresourcewinner "
        + decision["winner"]
        + " to57verifiednativeprocessesandworkerparents "
        + str(chosen),
    )
    try:
        for pid, command in identities.items():
            actual = [
                p.decode()
                for p in (Path("/proc") / str(pid) / "cmdline")
                .read_bytes()
                .split(b"\0")
                if p
            ]
            if actual != command:
                raise ValueError("Verifiedprocessidentitychangedbeforeaffinityapply")
            os.sched_setaffinity(pid, chosen)
    except Exception:
        for pid, mask in originals.items():
            actual = [
                p.decode()
                for p in (Path("/proc") / str(pid) / "cmdline")
                .read_bytes()
                .split(b"\0")
                if p
            ]
            if actual != identities[pid]:
                raise ValueError("Verifiedprocessidentitychangedbeforeresourcerollback")
            os.sched_setaffinity(pid, mask)
        benchmark.journal(
            "outcome",
            "Resourcepolicyapplyfailed; allverifiedoriginalaffinitiesrestored",
        )
        raise
    decision["applied"] = True
    (args.output / "decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    benchmark.journal(
        "outcome",
        "Measuredresourcewinnerapplied;57nativeMDstreamscontinue "
        + json.dumps(
            {
                k: v
                for k, v in decision.items()
                if k not in ("targets", "previous_affinities")
            }
        ),
    )


if __name__ == "__main__":
    main()
