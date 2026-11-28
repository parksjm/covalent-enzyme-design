"""Bounded final100 candidate watcher; never changes native sampling or primary aliases."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from nyl12_covalent.relay_final_analysis import analyze
from nyl12_covalent.relay_final_sources import (
    BASE,
    COHORTS,
    MIXING,
    final_data,
    readiness,
)
from nyl12_covalent.relay_pooled_analysis import write
from nyl12_covalent.relay_umbrella_analysis_io import digest


def journal(phase: str, action: str) -> None:
    """Append autonomous derived-monitor intent/outcome before dependent work."""
    subprocess.run(
        [
            sys.executable,
            "/data/jpc/covalent-enzyme-design/scripts/record_agent_action.py",
            "--log",
            "/data/jpc/covalent-enzyme-design/output_agent_actions/20261021.jsonl",
            "--agent",
            "relay_final100_watcher",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Authorized isolated final100 candidate; no native sampling changes",
            "--rollback",
            "Retain frozen source and derived evidence; stop only this derived monitor",
        ],
        check=True,
    )


def sync_commands(binding: dict) -> list[list[str]]:
    """Return six derived-only pulls for an explicitly reviewed host/prepared binding."""
    if (
        binding["host"] not in {"azalea", "geranium", "verbena"}
        or binding["reviewed"] is not True
    ):
        raise ValueError("Explicit reviewed science-host binding required")
    if digest(MIXING) != binding["prepared_sha256"]:
        raise ValueError("Reviewed host prepared identity mismatch")
    ssh = "ssh -i " + str(Path(binding["ssh_key"])) + " -o BatchMode=yes"
    commands = []
    for window in json.loads(MIXING.read_text())["windows"]:
        folder = MIXING.parent / window["label"]
        commands.append(
            [
                "rsync",
                "-a",
                "--prune-empty-dirs",
                "-e",
                ssh,
                "--include=*/",
                "--include=*.json",
                "--include=*.csv",
                "--include=umbrella_coordinate.dat",
                "--exclude=*",
                binding["host"]
                + ":"
                + str(Path(binding["remote_root"]) / window["label"])
                + "/",
                str(folder) + "/",
            ]
        )
    return commands


def watch(
    output: Path,
    source_manifest: Path,
    host_binding: Path,
    hours: float = 72,
    poll: float = 120,
) -> None:
    """Wait for actual60+30 accepted cohorts and run the candidate once."""
    if not 0 < hours <= 72 or poll < 60:
        raise ValueError("Require bounded72hours and >=60second polling")
    pinned = json.loads(source_manifest.read_text())
    for filename, sha in pinned["files"].items():
        if digest(source_manifest.parent / filename) != sha:
            raise ValueError("Frozen source changed:" + filename)
    for filename, sha in pinned["static_inputs"].items():
        if digest(Path(filename)) != sha:
            raise ValueError("Pinned projection/protocol changed:" + filename)
    identity = {str(path): digest(path) for path, *_ in COHORTS}
    if identity != pinned["external_inputs"]:
        raise ValueError("Prepared manifests do not match frozen final inventory")
    binding = json.loads(host_binding.read_text())
    binding_sha = digest(host_binding)
    commands = sync_commands(binding)
    write(
        output / "derived_sync_manifest.json",
        dict(
            commands=commands,
            native_writes=False,
            excluded="All trajectories/restarts/topologies/native MD outputs",
        ),
    )
    accepted = output / "main60_supplement43x30/acceptance.json"
    if accepted.exists():
        proof = json.loads(accepted.read_text())
        point = accepted.parent / "point.json"
        if (
            proof["states"] != 100
            or proof["actual_main_ps"] != 60
            or proof["actual_supplement_ps"] != 30
            or digest(point) != proof["point_sha256"]
        ):
            raise ValueError("Existing final candidate identity mismatch")
        return
    deadline = time.monotonic() + hours * 3600
    while time.monotonic() < deadline:
        if digest(host_binding) != binding_sha:
            raise ValueError("Reviewed host binding changed during running attempt")
        for filename, sha in identity.items():
            if digest(Path(filename)) != sha:
                raise ValueError("Pinned prepared input changed:" + filename)
        journal(
            "intent",
            "Read-only reviewed-host six-lineage derived sync and readiness check",
        )
        for command in commands:
            subprocess.run(command, check=True, timeout=120)
        status = readiness()
        journal(
            "outcome",
            "Derived sync complete; actual accepted cohort readiness recorded",
        )
        write(output / "monitor_state.json", dict(status="waiting", **status))
        if status["ready"] and not (BASE / ".compact_sync_in_progress").exists():
            for filename, sha in pinned["static_inputs"].items():
                if digest(Path(filename)) != sha:
                    raise ValueError(
                        "Pinned projection/protocol changed before final fitting:"
                        + filename
                    )
            journal(
                "intent",
                "Validate all100 actual lineages; fit final candidate and paired conditional uncertainty at max4CPU",
            )
            data = final_data()
            if data is None:
                raise ValueError("Readiness changed during final source validation")
            analyze(data, output / "main60_supplement43x30", workers=4)
            journal(
                "outcome",
                "Final100 candidate complete; stationarity unresolved and no primary promotion",
            )
            write(
                output / "monitor_state.json",
                dict(
                    status="candidate_complete",
                    **status,
                    primary_replacement=False,
                    sampling_converged=False,
                ),
            )
            return
        time.sleep(poll)
    raise TimeoutError("Final candidate monitor reached explicit72h bound")


def main() -> None:
    """Run the isolated final candidate monitor from frozen sources."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--host-binding", type=Path, required=True)
    parser.add_argument("--hours", type=float, default=72)
    parser.add_argument("--poll-seconds", type=float, default=120)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    watch(
        args.output,
        args.source_manifest,
        args.host_binding,
        args.hours,
        args.poll_seconds,
    )


if __name__ == "__main__":
    main()
