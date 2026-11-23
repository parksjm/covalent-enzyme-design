"""Wait for force-accepted w033, then verify full-shell local curvature once."""

from __future__ import annotations

import argparse
import logging
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from nyl12_covalent.relay_minimum import (
    freeze,
    geometry_screen,
    job_records,
    journal,
    prerequisites,
    proton_inventory,
)
from nyl12_covalent.relay_release import geometry
from nyl12_covalent.stationary_native import NativeState, write_json

LOGGER = logging.getLogger(__name__)


def resource_gate(output: Path, workers: int) -> dict[str, object]:
    """Require bounded CPU, memory and disk headroom before derivatives."""
    cpus = os.cpu_count()
    if cpus is None:
        raise RuntimeError("CPU count unavailable")
    load = os.getloadavg()[0]
    memory = {
        line.split(":")[0]: int(line.split()[1]) * 1024
        for line in Path("/proc/meminfo").read_text().splitlines()
    }
    disk = shutil.disk_usage(output).free
    return {
        "passed": cpus - load >= workers
        and memory["MemAvailable"] >= 16 * 2**30
        and disk >= 20 * 2**30,
        "logical_cpus": cpus,
        "load_1min": load,
        "workers": workers,
        "available_memory_bytes": memory["MemAvailable"],
        "free_disk_bytes": disk,
        "requirements": "load headroom >= workers; memory>=16GiB; disk>=20GiB",
    }


def run_owned(command: list[str], output: Path, seconds: float) -> int:
    """Run only the owned verifier process group within the remaining budget."""
    with (output / "verification.log").open("x") as log:
        process = subprocess.Popen(
            command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
        )
        write_json(
            output / "verification_process.json",
            {"pid": process.pid, "command": command, "budget_seconds": seconds},
        )
        try:
            return process.wait(timeout=seconds)
        except subprocess.TimeoutExpired:
            journal(
                "intent",
                f"Controller deadline: stop owned verifier group {process.pid}",
            )
            if process.poll() is None and os.getpgid(process.pid) == process.pid:
                os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                if os.getpgid(process.pid) != process.pid:
                    raise RuntimeError(
                        "Owned verifier group identity changed"
                    ) from None
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            journal(
                "outcome",
                "Owned verifier budget exhausted; original optimizers untouched",
            )
            return 124


def watch(
    campaign: Path,
    reference: Path,
    output: Path,
    interval: float,
    hours: float,
    workers: int,
) -> None:
    """Conditionally verify w033 without interrupting any original optimizer.

    Args:
        campaign: Existing three-release root.
        reference: Validated attempt05 source of original precision library.
        output: Single existing minimum_verification root with launcher subdirectory.
        interval: Poll period in seconds, normally 120.
        hours: Explicit total controller and owned verification budget, normally 168.
        workers: Maximum eight nice single-CPU gradient workers.
    """
    if (
        not all(math.isfinite(x) and x > 0 for x in (interval, hours))
        or not 1 <= workers <= 8
    ):
        raise ValueError("Require positive timing and 1..8 workers")
    deadline = time.monotonic() + hours * 3600
    journal(
        "intent",
        "Check conditional minimum prerequisites; preserve all original job sessions",
    )
    hashes = prerequisites(campaign, reference)
    write_json(output / "prerequisites.json", hashes)
    journal(
        "outcome",
        "All16 SKFs, original masks/inputs and validated precision library matched",
        (output / "prerequisites.json",),
    )
    while True:
        journal("intent", "Read three original job records and w033 conditional gate")
        records = job_records(campaign)
        record = records["ti1_window033"]
        write_json(output / "original_job_records.json", records)
        journal(
            "outcome",
            f"Original statuses: {[(k, v['process']['exit_status']) for k, v in records.items()]}",
        )
        if record["process"]["exit_status"] is not None:
            scientific = record["scientific"]
            if scientific is None:
                raise RuntimeError(
                    "Terminal w033 missing scientific acceptance; inspect retained logs"
                )
            if (
                record["process"]["exit_status"] != 0
                or not scientific.get("accepted")
                or not scientific.get("force_converged")
            ):
                write_json(
                    output / "scientific_acceptance.json",
                    {
                        "numerically_verified": False,
                        "status": "requires_review_or_continuation",
                        "reason": "Original native run not force-accepted; no automatic restart",
                    },
                )
                return
            seed = campaign / "ti1_window033"
            state = NativeState(seed, seed / "optimized.ncrst")
            screen = geometry_screen(geometry(state.xyz))
            screen["qm_proton_inventory"] = proton_inventory(state)
            write_json(output / "geometry_screen.json", screen)
            if screen["decision"] != "retained":
                write_json(
                    output / "scientific_acceptance.json",
                    {
                        "numerically_verified": False,
                        "status": screen["decision"],
                        "reason": "Declared geometry screen; inspect angles and proton inventory manually",
                    },
                )
                return
            resources = resource_gate(output, workers)
            write_json(output / "resource_preflight.json", resources)
            if resources["passed"]:
                break
        if time.monotonic() >= deadline:
            raise TimeoutError(
                "168h-style controller budget exhausted; original optimizers untouched"
            )
        time.sleep(min(interval, max(0, deadline - time.monotonic())))
    journal(
        "intent",
        "Freeze accepted w033 restart and exact original method/input/library hashes",
    )
    if prerequisites(campaign, reference) != hashes:
        raise ValueError("Prerequisites changed while controller waited")
    frozen = freeze(seed, reference, output)
    journal(
        "outcome",
        "Frozen accepted restart and full-shell original inputs",
        (frozen / "manifest.json",),
    )
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--verify",
        "--output",
        str(output),
        "--workers",
        str(workers),
    ]
    journal(
        "intent",
        "Launch owned full-shell verifier after accepted force and geometry gates",
    )
    code = run_owned(command, output, max(0.001, deadline - time.monotonic()))
    write_json(output / "process_acceptance.json", {"exit_status": code})
    journal(
        "outcome", f"Owned minimum verifier exit_status={code}; original jobs unchanged"
    )
    if code:
        raise RuntimeError(f"Minimum verifier failed or budget expired: {code}")


def main() -> None:
    """Run the persistent controller or its isolated derivative subprocess."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval", type=float, default=120)
    parser.add_argument("--hours", type=float, default=168)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    try:
        if args.verify:
            from nyl12_covalent.relay_minimum_run import verify_minimum

            verify_minimum(args.output / "01_frozen", args.output, args.workers)
        else:
            if args.campaign is None or args.reference is None:
                raise ValueError("Controller requires exact campaign and reference")
            watch(
                args.campaign,
                args.reference,
                args.output,
                args.interval,
                args.hours,
                args.workers,
            )
    except Exception as error:
        journal(
            "outcome", f"Conditional minimum failed: {type(error).__name__}: {error}"
        )
        write_json(
            args.output
            / (
                "verification_failure.json"
                if args.verify
                else "controller_failure.json"
            ),
            {
                "numerically_verified": False,
                "error": str(error),
                "reason": "Retain evidence; no automatic retries, threshold changes or new path work",
            },
        )
        raise


if __name__ == "__main__":
    main()
