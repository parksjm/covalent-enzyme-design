"""Wait for a completed frequency controller before conditional-curvature analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import subprocess
import sys
import time
from pathlib import Path

LOGGER = logging.getLogger(__name__)


def verify_sources(snapshot: Path) -> None:
    """Require every pinned analysis source to match its recorded digest."""
    identities = json.loads((snapshot / "source_identity.json").read_text())
    if not isinstance(identities, dict) or not identities:
        raise ValueError("Empty source identity")
    for name, expected in identities.items():
        path = snapshot / name
        if not path.is_relative_to(snapshot) or ".." in Path(name).parts:
            raise ValueError("Source identity escapes snapshot")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Pinned source changed: {name}")


def wait_for_completion(
    status_path: Path, wait_hours: float, poll_seconds: int
) -> None:
    """Wait for an explicit successful end record, or fail without analysis."""
    if wait_hours <= 0 or poll_seconds <= 0:
        raise ValueError("Positive wait duration and polling interval required")
    deadline = time.monotonic() + wait_hours * 3600
    while True:
        if status_path.is_file():
            status = json.loads(status_path.read_text())
            if status.get("ended_utc") is not None:
                if status.get("exit_status") != 0:
                    raise RuntimeError(
                        "Frequency controller did not finish successfully"
                    )
                return
        if time.monotonic() >= deadline:
            raise TimeoutError("Frequency controller wait expired")
        time.sleep(poll_seconds)


def run_analysis(
    case: Path,
    snapshot: Path,
    journal: Path,
    *,
    wait_hours: float,
    poll_seconds: int,
) -> int:
    """Run the pinned analyzer only after source and process prerequisites pass."""
    python = sys.executable
    recorder = snapshot / "scripts/record_agent_action.py"

    def record(phase: str, detail: str) -> None:
        subprocess.run(
            [
                python,
                str(recorder),
                "--log",
                str(journal),
                "--agent",
                "root",
                "--phase",
                phase,
                "--action",
                "conditional_curvature_analysis",
                "--permission-mode",
                "escalated",
                "--reason",
                detail,
                "--rollback",
                "Preserve analysis output; use a fresh directory for retries",
            ],
            check=True,
        )

    verify_sources(snapshot)
    record("intent", "Authorized post-frequency check; no fitting or engine launch")
    try:
        wait_for_completion(
            case / "frequency_controller_01/job/status.json", wait_hours, poll_seconds
        )
        verify_sources(snapshot)
        output = case / "conditional_curvature_01"
        if output.exists():
            raise FileExistsError(output)
        command = [
            python,
            str(snapshot / "scripts/analyze_nyl12_ti1_curvature.py"),
            "--frequency-dir",
            str(case / "frequency_01"),
            "--output",
            str(output),
            "--step-bohr",
            "0.0002",
            "--stationarity-tolerance",
            "0.00003",
            "--curvature-tolerance",
            "0.000001",
            "--derivative-relative-tolerance",
            "0.001",
        ]
        result = subprocess.run(command, check=False)
    except Exception as error:
        record("outcome", f"Analysis withheld/failed: {type(error).__name__}: {error}")
        raise
    record("outcome", f"Conditional-curvature analysis exit status {result.returncode}")
    LOGGER.info("Analysis exit status %d; output %s", result.returncode, output)
    return result.returncode


def main(argv: list[str] | None = None) -> int:
    """Parse explicit paths and wait limits, then run the gated analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--journal", type=Path, required=True)
    parser.add_argument("--wait-hours", type=float, required=True)
    parser.add_argument("--poll-seconds", type=int, required=True)
    args = parser.parse_args(argv)
    return run_analysis(
        args.case,
        args.snapshot,
        args.journal,
        wait_hours=args.wait_hours,
        poll_seconds=args.poll_seconds,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
