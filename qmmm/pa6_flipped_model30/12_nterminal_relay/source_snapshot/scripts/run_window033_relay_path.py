#!/usr/bin/env python3
"""Run one relay branch or orchestrate the two independently."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.stationary_native import digest, write_json  # noqa: E402
from nyl12_covalent.nterminal_relay_run import run_branch  # noqa: E402
from nyl12_covalent.relay_path_queue import verify_minimum  # noqa: E402


def require_dispatch(campaign: Path) -> None:
    """Reject direct invocation before verified inertia or after source changes."""
    config = json.loads((campaign / "config.json").read_text())
    if verify_minimum(Path(config["verification_root"])) is None:
        raise RuntimeError("Window033 is not yet minimum-verified; no native dispatch")
    for source, expected in config["source_sha256"].items():
        if digest(Path(source)) != expected:
            raise ValueError(f"Frozen relay source/input changed: {source}")


def main() -> int:
    """Run requested branch, or wait for both independent child processes."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--branch", choices=("reverse_activation", "forward_acylation"))
    args = parser.parse_args()
    campaign = args.campaign.resolve(strict=True)
    require_dispatch(campaign)
    if args.branch:
        return 0 if run_branch(campaign, args.branch)["status"] == "completed" else 1
    branches = ("reverse_activation", "forward_acylation")
    processes = {
        branch: subprocess.Popen(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--campaign",
                str(campaign),
                "--branch",
                branch,
            ]
        )
        for branch in branches
    }
    status = {
        "status": "running",
        "pids": {branch: proc.pid for branch, proc in processes.items()},
    }
    write_json(campaign / "orchestrator_status.json", status)
    codes = {branch: proc.wait() for branch, proc in processes.items()}
    status.update(
        status="completed" if all(code == 0 for code in codes.values()) else "failed",
        returncodes=codes,
    )
    write_json(campaign / "orchestrator_status.json", status)
    return 0 if status["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
