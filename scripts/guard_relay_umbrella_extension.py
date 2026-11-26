#!/usr/bin/env python3
"""Wait for the accepted 30-ps baseline, then launch the prepared 60-ps extension."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import json
import logging
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.jobs import launch_job  # noqa: E402
from nyl12_covalent.neb_umbrella_campaign import digest, write_json  # noqa: E402
from nyl12_covalent.relay_umbrella_analysis_io import load_block  # noqa: E402
from nyl12_covalent.relay_umbrella_runtime import verify_plan  # noqa: E402


def readiness(plan: dict, config: dict, campaign: Path) -> dict:
    """Check all baseline endpoints, lineage and strictly accepted sampling data.

    Args:
        plan: Immutable extension execution plan.
        config: Extended prepared configuration, with original stage identities.
        campaign: Existing campaign directory containing native outputs.

    Returns:
        Explicit readiness and missing-window list. Existing invalid or failed
        artifacts raise; absent future baseline outputs mean wait.
    """
    original = Path(plan["original_plan"])
    if digest(original) != plan["original_plan_sha256"]:
        raise ValueError("Original baseline plan changed")
    request = json.loads((campaign / "analysis_request.json").read_text())
    pending = []
    for row in config["windows"]:
        window = campaign / row["label"]
        status_path = window / "window_status.json"
        if (
            status_path.exists()
            and json.loads(status_path.read_text())["status"] == "failed"
        ):
            raise ValueError(f"Baseline window failed: {row['label']}")
        stage = window / "08_production_5ps"
        required = [
            stage / name
            for name in [
                "acceptance.json",
                "process.json",
                "restart.ncrst",
                "input_identity.json",
                "analysis_acceptance.json",
            ]
        ]
        segment = window / "segment_receipt_limit_06.json"
        if not all(path.is_file() for path in [*required, segment]):
            pending.append(row["label"])
            continue
        process = json.loads((stage / "process.json").read_text())
        if process["exit_status"] != 0:
            raise ValueError(f"Baseline native process failed: {row['label']}")
        receipt = json.loads(segment.read_text())
        if receipt["status"] != "accepted" or receipt["production_blocks"] != 6:
            raise ValueError("Baseline does not contain six accepted blocks")
        if receipt["checkpoint_sha256"] != digest(stage / "restart.ncrst"):
            raise ValueError("Baseline endpoint restart changed")
        identity = json.loads((stage / "input_identity.json").read_text())["sha256"]
        inputs = [
            original,
            campaign / config["topology"],
            window / "07_production_5ps/restart.ncrst",
            stage / "run.mdin",
            stage / "umbrella.RST",
        ]
        if identity != {str(path): digest(path) for path in inputs}:
            raise ValueError("Baseline final-stage input lineage mismatch")
        load_block(
            stage / "analysis_acceptance.json",
            temperature=config["temperature_K"],
            hamiltonian=config["physical_hamiltonian_sha256"],
            diagnostic_columns=request["diagnostic_columns"],
        )
    return dict(
        ready=not pending,
        windows_ready=len(config["windows"]) - len(pending),
        windows_total=len(config["windows"]),
        pending_windows=pending,
    )


def native_active(plan: dict, campaign: Path) -> list[int]:
    """Return verified Sander PIDs in this campaign, without signaling processes.

    Args:
        plan: Execution plan identifying the exact native executable.
        campaign: Working-directory ancestry identifying this campaign.

    Returns:
        Matching process IDs; unrelated scientific jobs are ignored.
    """
    active = []
    for command in Path("/proc").glob("[0-9]*/cmdline"):
        try:
            argv = command.read_bytes().split(b"\0")
            cwd = (command.parent / "cwd").resolve(strict=True)
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        if (
            argv
            and argv[0].decode(errors="replace") == plan["sander"]
            and cwd.is_relative_to(campaign)
        ):
            active.append(int(command.parent.name))
    return active


def main() -> None:
    """Launch only after root review, complete native lineage and zero overlapping workers."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--review", required=True, type=Path)
    parser.add_argument("--max-hours", required=True, type=float)
    args = parser.parse_args()
    if args.max_hours <= 0:
        raise ValueError("Positive guard duration required")
    plan, config, campaign = verify_plan(args.plan)
    review = json.loads(args.review.read_text())
    if review.get("accepted") is not True or review.get("plan_sha256") != digest(
        args.plan
    ):
        raise ValueError("Explicit reviewed extension plan required")
    source = args.plan.parent / "source_snapshot"
    status = args.plan.parent / "gate_status.json"
    journal = [
        sys.executable,
        str(Path("/data/jpc/covalent-enzyme-design/scripts/record_agent_action.py")),
        "--log",
        "/data/jpc/covalent-enzyme-design/output_agent_actions/20261021.jsonl",
        "--agent",
        "relay_extension_gate",
        "--permission-mode",
        "escalated",
        "--reason",
        "PIauthorizes30to60pscontinuation;57acceptedbaselineendpointsrequired",
        "--rollback",
        "Keepalloutputs;stoponlyextensionworkersafterverifyingtheiridentity",
    ]
    with (args.plan.parent / "extension_gate.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        deadline = time.monotonic() + args.max_hours * 3600
        while time.monotonic() < deadline:
            check = readiness(plan, config, campaign)
            active = native_active(plan, campaign) if check["ready"] else []
            write_json(
                status,
                dict(
                    **check,
                    native_campaign_pids=active,
                    status="waiting_for_baseline",
                    native_launched=False,
                    checked_utc=datetime.now(timezone.utc).isoformat(),
                ),
            )
            if check["ready"] and not active:
                break
            time.sleep(120)
        else:
            raise TimeoutError("Baseline did not complete within extension gate budget")
        if (args.plan.parent / "launch_claim").exists():
            raise FileExistsError(
                "Extension has already been claimed; review existing worker jobs"
            )
        subprocess.run(
            [
                *journal,
                "--phase",
                "intent",
                "--action",
                "Launch reviewed60psrelayextension",
            ],
            check=True,
        )
        subprocess.run(
            [
                plan["python"],
                str(source / "scripts/run_relay_umbrella.py"),
                "launch",
                "--plan",
                str(args.plan),
            ],
            check=True,
            env=None,
        )
        launch_job(
            "relay_us_cont60_diagnostics_20261126",
            args.plan.parent / "diagnostics_job",
            campaign,
            [
                plan["python"],
                str(source / "scripts/guard_relay_umbrella.py"),
                "guard",
                "--campaign",
                str(campaign),
                "--plan",
                str(args.plan),
                "--output",
                str(args.plan.parent / "derived_diagnostics"),
                "--duration-hours",
                "72",
                "--analysis-only",
            ],
            env_overrides=plan["environment"],
            input_paths=[args.plan, args.review],
        )
        write_json(
            status,
            dict(
                status="extension_workers_launched",
                workers=plan["workers"],
                native_launched=True,
                target_ps_per_window=60,
                sampling_converged=False,
            ),
        )
        subprocess.run(
            [
                *journal,
                "--phase",
                "outcome",
                "--action",
                "Launch reviewed60psrelayextension",
            ],
            check=True,
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
