#!/usr/bin/env python3
"""Bounded local watcher for compact accepted blocks synced onto Violet."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import tempfile
import time

from nyl12_covalent.relay_umbrella_analysis_io import analyze_request
from nyl12_covalent.relay_umbrella_analysis_reports import compare_reports, plot_report
from nyl12_covalent.relay_mbar_progression import plot_progression, select_horizon

ROOT = Path(__file__).resolve().parents[1]


def journal(phase: str, action: str, result: str = "", *, root: Path = ROOT) -> None:
    """Record watcher intent/outcome; logging failure stops dependent mutations."""
    subprocess.run(
        [
            "python3",
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_mbar_watcher",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "First full accepted five-ps coverage and subsequent block checkpoints",
            "--rollback",
            "Retain immutable input and report artifacts; stop only verified watcher process",
            "--result",
            result,
        ],
        check=True,
    )


def discover(
    request: dict, receipts_root: Path, pattern: str
) -> tuple[list[Path], bool]:
    """Read only manifest windows and direct production stages.

    Targeted diagnostics and derived trees are outside this campaign's source
    namespace. A mislabeled receipt inside a declared window remains an error.
    """
    receipts = sorted(
        path
        for window in request["windows"]
        for path in (receipts_root / window["label"]).glob(
            f"*_production_5ps/{pattern}"
        )
    )
    duration = dict.fromkeys((window["label"] for window in request["windows"]), 0.0)
    for path in receipts:
        proof = json.loads(path.read_text())
        label = proof["window"]
        if label not in duration:
            raise ValueError(
                f"Unexpected accepted-block window in watched tree: {label}"
            )
        if proof["accepted"] is not True:
            raise ValueError(
                f"Nonaccepted proof must not use accepted receipt name: {path}"
            )
        duration[label] += float(proof["sampled_duration_ps"])
    return receipts, bool(duration) and min(duration.values()) >= 5.0


def main() -> None:
    """Watch a bounded local compact tree; never launch native MD or sync data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    parser.add_argument("--receipts-root", required=True, type=Path)
    parser.add_argument("--receipt-pattern", default="analysis_acceptance.json")
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--max-minutes", required=True, type=float)
    parser.add_argument("--poll-seconds", default=60.0, type=float)
    parser.add_argument("--journal-root", type=Path, default=ROOT)
    parser.add_argument("--resume-output", action="store_true")
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    if args.max_minutes <= 0 or args.poll_seconds < 30:
        raise ValueError(
            "Positive runtime bound and poll interval at least 30 seconds required"
        )
    journal(
        "intent", "Start bounded local compact-block watcher", root=args.journal_root
    )
    try:
        watch(args)
    except Exception as exc:
        journal(
            "outcome",
            "Bounded watcher failed",
            f"{type(exc).__name__}: {exc}",
            root=args.journal_root,
        )
        if args.output_root.is_dir():
            (args.output_root / "watcher_error.json").write_text(
                json.dumps(
                    dict(
                        exception_type=type(exc).__name__,
                        message=str(exc),
                        process_completed=False,
                    ),
                    indent=2,
                )
                + "\n"
            )
        raise
    journal(
        "outcome",
        "Bounded watcher completed without launching MD",
        root=args.journal_root,
    )


def watch(args: argparse.Namespace) -> None:
    """Monitor local accepted blocks within the explicitly bounded request."""
    cfg = json.loads(args.request.read_text())
    args.output_root.mkdir(parents=True, exist_ok=args.resume_output)
    deadline = time.monotonic() + args.max_minutes * 60
    previous, readiness_identity = None, None
    checkpoints: dict[float, dict] = {}
    readiness_count = max(
        (
            int(path.stem.rsplit("_", 1)[1])
            for path in args.output_root.parent.glob(
                "readiness_update_[0-9][0-9][0-9].json"
            )
        ),
        default=0,
    )
    while time.monotonic() < deadline:
        if (args.receipts_root / ".compact_sync_in_progress").exists():
            logging.info(
                "Compact sync marker present; waiting for complete atomic snapshot"
            )
            if args.once:
                break
            time.sleep(min(args.poll_seconds, max(0, deadline - time.monotonic())))
            continue
        receipts, _ = discover(
            cfg, args.receipts_root.resolve(strict=True), args.receipt_pattern
        )
        identity = hashlib.sha256(
            "".join(str(p) + p.read_text() for p in receipts).encode()
        ).hexdigest()
        durations = dict.fromkeys((window["label"] for window in cfg["windows"]), 0.0)
        for receipt in receipts:
            proof = json.loads(receipt.read_text())
            durations[proof["window"]] += float(proof["sampled_duration_ps"])
        common_duration = min(durations.values(), default=0.0)
        maximum_horizon = int((common_duration + 1e-8) // 5.0) * 5
        if maximum_horizon < 5 and identity != readiness_identity:
            journal(
                "intent",
                "Refresh partial-coverage readiness without MBAR",
                identity,
                root=args.journal_root,
            )
            readiness_count += 1
            prefix = args.output_root.parent / f"readiness_update_{readiness_count:03d}"
            request_path = prefix.with_suffix(".request.json")
            request_path.write_text(
                json.dumps(
                    dict(
                        cfg,
                        accepted_block_receipts=[str(p.resolve()) for p in receipts],
                    ),
                    indent=2,
                )
                + "\n"
            )
            readiness = analyze_request(request_path)
            prefix.with_suffix(".json").write_text(
                json.dumps(readiness, indent=2, allow_nan=False) + "\n"
            )
            plot_report(readiness, prefix.with_suffix(".png"))
            readiness_identity = identity
            journal(
                "outcome",
                "Refreshed partial-coverage readiness without MBAR",
                readiness["status"],
                root=args.journal_root,
            )
        checkpoint_changed = False
        for horizon in (float(value) for value in range(5, maximum_horizon + 1, 5)):
            selection = select_horizon(receipts, cfg["windows"], horizon)
            if selection is None:
                continue
            selected, durations = selection
            checkpoint_id = hashlib.sha256(
                "".join(str(path) + path.read_text() for path in selected).encode()
            ).hexdigest()
            checkpoint = args.output_root / f"checkpoint_{int(horizon):02d}ps"
            if checkpoint.exists():
                source = json.loads((checkpoint / "source.json").read_text())
                if source["identity_sha256"] != checkpoint_id:
                    raise ValueError(
                        f"Cumulative checkpoint inputs changed: {checkpoint}"
                    )
                checkpoints[horizon] = json.loads(
                    (checkpoint / "analysis.json").read_text()
                )
                if checkpoints[horizon].get("pmf"):
                    previous = checkpoints[horizon]
                continue
            journal(
                "intent",
                f"Analyze cumulative {horizon:g} ps coverage",
                checkpoint_id,
                root=args.journal_root,
            )
            checkpoint.mkdir()
            request = dict(
                cfg, accepted_block_receipts=[str(path.resolve()) for path in selected]
            )
            request_path = checkpoint / "analysis_request.json"
            request_path.write_text(json.dumps(request, indent=2) + "\n")
            result = analyze_request(request_path)
            if previous is not None and result.get("pmf") and previous.get("pmf"):
                result["gates"]["block_stability"] = compare_reports(result, previous)
            result["checkpoint"] = {
                "actual_horizon_ps": horizon,
                "per_window_duration_ps": durations,
                "scope": "cumulative matched production data through this horizon",
                "plot_common_es_reference_s": next(
                    (
                        prior["pmf"]["reference_s"]
                        for _, prior in sorted(checkpoints.items())
                        if prior.get("status") == "preliminary_connected_not_converged"
                        and prior.get("pmf")
                    ),
                    result["pmf"].get("reference_s") if result.get("pmf") else None,
                ),
                "plot_normalization": "Each connected curve is shifted to the earliest connected checkpoint's sampled E-S reference bin; no uncertainty is shown",
                "converged": False,
            }
            (checkpoint / "analysis.json").write_text(
                json.dumps(result, indent=2, allow_nan=False) + "\n"
            )
            (checkpoint / "source.json").write_text(
                json.dumps(
                    {
                        "identity_sha256": checkpoint_id,
                        "accepted_receipts": [str(path.resolve()) for path in selected],
                        "receipt_sha256": {
                            str(path.resolve()): hashlib.sha256(
                                path.read_bytes()
                            ).hexdigest()
                            for path in selected
                        },
                        "actual_horizon_ps": horizon,
                        "per_window_duration_ps": durations,
                        "physical_hamiltonian_sha256": cfg[
                            "physical_hamiltonian_sha256"
                        ],
                        "es_search_bounds": cfg["es_search_bounds"],
                    },
                    indent=2,
                )
                + "\n"
            )
            plot_report(result, checkpoint / "analysis.png")
            checkpoints[horizon] = result
            previous = result if result.get("pmf") else previous
            plot_progression(
                sorted(checkpoints.items()), checkpoint / "progression.png"
            )
            checkpoint_changed = True
            journal(
                "outcome",
                f"Analyzed cumulative {horizon:g} ps coverage",
                result["status"],
                root=args.journal_root,
            )
        if checkpoint_changed:
            progression = sorted(checkpoints.items())
            temporary = args.output_root / ".latest_progression.tmp.png"
            temporary.unlink(missing_ok=True)
            with tempfile.NamedTemporaryFile(
                dir=args.output_root, suffix=".png", delete=False
            ) as temp:
                temporary = Path(temp.name)
            try:
                temporary.unlink()
                plot_progression(progression, temporary)
                os.replace(temporary, args.output_root / "latest_progression.png")
            finally:
                temporary.unlink(missing_ok=True)
            latest = max(checkpoints)
            latest_tmp = args.output_root / "latest_checkpoint.json.tmp"
            latest_result = checkpoints[latest]
            latest_tmp.write_text(
                json.dumps(
                    {
                        "actual_horizon_ps": latest,
                        "status": latest_result["status"],
                        "source_identity_sha256": json.loads(
                            (
                                args.output_root
                                / f"checkpoint_{int(latest):02d}ps"
                                / "source.json"
                            ).read_text()
                        )["identity_sha256"],
                        "analysis": f"checkpoint_{int(latest):02d}ps/analysis.json",
                        "figure": "latest_progression.png",
                        "converged": False,
                    },
                    indent=2,
                )
                + "\n"
            )
            os.replace(latest_tmp, args.output_root / "latest_checkpoint.json")
        else:
            logging.info(
                "Readiness: %d accepted receipts; full coverage=%s",
                len(receipts),
                bool(checkpoints),
            )
        if args.once:
            break
        time.sleep(min(args.poll_seconds, max(0, deadline - time.monotonic())))
    (args.output_root / "watcher_completion.json").write_text(
        json.dumps(
            dict(
                process_completed=True,
                reason="once" if args.once else "runtime_bound_reached",
                full_coverage_analyzed=previous is not None,
                converged=False,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
