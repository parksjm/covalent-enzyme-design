"""Validated source checks and publication helpers for direct PMF monitoring."""

from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

from nyl12_covalent.direct_us_continuation_io import (
    _trajectory_digest,
    validate_accepted_stage,
    validate_continuation_plan,
)
from nyl12_covalent.relay_umbrella_analysis_io import digest

LOGGER = logging.getLogger(__name__)
HORIZONS_PS = tuple(float(value) for value in range(35, 61, 5))
POLL_SECONDS = 120


def _repository_root(campaign: Path) -> Path:
    """Find the repository owning a campaign and its active action journal."""
    for parent in (campaign, *campaign.parents):
        if (parent / "scripts/record_agent_action.py").is_file() and (
            parent / "output_agent_actions/20261021.jsonl"
        ).parent.is_dir():
            return parent
    raise FileNotFoundError("Cannot locate project action journal from campaign path")


def _record_action(
    root: Path, phase: str, action: str, reason: str, file: Path, result: str = ""
) -> None:
    """Append one watcher action to the project journal."""
    command = [
        sys.executable,
        str(root / "scripts/record_agent_action.py"),
        "--log",
        str(root / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "direct_fixed_pmf_watcher",
        "--phase",
        phase,
        "--action",
        action,
        "--permission-mode",
        "escalated",
        "--reason",
        reason,
        "--rollback",
        "Retain derived evidence; native jobs are untouched",
        "--file",
        str(file),
    ]
    if result:
        command.extend(["--result", result])
    subprocess.run(command, check=True)


def _relay_request(campaign: Path) -> Path:
    """Recover the pinned relay analysis request from an existing direct report."""
    analysis = campaign / "analysis/fixed_burnin_20261126"
    for report in (
        analysis / "direct_30ps_analysis.json",
        analysis / "direct_25ps_analysis.json",
    ):
        if not report.is_file():
            continue
        data = json.loads(report.read_text())
        matches = [
            (Path(filename), expected)
            for filename, expected in data.get("source_sha256", {}).items()
            if Path(filename).name == "analysis_request.json"
        ]
        if len(matches) == 1 and digest(matches[0][0]) == matches[0][1]:
            return matches[0][0]
    raise ValueError(
        "No uniquely hash-verified relay analysis request in fixed reports"
    )


def _continuation_ready(
    campaign: Path, plan: Path, horizon: float
) -> tuple[bool, dict[str, list[dict[str, Any]]], dict[str, str], Path]:
    """Check whether all requested windows have accepted native and derived blocks."""
    plan_data = json.loads(plan.read_text())
    augmented = Path(plan_data["augmentation_root"]).resolve()
    records, sources = validate_continuation_plan(campaign, plan, horizon, augmented)
    for window_records in records.values():
        for record in window_records:
            stage = Path(record["stage"])
            derived = augmented / record["window"] / record["label"]
            required = [
                stage / "process.json",
                stage / "acceptance.json",
                stage / "observables/observables.json",
                stage / "observables/observables.csv",
                derived / "augmentation.json",
                derived / "observables.csv",
            ]
            if not all(path.is_file() for path in required):
                return False, records, sources, augmented
            sources.update(validate_accepted_stage(record, augmented))
    return True, records, sources, augmented


def _verify_report(
    report_path: Path,
    horizon_ps: float,
    required_sources: dict[str, str] | None = None,
    *,
    expected_windows: int = 63,
) -> dict[str, Any]:
    """Reject incomplete, mismatched, or stale fixed-PMF evidence."""
    report = json.loads(report_path.read_text())
    if (
        float(report.get("horizon_ps", -1)) != horizon_ps
        or report.get("windows") != expected_windows
        or report.get("protocol")
        != "fixed5ps_once_per_continuous_production_window_all_frames"
    ):
        raise ValueError(
            f"Saved report has the wrong fixed-analysis identity: {report_path}"
        )
    source_hashes = report.get("source_sha256")
    if not isinstance(source_hashes, dict) or not source_hashes:
        raise ValueError(f"Saved report lacks source hashes: {report_path}")
    for filename, expected in source_hashes.items():
        path = Path(filename)
        if not path.is_file():
            raise ValueError(f"Saved report source changed or is missing: {path}")
        actual = (
            _trajectory_digest(path) if path.name == "trajectory.nc" else digest(path)
        )
        if actual != expected:
            raise ValueError(f"Saved report source changed or is missing: {path}")
    for filename, expected in (required_sources or {}).items():
        if source_hashes.get(filename) != expected:
            raise ValueError(
                f"Saved report does not include current continuation source: {filename}"
            )
    return report


def _journaled_child(
    root: Path,
    evidence: Path,
    action: str,
    reason: str,
    command: list[str],
    log_path: Path,
    output_path: Path,
    timeout_seconds: float,
    source_root: Path,
) -> None:
    """Run one derived child, retaining logs/status and never retrying failures."""
    _record_action(root, "intent", action, reason, output_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_path.open("x") as stream:
            environment = os.environ.copy()
            package_root = str(source_root / "src")
            existing_pythonpath = environment.get("PYTHONPATH")
            environment["PYTHONPATH"] = os.pathsep.join(
                [package_root, *([existing_pythonpath] if existing_pythonpath else [])]
            )
            result = subprocess.run(
                command,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
                timeout=max(1.0, timeout_seconds),
                env=environment,
            )
    except (OSError, subprocess.TimeoutExpired) as error:
        status_path = evidence / f"{action.replace(' ', '_')}_failure.json"
        status_path.write_text(
            json.dumps(
                {"action": action, "error": str(error), "log": str(log_path)}, indent=2
            )
            + "\n"
        )
        _record_action(root, "outcome", action, reason, status_path, f"Failed: {error}")
        raise
    if result.returncode != 0:
        status_path = evidence / f"{action.replace(' ', '_')}_failure.json"
        status_path.write_text(
            json.dumps(
                {
                    "action": action,
                    "returncode": result.returncode,
                    "log": str(log_path),
                    "output": str(output_path),
                },
                indent=2,
            )
            + "\n"
        )
        _record_action(
            root,
            "outcome",
            action,
            reason,
            status_path,
            f"Child exited {result.returncode}",
        )
        raise subprocess.CalledProcessError(result.returncode, command)
    _record_action(
        root, "outcome", action, reason, output_path, f"Completed; log={log_path}"
    )


def _plot_metadata(output: Path, direct: Path, relay: Path) -> dict[str, Any] | None:
    """Return an existing plot identity only after validating both report hashes."""
    direct_h = float(json.loads(direct.read_text())["horizon_ps"])
    relay_h = float(json.loads(relay.read_text())["horizon_ps"])
    candidates = sorted(
        output.glob(f"comparison_{direct_h:g}ps_vs_relay_{relay_h:g}ps_*_source.json")
    )
    if not candidates:
        return None
    if len(candidates) != 1:
        raise ValueError(f"Expected one plot metadata file in {output}")
    metadata = json.loads(candidates[0].read_text())
    if (
        metadata.get("direct_report_sha256") != digest(direct)
        or metadata.get("relay_report_sha256") != digest(relay)
        or metadata.get("horizons_matched") is not True
    ):
        raise ValueError(f"Existing matched plot metadata is stale: {candidates[0]}")
    paths = metadata.get("plots", [])
    if not paths or any(not Path(path).is_file() for path in paths):
        raise ValueError(
            f"Existing plot metadata points to missing figures: {candidates[0]}"
        )
    return metadata


def _replace_latest(plots: Path, metadata: dict[str, Any]) -> dict[str, Any]:
    """Atomically replace only the known latest PNG/PDF symlinks, preserving targets."""
    files = [Path(path) for path in metadata["plots"]]
    aliases = {}
    for extension, source in zip(("png", "pdf"), files[:2], strict=True):
        alias = plots / f"pmf_latest.{extension}"
        old = None
        if alias.exists() or alias.is_symlink():
            if not alias.is_symlink():
                raise ValueError(
                    f"Refusing to replace non-symlink latest alias: {alias}"
                )
            old = os.readlink(alias)
            old_path = (alias.parent / old).resolve(strict=True)
            if plots.resolve() not in old_path.parents:
                raise ValueError(
                    f"Latest alias resolves outside campaign plots: {alias}"
                )
        target = source.resolve(strict=True)
        if plots.resolve() not in target.parents:
            raise ValueError(f"New latest target is outside campaign plots: {source}")
        pending = alias.with_name(alias.name + f".{os.getpid()}.pending")
        if pending.exists() or pending.is_symlink():
            raise FileExistsError(pending)
        os.symlink(os.path.relpath(target, alias.parent), pending)
        os.replace(pending, alias)
        aliases[str(alias)] = {
            "previous_target": old,
            "previous_sha256": hashlib.sha256(old_path.read_bytes()).hexdigest()
            if old
            else None,
            "new_target": os.readlink(alias),
            "new_sha256": digest(target),
        }
    return aliases


def _write_status(path: Path, value: dict[str, Any]) -> None:
    """Atomically retain watcher progress without replacing scientific reports."""
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(value, indent=2) + "\n")
    pending.replace(path)
