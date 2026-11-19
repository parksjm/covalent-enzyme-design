"""Monitor isolated24/48-image spacing trials using unchanged native capture/export."""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path
import sys
import tempfile
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_progress import parse_last_complete_neb_report  # noqa:E402
from nyl12_covalent.neb_xyz_snapshot import snapshot_complete_band_from_xyz  # noqa:E402
from nyl12_covalent.relay_neb_progress import (  # noqa:E402
    atomic_json,
    digest,
    publish_band,
    record_terminal_failure,
    retain_sources,
    source_paths,
)
from scripts.monitor_acylenzyme_candidates import find_sander_pid  # noqa:E402
from scripts.record_agent_action import append_action  # noqa:E402

LOGGER = logging.getLogger(__name__)
TRANSIENT = (
    "No complete NEB report found",
    "Latest XYZ band has no closed native save loop",
    "Incomplete ordered XYZ/native save loop",
    "Incomplete native XYZ final record",
    "Incomplete native XYZ frame",
    "XYZ snapshot inputs changed during byte capture",
    "XYZ data or writer changed during byte capture",
    "XYZ inputs changed before byte capture closed",
)


def audit(
    journal: Path, phase: str, action: str, files: list[Path] | None = None
) -> None:
    """Durably record monitor intent/outcome before dependent mutations.

    Args:
        journal: Active append-only action journal.
        phase: Required intent, outcome or note phase.
        action: Nonsecret operational description.
        files: Optional explicit files whose current hashes accompany the record.
    """
    append_action(
        journal,
        "neb_spacing_trial_monitor",
        phase,
        action,
        "tool",
        "Authorized bounded coherent native NEB progress export; no new calculations",
        "Preserve calculation files and immutable progress; monitor never signals processes",
        files=[str(path) for path in files] if files else None,
    )


def writer_state(campaign: Path, expected_executable: Path) -> tuple[int | None, bool]:
    """Identify the unique exact-cwd/executable writer or explicit terminal job.

    Args:
        campaign: NEB stage root containing calculation/ and job/status.json.
        expected_executable: Exact case Amber serial Sander executable.

    Returns:
        Live verified PID, if any, and launcher terminal status.
    """
    pid = find_sander_pid(campaign / "calculation")
    if pid is not None and (Path("/proc") / str(pid) / "exe").resolve(
        strict=True
    ) != expected_executable.resolve(strict=True):
        raise ValueError("NEB cwd writer executable differs from declared case Sander")
    path = campaign / "job/status.json"
    terminal = False
    if path.is_file():
        status = json.loads(path.read_text())
        terminal = (
            status["ended_utc"] is not None and type(status["exit_status"]) is int
        )
    if terminal and pid is not None:
        raise RuntimeError("Terminal launcher status conflicts with a live NEB writer")
    return pid, terminal


def deferred(output: Path, reason: str, terminal: bool) -> dict[str, Any]:
    """Preserve last coherent latest while publishing an explicit waiting reason."""
    record = {
        "terminal": terminal,
        "outcome": reason,
        "latest_preserved": (output / "latest").is_symlink(),
    }
    atomic_json(output / "monitor_status.json", record)
    return record


def refresh(config_path: Path, output: Path, journal: Path) -> dict[str, Any]:
    """Publish at most one new complete band; defer only recognized writer races.

    Args:
        config_path: Explicit root-prepared NEB config, not the original case.
        output: Progress destination with atomic latest pointer.
        journal: Active durable journal; failure stops dependent mutations.

    Returns:
        Terminal status and one of waiting, unchanged or published outcomes.
    """
    cfg = json.loads(config_path.read_text())
    if cfg["regular_images"] not in {24, 48} or cfg["active_atoms"] != 4662:
        raise ValueError(
            "Expected explicitly declared24 or48 regular images and4662 active atoms"
        )
    case_path = Path(cfg["case_path"])
    case = json.loads(case_path.read_text())
    if (
        cfg["topology"] != case["topology"]
        or case["qm_charge"] != -1
        or len(case["qm_atom_ids"]) != 102
    ):
        raise ValueError("NEB config differs from declared primary QM102/-1 case")
    pid, terminal = writer_state(
        config_path.parent, Path(case["amber_home"]) / "bin/sander"
    )
    output.mkdir(parents=True, exist_ok=True)
    native = config_path.parent / "calculation"
    mdout, xyz = native / "minimized.mdout", native / "nebpath.xyz"
    if pid is None and not terminal:
        return deferred(
            output,
            "waiting for unique declared native writer or terminal job status",
            False,
        )
    if not mdout.is_file() or not xyz.is_file():
        if terminal:
            raise FileNotFoundError(
                "Terminal NEB has no native mdout/XYZ band to export"
            )
        return deferred(output, "waiting for native mdout/XYZ", False)
    try:
        report = parse_last_complete_neb_report(
            mdout.read_bytes(), cfg["regular_images"]
        )
    except ValueError as error:
        if not terminal and str(error) in TRANSIENT:
            return deferred(output, str(error), False)
        raise
    paths = [Path(__file__).resolve(), *source_paths(), config_path, case_path]
    source_identity = {str(path): digest(path) for path in paths}
    identity = {
        "report_index": report.report_index,
        "cycle": report.cycle,
        "energies_hartree": [i.energy_hartree for i in report.images],
        "source_sha256": source_identity,
        "terminal": terminal,
    }
    progress = output / "progress.json"
    if progress.is_file() and json.loads(progress.read_text())["identity"] == identity:
        return deferred(
            output,
            "complete report unchanged; no reconstruction/export repeated",
            terminal,
        )
    with tempfile.TemporaryDirectory(prefix="relay_neb_verified_band_") as temporary:
        band = Path(temporary) / "band"
        try:
            captured = snapshot_complete_band_from_xyz(
                mdout,
                xyz,
                band,
                cfg,
                cfg["regular_images"],
                writer_pid=pid,
                completed_job_status=config_path.parent / "job/status.json"
                if terminal
                else None,
            )
        except FileNotFoundError as error:
            if pid is not None and str(error.filename).startswith(f"/proc/{pid}/"):
                return deferred(
                    output,
                    "verified writer exited during capture; await terminal job status",
                    False,
                )
            raise
        except (ValueError, RuntimeError) as error:
            if not terminal and str(error) in TRANSIENT:
                return deferred(output, str(error), False)
            raise
        # The writer may finish another report between preflight and capture.
        identity.update(
            report_index=captured.report_index,
            cycle=captured.cycle,
            energies_hartree=[i.energy_hartree for i in captured.images],
        )
        milestone = retain_sources(output, source_identity, paths)
        record = publish_band(
            cfg, band, captured, output, identity, terminal, milestone
        )
    record["outcome"] = "published coherent band"
    atomic_json(
        output / "monitor_status.json",
        {k: record[k] for k in ("terminal", "outcome", "cycle", "latest_update")},
    )
    return record


def main(argv: list[str] | None = None) -> None:
    """Run once or poll every30/60s until terminal export or a bounded deadline.

    Args:
        argv: Optional explicit CLI arguments for command-line tests.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--journal", type=Path, default=ROOT / "output_agent_actions/20261021.jsonl"
    )
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=int, choices=(30, 60), default=30)
    parser.add_argument("--hours", type=float, default=48)
    args = parser.parse_args(argv)
    if not math.isfinite(args.hours) or args.hours <= 0:
        parser.error("--hours must be finite and positive")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    output = (
        args.output.resolve()
        if args.output
        else args.config.resolve().parent / "progress"
    )
    deadline = time.monotonic() + args.hours * 3600
    while True:
        audit(
            args.journal,
            "intent",
            f"Inspect coherent relay NEB band: {args.config}",
            [args.config],
        )
        try:
            result = refresh(args.config.resolve(), output, args.journal)
        except Exception as error:
            audit(
                args.journal,
                "intent",
                "Check/write compact terminal failure acceptance record",
            )
            record_terminal_failure(args.config.resolve(), error)
            audit(
                args.journal,
                "outcome",
                f"NEB progress failed: {type(error).__name__}: {error}",
            )
            raise
        audit(
            args.journal,
            "outcome",
            f"NEB progress: {result['outcome']}; terminal={result['terminal']}",
            [output / "monitor_status.json"],
        )
        LOGGER.info("%s; terminal=%s", result["outcome"], result["terminal"])
        if result["terminal"] or not args.watch:
            return
        if time.monotonic() >= deadline:
            raise TimeoutError("Relay NEB monitor reached bounded deadline")
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
