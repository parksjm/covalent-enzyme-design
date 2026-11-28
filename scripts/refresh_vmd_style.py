#!/usr/bin/env python3
"""Refresh explicit curated VMD display styles without regenerating trajectories."""

import argparse
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import subprocess
import sys

SOURCE_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("NYL12_VMD_STYLE_REPOSITORY_ROOT", str(SOURCE_ROOT)))
sys.path.insert(0, str(SOURCE_ROOT / "src"))

from nyl12_covalent.vmd_style_manifests import (  # noqa: E402
    backup_file,
    checksum_updates,
    digest,
    verify_checksums,
)
from nyl12_covalent.vmd_style_refresh import prepare_correction  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Validate a JSON refresh plan and optionally apply reversible corrections."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--backup", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    rows = plan["movies"]
    report = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        plan=str(args.plan),
        applied=args.apply,
        movies=[],
        failures=[],
    )
    prepared = []
    # Validate every physical mapping before the first existing viewer changes.
    for row in rows:
        try:
            path, text = prepare_correction(row)
            prepared.append((row, path, text))
        except (ValueError, OSError, KeyError) as exc:
            report["failures"].append(dict(viewer=row.get("viewer"), error=str(exc)))
    if report["failures"]:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        raise ValueError(f"Mapping validation failed; see {args.report}")
    for row, viewer, text in prepared:
        changed = {viewer: text}
        changed.update(checksum_updates(viewer.parent, {viewer.name: text}))
        entry = dict(
            viewer=str(viewer),
            qm_count=len(row["qm_indices"]),
            qm_evidence=row["qm_evidence"],
            topology_validation=row.get("topology_validation"),
            before_checksum_failures=verify_checksums(viewer.parent),
            changes=[],
        )
        if args.apply:
            action = f"Refresh display-only viewer and derived checksums: {viewer}"
            command = [
                sys.executable,
                str(ROOT / "scripts/record_agent_action.py"),
                "--log",
                str(ROOT / "output_agent_actions/20261021.jsonl"),
                "--agent",
                "vmd_style_refresh",
                "--permission-mode",
                "escalated",
                "--reason",
                "Explicit curated movie style refresh; preserve scientific files",
                "--rollback",
                f"Restore only backed-up derived files from {args.backup}",
                "--action",
                action,
            ]
            subprocess.run(
                [*command, "--phase", "intent"], check=True, capture_output=True
            )
            for path, replacement in changed.items():
                if path.read_text() == replacement:
                    continue
                preserved = backup_file(path, ROOT, args.backup)
                path.write_text(replacement)
                preserved["after_sha256"] = digest(path)
                entry["changes"].append(preserved)
            subprocess.run(
                [
                    *command,
                    "--phase",
                    "outcome",
                    "--result",
                    json.dumps(entry["changes"]),
                ],
                check=True,
                capture_output=True,
            )
        entry["checksum_failures"] = (
            verify_checksums(viewer.parent) if args.apply else []
        )
        report["movies"].append(entry)
        LOGGER.info("%s %s", "Refreshed" if args.apply else "Validated", viewer)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    if any(row["checksum_failures"] for row in report["movies"]):
        raise ValueError(f"Portable checksum mismatches remain; see {args.report}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
