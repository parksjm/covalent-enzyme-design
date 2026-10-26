#!/usr/bin/env python3
"""Archive redundant QM/MM diagnostics after recording their identities."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import shutil
from typing import Any


ROOT = Path("/data/jpc/covalent-enzyme-design")
OLD_RUN = ROOT / "output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2"
REPORT = ROOT / "provenance/qmmm_cleanup_20261025/cleanup_manifest.json"
ARCHIVE = ROOT / "clutter/obsolete_qmmm_20261025"


def digest(path: Path) -> str:
    """Return a streaming SHA-256 digest for a potentially large file."""
    value = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            value.update(block)
    return value.hexdigest()


def candidates() -> list[Path]:
    """Return only redundant, regenerable files from superseded runs."""
    paths = list(OLD_RUN.rglob("path*.xyz"))
    paths.extend(OLD_RUN.rglob("optimization.nc"))
    for directory in (
        OLD_RUN / "live_energy_plot",
        OLD_RUN / "live_inspection",
        OLD_RUN / "vmd_progress_snapshot_20261025",
    ):
        if directory.exists():
            paths.extend(path for path in directory.rglob("*") if path.is_file())
    return sorted(set(path.resolve() for path in paths))


def write_report(payload: dict[str, Any]) -> None:
    """Persist the cleanup report after each file."""
    REPORT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def main() -> None:
    """Hash and archive exact redundant files while retaining core run records."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    selected = candidates()
    if not selected:
        raise FileNotFoundError(
            "No obsolete artifacts found; refusing ambiguous cleanup"
        )
    for path in selected:
        if OLD_RUN.resolve() not in path.parents:
            raise ValueError(f"Cleanup target escapes the obsolete run: {path}")
    report: dict[str, Any] = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "policy": (
            "Archive only redundant DL-FIND full-system XYZ files, superseded full "
            "optimizer NetCDF files, and regenerated live-publisher products. "
            "Retain restarts, mdout, inputs, acceptance, source snapshots, curated "
            "VMD, and append-only provenance."
        ),
        "archive": str(ARCHIVE),
        "rollback": "Move archived files back to their recorded original paths.",
        "files": [],
    }
    write_report(report)
    for index, path in enumerate(selected, start=1):
        record = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": digest(path),
            "archived": False,
        }
        report["files"].append(record)
        write_report(report)
        destination = ARCHIVE / path.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(path, destination)
        record["archive_path"] = str(destination)
        record["archived"] = True
        write_report(report)
        logging.info("Archived %d/%d: %s", index, len(selected), path)
    report["ended_utc"] = datetime.now(timezone.utc).isoformat()
    report["total_bytes"] = sum(item["bytes"] for item in report["files"])
    write_report(report)


if __name__ == "__main__":
    main()
