#!/usr/bin/env python3
"""Archive phase-II calculation directories and reset a validated run."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.unrestrained_relay import (  # noqa: E402
    PHASE2_CORRECTION_REASON,
)

LOGGER = logging.getLogger(__name__)
STAGE1 = [f"stage1_{index:02d}" for index in range(21)]


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a regular file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    """Read a JSON object, rejecting missing or malformed inputs."""
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def validate_run(run_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate completed phase I and its accepted reactive geometry."""
    status = _read_json(run_dir / "status.json")
    completed = status.get("completed_windows")
    if not isinstance(completed, list) or completed[: len(STAGE1)] != STAGE1:
        raise ValueError("status.json must begin with stage1_00 through stage1_20")
    if any(not name.startswith("stage2_") for name in completed[len(STAGE1) :]):
        raise ValueError("Only phase-II windows may follow completed phase I")
    missing = [name for name in STAGE1 if not (run_dir / name).is_dir()]
    if missing:
        raise ValueError(f"Missing completed phase-I directories: {missing}")
    acceptance = _read_json(run_dir / "stage1_20" / "acceptance.json")
    try:
        distances = acceptance["distances_A"]
        og_hg = float(distances["OG1_HG1"])
        asp_hg = float(distances["Asp272_OD1_HG1"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "stage1_20 acceptance lacks numeric required distances"
        ) from error
    if not og_hg > 1.3 or not asp_hg < 1.2:
        raise ValueError("stage1_20 acceptance does not meet the geometry gates")
    return status, acceptance


def _file_inventory(directory: Path) -> dict[str, str]:
    """Return relative file paths and hashes below a directory."""
    return {
        path.relative_to(directory).as_posix(): sha256(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def archive_phase2(run_dir: Path, clutter_path: Path, dry_run: bool = False) -> Path:
    """Archive all stage2 directories and reset run status to phase I.

    Args:
        run_dir: Calculation run containing status.json and stage directories.
        clutter_path: User-selected parent directory for the archive.
        dry_run: Validate and report planned changes without writing files.

    Returns:
        The archive directory path, including when running in dry-run mode.
    """
    run_dir = run_dir.resolve()
    clutter_path = clutter_path.resolve()
    status, _ = validate_run(run_dir)
    stages = sorted(path for path in run_dir.glob("stage2_*") if path.is_dir())
    archive_dir = clutter_path / f"{run_dir.name}_stage2_archive"
    if archive_dir.exists():
        raise FileExistsError(f"Archive destination already exists: {archive_dir}")
    if not stages:
        raise ValueError("No invalid stage2_* directories were found")
    LOGGER.info("Found %d phase-II directories", len(stages))
    if dry_run:
        LOGGER.info(
            "Dry run: would archive under %s and reset status.json", archive_dir
        )
        return archive_dir

    clutter_path.mkdir(parents=True, exist_ok=True)
    archive_dir.mkdir()
    status_before = archive_dir / "status_before.json"
    shutil.copy2(run_dir / "status.json", status_before)
    records = []
    for source in stages:
        before = _file_inventory(source)
        destination = archive_dir / source.name
        shutil.move(str(source), destination)
        after = _file_inventory(destination)
        if before != after:
            raise RuntimeError(f"Archive verification failed for {source}")
        records.append(
            {
                "source": str(source),
                "destination": str(destination),
                "sha256_before": before,
                "sha256_after": after,
            }
        )

    status["completed_windows"] = STAGE1
    status["state"] = "stopped"
    status["reason"] = PHASE2_CORRECTION_REASON
    status["updated_utc"] = datetime.now(timezone.utc).isoformat()
    status_path = run_dir / "status.json"
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    manifest = {
        "run_dir": str(run_dir),
        "status_before": str(status_before),
        "status_before_sha256": sha256(status_before),
        "status_after_sha256": sha256(status_path),
        "directories": records,
    }
    (archive_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info("Archived phase-II directories to %s", archive_dir)
    return archive_dir


def main() -> None:
    """Parse arguments and run the archive/reset operation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--clutter-path", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    archive_phase2(args.run_dir, args.clutter_path, args.dry_run)


if __name__ == "__main__":
    main()
