#!/usr/bin/env python3
"""Move superseded export and test-output directories into reversible clutter."""

from __future__ import annotations

import hashlib
import json
import logging
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[2]
DESTINATION_ROOT = ROOT / "clutter" / "superseded_exports_and_test_runs_20261025"
MANIFEST_PATH = Path(__file__).with_name("move_manifest.json")

EXACT_NAMES = (
    "output_laptop_inspection_before_dftb_review_20261025",
    "output_path_full_pytest_20261025",
    "output_path_inspection_allcov_20261025",
    "output_path_inspection_cov_retry_20261025",
    "output_path_inspection_pytest_20261025",
    "output_path_inspection_pytest_retry2_20261025",
    "output_path_inspection_pytest_retry_20261025",
    "output_qmmm_final_full_pytest_20261025",
    "output_qmmm_package_retry_20261025",
    "output_qmmm_package_tests_20261025",
    "output_vmd_inspection_20261024_failed_attempt2",
    "output_vmd_inspection_20261024_failed_attempt3",
    "output_vmd_inspection_20261024_failed_attempt4",
)
REVISION_PATTERN = re.compile(r"output_laptop_inspection_20261021_revision_\d{2}")
PROTECTED_NAMES = {
    "output_agent_actions",
    "output_laptop_inspection_20261021",
    "output_parameter_cache",
    "output_vmd_inspection_20261024",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _directory_inventory(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for candidate in sorted(path.rglob("*")):
        relative = candidate.relative_to(path).as_posix()
        if candidate.is_symlink():
            records.append(
                {
                    "kind": "symlink",
                    "path": relative,
                    "target": candidate.readlink().as_posix(),
                }
            )
        elif candidate.is_file():
            records.append(
                {
                    "kind": "file",
                    "path": relative,
                    "sha256": _sha256(candidate),
                    "size_bytes": candidate.stat().st_size,
                }
            )
    return records


def _write_manifest(manifest: dict[str, Any]) -> None:
    temporary = MANIFEST_PATH.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(MANIFEST_PATH)


def _selected_directories() -> list[Path]:
    revisions = sorted(
        path
        for path in ROOT.iterdir()
        if path.is_dir() and REVISION_PATTERN.fullmatch(path.name)
    )
    if len(revisions) != 37:
        raise RuntimeError(f"Expected 37 laptop revisions, found {len(revisions)}")
    exact = [ROOT / name for name in EXACT_NAMES]
    missing = [path.as_posix() for path in exact if not path.is_dir()]
    if missing:
        raise FileNotFoundError(f"Expected cleanup directories are missing: {missing}")
    selected = revisions + exact
    protected = sorted(path.name for path in selected if path.name in PROTECTED_NAMES)
    if protected:
        raise RuntimeError(f"Protected directories selected: {protected}")
    return selected


def move_superseded_outputs() -> dict[str, Any]:
    """Move the exact audited cleanup set and verify content hashes.

    Returns:
        Completed cleanup manifest.

    Raises:
        FileExistsError: If the holding directory already exists.
        FileNotFoundError: If an audited source directory is missing.
        RuntimeError: If the expected set or post-move hashes differ.
    """
    selected = _selected_directories()
    if DESTINATION_ROOT.exists():
        raise FileExistsError(f"Destination already exists: {DESTINATION_ROOT}")
    manifest: dict[str, Any] = {
        "created_utc": datetime.now(UTC).isoformat(),
        "destination_root": DESTINATION_ROOT.relative_to(ROOT).as_posix(),
        "entries": [],
        "policy": "reversible whole-directory moves; no deletion",
        "repository_root": ROOT.as_posix(),
        "status": "in_progress",
    }
    DESTINATION_ROOT.mkdir(parents=True)
    _write_manifest(manifest)
    for source in selected:
        inventory = _directory_inventory(source)
        destination = DESTINATION_ROOT / source.name
        entry = {
            "destination": destination.relative_to(ROOT).as_posix(),
            "files": inventory,
            "source": source.relative_to(ROOT).as_posix(),
            "status": "inventoried",
            "total_bytes": sum(item.get("size_bytes", 0) for item in inventory),
        }
        manifest["entries"].append(entry)
        _write_manifest(manifest)
        shutil.move(source.as_posix(), destination.as_posix())
        entry["status"] = "moved"
        _write_manifest(manifest)
        if _directory_inventory(destination) != inventory:
            raise RuntimeError(f"Post-move hash verification failed: {source.name}")
        entry["status"] = "verified"
        _write_manifest(manifest)
        LOGGER.info("Moved and verified %s", source.name)
    manifest["completed_utc"] = datetime.now(UTC).isoformat()
    manifest["directory_count"] = len(manifest["entries"])
    manifest["status"] = "complete"
    manifest["total_bytes"] = sum(entry["total_bytes"] for entry in manifest["entries"])
    manifest["rollback"] = (
        "Move each destination directory back to repository_root/source; "
        "verify against the recorded per-file hashes."
    )
    _write_manifest(manifest)
    return manifest


def main() -> int:
    """Run the reversible cleanup.

    Returns:
        Process exit status.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    manifest = move_superseded_outputs()
    LOGGER.info(
        "Moved %d directories (%d bytes) into %s",
        manifest["directory_count"],
        manifest["total_bytes"],
        manifest["destination_root"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
