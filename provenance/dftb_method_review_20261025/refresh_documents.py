"""Refresh the curated document export while preserving its previous revision."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.prepare_nyl12_inspection_download import build_package  # noqa: E402

DOCUMENTS = (
    "manuscript/INTRODUCTION_DRAFT.md",
    "manuscript/METHODS_DRAFT.md",
    "manuscript/REFERENCES.md",
    "manuscript/recent_design_and_relay_references.bib",
    "docs/asp272_relay_pilot_20261025.md",
    "docs/design_ranking_and_iteration_20261022.md",
    "docs/standalone_handoff_20261020.md",
)


def digest(path: Path) -> str:
    """Return the SHA-256 digest of one retained small file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Verify a fresh package before replacing the stable download location."""
    archive = Path(__file__).resolve().parent
    current = ROOT / "output_laptop_inspection_20261021"
    staging = ROOT / "output_laptop_inspection_dftb_review_staging"
    backup = ROOT / "output_laptop_inspection_before_dftb_review_20261025"
    for path in (staging, backup, archive / "snapshot"):
        if path.exists():
            raise FileExistsError(path)
    prior = json.loads((current / "provenance.json").read_text())
    files: dict[str, Path] = {}
    for name, row in prior["files"].items():
        source = current / name
        if digest(source) != row["sha256"]:
            raise ValueError(f"Previous package hash mismatch: {name}")
        files[name] = source
    for name in DOCUMENTS:
        target = archive / "snapshot" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
        if digest(target) != digest(ROOT / name):
            raise ValueError(f"Document changed during snapshot: {name}")
        files[name] = target
    # Retain existing provenance origins for all unchanged payload files.
    for name in set(files) - set(DOCUMENTS):
        files[name] = Path(prior["files"][name]["source"])
        if digest(files[name]) != prior["files"][name]["sha256"]:
            raise ValueError(f"Original source changed: {name}")
    report = build_package(files, staging)
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if digest(staging / name) != expected:
            raise ValueError(f"Fresh package hash mismatch: {name}")
    current.rename(backup)
    staging.rename(current)
    result = {
        "verified": True,
        "documents": list(DOCUMENTS),
        "retained_backup": str(backup),
        "files": len(report["files"]),
        "payload_bytes": report["payload_bytes"],
    }
    (archive / "export_verification.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    logging.info("Verified refreshed document export: %s files", result["files"])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
