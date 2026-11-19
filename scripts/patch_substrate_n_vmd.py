#!/usr/bin/env python3
"""Retain the substrate-N distance label in immutable-publisher VMD exports."""

from __future__ import annotations

import argparse
import hashlib
import logging
from pathlib import Path
import shutil
import time

from scripts.export_direct_pt_progress import PA6_SCISSILE_N, THR233_HG1
from scripts.record_agent_action import append_action

LOGGER = logging.getLogger(__name__)
LABEL = f"label add Bonds $molid/{THR233_HG1} $molid/{PA6_SCISSILE_N}"
ANCHOR = "label textsize 1.0"


def patch_bundle(bundle: Path, backup: Path, journal: Path) -> bool:
    """Add the missing label, preserving originals and the checksum manifest.

    Args:
        bundle: Published bundle containing VMD and its checksum manifest.
        backup: Unused directory for the bundle's original presentation files.
        journal: Append-only action journal.

    Returns:
        Whether a new label was installed.
    """
    vmd = bundle / "reaction_path.vmd"
    manifest = bundle / "MANIFEST.sha256"
    text = vmd.read_text()
    if LABEL in text:
        return False
    if text.count(ANCHOR) != 1:
        raise ValueError(f"Expected one label configuration anchor in {vmd}")
    old_digest = hashlib.sha256(vmd.read_bytes()).hexdigest()
    old_entry = f"{old_digest}  {vmd.name}\n"
    checksums = manifest.read_text()
    if checksums.count(old_entry) != 1:
        raise ValueError(f"VMD checksum does not match published manifest: {bundle}")
    arguments = dict(
        log_path=journal,
        agent="vmd_substrate_n",
        action=f"Add HG1-substrate N VMD distance label: {bundle}",
        permission_mode="escalated",
        reason="Correct inspection display without changing scientific inputs",
        rollback=f"Restore VMD and MANIFEST.sha256 from {backup}",
        files=[str(vmd), str(manifest)],
    )
    append_action(phase="intent", **arguments)
    backup.mkdir(parents=True, exist_ok=False)
    for source in (vmd, manifest):
        shutil.copy2(source, backup / source.name)
    updated = text.replace(ANCHOR, f"{LABEL}\n{ANCHOR}")
    temporary = bundle / ".substrate_n_label.tmp"
    temporary.write_text(updated)
    temporary.replace(vmd)
    new_digest = hashlib.sha256(vmd.read_bytes()).hexdigest()
    temporary.write_text(checksums.replace(old_entry, f"{new_digest}  {vmd.name}\n"))
    temporary.replace(manifest)
    append_action(phase="outcome", result="Label and checksum installed", **arguments)
    LOGGER.info("Added HG1-substrate N label to %s", bundle)
    return True


def patch_once(roots: list[Path], backup: Path, journal: Path) -> int:
    """Patch all fully published generations in the requested progress roots."""
    count = 0
    for root in roots:
        for bundle in sorted(root.glob("update_*")):
            if (bundle / "MANIFEST.sha256").is_file():
                # Keep each pathway's same-numbered generations distinct.
                destination = backup / root.parent.name / bundle.name
                count += patch_bundle(bundle, destination, journal)
    return count


def main() -> None:
    """Patch once, or monitor publication for an explicit bounded duration."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, action="append", required=True)
    parser.add_argument("--backup", type=Path, required=True)
    parser.add_argument("--journal", type=Path, required=True)
    parser.add_argument("--watch-seconds", type=float, default=0)
    parser.add_argument("--interval", type=float, default=15)
    args = parser.parse_args()
    if args.watch_seconds < 0 or args.interval <= 0:
        parser.error("watch-seconds must be nonnegative and interval positive")
    if not all(root.is_dir() for root in args.root):
        parser.error("Every root must be an existing VMD progress directory")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    deadline = time.monotonic() + args.watch_seconds
    while True:
        patch_once(args.root, args.backup, args.journal)
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return
        time.sleep(min(args.interval, remaining))


if __name__ == "__main__":
    main()
