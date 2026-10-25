"""Publish verified small path artifacts into the existing laptop package."""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path


def publish_download(generation: Path, destination: Path) -> None:
    """Replace live inspection files atomically, updating checksums last.

    Args:
        generation: Immutable, verified inspection generation.
        destination: Existing curated download directory with checksum manifest.
    """
    manifest = destination / "MANIFEST.sha256"
    original = manifest.read_text()
    (generation / "previous_download_manifest.sha256").write_text(original)
    entries = {}
    for line in original.splitlines():
        digest, name = line.split(maxsplit=1)
        entries[name.strip()] = digest
    names = (
        "reaction_path.dcd",
        "reaction_path.vmd",
        "reaction_distances.png",
        "reaction_distances.json",
        "reaction_distances.csv",
    )
    for name in names:
        source = generation / name
        pending = destination / ("." + name + ".pending")
        shutil.copy2(source, pending)
        entries[name] = hashlib.sha256(source.read_bytes()).hexdigest()
        if hashlib.sha256(pending.read_bytes()).hexdigest() != entries[name]:
            raise ValueError("Download copy digest mismatch: " + name)
        pending.replace(destination / name)
    pending_manifest = destination / ".MANIFEST.pending"
    pending_manifest.write_text(
        "".join(f"{digest}  {name}\n" for name, digest in entries.items())
    )
    pending_manifest.replace(manifest)
