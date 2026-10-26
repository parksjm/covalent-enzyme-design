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
    destination.mkdir(parents=True, exist_ok=True)
    manifest = destination / "MANIFEST.sha256"
    original = manifest.read_text() if manifest.exists() else ""
    (generation / "previous_download_manifest.sha256").write_text(original)
    entries = {}
    for line in original.splitlines():
        digest, name = line.split(maxsplit=1)
        entries[name.strip()] = digest
    files = {
        "reaction_path.dcd": "reaction_path.dcd",
        "reaction_path.vmd": "reaction_path.vmd",
        "stripped.parm7.gz": "stripped.parm7.gz",
        "qm_indices.tcl": "qm_indices.tcl",
        "nyl12_es_vmd_common.tcl": "nyl12_es_vmd_common.tcl",
        "nyl12_reactive_distances.tcl": "nyl12_reactive_distances.tcl",
        "reaction_distances.png": "distances.png",
        "reaction_distances.json": "distances.json",
        "reaction_distances.csv": "distances.csv",
        "reaction_energies.png": "energies.png",
        "reaction_energies.json": "energies.json",
        "reaction_energies.csv": "energies.csv",
    }
    for source_name, destination_name in files.items():
        source = generation / source_name
        pending = destination / ("." + destination_name + ".pending")
        shutil.copy2(source, pending)
        entries[destination_name] = hashlib.sha256(source.read_bytes()).hexdigest()
        if (
            hashlib.sha256(pending.read_bytes()).hexdigest()
            != entries[destination_name]
        ):
            raise ValueError("Download copy digest mismatch: " + destination_name)
        pending.replace(destination / destination_name)
    pending_manifest = destination / ".MANIFEST.pending"
    pending_manifest.write_text(
        "".join(f"{digest}  {name}\n" for name, digest in entries.items())
    )
    pending_manifest.replace(manifest)
