"""Build only explicitly inventoried student data; never package raw references."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import shutil
import tarfile
from typing import Any

from student_bundle_support import (
    RESERVED,
    SHA256,
    SOFTWARE_MAGIC,
    SOFTWARE_SUFFIXES,
    digest,
    regular_source,
    safe_relative,
    verify_contents,
)

LOGGER = logging.getLogger(__name__)


def selected_entries(
    manifest: dict[str, Any], profile: str
) -> tuple[list[dict[str, Any]], list[str]]:
    """Validate the explicit manifest and select available profile entries.

    Args:
        manifest: Object with entries containing the agreed source/destination contract.
        profile: Core alone, or inspection containing both core and inspection.

    Returns:
        Available selected entries and optional selected sources that are absent.
    """
    selected, skipped, destinations = [], [], set()
    for entry in manifest["entries"]:
        if entry["profile"] not in {"core", "inspection", "raw_reference_only"}:
            raise ValueError("Unknown manifest profile")
        if entry["profile"] == "raw_reference_only" or (
            profile == "core" and entry["profile"] == "inspection"
        ):
            continue
        name = str(safe_relative(entry["destination"]))
        if name in RESERVED or name in destinations:
            raise ValueError(f"Reserved or duplicate destination: {name}")
        if any(
            name.startswith(previous + "/") or previous.startswith(name + "/")
            for previous in destinations
        ):
            raise ValueError(f"Conflicting file/directory destinations: {name}")
        destinations.add(name)
        if (
            type(entry["required"]) is not bool
            or type(entry["bytes"]) is not int
            or entry["bytes"] < 0
        ):
            raise ValueError(
                "required must be bool and bytes must be a nonnegative int"
            )
        if not isinstance(entry["sha256"], str) or not SHA256.fullmatch(
            entry["sha256"]
        ):
            raise ValueError(f"Selected entry needs SHA-256: {name}")
        if not isinstance(entry["category"], str) or not entry["category"]:
            raise ValueError("Selected entry needs a category")
        source = Path(entry["source"])
        if not source.is_absolute():
            raise ValueError("Selected source must be absolute")
        if not source.exists() and not source.is_symlink() and not entry["required"]:
            skipped.append(str(source))
            continue
        regular_source(source)
        if source.suffix.lower() in SOFTWARE_SUFFIXES:
            raise ValueError(f"Software binaries are forbidden: {source}")
        selected.append(entry)
    if not selected:
        raise ValueError("No available explicit entries selected")
    return selected, skipped


def copy_verified(entry: dict[str, Any], destination: Path) -> None:
    """Stream one source into a new file and verify its exact size and SHA-256."""
    source = Path(entry["source"])
    regular_source(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    checksum, count = hashlib.sha256(), 0
    with source.open("rb") as incoming, destination.open("xb") as outgoing:
        prefix = incoming.read(4)
        if prefix in SOFTWARE_MAGIC or prefix.startswith(b"MZ"):
            raise ValueError(f"Software executable magic is forbidden: {source}")
        for block in (prefix,):
            outgoing.write(block)
            checksum.update(block)
            count += len(block)
        for block in iter(lambda: incoming.read(1024 * 1024), b""):
            outgoing.write(block)
            checksum.update(block)
            count += len(block)
    if count != entry["bytes"] or checksum.hexdigest() != entry["sha256"]:
        raise ValueError(f"Source size/SHA mismatch: {source}")
    shutil.copystat(source, destination, follow_symlinks=False)


def build(manifest_path: Path, profile: str, output: Path) -> Path:
    """Build a new directory/archive from only the selected manifest profile.

    Args:
        manifest_path: JSON object with explicit entries; no directory discovery.
        profile: core or inspection (the latter includes core).
        output: New nonexistent directory; siblings archive/sidecar must also be absent.

    Returns:
        Verified gzip tar archive, with its own adjacent checksum sidecar.
    """
    if profile not in {"core", "inspection"}:
        raise ValueError("Profile must be core or inspection")
    output = output.absolute()
    archive = output.with_name(output.name + ".tar.gz")
    sidecar = archive.with_name(archive.name + ".sha256")
    if any(path.exists() or path.is_symlink() for path in (output, archive, sidecar)):
        raise FileExistsError("Output directory/archive/checksum already exists")
    manifest = json.loads(manifest_path.read_text())
    entries, skipped = selected_entries(manifest, profile)
    output.mkdir(parents=True)
    for entry in entries:
        relative = safe_relative(entry["destination"])
        copy_verified(entry, output.joinpath(*relative.parts))
    record = {
        "profile": profile,
        "source_manifest_sha256": digest(manifest_path),
        "entries": entries,
        "skipped_optional_sources": skipped,
        "raw_references_included": False,
        "software_binaries_included": False,
    }
    (output / "export_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    index = (
        "# Student handoff bundle\n\n"
        + f"Profile: **{profile}**. Verify CONTENTS.sha256 before use.\n\n"
    )
    index += "This bundle contains only explicitly inventoried files; raw MD references and software binaries are excluded.\n\n"
    index += (
        "\n".join(
            f"- [{entry['destination']}]({entry['destination']})" for entry in entries
        )
        + "\n"
    )
    (output / "README.md").write_text(index)
    files = sorted(path for path in output.rglob("*") if path.is_file())
    (output / "CONTENTS.sha256").write_text(
        "".join(f"{digest(path)}  {path.relative_to(output)}\n" for path in files)
    )
    verify_contents(output)
    with tarfile.open(archive, "x:gz", dereference=False) as bundle:
        bundle.add(output, arcname=output.name, recursive=True)
    sidecar.write_text(f"{digest(archive)}  {archive.name}\n")
    LOGGER.info("Built %s: %d entries, profile=%s", archive, len(entries), profile)
    return archive


def main(argv: list[str] | None = None) -> None:
    """Parse explicit manifest/profile/output arguments and build a new bundle."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--profile", choices=("core", "inspection"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    build(args.manifest, args.profile, args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
