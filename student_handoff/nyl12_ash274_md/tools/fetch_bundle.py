"""Fetch one selected student archive and its own checksum, then optionally extract."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile

from student_bundle_support import SHA256, digest, safe_relative, verify_contents

LOGGER = logging.getLogger(__name__)


def extract_verified(archive: Path, destination: Path, package: str) -> Path:
    """Extract a checked archive without links, devices, traversal or overwrites.

    Args:
        archive: Archive whose transport checksum has already been verified.
        destination: Existing download directory.
        package: Expected single top-level archive directory.

    Returns:
        Extracted directory after verifying all inner CONTENTS checksums.
    """
    root = destination / package
    if root.exists() or root.is_symlink():
        raise FileExistsError(root)
    with tarfile.open(archive, "r:gz") as bundle:
        members = bundle.getmembers()
        names = set()
        for member in members:
            relative = safe_relative(member.name.rstrip("/"))
            if relative.parts[0] != package or not (member.isfile() or member.isdir()):
                raise ValueError(
                    "Archive must contain only regular files/directories under selected package"
                )
            if str(relative) in names:
                raise ValueError("Duplicate archive member")
            names.add(str(relative))
        for member in members:
            target = destination.joinpath(*safe_relative(member.name.rstrip("/")).parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                stream = bundle.extractfile(member)
                if stream is None:
                    raise ValueError("Regular archive member has no payload")
                with stream, target.open("xb") as output:
                    shutil.copyfileobj(stream, output, 1024 * 1024)
    verify_contents(root)
    return root


def fetch(
    source: str, remote_dir: str, package: str, destination: Path, extract: bool
) -> Path:
    """Download exactly one archive and its own checksum without remote refreshes.

    Args:
        source: Explicit safe SSH alias or user@host.
        remote_dir: Absolute POSIX directory containing published archives.
        package: core_data or inspection, matching the builder output basename.
        destination: Local download directory; existing archive files are rejected.
        extract: Whether to verify/extract the selected archive's payload too.

    Returns:
        Verified archive path, or verified extracted package directory.
    """
    if not re.fullmatch(
        r"(?:[A-Za-z0-9_][A-Za-z0-9_.-]*@)?[A-Za-z0-9_][A-Za-z0-9_.-]*", source
    ):
        raise ValueError("Source must be a safe SSH alias or user@host")
    remote = PurePosixPath(remote_dir)
    if (
        not remote.is_absolute()
        or not re.fullmatch(r"/[A-Za-z0-9_./-]+", remote_dir)
        or ".." in remote.parts
    ):
        raise ValueError("Remote directory must be a safe absolute POSIX path")
    if package not in {"core_data", "inspection"}:
        raise ValueError("Select core_data or inspection")
    destination = destination.absolute()
    if destination.is_symlink():
        raise ValueError("Symlink download directory is forbidden")
    archive = destination / f"{package}.tar.gz"
    sidecar = destination / f"{package}.tar.gz.sha256"
    if any(path.exists() or path.is_symlink() for path in (archive, sidecar)):
        raise FileExistsError("Selected archive/checksum already exists locally")
    destination.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "rsync",
            "-rt",
            "--protect-args",
            "--",
            f"{source}:{remote}/{archive.name}",
            f"{source}:{remote}/{sidecar.name}",
            str(destination) + "/",
        ],
        check=True,
    )
    lines = sidecar.read_text().splitlines()
    if len(lines) != 1 or "  " not in lines[0]:
        raise ValueError("Expected one selected-archive checksum record")
    checksum, name = lines[0].split("  ", 1)
    if (
        name != archive.name
        or not SHA256.fullmatch(checksum)
        or digest(archive) != checksum
    ):
        raise ValueError("Selected archive checksum/name mismatch")
    LOGGER.info("Verified selected archive %s", archive)
    return extract_verified(archive, destination, package) if extract else archive


def main(argv: list[str] | None = None) -> None:
    """Fetch only the chosen package; optional packages have independent checksums."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--remote-dir", required=True)
    parser.add_argument("--package", choices=("core_data", "inspection"), required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--extract", action="store_true")
    args = parser.parse_args(argv)
    fetch(args.source, args.remote_dir, args.package, args.destination, args.extract)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
