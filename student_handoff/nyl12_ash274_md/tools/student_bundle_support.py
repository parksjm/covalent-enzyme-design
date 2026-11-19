"""Small standard-library path and checksum primitives for student bundles."""

from __future__ import annotations

import hashlib
from pathlib import Path, PurePosixPath
import re

SHA256 = re.compile(r"[0-9a-f]{64}\Z")
RESERVED = {"README.md", "CONTENTS.sha256", "export_manifest.json"}
SOFTWARE_SUFFIXES = {".so", ".dll", ".dylib", ".exe", ".o", ".a", ".pyc"}
SOFTWARE_MAGIC = {
    b"\x7fELF",
    b"\xfe\xed\xfa\xce",
    b"\xce\xfa\xed\xfe",
    b"\xfe\xed\xfa\xcf",
    b"\xcf\xfa\xed\xfe",
}


def safe_relative(value: str) -> PurePosixPath:
    """Return a safe nonempty portable relative path or fail explicitly.

    Args:
        value: POSIX path from a manifest or archive.

    Returns:
        Validated relative path without parent traversal or control characters.
    """
    if not isinstance(value, str) or not value or any(ord(c) < 32 for c in value):
        raise ValueError("Destination must be a nonempty path without controls")
    if (
        "\\" in value
        or ":" in value
        or any(part in {"", ".", ".."} for part in value.split("/"))
    ):
        raise ValueError(f"Unsafe relative path: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute():
        raise ValueError(f"Absolute destination is forbidden: {value}")
    return path


def digest(path: Path) -> str:
    """Return streaming SHA-256 without loading a scientific file into memory."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def regular_source(path: Path) -> None:
    """Require an absolute regular source path without symlink components."""
    if not path.is_absolute():
        raise ValueError(f"Source must be absolute: {path}")
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise ValueError(f"Symlink source is forbidden: {path}")
    if not path.is_file():
        raise FileNotFoundError(path)


def verify_contents(directory: Path) -> None:
    """Verify every listed payload and reject unlisted files or symlinks."""
    names = set()
    for line in (directory / "CONTENTS.sha256").read_text().splitlines():
        checksum, name = line.split("  ", 1)
        relative = safe_relative(name)
        if not SHA256.fullmatch(checksum) or name in names:
            raise ValueError("Invalid or duplicate CONTENTS record")
        path = directory.joinpath(*relative.parts)
        regular_source(path)
        if digest(path) != checksum:
            raise ValueError(f"Payload checksum mismatch: {name}")
        names.add(name)
    actual = {
        str(path.relative_to(directory))
        for path in directory.rglob("*")
        if path.is_file()
    }
    if actual != names | {"CONTENTS.sha256"}:
        raise ValueError("Bundle contains missing or unlisted payload files")
