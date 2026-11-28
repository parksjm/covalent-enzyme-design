"""Preserve originals and refresh derived portable-viewer checksums."""

import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Any


def digest(path: Path) -> str:
    """Return a streaming SHA-256 digest of a file."""
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def backup_file(path: Path, root: Path, backup: Path) -> dict[str, str]:
    """Preserve an original derived file once, refusing conflicting backups.

    Args:
        path: Existing file inside root.
        root: Repository root.
        backup: Central backup directory retaining repository-relative paths.

    Returns:
        Original path, backup path and original SHA-256.

    Raises:
        ValueError: If an existing backup contains different bytes.
    """
    destination = backup / path.relative_to(root)
    old_hash = digest(path)
    if destination.exists():
        if digest(destination) != old_hash:
            raise ValueError(f"Conflicting preserved original: {destination}")
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    return dict(path=str(path), backup=str(destination), before_sha256=old_hash)


def checksum_updates(folder: Path, changes: dict[str, str]) -> dict[Path, str]:
    """Render updates to portable checksums and derived viewer hash maps.

    Args:
        folder: Movie bundle directory.
        changes: Changed filenames mapped to their new content, in memory.

    Returns:
        Derived JSON/checksum files mapped to replacement text. Immutable source
        identity and scientific acceptance records are never rewritten.
    """
    result = {}
    hashes = {
        name: hashlib.sha256(text.encode()).hexdigest()
        for name, text in changes.items()
    }
    for name in ("manifest.json", "export_manifest.json", "source_manifest.json"):
        path = folder / name
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        changed = False
        # These explicit keys describe current exported files. Source ancestry
        # keys (source_sha256, sources, source_identity, snapshots) stay intact.
        for key in ("files", "bundle_files", "export_files"):
            entries = data.get(key)
            if not isinstance(entries, dict):
                continue
            for filename, new_hash in hashes.items():
                if filename in entries:
                    value = entries[filename]
                    if isinstance(value, str):
                        entries[filename] = new_hash
                    elif isinstance(value, dict) and "sha256" in value:
                        value["sha256"] = new_hash
                    else:
                        raise ValueError(
                            f"Unsupported derived hash entry {path}:{filename}"
                        )
                    changed = True
        if changed:
            text = json.dumps(data, indent=2) + "\n"
            result[path] = text
            hashes[name] = hashlib.sha256(text.encode()).hexdigest()
    for name in ("MANIFEST.sha256", "SHA256SUMS.txt", "SHA256SUMS"):
        path = folder / name
        if not path.exists():
            continue
        lines = []
        for line in path.read_text().splitlines():
            match = re.fullmatch(r"([0-9a-f]{64})(\s+)(\*?)(.+)", line)
            if not match:
                raise ValueError(f"Unsupported checksum record in {path}: {line!r}")
            old, space, star, filename = match.groups()
            normalized = filename.removeprefix("./")
            lines.append(f"{hashes.get(normalized, old)}{space}{star}{filename}")
        result[path] = "\n".join(lines) + "\n"
    return result


def verify_checksums(folder: Path) -> list[dict[str, Any]]:
    """Return mismatching current portable checksum entries without rewriting.

    Args:
        folder: Movie bundle containing portable SHA-256 manifests.

    Returns:
        Mismatch records; an empty list establishes complete listed integrity.
    """
    failures = []
    for name in ("MANIFEST.sha256", "SHA256SUMS.txt", "SHA256SUMS"):
        manifest = folder / name
        if not manifest.exists():
            continue
        for line in manifest.read_text().splitlines():
            expected, filename = line.split(maxsplit=1)
            target = folder / filename.lstrip("*")
            if not target.is_file() or digest(target) != expected:
                failures.append(
                    dict(manifest=str(manifest), file=str(target), expected=expected)
                )
    return failures
