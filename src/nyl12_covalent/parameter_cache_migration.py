"""Explicit migrations for legacy parameter-cache directory names."""

from __future__ import annotations

import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .parameter_cache import _entry_label

if TYPE_CHECKING:
    from .parameter_cache import ParameterArtifactCache


def migrate_entry(cache: ParameterArtifactCache, identity: Any, label: str) -> Path:
    """Rename and record a verified hash-named cache entry.

    Args:
        cache: Cache containing the legacy entry.
        identity: Complete JSON-serializable identity used for lookup.
        label: Safe readable directory name for the migrated entry.

    Returns:
        New cache-entry directory.

    Raises:
        FileNotFoundError: No entry matches the identity.
        FileExistsError: The destination label or record already exists.
        RuntimeError: The entry manifest or artifacts are corrupt.
        ValueError: The entry is already named or the label is unsafe.
    """
    identity = cache._canonical_identity(identity)
    key = cache.key_for(identity)
    directory_name = _entry_label(label)
    records = cache.root / ".migrations"
    record_path = records / f"{key}.json"
    with cache._locked():
        located = cache._find_entry(identity)
        if located is None:
            raise FileNotFoundError(f"No parameter cache entry for {key}")
        source, _ = located
        if source.name != key:
            raise ValueError(f"Only hash-named entries can be migrated: {source}")
        destination = cache.root / directory_name
        if destination.exists():
            raise FileExistsError(
                f"Parameter cache label is already in use: {destination}"
            )
        if record_path.exists():
            raise FileExistsError(f"Migration record already exists: {record_path}")
        records.mkdir(exist_ok=True)
        record = {
            "key": key,
            "old_directory": key,
            "new_directory": directory_name,
            "migrated_utc": datetime.now(timezone.utc).isoformat(),
        }
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=records, delete=False
        ) as stream:
            json.dump(record, stream, indent=2, sort_keys=True)
            stream.write("\n")
            temporary_record = Path(stream.name)
        try:
            os.replace(source, destination)
            os.replace(temporary_record, record_path)
        except Exception:
            if destination.exists() and not source.exists():
                os.replace(destination, source)
            raise
        finally:
            temporary_record.unlink(missing_ok=True)
        return destination
