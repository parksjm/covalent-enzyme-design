"""Content-addressed storage for reproducible Nyl12 parameter artifacts."""

from __future__ import annotations

import fcntl
import hashlib
import json
import logging
import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from typing import Any, Iterator, Mapping

LOGGER = logging.getLogger(__name__)
MANIFEST = "manifest.json"


def _digest(path: Path) -> str:
    """Return the SHA256 digest of a file, streaming large artifacts."""
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def _artifact_name(name: str) -> str:
    """Validate a portable, relative artifact name."""
    candidate = PurePosixPath(name)
    if (
        not name
        or candidate.is_absolute()
        or ".." in candidate.parts
        or name in {".", ".."}
        or "." in candidate.parts
        or name == MANIFEST
        or "\\" in name
    ):
        raise ValueError(f"Unsafe or reserved artifact name: {name!r}")
    return candidate.as_posix()


class ParameterArtifactCache:
    """Store and materialize parameter files using scientific identity hashes.

    Args:
        root: Directory that holds immutable cache entries.

    An identity should include all scientific inputs that can change the
    parameters, such as graph, state, protonation, conformer, charge method,
    tool versions, and force-field fingerprints.
    """

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self._lock_path = self.root / ".cache.lock"

    @staticmethod
    def key_for(identity: Any) -> str:
        """Return a stable SHA256 key for a JSON-serializable identity."""
        encoded = json.dumps(
            identity, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def _locked(self) -> Any:
        """Open an exclusive interprocess lock for cache mutations and reads."""
        lock = self._lock_path.open("a+b")
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        return lock

    @contextmanager
    def calculation_lock(self, identity: Any) -> Iterator[None]:
        """Serialize expensive calculations for one identity across processes."""
        key = self.key_for(identity)
        lock_dir = self.root / ".locks"
        lock_dir.mkdir(exist_ok=True)
        with (lock_dir / f"{key}.lock").open("a+b") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            yield
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)

    @staticmethod
    def _canonical_identity(identity: Any) -> Any:
        """Normalize JSON values so tuples and lists compare consistently."""
        encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"),
                             ensure_ascii=False, allow_nan=False)
        return json.loads(encoded)

    def _read_entry(self, identity: Any) -> dict[str, Any] | None:
        identity = self._canonical_identity(identity)
        key = self.key_for(identity)
        entry = self.root / key
        if not entry.exists():
            return None
        manifest_path = entry / MANIFEST
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(f"Corrupt parameter cache manifest: {manifest_path}") from error
        if manifest.get("key") != key or manifest.get("identity") != identity:
            raise RuntimeError(f"Parameter cache identity mismatch: {entry}")
        files = manifest.get("files")
        if not isinstance(files, dict) or not files:
            raise RuntimeError(f"Invalid parameter cache file list: {manifest_path}")
        for name, expected in files.items():
            safe_name = _artifact_name(name)
            path = entry / safe_name
            if not path.is_file() or _digest(path) != expected:
                raise RuntimeError(f"Parameter cache artifact is corrupt: {path}")
        return manifest

    def store(self, identity: Any, artifacts: Mapping[str, Path]) -> str:
        """Store named existing files and return their content-addressed key.

        Args:
            identity: Complete JSON-serializable scientific provenance.
            artifacts: Mapping from stable artifact names to existing files.

        Raises:
            FileNotFoundError: An input artifact does not exist.
            RuntimeError: An existing cache entry is corrupt.
            ValueError: Identity or artifact names are invalid.
        """
        identity = self._canonical_identity(identity)
        key = self.key_for(identity)
        if not artifacts:
            raise ValueError("At least one parameter artifact is required")
        normalized = {_artifact_name(name): Path(source) for name, source in artifacts.items()}
        if len(normalized) != len(artifacts):
            raise ValueError("Artifact names must be unique after normalization")
        names = set(normalized)
        if any(any(other.startswith(f"{name}/") for other in names) for name in names):
            raise ValueError("Artifact names cannot be both files and directories")
        for source in normalized.values():
            if not source.is_file():
                raise FileNotFoundError(source)
        source_hashes = {name: _digest(source) for name, source in normalized.items()}
        with self._locked() as lock:
            existing = self._read_entry(identity)
            if existing is not None:
                if existing["files"] != source_hashes:
                    raise RuntimeError("Existing parameter cache entry has different artifacts")
                LOGGER.info("Reusing parameter cache entry %s", key)
                return key
            staging = Path(tempfile.mkdtemp(prefix=f".{key}.", dir=self.root))
            try:
                hashes: dict[str, str] = {}
                for name, source in normalized.items():
                    target = staging / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source, target)
                    hashes[name] = _digest(target)
                manifest = {"key": key, "identity": identity, "files": hashes}
                (staging / MANIFEST).write_text(
                    json.dumps(manifest, sort_keys=True, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n",
                    encoding="utf-8",
                )
                os.replace(staging, self.root / key)
            finally:
                if staging.exists():
                    shutil.rmtree(staging)
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        return key

    def materialize(self, identity: Any, target_dir: Path) -> list[Path] | None:
        """Verify and copy a cached entry into an empty target directory.

        Returns ``None`` on a cache miss; corruption and target collisions raise.
        Existing target files are never overwritten.
        """
        target_dir = Path(target_dir)
        with self._locked():
            manifest = self._read_entry(identity)
            if manifest is None:
                return None
            names = [_artifact_name(name) for name in manifest["files"]]
            collisions = [target_dir / name for name in names if (target_dir / name).exists()]
            if collisions:
                raise FileExistsError(f"Refusing to overwrite cached artifacts: {collisions}")
            target_dir.parent.mkdir(parents=True, exist_ok=True)
            staging = Path(tempfile.mkdtemp(prefix=".parameter-cache.", dir=target_dir.parent))
            try:
                entry = self.root / manifest["key"]
                for name in names:
                    staged = staging / name
                    staged.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(entry / name, staged)
                target_dir.mkdir(parents=True, exist_ok=True)
                paths = []
                for name in names:
                    destination = target_dir / name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(staging / name, destination)
                    paths.append(destination)
                return paths
            finally:
                if staging.exists():
                    shutil.rmtree(staging)
