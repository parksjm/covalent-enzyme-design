"""Hash and copy helpers for frozen relay umbrella preparation bundles."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
PARAMETER_RELEASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/umbrella_sampling/nterminal_relay_releases_20261123"
)
AMBERHOME = Path("/data/jpc/code/ambertools26/venv")
SANDER = AMBERHOME / "bin/sander"
SK_PARAMETER_DIR = AMBERHOME / "dat/slko/3ob-3-1"


def digest(path: Path) -> str:
    """Return streaming SHA-256 for an immutable input file."""
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def validated_parameters() -> dict[str, str]:
    """Require all installed and released SKF files match the validated 3ob set."""
    release = json.loads((PARAMETER_RELEASE / "prepared.json").read_text())
    expected = release.get("parameter_sha256", {})
    names = {f"{first}-{second}.skf" for first in "HCNO" for second in "HCNO"}
    if set(expected) != names:
        raise ValueError("Validated release must declare all sixteen 3ob SKF hashes")
    for name, sha in expected.items():
        installed = SK_PARAMETER_DIR / name
        released = PARAMETER_RELEASE / "parameters_3ob" / name
        if digest(installed) != sha or digest(released) != sha:
            raise ValueError(f"3ob parameter differs from validated release: {name}")
    return expected


def copy_checked(source: Path, destination: Path) -> str:
    """Copy one immutable input and verify byte identity."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    source_hash = digest(source)
    if digest(destination) != source_hash:
        raise ValueError(f"Copied input differs from source: {destination}")
    return source_hash


def verify_prepared(directory: Path) -> dict[str, Any]:
    """Validate only the bundle's copied inputs, source and evidence snapshots."""
    directory = directory.resolve(strict=True)
    prepared = directory / "prepared.json"
    actual_hash = digest(prepared)
    if (directory / "prepared.sha256").read_text().strip() != actual_hash:
        raise ValueError("Prepared manifest digest mismatch")
    payload = json.loads(prepared.read_text())
    for relative, expected in payload["frozen_files_sha256"].items():
        if digest(directory / relative) != expected:
            raise ValueError(f"Frozen pilot file changed: {relative}")
    if payload["native_calculations_launched"] is not False:
        raise ValueError("Preparation manifest incorrectly claims a launch")
    return payload
