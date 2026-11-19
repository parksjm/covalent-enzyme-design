"""Inventory existing student handoff data; never run or alter scientific inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path("/data/jpc/covalent-enzyme-design")
FF = Path("/data/jpc/polyamide-amber-ff")
LEGACY = Path("/data/jpc/enzyme-substrate-modeling")
OUT = ROOT / "student_handoff/nyl12_ash274_md"
BASE = ROOT / "qmmm/pa6_flipped_model30"
RELEASE = FF / "parameters/releases/pa6_l4_modular_resp_gaff2_v1_20261031"
ENTRIES: dict[str, dict] = {}
MISSING_OPTIONAL: list[str] = []


def sha(path: Path) -> str:
    """Return a streaming source file SHA256."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write(path: Path, value: object) -> None:
    """Create a new finite authoritative JSON inventory artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise FileExistsError(f"Existing generated receipt differs: {path}")
        return
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def add(
    path: Path,
    destination: str,
    category: str,
    *,
    profile: str = "core",
    required: bool = True,
) -> dict:
    """Add one explicit source file; never expand a recursive export."""
    path = path.resolve(strict=True)
    record = dict(
        source=str(path),
        destination=destination,
        profile=profile,
        required=required,
        bytes=path.stat().st_size,
        sha256=sha(path),
        category=category,
    )
    if destination in ENTRIES and ENTRIES[destination] != record:
        raise ValueError(f"Conflicting portable destination: {destination}")
    ENTRIES[destination] = record
    return dict(path=destination.removeprefix("datasets/"), sha256=record["sha256"])
