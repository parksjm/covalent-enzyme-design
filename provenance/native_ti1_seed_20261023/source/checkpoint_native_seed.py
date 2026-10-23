"""Preserve native seed evidence and refresh the reversible laptop export."""

from __future__ import annotations

import fcntl
import hashlib
import json
import logging
import shutil
import sys
from pathlib import Path

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(ROOT))
from scripts.prepare_nyl12_inspection_download import build_package  # noqa: E402

BASE = ROOT / "output_nyl12_native_reference_20261023"
ARCHIVE = ROOT / "provenance/native_ti1_seed_20261023"


def digest(path: Path) -> str:
    """Return the complete file SHA256."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Archive small immutable evidence and verify a new inspection package."""
    current = ROOT / "output_laptop_inspection_20261021"
    backup = ROOT / "output_laptop_inspection_20261021_revision_17"
    staging = ROOT / "output_laptop_inspection_20261021_native_staging"
    if any(path.exists() for path in (ARCHIVE, backup, staging)):
        raise FileExistsError("Checkpoint requires fresh archive/staging/backup")
    seed = BASE / "native_seed_01"
    if not (seed / "preflight.json").is_file():
        raise FileNotFoundError("Retain actual seed assessment before checkpoint")
    ARCHIVE.mkdir()
    for name in ("source_audit_01", "native_seed_01", "freeze_flags_preflight_01"):
        for source in (BASE / name).rglob("*"):
            if (
                source.is_file()
                and source.stat().st_size < 1_000_000
                and source.suffix not in {".chk", ".fchk", ".pyc"}
            ):
                target = ARCHIVE / name / source.relative_to(BASE / name)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                if digest(source) != digest(target):
                    raise ValueError(f"Copy-time hash mismatch: {source}")
    pilot = BASE / "native_frozen_scaffold_pilot_01"
    if pilot.is_dir():
        for name in (
            "input.gjf",
            "run_gaussian.sh",
            "pilot_policy.json",
            "staging_manifest.json",
        ):
            source = pilot / name
            target = ARCHIVE / "pilot_inputs" / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    for source in [*BASE.glob("*.py"), BASE / "PILOT_PROTOCOL.md"]:
        target = ARCHIVE / "source" / source.name
        target.parent.mkdir(exist_ok=True)
        shutil.copy2(source, target)
    for source in (ROOT / "docs").glob("*.md"):
        if source.name in {
            "native_ti1_reference_20261023.md",
            "standalone_handoff_20261020.md",
            "laptop_inspection_20261021.md",
        }:
            target = ARCHIVE / "docs" / source.name
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(source, target)
    journal = ARCHIVE / "agent_actions_at_checkpoint.jsonl"
    with (ROOT / "output_agent_actions/20261021.jsonl").open("rb") as source:
        fcntl.flock(source, fcntl.LOCK_SH)
        with journal.open("xb") as target:
            shutil.copyfileobj(source, target)
        fcntl.flock(source, fcntl.LOCK_UN)
    old = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in old["files"].items()}
    changed = {"README.md", "provenance/agent_actions_at_export.jsonl"}
    files["README.md"] = ROOT / "docs/laptop_inspection_20261021.md"
    files["provenance/agent_actions_at_export.jsonl"] = journal
    for source in (ROOT / "docs").glob("*.md"):
        name = f"docs/{source.name}"
        changed.add(name)
        files[name] = source
    for source in seed.iterdir():
        if source.is_file() and source.suffix in {".sdf", ".xyz", ".json", ".pml"}:
            files[f"qm/native_ti1_seed/{source.name}"] = source
    result = build_package(files, staging)
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if digest(staging / name) != expected:
            raise ValueError(f"Package hash mismatch: {name}")
    for name in old["files"]:
        if (
            name not in changed
            and (current / name).read_bytes() != (staging / name).read_bytes()
        ):
            raise ValueError(f"Old nondocument payload changed: {name}")
    current.rename(backup)
    staging.rename(current)
    report = {
        "verified": True,
        "payload_files": len(result["files"]),
        "payload_bytes": result["payload_bytes"],
        "retained_backup": str(backup),
        "forcefield_validated": False,
        "old_nondocument_payloads_preserved": True,
    }
    (ARCHIVE / "package_verification.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    hashes = {
        str(p.relative_to(ARCHIVE)): digest(p)
        for p in ARCHIVE.rglob("*")
        if p.is_file()
    }
    (ARCHIVE / "source_manifest.json").write_text(json.dumps(hashes, indent=2) + "\n")
    logging.info("Native checkpoint: %s", report)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
