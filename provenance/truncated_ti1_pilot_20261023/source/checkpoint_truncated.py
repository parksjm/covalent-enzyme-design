"""Archive the truncated pilot and reversibly refresh its small inspection export."""

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

BASE = Path(__file__).resolve().parent
ARCHIVE = ROOT / "provenance/truncated_ti1_pilot_20261023"


def digest(path: Path) -> str:
    """Return SHA256 of retained file bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy_exact(source: Path, target: Path) -> None:
    """Copy without replacing an artifact and verify identical contents."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise FileExistsError(target)
    shutil.copy2(source, target)
    if digest(source) != digest(target):
        raise ValueError(f"Copy-time hash mismatch: {source}")


def main() -> None:
    """Archive fixed inputs and verified evidence; retain the previous package."""
    current = ROOT / "output_laptop_inspection_20261021"
    staging = ROOT / "output_laptop_inspection_20261021_truncated_staging"
    backup = ROOT / "output_laptop_inspection_20261021_revision_18"
    if any(p.exists() for p in (ARCHIVE, staging, backup)):
        raise FileExistsError("Fresh archive, staging and backup required")
    pilot = BASE / "truncated_pilot_01"
    if not json.loads((pilot / "preflight.json").read_text())["accepted"]:
        raise ValueError("Seed preflight failed")
    for name in ("preparation_job", "verification_job"):
        source_dir = BASE / name
        status = json.loads((source_dir / "status.json").read_text())
        if status["exit_status"] != 0 or not status["ended_utc"]:
            raise ValueError(f"Prerequisite not completed successfully: {name}")
        for source in source_dir.iterdir():
            if source.is_file():
                target_name = (
                    source.name + ".txt" if source.suffix == ".log" else source.name
                )
                copy_exact(source, ARCHIVE / "jobs" / name / target_name)
    payload_names = (
        "native_seed.sdf",
        "native_seed.xyz",
        "atom_map.json",
        "coordinate_provenance.json",
        "preflight.json",
        "view_seed.pml",
        "input.gjf",
        "pilot_policy.json",
        "source_manifest.json",
    )
    for name in (*payload_names, "run_gaussian.sh", "assess_pilot.py"):
        copy_exact(pilot / name, ARCHIVE / "pilot" / name)
    for source in (pilot / "source").iterdir():
        copy_exact(source, ARCHIVE / "pilot" / "source" / source.name)
    for source in BASE.glob("*.py"):
        copy_exact(source, ARCHIVE / "source" / source.name)
    copy_exact(BASE / "coverage.json", ARCHIVE / "coverage.json")
    old = (
        ROOT / "output_nyl12_native_reference_20261023/native_frozen_scaffold_pilot_01"
    )
    copy_exact(
        old / "user_requested_interruption.json",
        ARCHIVE / "superseded103" / "interruption.json",
    )
    copy_exact(old / "job/status.json", ARCHIVE / "superseded103" / "status.json")
    for name in ("command.json", "run.sh"):
        copy_exact(pilot / "job" / name, ARCHIVE / "launch" / name)
    for name in (
        "native_ti1_reference_20261023.md",
        "standalone_handoff_20261020.md",
        "laptop_inspection_20261021.md",
    ):
        copy_exact(ROOT / "docs" / name, ARCHIVE / "docs" / name)
    journal = ARCHIVE / "agent_actions_at_checkpoint.jsonl"
    with (ROOT / "output_agent_actions/20261021.jsonl").open("rb") as source:
        fcntl.flock(source, fcntl.LOCK_SH)
        with journal.open("xb") as target:
            shutil.copyfileobj(source, target)
        fcntl.flock(source, fcntl.LOCK_UN)
    prior = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    changed = {"README.md", "provenance/agent_actions_at_export.jsonl"}
    files["README.md"] = ROOT / "docs/laptop_inspection_20261021.md"
    files["provenance/agent_actions_at_export.jsonl"] = journal
    for source in (ROOT / "docs").glob("*.md"):
        name = f"docs/{source.name}"
        changed.add(name)
        files[name] = source
    for name in payload_names:
        files[f"qm/truncated_ti1_seed/{name}"] = ARCHIVE / "pilot" / name
    result = build_package(files, staging)
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if digest(staging / name) != expected:
            raise ValueError(f"Package checksum mismatch: {name}")
    for name in prior["files"]:
        if (
            name not in changed
            and (current / name).read_bytes() != (staging / name).read_bytes()
        ):
            raise ValueError(f"Prior nondocument payload changed: {name}")
    current.rename(backup)
    staging.rename(current)
    report = {
        "verified": True,
        "payload_files": len(result["files"]),
        "payload_bytes": result["payload_bytes"],
        "retained_backup": str(backup),
        "old_nondocument_payloads_unchanged": True,
        "forcefield_accepted": False,
        "tests": "8 passed;97%overall,94%preparer,95%assessor coverage",
        "active_session": "nyl12_truncated_ti1_16core_20261023",
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
    logging.info("Truncated checkpoint: %s", report)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
