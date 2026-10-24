"""Freeze the mechanistic-study proposal and verify its small-file inspection export."""

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


def main() -> None:
    """Preserve source documents, journal snapshot and a reversible download."""
    archive = ROOT / "provenance/mechanistic_project_scope_20261024"
    archive.mkdir(exist_ok=False)
    sources = {
        Path(__file__): "source/checkpoint.py",
        Path(
            "/tmp/finalize_mechanistic_scope.py"
        ): "source/finalize_mechanistic_scope.py",
    }
    for name in (
        "standalone_handoff_20261020.md",
        "nyl12_mechanistic_project_scope_20261024.md",
    ):
        sources[ROOT / "docs" / name] = "docs/" + name
    for p in (ROOT / "manuscript").iterdir():
        if p.is_file() and p.suffix in {".md", ".bib"}:
            sources[p] = "manuscript/" + p.name
    for source, name in sources.items():
        target = archive / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if source.read_bytes() != target.read_bytes():
            raise ValueError("Copy-time digest mismatch")
    journal = archive / "agent_actions_at_checkpoint.jsonl"
    with (ROOT / "output_agent_actions/20261021.jsonl").open("rb") as source:
        fcntl.flock(source, fcntl.LOCK_SH)
        with journal.open("xb") as target:
            shutil.copyfileobj(source, target)
        fcntl.flock(source, fcntl.LOCK_UN)
    current = ROOT / "output_laptop_inspection_20261021"
    prior = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    changed = {"README.md", "provenance/agent_actions_at_export.jsonl"}
    for category in ("docs", "manuscript"):
        for p in (archive / category).iterdir():
            name = category + "/" + p.name
            files[name] = p
            changed.add(name)
    guide = archive / "download_README.md"
    guide.write_text(
        "# Proposed PA6 mechanistic study, 2026-10-24\n\nSee docs/nyl12_mechanistic_project_scope_20261024.md for the proposed staged DFTB QM/MM study, prior PA6 MD provenance and NylB comparison. No new simulation or enzyme pivot. Updated references R34/R35 are included.\n\n"
        + (current / "README.md").read_text()
    )
    files["README.md"] = guide
    files["provenance/agent_actions_at_export.jsonl"] = journal
    staging = ROOT / "output_laptop_inspection_20261021_mechanistic_scope_staging"
    backup = ROOT / "output_laptop_inspection_20261021_revision_29"
    if staging.exists() or backup.exists():
        raise FileExistsError("Fresh package paths required")
    report = build_package(files, staging)
    for name in prior["files"]:
        if (
            name not in changed
            and (current / name).read_bytes() != (staging / name).read_bytes()
        ):
            raise ValueError(f"Unrelated artifact changed:{name}")
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        if hashlib.sha256((staging / name).read_bytes()).hexdigest() != digest:
            raise ValueError("Export digest mismatch")
    current.rename(backup)
    staging.rename(current)
    (archive / "package_verification.json").write_text(
        json.dumps(
            {
                "verified": True,
                "payload_files": len(report["files"]),
                "bytes": report["payload_bytes"],
                "retained_backup": str(backup),
            },
            indent=2,
        )
        + "\n"
    )
    (archive / "source_manifest.json").write_text(
        json.dumps(
            {
                str(p.relative_to(archive)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in archive.rglob("*")
                if p.is_file()
            },
            indent=2,
        )
        + "\n"
    )
    logging.info(
        "Mechanistic-scope snapshot and unchanged unrelatedpayloads verified"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
