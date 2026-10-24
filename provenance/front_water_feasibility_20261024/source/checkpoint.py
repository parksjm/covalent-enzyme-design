"""Checkpoint full-pocket geometry evidence and the mechanism literature review."""

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
    """Retain immutable small-file evidence and verify the curated download."""
    base = ROOT / "output_nyl12_front_water_20261024"
    archive = ROOT / "provenance/front_water_feasibility_20261024"
    archive.mkdir(exist_ok=False)
    sources = {Path(__file__): "source/checkpoint.py"}
    for p in base.glob("*.py"):
        sources[p] = "source/" + p.name
    for p in base.glob("*.json"):
        sources[p] = "reports/" + p.name
    for p in (base / "inspection").iterdir():
        if p.is_file():
            sources[p] = "inspection/" + p.name
    for job in base.glob("*job"):
        for name in ("status.json", "command.json", "run.sh", "job.log"):
            p = job / name
            if p.is_file():
                sources[p] = (
                    "jobs/"
                    + job.name
                    + "/"
                    + ("stdout.txt" if name == "job.log" else name)
                )
    for name in (
        "standalone_handoff_20261020.md",
        "front_water_feasibility_20261024.md",
        "water_relay_mechanism_evidence_20261024.md",
    ):
        sources[ROOT / "docs" / name] = "docs/" + name
    for source, name in sources.items():
        target = archive / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if source.read_bytes() != target.read_bytes():
            raise ValueError("Copy mismatch")
    journal = archive / "agent_actions_at_checkpoint.jsonl"
    with (ROOT / "output_agent_actions/20261021.jsonl").open("rb") as f:
        fcntl.flock(f, fcntl.LOCK_SH)
        with journal.open("xb") as out:
            shutil.copyfileobj(f, out)
        fcntl.flock(f, fcntl.LOCK_UN)
    current = ROOT / "output_laptop_inspection_20261021"
    prior = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    readme = archive / "download_README.md"
    readme.write_text(
        "# Latest: front-water feasibility and mechanism review\n\nNo acceptable front-side water found in the fixed full pocket. No new Gaussian job or parameter fit was launched. Inspect qm/front_water_feasibility/view_front_search.pml; both probes are rejected, not seeds. See docs/front_water_feasibility_20261024.md and docs/water_relay_mechanism_evidence_20261024.md.\n\n"
        + (current / "README.md").read_text()
    )
    files["README.md"] = readme
    files["provenance/agent_actions_at_export.jsonl"] = journal
    for source in (archive / "docs").iterdir():
        files["docs/" + source.name] = source
    for source in (archive / "inspection").iterdir():
        files["qm/front_water_feasibility/" + source.name] = source
    for source in (archive / "reports").iterdir():
        files["qm/front_water_feasibility/reports/" + source.name] = source
    staging = ROOT / "output_laptop_inspection_20261021_front_water_staging"
    backup = ROOT / "output_laptop_inspection_20261021_revision_26"
    if staging.exists() or backup.exists():
        raise FileExistsError("Fresh export paths required")
    report = build_package(files, staging)
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        if hashlib.sha256((staging / name).read_bytes()).hexdigest() != digest:
            raise ValueError("Package digest failure")
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
    logging.info("Front-water reports and literature export verified")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
