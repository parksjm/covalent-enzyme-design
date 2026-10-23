"""Archive the reviewed donor cluster and update the curated laptop payload."""

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
    """Preserve source records and verify an additive inspection export."""
    base = ROOT / "output_nyl12_oxyanion_cluster_20261023"
    seed = base / "seed_01"
    archive = ROOT / "provenance/ti1_oxyanion_cluster_20261023"
    current = ROOT / "output_laptop_inspection_20261021"
    backup = ROOT / "output_laptop_inspection_20261021_revision_22"
    staging = ROOT / "output_laptop_inspection_20261021_oxyanion_staging"
    if any(p.exists() for p in (archive, backup, staging)):
        raise FileExistsError("Fresh checkpoint and package paths required")
    if not seed.is_dir():
        raise FileNotFoundError(seed)
    archive.mkdir()
    sources = {Path(__file__): "source/checkpoint.py"}
    for p in base.glob("*.py"):
        sources[p] = "source/" + p.name
    for p in seed.iterdir():
        if p.is_file() and p.suffix in {".json", ".sdf", ".xyz", ".pml"}:
            sources[p] = "seed/" + p.name
    for job in [*base.glob("*job"), *(base / "jobs").iterdir()]:
        for name in ("status.json", "command.json", "run.sh", "job.log"):
            p = job / name
            if p.is_file():
                sources[p] = (
                    "jobs/"
                    + job.name
                    + "/"
                    + ("stdout.txt" if name == "job.log" else name)
                )
    sources[base / "final_coverage.json"] = "final_coverage.json"
    pilot = base / "svp_pilot_01"
    for name in (
        "input.gjf",
        "pilot_policy.json",
        "run_gaussian.sh",
        "source_manifest.json",
    ):
        sources[pilot / name] = "pilot/" + name
    for name in ("status.json", "command.json", "run.sh"):
        sources[pilot / "job" / name] = "pilot/job/" + name
    for name in (
        "standalone_handoff_20261020.md",
        "laptop_inspection_20261021.md",
        "ti1_state_loss_and_qm_options_20261023.md",
    ):
        sources[ROOT / "docs" / name] = "docs/" + name
    for source, name in sources.items():
        target = archive / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if source.read_bytes() != target.read_bytes():
            raise ValueError(f"Copy mismatch: {source}")
    journal = archive / "agent_actions_at_checkpoint.jsonl"
    with (ROOT / "output_agent_actions/20261021.jsonl").open("rb") as source:
        fcntl.flock(source, fcntl.LOCK_SH)
        with journal.open("xb") as target:
            shutil.copyfileobj(source, target)
        fcntl.flock(source, fcntl.LOCK_UN)
    prior = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    changed = {"README.md", "provenance/agent_actions_at_export.jsonl"}
    files["README.md"] = archive / "docs/laptop_inspection_20261021.md"
    files["provenance/agent_actions_at_export.jsonl"] = journal
    for source in (archive / "docs").iterdir():
        name = "docs/" + source.name
        changed.add(name)
        files[name] = source
    for source in (archive / "seed").iterdir():
        files["qm/ti1_oxyanion_cluster/" + source.name] = source
    for source in (archive / "pilot").rglob("*"):
        if source.is_file():
            name = "qm/ti1_oxyanion_cluster/pilot/" + str(
                source.relative_to(archive / "pilot")
            )
            files[name] = source
    report = build_package(files, staging)
    for name in prior["files"]:
        if (
            name not in changed
            and (current / name).read_bytes() != (staging / name).read_bytes()
        ):
            raise ValueError(f"Unrelated payload changed: {name}")
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if hashlib.sha256((staging / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Export checksum mismatch: {name}")
    current.rename(backup)
    staging.rename(current)
    result = {
        "verified": True,
        "payload_files": len(report["files"]),
        "payload_bytes": report["payload_bytes"],
        "retained_backup": str(backup),
        "new_QM_launched": True,
        "job_status_is_snapshot": True,
    }
    (archive / "package_verification.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    manifest = {
        str(p.relative_to(archive)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in archive.rglob("*")
        if p.is_file()
    }
    (archive / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    logging.info("Donor-cluster checkpoint verified: %s", result)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
