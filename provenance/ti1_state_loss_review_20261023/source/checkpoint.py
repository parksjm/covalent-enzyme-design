"""Preserve the stopped-state review and refresh small inspection documents."""

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
    """Verify an additive, reversible update containing no new QM calculation."""
    review = Path(__file__).resolve().parent
    archive = ROOT / "provenance/ti1_state_loss_review_20261023"
    current = ROOT / "output_laptop_inspection_20261021"
    backup = ROOT / "output_laptop_inspection_20261021_revision_21"
    staging = ROOT / "output_laptop_inspection_20261021_state_loss_staging"
    if any(p.exists() for p in (archive, backup, staging)):
        raise FileExistsError("Fresh checkpoint directories required")
    archive.mkdir()
    sources = {
        review / "report.json": "report.json",
        review / "review.py": "source/review.py",
        Path(__file__): "source/checkpoint.py",
        review / "job/status.json": "analysis_job/status.json",
        review / "job/command.json": "analysis_job/command.json",
        review / "job/run.sh": "analysis_job/run.sh",
        review / "job/job.log": "analysis_job/stdout.txt",
    }
    pilot = ROOT / "output_nyl12_waterfree_reference_20261023/waterfree_pilot_01"
    sources[pilot / "user_requested_interruption.json"] = (
        "stopped_job/interruption.json"
    )
    sources[pilot / "job/status.json"] = "stopped_job/status.json"
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
            raise ValueError("Copy-time mismatch")
    recorded_inputs = [
        pilot / "optimization.log",
        pilot / "atom_map.json",
        ROOT
        / "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation/model.prmtop",
        ROOT
        / "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation/prepared.rst7",
        ROOT / "scripts/gaussian_log_to_xyz.py",
        ROOT / "src/nyl12_covalent/active_site_geometry.py",
        ROOT / "src/nyl12_covalent/local_repair_model.py",
    ]
    (archive / "source_inputs_at_checkpoint.json").write_text(
        json.dumps(
            {
                str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in recorded_inputs
            },
            indent=2,
        )
        + "\n"
    )
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
    files["reports/parameters/ti1_state_loss_review.json"] = archive / "report.json"
    report = build_package(files, staging)
    for name in prior["files"]:
        if (
            name not in changed
            and (current / name).read_bytes() != (staging / name).read_bytes()
        ):
            raise ValueError(f"Prior payload changed: {name}")
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if hashlib.sha256((staging / name).read_bytes()).hexdigest() != expected:
            raise ValueError("Package hash mismatch")
    current.rename(backup)
    staging.rename(current)
    result = {
        "verified": True,
        "payload_files": len(report["files"]),
        "payload_bytes": report["payload_bytes"],
        "retained_backup": str(backup),
        "new_QM_launched": False,
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
    logging.info("Stopped-state checkpoint verified: %s", result)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
