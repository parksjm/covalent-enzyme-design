"""Promote a verified trajectory snapshot to the reversible laptop package."""

from __future__ import annotations

import argparse
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


def digest(path: Path) -> str:
    """Return complete SHA256 of file contents."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Archive a new snapshot and promote it after byte-for-byte verification."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    snapshot = args.snapshot.resolve(strict=True)
    report = json.loads((snapshot / "trajectory_report.json").read_text())
    if (
        report["frames"] < 1
        or report["atoms_per_frame"] != 39
        or report["explicit_H_per_frame"] != 23
    ):
        raise ValueError("Expected a nonempty, explicit-H39 atom trajectory")
    if digest(snapshot / "optimization.log.txt") != report["source_sha256"]:
        raise ValueError("Snapshot source hash mismatch")
    archive = args.archive.absolute()
    if archive.exists():
        raise FileExistsError(archive)
    current = ROOT / "output_laptop_inspection_20261021"
    number = 1 + max(
        int(p.name.rsplit("_", 1)[1])
        for p in ROOT.glob("output_laptop_inspection_20261021_revision_*")
        if p.is_dir()
    )
    backup = ROOT / f"output_laptop_inspection_20261021_revision_{number}"
    staging = ROOT / "output_laptop_inspection_20261021_trajectory_staging"
    if staging.exists() or backup.exists():
        raise FileExistsError("Fresh staging and backup required")
    shutil.copytree(snapshot, archive / "trajectory")
    for source in snapshot.rglob("*"):
        if source.is_file() and digest(source) != digest(
            archive / "trajectory" / source.relative_to(snapshot)
        ):
            raise ValueError(f"Copy-time hash mismatch: {source}")
    code = archive / "source"
    code.mkdir()
    for name in (
        "export_optimization_trajectory.py",
        "test_trajectory_export.py",
        "package_trajectory_snapshot.py",
    ):
        shutil.copy2(Path(__file__).parent / name, code / name)
    shutil.copy2(
        ROOT / "scripts/gaussian_log_to_xyz.py", code / "gaussian_log_to_xyz.py"
    )
    docs = archive / "docs"
    docs.mkdir()
    for name in ("laptop_inspection_20261021.md", "standalone_handoff_20261020.md"):
        shutil.copy2(ROOT / "docs" / name, docs / name)
    journal = archive / "agent_actions_at_checkpoint.jsonl"
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
    prefix = "qm/waterfree_ti1_optimization"
    for name in (
        "optimization_trajectory.xyz",
        "optimization_trajectory.pdb",
        "trajectory_metrics.csv",
        "trajectory_report.json",
        "view_optimization.pml",
        "README.md",
    ):
        files[f"{prefix}/{name}"] = archive / "trajectory" / name
    package = build_package(files, staging)
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if digest(staging / name) != expected:
            raise ValueError(f"Package hash mismatch: {name}")
    for name in prior["files"]:
        if (
            name not in changed
            and not name.startswith(prefix + "/")
            and (current / name).read_bytes() != (staging / name).read_bytes()
        ):
            raise ValueError(f"Prior payload changed: {name}")
    current.rename(backup)
    staging.rename(current)
    verification = {
        "verified": True,
        "frames": report["frames"],
        "payload_files": len(package["files"]),
        "payload_bytes": package["payload_bytes"],
        "retained_backup": str(backup),
        "old_nondocument_payloads_preserved": True,
        "geometry_acceptance_claimed": False,
    }
    (archive / "package_verification.json").write_text(
        json.dumps(verification, indent=2) + "\n"
    )
    hashes = {
        str(p.relative_to(archive)): digest(p)
        for p in archive.rglob("*")
        if p.is_file()
    }
    (archive / "source_manifest.json").write_text(json.dumps(hashes, indent=2) + "\n")
    logging.info("Trajectory inspection package: %s", verification)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
