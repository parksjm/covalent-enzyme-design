"""Freeze the PA6 ES reuse audit and rebuild its verified small-file export."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import logging
import shutil
import sys
from pathlib import Path
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
MAX_BYTES = 25_000_000
ENGINE = ROOT / "output_nyl12_dftb_engine_checks_20261024"
RUN = ROOT / "output_nyl12_pa6_es_audit_20261024/reuse_audit_attempt1"
HYDRATION = ROOT / "output_nyl12_pa6_hydration_20261024"
sys.path.insert(0, str(ROOT))
from scripts.prepare_nyl12_inspection_download import build_package  # noqa: E402

SOURCE_FILES = (
    "configs/pa6_wt_es_reuse_20261024.json",
    "scripts/audit_nyl12_pa6_md_reuse.py",
    "scripts/inspect_nyl12_pa6_es_seeds.py",
    "src/nyl12_covalent/es_snapshot_selection.py",
    "src/nyl12_covalent/es_topology_audit.py",
    "src/nyl12_covalent/es_water_inventory.py",
    "src/nyl12_covalent/es_seed_display.py",
    "src/nyl12_covalent/es_seed_legacy_geometry.py",
    "tests/test_audit_nyl12_pa6_md_reuse.py",
    "tests/test_inspect_nyl12_pa6_es_seeds.py",
    "tests/test_es_snapshot_selection.py",
    "tests/test_es_topology_audit.py",
    "tests/test_es_water_inventory.py",
)
RUN_SOURCE_FILES = (
    "audit_nyl12_pa6_md_reuse.py",
    "es_topology_audit.py",
    "es_snapshot_selection.py",
    "pa6_wt_es_reuse_20261024.json",
)


def digest(path: Path) -> str:
    """Return a small artifact's SHA256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    """Write an exclusively new JSON record."""
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")


def copy_checked(source: Path, target: Path) -> None:
    """Copy a small source into a fresh destination and verify its digest."""
    if target.exists():
        raise FileExistsError(target)
    if not source.is_file() or source.stat().st_size > MAX_BYTES:
        raise ValueError(f"Missing or oversized source: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    before = digest(source)
    shutil.copy2(source, target)
    if digest(target) != before or digest(source) != before:
        raise ValueError(f"Copy-time digest mismatch: {source}")


def main() -> None:
    """Freeze actual and current sources and preserve the previous export."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inspection", required=True, type=Path)
    args = parser.parse_args()
    inspection = args.inspection.resolve(strict=True)
    expected_inspection = (
        ROOT / "output_nyl12_pa6_es_audit_20261024/seed_inspection_attempt4"
    ).resolve(strict=True)
    if inspection != expected_inspection:
        raise ValueError("This checkpoint requires the accepted attempt4 inspection")
    archive = ROOT / "provenance/pa6_es_reuse_checkpoint_20261024"
    current = ROOT / "output_laptop_inspection_20261021"
    staging = ROOT / "output_laptop_inspection_20261021_pa6_es_reuse_staging"
    backup = ROOT / "output_laptop_inspection_20261021_revision_30"
    for fresh in (archive, staging, backup):
        if fresh.exists():
            raise FileExistsError(f"Fresh path required: {fresh}")

    inspection_report = json.loads((inspection / "inspection.json").read_text())
    names = [frame["name"] for frame in inspection_report["frames"]]
    if len(names) != 9 or len(set(names)) != 9:
        raise ValueError("Exactly nine distinct inspected ES frames required")
    inspection_files = [f"{name}.pdb" for name in names] + [
        "fullprotein_substrate.pdb",
        "view_es_seeds.py",
        "view_es_seeds.pml",
        "README.txt",
        "inspection.json",
        "output_hashes.json",
    ]
    # Pocket mappings are useful; the large context mapping is excluded.
    inspection_files += [f"{name}_mapping.json" for name in names]
    output_records = json.loads((inspection / "output_hashes.json").read_text())
    recorded = {record["path"] for record in output_records}
    if not (set(inspection_files) - {"output_hashes.json"}).issubset(recorded):
        raise ValueError("Inspection manifest omits requested artifacts")
    for record in output_records:
        name = record["path"]
        if name in inspection_files and digest(inspection / name) != record["sha256"]:
            raise ValueError(f"Inspection output digest mismatch: {name}")
    summary_records = json.loads((inspection / "summary_provenance.json").read_text())
    summary_required = {
        str(ROOT / "output_nyl12_pa6_es_audit_20261024/export_seed_summary.py"),
        str(inspection / "inspection.json"),
        str(inspection / "seed_summary.csv"),
        str(inspection / "seed_summary.md"),
    }
    if {record["path"] for record in summary_records} != summary_required:
        raise ValueError("Summary provenance must identify exact source and summaries")
    for record in summary_records:
        path = Path(record["path"])
        if not path.is_absolute() or digest(path) != record["sha256"]:
            raise ValueError(f"Summary provenance mismatch: {path}")
    acceptance = json.loads((inspection / "inspection_acceptance.json").read_text())
    if (
        acceptance["process_exit_status"] != 0
        or acceptance["nine_seed_geometry_source_precision_checks_passed"] is not True
        or acceptance["source_output_hashes_verified"] is not True
        or acceptance["inspection_sha256"] != digest(inspection / "inspection.json")
    ):
        raise ValueError("Inspection acceptance or current inspection hash failed")
    inspection_files += [
        "seed_summary.csv",
        "seed_summary.md",
        "summary_provenance.json",
        "inspection_acceptance.json",
    ]

    sources = {ROOT / name: "source/" + name for name in SOURCE_FILES}
    sources[Path(__file__).resolve()] = "source/checkpoint_pa6_es_reuse.py"
    sources[ROOT / "output_nyl12_pa6_es_audit_20261024/export_seed_summary.py"] = (
        "source/export_seed_summary.py"
    )
    for name in (
        "pa6_es_reuse_audit_20261024.md",
        "pa6_dftb_validation_pilot_20261024.md",
        "pa6_es_hydration_20261024.md",
        "standalone_handoff_20261020.md",
    ):
        sources[ROOT / "docs" / name] = "docs/" + name
    for path in sorted((ROOT / "manuscript").iterdir()):
        if path.is_file() and path.suffix in {".md", ".bib"}:
            sources[path] = "manuscript/" + path.name
    for name in RUN_SOURCE_FILES:
        sources[RUN / "source" / name] = "run_source/" + name
    sources[RUN / "reuse_audit.json"] = "reports/reuse_audit.json"
    for name in ("comparison.json", "engine_parameters.json", "source_manifest.json"):
        sources[ENGINE / name] = "engine/" + name
    for path in sorted((ENGINE / "source").rglob("*")):
        if path.is_file():
            sources[path] = "engine_source/" + str(path.relative_to(ENGINE / "source"))
    inspection_source = inspection / "source"
    if not inspection_source.is_dir():
        raise FileNotFoundError(
            f"Inspection source directory required: {inspection_source}"
        )
    for path in sorted(inspection_source.rglob("*")):
        if path.is_file():
            sources[path] = "inspection_source/" + str(
                path.relative_to(inspection_source)
            )
    for name in inspection_files:
        sources[inspection / name] = "inspection/" + name
    hydration_acceptance = json.loads(
        (HYDRATION / "scientific_acceptance.json").read_text()
    )
    positive_control = json.loads((HYDRATION / "positive_control.json").read_text())
    if (
        hydration_acceptance["accepted"] is not True
        or positive_control["accepted"] is not True
        or positive_control["process_exit_status"] != 0
        or not all(
            item["passed"] is True for item in positive_control["comparison"].values()
        )
    ):
        raise ValueError("Hydration analysis and positive control must be accepted")
    hydration_reports = (
        "hydration_summary.json",
        "scientific_acceptance.json",
        "source_manifest.json",
        "positive_control.json",
    )
    for name in hydration_reports:
        sources[HYDRATION / name] = "hydration/" + name
    for name in (
        "hydration_audit.py",
        "inspection.json",
        "pa6_wt_es_reuse_20261024.json",
        "positive_control.py",
        "positive_control_triclinic.py",
        "reuse_audit.json",
    ):
        sources[HYDRATION / "source" / name] = "hydration_source/" + name
    for replica in ("replica_01", "replica_02", "replica_03"):
        for name in (
            "watershell.cpptraj",
            "cpptraj.log",
            "n35.dat",
            "og35.dat",
            "union35.dat",
            "union4.dat",
            "hydration_counts.csv",
            "process_status.json",
        ):
            relative = f"attempts/{replica}/{name}"
            sources[HYDRATION / relative] = "hydration_evidence/" + relative
    for name in (
        "positive_control.cpptraj",
        "cpptraj.log",
        "lower35.dat",
        "lower4.dat",
    ):
        relative = "positive_control_attempt2/" + name
        sources[HYDRATION / relative] = "hydration_evidence/" + relative
    for job in (
        "full_hydration_attempt1",
        "positive_control_attempt1",
        "positive_control_attempt2",
    ):
        for name in ("status.json", "command.json", "run.sh", "job.log"):
            relative = f"jobs/{job}/{name}"
            sources[HYDRATION / relative] = "hydration_evidence/" + relative
    for name in ("status.json", "command.json", "job.log", "run.sh"):
        sources[
            ROOT / "output_nyl12_pa6_es_audit_20261024/jobs/full_checks_complete" / name
        ] = "verification/full_checks_complete/" + name
    for source in sources:
        if not source.is_file() or source.stat().st_size > MAX_BYTES:
            raise ValueError(f"Missing or oversized checkpoint source: {source}")

    prior = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    prior_digests = {name: digest(current / name) for name in prior["files"]}
    archive.mkdir()
    for source, name in sources.items():
        copy_checked(source, archive / name)
    journal = archive / "agent_actions_at_checkpoint.jsonl"
    with (ROOT / "output_agent_actions/20261021.jsonl").open("rb") as source:
        fcntl.flock(source, fcntl.LOCK_SH)
        try:
            with journal.open("xb") as target:
                shutil.copyfileobj(source, target)
        finally:
            fcntl.flock(source, fcntl.LOCK_UN)

    changed = {"README.md", "provenance/agent_actions_at_export.jsonl"}
    for category in ("docs", "manuscript"):
        for path in (archive / category).iterdir():
            name = category + "/" + path.name
            files[name] = path
            changed.add(name)
    for name in inspection_files:
        destination = "md/pa6_wt_es_seeds/" + name
        files[destination] = archive / "inspection" / name
        changed.add(destination)
    for path in (archive / "reports").iterdir():
        destination = "md/pa6_wt_es_seeds/" + path.name
        files[destination] = path
        changed.add(destination)
    for path in (archive / "engine").iterdir():
        destination = "md/pa6_wt_es_seeds/engine/" + path.name
        files[destination] = path
        changed.add(destination)
    for name in hydration_reports:
        destination = "md/pa6_wt_es_seeds/hydration/" + name
        files[destination] = archive / "hydration" / name
        changed.add(destination)
    for name in ("lower35.dat", "lower4.dat"):
        destination = "md/pa6_wt_es_seeds/hydration/positive_control/" + name
        files[destination] = (
            archive / "hydration_evidence/positive_control_attempt2" / name
        )
        changed.add(destination)
    guide = archive / "download_README.md"
    guide.write_text(
        "# PA6 ES reuse checkpoint, 2026-10-24\n\n"
        "Nine retained classical ES snapshots are actual trajectory frames, "
        "not optimized TI/TS structures. See docs/pa6_es_reuse_audit_20261024.md.\n"
        "In PyMOL: run /absolute/path/to/md/pa6_wt_es_seeds/view_es_seeds.py\n"
        "Replace /absolute/path/to with your downloaded package directory.\n\n"
        "All 6000 saved frames have no water oxygen within 4 A of the tested "
        "catalytic-site selections in this fixed ES model at 50 ps cadence; "
        "this is not proof of a catalytic mechanism.\n"
        "All original hydrogenated solvent coordinates remain retained on "
        "the workstation.\n\n" + (current / "README.md").read_text().rstrip() + "\n"
    )
    files["README.md"] = guide
    files["provenance/agent_actions_at_export.jsonl"] = journal
    if sum(path.stat().st_size for path in files.values()) > MAX_BYTES:
        raise ValueError("Combined export payload exceeds 25 MB")
    report = build_package(files, staging)
    for name, before in prior_digests.items():
        if digest(current / name) != before:
            raise ValueError(f"Existing package changed during rebuild: {name}")
        if name not in changed and digest(staging / name) != before:
            raise ValueError(f"Unrelated artifact changed: {name}")
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if digest(staging / name) != expected:
            raise ValueError(f"Export digest mismatch: {name}")
    actual_bytes = sum(
        path.stat().st_size for path in staging.rglob("*") if path.is_file()
    )
    if actual_bytes > MAX_BYTES:
        raise ValueError("Complete export including manifests exceeds 25 MB")
    current.rename(backup)
    staging.rename(current)
    write_json(
        archive / "package_verification.json",
        {
            "verified": True,
            "payload_files": len(report["files"]),
            "payload_bytes": report["payload_bytes"],
            "complete_package_bytes": actual_bytes,
            "retained_backup": str(backup),
            "inspection_source": str(inspection),
            "checkpoint_script_sha256": digest(Path(__file__)),
            "prior_artifacts_verified": len(prior_digests),
        },
    )
    write_json(
        archive / "source_manifest.json",
        {
            str(path.relative_to(archive)): digest(path)
            for path in sorted(archive.rglob("*"))
            if path.is_file()
        },
    )
    for name, expected in json.loads(
        (archive / "source_manifest.json").read_text()
    ).items():
        if digest(archive / name) != expected:
            raise ValueError(f"Archive digest mismatch: {name}")
    logging.info("ES checkpoint verified; retained previous package at %s", backup)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
