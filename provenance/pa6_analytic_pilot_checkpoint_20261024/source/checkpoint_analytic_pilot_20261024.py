"""Freeze analytic-pilot proof separately from local unpublished coordinates."""

from __future__ import annotations
import argparse
import fcntl
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
PILOT = ROOT / "output_nyl12_pa6_dftb_pilot_20261024"
LIMIT = 25_000_000
sys.path.insert(0, str(ROOT))
from scripts.prepare_nyl12_inspection_download import build_package  # noqa: E402

DOCS = (
    "ti1_parameterization_status_20261024.md",
    "local_agent_mechanism_project_20261024.md",
    "standalone_handoff_20261020.md",
    "pa6_dftb_validation_pilot_20261024.md",
    "pa6_analytic_pilot_results_20261024.md",
)
MODULES = (
    "es_dftb_analytic_pilot",
    "es_dftb_checks",
    "es_dftb_engine",
    "es_dftb_force_evidence",
    "es_dftb_relaxation",
    "es_dftb_relaxation_gate",
    "es_qm_region",
)
TESTS = (
    "test_es_dftb_engine.py",
    "test_es_qm_region.py",
    "test_run_nyl12_es_dftb_relaxation.py",
)
CONFIGS = ("pa6_es_analytic_pilot_20261024.json", "pa6_es_dftb_pilot_20261024.json")
BASE = (
    "atom_mapping.json.gz",
    "export_report.json",
    "initial_fullprotein_ligand.pdb.gz",
    "optimization_evaluations.pdb.gz",
    "pocket_endpoint.pdb",
    "pocket_initial.pdb",
    "source_manifest.json",
    "view_analytic_pilot.py",
    "full_enzyme_substrate_endpoint.pdb.gz",
    "view_full_enzyme_endpoint.py",
)
QM = (
    "qm_subsystem_initial.xyz",
    "qm_subsystem_endpoint.xyz",
    "qm_export_manifest.json",
    "qm_subsystem_atom_map.json",
)


def sha(path: Path) -> str:
    """Hash a small artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data: Any) -> None:
    """Write an exclusively new proof record."""
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def copy_checked(source: Path, target: Path) -> None:
    """Freeze one small artifact and verify source and destination hashes."""
    if not source.is_file() or source.stat().st_size > LIMIT:
        raise ValueError(f"Missing or oversized input: {source}")
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    before = sha(source)
    shutil.copy2(source, target)
    if sha(source) != before or sha(target) != before:
        raise ValueError(f"Copy-time source change: {source}")


def verify_package(directory: Path) -> None:
    """Verify every package manifest entry."""
    for line in (directory / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if sha(directory / name) != expected:
            raise ValueError(f"Package digest mismatch: {name}")


def main() -> None:
    """Create reviewed source proof and a reversible small inspection export."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inspection", required=True, type=Path)
    parser.add_argument("--qm-export-directory", required=True, type=Path)
    parser.add_argument("--endpoint-review", required=True, type=Path)
    args = parser.parse_args()
    inspection = args.inspection.resolve(strict=True)
    qm_export = args.qm_export_directory.resolve(strict=True)
    review = args.endpoint_review.resolve(strict=True)
    proof = ROOT / "provenance/pa6_analytic_pilot_checkpoint_20261024"
    local = PILOT / "analytic_checkpoint_payload_20261024"
    current = ROOT / "output_laptop_inspection_20261021"
    staging = ROOT / "output_laptop_inspection_20261021_analytic_pilot_staging"
    backup = ROOT / "output_laptop_inspection_20261021_revision_32"
    for fresh in (proof, local, staging, backup):
        if fresh.exists():
            raise FileExistsError(f"Fresh path required: {fresh}")
    sources = {ROOT / "docs" / name: "docs/" + name for name in DOCS}
    sources.update({ROOT / "configs" / name: "configs/" + name for name in CONFIGS})
    sources.update(
        {
            ROOT / "src/nyl12_covalent" / (name + ".py"): "source/src/nyl12_covalent/"
            + name
            + ".py"
            for name in MODULES
        }
    )
    sources.update({ROOT / "tests" / name: "source/tests/" + name for name in TESTS})
    sources[ROOT / "scripts/run_nyl12_es_dftb_relaxation.py"] = (
        "source/scripts/run_nyl12_es_dftb_relaxation.py"
    )
    sources[Path(__file__).resolve()] = "source/checkpoint_analytic_pilot_20261024.py"
    for path in sorted((ROOT / "manuscript").iterdir()):
        if path.is_file() and path.suffix in {".md", ".bib"}:
            sources[path] = "manuscript/" + path.name
    for attempt in ("analytic_relaxation_attempt1", "analytic_relaxation_attempt2"):
        folder = PILOT / attempt / "source"
        if not folder.is_dir():
            raise FileNotFoundError(folder)
        for path in sorted(folder.rglob("*")):
            if path.is_file() and path.suffix in {".py", ".json", ".sh", ".in"}:
                sources[path] = f"run_source/{attempt}/" + str(path.relative_to(folder))
        for name in (
            "source_manifest.json",
            "process.json",
            "geometry_acceptance.json",
        ):
            sources[PILOT / attempt / name] = f"run_proof/{attempt}/" + name
    folder = inspection / "source"
    if not folder.is_dir():
        raise FileNotFoundError(folder)
    for path in sorted(folder.rglob("*")):
        if path.is_file() and path.suffix in {".py", ".json"}:
            sources[path] = "inspection_source/" + str(path.relative_to(folder))
    payload = {inspection / name: name for name in BASE}
    payload.update({qm_export / name: name for name in QM})
    payload[review] = "endpoint_review.json"
    for path in (*sources, *payload):
        if not path.is_file() or path.stat().st_size > LIMIT:
            raise ValueError(f"Missing or oversized checkpoint input: {path}")
    prior = json.loads((current / "provenance.json").read_text())
    if len(prior["files"]) != 306:
        raise ValueError("Prior inventory differs from reviewed306 artifacts")
    verify_package(current)
    prior_hashes = {name: sha(current / name) for name in prior["files"]}
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    proof.mkdir()
    local.mkdir()
    for path, name in sources.items():
        copy_checked(path, proof / name)
    for path, name in payload.items():
        copy_checked(path, local / name)
    write_json(
        proof / "local_payload_manifest.json",
        [
            {
                "source": str(path),
                "retained_copy": str(local / name),
                "sha256": sha(local / name),
                "bytes": (local / name).stat().st_size,
            }
            for path, name in payload.items()
        ],
    )
    journal = proof / "agent_actions_at_checkpoint.jsonl"
    with (ROOT / "output_agent_actions/20261021.jsonl").open("rb") as handle:
        fcntl.flock(handle, fcntl.LOCK_SH)
        try:
            with journal.open("xb") as target:
                shutil.copyfileobj(handle, target)
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)
    changed = {"README.md", "provenance/agent_actions_at_export.jsonl"}
    for category in ("docs", "manuscript"):
        for path in (proof / category).iterdir():
            name = category + "/" + path.name
            files[name] = path
            changed.add(name)
    for name in payload.values():
        destination = "md/pa6_es_analytic_pilot/" + name
        files[destination] = local / name
        changed.add(destination)
    for name in CONFIGS:
        destination = "md/pa6_es_analytic_pilot/reference_configs/" + name
        files[destination] = proof / "configs" / name
        changed.add(destination)
    guide = proof / "download_README.md"
    guide.write_text(
        "# Conditional analytic ES pilot, 2026-10-24\n\nSee docs/pa6_analytic_pilot_results_20261024.md and md/pa6_es_analytic_pilot/endpoint_review.json for final interpretation. This bounded analytic-gradient pilot does not establish an optimized TI/TS or validate the full mechanism.\nIn PyMOL: run /absolute/path/to/md/pa6_es_analytic_pilot/view_analytic_pilot.py\n\n"
        + (current / "README.md").read_text().rstrip()
        + "\n"
    )
    files["README.md"] = guide
    files["provenance/agent_actions_at_export.jsonl"] = journal
    if sum(path.stat().st_size for path in files.values()) > LIMIT:
        raise ValueError("Export payload exceeds25 MB")
    report = build_package(files, staging)
    verify_package(staging)
    for name, expected in prior_hashes.items():
        if sha(current / name) != expected:
            raise ValueError(f"Prior artifact changed during build: {name}")
        if name not in changed and sha(staging / name) != expected:
            raise ValueError(f"Unrelated artifact changed: {name}")
    size = sum(path.stat().st_size for path in staging.rglob("*") if path.is_file())
    if size > LIMIT:
        raise ValueError("Complete package exceeds25 MB")
    current.rename(backup)
    staging.rename(current)
    write_json(
        proof / "package_verification.json",
        {
            "verified": True,
            "payload_files": len(report["files"]),
            "complete_package_bytes": size,
            "prior_payload_count": len(prior_hashes),
            "retained_backup": str(backup),
            "local_coordinate_payload": str(local),
            "checkpoint_script_sha256": sha(Path(__file__)),
            "declared_updated_prior_paths": sorted(changed & set(prior_hashes)),
        },
    )
    manifest = {
        str(path.relative_to(proof)): sha(path)
        for path in sorted(proof.rglob("*"))
        if path.is_file()
    }
    write_json(proof / "source_manifest.json", manifest)
    if any(sha(proof / name) != expected for name, expected in manifest.items()):
        raise ValueError("Proof archive digest mismatch")
    logging.info("Analytic pilot checkpoint verified; backup retained at %s", backup)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
