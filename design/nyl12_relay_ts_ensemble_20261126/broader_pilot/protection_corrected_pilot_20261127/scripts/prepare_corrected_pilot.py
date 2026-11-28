#!/usr/bin/env python3
"""Freeze a new 24-sample pilot using the current sequence-protection census."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(ROOT / "src"))
CAMPAIGN = ROOT / "design/nyl12_relay_ts_ensemble_20261126"
PILOT = CAMPAIGN / "broader_pilot/protection_corrected_pilot_20261127"
LEGACY = CAMPAIGN / "broader_pilot/four_template_20261127/generation_run01"
LEGACY_MANIFEST = LEGACY / "manifest.json"
CENSUS_PARENTS = [
    CAMPAIGN / "templates/medoid_01_delivery_full_010_frame_068/complex_heavy.pdb",
    CAMPAIGN / "templates/medoid_02_delivery_full_011_frame_006/complex_heavy.pdb",
    CAMPAIGN / "templates/medoid_03_delivery_full_007_frame_095/complex_heavy.pdb",
]


def sha256(path: Path) -> str:
    """Return the digest of an exact file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path: Path, value: Any) -> None:
    """Write one immutable JSON artifact and refuse replacement."""
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def copy_source_tree(source: Path, destination: Path) -> None:
    """Copy a validated source tree without following hidden defaults."""
    if not source.is_dir():
        raise FileNotFoundError(source)
    shutil.copytree(source, destination)


def validate_census(protection: dict[str, Any]) -> None:
    """Require the PI-directed K156 release and conservative L235 hold."""
    if protection["counts"]["fixed_total"] != 100:
        raise ValueError("Expected 100 fixed positions in the current policy")
    if protection["counts"]["mutable_total"] != 228:
        raise ValueError("Expected 228 mutable positions in the current policy")
    if 156 not in protection["mutable"]["alpha"]:
        raise ValueError("Canonical Lys156 is not mutable")
    if 3 in protection["mutable"]["beta"]:
        raise ValueError("Canonical Leu235 must remain fixed as beta-local 3")


def current_source_snapshot() -> dict[str, Path]:
    """Copy the complete project package and validated generation/scoring code."""
    snapshot = PILOT / "source_snapshot"
    snapshot.mkdir(exist_ok=False)
    sources = {
        "project_package": (
            ROOT / "src/nyl12_covalent",
            snapshot / "src/nyl12_covalent",
        ),
        "generation_runtime": (
            CAMPAIGN / "preparation/source_snapshot_final",
            snapshot / "generation",
        ),
        "scoring_runtime": (
            CAMPAIGN / "preparation/source_snapshot_scoring",
            snapshot / "scoring",
        ),
        "runtime_support": (
            ROOT / "design/nyl12_ligandmpnn_pilot/source_snapshot/runtime_support",
            snapshot / "runtime_support",
        ),
    }
    for source, destination in sources.values():
        copy_source_tree(source, destination)
    (snapshot / "scripts").mkdir()
    runner_source = ROOT / "scripts/run_relay_broad_pilot.py"
    runner_dest = snapshot / "scripts/run_relay_broad_pilot.py"
    runner = runner_source.read_text(encoding="utf-8")
    if "2026112620 + index" not in runner:
        raise ValueError("Expected validated legacy seed schedule to patch")
    runner_dest.write_text(
        runner.replace("2026112620 + index", "2026112720 + index"),
        encoding="utf-8",
    )
    return {key: dest for key, (_, dest) in sources.items()} | {
        "runner": runner_dest,
    }


def build_manifest(
    protection: dict[str, Any], snapshot_paths: dict[str, Path]
) -> dict[str, Any]:
    """Derive a namespaced, four-parent request from the completed legacy manifest."""
    manifest = copy.deepcopy(json.loads(LEGACY_MANIFEST.read_text(encoding="utf-8")))
    request = manifest["request"]
    request.update(
        {
            "status": "prepared_waiting_root_wt_acceptance_and_resource_release",
            "campaign": "nyl12_relay_protection_corrected_broad_pilot_20261127",
            "protection_receipt": str(PILOT / "canonical_protection.json"),
            "protection_sha256": sha256(PILOT / "canonical_protection.json"),
        }
    )
    generation = request["generation"]
    generation["mask"] = protection["mutable"]
    generation["original_samples_per_parent"] = 6
    generation["remodeled_samples_per_parent"] = 6
    generation["allocation"] = "two_original_two_remodeled_six_each"
    generation["symmetry_residues"] = protection["symmetry_residues"]
    generation["symmetry_weights"] = protection["symmetry_weights"]
    generation["redesigned_residues"] = protection["redesigned_residues"]
    manifest["allocation_protocol"] = "two_original_two_remodeled_six_each_20261127"

    native = manifest["native"]
    snapshot = PILOT / "source_snapshot"
    native.update(
        {
            "approved_for_launch": False,
            "samples_per_template": 6,
            "runner": str(snapshot / "scripts/run_relay_broad_pilot.py"),
            "generation_engine": str(snapshot / "generation/ligandmpnn/run.py"),
            "score_engine": str(snapshot / "scoring/ligandmpnn/score.py"),
            "score_wrapper": str(
                snapshot / "scoring/scripts/score_relay_ensemble_sequence.py"
            ),
            "support": str(snapshot / "runtime_support"),
            "mask": protection["mutable"],
            "symmetry_residues": generation["symmetry_residues"],
            "symmetry_weights": generation["symmetry_weights"],
            "redesigned_residues": generation["redesigned_residues"],
            "environment": {
                **native["environment"],
                "PYTHONPATH": str(snapshot / "runtime_support"),
                "CUDA_VISIBLE_DEVICES": "ROOT_RESOURCE_RELEASE_REQUIRED",
            },
        }
    )
    accepted = manifest["accepted_manifests"]
    if len(accepted) != 1 or len(accepted[0]["parents"]) != 2:
        raise ValueError("Expected the two independently accepted 5A remodel parents")
    parent_sources = {p["parent_source_id"] for p in accepted[0]["parents"]}
    if parent_sources != {"medoid01", "medoid03"}:
        raise ValueError("Accepted remodels must retain both independent originals")
    source_files: set[Path] = set()
    for _, destination in snapshot_paths.items():
        if destination.is_file():
            source_files.add(destination)
        else:
            source_files.update(p for p in destination.rglob("*") if p.is_file())
    source_files.update(CENSUS_PARENTS)
    source_files.update(Path(p) for p in request["original_parent_pdbs"].values())
    source_files.add(Path(native["ES_template"]["pdb"]))
    for parent in accepted[0]["parents"]:
        source_files.update(
            Path(parent[field])
            for field in (
                "pdb_path",
                "canonical_mapping",
                "pa6_graph_receipt",
                "scientific_acceptance",
                "authoritative_coordinates",
                "explicit_graph_cif",
            )
        )
    source_files.update(
        {
            Path(native["python"]),
            Path(native["python"]).parent.parent / "pyvenv.cfg",
            Path(native["checkpoint"]),
            ROOT / "scripts/record_agent_action.py",
            CAMPAIGN
            / "backbone_qualification/accepted_parent_manifest_adapter_20261127.json",
            CAMPAIGN
            / "backbone_qualification/qualification_terminal_adapter_20261127.json",
            PILOT / "canonical_protection.json",
        }
    )
    manifest["source_sha256"] = {
        str(path): sha256(path) for path in sorted(source_files)
    }
    native["source_sha256"] = manifest["source_sha256"]
    return manifest


def prepare() -> Path:
    """Create the immutable new census, snapshots, manifest, and launch instructions."""
    PILOT.mkdir(parents=True, exist_ok=True)
    if any(
        (PILOT / name).exists()
        for name in ("canonical_protection.json", "manifest.json", "source_snapshot")
    ):
        raise FileExistsError("Pilot artifacts already exist; refusing overwrite")
    (PILOT / "inputs").mkdir(exist_ok=True)
    from nyl12_covalent.relay_broad_pilot import protection_census, symmetry_manifest

    protection = protection_census(CENSUS_PARENTS)
    protection["source_sha256"] = {str(path): sha256(path) for path in CENSUS_PARENTS}
    protection["method"] = {
        "rule": "current shared hard-protection policy plus recurring inter-repeat contacts",
        "interface_distance_angstrom": 4.0,
        "recurrence": "at least 3 of 4 equivalent copies in each of all 3 original TS parents",
        "interface_protection_is_operational_not_proven_essential": True,
        "Leu235": "retained conservatively as beta-local3 mature N-terminal context; not claimed essential",
        "Lys156": "mutable under PI direction; no catalytic-essentiality claim",
    }
    ties = symmetry_manifest(protection["mutable"])
    protection.update(ties)
    validate_census(protection)
    protection_path = PILOT / "canonical_protection.json"
    write_new(protection_path, protection)
    snapshot_paths = current_source_snapshot()
    manifest = build_manifest(protection, snapshot_paths)
    manifest_path = PILOT / "manifest.json"
    write_new(manifest_path, manifest)
    binding = {
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "campaign_id_prefix": "relay_corr01_",
        "raw_sample_id_prefix": "corr01_",
        "protection_sha256": sha256(protection_path),
        "source_sha256_count": len(manifest["source_sha256"]),
        "root_wt_release": str(PILOT / "root_wt_acceptance_release.json"),
        "resource_preflight": str(PILOT / "root_resource_preflight.json"),
        "launch_command": f"uv run python {PILOT / 'scripts/run_corrected_pilot.py'} {manifest_path} {PILOT / 'run01'}",
    }
    write_new(PILOT / "source_binding.json", binding)
    (PILOT / "README.md").write_text(
        "# Protection-corrected broad pilot\n\n"
        "Prepared only. Root WT and resource receipts are required before launch.\n"
        "This run uses four parent backbones, six raw samples per parent (24 total), "
        "current 100-fixed/228-mutable census, four-copy symmetry, and fresh IDs.\n"
        "Lys156 is mutable. Leu235 (beta local 3) remains conservatively fixed as "
        "mature beta N-terminal context; this is not evidence of essentiality.\n\n"
        f"Manifest SHA256: `{binding['manifest_sha256']}`\n\n"
        f"Launch command (gated): `{binding['launch_command']}`\n",
        encoding="utf-8",
    )
    return manifest_path


if __name__ == "__main__":
    print(prepare())
