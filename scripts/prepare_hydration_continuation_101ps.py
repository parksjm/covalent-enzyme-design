#!/usr/bin/env python3
"""Prepare a reviewed, conditional 51-to-101 ps PA6 hydration continuation."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.hydration_continuation_101ps import (  # noqa: E402
    candidate_found,
    validate_prefix,
)
from nyl12_covalent.hydration_continuation_sources_101ps import (  # noqa: E402
    patch_native_sources,
    patch_postprocess_sources,
)

BASE_RELATIVE = Path(
    "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
)
SEGMENTS = tuple(f"extend_{i:02d}" for i in range(10, 20))
SEEDS = (
    2026112610,
    2026112611,
    2026112612,
    2026112613,
    2026112614,
    2026112615,
    2026112616,
    2026112617,
    2026112618,
    2026112619,
)
LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare(base: Path) -> dict[str, Any]:
    """Create pinned source snapshots and a review-only 51-to-101 ps plan."""
    base = base.resolve(strict=True)
    outputs = [
        base / name
        for name in (
            "continuation_101ps_source",
            "continuation_101ps_postprocess_source",
            "continuation_101ps_template.json",
            "continuation_101ps_config.json",
            "continuation_101ps_plan.json",
            "continuation_101ps_plan.sha256",
        )
    ]
    if any(path.exists() for path in outputs):
        raise FileExistsError(
            "A continuation_101ps preparation artifact already exists"
        )
    if validate_prefix(base, 10, 19) != 0:
        raise ValueError("Continuation blocks already attempted; use a new version")
    if candidate_found(base, 9):
        raise ValueError("A near-attack water exists; stop for PI review")
    native_hashes = patch_native_sources(
        base / "followon_source_snapshot", base / "continuation_101ps_source"
    )
    post_hashes = patch_postprocess_sources(
        base / "followon_postprocess_source",
        base / "continuation_101ps_postprocess_source",
    )
    # Align the retained receipt copies with the reviewed postprocess snapshot.
    retained = base / "continuation_101ps_source/postprocess"
    for rel in post_hashes:
        if rel.startswith("src/nyl12_covalent/"):
            shutil.copy2(
                base / "continuation_101ps_postprocess_source" / rel,
                retained / Path(rel).name,
            )
        elif rel.startswith("scripts/"):
            shutil.copy2(
                base / "continuation_101ps_postprocess_source" / rel,
                retained / Path(rel).name,
            )
    native_hashes = {
        str(path.relative_to(base / "continuation_101ps_source")): sha256(path)
        for path in sorted((base / "continuation_101ps_source").rglob("*.py"))
    }
    source_hashes = {str(p): sha256(p) for p in sorted(retained.glob("*.py"))}
    post_paths = {
        str((base / "continuation_101ps_postprocess_source" / rel).resolve()): digest
        for rel, digest in post_hashes.items()
        if rel.startswith("src/nyl12_covalent/") or rel.startswith("scripts/")
    }
    cfg = json.loads((base / "followon_runtime_config.json").read_text())
    restart = base / "extend_09/restart.ncrst"
    cfg.update(
        endpoint_restart_path=str(restart.resolve()),
        endpoint_restart_sha256=sha256(restart),
        segment_names=list(SEGMENTS),
        seeds=dict(zip(SEGMENTS, SEEDS, strict=True)),
        source_snapshot_sha256=native_hashes,
        postprocess_source_sha256=post_paths,
        postprocess_argv=[
            sys.executable,
            str(
                (
                    base
                    / "continuation_101ps_postprocess_source/scripts/refresh_deacylation_hydration.py"
                ).resolve()
            ),
            "--base",
            str(base.resolve()),
            "--project-root",
            str((base / "continuation_101ps_postprocess_source").resolve()),
        ],
        postprocess_environment={
            "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        },
    )
    (base / "continuation_101ps_template.json").write_text(
        json.dumps(cfg, indent=2) + "\n"
    )
    plan = {
        "status": "PREPARED_ONLY_NOT_LAUNCHED",
        "requires_root_review": True,
        "root": str(ROOT),
        "base": str(base),
        "initial_required_ps": 51,
        "target_total_ps": 101,
        "segment_names": list(SEGMENTS),
        "seeds": cfg["seeds"],
        "poll_seconds": 120,
        "wait_seconds": 1,
        "scope": "Unbiased NVT; one CPU; no reheating; original protonation, QM partition, Hamiltonian, active region, topology and MM-water model retained.",
        "water_handling": "No water insertion or promotion; all solvent waters remain MM. Any near-attack result is only a position candidate and requires separate QM promotion and reoptimization.",
        "candidate_gate": {
            "O_C_A_max": 3.5,
            "attack_angle_deg": [95, 125],
            "action": "stop and request PI review after refreshing analysis",
        },
        "native_source_sha256": native_hashes,
        "postprocess_source_sha256": post_paths,
        "retained_postprocess_copy_sha256": source_hashes,
        "pinned_sha256": {},
        "postprocess_argv": cfg["postprocess_argv"],
        "postprocess_environment": cfg["postprocess_environment"],
        "recovery_source_sha256": post_paths,
    }
    pinned = [
        base / "config.json",
        base / "hydration_endpoint_status_20261124.json",
        base / "extension_20ps_config.json",
        base / "extension_20ps_status.json",
        base / "followon_runtime_config.json",
        base / "continuation_101ps_template.json",
        base / "followon_completion.json",
        base / "followon_queue_resume_job/execution/status.json",
        restart,
        base / "inputs/source.mdin",
        base / "inputs/system.prmtop",
    ]
    for index in range(1, 10):
        block = base / f"extend_{index:02d}"
        pinned.extend(
            block / rel
            for rel in (
                "status.json",
                "process.json",
                "scientific_acceptance.json",
                "restart.ncrst",
                "inputs/segment.mdin",
                "inputs/source.mdin",
                "inputs/system.prmtop",
                "inputs/restart.ncrst",
            )
        )
    pinned.extend(path for path in (base / "continuation_101ps_source").rglob("*.py"))
    pinned.extend(
        path for path in (base / "continuation_101ps_postprocess_source").rglob("*.py")
    )
    plan["pinned_sha256"] = {str(path.resolve()): sha256(path) for path in pinned}
    plan["initial_analysis_manifest_sha256"] = sha256(
        base / "analysis/source_manifest.json"
    )
    plan_path = base / "continuation_101ps_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    plan_hash = sha256(plan_path)
    plan_path.with_suffix(".sha256").write_text(plan_hash + "\n")
    return {
        "status": plan["status"],
        "plan": str(plan_path),
        "plan_sha256": plan_hash,
        "segments": list(SEGMENTS),
    }


def run(plan_path: Path, reviewed_sha256: str) -> None:
    """Run the frozen queue only when the caller supplies its reviewed hash."""
    plan_path = plan_path.resolve(strict=True)
    actual = sha256(plan_path)
    expected = plan_path.with_suffix(".sha256").read_text().strip()
    if actual != expected or reviewed_sha256 != actual:
        raise ValueError("Root-reviewed plan SHA-256 is required")
    plan = json.loads(plan_path.read_text())
    if plan.get("status") != "PREPARED_ONLY_NOT_LAUNCHED":
        raise ValueError("Plan is not in review-only state")
    snapshot_runner = (
        Path(plan["base"])
        / "continuation_101ps_source/scripts/run_hydration_followon_queue.py"
    )
    subprocess.run(
        [sys.executable, str(snapshot_runner), "--plan", str(plan_path)],
        cwd=plan["root"],
        check=True,
    )


def main() -> None:
    """Prepare a plan or run it after explicit hash review."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / BASE_RELATIVE)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-plan", type=Path)
    parser.add_argument("--reviewed-sha256")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.prepare == (args.run_plan is not None):
        parser.error("Choose exactly one of --prepare or --run-plan")
    if args.prepare:
        LOGGER.info(
            "Prepared review-only continuation: %s", json.dumps(prepare(args.base))
        )
    else:
        if not args.reviewed_sha256:
            parser.error("--reviewed-sha256 is required with --run-plan")
        run(args.run_plan, args.reviewed_sha256)


if __name__ == "__main__":
    main()
