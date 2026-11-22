#!/usr/bin/env python3
"""Collect disjoint predictions, publish hashed relax requests, then compare models."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.design_component_normalization import normalize_model  # noqa: E402
from nyl12_covalent.design_geometry_comparison import compare_models  # noqa: E402
from nyl12_covalent.design_geometry_plots import plot_comparison  # noqa: E402
from nyl12_covalent.ligandmpnn_pilot_audit import _pdb_sequences  # noqa: E402


def write_json(path: Path, payload: dict) -> None:
    """Publish an atomic progress or request document."""
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temp, path)


def digest(path: Path) -> str:
    """Return the SHA256 digest of an immutable model input."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sync_tree(source: str, destination: str, options: list[str], identity: str) -> None:
    """Run one authorized, nondeleting transfer and fail explicitly on an error."""
    subprocess.run(
        [
            "rsync",
            "-rtz",
            "--safe-links",
            "-e",
            f"ssh -i {identity} -o BatchMode=yes",
            *options,
            source,
            destination,
        ],
        check=True,
        timeout=600,
    )


def component_path(base: Path, job: dict) -> Path:
    """Locate the native chemical component, including the preserved WT sentinel."""
    variant, index = job["variant"], job["model_index"]
    if variant == "WT" and index == 0:
        run = base / "predictions/WT/seed_01"
    else:
        run = base / "predictions/runs" / variant / f"model_{index}" / "output"
    return run / f"boltz_results_{variant}" / "processed/mols" / f"{variant}.pkl"


def publish_requests(base: Path, manifest: dict, config: dict) -> int:
    """Normalize ready models using their own molecular graphs and verify sequences."""
    count = 0
    for job in manifest["jobs"]:
        variant, index = job["variant"], job["model_index"]
        case = base / "models" / variant / f"seed_{index + 1:02d}"
        request = case / "request.json"
        if request.exists():
            payload = json.loads(request.read_text())
            if digest(Path(payload["source"])) != payload["sha256"]:
                raise ValueError(f"Published normalized model changed: {request}")
            count += 1
            continue
        original = (
            base / "predictions/canonical" / variant / f"{variant}_model_{index}.pdb"
        )
        confidence = original.with_name(f"confidence_{variant}_model_{index}.json")
        component = component_path(base, job)
        if not all(p.is_file() for p in (original, confidence, component)):
            continue
        case.mkdir(parents=True, exist_ok=True)
        normalization = normalize_model(
            original, component, Path(config["mol2"]), case / "raw"
        )
        roles = normalization["role_mapping"]
        if roles["C54"] != "C54" or roles["O49"] != "O49":
            raise ValueError("Generation contacts address a different amide")
        normalized = case / "raw/normalized.pdb"
        sequences = _pdb_sequences(normalized)
        for chain in "ABCD":
            if sequences[chain] != job["expected_alpha_sequence"]:
                raise ValueError(f"Predicted alpha sequence changed in {case}")
        for chain in "EFGH":
            if sequences[chain] != job["expected_beta_sequence"]:
                raise ValueError(f"Predicted beta sequence changed in {case}")
        write_json(
            request,
            {
                "source": str(normalized),
                "sha256": digest(normalized),
                "variant": variant,
                "seed": job["seed"],
                "model_index": index,
                "ligand_count": 4,
            },
        )
        logging.info(
            "Published graph/sequence-checked request %s seed%d", variant, index + 1
        )
        count += 1
    return count


def final_comparison(base: Path, manifest: dict, config: dict) -> None:
    """Write complete before/after metrics and figures with all planned seeds retained."""
    entries = []
    for job in manifest["jobs"]:
        case = base / "models" / job["variant"] / f"seed_{job['model_index'] + 1:02d}"
        for stage, path in [
            ("raw", case / "raw/normalized.pdb"),
            ("relaxed", case / "relaxation/relaxed.pdb"),
        ]:
            entries.append(
                {
                    "variant": job["variant"],
                    "seed": job["seed"],
                    "stage": stage,
                    "path": str(path),
                    "ligand_count": 4,
                }
            )
    result = compare_models(entries, Path(config["reference_cif"]), base / "analysis")
    plot_comparison(
        base / "analysis/occupied_sites.csv",
        base / "analysis/model_summary.csv",
        base / "analysis/plots",
        planned_models=len(manifest["jobs"]),
    )
    with (base / "analysis/relaxation_acceptance.csv").open("w", newline="") as f:
        fields = [
            "variant",
            "seed",
            "force_converged",
            "final_rms_gradient",
            "acceptance_json",
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for j in manifest["jobs"]:
            p = (
                base
                / "models"
                / j["variant"]
                / f"seed_{j['model_index'] + 1:02d}"
                / "relaxation/acceptance.json"
            )
            acceptance = json.loads(p.read_text())
            writer.writerow(
                {
                    "variant": j["variant"],
                    "seed": j["seed"],
                    "force_converged": acceptance["force_converged"],
                    "final_rms_gradient": acceptance["final_rms_gradient"],
                    "acceptance_json": str(p),
                }
            )
    write_json(
        base / "completion.json",
        {
            "process_complete": True,
            "planned_models": 45,
            "geometry_report": result,
            "interpretation": "Model seeds are replicates; four sites are nested; geometry and minimization do not validate catalytic activity.",
        },
    )


def main() -> None:
    """Keep the finite, PI-authorized three-host campaign progressing in tmux."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True, type=Path)
    args = parser.parse_args()
    base = args.campaign.resolve()
    manifest = json.loads((base / "prepared/candidate_manifest.json").read_text())
    config = json.loads((base / "configuration/relaxation.json").read_text())
    identity = config["ssh_identity"]
    if len(manifest["jobs"]) != 45:
        raise ValueError("Expected45models")
    deadline = time.monotonic() + 43200
    status_path = base / "status.json"
    write_json(
        status_path,
        {
            "status": "running",
            "planned_models": 45,
            "normalized_models": 0,
            "relaxed_models": 0,
        },
    )
    try:
        while time.monotonic() < deadline:
            sync_tree(
                f"jpc@tulip:{base}/predictions/",
                f"{base}/predictions/",
                [
                    "--include=*/",
                    "--include=*.pdb",
                    "--include=confidence_*.json",
                    "--include=*.pkl",
                    "--exclude=*",
                ],
                identity,
            )
            requested = publish_requests(base, manifest, config)
            if requested:
                sync_tree(
                    f"{base}/models/",
                    f"jpc@geranium:{base}/models/",
                    [
                        "--exclude=relaxation/",
                        "--exclude=relaxation_worker_status.json",
                        "--exclude=request.json",
                    ],
                    identity,
                )
                # Readiness is transferred only after all referenced model files.
                sync_tree(
                    f"{base}/models/",
                    f"jpc@geranium:{base}/models/",
                    ["--include=*/", "--include=request.json", "--exclude=*"],
                    identity,
                )
                sync_tree(
                    f"jpc@geranium:{base}/models/",
                    f"{base}/models/",
                    ["--exclude=raw/", "--exclude=request.json"],
                    identity,
                )
            completed = 0
            for j in manifest["jobs"]:
                case = (
                    base / "models" / j["variant"] / f"seed_{j['model_index'] + 1:02d}"
                )
                p = case / "relaxation_worker_status.json"
                if p.exists():
                    value = json.loads(p.read_text())
                    if value["status"] == "failed":
                        raise RuntimeError(f"Geranium model failed; inspect {p}")
                    if value["status"] == "completed":
                        completed += 1
            write_json(
                status_path,
                {
                    "status": "running",
                    "planned_models": 45,
                    "normalized_models": requested,
                    "relaxed_models": completed,
                },
            )
            logging.info(
                "Progress: %d/45 graph-normalized, %d/45 minimized",
                requested,
                completed,
            )
            if completed == 45:
                final_comparison(base, manifest, config)
                write_json(
                    status_path,
                    {
                        "status": "completed",
                        "planned_models": 45,
                        "normalized_models": 45,
                        "relaxed_models": 45,
                    },
                )
                return
            time.sleep(60)
        raise TimeoutError(
            "Campaign did not finish within12h; inspect retained per-model records"
        )
    except Exception as exc:
        write_json(
            status_path,
            {
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
                "planned_models": 45,
            },
        )
        raise


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
