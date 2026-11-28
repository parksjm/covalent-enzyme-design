#!/usr/bin/env python3
"""Resume matched scores from immutable generated sequences; never regenerate."""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
from pathlib import Path
import sys
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(HERE / "scoring/src"))

from nyl12_covalent.ligandmpnn_pilot_audit import _pdb_sequences  # noqa: E402
from nyl12_covalent.relay_broad_runtime import (  # noqa: E402
    digest,
    native_run,
    rank_broad,
    score_summary,
    write_new,
)
from nyl12_covalent.relay_ensemble_sequences import audit_sequence  # noqa: E402


def verify_inputs(config: dict[str, Any]) -> None:
    """Require all bound input/source bytes before any score process starts."""
    for path, expected in config["sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Score continuation input/source changed: {path}")
    if config["raw_count"] != 24 or config["unique_count"] != 24:
        raise ValueError("Score-only recovery requires this exact 24-sequence pool")


def load_runner(path: Path) -> Any:
    """Load the retained broad runner solely for its unchanged flag builder."""
    spec = importlib.util.spec_from_file_location("retained_broad_flags", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot import the retained scoring flag builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def describe_stats(
    path: Path, sequence: str, mask: dict[str, Any], local: dict[str, Any]
) -> tuple[dict[str, float], np.ndarray]:
    """Validate native stats and return existing score summaries and orders."""
    import torch

    stats = torch.load(path, map_location="cpu", weights_only=False)
    order = np.asarray(stats["decoding_order"])
    if order.shape != (2, 1312):
        raise ValueError("Expected the same two full-length native decoding orders")
    return score_summary(stats, sequence, mask, local), order


def run(config: dict[str, Any], output: Path) -> None:
    """Cross-score the existing pool with matched retained WT decoding orders.

    Args:
        config: Exact source, pool, template and retained-WT bindings.
        output: New output directory; previous attempts cannot be overwritten.
    """
    verify_inputs(config)
    manifest = json.loads(Path(config["generation_manifest"]).read_text())
    original = Path(config["generation_output"])
    candidates = json.loads((original / "deduplicated_pool.json").read_text())
    templates = config["templates"]
    if len(candidates) != 24 or len(templates) != 5:
        raise ValueError("Expected 24 existing designs and four TS plus one ES")
    if [item["state"] for item in templates] != ["TS"] * 4 + ["ES"]:
        raise ValueError("Template state/order mismatch")
    native = manifest["native"]
    native["score_wrapper"] = config["coherent_score_wrapper"]
    mask = manifest["request"]["generation"]["mask"]
    local = manifest["local_pocket_mask"]
    reference = _pdb_sequences(Path(templates[0]["pdb"]))
    ids = [item["candidate_id"] for item in candidates]
    if len(set(ids)) != 24 or not all(i.startswith("relay_corr01_") for i in ids):
        raise ValueError("Unexpected or duplicate corrected candidate IDs")
    for candidate in candidates:
        checked = audit_sequence(candidate["sequence"], reference, mask)
        if checked["mutation_count"] != candidate["mutation_count"]:
            raise ValueError("Retained candidate mutation audit changed")
    output.mkdir(exist_ok=False)
    write_new(output / "resume_manifest.json", config)
    baseline = {
        **audit_sequence(":".join(reference.values()), reference, mask),
        "candidate_id": "WT",
        "reference": True,
        "scores": [],
    }
    orders = {}
    for index, template in enumerate(templates, 1):
        receipt = json.loads(
            (original / f"score_WT_{index:02d}.process.json").read_text()
        )
        if receipt["exit_code"] != 0:
            raise ValueError("Cannot reuse a failed WT score")
        scores, order = describe_stats(
            Path(template["WT_stats"]), baseline["sequence"], mask, local
        )
        orders[template["id"]] = order
        baseline["scores"].append(
            {"template": template["id"], "state": template["state"], **scores}
        )
    write_new(output / "WT_reference_scores.json", baseline)
    runner = load_runner(Path(config["broad_runner"]))
    for candidate in candidates:
        candidate["scores"] = []
        for index, template in enumerate(templates, 1):
            tag = f"{candidate['candidate_id']}_{index:02d}"
            target = output / f"scores_{tag}"
            request = {
                "sequence": candidate["sequence"],
                "mask": mask,
                "score_engine": native["score_engine"],
                "pdb": template["pdb"],
                "argv": [
                    *runner.flags(manifest, template["pdb"], target),
                    "--seed",
                    "2026112699",
                    "--number_of_batches",
                    "2",
                    "--use_sequence",
                    "1",
                    "--autoregressive_score",
                    "1",
                    "--single_aa_score",
                    "0",
                    "--verbose",
                    "0",
                ],
            }
            path = output / f"request_{tag}.json"
            write_new(path, request)
            code = native_run(
                ROOT,
                [native["python"], native["score_wrapper"], str(path)],
                output / f"score_{tag}",
                native["environment"],
            )
            if code:
                raise RuntimeError(
                    f"Score-only native evaluation failed: {tag}, {code}"
                )
            scores, order = describe_stats(
                target / (Path(template["pdb"]).stem + ".pt"),
                candidate["sequence"],
                mask,
                local,
            )
            if not np.array_equal(order, orders[template["id"]]):
                raise ValueError("WT/design decoding orders differ")
            candidate["scores"].append(
                {"template": template["id"], "state": template["state"], **scores}
            )
        write_new(output / f"completed_{candidate['candidate_id']}.json", candidate)
        logging.info("Completed matched scores for %s", candidate["candidate_id"])
    ranked = rank_broad(candidates)
    write_new(output / "ranked_candidates.json", ranked)
    write_new(output / "preselected_candidate.json", ranked[0])
    write_new(
        output / "scientific_acceptance.json",
        {
            "accepted": True,
            "raw_count": 24,
            "unique_count": 24,
            "retained_WT_scores_reused": 5,
            "native_candidate_scores": 120,
            "ties_protections_and_matched_crossscores_valid": True,
            "activity_validated": False,
            "generation_was_not_repeated": True,
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(json.loads(args.config.read_text()), args.output)
