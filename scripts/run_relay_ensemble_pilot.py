#!/usr/bin/env python3
"""Execute a frozen 3-by-8 exploratory LigandMPNN pilot and crossscore its pool."""

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

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.ligandmpnn_pilot_audit import (  # noqa: E402
    _fasta,
    _pdb_sequences,
)
from nyl12_covalent.relay_ensemble_sequences import (  # noqa: E402
    audit_sequence,
    pool_samples,
    rank_candidates,
    summarize_score,
)

LOGGER = logging.getLogger(__name__)
ENVIRONMENT = {
    "CUDA_VISIBLE_DEVICES": "0",
    "OMP_NUM_THREADS": "2",
    "OPENBLAS_NUM_THREADS": "2",
    "MKL_NUM_THREADS": "2",
}


def write_json(path: Path, data: object) -> None:
    """Write a new JSON artifact, refusing replacement."""
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def native_run(argv: list[str], output: Path, env: dict[str, str]) -> None:
    """Retain an exact native command, log and independent process exit record."""
    write_json(output.with_suffix(".command.json"), argv)
    with output.with_suffix(".log").open("x") as handle:
        result = subprocess.run(
            argv, env=env, stdout=handle, stderr=subprocess.STDOUT, check=False
        )
    write_json(output.with_suffix(".process.json"), {"exit_code": result.returncode})
    result.check_returncode()


def common_flags(manifest: dict, template: dict, output: Path) -> list[str]:
    """Return explicit shared native flags for this unchanged atomic parent."""
    return [
        "--model_type",
        "ligand_mpnn",
        "--checkpoint_ligand_mpnn",
        manifest["checkpoint"],
        "--pdb_path",
        template["pdb"],
        "--out_folder",
        str(output),
        "--redesigned_residues",
        manifest["redesigned_residues"],
        "--symmetry_residues",
        manifest["symmetry_residues"],
        "--homo_oligomer",
        "0",
        "--ligand_mpnn_use_atom_context",
        "1",
        "--ligand_mpnn_use_side_chain_context",
        "1",
        "--parse_atoms_with_zero_occupancy",
        "1",
        "--batch_size",
        "1",
    ]


def validate(manifest: dict) -> None:
    """Check the immutable budget, approved mask, environment and content hashes."""
    if manifest["approved_for_launch"] is not True:
        raise ValueError("Final canonical mask has not been approved for launch")
    if len(manifest["templates"]) != 3 or manifest["samples_per_template"] != 8:
        raise ValueError("Only three templates with eight samples each are authorized")
    if manifest["environment"] != {**ENVIRONMENT, "PYTHONPATH": manifest["support"]}:
        raise ValueError("Unexpected GPU0/two-thread runtime environment")
    for path, expected in manifest["source_sha256"].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source/input hash mismatch: {path}")
    required = {
        manifest["checkpoint"],
        manifest["generation_engine"],
        manifest["score_engine"],
        manifest["score_wrapper"],
        str(Path(__file__)),
    }
    required.update(t["pdb"] for t in manifest["templates"])
    if not required <= manifest["source_sha256"].keys():
        raise ValueError("Manifest omits required source/input hashes")
    groups = manifest["symmetry_residues"].split("|")
    if manifest["symmetry_weights"].split("|") != ["1,0,0,0"] * len(groups):
        raise ValueError("Generation weights must select only active A/E copies")
    expected = [
        ",".join(f"{c}{p}" for c in chains)
        for family, chains in (("alpha", "ABCD"), ("beta", "EFGH"))
        for p in manifest["mask"][family]
    ]
    if groups != expected:
        raise ValueError("Symmetry groups differ from approved canonical mask")
    if set(manifest["redesigned_residues"].split()) != {
        residue for group in expected for residue in group.split(",")
    }:
        raise ValueError("Designed positions differ from tied mask")
    if not expected:
        raise ValueError("Empty mutation mask")
    native = _pdb_sequences(Path(manifest["templates"][0]["pdb"]))
    parents = [*manifest["templates"]]
    if manifest.get("ES_template") is not None:
        if manifest["ES_control"]["status"] != "accepted_exact_mapping":
            raise ValueError("ES control is not accepted with exact mapping")
        parents.append(manifest["ES_template"])
        if manifest["ES_template"]["pdb"] not in manifest["source_sha256"]:
            raise ValueError("ES control lacks input hash")
    for template in parents:
        if _pdb_sequences(Path(template["pdb"])) != native:
            raise ValueError("Templates have different canonical WT sequences")
        audit_sequence(":".join(native.values()), native, manifest["mask"])
    if not manifest["parser_audit"]["accepted"]:
        raise ValueError("Actual native parsing/catalytic-context audit failed")


def run(
    manifest: dict,
    output: Path,
    generation_only: bool = False,
    generation_source: Path | None = None,
) -> list[dict]:
    """Generate at most 24 unique candidates, then score each complete sequence."""
    validate(manifest)
    output.mkdir(exist_ok=False)
    write_json(output / "manifest.json", manifest)
    env = {**os.environ, **manifest["environment"]}
    # No automatic CPU fallback when the scientific runtime loses GPU access.
    native_run(
        [
            manifest["python"],
            "-c",
            "import torch; assert torch.cuda.is_available(), 'GPU required'",
        ],
        output / "gpu_prerequisite",
        env,
    )
    samples = []
    if generation_source is not None:
        samples = json.loads((generation_source / "raw_samples.json").read_text())
        if len(samples) != 24:
            raise ValueError("Retained generation must contain exactly24raw samples")
        native = _pdb_sequences(Path(manifest["templates"][0]["pdb"]))
        for sample in samples:
            audit_sequence(sample["sequence"], native, manifest["mask"])
    templates = [] if generation_source is not None else manifest["templates"]
    for index, template in enumerate(templates):
        job = output / f"generation_{index + 1:02d}"
        argv = [
            manifest["python"],
            manifest["generation_engine"],
            *common_flags(manifest, template, job),
            "--symmetry_weights",
            manifest["symmetry_weights"],
            "--temperature",
            "0.1",
            "--seed",
            str(2026112601 + index),
            "--number_of_batches",
            "8",
            "--save_stats",
            "1",
            "--pack_side_chains",
            "0",
        ]
        native_run(argv, output / f"generation_{index + 1:02d}_native", env)
        records = _fasta(job / "seqs" / (Path(template["pdb"]).stem + ".fa"))
        native = _pdb_sequences(Path(template["pdb"]))
        if len(records) != 9 or records[0][1] != ":".join(native.values()):
            raise ValueError("Native record or eight-sample count mismatch")
        for sample_index, (header, sequence) in enumerate(records[1:], 1):
            sample = audit_sequence(sequence, native, manifest["mask"])
            samples.append(
                {
                    **sample,
                    "sample_id": f"TS{index + 1}_{sample_index}",
                    "native_header": header,
                    "parent": template["id"],
                }
            )
    candidates = pool_samples(samples)
    write_json(output / "raw_samples.json", samples)
    write_json(output / "deduplicated_pool.json", candidates)
    if generation_only:
        write_json(
            output / "generation_scientific_acceptance.json",
            {
                "accepted": bool(candidates),
                "raw_count": len(samples),
                "unique_count": len(candidates),
                "ties_and_protections_valid": True,
                "crossscoring": "pending",
                "empty_pool": not candidates,
            },
        )
        return candidates
    if not candidates:
        write_json(
            output / "scientific_acceptance.json",
            {"accepted": False, "outcome": "empty_pool", "raw_count": len(samples)},
        )
        return []
    import torch

    native = _pdb_sequences(Path(manifest["templates"][0]["pdb"]))
    baseline = {
        **audit_sequence(":".join(native.values()), native, manifest["mask"]),
        "candidate_id": "WT",
        "reference": True,
        "lineage": [],
    }
    orders = {}
    for candidate in [baseline, *candidates]:
        candidate["scores"] = []
        score_parents = [*manifest["templates"]]
        if manifest.get("ES_template") is not None:
            score_parents.append(manifest["ES_template"])
        for template in score_parents:
            tag = f"{candidate['candidate_id']}_{template['id']}"
            native_output = output / f"scores_{tag}"
            request = {
                "sequence": candidate["sequence"],
                "mask": manifest["mask"],
                "score_engine": manifest["score_engine"],
                "pdb": template["pdb"],
                "argv": [
                    *common_flags(manifest, template, native_output),
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
            write_json(path, request)
            native_run(
                [manifest["python"], manifest["score_wrapper"], str(path)],
                output / f"score_{tag}_native",
                env,
            )
            stats = torch.load(
                native_output / (Path(template["pdb"]).stem + ".pt"),
                map_location="cpu",
                weights_only=False,
            )
            order = np.asarray(stats["decoding_order"])
            if order.shape != (2, 1312):
                raise ValueError("Crossscore must retain two complete decoding orders")
            if template["id"] not in orders:
                orders[template["id"]] = order.copy()
            elif not np.array_equal(order, orders[template["id"]]):
                raise ValueError("Candidate and WT decoding orders are unmatched")
            candidate["scores"].append(
                {
                    "template": template["id"],
                    "state": "ES" if template["id"] == "ES" else "TS",
                    **summarize_score(stats, candidate["sequence"], manifest["mask"]),
                }
            )
    write_json(output / "WT_reference_scores.json", baseline)
    ranked = rank_candidates(candidates)
    write_json(output / "ranked_candidates.json", ranked)
    with (output / "ranking.csv").open("x", newline="") as handle:
        keys = [
            "candidate_id",
            "rank",
            "mutation_count",
            "ts_mean_nll",
            "ts_worst_nll",
            "pareto",
        ]
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows({k: c[k] for k in keys} for c in ranked)
    write_json(
        output / "scientific_acceptance.json",
        {
            "accepted": True,
            "raw_count": len(samples),
            "unique_count": len(ranked),
            "sequence_protection_and_ties_valid": True,
            "finite_crossscores": True,
            "ES_control": manifest["ES_control"],
            "catalytic_improvement_validated": False,
        },
    )
    LOGGER.info(
        "Audited %s raw samples, %s unique complete sequences",
        len(samples),
        len(ranked),
    )
    return ranked


def main() -> None:
    """Read the approved frozen manifest and create one retained run attempt."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--generation-only", action="store_true")
    parser.add_argument("--score-existing", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(
        json.loads(args.manifest.read_text()),
        args.output,
        args.generation_only,
        args.score_existing,
    )


if __name__ == "__main__":
    main()
