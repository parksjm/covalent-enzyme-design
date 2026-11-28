#!/usr/bin/env python3
"""Execute the qualified, bounded broad LigandMPNN pilot and matched scores."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.ligandmpnn_pilot_audit import _fasta, _pdb_sequences  # noqa: E402
from nyl12_covalent.relay_broad_pilot import sample_schedule  # noqa: E402
from nyl12_covalent.relay_broad_release import (  # noqa: E402
    FOUR_TEMPLATE_PROTOCOL,
    four_template_schedule,
)
from nyl12_covalent.relay_broad_runtime import (  # noqa: E402
    digest,
    export_local,
    journal,
    native_run,
    parser_audit,
    rank_broad,
    score_summary,
    verify_parent,
    write_new,
)
from nyl12_covalent.relay_ensemble_sequences import audit_sequence, pool_samples  # noqa: E402


def flags(manifest: dict[str, Any], pdb: str, out: Path) -> list[str]:
    """Return the frozen native atomic-context flags."""
    generation = manifest["request"]["generation"]
    return [
        "--model_type",
        "ligand_mpnn",
        "--checkpoint_ligand_mpnn",
        manifest["native"]["checkpoint"],
        "--pdb_path",
        pdb,
        "--out_folder",
        str(out),
        "--redesigned_residues",
        generation["redesigned_residues"],
        "--symmetry_residues",
        generation["symmetry_residues"],
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


def selected_remodels(manifests: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Select one accepted native output per original-parent/noise pair."""
    selected = {}
    for parent in sorted(
        [p for m in manifests for p in m["parents"]], key=lambda p: p["parent_id"]
    ):
        verify_parent(parent)
        key = (parent["parent_source_id"], parent["noise_A"])
        if key[0] not in {"medoid01", "medoid03"} or key[1] not in {2, 5}:
            raise ValueError("Unexpected original/noise parent class")
        selected.setdefault(key, parent)
    return list(selected.values())


def prepare_templates(manifest: dict[str, Any], output: Path) -> list[dict[str, Any]]:
    """Export native inputs and require actual native parser acceptance."""
    original = manifest["request"]["original_parent_ids"]
    selected = selected_remodels(manifest["accepted_manifests"])
    remodeled = [p["parent_id"] for p in selected]
    if "allocation_protocol" in manifest:
        if manifest["allocation_protocol"] != FOUR_TEMPLATE_PROTOCOL:
            raise ValueError("Unknown explicit broad allocation protocol")
        schedule = four_template_schedule(original, remodeled)
    else:
        schedule = sample_schedule(original, remodeled)
    parents = [
        {
            "parent_id": name,
            "pdb_path": manifest["request"]["original_parent_pdbs"][name],
        }
        for name in original
    ]
    parents.extend(selected)
    templates = []
    for parent, allocation in zip(parents, schedule, strict=True):
        job = output / parent["parent_id"]
        job.mkdir()
        pdb = job / "ligandmpnn_input.pdb"
        export_local(Path(parent["pdb_path"]), pdb)
        audit = parser_audit(pdb, Path(manifest["native"]["generation_engine"]))
        write_new(job / "native_parser_acceptance.json", audit)
        templates.append(
            {"id": parent["parent_id"], "pdb": str(pdb), "state": "TS", **allocation}
        )
    es = manifest["native"]["ES_template"]
    audit = parser_audit(Path(es["pdb"]), Path(manifest["native"]["generation_engine"]))
    write_new(output / "ES_native_parser_acceptance.json", audit)
    templates.append({**es, "state": "ES"})
    return templates


def run(manifest: dict[str, Any], output: Path) -> None:
    """Generate the fixed raw budget, audit whole sequences and crossscore."""
    for path, expected in manifest["source_sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Frozen source hash mismatch: {path}")
    request = manifest["request"]
    generation = request["generation"]
    if (
        generation["temperature"],
        generation["native_identity_bias"],
        generation["max_raw_samples"],
    ) != (0.1, 0, 24):
        raise ValueError("Unauthorized broad generation settings")
    output.mkdir(exist_ok=False)
    write_new(output / "manifest.json", manifest)
    journal(ROOT, "intent", "Prepare broad parent exports", str(output))
    templates = prepare_templates(manifest, output)
    native = manifest["native"]
    env = native["environment"]
    reference = _pdb_sequences(Path(templates[0]["pdb"]))
    mask = generation["mask"]
    for template in templates:
        if _pdb_sequences(Path(template["pdb"])) != reference:
            raise ValueError("Crossscore template sequences differ from WT")
    journal(
        ROOT, "outcome", "Prepare broad parent exports", "All native parsers accepted"
    )
    code = native_run(
        ROOT,
        [
            native["python"],
            "-c",
            "import torch; assert torch.cuda.is_available(), 'GPU required'",
        ],
        output / "GPU_prerequisite",
        env,
    )
    if code:
        raise RuntimeError("Authorized GPU runtime is unavailable")
    samples = []
    for index, template in enumerate(templates[:-1]):
        count = template["raw_samples"]
        job = output / f"generation_{index + 1:02d}"
        argv = [
            native["python"],
            native["generation_engine"],
            *flags(manifest, template["pdb"], job),
            "--symmetry_weights",
            generation["symmetry_weights"],
            "--temperature",
            "0.1",
            "--seed",
            str(2026112620 + index),
            "--number_of_batches",
            str(count),
            "--save_stats",
            "1",
            "--pack_side_chains",
            "0",
        ]
        code = native_run(ROOT, argv, output / f"generate_{index + 1:02d}", env)
        if code:
            raise RuntimeError(f"Native generation failed: exit={code}")
        records = _fasta(job / "seqs/ligandmpnn_input.fa")
        if len(records) != count + 1 or records[0][1] != ":".join(reference.values()):
            raise ValueError("Unexpected native identity or raw sample count")
        for number, (header, sequence) in enumerate(records[1:], 1):
            samples.append(
                {
                    **audit_sequence(sequence, reference, mask),
                    "sample_id": f"{template['id']}_{number}",
                    "parent": template["id"],
                    "native_header": header,
                }
            )
    candidates = pool_samples(samples)
    write_new(output / "raw_samples.json", samples)
    write_new(output / "deduplicated_pool.json", candidates)
    if not candidates:
        raise RuntimeError("Native generation returned an empty candidate pool")
    baseline = {
        **audit_sequence(":".join(reference.values()), reference, mask),
        "candidate_id": "WT",
        "reference": True,
    }
    import torch

    orders = {}
    for candidate in [baseline, *candidates]:
        candidate["scores"] = []
        for index, template in enumerate(templates):
            tag = f"{candidate['candidate_id']}_{index + 1:02d}"
            target = output / f"scores_{tag}"
            score_request = {
                "sequence": candidate["sequence"],
                "mask": mask,
                "score_engine": native["score_engine"],
                "pdb": template["pdb"],
                "argv": [
                    *flags(manifest, template["pdb"], target),
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
            write_new(path, score_request)
            code = native_run(
                ROOT,
                [native["python"], native["score_wrapper"], str(path)],
                output / f"score_{tag}",
                env,
            )
            if code:
                raise RuntimeError(f"Native crossscore failed: exit={code}")
            stats = torch.load(
                target / (Path(template["pdb"]).stem + ".pt"),
                map_location="cpu",
                weights_only=False,
            )
            order = np.asarray(stats["decoding_order"])
            if order.shape != (2, 1312):
                raise ValueError("Crossscore requires two full decoding orders")
            if template["id"] in orders and not np.array_equal(
                order, orders[template["id"]]
            ):
                raise ValueError("WT/design decoding orders differ")
            orders[template["id"]] = order.copy()
            candidate["scores"].append(
                {
                    "template": template["id"],
                    "state": template["state"],
                    **score_summary(
                        stats,
                        candidate["sequence"],
                        mask,
                        manifest["local_pocket_mask"],
                    ),
                }
            )
    ranked = rank_broad(candidates)
    write_new(output / "WT_reference_scores.json", baseline)
    write_new(output / "ranked_candidates.json", ranked)
    write_new(
        output / "preselected_candidate.json",
        {
            **ranked[0],
            "selection_basis": "Frozen local-pocket TS mean NLL, then worst NLL, mutation count, candidate ID; diagnostic only",
        },
    )
    write_new(
        output / "scientific_acceptance.json",
        {
            "accepted": True,
            "raw_count": len(samples),
            "unique_count": len(candidates),
            "ties_protections_and_matched_crossscores_valid": True,
            "activity_validated": False,
            "prediction_candidate_id": ranked[0]["candidate_id"],
        },
    )


def main() -> None:
    """Run one immutable broad attempt from its concrete frozen manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(json.loads(args.manifest.read_text()), args.output)


if __name__ == "__main__":
    main()
