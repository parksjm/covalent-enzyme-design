"""Verification and provenance records for reusable WT predictions."""

from __future__ import annotations
import csv
from pathlib import Path
from typing import Any
from nyl12_covalent.design_mask_validation_inputs import (
    SEEDS,
    _read_json,
    _sha,
    _yaml_identity,
)


def _existing_wt_reuse(
    previous: Path, identity: dict[str, Any], native_alpha: str, native_beta: str
) -> list[dict[str, Any]]:
    prior = _read_json(previous / "prepared/candidate_manifest.json")
    if prior.get("ligand_count") != 4 or prior.get("seeds") != SEEDS:
        raise ValueError(
            "Previous validation campaign has incompatible ligand count or seed set"
        )
    old_jobs = {
        (row["model_index"], row["seed"]): row
        for row in prior["jobs"]
        if row["variant"] == "WT"
    }
    if len(old_jobs) != 5:
        raise ValueError("Previous campaign must contain all five WT seed jobs")
    hosts = list(
        csv.DictReader(
            (previous / "dispatch/prediction_host_assignments.csv").open(newline="")
        )
    )
    host_by_key = {
        (int(row["model_index"]), int(row["seed"])): row["prediction_host"]
        for row in hosts
        if row["variant"] == "WT"
    }
    runtime = _read_json(previous / "dispatch/runtime_versions.json")
    if not runtime.get("same_boltz_source") or not runtime.get("checkpoint_sha256"):
        raise ValueError(
            "Previous runtime lacks a pinned Boltz source/checkpoint identity"
        )
    reused, new_jobs = [], []
    for index, seed in enumerate(SEEDS):
        old = old_jobs[(index, seed)]
        if (old.get("expected_alpha_sequence"), old.get("expected_beta_sequence")) != (
            native_alpha,
            native_beta,
        ):
            raise ValueError(
                "Previous WT job does not match the accepted native sequence"
            )
        if _yaml_identity(Path(old["yaml_path"])) != identity:
            raise ValueError(
                "Previous WT YAML/MSA/template scientific identity differs"
            )
        if index not in {0, 2, 4}:
            new_jobs.append({"model_index": index, "seed": seed})
            continue
        if host_by_key.get((index, seed)) != "violet":
            raise ValueError(f"WT index {index} was not predicted on Violet")
        base = previous / "models/WT" / f"seed_{index + 1:02d}"
        request = _read_json(base / "request.json")
        mapping = _read_json(base / "raw/mapping.json")
        acceptance_path = base / "relaxation/acceptance.json"
        acceptance = _read_json(acceptance_path)
        raw_path, relaxed_path = (
            base / "raw/normalized.pdb",
            base / "relaxation/relaxed.pdb",
        )
        raw_hash, relaxed_hash = _sha(raw_path), _sha(relaxed_path)
        if (
            request.get("variant"),
            request.get("seed"),
            request.get("model_index"),
            request.get("ligand_count"),
        ) != ("WT", seed, index, 4):
            raise ValueError(
                "Reusable WT model request does not match exact sequence/seed/index"
            )
        if (
            mapping.get("status") != "accepted_graph_isomorphic_component_normalization"
            or mapping.get("output_sha256") != raw_hash
            or acceptance.get("heavy_export", {}).get("output_sha256") != relaxed_hash
        ):
            raise ValueError(
                f"Reusable WT artifacts changed or lack matching hashes: {base}"
            )
        if (
            acceptance.get("status") != "minimization_complete"
            or acceptance.get("process_exit_code") != 0
        ):
            raise ValueError(f"Reusable WT Amber relaxation is incomplete: {base}")
        reused.append(
            {
                "variant": "WT",
                "model_index": index,
                "seed": seed,
                "raw_path": str(raw_path.resolve()),
                "relaxed_path": str(relaxed_path.resolve()),
                "acceptance_json": str(acceptance_path.resolve()),
                "source_sha256": {
                    "raw": raw_hash,
                    "relaxed": relaxed_hash,
                    "mapping": _sha(base / "raw/mapping.json"),
                    "acceptance": _sha(acceptance_path),
                    "yaml": _sha(Path(old["yaml_path"])),
                },
            }
        )
    if len(reused) != 3 or [row["model_index"] for row in reused] != [0, 2, 4]:
        raise ValueError("Exactly Violet WT model indices 0, 2 and 4 must be reused")
    if sorted((row["model_index"], row["seed"]) for row in new_jobs) != [
        (1, SEEDS[1]),
        (3, SEEDS[3]),
    ]:
        raise ValueError("WT indices 1 and 3 must be predicted anew")
    return reused
