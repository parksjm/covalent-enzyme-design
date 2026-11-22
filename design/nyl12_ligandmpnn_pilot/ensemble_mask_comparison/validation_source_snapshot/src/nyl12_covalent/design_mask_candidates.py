"""Prepare diverse audited mask designs for matched Boltz validation."""

from __future__ import annotations
import csv
import json
import logging
import shutil
from pathlib import Path
from typing import Any
from nyl12_covalent.design_boltz_inputs import (
    _native_sequences,
    _prepare_yaml,
    _write_msa,
)
from nyl12_covalent.design_mask_selection import select_candidates
from nyl12_covalent.design_mask_validation_inputs import (
    ARMS,
    SEEDS,
    _manifest_masks,
    _msa_hash,
    _read_json,
    _sha,
    _yaml_identity,
    _yaml_identity_data,
)
from nyl12_covalent.design_mask_validation_records import (
    _audit_samples_and_build_records,
)
from nyl12_covalent.design_mask_validation_reuse import _existing_wt_reuse

LOGGER = logging.getLogger(__name__)


def prepare_validation(
    report_path: Path, previous_campaign: Path, output: Path, seeds: list[int]
) -> dict[str, Any]:
    """Prepare a bounded four-ligand Boltz cohort and reusable WT controls.

    Args:
        report_path: Completed 48-sample design-mask comparison report.
        previous_campaign: Completed 45-model ``boltz_validation`` directory.
        output: New validation root for prepared inputs and manifests.
        seeds: Fixed five Boltz prediction seeds authorized by the protocol.

    Returns:
        Candidate/job manifest, including three exact Violet WT model reuses.

    Raises:
        FileNotFoundError: A declared input/source is missing.
        FileExistsError: The output root already exists.
        ValueError: Design, mask, protection, sequence, identity, or reuse gates fail.
    """
    if seeds != SEEDS:
        raise ValueError("Prediction seeds must match the fixed five-seed protocol")
    if output.exists():
        raise FileExistsError(output)
    report = _read_json(report_path)
    samples = report.get("samples")
    if not isinstance(samples, list) or len(samples) != 48:
        raise ValueError("Preparation requires all 48 audited raw samples")
    acceptance_path = report_path.parent / "acceptance.json"
    acceptance = _read_json(acceptance_path)
    if acceptance.get("accepted") is not True or acceptance.get(
        "report_sha256"
    ) != _sha(report_path):
        raise ValueError(
            "Generation report lacks an accepted, hash-matched acceptance record"
        )
    manifest_path = Path(report["manifest"])
    source_hashes = report.get("source_sha256", {})
    if not source_hashes:
        raise ValueError("Generation report has no source hashes")
    for raw_path, expected_hash in source_hashes.items():
        path = Path(raw_path)
        if not path.is_file() or _sha(path) != expected_hash:
            raise ValueError(f"Generation source hash mismatch: {path}")
    masks, protected, generation, parent_paths = _manifest_masks(manifest_path)
    selected = select_candidates(samples)
    prior = _read_json(previous_campaign / "prepared/candidate_manifest.json")
    native_pdb = Path(prior["native_pdb"])
    native = _native_sequences(native_pdb)
    native_alpha, native_beta = native["A"], native["E"]
    for parent_path in parent_paths.values():
        parent_native = _native_sequences(parent_path)
        if (parent_native["A"], parent_native["E"]) != (native_alpha, native_beta):
            raise ValueError(
                f"Accepted WT parent sequence differs from previous native: {parent_path}"
            )
    template_path = previous_campaign / "configuration/template.json"
    template = _read_json(template_path)
    proteins = [
        entry["protein"] for entry in template["sequences"] if "protein" in entry
    ]
    if len(proteins) != 2 or (proteins[0]["sequence"], proteins[1]["sequence"]) != (
        native_alpha,
        native_beta,
    ):
        raise ValueError("Template sequences do not match the previous native PDB")
    if len(template.get("constraints", [])) != 2 or not template.get("templates"):
        raise ValueError(
            "Validated template must preserve both PA6 contacts and N12C template"
        )
    alpha_source, beta_source = Path(proteins[0]["msa"]), Path(proteins[1]["msa"])
    source_strata = {
        (
            job["arm"],
            str(job["parent"]),
            float(job["argv"][job["argv"].index("--temperature") + 1]),
        )
        for job in generation["jobs"]
    }
    records, union = _audit_samples_and_build_records(
        samples,
        selected,
        masks,
        protected,
        native_alpha,
        native_beta,
        source_strata,
    )

    input_root = output / "prepared"
    msa_dir = input_root / "msas"
    for record in records:
        cid = record["id"]
        amsa, bmsa = msa_dir / cid / "alpha.csv", msa_dir / cid / "beta.csv"
        if cid == "WT":
            data = _prepare_yaml(template, native_alpha, native_beta, amsa, bmsa, 4)
            identity = _yaml_identity_data(
                data,
                _msa_hash(alpha_source, native_alpha, native_alpha),
                _msa_hash(beta_source, native_beta, native_beta),
            )
    old_jobs = [job for job in prior["jobs"] if job["variant"] == "WT"]
    if len(old_jobs) != 5:
        raise ValueError("Previous candidate manifest does not contain five WT jobs")
    if any(_yaml_identity(Path(job["yaml_path"])) != identity for job in old_jobs):
        raise ValueError("Previous WT YAML/MSA/template scientific identity differs")
    reused = _existing_wt_reuse(previous_campaign, identity, native_alpha, native_beta)

    relaxation_source = previous_campaign / "configuration/relaxation.json"
    relaxation = _read_json(relaxation_source)
    config_sources = {
        "whole_l4_source.mol2": Path(relaxation["mol2"]),
        "whole_l4_context_mapping.json": Path(relaxation["mapping"]),
    }
    for path in (
        Path(relaxation["release"]),
        *config_sources.values(),
        Path(relaxation["reference_cif"]),
    ):
        if not path.exists():
            raise FileNotFoundError(path)
    for entry in template.get("templates", []):
        if _sha(Path(entry["cif"])) != _sha(Path(relaxation["reference_cif"])):
            raise ValueError("Boltz template CIF differs from Amber reference CIF")

    output.mkdir(parents=True)
    (output / "prepared/inputs").mkdir(parents=True)
    (output / "prepared/msas").mkdir(parents=True)
    (output / "configuration").mkdir()
    shutil.copy2(relaxation_source, output / "configuration/relaxation.json")
    for name, source in config_sources.items():
        shutil.copy2(source, output / "configuration" / name)
    jobs: list[dict[str, Any]] = []
    reused_keys = {(row["model_index"], row["seed"]) for row in reused}
    for variant_index, record in enumerate(records):
        cid = record["id"]
        alpha_msa, beta_msa = msa_dir / cid / "alpha.csv", msa_dir / cid / "beta.csv"
        _write_msa(alpha_source, alpha_msa, native_alpha, record["alpha_sequence"])
        _write_msa(beta_source, beta_msa, native_beta, record["beta_sequence"])
        data = _prepare_yaml(
            template,
            record["alpha_sequence"],
            record["beta_sequence"],
            alpha_msa,
            beta_msa,
            4,
        )
        yaml_path = input_root / "inputs" / f"{cid}.yaml"
        yaml_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        for model_index, seed in enumerate(seeds):
            if cid == "WT" and (model_index, seed) in reused_keys:
                continue
            jobs.append(
                {
                    "campaign": "nyl12_model23_mask_validation",
                    "sequence_name": cid,
                    "variant": cid,
                    "model_index": model_index,
                    "seed": seed,
                    "yaml_path": str(yaml_path.resolve()),
                    "expected_alpha_sequence": record["alpha_sequence"],
                    "expected_beta_sequence": record["beta_sequence"],
                }
            )
    with (input_root / "jobs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "campaign",
                "sequence_name",
                "model_index",
                "seed",
                "yaml_path",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows({key: job[key] for key in writer.fieldnames} for job in jobs)
    hashes = {
        str(path.resolve()): _sha(path)
        for path in (
            report_path,
            manifest_path,
            native_pdb,
            template_path,
            alpha_source,
            beta_source,
            relaxation_source,
            *config_sources.values(),
            Path(relaxation["reference_cif"]),
        )
    }
    hashes.update({row["raw_path"]: row["source_sha256"]["raw"] for row in reused})
    hashes.update(
        {row["relaxed_path"]: row["source_sha256"]["relaxed"] for row in reused}
    )
    result = {
        "schema_version": 1,
        "report": str(report_path.resolve()),
        "report_sha256": _sha(report_path),
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": _sha(manifest_path),
        "native_pdb": str(native_pdb.resolve()),
        "template_yaml_json_source": str(template_path.resolve()),
        "ligand_count": 4,
        "seeds": seeds,
        "candidates": records,
        "jobs": jobs,
        "reused_models": reused,
        "source_sha256": hashes,
        "source_masks": {
            arm: {family: sorted(pos) for family, pos in masks[arm].items()}
            for arm in ARMS
        },
        "source_mask_union": {family: sorted(pos) for family, pos in union.items()},
        "hard_protected_positions": {
            family: sorted(pos) for family, pos in protected.items()
        },
        "raw_sample_count": len(samples),
        "runtime_requirement": {
            "host": "violet",
            "source_runtime": _read_json(
                previous_campaign / "dispatch/runtime_versions.json"
            )["violet"],
            "checkpoint_sha256": _read_json(
                previous_campaign / "dispatch/runtime_versions.json"
            )["checkpoint_sha256"],
            "same_boltz_source": True,
        },
        "new_prediction_job_count": len(jobs),
        "calculation_launched": False,
    }
    manifest_out = input_root / "candidate_manifest.json"
    manifest_out.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for job in jobs:
        job["yaml_path"] = str(Path(job["yaml_path"]).resolve())
    LOGGER.info(
        "Prepared %d candidates and %d new prediction jobs", len(records), len(jobs)
    )
    return result
