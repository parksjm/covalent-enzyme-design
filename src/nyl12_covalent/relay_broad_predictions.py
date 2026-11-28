"""Concrete eight-model unguided AF3/Boltz comparison preparation and execution."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from .relay_broad_graph import ROOT, COMPONENT, SMILES, digest, prepare_graph
from .relay_broad_predictor_runtime import (
    AF3,
    BOLTZ_WRAPPER,
    commands,
    run_comparison as run_comparison,
)

LOGGER = logging.getLogger(__name__)
SEEDS = (20261126, 20261127)
DEFAULT_PROTECTION = (
    ROOT
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/canonical_protection.json"
)
DEFAULT_ATOMS = (
    ROOT
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/prediction_preparation/boltz_actual_ligand_atoms.json"
)


def _write(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def _sequences(data: dict) -> tuple[str, str]:
    alpha, beta = data["alpha_sequence"], data["beta_sequence"]
    if (
        len(alpha) != 232
        or len(beta) != 96
        or set(alpha + beta) - set("ACDEFGHIKLMNPQRSTVWY")
    ):
        raise ValueError("Expectedcanonical232alphaand96betaproteinsequences")
    if not beta.startswith("T"):
        raise ValueError("MaturebetachainmuststartwithThr233")
    return alpha, beta


def prepare_comparison(
    campaign: Path,
    design_json: Path,
    protection: Path = DEFAULT_PROTECTION,
    native_atoms: Path = DEFAULT_ATOMS,
) -> dict:
    """Prepare exactly WT/design × two seeds × two unguided architectures.

    Args:
        campaign: Dedicated new prediction directory; existing manifest is verified.
        design_json: Single preselected broad design with alpha/beta sequences.
        protection: Canonical protection manifest containing authoritative WT sequences.
        native_atoms: Actual installed Boltz parser's ligand atom table.

    Returns:
        Eight-job manifest with query-only MSAs, fresh E-S graph and expected CIFs.
    """
    inputs = {
        str(p.resolve()): digest(p) for p in [design_json, protection, native_atoms]
    }
    manifest_path = campaign / "prediction_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest["input_source_sha256"] != inputs:
            raise ValueError("Preparedpredictioninputsourceschanged")
        for filename, expected in manifest["prepared_sha256"].items():
            if digest(Path(filename)) != expected:
                raise ValueError("Preparedpredictorinputchanged:" + filename)
        return manifest
    if campaign.exists():
        raise FileExistsError(campaign)
    native = json.loads(protection.read_text())["source_sequences"]
    wt = dict(alpha_sequence=native["alpha"], beta_sequence=native["beta"])
    design = json.loads(design_json.read_text())
    candidates = {"WT": _sequences(wt), "DESIGN": _sequences(design)}
    if candidates["WT"] == candidates["DESIGN"]:
        raise ValueError("PreselecteddesignmustdifferfromWT")
    campaign.mkdir(parents=True)
    graph = prepare_graph(campaign, native_atoms)
    jobs = []
    for role, (alpha, beta) in candidates.items():
        msa = campaign / "query_only_msas" / role
        msa.mkdir(parents=True)
        alpha_msa, beta_msa = msa / "alpha.a3m", msa / "beta.a3m"
        alpha_msa.write_text(f">query\n{alpha}\n")
        beta_msa.write_text(f">query\n{beta}\n")
        for predictor in ("af3", "boltz"):
            for seed in SEEDS:
                name = f"{role.lower()}_{predictor}_seed{seed}"
                work = campaign / "models" / name
                work.mkdir(parents=True)
                output = work / "prediction"
                if predictor == "af3":
                    path = work / f"{name}.json"
                    data = dict(
                        name=name,
                        dialect="alphafold3",
                        version=3,
                        modelSeeds=[seed],
                        sequences=[
                            dict(
                                protein=dict(
                                    id=list("ABCD"),
                                    sequence=alpha,
                                    unpairedMsa=alpha_msa.read_text(),
                                    pairedMsa="",
                                    templates=[],
                                )
                            ),
                            dict(
                                protein=dict(
                                    id=list("EFGH"),
                                    sequence=beta,
                                    unpairedMsa=beta_msa.read_text(),
                                    pairedMsa="",
                                    templates=[],
                                )
                            ),
                            dict(ligand=dict(id=list("IJKL"), ccdCodes=[COMPONENT])),
                        ],
                        userCCDPath=graph["user_ccd_path"],
                    )
                    expected = output / name / f"{name}_model.cif"
                else:
                    path = work / f"{name}.yaml"
                    data = dict(
                        version=1,
                        sequences=[
                            dict(
                                protein=dict(
                                    id=list("ABCD"),
                                    sequence=alpha,
                                    msa=str(alpha_msa.resolve()),
                                )
                            ),
                            dict(
                                protein=dict(
                                    id=list("EFGH"),
                                    sequence=beta,
                                    msa=str(beta_msa.resolve()),
                                )
                            ),
                            dict(ligand=dict(id=list("IJKL"), smiles=SMILES)),
                        ],
                        templates=[],
                        constraints=[],
                    )
                    expected = (
                        output
                        / f"boltz_results_{name}"
                        / "predictions"
                        / name
                        / f"{name}_model_0.cif"
                    )
                _write(path, data)
                job = dict(
                    job_id=name,
                    role=role,
                    predictor=predictor,
                    seed=seed,
                    input=str(path.resolve()),
                    output=str(output.resolve()),
                    expected_cif=str(expected.resolve()),
                    alpha_sequence=alpha,
                    beta_sequence=beta,
                    samples=1,
                    ligand_atom_maps=graph["ligand_atom_maps"][predictor],
                )
                cmd, env = commands(job)
                job.update(command=cmd, environment=env)
                _write(work / "command.json", job)
                jobs.append(job)
    manifest = dict(
        schema_version=1,
        jobs=jobs,
        model_count=8,
        seeds=list(SEEDS),
        selected_candidate_id=design["candidate_id"],
        input_source_sha256=inputs,
        prepared_sha256={
            str(p.resolve()): digest(p) for p in campaign.rglob("*") if p.is_file()
        },
        graph=graph,
        MSA="Query-onlyA3M;nohomologs;AF3pairedMsaempty",
        templates="emptyforWTanddesignbotharchitectures",
        guidance="None;nocontactconstraints,nocovalentE-Sbond,noprotein/TStemplate;Boltzuse_potentialsdisabled",
        source_sha256={
            str(p): digest(p)
            for p in [
                Path(__file__),
                Path(__file__).with_name("relay_broad_predictor_runtime.py"),
                BOLTZ_WRAPPER,
                AF3 / "run_alphafold.py",
            ]
        },
    )
    _write(manifest_path, manifest)
    LOGGER.info("Preparedexactly8unguidedcomparisonmodels:%s", manifest_path)
    return manifest
