"""Sequential native runtime for the gated eight-model predictor comparison."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from .relay_broad_graph import ROOT, digest
from .relay_broad_prediction_io import audit_prediction, export_geometry

AF3 = Path("/data/jpc/code/alphafold3-contact-restraints")
AF3_PYTHON = Path("/data/jpc/code/.conda/envs/af3/bin/python")
BOLTZ_PYTHON = Path("/data/jpc/enzyme-substrate-modeling/.pixi/envs/boltz/bin/python")
BOLTZ_WRAPPER = (
    ROOT
    / "design/nyl12_ligandmpnn_pilot/boltz_validation/source_snapshot/run_boltz_cli.py"
)


def _write(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def commands(job: dict) -> tuple[list[str], dict[str, str]]:
    """Return exact one-sample native CLI and explicit two-thread GPU1 environment.

    Args:
        job: Prepared single architecture/sequence/seed manifest entry.

    Returns:
        Argument vector and nonsecret environment overrides; no execution.
    """
    env = dict(
        CUDA_VISIBLE_DEVICES="1",
        OMP_NUM_THREADS="2",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="2",
        NUMEXPR_NUM_THREADS="2",
    )
    if job["predictor"] == "af3":
        purelib = "/data/jpc/code/.conda/envs/af3/lib/python3.11/site-packages"
        env.update(
            PYTHONPATH=f"{AF3}/src:{purelib}",
            XLA_PYTHON_CLIENT_PREALLOCATE="false",
            PATH="/data/jpc/code/.conda/envs/af3/bin:/usr/local/bin:/usr/bin:/bin",
        )
        cmd = [
            str(AF3_PYTHON),
            "-S",
            str(AF3 / "run_alphafold.py"),
            "--json_path",
            job["input"],
            "--output_dir",
            job["output"],
            "--model_dir",
            "/data/jpc/code/alphafold3/models",
            "--num_diffusion_samples=1",
            "--num_recycles=3",
            "--run_data_pipeline=false",
            "--buckets=1536",
            "--gpu_device=0",
            "--contact_guidance_weight=0",
            "--jax_compilation_cache_dir",
            str(Path(job["output"]).parent / "jax_cache"),
        ]
    elif job["predictor"] == "boltz":
        env["PYTHONPATH"] = str(BOLTZ_WRAPPER.parent)
        cmd = [
            str(BOLTZ_PYTHON),
            str(BOLTZ_WRAPPER),
            "predict",
            job["input"],
            "--model",
            "boltz2",
            "--output_format",
            "mmcif",
            "--cache",
            "/data/jpc/code/boltz",
            "--recycling_steps",
            "3",
            "--sampling_steps",
            "200",
            "--diffusion_samples",
            "1",
            "--seed",
            str(job["seed"]),
            "--out_dir",
            job["output"],
            "--no_kernels",
            "--num_workers",
            "0",
        ]
    else:
        raise ValueError("Unsupportedpredictor")
    return cmd, env


def run_comparison(manifest: dict, campaign: Path) -> dict:
    """Run the externally gated eight-model comparison sequentially on GPU1.

    Args:
        manifest: Source-verified eight-job preparation manifest.
        campaign: Matching prepared prediction directory.

    Returns:
        Aggregate software/graph acceptance; geometry utility is evaluated separately.
    """
    if len(manifest["jobs"]) != 8 or any(j["samples"] != 1 for j in manifest["jobs"]):
        raise ValueError("Exact8one-samplebudgetrequired")
    expected = {
        (role, predictor, seed)
        for role in ("WT", "DESIGN")
        for predictor in ("af3", "boltz")
        for seed in (20261126, 20261127)
    }
    if {(j["role"], j["predictor"], j["seed"]) for j in manifest["jobs"]} != expected:
        raise ValueError("ExactWT-design-twoarchitecture-twoseedbudgetrequired")
    for job in manifest["jobs"]:
        command, environment = commands(job)
        if job["command"] != command or job["environment"] != environment:
            raise ValueError("Declarednativecommandorenvironmentchanged")
    results = []
    for job in manifest["jobs"]:
        work = Path(job["input"]).parent
        acceptance = work / "software_acceptance.json"
        if acceptance.exists():
            record = json.loads(acceptance.read_text())
            if record["cif_sha256"] != digest(Path(job["expected_cif"])):
                raise ValueError("Acceptedpredictionchanged")
            geometry = export_geometry(
                Path(job["expected_cif"]), job["ligand_atom_maps"], work
            )
            record["geometry"] = geometry
            results.append(record)
            continue
        journal = [
            str(ROOT / ".venv/bin/python"),
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_broad_predictions",
            "--permission-mode",
            "escalated",
            "--reason",
            "Externallygated8modelunguidedcomparison",
            "--rollback",
            "Retainnativepredictoroutputsandlogs;stoponlyverifiedpredictionprocess",
            "--action",
            job["job_id"],
        ]
        subprocess.run([*journal, "--phase", "intent"], check=True)
        env = os.environ.copy()
        env.update(job["environment"])
        with (work / "native.log").open("xb") as log:
            result = subprocess.run(
                job["command"],
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        _write(
            work / "process.json",
            dict(
                exit_status=result.returncode,
                command=job["command"],
                environment=job["environment"],
            ),
        )
        subprocess.run(
            [
                *journal,
                "--phase",
                "outcome",
                "--result",
                f"Nativepredictorexit{result.returncode}",
            ],
            check=True,
        )
        if result.returncode != 0:
            raise RuntimeError("Nativepredictionfailed:" + job["job_id"])
        cif = Path(job["expected_cif"])
        record = audit_prediction(
            cif,
            job["alpha_sequence"],
            job["beta_sequence"],
            job["ligand_atom_maps"]["I"],
        )
        record.update(
            job_id=job["job_id"],
            predictor=job["predictor"],
            role=job["role"],
            seed=job["seed"],
            cif=str(cif),
            cif_sha256=digest(cif),
        )
        _write(acceptance, record)
        record["geometry"] = export_geometry(cif, job["ligand_atom_maps"], work)
        results.append(record)
    final = dict(
        model_count=8,
        software_acceptance=True,
        models=results,
        geometric_utility="Separateany-amidemetricsJSONand48rowCSV/model;noWTcontactgate",
    )
    _write(campaign / "prediction_acceptance.json", final)
    return final
