"""Source-pinned local stock Boltz audits and per-case native completion receipts."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
from typing import Any

from nyl12_covalent.matched_original_boltz import (
    candidate,
    journal,
    require_idle,
    sha,
    write_json,
)

NATIVE_FLAGS = [
    "--model",
    "boltz2",
    "--recycling_steps",
    "3",
    "--diffusion_samples",
    "1",
    "--sampling_steps",
    "200",
    "--devices",
    "1",
    "--num_workers",
    "0",
    "--no_kernels",
    "--output_format",
    "mmcif",
]


def audit(plan: dict[str, Any], work: Path) -> dict[str, Any]:
    """Run the full stock native dataset audit in the explicitly selected environment."""
    engine = plan["engine"]
    argv = [
        engine["python"],
        "-m",
        "nyl12_covalent.matched_original_native",
        "--audit",
        str(work / "case.json"),
        engine["cache"],
    ]
    overrides = {
        **engine["env_overrides"],
        "PYTHONPATH": str(Path(__file__).resolve().parents[1]),
        "CUDA_VISIBLE_DEVICES": "",
    }
    with (work / "native_cpu_audit.log").open("w") as handle:
        result = subprocess.run(
            argv,
            env={**os.environ, **overrides},
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    write_json(
        work / "native_cpu_process.json",
        {
            "command": argv,
            "environment_overrides": overrides,
            "exit_status": result.returncode,
        },
    )
    if result.returncode:
        raise RuntimeError(f"Full native dataset CPU audit failed: {work.name}")
    acceptance = json.loads((work / "cpu_acceptance.json").read_text())
    if acceptance["status"] != "accepted" or not acceptance["stock_native_entry_only"]:
        raise ValueError("Full stock native dataset acceptance required")
    return acceptance


def predict_candidate(
    plan: dict[str, Any], output: Path, name: str, seed: int | None = None
) -> None:
    """Predict each fixed matched seed only after an unchanged native dataset audit."""
    item = candidate(plan, name)
    engine = plan["engine"]
    executable, checkpoint = Path(engine["executable"]), Path(engine["checkpoint"])
    if executable.name != "boltz" or sha(checkpoint) != engine["checkpoint_sha256"]:
        raise ValueError("Pinned stock native Boltz executable/checkpoint required")
    if seed is not None and seed not in plan["seeds"]:
        raise ValueError("Selected seed must belong to the explicit matched pair")
    selected_seeds = plan["seeds"] if seed is None else [seed]
    for seed in selected_seeds:
        work = output / "cases" / f"{name}_seed_{seed}"
        case = json.loads((work / "case.json").read_text())
        journal("intent", f"Audit and predict {case['case']} on {item['gpu_uuid']}")
        acceptance = audit(plan, work)
        inputs = [
            Path(case["input"]),
            work / "case.json",
            work / "cpu_acceptance.json",
            Path(case["msas"]) / "alpha_paired.csv",
            Path(case["msas"]) / "beta_paired.csv",
            Path(case["msas"]) / "precursor_pairing_receipt.json",
            Path(acceptance["native_featurizer_source_path"]),
            Path(__file__),
            Path(__file__).with_name("matched_original_native_audit.py"),
        ]
        if sha(Path(case["input"])) != acceptance["input_sha256"]:
            raise ValueError("CPU-qualified input changed")
        if (
            sha(Path(acceptance["native_featurizer_source_path"]))
            != acceptance["native_featurizer_source_sha256"]
        ):
            raise ValueError("CPU-qualified native feature source changed")
        record = (
            Path(case["remote_output"])
            / f"boltz_results_{case['case']}"
            / "processed/records"
            / f"{case['case']}.json"
        )
        if sha(record) != acceptance["record_sha256"]:
            raise ValueError("CPU-qualified template-free record changed")
        if (
            sha(Path(case["msas"]) / "precursor_pairing_receipt.json")
            != acceptance["paired_source_receipt_sha256"]
        ):
            raise ValueError("CPU-qualified paired source receipt changed")
        for filename, digest in case["source_sha256"].items():
            if sha(Path(filename)) != digest:
                raise ValueError(f"Staged scientific source changed: {filename}")
        inputs.append(record)
        require_idle(item["gpu_uuid"])
        argv = [
            str(executable),
            "predict",
            case["input"],
            *NATIVE_FLAGS,
            "--cache",
            engine["cache"],
            "--checkpoint",
            str(checkpoint),
            "--seed",
            str(seed),
            "--out_dir",
            case["remote_output"],
        ]
        overrides = {
            **engine["env_overrides"],
            "CUDA_VISIBLE_DEVICES": item["gpu_uuid"],
        }
        receipt = {
            "status": "native_running",
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "case": case["case"],
            "candidate": name,
            "seed": seed,
            "command": argv,
            "environment_overrides": overrides,
            "gpu_uuid": item["gpu_uuid"],
            "engine_executable_sha256": sha(executable),
            "checkpoint_sha256": engine["checkpoint_sha256"],
            "input_sha256": {str(p): sha(p) for p in inputs},
            "ligand_graph_manifest": plan["ligand_graph_manifest"],
            "approved_reference_manifest": plan["approved_reference_manifest"],
            "target_phi_sign": 1,
            "scientific_acceptance": "Not evaluated; process completion only",
        }
        receipt_path = work / "native_prediction.json"
        if receipt_path.exists():
            raise FileExistsError(
                "Native prediction receipt already exists; preserve attempt"
            )
        write_json(receipt_path, receipt)
        with (work / "native_prediction.log").open("w") as handle:
            result = subprocess.run(
                argv,
                env={**os.environ, **overrides},
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=False,
            )
        receipt.update(
            exit_status=result.returncode,
            ended_utc=datetime.now(timezone.utc).isoformat(),
            status="native_process_completed"
            if result.returncode == 0
            else "native_failed",
        )
        write_json(receipt_path, receipt)
        journal("outcome", f"Native prediction {case['case']} exit {result.returncode}")
        if result.returncode:
            raise RuntimeError(f"Native prediction failed: {case['case']}")
        cif = (
            Path(case["remote_output"])
            / f"boltz_results_{case['case']}"
            / "predictions"
            / case["case"]
            / f"{case['case']}_model_0.cif"
        )
        receipt.update(
            prediction_cif=str(cif),
            prediction_cif_sha256=sha(cif),
            status="native_prediction_completed",
        )
        write_json(receipt_path, receipt)


def main() -> None:
    """Expose only the internal CPU audit subprocess used by the stable runner."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("cache", type=Path)
    args = parser.parse_args()
    from nyl12_covalent.matched_original_native_audit import audit_case

    result = audit_case(json.loads(args.audit.read_text()), args.cache)
    write_json(args.audit.parent / "cpu_acceptance.json", result)


if __name__ == "__main__":
    main()
