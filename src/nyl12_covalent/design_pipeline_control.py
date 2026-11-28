"""Connect accepted sequences and fresh Boltz predictions to OpenMM controls."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np
import parmed

from .design_component_normalization import normalize_model
from .design_geometry_comparison import compare_models
from .design_model_relaxation import relax_model
from .design_pa6_mapping import export_canonical_heavy
from .ligandmpnn_pilot_audit import _pdb_sequences

LOGGER = logging.getLogger(__name__)


def verify_inputs(plan: dict[str, Any]) -> None:
    """Reject changed inputs or an expanded two-model control budget.

    Args:
        plan: Frozen control plan including SHA-256 input identities.
    """
    names = [candidate["id"] for candidate in plan["candidates"]]
    if names != ["WT", "variant_01"] or plan["candidate_count"] != 2:
        raise ValueError("Control requires exactly WT and variant_01")
    if plan["physical_gpu"] != 0 or plan["new_sequence_samples"] != 0:
        raise ValueError("Control is GPU0 only and reuses accepted sequences")
    for raw, expected in plan["input_sha256"].items():
        path = Path(raw)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Input hash changed: {path}")


def wait_for_predictions(root: Path, timeout_seconds: float = 21600) -> None:
    """Wait for the specific submitted prediction job, failing on its exit.

    Args:
        root: Logical control directory containing prediction_job/status.json.
        timeout_seconds: Maximum wall time before a recorded failure.
    """
    status_path = root / "prediction_job/status.json"
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        status = json.loads(status_path.read_text())
        code = status["exit_status"]
        if code is not None:
            if code != 0:
                raise RuntimeError(f"Prediction job exited {code}: {status_path}")
            marker = json.loads((root / "predictions/WORKER_COMPLETE.json").read_text())
            if marker.get("jobs") != 2 or marker.get("complete") != 2:
                raise ValueError("Two-model prediction completeness mismatch")
            return
        LOGGER.info("Waiting for the two submitted Boltz predictions")
        time.sleep(60)
    raise TimeoutError("Two-model prediction job exceeded the control deadline")


def run_control(root: Path, execution_config: Path) -> dict[str, Any]:
    """Prepare, minimize and compare the two explicitly authorized controls.

    Args:
        root: Absolute logical control directory with frozen execution config.
        execution_config: Explicit, immutable attempt-specific execution config.

    Returns:
        Process, geometry and independent force-convergence summary.
    """
    plan = json.loads((root / "plan.json").read_text())
    config = json.loads(execution_config.read_text())
    execution_plan = dict(plan)
    execution_plan["input_sha256"] = {
        **plan["input_sha256"],
        **config["dependency_sha256"],
    }
    verify_inputs(execution_plan)
    wait_for_predictions(root)
    relaxation = json.loads(Path(config["relaxation_config"]).read_text())
    relaxation["gpu"] = "0"
    entries, cases = [], []
    for candidate in plan["candidates"]:
        name = candidate["id"]
        LOGGER.info("Preparing and minimizing %s with OpenMM", name)
        case = root / "models" / name
        case.mkdir(parents=True, exist_ok=False)
        source = root / f"predictions/canonical/{name}/{name}_model_0.pdb"
        component = root / (
            f"predictions/runs/{name}/model_0/output/boltz_results_{name}/"
            f"processed/mols/{name}.pkl"
        )
        normalize_model(source, component, Path(relaxation["mol2"]), case / "raw")
        raw = case / "raw/normalized.pdb"
        sequences = _pdb_sequences(raw)
        for chains, key in (("ABCD", "alpha_sequence"), ("EFGH", "beta_sequence")):
            if any(sequences[chain] != candidate[key] for chain in chains):
                raise ValueError(f"Predicted protein sequence mismatch: {name}")
        prepared = relax_model(
            raw, case / "preparation", relaxation, preparation_only=True
        )
        output = case / "openmm"
        argv = [
            config["openmm_python"],
            "-m",
            "nyl12_covalent.design_openmm_control",
            "--topology",
            prepared["topology"],
            "--coordinates",
            prepared["coordinates"],
            "--output",
            str(output),
        ]
        env = os.environ.copy()
        env.update(config["environment"])
        (case / "openmm_invocation.json").write_text(
            json.dumps({"argv": argv, "environment": config["environment"]}, indent=2)
            + "\n"
        )
        with (case / "openmm.log").open("x") as log:
            subprocess.run(
                argv,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                timeout=21600,
            )
        acceptance = json.loads((output / "acceptance.json").read_text())
        endpoint = parmed.load_file(prepared["topology"], xyz=prepared["coordinates"])
        endpoint.coordinates = np.load(output / "positions_nm.npy") * 10.0
        endpoint.strip(":WAT,Na+,Cl-")
        top, coords = case / "endpoint.prmtop", case / "endpoint.rst7"
        endpoint.save(str(top), overwrite=False)
        endpoint.save(str(coords), overwrite=False)
        relaxed = case / "relaxed.pdb"
        export_canonical_heavy(top, coords, Path(prepared["conversion"]), relaxed)
        for stage, path in (("raw", raw), ("relaxed", relaxed)):
            entries.append(
                dict(
                    variant=name,
                    seed=plan["matched_seed"],
                    stage=stage,
                    path=str(path),
                    ligand_count=4,
                )
            )
        cases.append({"variant": name, "openmm_acceptance": acceptance})
    geometry = compare_models(entries, Path(config["reference_cif"]), root / "analysis")
    summary = {
        "process_completed": True,
        "cases": cases,
        "geometry": geometry,
        "scientific_scope": plan["scientific_scope"],
        "claim_of_improved_activity": False,
        "new_md_or_quantum_calculations": False,
    }
    (root / "COMPLETED.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    """Run the explicitly frozen two-model control from its logical directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--execution-config", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    run_control(args.root.resolve(), args.execution_config.resolve())


if __name__ == "__main__":
    main()
