#!/usr/bin/env python3
"""Run an unrestrained E-S continuation followed by a bounded relay path."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.unrestrained_relay import (  # noqa: E402
    coordinates_valid,
    distances,
    driven_input,
    es_input,
    native_converged,
    PHASE2_CORRECTION_REASON,
    restraint_text,
    schedule,
    sha256,
    stage2_schedule,
    write_json,
)

LOGGER = logging.getLogger(__name__)
SOURCE = (
    ROOT
    / "output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/phase1_unrestrained"
)
BASE = ROOT / "qmmm/pa6_asp272"
ES = BASE / "00_es_unrestrained"
PATH = BASE / "01_coordinate_drive_to_ti1"
TI1 = BASE / "02_ti1_unrestrained_validation"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _command(folder: Path, topology: Path, start: Path) -> list[str]:
    return [
        "/data/jpc/code/ambertools26/venv/bin/sander",
        "-O",
        "-i",
        str(folder / "minimize.mdin"),
        "-p",
        str(topology),
        "-c",
        str(start),
        "-o",
        str(folder / "minimized.mdout"),
        "-r",
        str(folder / "optimized.ncrst"),
        "-x",
        str(folder / "optimization.nc"),
        "-inf",
        str(folder / "mdinfo"),
    ]


def _run_engine(folder: Path, command: list[str], environment: dict[str, str]) -> int:
    record: dict[str, Any] = {"command": command, "started_utc": _now()}
    write_json(folder / "process.json", record)
    env = os.environ.copy()
    env.update(environment)
    with (folder / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=folder,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(exit_status=result.returncode, ended_utc=_now())
    write_json(folder / "process.json", record)
    return result.returncode


def prepare() -> None:
    """Freeze the prior checkpoint and prepare the corrected sequential run."""
    if ES.exists() or PATH.exists() or TI1.exists():
        raise FileExistsError("Corrected run directories already exist")
    preparation = json.loads((SOURCE / "preparation.json").read_text())
    topology_source = Path(preparation["topology"])
    start_source = SOURCE / "optimized.ncrst"
    input_source = SOURCE / "minimize.mdin"
    if not all(path.exists() for path in (topology_source, start_source, input_source)):
        raise FileNotFoundError("Required prior checkpoint is missing")
    system = BASE / "system"
    system.mkdir(parents=True, exist_ok=True)
    topology = system / "md_solvated.prmtop"
    if topology.exists() and sha256(topology) != sha256(topology_source):
        raise ValueError("Existing system topology has a different identity")
    if not topology.exists():
        shutil.copy2(topology_source, topology)
    ES.mkdir()
    shutil.copy2(start_source, ES / "starting.ncrst")
    base = input_source.read_text()
    (ES / "minimize.mdin").write_text(es_input(base))
    environment = preparation["environment"]
    initial = distances(topology, ES / "starting.ncrst")
    write_json(
        ES / "inputs.json",
        {
            "topology": str(topology),
            "starting_coordinates": str(ES / "starting.ncrst"),
            "source_checkpoint": str(start_source),
            "source_sha256": {
                str(topology_source): sha256(topology_source),
                str(start_source): sha256(start_source),
                str(input_source): sha256(input_source),
                str(Path(__file__)): sha256(Path(__file__)),
            },
            "environment": environment,
            "box": preparation["box"],
            "initial_distances_A": initial,
            "restraints": "none; ntr=0, nmropt=0, all atoms active",
        },
    )
    PATH.mkdir()
    write_json(
        PATH / "inputs.json",
        {
            "topology": str(topology),
            "initial": str(ES / "optimized.ncrst"),
            "box": preparation["box"],
            "environment": environment,
            "windows_per_stage": 21,
            "positional_restraints": "none",
            "driving_biases": [
                "d(OG1-C)",
                "stage1: d(OG1-HG1)-d(Asp272_OD1-HG1)",
                "stage2: d(Asp272_OD1-HG1)-d(Thr_N-HG1)",
            ],
        },
    )
    write_json(
        PATH / "status.json",
        {"state": "waiting_for_es", "completed_windows": [], "updated_utc": _now()},
    )
    TI1.mkdir()
    write_json(
        TI1 / "inputs.json",
        {
            "topology": str(topology),
            "source": "Final stage-2 coordinate-driving endpoint",
            "environment": environment,
            "restraints": "none; TI1 must persist after every bias is released",
        },
    )


def _validate_es(topology: Path, endpoint: Path) -> dict[str, float]:
    values = distances(topology, endpoint)
    valid = (
        values["OG1_C"] > 2.5
        and values["OG1_HG1"] < 1.2
        and values["Asp272_OD1_HG1"] > 1.3
        and values["Thr_N_H1"] < 1.25
        and values["Thr_N_H2"] < 1.25
        and 1.15 < values["C_O"] < 1.4
        and 1.2 < values["C_N"] < 1.5
    )
    if not valid:
        raise ValueError("Converged endpoint does not pass the noncovalent E-S gate")
    return values


def _run_es(inputs: dict[str, Any]) -> dict[str, float]:
    topology = Path(inputs["topology"])
    command = _command(ES, topology, ES / "starting.ncrst")
    if _run_engine(ES, command, inputs["environment"]):
        raise RuntimeError("Unrestrained E-S continuation failed")
    converged = native_converged(ES / "minimized.mdout")
    valid = coordinates_valid(topology, ES / "optimized.ncrst", inputs["box"])
    values = distances(topology, ES / "optimized.ncrst")
    write_json(
        ES / "acceptance.json",
        {
            "native_converged": converged,
            "finite_coordinates_and_fixed_box": valid,
            "distances_A": values,
            "restraints": "none",
        },
    )
    if not (converged and valid):
        raise ValueError("Unrestrained E-S continuation did not converge")
    return _validate_es(topology, ES / "optimized.ncrst")


def _run_window(
    row: dict[str, Any], restart: Path, topology: Path, inputs: dict[str, Any]
) -> Path:
    folder = PATH / row["label"]
    folder.mkdir()
    (folder / "minimize.mdin").write_text(
        driven_input((ES / "minimize.mdin").read_text())
    )
    (folder / "restraints.RST").write_text(
        restraint_text(
            row["stage"],
            row["attack_target_A"],
            row["proton_difference_target_A"],
            (row["relay_retention_target_A"] if row["stage"] == 2 else None),
        )
    )
    command = _command(folder, topology, restart)
    process = {**row, "command": command, "input_restart": str(restart)}
    write_json(folder / "process.json", process)
    if _run_engine(folder, command, inputs["environment"]):
        raise RuntimeError(f"Engine failed in {row['label']}")
    process.update(json.loads((folder / "process.json").read_text()))
    process.update(row)
    process["input_restart"] = str(restart)
    write_json(folder / "process.json", process)
    converged = native_converged(folder / "minimized.mdout")
    valid = coordinates_valid(topology, folder / "optimized.ncrst", inputs["box"])
    values = distances(topology, folder / "optimized.ncrst")
    actual = (
        values["OG1_HG1"] - values["Asp272_OD1_HG1"]
        if row["stage"] == 1
        else values["Asp272_OD1_HG1"] - values["Thr_N_HG1"]
    )
    relay_retention_actual = values["OG1_HG1"] - values["Asp272_OD1_HG1"]
    proton_ownership_valid = row["stage"] == 1 or (
        values["OG1_HG1"] >= row["og1_hg1_minimum_A"]
    )
    write_json(
        folder / "acceptance.json",
        {
            "native_converged": converged,
            "finite_coordinates_and_fixed_box": valid,
            "tracking_is_diagnostic_only": True,
            "attack_error_A": abs(values["OG1_C"] - row["attack_target_A"]),
            "proton_coordinate_error_A": abs(
                actual - row["proton_difference_target_A"]
            ),
            "relay_retention_error_A": (
                abs(relay_retention_actual - row["relay_retention_target_A"])
                if row["stage"] == 2
                else None
            ),
            "proton_ownership_valid": proton_ownership_valid,
            "distances_A": values,
            "positional_restraints": "none",
        },
    )
    if not (converged and valid and proton_ownership_valid):
        raise ValueError(
            f"Window {row['label']} failed convergence, coordinates, or "
            "proton-ownership acceptance"
        )
    return folder / "optimized.ncrst"


def _run_ti1_release(
    restart: Path, topology: Path, inputs: dict[str, Any]
) -> dict[str, Any]:
    """Release every construction bias and test whether TI1 is a minimum."""
    shutil.copy2(restart, TI1 / "starting.ncrst")
    (TI1 / "minimize.mdin").write_text(es_input((ES / "minimize.mdin").read_text()))
    command = _command(TI1, topology, TI1 / "starting.ncrst")
    if _run_engine(TI1, command, inputs["environment"]):
        raise RuntimeError("Unrestrained TI1 validation failed")
    converged = native_converged(TI1 / "minimized.mdout")
    valid = coordinates_valid(topology, TI1 / "optimized.ncrst", inputs["box"])
    values = distances(topology, TI1 / "optimized.ncrst")
    retained = (
        1.3 < values["OG1_C"] < 1.9
        and values["OG1_HG1"] > 1.3
        and values["Asp272_OD1_HG1"] > 1.3
        and values["Thr_N_HG1"] < 1.25
        and 1.2 < values["C_O"] < 1.6
        and 1.2 < values["C_N"] < 1.7
    )
    acceptance = {
        "native_converged": converged,
        "finite_coordinates_and_fixed_box": valid,
        "ti1_geometry_retained": retained,
        "distances_A": values,
        "restraints": "none",
        "interpretation": (
            "A retained, converged tetrahedral state is a TI1 minimum candidate; "
            "the driven path itself is not a barrier or TS calculation."
        ),
    }
    write_json(TI1 / "acceptance.json", acceptance)
    if not (converged and valid):
        raise ValueError("Unrestrained TI1 validation did not converge")
    return acceptance


def _run_path_from_endpoint(endpoint: Path, og1_hg1_minimum: float) -> None:
    """Run coordinate driving from one accepted, converged E-S endpoint."""
    path_inputs = json.loads((PATH / "inputs.json").read_text())
    status = json.loads((PATH / "status.json").read_text())
    resumable_stopped = (
        status["state"] == "stopped"
        and status.get("reason")
        in {
            "Unrestrained E-S continuation did not converge",
            "Engine failed in stage1_00",
        }
        and not status["completed_windows"]
    )
    resumable_running = status["state"] == "running" and bool(
        status["completed_windows"]
    )
    resumable_corrected_phase2 = (
        status["state"] == "stopped"
        and status.get("reason") == PHASE2_CORRECTION_REASON
        and status["completed_windows"]
        == [f"stage1_{index:02d}" for index in range(21)]
    )
    if (
        status["state"] != "waiting_for_es"
        and not resumable_stopped
        and not resumable_running
        and not resumable_corrected_phase2
    ):
        raise ValueError("Coordinate-driving workflow has already been started")
    topology = Path(path_inputs["topology"])
    if not coordinates_valid(topology, endpoint, path_inputs["box"]):
        raise ValueError("Accepted E-S endpoint has invalid coordinates or cell")
    initial = _validate_es(topology, endpoint)
    try:
        status.update(
            state="running",
            reason=None,
            es_endpoint=str(endpoint),
            es_endpoint_sha256=sha256(endpoint),
            es_endpoint_distances_A=initial,
            updated_utc=_now(),
        )
        write_json(PATH / "status.json", status)
        restart = endpoint
        rows = schedule(initial, path_inputs["windows_per_stage"])
        for row in rows:
            if row["label"] in status["completed_windows"]:
                restart = PATH / row["label"] / "optimized.ncrst"
                continue
            restart = _run_window(row, restart, topology, path_inputs)
            status["completed_windows"].append(row["label"])
            status["updated_utc"] = _now()
            write_json(PATH / "status.json", status)
        stage1_final = distances(topology, restart)
        q2 = stage1_final["Asp272_OD1_HG1"] - stage1_final["Thr_N_HG1"]
        relay_retention = stage1_final["OG1_HG1"] - stage1_final["Asp272_OD1_HG1"]
        for row in stage2_schedule(
            q2,
            relay_retention,
            og1_hg1_minimum,
            path_inputs["windows_per_stage"],
        ):
            if row["label"] in status["completed_windows"]:
                restart = PATH / row["label"] / "optimized.ncrst"
                continue
            restart = _run_window(row, restart, topology, path_inputs)
            status["completed_windows"].append(row["label"])
            status["updated_utc"] = _now()
            write_json(PATH / "status.json", status)
        status.update(state="validating_unrestrained_ti1", updated_utc=_now())
        write_json(PATH / "status.json", status)
        ti1_acceptance = _run_ti1_release(restart, topology, path_inputs)
        status.update(
            state="completed",
            ti1_validation=str(TI1 / "acceptance.json"),
            ti1_geometry_retained=ti1_acceptance["ti1_geometry_retained"],
            updated_utc=_now(),
        )
    except Exception as error:
        status.update(state="stopped", reason=str(error), updated_utc=_now())
        raise
    finally:
        write_json(PATH / "status.json", status)


def run(og1_hg1_minimum: float) -> None:
    """Run E-S convergence, then two bounded coordinate-driving stages."""
    es_inputs = json.loads((ES / "inputs.json").read_text())
    _run_es(es_inputs)
    _run_path_from_endpoint(ES / "optimized.ncrst", og1_hg1_minimum)


def main() -> None:
    """Prepare or execute the corrected sequential workflow."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "run", "resume-path"))
    parser.add_argument("--endpoint", type=Path)
    parser.add_argument("--stage2-og1-hg1-minimum", type=float)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.mode == "prepare":
        prepare()
    else:
        if args.stage2_og1_hg1_minimum is None:
            parser.error("run and resume-path require --stage2-og1-hg1-minimum")
        if args.mode == "run":
            run(args.stage2_og1_hg1_minimum)
            return
        if args.endpoint is None:
            parser.error("resume-path requires --endpoint")
        _run_path_from_endpoint(
            args.endpoint.resolve(),
            args.stage2_og1_hg1_minimum,
        )


if __name__ == "__main__":
    main()
