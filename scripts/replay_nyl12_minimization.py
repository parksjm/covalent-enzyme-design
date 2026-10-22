"""Replay a retained local QM/MM stage with actual trajectory output enabled."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.local_repair_engine import run_checked_engine  # noqa: E402
from nyl12_covalent.local_repair_model import export_pdb  # noqa: E402
from nyl12_covalent.local_repair_checks import audit_local_model  # noqa: E402
from nyl12_covalent.ti1_attempt import _sha, _snapshot  # noqa: E402


def main() -> None:
    """Retain exact source inputs and rerun with coordinate output only changed."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calculation", type=Path, required=True)
    parser.add_argument("--amber-bin", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    source, output = args.calculation.absolute(), args.output.absolute()
    stage = source / "01_hydrogens_water"
    output.mkdir(parents=True, exist_ok=False)
    inventory = {}
    for destination, path in {
        "model.prmtop": source / "model.prmtop",
        "start.rst7": stage / "start.rst7",
        "reference.rst7": source / "prepared.rst7",
        "bridge.RST": stage / "bridge.RST",
        "original_optimized.rst7": stage / "optimized.rst7",
        "original.mdin": stage / "minimize.mdin",
        "original_command.json": stage / "command.json",
        "original_manifest.json": source / "manifest.json",
    }.items():
        shutil.copy2(path, output / destination)
        inventory[destination] = {
            "source": str(path),
            "sha256": _sha(output / destination),
        }
    mdin = (output / "original.mdin").read_text()
    if mdin.count("ntwr=50,") != 1 or "ntwx" in mdin or "ioutfm" in mdin:
        raise ValueError("Unexpected source output settings; inspect before replay")
    mdin = mdin.replace("ntwr=50,", "ntwr=50, ntwx=1, ioutfm=0,")
    (output / "minimize.mdin").write_text(mdin)
    command = [
        str((args.amber_bin / "sander").absolute()),
        "-i",
        str(output / "minimize.mdin"),
        "-o",
        str(output / "minimize.mdout"),
        "-p",
        str(output / "model.prmtop"),
        "-c",
        str(output / "start.rst7"),
        "-ref",
        str(output / "reference.rst7"),
        "-r",
        str(output / "optimized.rst7"),
        "-x",
        str(output / "optimization.mdcrd"),
        "-inf",
        str(output / "minimize.mdinfo"),
    ]
    original_command = json.loads((stage / "command.json").read_text())
    original_manifest = json.loads((source / "manifest.json").read_text())
    if _sha(Path(command[0])) != original_manifest["sander_sha256"]:
        raise ValueError("Sander executable changed since original run")
    for path, digest in original_manifest["engine_data_sha256"].items():
        if _sha(Path(path)) != digest:
            raise ValueError("DFTB parameter data changed since original run")
    overrides = original_command["environment_overrides"]
    environment = os.environ.copy()
    environment.update(overrides)
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": inventory,
                "source_sha256": _snapshot(output),
                "argv": command,
                "environment_overrides": overrides,
                "changes": {"ntwx": 1, "ioutfm": 0, "trajectory": "optimization.mdcrd"},
                "sampling": "Sander minimizer force-call snapshots; may include line-search/trial points; not MD time",
            },
            indent=2,
        )
        + "\n"
    )
    process = run_checked_engine(
        command,
        cwd=output,
        env=environment,
        log_path=output / "sander.log",
        mdout=output / "minimize.mdout",
    )
    (output / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    if process["exit_status"] or process["stopped_on_scf_failure"]:
        raise RuntimeError("Trajectory replay engine failed; retained logs")
    model = AmberParm(str(output / "model.prmtop"), xyz=str(output / "optimized.rst7"))
    original = AmberParm(
        str(output / "model.prmtop"), xyz=str(output / "original_optimized.rst7")
    )
    difference = np.linalg.norm(model.coordinates - original.coordinates, axis=1)
    report = audit_local_model(model)
    report["replay_endpoint_max_difference_A"] = float(difference.max())
    report["replay_endpoint_coordinate_rms_difference_A"] = float(
        np.sqrt(np.mean(difference**2))
    )
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    export_pdb(model, output / "optimized.pdb")
    logging.info(
        "Replay complete; maximum endpoint difference %.8f A", difference.max()
    )


if __name__ == "__main__":
    main()
