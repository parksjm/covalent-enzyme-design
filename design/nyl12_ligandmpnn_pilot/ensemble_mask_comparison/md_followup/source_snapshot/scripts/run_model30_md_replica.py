#!/usr/bin/env python3
"""Execute one matched10ns PA6 pilot from an accepted restart."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.pilot_md_protocol import mdin_files, protocol_manifest  # noqa: E402

LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Write an explicit run record to a stage-owned JSON file."""
    path.write_text(json.dumps(data, indent=2) + "\n")


def stage_acceptance(path: Path, steps: int) -> dict[str, Any]:
    """Require complete finite MD output with a reasonable final temperature.

    Args:
        path: Amber output file.
        steps: Required complete timestep count.

    Returns:
        Stage completion and final temperature observations.
    """
    text = path.read_text()
    pattern = (
        r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([-+0-9.E]+)"
        r"\s+TEMP\(K\)\s*=\s*([-+0-9.E]+)"
    )
    instantaneous = text.split("A V E R A G E S", 1)[0]
    blocks = re.findall(pattern, instantaneous)
    if not blocks or int(blocks[-1][0]) != steps:
        raise ValueError(f"Missing final step {steps}: {path}")
    temperature = float(blocks[-1][2])
    if not 250 < temperature < 350:
        raise ValueError(f"Final temperature outside250-350K: {temperature}")
    if re.search(
        r"\b(?:NaN|Infinity)\b|SHAKE cannot be accomplished|Error:|ERROR:",
        text,
        re.IGNORECASE,
    ):
        raise ValueError(f"Numeric/engine error: {path}")
    return {
        "completed_steps": steps,
        "final_time_ps": float(blocks[-1][1]),
        "final_temperature_K": temperature,
        "temperature_gate_pass": True,
        "scientific_scope": "Stage completion does not prove binding stability or equilibration",
    }


def main() -> None:
    """Run a new independent replica and retain every stage and acceptance record."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("topology", "coordinates", "engine", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument(
        "--resume-after", help="Explicit completed stage to continue after"
    )
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.resume_after is None:
        args.output.mkdir(parents=True, exist_ok=False)
    elif not args.output.is_dir():
        raise FileNotFoundError(args.output)
    inputs = mdin_files(args.seed)
    manifest = protocol_manifest(args.seed)
    sources = (
        args.topology,
        args.coordinates,
        args.engine,
        Path(__file__),
        Path(__file__).resolve().parents[1] / "src/nyl12_covalent/pilot_md_protocol.py",
    )
    manifest.update(
        label=args.label,
        executable=str(args.engine),
        input_sha256={str(p): digest(p) for p in sources},
    )
    if args.resume_after is None:
        write_json(args.output / "protocol.json", manifest)
    else:
        original = json.loads((args.output / "protocol.json").read_text())
        if original["stages"] != manifest["stages"] or original["label"] != args.label:
            raise ValueError("Continuation protocol or replica identity mismatch")
    topology = args.output / "system.prmtop"
    reference = args.output / "starting_minimized.rst7"
    if args.resume_after is None:
        shutil.copy2(args.topology, topology)
        shutil.copy2(args.coordinates, reference)
    elif digest(topology) != digest(args.topology) or digest(reference) != digest(
        args.coordinates
    ):
        raise ValueError("Continuation input topology or original reference changed")
    coordinates = reference
    stages = manifest["stages"]
    if args.resume_after is not None:
        names = [stage["name"] for stage in stages]
        completed = names.index(args.resume_after)
        previous = stages[completed]
        acceptance = stage_acceptance(
            args.output / f"{args.resume_after}.out", previous["steps"]
        )
        coordinates = args.output / f"{args.resume_after}.rst7"
        acceptance["restart_sha256"] = digest(coordinates)
        continuation = {
            "checkpoint_stage": args.resume_after,
            "checkpoint_acceptance": acceptance,
            "source_sha256": manifest["input_sha256"],
            "started_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_json(
            args.output / f"continuation_after_{args.resume_after}.json", continuation
        )
        stages = stages[completed + 1 :]
    for stage in stages:
        if (args.output / f"{stage['name']}.out").exists():
            raise FileExistsError("Refusing to overwrite existing MD stage output")
        name = stage["name"]
        mdin = args.output / f"{name}.in"
        mdin.write_text(inputs[name])
        restart = args.output / f"{name}.rst7"
        command = [
            str(args.engine),
            "-O",
            "-i",
            str(mdin),
            "-o",
            str(args.output / f"{name}.out"),
            "-p",
            str(topology),
            "-c",
            str(coordinates),
            "-r",
            str(restart),
            "-x",
            str(args.output / f"{name}.nc"),
            "-inf",
            str(args.output / f"{name}.mdinfo"),
        ]
        if stage["heavy_restraint_kcal_mol_A2"] is not None:
            command.extend(("-ref", str(reference)))
        record = {
            "command": command,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "input_sha256": {str(p): digest(p) for p in (mdin, topology, coordinates)},
            "reference_sha256": digest(reference),
        }
        write_json(args.output / f"{name}.process.json", record)
        LOGGER.info("Starting %s %s", args.label, name)
        with (args.output / f"{name}.engine.log").open("w") as handle:
            result = subprocess.run(
                command,
                cwd=args.output,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=False,
            )
        record.update(
            returncode=result.returncode,
            ended_utc=datetime.now(timezone.utc).isoformat(),
        )
        write_json(args.output / f"{name}.process.json", record)
        if result.returncode:
            raise RuntimeError(f"{args.label} {name} exited {result.returncode}")
        acceptance = stage_acceptance(args.output / f"{name}.out", stage["steps"])
        acceptance["restart_sha256"] = digest(restart)
        write_json(args.output / f"{name}.acceptance.json", acceptance)
        coordinates = restart
        LOGGER.info(
            "Completed %s %s finalT=%.2fK",
            args.label,
            name,
            acceptance["final_temperature_K"],
        )
    write_json(
        args.output / "completion.json",
        {
            "status": "10ns_pilot_complete_requires_projected_phi_native_analysis",
            "last_restart_sha256": digest(coordinates),
        },
    )


if __name__ == "__main__":
    main()
