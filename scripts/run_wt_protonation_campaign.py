#!/usr/bin/env python3
"""Prepare or execute accepted WT minimizations, then delegate twelve MD jobs."""

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

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "student_handoff/nyl12_ash274_md/tools"))
from nyl12_covalent.wt_protonation_minimization import (  # noqa:E402
    MIN_INPUT,
    STATES,
    configs,
    coordinate_acceptance,
    digest,
    native_acceptance,
    protected_rows,
    write_json,
)
from prepared_md_support import require_gpu  # noqa:E402

SOURCES = (
    "scripts/run_model30_md_replica.py",
    "scripts/continue_model30_md.py",
    "scripts/queue_prepared_md.py",
    "scripts/run_wt_protonation_campaign.py",
    "src/nyl12_covalent/__init__.py",
    "src/nyl12_covalent/pilot_md_protocol.py",
    "src/nyl12_covalent/jobs.py",
    "src/nyl12_covalent/prepared_md_queue_plan.py",
    "src/nyl12_covalent/wt_protonation_minimization.py",
    *(
        "student_handoff/nyl12_ash274_md/tools/" + n
        for n in (
            "run_prepared_md.py",
            "prepared_md_support.py",
            "prepared_terminal_state.py",
        )
    ),
)
ENGINE = Path("/data/jpc/code/pmemd26/bin/pmemd.cuda")
LOGGER = logging.getLogger(__name__)


def prepare(root: Path, repo: Path) -> None:
    """Freeze the exact minimal executable tree and minimization commands only."""
    manifest = json.loads((root / "manifest.json").read_text())
    if set(manifest["states"]) != set(STATES):
        raise ValueError("Three accepted states required")
    frozen = root / "frozen_protocol"
    frozen.mkdir()
    source_hashes = {}
    for name in SOURCES:
        source = repo / name
        target = frozen / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        source_hashes[str(target)] = digest(target)
    plans = []
    engine = ENGINE
    if not os.access(engine, os.X_OK):
        raise FileNotFoundError(engine)
    for i, state in enumerate(STATES):
        folder = root / state
        accepted = json.loads((folder / "preparation_acceptance.json").read_text())
        if not accepted["preparation_accepted"] or not accepted["process_complete"]:
            raise ValueError(state)
        for name, key in [
            ("topology.prmtop", "topology_sha256"),
            ("start.rst7", "coordinates_sha256"),
        ]:
            if digest(folder / name) != accepted[key]:
                raise ValueError(f"Prepared source changed: {state}/{name}")
        output = folder / "unrestrained_minimization"
        output.mkdir()
        mdin = output / "minimize.mdin"
        mdin.write_text(MIN_INPUT)
        command = [
            str(engine),
            "-O",
            "-i",
            str(mdin),
            "-p",
            str(folder / "topology.prmtop"),
            "-c",
            str(folder / "start.rst7"),
            "-o",
            str(output / "minimized.mdout"),
            "-r",
            str(output / "minimized.ncrst"),
            "-x",
            str(output / "minimization.nc"),
            "-inf",
            str(output / "mdinfo"),
        ]
        plans.append(
            dict(
                state=state,
                gpu=str(i % 2),
                cwd=str(output),
                command=command,
                topology_sha256=accepted["topology_sha256"],
                coordinates_sha256=accepted["coordinates_sha256"],
                mdin_sha256=digest(mdin),
            )
        )
    write_json(
        root / "orchestration_plan.json",
        dict(
            states=plans,
            frozen_repo=str(frozen),
            source_sha256=source_hashes,
            engine_sha256=digest(engine),
            created_utc=datetime.now(timezone.utc).isoformat(),
            policy="All three serial unrestrained minimizations accepted before any MD; no fallback",
        ),
    )


def execute(root: Path) -> None:
    """Run three gated minimizations, then plan/launch existing two-GPU queues."""
    plan = json.loads((root / "orchestration_plan.json").read_text())
    frozen = Path(plan["frozen_repo"])
    for path, sha in plan["source_sha256"].items():
        if digest(Path(path)) != sha:
            raise ValueError(f"Frozen source changed: {path}")
    datasets = {}
    for row in plan["states"]:
        state = row["state"]
        folder = root / state
        output = Path(row["cwd"])
        if (output / "process.json").exists():
            raise FileExistsError(
                "Existing minimization attempt; explicit reviewed continuation required"
            )
        for name, key in [
            ("topology.prmtop", "topology_sha256"),
            ("start.rst7", "coordinates_sha256"),
        ]:
            if digest(folder / name) != row[key]:
                raise ValueError(f"Prepared source changed: {state}")
        if (
            digest(output / "minimize.mdin") != row["mdin_sha256"]
            or digest(Path(row["command"][0])) != plan["engine_sha256"]
        ):
            raise ValueError("Minimization source changed")
        gpu = require_gpu(row["gpu"])
        env = os.environ.copy()
        env.update(
            CUDA_VISIBLE_DEVICES=row["gpu"],
            AMBERHOME="/data/jpc/code/ambertools26/venv",
            OMP_NUM_THREADS="1",
            LD_LIBRARY_PATH="/data/jpc/code/ambertools26/venv/lib:"
            + env.get("LD_LIBRARY_PATH", ""),
        )
        process = dict(
            command=row["command"],
            gpu=gpu,
            started_utc=datetime.now(timezone.utc).isoformat(),
        )
        write_json(output / "process.json", process)
        with (output / "launcher.log").open("x") as log:
            result = subprocess.run(
                row["command"],
                cwd=output,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        write_json(
            output / "process_outcome.json",
            dict(
                returncode=result.returncode,
                ended_utc=datetime.now(timezone.utc).isoformat(),
            ),
        )
        try:
            acceptance = native_acceptance(
                (output / "minimized.mdout").read_text(), result.returncode
            )
            acceptance.update(
                coordinate_acceptance(
                    folder / "topology.prmtop", output / "minimized.ncrst", state
                )
            )
        except (ValueError, OSError) as error:
            write_json(
                output / "scientific_rejection.json",
                dict(accepted=False, reason=str(error), exception=type(error).__name__),
            )
            raise
        acceptance.update(
            state=state,
            native_output_sha256=digest(output / "minimized.mdout"),
            mdin_sha256=row["mdin_sha256"],
            engine_sha256=plan["engine_sha256"],
            prepared_coordinates_sha256=row["coordinates_sha256"],
            frozen_source_sha256=plan["source_sha256"],
        )
        receipt = output / "minimization_acceptance.json"
        write_json(receipt, acceptance)

        def artifact(path: Path) -> dict[str, str]:
            return dict(path=str(path), sha256=digest(path))

        datasets[state] = dict(
            microstate="ASP272/ASH274" if state.endswith("ash274") else "ASP272/ASP274",
            terminal_protonation="NH3+" if state.startswith("nh3") else "NH2",
            restart_stage="unrestrained_minimized",
            protected_residues=protected_rows(state),
            topology=artifact(folder / "topology.prmtop"),
            restart=artifact(output / "minimized.ncrst"),
            preparation_acceptance=artifact(folder / "preparation_acceptance.json"),
            minimization_acceptance=artifact(receipt),
        )
        LOGGER.info("Accepted %s RMS=%g", state, acceptance["rms_gradient"])
    write_json(root / "datasets.json", dict(datasets=datasets))
    paths = configs(root)
    queue = root / "queue_config.json"
    write_json(
        queue,
        dict(
            configs=[str(p) for p in paths],
            gpus=["0", "1"],
            output="queue",
            session_prefix="wt_matched_3state",
        ),
    )
    for action in ("plan", "launch"):
        subprocess.run(
            [
                sys.executable,
                str(frozen / "scripts/queue_prepared_md.py"),
                action,
                "--config",
                str(queue),
                "--repo",
                str(frozen),
            ],
            check=True,
        )


def main(argv: list[str] | None = None) -> None:
    """Prepare reviewable inputs, or explicitly execute the frozen preparation."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "execute"))
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--repo", type=Path, default=REPO)
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    if args.action == "prepare":
        prepare(args.root.resolve(), args.repo.resolve())
    else:
        execute(args.root.resolve())


if __name__ == "__main__":
    main()
