#!/usr/bin/env python3
"""Drive Thr233 HG1 directly to its N terminus, then release the bias."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.direct_proton_transfer import (  # noqa: E402
    driven_input,
    release_input,
    restraint_text,
    schedule,
)
from nyl12_covalent.unrestrained_relay import native_converged  # noqa: E402
from nyl12_covalent.coupled_ti1_drive import coupled_geometry as proton_distances  # noqa: E402
from nyl12_covalent.path_energies import (  # noqa: E402
    driven_path_energies,
    write_driven_path_energy_data,
)
from parmed.amber import AmberParm  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def _command(sander: Path, folder: Path, topology: Path, start: Path) -> list[str]:
    folder = folder.resolve()
    return [
        str(sander.resolve()),
        "-O",
        "-i",
        str(folder / "minimize.mdin"),
        "-p",
        str(topology.resolve()),
        "-c",
        str(start.resolve()),
        "-o",
        str(folder / "minimized.mdout"),
        "-r",
        str(folder / "optimized.ncrst"),
        "-x",
        str(folder / "optimization.nc"),
        "-inf",
        str(folder / "mdinfo"),
    ]


def _run(command: list[str], folder: Path, environment: dict[str, str]) -> None:
    record: dict[str, Any] = {"command": command, "started_utc": _now()}
    _write_json(folder / "process.json", record)
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
    _write_json(folder / "process.json", record)
    if result.returncode:
        raise RuntimeError(f"Sander failed in {folder.name}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--start", required=True, type=Path)
    parser.add_argument("--base-mdin", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--final-coordinate", type=float, default=0.75)
    parser.add_argument("--maximum-increment", type=float, default=0.10)
    parser.add_argument("--cpptraj", type=Path)
    return parser


def run(args: argparse.Namespace) -> None:
    """Execute the sequential drive and unbiased endpoint optimization."""
    paths = (args.topology, args.start, args.base_mdin, args.sander)
    if not all(path.is_file() for path in paths):
        raise FileNotFoundError("One or more required source files are missing")
    system = AmberParm(str(args.topology))
    expected = {
        13905: (929, "N"),
        13906: (929, "H1"),
        13907: (929, "H2"),
        13916: (929, "OG1"),
        13917: (929, "HG1"),
        19916: (1314, "C"),
        19917: (1314, "O"),
        19929: (1315, "N"),
    }
    for index, identity in expected.items():
        atom = system.atoms[index - 1]
        if (atom.residue.idx + 1, atom.name) != identity:
            raise ValueError(f"Reactive atom map mismatch at {index}: {atom}")
    if args.cpptraj is not None and not args.cpptraj.is_file():
        raise FileNotFoundError(args.cpptraj)
    args.output.mkdir(parents=True, exist_ok=False)
    windows = args.output / "01_coordinate_drive"
    release = args.output / "02_unrestrained_endpoint"
    windows.mkdir()
    release.mkdir()
    base = args.base_mdin.read_text()
    initial = proton_distances(args.topology, args.start)
    rows = schedule(
        initial["proton_difference_A"],
        args.final_coordinate,
        args.maximum_increment,
    )
    environment = {
        "AMBERHOME": str(args.sander.parent.parent),
        "LD_LIBRARY_PATH": str(args.sander.parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    _write_json(
        args.output / "inputs.json",
        {
            "scientific_intent": "Direct Thr233 HG1 transfer to neutral N terminus",
            "reaction_coordinate": "d(OG1,HG1)-d(N,HG1)",
            "atom_indices_one_based": {"N": 13905, "OG1": 13916, "HG1": 13917},
            "initial_distances_A": initial,
            "schedule": rows,
            "positional_restraints": "none",
            "active_region": "inherited unchanged from accepted 12 A endpoint input",
            "sources": [
                {"path": str(path.resolve()), "sha256": _sha256(path)} for path in paths
            ],
            "environment": environment,
        },
    )
    status: dict[str, Any] = {
        "state": "running_coordinate_drive",
        "completed_windows": [],
        "started_utc": _now(),
    }
    _write_json(args.output / "status.json", status)
    restart = args.start.resolve()
    try:
        for row in rows:
            folder = windows / row["label"]
            folder.mkdir()
            (folder / "minimize.mdin").write_text(driven_input(base))
            (folder / "restraints.RST").write_text(restraint_text(row["target_A"]))
            _run(
                _command(args.sander, folder, args.topology, restart),
                folder,
                environment,
            )
            endpoint = folder / "optimized.ncrst"
            values = proton_distances(args.topology, endpoint)
            accepted = native_converged(folder / "minimized.mdout")
            _write_json(
                folder / "acceptance.json",
                {
                    "native_converged": accepted,
                    "target_A": row["target_A"],
                    "tracking_error_A": abs(
                        values["proton_difference_A"] - row["target_A"]
                    ),
                    "tracking_is_diagnostic_only": True,
                    "distances_A": values,
                    "positional_restraints": "none",
                },
            )
            if values["N_H1_A"] > 1.25 or values["N_H2_A"] > 1.25:
                raise RuntimeError(f"Terminal H1/H2 identity lost in {row['label']}")
            if not accepted:
                raise RuntimeError(f"DL-FIND did not converge in {row['label']}")
            restart = endpoint
            status["completed_windows"].append(row["label"])
            status["latest_distances_A"] = values
            status["updated_utc"] = _now()
            _write_json(args.output / "status.json", status)
            _publish(args)
        status.update(state="running_unrestrained_endpoint", updated_utc=_now())
        _write_json(args.output / "status.json", status)
        (release / "minimize.mdin").write_text(release_input(base))
        _run(
            _command(args.sander, release, args.topology, restart), release, environment
        )
        final = proton_distances(args.topology, release / "optimized.ncrst")
        converged = native_converged(release / "minimized.mdout")
        retained = (
            final["N_HG1_A"] < 1.25
            and final["OG1_HG1_A"] > 1.30
            and final["N_H1_A"] < 1.25
            and final["N_H2_A"] < 1.25
        )
        acceptance = {
            "native_converged": converged,
            "proton_transferred_state_retained": retained,
            "distances_A": final,
            "reaction_coordinate_restraints": "none",
            "positional_restraints": "none",
            "active_region": "same 12 A DL-FIND active shell as the accepted E-S input",
        }
        _write_json(release / "acceptance.json", acceptance)
        status.update(
            state="completed", endpoint_acceptance=acceptance, updated_utc=_now()
        )
    except Exception as error:
        status.update(state="stopped", reason=str(error), updated_utc=_now())
        raise
    finally:
        _write_json(args.output / "status.json", status)
    _publish(args)


def _publish(args: argparse.Namespace) -> None:
    if args.cpptraj is None:
        return
    status = json.loads((args.output / "status.json").read_text())
    rows = driven_path_energies(
        args.output / "01_coordinate_drive", status["completed_windows"]
    )
    data, _ = write_driven_path_energy_data(rows, args.output)
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/plot_reaction_energies.py"),
            str(data),
            str(args.output / "potential_energy.png"),
            "--title",
            "Neutral flipped-model30 direct proton transfer",
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/export_direct_pt_progress.py"),
            "--run",
            str(args.output),
            "--topology",
            str(args.topology),
            "--cpptraj",
            str(args.cpptraj),
            "--output",
            str(args.output / "vmd_progress"),
        ],
        check=True,
    )


def main() -> None:
    """Parse arguments and execute the direct proton-transfer workflow."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    run(_parser().parse_args())


if __name__ == "__main__":
    main()
