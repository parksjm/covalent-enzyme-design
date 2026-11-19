#!/usr/bin/env python3
"""Drive direct OG1 proton transfer to substrate N and acyl formation."""

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
sys.path.insert(0, str(ROOT))

from nyl12_covalent.coupled_ti1_drive import (  # noqa: E402
    asynchronous_schedule,
)
from nyl12_covalent.direct_proton_transfer import (  # noqa: E402
    driven_input,
    release_input,
)
from nyl12_covalent.unrestrained_relay import native_converged  # noqa: E402

from nyl12_covalent.direct_acyl_drive import (  # noqa: E402
    direct_geometry as coupled_geometry,
    direct_restraint as asynchronous_restraint,
    composite,
    candidate_geometry,
    chemistry_retained,
)
from scripts.concurrent_path_publication import publish_progress  # noqa: E402
from nyl12_covalent.direct_acyl_resume import prepare_resume  # noqa: E402

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


def _run_engine(command: list[str], folder: Path, environment: dict[str, str]) -> None:
    record: dict[str, Any] = {"command": command, "started_utc": _now()}
    _write_json(folder / "process.json", record)
    with (folder / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=folder,
            env={**os.environ, **environment},
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(exit_status=result.returncode, ended_utc=_now())
    _write_json(folder / "process.json", record)
    if result.returncode:
        raise RuntimeError(f"Sander failed in {folder.name}")


def run(args: argparse.Namespace) -> None:
    """Run the coupled coordinate path and unbiased TI1 endpoint test."""
    sources = (args.topology, args.start, args.base_mdin, args.sander)
    if not all(path.is_file() for path in sources):
        raise FileNotFoundError("Required coupled-path source is missing")
    if not args.resume_after:
        args.output.mkdir(parents=True, exist_ok=False)
    windows = args.output / "01_coupled_coordinate_drive"
    release = args.output / "02_unrestrained_ti1_endpoint"
    if not args.resume_after:
        windows.mkdir()
        release.mkdir()
    base = args.base_mdin.read_text()
    initial = coupled_geometry(args.topology, args.start)
    rows = asynchronous_schedule(
        initial["proton_difference_A"] + initial["carbonyl_C_N_A"],
        initial["OG1_carbonyl_C_A"],
        final_composite=args.final_composite,
        maximum_increment=args.maximum_increment,
    )
    environment = {
        "AMBERHOME": str(args.sander.resolve().parent.parent),
        "LD_LIBRARY_PATH": str(args.sander.resolve().parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    restart = args.start.resolve()
    if args.resume_after:
        rows, status = prepare_resume(args)
    else:
        _write_json(
            args.output / "inputs.json",
            {
                "scientific_intent": "Direct HG1 transfer to substrate N plus ester formation and scissile cleavage",
                "driven_coordinate": "d(OG1,HG1)-d(HG1,substrate_N)-d(OG1,carbonyl_C)+d(carbonyl_C,substrate_N)",
                "monitored_components": [
                    "d(OG1,HG1)-d(N,HG1)",
                    "d(OG1,carbonyl_C)",
                ],
                "atom_indices_one_based": {
                    "Thr_N": 13905,
                    "substrate_N": 19929,
                    "substrate_HN": 19937,
                    "OG1": 13916,
                    "HG1": 13917,
                    "carbonyl_C": 19916,
                },
                "initial_geometry_A": initial,
                "diagnostic_only_not_pmf": True,
                "endpoint_label": "acyl-enzyme-like candidate; no Hessian validation",
                "early_release": "Converged mapped candidate; all reaction bias removed",
                "schedule": rows,
                "positional_restraints": "none",
                "active_region": "unchanged accepted 12 A DL-FIND shell",
                "sources": [
                    {"path": str(path.resolve()), "sha256": _sha256(path)}
                    for path in sources
                ],
            },
        )
        status: dict[str, Any] = {
            "state": "running",
            "completed_windows": [],
            "started_utc": _now(),
        }
        _write_json(args.output / "status.json", status)
    try:
        for row in rows:
            folder = windows / row["label"]
            folder.mkdir()
            (folder / "minimize.mdin").write_text(driven_input(base))
            (folder / "restraints.RST").write_text(
                asynchronous_restraint(
                    row["composite_target_A"],
                    substrate_proton_guard=args.substrate_proton_guard,
                )
            )
            _run_engine(
                _command(args.sander, folder, args.topology, restart),
                folder,
                environment,
            )
            endpoint = folder / "optimized.ncrst"
            values = coupled_geometry(args.topology, endpoint)
            coordinate_value = composite(values)
            converged = native_converged(folder / "minimized.mdout")
            _write_json(
                folder / "acceptance.json",
                {
                    "native_converged": converged,
                    "chemistry_retained": chemistry_retained(values),
                    "geometry_A": values,
                    "composite_coordinate_A": coordinate_value,
                    "composite_tracking_error_A": abs(
                        coordinate_value - row["composite_target_A"]
                    ),
                    "tracking_is_diagnostic_only": True,
                    "positional_restraints": "none",
                },
            )
            status["completed_windows"].append(row["label"])
            status.update(latest_geometry_A=values, updated_utc=_now())
            _write_json(args.output / "status.json", status)
            publish_progress(args, ROOT)
            if not converged or not chemistry_retained(values):
                raise RuntimeError(
                    f"Convergence or mapped-chemistry gate failed in {row['label']}"
                )
            restart = endpoint
            if candidate_geometry(values):
                break
        if not candidate_geometry(values):
            raise RuntimeError(
                "Bounded diagnostic schedule exhausted without acyl-enzyme-like candidate"
            )
        status.update(state="running_unrestrained_ti1_endpoint", updated_utc=_now())
        _write_json(args.output / "status.json", status)
        (release / "minimize.mdin").write_text(release_input(base))
        _run_engine(
            _command(args.sander, release, args.topology, restart), release, environment
        )
        final = coupled_geometry(args.topology, release / "optimized.ncrst")
        converged = native_converged(release / "minimized.mdout")
        retained = converged and candidate_geometry(final)
        acceptance = {
            "native_converged": converged,
            "acyl_candidate_retained": retained,
            "chemistry_retained": chemistry_retained(final),
            "geometry_A": final,
            "reaction_coordinate_restraints": "none",
            "positional_restraints": "none",
        }
        _write_json(release / "acceptance.json", acceptance)
        status.update(
            state="completed", endpoint_acceptance=acceptance, updated_utc=_now()
        )
        _write_json(args.output / "status.json", status)
        publish_progress(args, ROOT)
        if not retained:
            raise RuntimeError(
                "Unbiased endpoint failed convergence or candidate retention"
            )
    except Exception as error:
        status.update(state="stopped", reason=str(error), updated_utc=_now())
        raise
    finally:
        _write_json(args.output / "status.json", status)


def main() -> None:
    """Parse inputs and run the asynchronous E-S-to-TI1 path."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--start", required=True, type=Path)
    parser.add_argument("--base-mdin", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--plot-python", required=True, type=Path)
    parser.add_argument("--resume-after", help="Explicit last accepted existing window")
    parser.add_argument("--checkpoint-sha256", help="Required digest for resumption")
    parser.add_argument("--substrate-proton-guard", action="store_true")
    parser.add_argument("--final-composite", type=float, default=2.5)
    parser.add_argument("--maximum-increment", type=float, default=0.15)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
