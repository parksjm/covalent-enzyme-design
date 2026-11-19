#!/usr/bin/env python3
"""Drive Thr233 OG1 attack from the valid proton-transferred window 24."""

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
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.coupled_ti1_drive import (  # noqa: E402
    attack_only_schedule,
    attack_with_proton_identity_restraints,
    coupled_geometry,
    ti1_like_geometry,
)
from nyl12_covalent.direct_proton_transfer import driven_input, release_input  # noqa: E402
from nyl12_covalent.unrestrained_relay import native_converged  # noqa: E402
from scripts.export_direct_pt_progress import export  # noqa: E402
from scripts.concurrent_path_publication import publish_progress  # noqa: E402

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


def _proton_identity_retained(values: dict[str, float]) -> bool:
    return all(values[f"N_{name}_A"] <= 1.25 for name in ("HG1", "H1", "H2")) and all(
        values[f"OG1_{name}_A"] >= 1.30 for name in ("HG1", "H1", "H2")
    )


def _resume_state(
    output: Path, rows: list[dict[str, Any]]
) -> tuple[dict[str, Any], Path, int]:
    """Promote a converged identity-safe tracking-only stop for continuation."""
    status_path = output / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError(status_path)
    status = json.loads(status_path.read_text())
    if status["state"] != "stopped":
        raise ValueError("Resume requires a stopped corrected-attack calculation")
    completed = list(status["completed_windows"])
    if completed != [row["label"] for row in rows[: len(completed)]]:
        raise ValueError("Completed-window sequence disagrees with the schedule")
    if len(completed) >= len(rows):
        raise ValueError("No corrected-attack windows remain")
    next_folder = output / "windows" / rows[len(completed)]["label"]
    if completed and not next_folder.exists():
        folder = output / "windows" / completed[-1]
        acceptance = json.loads((folder / "acceptance.json").read_text())
        endpoint = folder / "optimized.ncrst"
        if not endpoint.is_file() or not all(
            acceptance[key]
            for key in (
                "accepted",
                "native_converged",
                "proton_identity_retained",
                "scissile_connectivity_retained",
            )
        ):
            raise ValueError("Last completed window is not safe for continuation")
        history = list(status.get("resume_history", []))
        history.append(
            {
                "resumed_utc": _now(),
                "accepted_checkpoint": completed[-1],
                "previous_reason": status.get("reason"),
                "checkpoint_sha256": _sha256(endpoint),
            }
        )
        status.update(state="running", resume_history=history, updated_utc=_now())
        status.pop("reason", None)
        _write_json(status_path, status)
        return status, endpoint, len(completed)
    row = rows[len(completed)]
    folder = output / "windows" / row["label"]
    acceptance = json.loads((folder / "acceptance.json").read_text())
    endpoint = folder / "optimized.ncrst"
    if (
        not endpoint.is_file()
        or not acceptance["native_converged"]
        or not acceptance["proton_identity_retained"]
    ):
        raise ValueError("Stopped window is not safe for in-place continuation")
    revision = {
        "accepted_for_continuation": True,
        "reason": "Target tracking is diagnostic; convergence and proton identity pass",
        "original_acceptance_sha256": _sha256(folder / "acceptance.json"),
        "attack_tracking_error_A": acceptance["attack_tracking_error_A"],
        "recorded_utc": _now(),
    }
    _write_json(folder / "acceptance_revision.json", revision)
    completed.append(row["label"])
    history = list(status.get("resume_history", []))
    history.append(
        {
            "resumed_utc": _now(),
            "promoted_window": row["label"],
            "previous_reason": status.get("reason"),
        }
    )
    status.update(
        state="running",
        completed_windows=completed,
        resume_history=history,
        updated_utc=_now(),
    )
    status.pop("reason", None)
    _write_json(status_path, status)
    return status, endpoint, len(completed)


def run(args: argparse.Namespace) -> None:
    """Run the corrected attack coordinate with explicit proton-identity gates."""
    sources = (
        args.topology,
        args.start,
        args.base_mdin,
        args.sander,
        args.cpptraj,
    )
    if not all(path.is_file() for path in sources):
        raise FileNotFoundError("Required corrected-attack source is missing")
    attack_only = getattr(args, "attack_only", False)
    release_endpoint = getattr(args, "release_endpoint", False)
    if release_endpoint and not attack_only:
        raise ValueError("Endpoint release requires explicit attack-only mode")
    windows = args.output / "windows"
    initial = coupled_geometry(args.topology, args.start)
    if not _proton_identity_retained(initial):
        raise ValueError("Starting window does not retain all three N-H identities")
    rows = attack_only_schedule(
        initial["OG1_carbonyl_C_A"],
        final_attack=args.final_attack,
        maximum_increment=args.maximum_increment,
    )
    base = args.base_mdin.read_text()
    environment = {
        "AMBERHOME": str(args.sander.resolve().parent.parent),
        "LD_LIBRARY_PATH": str(args.sander.resolve().parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    if getattr(args, "resume", False):
        inputs = json.loads((args.output / "inputs.json").read_text())
        if inputs["schedule"] != rows:
            raise ValueError("Resume arguments disagree with the recorded schedule")
        status, restart, next_index = _resume_state(args.output, rows)
    else:
        if args.output.exists() and not (args.output / "preparation.json").is_file():
            raise FileExistsError(args.output)
        args.output.mkdir(parents=True, exist_ok=True)
        windows.mkdir()
        _write_json(
            args.output / "inputs.json",
            {
                "scientific_intent": (
                    "Sequential nucleophilic attack from the recorded proton-transferred "
                    "restart; proton identities are monitored after every window"
                ),
                "starting_restart": str(args.start.resolve()),
                "initial_geometry_A": initial,
                "driven_coordinate": "d(Thr233_OG1, substrate_carbonyl_C)",
                "proton_identity_gates": {
                    "restraints_applied": not attack_only,
                    "retained_acceptance_N_H_max_A": 1.25,
                    "retained_acceptance_OG1_H_min_A": 1.30,
                },
                "schedule": rows,
                "reaction_coordinate_positional_restraints": "none",
                "active_region": "unchanged accepted 12 A DL-FIND shell",
                "endpoint_release": release_endpoint,
                "attack_only": attack_only,
                "sources": [
                    {"path": str(path.resolve()), "sha256": _sha256(path)}
                    for path in sources
                ],
            },
        )
        status = {
            "state": "running",
            "starting_restart": str(args.start.resolve()),
            "completed_windows": [],
            "started_utc": _now(),
        }
        _write_json(args.output / "status.json", status)
        restart = args.start.resolve()
        next_index = 0
    try:
        for row in rows[next_index:]:
            folder = windows / row["label"]
            folder.mkdir()
            (folder / "minimize.mdin").write_text(
                driven_input(base, max_cycles=args.max_cycles)
            )
            (folder / "restraints.RST").write_text(
                _restraints(row["attack_target_A"], attack_only)
            )
            _run_engine(
                _command(args.sander, folder, args.topology, restart),
                folder,
                environment,
            )
            endpoint = folder / "optimized.ncrst"
            values = coupled_geometry(args.topology, endpoint)
            converged = native_converged(folder / "minimized.mdout")
            tracking_error = abs(values["OG1_carbonyl_C_A"] - row["attack_target_A"])
            retained = _proton_identity_retained(values)
            connectivity = not attack_only or values["carbonyl_C_N_A"] <= 1.80
            accepted = converged and retained and connectivity
            acceptance = {
                "accepted": accepted,
                "native_converged": converged,
                "attack_tracking_error_A": tracking_error,
                "tracking_is_diagnostic_only": True,
                "proton_identity_retained": retained,
                "scissile_connectivity_retained": connectivity,
                "geometry_A": values,
                "ti1_like_geometry": ti1_like_geometry(values),
                "positional_restraints": "none",
            }
            _write_json(folder / "acceptance.json", acceptance)
            if not accepted:
                raise RuntimeError(f"Scientific gate failed in {row['label']}")
            restart = endpoint
            status["completed_windows"].append(row["label"])
            status.update(latest_geometry_A=values, updated_utc=_now())
            _write_json(args.output / "status.json", status)
            bundle = _publish(args)
            LOGGER.info("Accepted %s and published %s", row["label"], bundle)
        if release_endpoint:
            restart = _release(args, base, restart, environment)
        status.update(
            state="completed",
            endpoint_restart=str(restart.resolve()),
            endpoint_is_benchmark_seed=not release_endpoint,
            endpoint_is_stationary_point=False,
            updated_utc=_now(),
        )
        _write_json(args.output / "status.json", status)
        _publish(args)
    except Exception as error:
        status.update(state="stopped", reason=str(error), updated_utc=_now())
        raise
    finally:
        _write_json(args.output / "status.json", status)


def _publish(args: argparse.Namespace) -> Path:
    """Publish each finalized state once, including plots when requested."""
    if getattr(args, "plot_python", None):
        return publish_progress(args, ROOT, "attack_only", "windows")
    return export(
        args.output,
        args.topology,
        args.cpptraj,
        args.output / "vmd_progress",
        "attack_only",
    )


def _restraints(target: float, attack_only: bool) -> str:
    records = attack_with_proton_identity_restraints(target)
    return records.splitlines()[0] + "\n" if attack_only else records


def _release(
    args: argparse.Namespace, base: str, start: Path, environment: dict[str, str]
) -> Path:
    folder = args.output / "02_unrestrained_endpoint"
    folder.mkdir()
    (folder / "minimize.mdin").write_text(release_input(base))
    _run_engine(
        _command(args.sander, folder, args.topology, start), folder, environment
    )
    endpoint = folder / "optimized.ncrst"
    values = coupled_geometry(args.topology, endpoint)
    converged = native_converged(folder / "minimized.mdout")
    retained = _proton_identity_retained(values)
    connectivity = values["carbonyl_C_N_A"] <= 1.80
    accepted = converged and retained and connectivity
    _write_json(
        folder / "acceptance.json",
        {
            "accepted": accepted,
            "native_converged": converged,
            "proton_identity_retained": retained,
            "scissile_connectivity_retained": connectivity,
            "geometry_A": values,
            "ti1_like_geometry": ti1_like_geometry(values),
            "reaction_restraints": "none",
            "stationary_point_validated": False,
        },
    )
    if not accepted:
        raise RuntimeError("Scientific gate failed in reaction-unrestrained endpoint")
    return endpoint


def main() -> None:
    """Parse inputs and run the corrected attack branch."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--start", required=True, type=Path)
    parser.add_argument("--base-mdin", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--final-attack", type=float, default=1.50)
    parser.add_argument("--maximum-increment", type=float, default=0.075)
    parser.add_argument("--max-cycles", type=int, default=2000)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--attack-only", action="store_true")
    parser.add_argument("--release-endpoint", action="store_true")
    parser.add_argument("--plot-python", type=Path)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
