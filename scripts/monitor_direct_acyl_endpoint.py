"""Monitor direct endpoint optimization and export immutable VMD snapshots."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from parmed.amber import Rst7

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from export_direct_pt_progress import (  # noqa: E402
    STRIP_MASK,
    run_cpptraj,
    sha256,
    vmd_script,
    write_checksum_manifest,
)
from nyl12_covalent.path_inspection import dcd_counts  # noqa: E402
from nyl12_covalent.direct_acyl_drive import direct_geometry, candidate_geometry  # noqa: E402
from nyl12_covalent.unrestrained_relay import native_converged  # noqa: E402

LOGGER = logging.getLogger(__name__)
PAIRS = {
    "OG1_C": (13916, 19916),
    "C_O": (19916, 19917),
    "C_N": (19916, 19929),
    "Nsub_HG1": (19929, 13917),
    "OG1_HG1": (13916, 13917),
    "ThrN_HG1": (13905, 13917),
    "ThrN_H1": (13905, 13906),
    "ThrN_H2": (13905, 13907),
    "Nsub_HN": (19929, 19937),
}


def journal(phase: str, action: str) -> None:
    """Record a monitor action in the project append-only audit trail."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "acyl_endpoint_monitor",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Observe authorized reaction-unrestrained endpoint",
            "--rollback",
            "Retain immutable snapshots; original calculation untouched",
        ],
        check=True,
    )


def snapshot(base: Path) -> None:
    """Export starting seed, saved optimizer frames and latest restart."""
    journal(
        "intent",
        "Snapshot endpoint restart/trajectory; evaluate chemistry and export VMD",
    )
    work = base / "unrestrained_optimization_attempt2"
    generation = (
        base
        / "visualization"
        / datetime.now(timezone.utc).strftime("checkpoint_%Y%m%dT%H%M%S")
    )
    generation.mkdir(parents=True)
    source_paths = [
        work / "optimized.ncrst",
        work / "trajectory.nc",
        work / "minimize.mdout",
    ]
    for path in source_paths:
        shutil.copy2(path, generation / path.name)
    xyz = Rst7.open(str(generation / "optimized.ncrst")).coordinates.reshape(-1, 3)
    metrics = {
        name: float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))
        for name, (a, b) in PAIRS.items()
    }
    full_geometry = direct_geometry(
        base / "inputs/system.prmtop", generation / "optimized.ncrst"
    )
    complete = json.loads((base / "job_attempt2/status.json").read_text())
    report = {
        "distances_A": metrics,
        "process_exit_status": complete["exit_status"],
        "native_convergence_reported": native_converged(generation / "minimize.mdout"),
        "chemical_retention_gate": candidate_geometry(full_geometry),
        "full_geometry_A": full_geometry,
        "limitations": "Chemistry gate is independent of convergence; outer shell frozen; no Hessian or NEB; source topology is not a classical product topology",
        "snapshot_sha256": {p.name: sha256(generation / p.name) for p in source_paths},
    }
    cpptraj = Path("/data/jpc/code/ambertools26/venv/bin/cpptraj")
    top = (base / "inputs/system.prmtop").resolve()
    export_input = generation / "trajectory.cpptraj"
    export_input.write_text(
        f"parm {top}\ntrajin {(base / 'inputs/acylenzyme_seed.ncrst').resolve()}\n"
        f"trajin {(generation / 'trajectory.nc').resolve()}\n"
        f"trajin {(generation / 'optimized.ncrst').resolve()}\n"
        f"strip {STRIP_MASK}\ntrajout {(generation / 'reaction_path.dcd').resolve()} dcd\nrun\nquit\n"
    )
    run_cpptraj(cpptraj, export_input, generation, generation / "trajectory_export.log")
    top_input = generation / "topology.cpptraj"
    top_input.write_text(
        f"parm {top}\nparmstrip {STRIP_MASK}\nparmwrite out {(generation / 'stripped.parm7').resolve()}\nquit\n"
    )
    run_cpptraj(cpptraj, top_input, generation, generation / "topology_export.log")
    frames, atoms = dcd_counts(generation / "reaction_path.dcd")
    report.update(display_frames=frames, display_atoms=atoms)
    viewer = vmd_script(
        frames, description="acylenzyme seed, saved unbiased optimizer, latest restart"
    )
    viewer += (
        "\n# Additional endpoint chemistry distances: C=O, C-N, original N-H bonds.\n"
    )
    for name in ("C_O", "C_N", "ThrN_H1", "ThrN_H2", "Nsub_HN"):
        a, b = PAIRS[name]
        viewer += f"label add Bonds $molid/{a - 1} $molid/{b - 1}\n"
    (generation / "reaction_path.vmd").write_text(viewer)
    (generation / "scientific_checkpoint.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    write_checksum_manifest(generation)
    latest = base / "visualization/latest"
    pending = base / "visualization/.latest.pending"
    pending.symlink_to(generation.name)
    pending.replace(latest)
    journal(
        "outcome",
        f"Exported {frames} frames to {generation}; chemistry_gate={report['chemical_retention_gate']}; native_convergence={report['native_convergence_reported']}",
    )
    LOGGER.info("Endpoint checkpoint: %s", json.dumps(report))


def main() -> None:
    """Run one snapshot or monitor until process exit, with a 24-hour bound."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    for _ in range(1440 if args.watch else 1):
        snapshot(args.base.resolve())
        status = json.loads((args.base / "job_attempt2/status.json").read_text())
        if not args.watch or status["exit_status"] is not None:
            break
        time.sleep(60)


if __name__ == "__main__":
    main()
