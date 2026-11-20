"""Native stage execution and uncensored proton observables for NEB umbrellas."""

from __future__ import annotations

from datetime import datetime, timezone
import math
import os
import subprocess
from pathlib import Path
import re
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from .neb_umbrella_campaign import digest, frozen, write_json
from .neb_umbrella_inputs import STAGE_STEPS


def observables(
    xyz: np.ndarray, topology: Any, roles: dict[str, int]
) -> dict[str, Any]:
    """Retain mapped bond/proton distances and unexpected-parent flags without censoring.

    Args:
        xyz: Finite complete coordinates in Angstrom.
        topology: Source topology for element identities only, not reacting bond inference.
        roles: Audited one-based old-band atom identities.

    Returns:
        Both CVs, reactive/Nalpha distances, all three tracked H owners and branch flags.
    """
    r = {k: int(v) - 1 for k, v in roles.items()}

    def distance(a: str, b: str) -> float:
        return float(np.linalg.norm(xyz[r[a]] - xyz[r[b]]))

    pairs = {
        "C_Nleave": ("carbon", "leaving_n"),
        "C_OG1": ("carbon", "thr_og"),
        "OG1_HG1": ("thr_og", "proton"),
        "Nleave_HG1": ("leaving_n", "proton"),
        "Nalpha_HG1": ("thr_n", "proton"),
        "Nalpha_H1": ("thr_n", "nalpha_h1"),
        "Nalpha_H2": ("thr_n", "nalpha_h2"),
    }
    distances = {name: distance(a, b) for name, (a, b) in pairs.items()}
    heavy = np.array([a.idx for a in topology.atoms if a.atomic_number > 1])
    owners = []
    unexpected = False
    for name, allowed in [
        ("proton", {r["thr_og"], r["leaving_n"]}),
        ("nalpha_h1", {r["thr_n"]}),
        ("nalpha_h2", {r["thr_n"]}),
    ]:
        lengths = np.linalg.norm(xyz[heavy] - xyz[r[name]], axis=1)
        order = np.argsort(lengths)[:2]
        nearest = int(heavy[order[0]])
        minimum = float(lengths[order[0]])
        switch = nearest not in allowed and minimum <= 1.45
        unexpected |= switch
        owners.append(
            {
                "identity": name,
                "atom_id_1based": r[name] + 1,
                "nearest_heavy_id_1based": nearest + 1,
                "distance_A": minimum,
                "second_distance_A": float(lengths[order[1]]),
                "shared": float(lengths[order[1]]) <= 1.35,
                "unassigned": minimum > 1.45,
                "unexpected_parent": switch,
            }
        )
    return {
        "distances_A": distances,
        "cv_values": [
            distances["C_Nleave"] - distances["C_OG1"],
            distances["OG1_HG1"] - distances["Nleave_HG1"],
        ],
        "tracked_protons": owners,
        "unexpected_proton_switch_flag": unexpected,
        "direct_route_interpretation_blocked": unexpected,
        "branch_analysis_required": unexpected,
        "sample_retention": "AllfinitevalidMDretained; shared/unassigned/changedprotonsarenotcensored",
    }


def stage_acceptance(stage: Path, kind: str, gate: list[float]) -> dict[str, Any]:
    """Require finite completed native MD and temperature; never discard by proton pathway."""
    text = (stage / "run.mdout").read_text()
    if re.search(
        r"\b(?:NaN|Infinity)\b|SHAKE cannot be accomplished|FATAL|ERROR:", text, re.I
    ):
        raise ValueError("Native numerical/engine failure")
    if re.search(
        r"(?:SCC|SCF)[^\n]*(?:failed|failure|not converg|did not converge)", text, re.I
    ):
        raise ValueError("Native SCC/SCF failure")
    rows = re.findall(
        r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([-+\d.EeDd]+)\s+TEMP\(K\)\s*=\s*([-+\d.EeDd]+)",
        text.split("A V E R A G E S", 1)[0],
    )
    if not rows or int(rows[-1][0]) != STAGE_STEPS[kind]:
        raise ValueError("Missing complete stage timestep count")
    temperatures = [float(row[2].replace("D", "E")) for row in rows]
    energies = [
        float(value.replace("D", "E"))
        for value in re.findall(r"EPtot\s*=\s*(\S+)", text)
    ]
    if not energies or not all(math.isfinite(v) for v in [*energies, *temperatures]):
        raise ValueError("Nonfinite/missing native energies or temperatures")
    if not gate[0] <= temperatures[-1] <= gate[1]:
        raise ValueError("Final stage temperature outside explicit gate")
    return {
        "completed_steps": int(rows[-1][0]),
        "final_temperature_K": temperatures[-1],
        "temperature_gate_K": gate,
        "finite_energy_temperature": True,
        "SCC_error_detected": False,
    }


def _run_stage(
    stage: Path, sander: Path, topology: Path, start: Path, environment: dict[str, str]
) -> Path:
    """Run native Sander with short local names to avoid Fortran filename truncation."""
    for name, source in (("topology.prmtop", topology), ("input.ncrst", start)):
        (stage / name).symlink_to(source.resolve(strict=True))
    command = [
        str(sander),
        "-O",
        "-i",
        "run.mdin",
        "-p",
        "topology.prmtop",
        "-c",
        "input.ncrst",
        "-o",
        "run.mdout",
        "-r",
        "restart.ncrst",
        "-x",
        "trajectory.nc",
        "-inf",
        "mdinfo",
    ]
    record = {
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "cwd": str(stage),
        "short_alias_sources": {
            "topology.prmtop": str(topology),
            "input.ncrst": str(start),
        },
    }
    write_json(stage / "process.json", record)
    with (stage / "engine.log").open("w") as stream:
        result = subprocess.run(
            command,
            cwd=stage,
            env={**os.environ, **environment},
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(
        exit_status=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    write_json(stage / "process.json", record)
    if result.returncode:
        raise RuntimeError(f"Sander failed in {stage}")
    return stage / "restart.ncrst"


def run_window(prepared: Path, label: str) -> None:
    """Run one complete window using existing native process helper and checkpoint gates."""
    from .neb_umbrella_diagnostics import block_diagnostics, export_view

    cfg = frozen(prepared)
    output = prepared.parent
    row = next(r for r in cfg["windows"] if r["label"] == label)
    folder = output / label
    status = folder / "window_status.json"
    if status.exists():
        raise FileExistsError("Windowalreadyattempted; no silent retry")
    topology = parmed.load_file(cfg["topology"])
    source = Rst7.open(row["seed"])
    baseline = np.asarray(source.coordinates).reshape(-1, 3).copy()
    active = np.array(list(AmberMask(topology, cfg["active_mask"]).Selected()))
    if len(active) != 4165:
        raise ValueError("Frozen-boundaryactiveatomcountchanged")
    movable = np.zeros(len(baseline), dtype=bool)
    movable[active] = True
    environment = {
        "AMBERHOME": str(Path(cfg["sander"]).parent.parent),
        "LD_LIBRARY_PATH": str(Path(cfg["sander"]).parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    restart = Path(row["seed"])
    completed = []
    try:
        for spec in row["stages"]:
            stage = folder / spec["label"]
            if (stage / "process.json").exists():
                raise FileExistsError("Stagealreadyattempted")
            write_json(
                status,
                {"status": "running", "stage": spec["label"], "completed": completed},
            )
            write_json(
                stage / "input_identity.json",
                {
                    "source_sha256": {
                        str(p): digest(p)
                        for p in [
                            prepared,
                            Path(cfg["topology"]),
                            restart,
                            stage / "run.mdin",
                            stage / "umbrella.RST",
                        ]
                    },
                    "environment": environment,
                    "started_utc": datetime.now(timezone.utc).isoformat(),
                },
            )
            restart = _run_stage(
                stage, Path(cfg["sander"]), Path(cfg["topology"]), restart, environment
            )
            acceptance = stage_acceptance(
                stage,
                spec["kind"],
                cfg["heating_temperature_gate_K"]
                if spec["kind"] == "heat"
                else cfg["temperature_gate_K"],
            )
            checkpoint = Rst7.open(str(restart))
            xyz = np.asarray(checkpoint.coordinates).reshape(-1, 3)
            velocities = np.asarray(checkpoint.velocities)
            if (
                xyz.shape != baseline.shape
                or not np.isfinite(xyz).all()
                or not checkpoint.hasvels
                or not np.isfinite(velocities).all()
            ):
                raise ValueError("Invalid/nonfinite checkpoint coordinates/velocities")
            if not np.allclose(checkpoint.box, source.box, rtol=0, atol=1e-5):
                raise ValueError("Fixed periodic box changed")
            drift = float(np.max(np.abs(xyz[~movable] - baseline[~movable])))
            if drift > 1e-6:
                raise ValueError("Frozen exterior moved beyond native output precision")
            acceptance.update(
                checkpoint_sha256=digest(restart),
                frozen_outer_max_displacement_A=drift,
                observables=observables(xyz, topology, cfg["roles_1based"]),
            )
            write_json(stage / "acceptance.json", acceptance)
            completed.append(spec["label"])
            if spec["kind"] == "production":
                block_diagnostics(cfg, row, stage, topology)
                export_view(cfg, row, folder, topology)
        write_json(
            status,
            {
                "status": "completed",
                "completed": completed,
                "production_ps": 30,
                "sampling_convergence": "Not established by completion",
            },
        )
    except BaseException as error:
        write_json(
            status,
            {
                "status": "failed",
                "completed": completed,
                "error": f"{type(error).__name__}: {error}",
            },
        )
        raise
