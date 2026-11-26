"""Prepare and run a bounded frozen-exterior acyl-enzyme hydration probe."""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from .dftb_scc_validation import audit_dftb_scc
from .neb_umbrella_inputs import active_group_cards, audited_qmmm_block


def digest(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: dict[str, Any]) -> None:
    """Atomically save finite JSON evidence."""
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def build_input(source: str, active: str, stage: str, seed: int) -> str:
    """Build one unbiased NVT heating or equilibration input.

    Args:
        source: Original native input containing the audited QM block.
        active: Exact original active mask, checked against source.
        stage: Either heat or equil.
        seed: Positive signed-32-bit Langevin random seed.

    Returns:
        Complete Sander input with no reaction or positional bias.
    """
    matches = re.findall(r"\bactive\s*=\s*'([^']+)'", source, re.I)
    if matches != [active] or stage not in {"heat", "equil"}:
        raise ValueError("Invalid stage or changed active mask")
    if type(seed) is not int or not 0 < seed <= 2147483647:
        raise ValueError("Invalid Langevin seed")
    heat = stage == "heat"
    steps = 2000 if heat else 10000
    initial = "irest=0, ntx=1, tempi=25.0" if heat else "irest=1, ntx=5"
    ramp = (
        "&wt type='TEMP0', istep1=0, istep2=2000, value1=25.0, value2=300.0 /\n&wt type='END' /\n"
        if heat
        else ""
    )
    return f"""Acyl enzyme hydration only: {stage}
&cntrl
 imin=0, {initial}, nstlim={steps}, dt=0.0005,
 ntb=1, ntp=0, igb=0, cut=10.0, iwrap=0,
 ntt=3, gamma_ln=2.0, temp0=300.0, ig={seed},
 ntc=1, ntf=1, ifqnt=1, jfastw=4, ntr=0, nmropt={int(heat)},
 ibelly=1, bellymask='', ntpr=100, ntwx=100, ntwr=500, ioutfm=1, ntxo=2,
/
{audited_qmmm_block(source)}
{active_group_cards(active)}{ramp}"""


def audit_output(text: str, steps: int) -> dict[str, Any]:
    """Validate native completed steps, finite thermodynamics and SCC history.

    Args:
        text: Complete native output.
        steps: Required exact last reported step.

    Returns:
        Separate numeric and electronic gate evidence.
    """
    electronic = audit_dftb_scc(text)
    if not electronic["verified_scc_history"]:
        raise ValueError("SCC history failed")
    if re.search(
        r"\b(?:NaN|Infinity)\b|SHAKE cannot be accomplished|FATAL|ERROR:", text, re.I
    ):
        raise ValueError("Native numeric or engine failure")
    data = text.split("A V E R A G E S", 1)[0]
    rows = re.findall(
        r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([-+\d.EeDd]+)\s+TEMP\(K\)\s*=\s*([-+\d.EeDd]+)",
        data,
    )
    energies = re.findall(r"EPtot\s*=\s*(\S+)", data)
    if not rows or int(rows[-1][0]) != steps or not energies:
        raise ValueError("Missing exact completed steps or energies")
    values = [
        float(v.replace("D", "E").replace("d", "e")) for row in rows for v in row[1:]
    ] + [float(v.replace("D", "E").replace("d", "e")) for v in energies]
    if not all(math.isfinite(v) for v in values):
        raise ValueError("Nonfinite energy, time or temperature")
    temperature = float(rows[-1][2].replace("D", "E"))
    lower = {2000: 150, 10000: 250}.get(steps)
    if lower is None or not lower <= temperature <= 350:
        raise ValueError("Temperature outside declared range")
    return {
        "completed_steps": steps,
        "final_temperature_K": temperature,
        "electronic_history": electronic,
    }


def audit_restart(
    path: Path, baseline: np.ndarray, active: np.ndarray, box: np.ndarray
) -> tuple[np.ndarray, dict[str, Any]]:
    """Check finite complete coordinates, velocities, fixed exterior and box.

    Args:
        path: Native full-system restart.
        baseline: Frozen original endpoint coordinates.
        active: Boolean original movable atom mask.
        box: Original six periodic box values.

    Returns:
        Coordinates and measured boundary evidence.
    """
    restart = Rst7.open(str(path))
    xyz = np.asarray(restart.coordinates).reshape(-1, 3)
    vel = np.asarray(restart.vels)
    newbox = np.asarray(restart.box)
    if xyz.shape != baseline.shape or vel.size != xyz.size or newbox.shape != (6,):
        raise ValueError("Restart shape or velocities invalid")
    if not all(np.isfinite(v).all() for v in (xyz, vel, newbox)):
        raise ValueError("Nonfinite restart")
    displacement = float(np.max(np.abs(xyz[~active] - baseline[~active])))
    if displacement > 1e-6 or not np.allclose(newbox, box, atol=1e-5, rtol=0):
        raise ValueError("Frozen exterior or periodic box changed")
    return xyz, {"frozen_max_displacement_A": displacement, "box_unchanged": True}


def _journal(repo: Path, phase: str, message: str) -> None:
    subprocess.run(
        [
            sys.executable,
            str(repo / "scripts/record_agent_action.py"),
            "--log",
            str(repo / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "deacylation_hydration",
            "--phase",
            phase,
            "--action",
            message,
            "--permission-mode",
            "escalated",
            "--reason",
            "PI authorized short original-QM hydration reconnaissance",
            "--rollback",
            "Preserve completed stage artifacts and stop only verified pilot process",
        ],
        check=True,
    )


def run(campaign: Path) -> dict[str, Any]:
    """Run the two prepared native stages, preserving separate scientific gates.

    Args:
        campaign: Previously prepared immutable campaign directory.

    Returns:
        Completed status or chemistry-stopped status.
    """
    campaign = campaign.resolve(strict=True)
    cfg = json.loads((campaign / "config.json").read_text())
    repo = Path(cfg["journal_repo"])
    if (
        digest(Path(cfg["audit_path"])) != cfg["audit_sha256"]
        or digest(Path(cfg["sander"]["path"])) != cfg["sander"]["sha256"]
    ):
        raise ValueError("Audit or executable identity changed")
    for name, h in cfg["input_sha256"].items():
        if digest(campaign / name) != h:
            raise ValueError("Prepared input hash changed")
    for name, h in cfg["source_snapshot_sha256"].items():
        if digest(campaign / "source_snapshot" / name) != h:
            raise ValueError("Frozen source hash changed")
    for name, expected in cfg["dftb_parameters_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError("DFTB parameter identity changed")
    top = campaign / "inputs/system.prmtop"
    start = campaign / "inputs/endpoint.ncrst"
    parm = AmberParm(str(top))
    active = np.array(AmberMask(parm, cfg["active_mask"]).Selection(), dtype=bool)
    qm = np.array(AmberMask(parm, cfg["qmmask"]).Selection(), dtype=bool)
    if qm.shape != active.shape or qm.sum() != 98 or not np.all(active[qm]):
        raise ValueError("Original QM/active membership changed")
    links = sum(bool(qm[b.atom1.idx]) != bool(qm[b.atom2.idx]) for b in parm.bonds)
    if links != 8:
        raise ValueError("Original eight link bonds changed")
    initial = Rst7.open(str(start))
    base = np.asarray(initial.coordinates).reshape(-1, 3)
    box = np.asarray(initial.box)
    if (
        base.shape != (104922, 3)
        or active.sum() != 4165
        or box.shape != (6,)
        or not np.isfinite(base).all()
        or not np.isfinite(box).all()
    ):
        raise ValueError("Invalid original endpoint or active boundary")
    status = {"status": "running", "scope": cfg["scope"]}
    write_json(campaign / "status.json", status)
    for stage, steps in [("heat", 2000), ("equil", 10000)]:
        folder = campaign / stage
        _journal(repo, "intent", f"Run hydration {stage} at {folder}")
        folder.mkdir(exist_ok=False)
        cmd = [
            cfg["sander"]["path"],
            "-O",
            "-i",
            str(campaign / "inputs" / f"{stage}.mdin"),
            "-p",
            str(top),
            "-c",
            str(start),
            "-o",
            str(folder / "run.mdout"),
            "-r",
            str(folder / "restart.ncrst"),
            "-x",
            str(folder / "trajectory.nc"),
            "-inf",
            str(folder / "mdinfo"),
        ]
        env = dict(
            os.environ,
            AMBERHOME="/data/jpc/code/ambertools26/venv",
            LD_LIBRARY_PATH="/data/jpc/code/ambertools26/venv/lib",
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
        )
        process = {
            "status": "running",
            "command": cmd,
            "start_sha256": digest(start),
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "env_overrides": {
                key: env[key]
                for key in (
                    "AMBERHOME",
                    "LD_LIBRARY_PATH",
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                )
            },
        }
        write_json(folder / "process.json", process)
        with (folder / "engine.log").open("x") as log:
            result = subprocess.run(
                cmd,
                cwd=folder,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        process.update(
            status="completed" if result.returncode == 0 else "failed",
            exit_code=result.returncode,
            ended_utc=datetime.now(timezone.utc).isoformat(),
        )
        write_json(folder / "process.json", process)
        try:
            if result.returncode:
                raise RuntimeError("Native process failed")
            numeric = audit_output((folder / "run.mdout").read_text(), steps)
            xyz, boundary = audit_restart(folder / "restart.ncrst", base, active, box)

            def d(a: int, b: int) -> float:
                return float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))

            chemical = {
                "C_OG_A": d(19916, 13916),
                "C_Nleave_A": d(19916, 19929),
                "C_O_A": d(19916, 19917),
                "Nalpha_H1_A": d(13905, 13906),
                "Nalpha_H2_A": d(13905, 13907),
                "Nleave_HG1_A": d(19929, 13917),
                "Nleave_HN_A": d(19929, 19937),
            }
            chemical["acyl_bond_retained"] = chemical["C_OG_A"] < 1.8
            write_json(
                folder / "scientific_acceptance.json",
                {
                    "numeric": numeric,
                    "boundary": boundary,
                    "chemistry": chemical,
                    "sampling_converged": False,
                },
            )
            _journal(
                repo,
                "outcome",
                f"Hydration {stage} complete; acyl retention {chemical['acyl_bond_retained']}",
            )
            if not chemical["acyl_bond_retained"]:
                status.update(
                    status="stopped", reason="acyl_bond_not_retained", stage=stage
                )
                write_json(campaign / "status.json", status)
                return status
            start = folder / "restart.ncrst"
        except Exception as exc:
            status.update(
                status="failed", stage=stage, error=f"{type(exc).__name__}: {exc}"
            )
            write_json(campaign / "status.json", status)
            _journal(
                repo,
                "outcome",
                f"Hydration {stage} failed: {type(exc).__name__}: {exc}",
            )
            raise
    status.update(status="completed", sampling_converged=False)
    write_json(campaign / "status.json", status)
    return status
