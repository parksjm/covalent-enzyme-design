"""Native Amber execution and explicit acceptance for the window033 relay."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
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
from .nterminal_relay_coordinates import geometry
from .stationary_native import NativeState, digest, write_json

LOGGER = logging.getLogger(__name__)


def action(phase: str, message: str, files: tuple[Path, ...] = ()) -> None:
    """Record intent/outcome in the mandatory project action journal."""
    repo = Path(os.environ["NYL12_REPOSITORY"])
    command = [
        sys.executable,
        str(repo / "scripts/record_agent_action.py"),
        "--log",
        str(repo / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "nterminal_relay_path",
        "--phase",
        phase,
        "--action",
        message,
        "--permission-mode",
        "escalated",
        "--reason",
        "PI conditional relay pathway/TS authorization after window033 minimum gate",
        "--rollback",
        "Preserve native attempts; stop only verified new relay jobs",
    ]
    for path in files:
        command.extend(("--file", str(path)))
    subprocess.run(command, check=True)


def native_energies(text: str, biased: bool) -> dict[str, float]:
    """Read final native total/bias energies without inventing a missing bias.

    Args:
        text: Complete native output.
        biased: Whether the input explicitly contains the reaction bias.

    Returns:
        Total, applied bias and bias-removed potential energies in kcal/mol.
    """

    def last(pattern: str) -> float:
        values = re.findall(pattern, text)
        if not values:
            raise ValueError(f"Missing native energy record: {pattern}")
        value = float(values[-1].replace("D", "E").replace("d", "e"))
        if not np.isfinite(value):
            raise ValueError("Nonfinite native energy")
        return value

    total = last(r"Etot\s*=\s*([-+\d.EeDd]+)")
    bias = last(r"RESTRAINT\s*=\s*([-+\d.EeDd]+)") if biased else 0.0
    if bias < -1e-6:
        raise ValueError("Negative harmonic driving energy")
    return dict(
        total_kcal_mol=total, bias_kcal_mol=bias, physical_kcal_mol=total - bias
    )


def unchanged_exterior(
    initial: np.ndarray, final: np.ndarray, active: np.ndarray
) -> float:
    """Require finite matching coordinates and unchanged original fixed exterior."""
    if final.shape != initial.shape or not np.isfinite(final).all():
        raise ValueError("Invalid native restart shape/coordinates")
    maximum = float(np.max(np.abs(final[~active] - initial[~active])))
    if maximum > 1e-10:
        raise ValueError("Native operation changed the original fixed exterior")
    return maximum


def hydrogen_audit(model: Path, xyz: np.ndarray) -> dict[str, Any]:
    """Monitor other QM protons without restraining the three relay identities."""
    parm = AmberParm(str(model / "inputs/system.prmtop"))
    source = (model / "inputs/dimer.mdin").read_text()
    mask = re.findall(r"\bqmmask\s*=\s*'([^']+)'", source, re.I)
    if len(mask) != 1:
        raise ValueError("Missing unique QM mask")
    qm = AmberMask(parm, mask[0]).Selection()
    excluded = {13906, 13907, 13917, 19937}
    distances = []
    for atom, selected in zip(parm.atoms, qm, strict=True):
        if not selected or atom.atomic_number != 1 or atom.idx + 1 in excluded:
            continue
        if len(atom.bond_partners) != 1:
            raise ValueError("Expected one original parent per non-relay QM proton")
        parent = atom.bond_partners[0]
        distance = float(np.linalg.norm(xyz[atom.idx] - xyz[parent.idx]))
        distances.append(
            dict(atom=atom.idx + 1, parent=parent.idx + 1, distance_A=distance)
        )
    if not distances:
        raise ValueError("No monitored non-relay QM hydrogens")
    passed = all(0.65 <= row["distance_A"] <= 1.50 for row in distances)
    return dict(
        passed=passed,
        original_parent_distances=distances,
        interpretation="Monitoring only; no proton-identity restraint applied",
    )


def run_native(
    campaign: Path,
    folder: Path,
    start: Path,
    mdin: str,
    restraints: str | None = None,
) -> dict[str, Any]:
    """Run one immutable native stage and retain separate scientific gates.

    Args:
        campaign: Prepared and minimum-gated campaign.
        folder: Fresh native stage directory.
        start: Explicit preceding accepted full-system restart.
        mdin: Complete native input retaining the exact model masks.
        restraints: The declared composite reaction potential, or None.

    Returns:
        Numeric/electronic/exterior acceptance and measured chemistry/energies.
    """
    config = json.loads((campaign / "config.json").read_text())
    model = campaign / config["model"]
    state = NativeState(model)
    if digest(state.topology) != config["source_topology_sha256"]:
        raise ValueError("Frozen topology differs from declared source")
    for name in ("active", "qmmask", "qm_theory"):
        original = re.findall(rf"\b{name}\s*=\s*'([^']+)'", state.mdin_source, re.I)
        requested = re.findall(rf"\b{name}\s*=\s*'([^']+)'", mdin, re.I)
        if len(original) != 1 or requested != original:
            raise ValueError(f"Native input changed the original {name}")
    for name in (
        "qmcharge",
        "scfconv",
        "dftb_telec",
        "qmcut",
        "qm_ewald",
        "qm_pme",
        "qmmm_int",
        "adjust_q",
        "lnk_dis",
        "lnk_atomic_no",
        "lnk_method",
        "tight_p_conv",
        "ntb",
        "igb",
        "ntc",
        "ntf",
        "cut",
        "ifqnt",
        "jfastw",
    ):
        pattern = rf"\b{name}\s*=\s*([-+\d.EeDd]+)"
        original = re.findall(pattern, state.mdin_source, re.I)
        requested = re.findall(pattern, mdin, re.I)
        if len(original) > 1 or requested != original:
            raise ValueError(f"Native input changed original QM control {name}")
    if re.findall(r"\bntr\s*=\s*(\d+)", mdin, re.I) != ["0"]:
        raise ValueError("No positional restraint may be introduced")
    if re.findall(r"\bnmropt\s*=\s*(\d+)", mdin, re.I) != [
        "1" if restraints is not None else "0"
    ]:
        raise ValueError("Reaction-bias declaration/input mismatch")
    if not np.array_equal(
        np.asarray(Rst7.open(str(start)).coordinates).reshape((-1, 3))[~state.active],
        state.xyz[~state.active],
    ):
        raise ValueError("Start uses a different fixed environment")
    action("intent", f"Run native relay stage {folder}")
    folder.mkdir(parents=True, exist_ok=False)
    input_path = folder / "minimize.mdin"
    input_path.write_text(mdin)
    sources = [
        input_path,
        state.topology,
        start,
        state.sander,
        model / "library_build/libdlfind.so",
    ]
    if restraints is not None:
        path = folder / "restraints.RST"
        path.write_text(restraints)
        sources.append(path)
    command = [
        str(state.sander),
        "-O",
        "-i",
        str(input_path),
        "-p",
        str(state.topology),
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
    record = dict(
        command=command,
        started_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256={str(p): digest(p) for p in sources},
    )
    write_json(folder / "process.json", record)
    env = dict(
        os.environ,
        AMBERHOME="/data/jpc/code/ambertools26/venv",
        LD_LIBRARY_PATH=f"{model / 'library_build'}:/data/jpc/code/ambertools26/venv/lib",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    with (folder / "engine.log").open("x") as stream:
        result = subprocess.run(
            command,
            cwd=folder,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(
        exit_status=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    write_json(folder / "process.json", record)
    if result.returncode:
        action(
            "outcome", f"Native stage failed with exit {result.returncode}: {folder}"
        )
        raise RuntimeError(f"Native engine failed: {folder}")
    text = (folder / "minimized.mdout").read_text()
    electronic = audit_dftb_scc(text)
    if not electronic["passed"]:
        write_json(folder / "electronic_failure.json", electronic)
        raise RuntimeError("Unrecovered native SCC error; preserve this stage")
    xyz = np.asarray(Rst7.open(str(folder / "optimized.ncrst")).coordinates).reshape(
        (-1, 3)
    )
    exterior = unchanged_exterior(state.xyz, xyz, state.active)
    hydrogens = hydrogen_audit(model, xyz)
    converged = "Optimisation converged" in text
    accepted = converged and hydrogens["passed"]
    response = dict(
        restart=str(folder / "optimized.ncrst"),
        native_converged=converged,
        electronic_history=electronic,
        geometry=geometry(xyz),
        accepted=accepted,
        frozen_max_abs_change_A=exterior,
        nonrelay_hydrogen_audit=hydrogens,
        **native_energies(text, restraints is not None),
    )
    if (
        accepted
        and restraints is None
        and not re.search(r"dimer\s*=\s*\.true\.", mdin, re.I)
    ):
        energy, gradient = state.gradient(
            xyz, folder / "independent_gradient", precision=True
        )
        maximum = float(np.max(np.abs(gradient)) * 0.529177)
        response.update(
            max_gradient_hartree_per_bohr=maximum,
            independent_energy_hartree=energy,
            accepted=maximum <= 1e-5,
        )
    write_json(folder / "acceptance.json", response)
    action(
        "outcome",
        f"Native relay stage {folder}; accepted={response['accepted']}",
        (folder / "acceptance.json",),
    )
    LOGGER.info("Native %s accepted=%s", folder, response["accepted"])
    return response
