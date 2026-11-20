"""Unrestrained MM preparation gates for the matched three-state WT campaign."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

STATES = ("nh2_asp272_asp274", "nh2_asp272_ash274", "nh3_asp272_asp274")
MIN_INPUT = """Matched WT fully unrestrained MM minimization
&cntrl
 imin=1, maxcyc=10000, ncyc=1000, drms=0.01,
 ntb=1, cut=10.0, ntr=0, nmropt=0, ntc=1, ntf=1,
 ntpr=100, ntwx=100, ntwr=1000, ioutfm=1, ntxo=2,
/
"""


def digest(path: Path) -> str:
    """Return a streamed SHA-256 identity."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def write_json(path: Path, value: Any) -> None:
    """Write a new immutable JSON artifact, refusing replacement."""
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def protected_rows(state: str) -> list[dict[str, Any]]:
    """Return explicit full-sequence/Amber mapping for all24 protected residues."""
    if state not in STATES:
        raise ValueError("Unreviewed preparation state")
    rows = []
    for i, site in enumerate(("AE", "BF", "CG", "DH")):
        for sequence, amber, name in [
            (113, 113 + 232 * i, "PHE"),
            (156, 156 + 232 * i, "LYS"),
            (186, 186 + 232 * i, "ASN"),
            (233, 929 + 96 * i, "THR"),
            (272, 968 + 96 * i, "ASP"),
            (274, 970 + 96 * i, "ASP"),
        ]:
            rows.append(
                dict(
                    site=site,
                    sequence_residue=sequence,
                    amber_residue=amber,
                    chain=site[0] if sequence < 233 else site[1],
                    resname="ASH"
                    if sequence == 274 and state.endswith("ash274")
                    else name,
                )
            )
    return rows


def native_acceptance(text: str, returncode: int) -> dict[str, Any]:
    """Require successful finite native minimization and RMS<=0.1 for MD prep.

    The requested0.01 force target is recorded separately from the0.1 MD gate.
    """
    if returncode or not re.search(r"TIMINGS", text, re.I):
        raise ValueError("Minimization process incomplete or failed")
    if re.search(r"\b(?:NaN|Infinity)\b|ERROR|FATAL|CUDA error", text, re.I):
        raise ValueError("Native numeric or engine failure")
    rows = re.findall(
        r"NSTEP\s+ENERGY\s+RMS\s+GMAX\s+NAME\s+NUMBER\s*\n\s*"
        r"(\d+)\s+([-+\d.EeDd]+)\s+([-+\d.EeDd]+)\s+([-+\d.EeDd]+)",
        text,
    )
    if not rows:
        raise ValueError("Missing native minimization force summary")
    step, energy, rms, gmax = rows[-1]
    energy, rms, gmax = [
        float(v.replace("D", "E").replace("d", "e")) for v in (energy, rms, gmax)
    ]
    if not all(math.isfinite(v) for v in (energy, rms, gmax)) or rms < 0 or rms > 0.1:
        raise ValueError(f"Minimization RMS {rms} exceeds MD-preparation0.1 gate")
    return dict(
        process_complete=True,
        process_exit_status=returncode,
        final_step=int(step),
        energy_kcal_mol=energy,
        rms_gradient=rms,
        max_gradient=gmax,
        force_converged_target_0_01=rms <= 0.01,
        force_converged_target=0.01,
        force_converged=rms <= 0.01,
        MD_preparation_RMS_limit=0.1,
        MD_preparation_accepted=True,
        restraints="ntr=0 nmropt=0; all atoms movable",
    )


def coordinate_acceptance(topology: Path, restart: Path, state: str) -> dict[str, Any]:
    """Validate finite matching coordinates, periodic box and actual proton bonds."""
    import numpy as np
    from parmed.amber import AmberParm
    from prepared_md_support import validate_protected

    top = AmberParm(str(topology), xyz=str(restart))
    xyz, box = np.asarray(top.coordinates), np.asarray(top.box)
    if (
        xyz.shape != (len(top.atoms), 3)
        or not np.isfinite(xyz).all()
        or box.shape != (6,)
        or not np.isfinite(box).all()
        or np.any(box <= 0)
    ):
        raise ValueError("Invalid topology/restart coordinates or periodic box")
    charged = state.startswith("nh3")
    rows = protected_rows(state)
    validate_protected(top, rows, "NH3+" if charged else "NH2", 1.0 if charged else 0.0)
    observed = []
    for row in rows:
        residue = top.residues[row["amber_residue"] - 1]
        atoms = {a.name: a for a in residue.atoms}
        pairs = (
            [("N", h) for h in ("H1", "H2", "H3") if h in atoms]
            if row["sequence_residue"] == 233
            else []
        )
        if residue.name == "ASH":
            pairs.append(("OD2", "HD2"))
        if row["sequence_residue"] == 233:
            pairs.append(("OG1", "HG1"))
        for heavy, hydrogen in pairs:
            length = float(
                np.linalg.norm(xyz[atoms[heavy].idx] - xyz[atoms[hydrogen].idx])
            )
            if not 0.7 <= length <= 1.3:
                raise ValueError(
                    f"Proton bond geometry changed: {row}/{hydrogen}={length}"
                )
            observed.append(
                dict(site=row["site"], hydrogen=hydrogen, distance_A=length)
            )
    return dict(
        restart_readable=True,
        coordinate_finite=True,
        topology_atom_match=True,
        proton_inventory_accepted=True,
        proton_bonds=observed,
        topology_sha256=digest(topology),
        coordinates_sha256=digest(restart),
    )


def configs(root: Path) -> list[Path]:
    """Create three explicit four-replica configs with72distinct stage seeds."""
    folder = root / "md_configs"
    folder.mkdir()
    paths = []
    for state_index, state in enumerate(STATES):
        path = folder / f"{state}.json"
        replicas = [
            dict(
                name=f"replica_{i + 1:02d}",
                seed=20270001 + state_index * 400000 + i * 60000,
                gpu=str(i % 2),
            )
            for i in range(4)
        ]
        write_json(
            path,
            dict(
                dataset_manifest="../datasets.json",
                dataset=state,
                amber_home="/data/jpc/code/ambertools26/venv",
                cuda_engine="/data/jpc/code/pmemd26/bin/pmemd.cuda",
                output=f"../md/{state}",
                session_prefix=f"wt_{state}",
                replicas=replicas,
            ),
        )
        paths.append(path)
    return paths
