"""Immutable inputs and conservative screens for conditional relay minima."""

from __future__ import annotations

from itertools import product
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import numpy as np
from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask
import re

from nyl12_covalent.stationary_native import NativeState, digest, write_json

NAMES = ("ti1_window032", "ti1_window033", "early_window028")
AMBER = Path("/data/jpc/code/ambertools26/venv")


def journal(phase: str, action: str, files: tuple[Path, ...] = ()) -> None:
    """Append a mandatory action record, failing before dependent mutations.

    Args:
        phase: Intent or outcome.
        action: Exact bounded action or result.
        files: Existing affected artifacts to hash.
    """
    repository = Path(os.environ["NYL12_REPOSITORY"])
    command = [
        sys.executable,
        str(repository / "scripts/record_agent_action.py"),
        "--log",
        str(repository / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "relay_minimum_controller",
        "--phase",
        phase,
        "--action",
        action,
        "--permission-mode",
        "escalated",
        "--reason",
        "Authorized conditional full-shell local minimum verification",
        "--rollback",
        "Preserve all stages and original optimizers; stop only owned controller",
    ]
    for path in files:
        command.extend(("--file", str(path)))
    subprocess.run(command, check=True)


def geometry_screen(geometry: dict[str, Any], margin: float = 0.001) -> dict[str, Any]:
    """Screen tetrahedral connectivity without requiring a proton assignment.

    Args:
        geometry: Full geometry report including proton alternatives and angles.
        margin: Angstrom band around any screening boundary needing manual review.

    Returns:
        Retained, outside_screen or ambiguous_boundary; never a chemical proof.
    """
    if not np.isfinite(margin) or margin < 0:
        raise ValueError("Nonnegative finite ambiguity margin required")
    distances = geometry["distances_A"]
    gaps = np.asarray(
        [1.8 - distances["C_OG"], 1.7 - distances["C_NL"], distances["C_O"] - 1.28]
    )
    if not np.isfinite(gaps).all():
        raise ValueError("Nonfinite geometry screen")
    decision = "outside_screen" if np.any(gaps < -margin) else "retained"
    if np.any(np.abs(gaps) <= margin):
        decision = "ambiguous_boundary"
    return {
        "decision": decision,
        "geometry": geometry,
        "ambiguity_margin_A": margin,
        "thresholds_A": {"C_OG_max": 1.8, "C_NL_max": 1.7, "C_O_min": 1.28},
        "meaning": "Declared connectivity screen only; proton transfer unrestricted",
    }


def job_records(campaign: Path) -> dict[str, Any]:
    """Read independent launcher and scientific records for all original jobs.

    Args:
        campaign: Existing unchanged three-seed campaign.

    Returns:
        Records with missing scientific acceptance explicit as null.
    """
    records = {}
    for name in NAMES:
        seed = campaign / name
        acceptance = seed / "scientific_acceptance.json"
        process = json.loads((seed / "job/status.json").read_text())
        records[name] = {
            "process": process,
            "scientific": json.loads(acceptance.read_text())
            if process["exit_status"] is not None and acceptance.exists()
            else None,
        }
    return records


def prerequisites(campaign: Path, reference: Path) -> dict[str, str]:
    """Verify exact method, parameter, input and serialization-library identity.

    Args:
        campaign: Existing relay release directory.
        reference: Validated attempt05 root supplying the serialization library.

    Returns:
        Absolute prerequisite file paths mapped to immutable SHA-256 digests.
    """
    if os.environ["AMBERHOME"] != str(AMBER):
        raise ValueError("Require explicit original AMBERHOME")
    if not os.access(AMBER / "bin/sander", os.X_OK):
        raise FileNotFoundError("Original native sander executable unavailable")
    prepared = json.loads((campaign / "prepared.json").read_text())
    parameters = prepared["parameter_sha256"]
    expected_names = {f"{a}-{b}.skf" for a, b in product("HCNO", repeat=2)}
    if set(parameters) != expected_names:
        raise ValueError("Require all sixteen original 3ob parameter pairs")
    files = [
        campaign / "prepared.json",
        reference.parent / "build_receipt.json",
        reference / "library_build/libdlfind.so",
        AMBER / "bin/sander",
    ]
    receipt = json.loads(files[1].read_text())
    if receipt["change"] != "write_xyz format only: f15.7 to es25.16":
        raise ValueError("Unqualified derivative library")
    if digest(files[2]) != receipt["library_sha256"]:
        raise ValueError("Serialization library differs from validated build")
    seed = campaign / "ti1_window033"
    identity = json.loads((seed / "source_identity.json").read_text())
    for key, value in {
        "active_atom_count": 4165,
        "qm_atom_count": 98,
        "qm_mm_link_count": 8,
        "qm_charge": -2,
    }.items():
        if identity[key] != value:
            raise ValueError(f"Original {key} changed")
    files.append(seed / "source_identity.json")
    for name in ("system.prmtop", "dimer.mdin", "scientific_acceptance.json"):
        path = seed / "inputs" / name
        if digest(path) != digest(reference / "inputs" / name):
            raise ValueError(f"Original reference input differs: {name}")
        files.append(path)
    for name, expected in parameters.items():
        for path in (
            AMBER / "dat/slko/3ob-3-1" / name,
            campaign / "parameters_3ob" / name,
        ):
            if digest(path) != expected:
                raise ValueError(f"Original parameter changed: {path}")
            files.append(path)
    return {str(path): digest(path) for path in files}


def freeze(seed: Path, reference: Path, output: Path) -> Path:
    """Freeze an accepted restart and exact original inputs without overwriting.

    Args:
        seed: Force-accepted w033 root.
        reference: Original validated precision library reference.
        output: One minimum_verification directory.

    Returns:
        Frozen NativeState-compatible root.
    """
    frozen = output / "01_frozen"
    frozen.mkdir()
    (frozen / "inputs").mkdir()
    (frozen / "library_build").mkdir()
    copies = [
        (seed / "inputs" / n, frozen / "inputs" / n)
        for n in ("system.prmtop", "dimer.mdin", "scientific_acceptance.json")
    ]
    copies += [
        (seed / "optimized.ncrst", frozen / "inputs/optimized.ncrst"),
        (
            reference / "library_build/libdlfind.so",
            frozen / "library_build/libdlfind.so",
        ),
        (
            reference.parent / "build_receipt.json",
            frozen / "library_build/build_receipt.json",
        ),
        (seed / "scientific_acceptance.json", frozen / "release_acceptance.json"),
        (seed / "source_identity.json", frozen / "source_identity.json"),
    ]
    for source, destination in copies:
        shutil.copy2(source, destination)
        if digest(source) != digest(destination):
            raise RuntimeError("Frozen copy digest mismatch")
    write_json(
        frozen / "manifest.json",
        {
            str(destination.relative_to(frozen)): digest(destination)
            for _, destination in copies
        },
    )
    return frozen


def verify_frozen(frozen: Path) -> None:
    """Reject changed frozen inputs before a new derivative stage."""
    for name, expected in json.loads((frozen / "manifest.json").read_text()).items():
        if digest(frozen / name) != expected:
            raise ValueError(f"Frozen input changed: {name}")


def probe_columns(state: NativeState) -> list[int]:
    """Select reactive C, neighboring Thr930 HB and first active MM xyz columns."""
    parm = AmberParm(str(state.topology))
    qm_mask = re.findall(r"qmmask='([^']+)'", state.mdin_source)
    if len(qm_mask) != 1:
        raise ValueError("Require one unchanged QM mask")
    qm = np.asarray(AmberMask(parm, qm_mask[0]).Selection(), dtype=bool)
    mm = np.flatnonzero(state.active & ~qm)
    if not len(mm):
        raise ValueError("No representative movable MM atom")
    atoms = [19916 - 1, 13925 - 1, int(mm[0])]
    if parm.atoms[13925 - 1].name != "HB":
        raise ValueError("Thr930 HB atom mapping changed")
    columns = []
    for atom in atoms:
        locations = np.flatnonzero(state.indices == atom)
        if len(locations) != 1:
            raise ValueError("Required probe atom outside unchanged active shell")
        columns.extend(int(locations[0]) * 3 + axis for axis in range(3))
    return sorted(set(columns))


def proton_inventory(state: NativeState) -> dict[str, Any]:
    """Report every QM proton's nearest heavy atoms and possible bond partners."""
    parm = AmberParm(str(state.topology))
    mask = re.findall(r"qmmask='([^']+)'", state.mdin_source)[0]
    qm = np.asarray(AmberMask(parm, mask).Selection(), dtype=bool)
    heavy = np.asarray([a.idx for a in parm.atoms if a.atomic_number > 1])
    hydrogens = [a for a in parm.atoms if qm[a.idx] and a.atomic_number == 1]
    inventory = []
    for atom in hydrogens:
        distances = np.linalg.norm(state.xyz[heavy] - state.xyz[atom.idx], axis=1)
        nearest = np.argsort(distances)[:3]
        inventory.append(
            {
                "hydrogen_id_1based": atom.idx + 1,
                "name": atom.name,
                "topology_partners_1based": [a.idx + 1 for a in atom.bond_partners],
                "nearest_heavy": [
                    {
                        "atom_id_1based": int(heavy[i]) + 1,
                        "distance_A": float(distances[i]),
                    }
                    for i in nearest
                ],
                "partners_within_1p2A": (heavy[distances < 1.2] + 1).tolist(),
                "ambiguous_cutoff": bool(np.any(np.abs(distances - 1.2) <= 0.02)),
            }
        )
    return {
        "qm_hydrogen_count": len(hydrogens),
        "protons": inventory,
        "meaning": "Geometric inventory; no required N-terminal occupancy",
    }
