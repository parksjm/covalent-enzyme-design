#!/usr/bin/env python3
"""Audit explicit WT PA6 MD sources and extract reproducible solvated seeds."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

import numpy as np
from parmed.amber import AmberParm
from parmed.amber.mdin import Mdin
from parmed.formats import CIFFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.es_snapshot_selection import select_es_frames  # noqa: E402
from nyl12_covalent.es_topology_audit import audit_es_topology  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _source_record(path: Path) -> dict[str, Any]:
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source changed during hashing: {path}")
    return {
        "path": str(path.resolve()),
        "bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256": digest.hexdigest(),
    }


def _production_record(directory: Path) -> dict[str, Any]:
    control = Mdin()
    control.read(str(directory / "07_prod_100ns.in"))
    values = control.cntrl_nml
    if (
        not math.isfinite(values["dt"])
        or values["dt"] <= 0
        or not isinstance(values["ntwx"], int)
        or values["ntwx"] <= 0
        or not isinstance(values["nstlim"], int)
        or values["nstlim"] <= 0
    ):
        raise ValueError(f"Invalid production duration/output cadence: {directory}")
    text = (directory / "07_prod_100ns.out").read_text()
    steps = [int(n) for n in re.findall(r"NSTEP\s*=\s*(\d+)", text)]
    if "Final Performance Info" not in text or not steps:
        raise ValueError(f"No normal completion evidence: {directory}")
    if max(steps) != values["nstlim"] or values["ntr"] != 0:
        raise ValueError(f"Incomplete or restrained production: {directory}")
    if values["imin"] != 0 or values["irest"] != 1 or values["ntb"] != 2:
        raise ValueError(f"Unexpected production mode: {directory}")
    return {
        "directory": str(directory),
        "last_step": max(steps),
        "declared": {
            k: values[k]
            for k in ("nstlim", "dt", "ntwx", "ig", "ntr", "ntc", "ntf", "ntb")
        },
        "duration_ns": values["nstlim"] * values["dt"] / 1000,
        "frame_dt_ns": values["ntwx"] * values["dt"] / 1000,
        "normal_completion": True,
    }


def _native_mapping(path: Path, rows: list[dict[str, Any]]) -> dict[str, Any]:
    crystal = CIFFile.parse(str(path), skip_bonds=True)
    if crystal.get_coordinates().shape[0] != 1:
        raise ValueError("Native reference must have exactly one model")
    observed, positions = [], {}
    for row in rows:
        matches = [
            r
            for r in crystal.residues
            if r.chain == row["chain"] and r.number == row["number"]
        ]
        if len(matches) != 1 or matches[0].name != row["resname"]:
            raise ValueError(f"Native residue identity mismatch: {row}")
        atoms = {a.name: a for a in matches[0].atoms}
        if len(atoms) != len(matches[0].atoms) or any(
            a.other_locations for a in atoms.values()
        ):
            raise ValueError(
                "Native alternate/duplicate atoms require explicit selection"
            )
        xyz = {name: np.array([a.xx, a.xy, a.xz]) for name, a in atoms.items()}
        if not all(np.isfinite(v).all() for v in xyz.values()):
            raise ValueError("Nonfinite native coordinates")
        positions[row["label"]] = xyz
        observed.append({**row, "native_identity_verified": True})
    distances = {}
    for label in ("asp272", "asp274"):
        for oxygen in ("OD1", "OD2"):
            for target in ("N", "OG1"):
                distances[f"{label}:{oxygen}-thr233:{target}"] = float(
                    np.linalg.norm(
                        positions[label][oxygen] - positions["thr233"][target]
                    )
                )
    return {
        "rows": observed,
        "distances_angstrom": distances,
        "literature_residue_correspondence": "Not established by this identity check",
    }


def _extract_frame(
    system: Path,
    replica: str,
    frame: dict[str, Any],
    output: Path,
    amberhome: Path,
) -> dict[str, Any]:
    target = output / f"{replica}_frame_{frame['full_frame']:04d}"
    target.mkdir()
    trajectory = target / "solvated.nc"
    script = target / "extract.cpptraj"
    script.write_text(
        f"parm {system / 'md_solvated.prmtop'}\n"
        f"trajin {system / replica / '07_prod_100ns.nc'} "
        f"{frame['full_frame']} {frame['full_frame']}\n"
        "autoimage anchor :1313\n"
        f"trajout {trajectory} netcdf\nrun\n"
    )
    env = os.environ.copy()
    env.update(AMBERHOME=str(amberhome), OMP_NUM_THREADS="1")
    env["LD_LIBRARY_PATH"] = str(amberhome / "lib") + (
        ":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else ""
    )
    with (target / "cpptraj.log").open("w") as log:
        subprocess.run(
            [str(amberhome / "bin/cpptraj"), "-i", str(script)],
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    if not trajectory.is_file():
        raise ValueError(f"Expected cpptraj coordinates missing: {trajectory}")
    structure = AmberParm(str(system / "md_solvated.prmtop"), xyz=str(trajectory))
    if (
        structure.get_coordinates().shape[0] != 1
        or not np.isfinite(structure.coordinates).all()
    ):
        raise ValueError("Expected exactly one finite coordinate frame")
    structure.save(str(target / "solvated.rst7"), overwrite=False)
    structure[":1-1313"].save(str(target / "protein_substrate.pdb"), overwrite=False)
    return {
        **frame,
        "replica": replica,
        "directory": str(target),
        "solvated_coordinates": _source_record(trajectory),
        "atom_count": len(structure.atoms),
        "water_count": sum(r.name == "WAT" for r in structure.residues),
    }


def main() -> None:
    """Audit explicitly configured sources in a fresh output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.manifest.read_text())
    args.output.mkdir(exist_ok=False)
    output = args.output.resolve()
    system, amberhome = Path(config["system"]), Path(config["amberhome"])
    records = [_production_record(system / replica) for replica in config["replicas"]]
    seeds = [row["declared"]["ig"] for row in records]
    if len(set(seeds)) != len(seeds) or any(seed <= 0 for seed in seeds):
        raise ValueError("Replicas must declare distinct positive thermostat seeds")
    audit = audit_es_topology(system / "md_dry.prmtop", config["atom_map"], 1313)
    native = _native_mapping(Path(config["native_reference"]), config["native_mapping"])
    topology = AmberParm(str(system / "md_solvated.prmtop"))
    for row in config["native_mapping"]:
        residue = topology.residues[row["topology_residue"] - 1]
        if residue.name != row["resname"]:
            raise ValueError(f"Mapped topology residue mismatch: {row}")
    source_paths = [
        args.manifest,
        Path(__file__),
        Path(config["geometry_csv"]),
        Path(config["native_reference"]),
        ROOT / "src/nyl12_covalent/es_topology_audit.py",
        ROOT / "src/nyl12_covalent/es_snapshot_selection.py",
    ]
    source_paths += [system / name for name in config["source_files"]]
    source_paths += [
        system / replica / name
        for replica in config["replicas"]
        for name in config["replica_files"]
    ]
    LOGGER.info("Hashing original sources; large trajectories remain in place")
    sources = [_source_record(path) for path in source_paths]
    frozen = output / "source"
    frozen.mkdir()
    for path in source_paths[:6]:
        shutil.copy2(path, frozen / path.name)
    preliminary = {
        "sources": sources,
        "production": records,
        "topology": audit,
        "native_mapping": native,
    }
    (output / "input_audit.json").write_text(json.dumps(preliminary, indent=2) + "\n")
    frames = []
    for record in records:
        replica = Path(record["directory"]).name
        selected = select_es_frames(
            Path(config["geometry_csv"]),
            config["model"],
            replica,
            [tuple(w) for w in config["windows_ns"]],
            config["geometry_stride"],
            record["frame_dt_ns"],
        )
        for frame in selected:
            LOGGER.info("Extract %s full frame %s", replica, frame["full_frame"])
            frames.append(_extract_frame(system, replica, frame, output, amberhome))
    result = {
        **preliminary,
        "selected_frames": frames,
        "acceptance": {
            "reuse_audit_complete": True,
            "qmmm_hamiltonian_validated": False,
            "production_free_energy_authorized_by_audit": False,
        },
        "limitations": [
            "Independent seeds do not prove ensemble independence.",
            "Selected frames are conditional seeds, not occupancy estimates.",
            "AM1-BCC provenance requires the original charge-generation record.",
        ],
    }
    (output / "reuse_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    LOGGER.info("Completed source audit and %d solvated seed exports", len(frames))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
