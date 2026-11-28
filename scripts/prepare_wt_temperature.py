#!/usr/bin/env python3
"""Freeze accepted WT E-S endpoints and bounded 348.15 K equilibration inputs."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Any

import numpy as np
from parmed.amber import AmberParm, Rst7

from nyl12_covalent.hydration_followon_gates import digest, read
from nyl12_covalent.wt_temperature_protocol import STAGES, mdin_files

ROOT = Path(__file__).resolve().parents[1]
FOLLOWUP = ROOT / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup"
BASE = FOLLOWUP / "temperature_348K"
TOPOLOGY_SHA = "2a82bfac3c467f12db72e844b7a6c434ae4dfff20d351cce62670d6c1b37e3d5"
HOSTS = {
    "violet": {
        "engine": "/data/jpc/code/pmemd26/bin/pmemd.cuda_SPFP",
        "sha256": "7ed3f439d09d6beff94210bb3f356585e6918cd3f500dcc58046067ac8fd9602",
        "LD_LIBRARY_PATH": "/data/jpc/code/pmemd26/lib:/data/jpc/code/ambertools26/venv/lib",
    },
    "azalea": {
        "engine": "/data/jpc/code/pmemd26_blackwell_build_20261005/pmemd26_src/build_sysgcc_cuda128_v2/src/pmemd/src/pmemd.cuda_SPFP",
        "sha256": "f4b0b32eda07d1e6b4fbe597e8c3bd4d39cc492a2b01c19c13db7ebe92ff52b1",
        "LD_LIBRARY_PATH": "/data/jpc/code/cuda128_findcuda_compat/lib64:/data/jpc/code/pmemd26/lib:/data/jpc/code/ambertools26/venv/lib",
    },
}
LOGGER = logging.getLogger(__name__)


def write(path: Path, value: dict[str, Any]) -> None:
    """Save exclusive finite JSON preparation evidence."""
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def collect(host: str, source: Path, target: Path) -> None:
    """Copy a local or remote immutable artifact without changing its source."""
    if target.exists():
        raise FileExistsError(target)
    if host == "violet":
        shutil.copyfile(source, target)
    else:
        subprocess.run(["scp", f"{host}:{source}", str(target)], check=True)


def endpoint_identity(topology: Path, restart: Path) -> dict[str, Any]:
    """Verify a complete finite WT E-S restart with inherited velocities and box."""
    if digest(topology) != TOPOLOGY_SHA:
        raise ValueError("WT topology identity changed")
    parm = AmberParm(str(topology))
    frame = Rst7.open(str(restart))
    xyz = np.asarray(frame.coordinates).reshape(-1, 3)
    velocities = np.asarray(frame.velocities).reshape(-1, 3)
    box = np.asarray(frame.box)
    if (
        xyz.shape != (len(parm.atoms), 3)
        or velocities.shape != xyz.shape
        or box.shape != (6,)
    ):
        raise ValueError("Full-system coordinates, velocities or box missing")
    if (
        not all(np.isfinite(array).all() for array in (xyz, velocities, box))
        or not np.any(velocities)
        or not np.all(box[:3] > 0)
    ):
        raise ValueError("Invalid full-system restart")
    # Native 2fs time accumulation differs by0.000108ps after100ns.
    # Bound the discrepancy to half a timestep; hashes establish exact identity.
    if abs(frame.time - 100500.0) > 0.001:
        raise ValueError("Not the final accepted100ns production endpoint")
    ligand_names = {"P6N", "P6I", "P6C"}
    solute_ligands = [res for res in parm.residues if res.name in ligand_names]
    if not solute_ligands:
        raise ValueError("Expected parameterized PA6 E-S substrates missing")
    cross_bonds = [
        (b.atom1.idx + 1, b.atom2.idx + 1)
        for b in parm.bonds
        if (b.atom1.residue.name in ligand_names)
        != (b.atom2.residue.name in ligand_names)
    ]
    if cross_bonds:
        raise ValueError("Unexpected enzyme-substrate covalent bond")
    return {
        "atoms": len(parm.atoms),
        "restart_time_ps": frame.time,
        "coordinates_velocities_box_finite": True,
        "protein_substrate_cross_bonds": cross_bonds,
        "system": "Full-system classical WT E-S; unchanged topology and explicit solvent",
        "terminal_protonation": "NH2",
        "microstate": "ASP272/ASP274",
    }


def prepare() -> None:
    """Stage two independent thermal starts and reviewed host configurations only."""
    if (BASE / "plan.json").exists() or (BASE / "inputs").exists():
        raise FileExistsError("Thermal preparation already attempted; do not overwrite")
    BASE.mkdir(exist_ok=True)
    (BASE / "inputs").mkdir()
    (BASE / "provenance").mkdir()
    source_top = FOLLOWUP / "md/WT/replica_01/dynamics/system.prmtop"
    topology = BASE / "inputs/system.prmtop"
    collect("violet", source_top, topology)
    replicas = []
    for index, host in ((1, "violet"), (2, "azalea")):
        name = f"replica_{index:02d}"
        source = FOLLOWUP / "md/WT" / name
        provenance = BASE / "provenance" / name
        provenance.mkdir()
        start = BASE / "inputs" / f"{name}_accepted300K_100ns.rst7"
        files = {
            "dynamics/06_production_90ns.acceptance.json": "acceptance.json",
            "dynamics/extension_100ns.json": "extension_100ns.json",
            "job/status.json": "job_status.json",
            "dynamics/06_production_90ns.in": "native_300K.in",
            "dynamics/06_production_90ns.out": "native_300K.out",
            "100ns_process_complete.json": "100ns_process_complete.json",
        }
        collect(host, source / "dynamics/06_production_90ns.rst7", start)
        for relative, filename in files.items():
            collect(host, source / relative, provenance / filename)
        acceptance = read(provenance / "acceptance.json")
        job = read(provenance / "job_status.json")
        extension = read(provenance / "extension_100ns.json")
        if (
            acceptance["completed_steps"] != 45_000_000
            or acceptance["total_production_ns"] != 100
            or acceptance["temperature_gate_pass"] is not True
            or acceptance["restart_sha256"] != digest(start)
            or job["exit_status"] != 0
            or job["ended_utc"] is None
            or extension["returncode"] != 0
            or extension["ended_utc"] is None
            or read(provenance / "100ns_process_complete.json")["completed"] is not True
        ):
            raise ValueError(f"AcceptedWT100ns endpoint gates failed: {name}")
        source_top_hash = extension["input_sha256"][
            str(source / "dynamics/system.prmtop")
        ]
        if (
            source_top_hash != TOPOLOGY_SHA
            or job["input_sha256"][HOSTS[host]["engine"]] != HOSTS[host]["sha256"]
        ):
            raise ValueError("Source topology or host engine provenance mismatch")
        identity = endpoint_identity(topology, start)
        seed = 34800000 + index * 1000
        target = BASE / name
        target.mkdir()
        (target / "inputs").mkdir()
        for stage, text in mdin_files(seed).items():
            (target / "inputs" / f"{stage}.in").write_text(text)
        config = {
            "status": "PREPARED_ONLY_ROOT_REVIEW_BEFORE_NATIVE_LAUNCH",
            "replica": name,
            "source_host": host,
            "source_directory": str(source),
            "source_restart": str(source / "dynamics/06_production_90ns.rst7"),
            "source_acceptance": acceptance,
            "endpoint_identity": identity,
            "start_restart": str(start),
            "topology": str(topology),
            "start_sha256": digest(start),
            "topology_sha256": digest(topology),
            "seed_by_stage": {stage: seed + i for i, (stage, _) in enumerate(STAGES)},
            "stages": [
                {
                    "name": stage,
                    "steps": steps,
                    "dt_ps": 0.002,
                    "duration_ps": steps * 0.002,
                    "input": str(target / "inputs" / f"{stage}.in"),
                    "ntr": 0,
                }
                for stage, steps in STAGES
            ],
            "preferred_host": host,
            "host_runtime": HOSTS[host],
            "environment": {
                **job["env_overrides"],
                "CUDA_VISIBLE_DEVICES": "UNASSIGNED_REQUIRES_PREFLIGHT",
            },
            "output": str(target / "dynamics"),
            "native_launch_authorized": False,
        }
        write(target / "config.json", config)
        replicas.append(config)
    write(
        BASE / "plan.json",
        {
            "schema_version": 1,
            "status": "PREPARED_ONLY",
            "scope": "Two independent WT E-S starts; fullsystem200ps300->348.15KNVT+1nsNPT+2nsunrestrainedNPT; stop forreview, no100nsproduction",
            "replicas": replicas,
            "temperature_K": 348.15,
            "duration_per_replica_ps": 3200,
            "same_force_field_protonation_topology": True,
            "native_launch_authorized": False,
            "temperature_effect_claims": "Equilibration alone cannot establish a temperature effect; matched300Kcontrols and longer sampling need separate design",
        },
    )
    snapshot = BASE / "source_snapshot"
    snapshot.mkdir()
    sources = [
        ROOT / "scripts/prepare_wt_temperature.py",
        ROOT / "src/nyl12_covalent/wt_temperature_protocol.py",
        ROOT / "src/nyl12_covalent/hydration_followon_gates.py",
        ROOT / "docs/AMBER_WORKSTATION_RUNTIME.md",
    ]
    for source in sources:
        shutil.copyfile(source, snapshot / source.name)
    write(
        BASE / "source_manifest.json",
        {
            "sha256": {str(p): digest(p) for p in BASE.rglob("*") if p.is_file()},
            "original_sources": {str(p): digest(p) for p in sources},
            "preservation": "Original local and remote endpoints/provenance are unchanged; preparation copies retained",
        },
    )
    LOGGER.info("Prepared two348.15K E-S branches; no native GPU command executed")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    prepare()
