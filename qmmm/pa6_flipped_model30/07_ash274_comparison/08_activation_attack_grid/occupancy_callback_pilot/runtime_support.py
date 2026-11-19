"""Explicit scientific inputs and native-evaluation gates for occupancy pilot."""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import parmed
import scipy
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from nyl12_covalent.native_sander_cli import (
    NativeSanderEvaluator,
    load_calibration,
    source_sha256,
)
from nyl12_covalent.occupancy_bias import OccupancyBiasConfig, evaluate_bias
from scripts.prepare_dimer_continuation import read_xyz_frames

ROOT = Path("/data/jpc/covalent-enzyme-design")
GRID = ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid"
SOURCE = (
    GRID.parent / "06_qmmm_optimization/replica_1_tier_A_seed/water_qm_reoptimization"
)
AMBER = Path("/data/jpc/code/ambertools26/venv")
CALIBRATION = ROOT / (
    "qmmm/pa6_flipped_model30/05_neutral_reaction_path/"
    "07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113/"
    "native_gradient_0K_onecycle_20261113_attempt2/scientific_acceptance.json"
)
POINTS = [
    ("A", 2.10, 3.125178400033067),
    ("B", 2.10, 2.60),
    ("C", 2.10, 2.10),
    ("D", 2.10, 1.60),
]


@dataclass
class Model:
    """Validated explicit model, source provenance and native settings."""

    topology: Any
    topology_path: Path
    restart: Rst7
    coordinates: np.ndarray
    active: np.ndarray
    qm: np.ndarray
    hydrogens: np.ndarray
    mdin: str
    source_paths: list[Path]
    scale: float


def configuration(
    model: Model, q_target: float, r_target: float
) -> OccupancyBiasConfig:
    """Return the declared transfer and independent attack bias definition."""
    return OccupancyBiasConfig(
        13915,
        13904,
        model.hydrogens,
        (13915, 19919),
        q_target,
        r_target,
        100.0,
        200.0,
        1.25,
    )


def save(path: Path, data: dict[str, Any]) -> None:
    """Write a new authoritative JSON record without overwriting evidence."""
    with path.open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def load_model() -> Model:
    """Validate mapping, accepted source, calibration hashes and input geometry."""
    topology_path = GRID.parent / "00_system/ash274_solvated.prmtop"
    source_mdin = SOURCE / "minimize.mdin"
    source_acceptance = SOURCE / "acceptance.json"
    accepted = json.loads(source_acceptance.read_text())
    if accepted["status"] != "accepted" or not accepted["acceptance"]["accepted"]:
        raise ValueError("Original water-QM model is not accepted")
    points = json.loads((GRID / "grid_points.json").read_text())
    point = points["points"]["q_004_attack_001"]
    restart_path = GRID / point["restart"]
    if source_sha256(restart_path) != point["restart_sha256"]:
        raise ValueError("Initial grid restart hash mismatch")
    restart = Rst7.open(str(restart_path))
    xyz = np.asarray(restart.coordinates, dtype=float).reshape(-1, 3).copy()
    if restart.box is None or not np.isfinite(restart.box).all():
        raise ValueError("Source periodic box missing/invalid")
    topology = parmed.load_file(str(topology_path))
    if len(xyz) != len(topology.atoms):
        raise ValueError("Topology/restart atom count mismatch")
    text = source_mdin.read_text()
    qmmask = re.findall(r"\bqmmask='([^']+)'", text)
    masks = re.findall(r"\bactive='([^']+)'", text)
    if len(qmmask) != 1 or len(masks) != 1:
        raise ValueError("Expected one explicit QM and active mask")
    qm = np.array(list(AmberMask(topology, qmmask[0]).Selected()), dtype=int)
    active = np.array(AmberMask(topology, masks[0]).Selection(), dtype=bool)
    hydrogens = np.array([i for i in qm if topology.atoms[i].atomic_number == 1])
    protocol = json.loads((GRID / "protocol.json").read_text())
    if qm.tolist() != [a["amber_atom_id"] - 1 for a in protocol["qm_atoms_ordered"]]:
        raise ValueError("QM map differs from declared grid source")
    if (len(qm), len(hydrogens), int(active.sum())) != (102, 54, 4662):
        raise ValueError("Expected 102 real QM atoms, 54 H, 4662 active atoms")
    if not active[qm].all() or not re.search(r"qmcharge\s*=\s*-1\s*,", text):
        raise ValueError("QM active selection or charge is inconsistent")
    # Preserve the complete source Hamiltonian prefix; change only cycle count
    # and explicitly set the previously implicit zero electronic temperature.
    prefix, separator, _ = text.partition("&dlfind")
    if not separator or "dftb_telec" in prefix.lower():
        raise ValueError("Unexpected source DL-FIND/temperature configuration")
    prefix, count = re.subn(r"\bmaxcyc=2000\b", "maxcyc=1", prefix)
    if count != 1:
        raise ValueError("Unexpected source cycle count")
    prefix = prefix.replace("qm_theory='DFTB3',", "qm_theory='DFTB3', dftb_telec=0.0,")
    if prefix.count("dftb_telec=0.0") != 1:
        raise ValueError("Source method must explicitly be DFTB3")
    mdin = prefix + (
        "&dlfind\n optalg='LBFGS', crdrep='CART', dimer=.false.,\n"
        f" active='{masks[0]}', maxcycle=1, tol=1.d-4, tole=1.d-5,\n"
        " maxstep=0.1, verbosity=2,\n/\n"
    )
    scale, _ = load_calibration(CALIBRATION)
    paths = [
        topology_path,
        restart_path,
        source_mdin,
        source_acceptance,
        SOURCE / "optimized.ncrst",
        GRID / "grid_points.json",
        GRID / "protocol.json",
        CALIBRATION,
        AMBER / "bin/sander",
        AMBER / "lib/libdlfind.so",
        AMBER / "lib/libsander.so",
    ]
    model = Model(
        topology, topology_path, restart, xyz, active, qm, hydrogens, mdin, paths, scale
    )
    evaluate_bias(xyz, configuration(model, 2.1, POINTS[0][2]))
    return model


def environment() -> dict[str, str]:
    """Set and return the explicit single-threaded native runtime environment."""
    settings = {
        "AMBERHOME": str(AMBER),
        "LD_LIBRARY_PATH": str(AMBER / "lib"),
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    os.environ.update(settings)
    return settings


def snapshot_sources(destination: Path, paths: list[Path]) -> dict[str, str]:
    """Copy explicitly indexed sources and return authoritative source hashes."""
    destination.mkdir()
    identities = {}
    for index, path in enumerate(paths):
        digest = source_sha256(path)
        copied = destination / f"{index:02d}_{path.name}"
        shutil.copy2(path, copied)
        if source_sha256(copied) != digest:
            raise ValueError(f"Snapshot hash mismatch: {path}")
        identities[str(path)] = digest
    return identities


class CheckedPhysical:
    """Native callable with first-call geometry/map and per-call SCC gates."""

    def __init__(self, native: NativeSanderEvaluator, active: np.ndarray) -> None:
        """Retain native evaluator and exact boolean active selection."""
        self.native = native
        self.active = active
        self.counter: int = 0
        self.initial_geometry_error_A: float | None = None

    def __call__(self, coordinates: np.ndarray) -> tuple[float, np.ndarray]:
        """Evaluate native geometry and fail explicitly on SCC/map failure."""
        self.counter += 1
        energy, gradient = self.native(coordinates)
        tag = f"eval_{self.native.counter:05d}"
        mdout = (self.native.branch_dir / f"{tag}.mdout").read_text()
        if re.search(r"(?:SCC|SCF)[^\n]*(?:failed|not converg|failure)", mdout, re.I):
            raise ValueError("Native SCC/SCF failure detected")
        energies = re.findall(r"Energy calculation finished, energy:", mdout)
        forces = read_xyz_frames(
            self.native.branch_dir / f"{tag}.path_force.xyz", int(self.active.sum())
        )
        if len(energies) != 1 or len(forces) != 1:
            raise ValueError("Native call must contain exactly one energy/gradient")
        if self.counter == 1:
            frames = read_xyz_frames(
                self.native.branch_dir / f"{tag}.path_active.xyz",
                int(self.active.sum()),
            )
            if len(frames) != 1:
                raise ValueError("Native call must contain exactly one active geometry")
            self.initial_geometry_error_A = float(
                np.max(np.abs(frames[0] - coordinates[self.active]))
            )
            if self.initial_geometry_error_A > 2e-6:
                raise ValueError("Native path_active mapping/initial geometry mismatch")
        with (self.native.branch_dir / "native_evidence.jsonl").open("a") as stream:
            stream.write(json.dumps(self.native.evidence[-1]) + "\n")
        return energy, gradient


def state_record(model: Model, xyz: np.ndarray) -> dict[str, Any]:
    """Record both occupancies, original QM bonds and all real-H ownership."""
    bias = evaluate_bias(xyz, configuration(model, 2.1, POINTS[0][2]))
    heavy = np.array(
        [i for i in model.qm if model.topology.atoms[i].atomic_number != 1]
    )
    owners = []
    for h in model.hydrogens:
        distances = np.linalg.norm(xyz[heavy] - xyz[h], axis=1)
        closest = np.argsort(distances)[:2]
        d1, d2 = distances[closest]
        owners.append(
            {
                "hydrogen_1based": int(h + 1),
                "nearest_qm_heavy_1based": int(heavy[closest[0]] + 1),
                "nearest_distance_A": float(d1),
                "second_qm_heavy_1based": int(heavy[closest[1]] + 1),
                "second_distance_A": float(d2),
                "classification": "shared_ambiguous"
                if d2 <= 1.35
                else ("covalent" if d1 <= 1.25 else "unassigned"),
            }
        )
    qm_set = set(model.qm.tolist())
    bonds = [
        {
            "atoms_1based": [b.atom1.idx + 1, b.atom2.idx + 1],
            "distance_A": float(np.linalg.norm(xyz[b.atom1.idx] - xyz[b.atom2.idx])),
        }
        for b in model.topology.bonds
        if b.atom1.idx in qm_set and b.atom2.idx in qm_set
    ]
    return {
        "CN_N": bias.cn_acceptor,
        "CN_OG1": bias.cn_donor,
        "q": bias.q,
        "attack_A": bias.r,
        "original_qm_bond_distances": bonds,
        "proton_owners": owners,
    }


def preflight(model: Model) -> dict[str, Any]:
    """Validate declared units and model without invoking the native evaluator."""
    config = configuration(model, 2.1, POINTS[0][2])
    initial = evaluate_bias(model.coordinates, config)
    config = configuration(model, initial.q - 0.2, initial.r)
    known = evaluate_bias(model.coordinates, config)
    if not np.isclose(known.energy, 2.0, atol=1e-10):
        raise ValueError("kq=100 and dq=.2 must produce 2 kcal/mol")
    from nyl12_covalent.native_sander_cli import source_sha256

    hashes = {str(path): source_sha256(path) for path in model.source_paths}
    return {
        "status": "passed_no_engine_execution",
        "real_qm_atoms": len(model.qm),
        "real_qm_hydrogens": len(model.hydrogens),
        "qm_charge": -1,
        "active_atoms": int(model.active.sum()),
        "source_sha256": hashes,
        "disk_free_bytes": shutil.disk_usage(ROOT).free,
        "scipy_version": scipy.__version__,
        "python": sys.executable,
        "bias_unit_check_kcal_mol": known.energy,
        "initial_state": state_record(model, model.coordinates),
    }
