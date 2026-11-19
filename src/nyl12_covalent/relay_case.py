"""Validate explicit relay cases and audit native evaluations and chemistry."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
import shutil
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from .native_sander_cli import NativeSanderEvaluator, load_calibration, source_sha256
from .proton_coordination import smooth_coordination
from scripts.prepare_dimer_continuation import read_xyz_frames


@dataclass
class RelayCase:
    """Validated configuration with explicit zero-based computational selections."""

    config: dict[str, Any]
    topology: Any
    restart: Rst7
    coordinates: np.ndarray
    active: np.ndarray
    qm: np.ndarray
    hydrogens: np.ndarray
    atom_map: dict[str, int]
    scale: float


def write_json(path: Path, record: dict[str, Any]) -> None:
    """Write a new immutable authoritative JSON record."""
    with path.open("x") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")


def load_case(case_path: Path) -> RelayCase:
    """Load and validate a case, Hamiltonian masks, PBC and electron parity.

    Args:
        case_path: Root-prepared JSON specifying all scientific inputs.

    Returns:
        Model with immutable-source coordinates and validated native force scale.

    Raises:
        ValueError: If atom counts, masks, charges, method or calibration disagree.
    """
    cfg = json.loads(case_path.read_text())
    topology = parmed.load_file(cfg["topology"])
    restart = Rst7.open(cfg["source_restart"])
    xyz = np.asarray(restart.coordinates, dtype=float).reshape(-1, 3).copy()
    if xyz.shape != (len(topology.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Topology/source coordinate count or finiteness mismatch")
    box = np.asarray(restart.box, dtype=float)
    if box.shape != (6,) or not np.isfinite(box).all() or np.any(box <= 0):
        raise ValueError("Expected finite positive six-component periodic box")
    selections = []
    for key in ("qm_atom_ids", "active_atom_ids"):
        values = np.asarray(cfg[key])
        if (
            values.ndim != 1
            or not len(values)
            or values.dtype.kind not in "iu"
            or len(set(values.tolist())) != len(values)
            or values.min() < 1
            or values.max() > len(xyz)
        ):
            raise ValueError(f"Invalid distinct one-based selection: {key}")
        selections.append(values.astype(int) - 1)
    qm, active_ids = selections
    active = np.zeros(len(xyz), dtype=bool)
    active[active_ids] = True
    hydrogens = np.array([i for i in qm if topology.atoms[i].atomic_number == 1])
    if len(qm) != 102 or len(hydrogens) != 54 or not active[qm].all():
        raise ValueError("Expected 102 active QM atoms including 54 real hydrogens")
    text = cfg["mdin"]
    for key, expected in (("qmmask", qm), ("active", np.sort(active_ids))):
        masks = re.findall(rf"\b{key}\s*=\s*'([^']+)'", text, re.I)
        if len(masks) != 1 or not np.array_equal(
            np.array(list(AmberMask(topology, masks[0]).Selected())), expected
        ):
            raise ValueError(f"Declared {key} does not match native mask")
    for key, value in (
        ("maxcyc", 1),
        ("maxcycle", 1),
        ("dftb_telec", 0),
        ("qmcharge", cfg["qm_charge"]),
    ):
        matches = re.findall(rf"\b{key}\s*=\s*([-+\d.EeDd]+)", text, re.I)
        if (
            len(matches) != 1
            or float(matches[0].replace("d", "e").replace("D", "E")) != value
        ):
            raise ValueError(f"Expected explicit native setting {key}={value}")
    if type(cfg["qm_charge"]) is not int or cfg["qm_charge"] != -1:
        raise ValueError("Primary relay case requires unchanged formal QM charge -1")
    if not re.search(r"qm_theory\s*=\s*'DFTB3'", text, re.I):
        raise ValueError("Expected explicit DFTB3 method")
    qm_set = set(qm.tolist())
    links = sum(
        (b.atom1.idx in qm_set) != (b.atom2.idx in qm_set) for b in topology.bonds
    )
    electrons = (
        sum(topology.atoms[i].atomic_number for i in qm) + links - cfg["qm_charge"]
    )
    if electrons % 2:
        raise ValueError("QM region plus link-H electron parity is odd")
    mapping = cfg["atom_map"]
    if any(type(i) is not int or i - 1 not in qm_set for i in mapping.values()):
        raise ValueError("Every explicit reactive atom must be in the QM region")
    expected_elements = {
        "thr_n": 7,
        "thr_h1": 1,
        "thr_h2": 1,
        "thr_cb": 6,
        "thr_og": 8,
        "proton": 1,
        "carbon": 6,
        "oxygen": 8,
        "leaving_n": 7,
        "leaving_h": 1,
        "water_o": 8,
        "water_h1": 1,
        "water_h2": 1,
        "ash_oxygen": 8,
        "ash_proton": 1,
        "asp_o1": 8,
        "asp_o2": 8,
        "phe_n": 7,
        "asn_nd2": 7,
    }
    if any(
        topology.atoms[mapping[key] - 1].atomic_number != z
        for key, z in expected_elements.items()
    ):
        raise ValueError("Reactive atom mapping element mismatch")
    for key in ("max_atom_step_A", "force_gate"):
        if not np.isfinite(cfg[key]) or cfg[key] <= 0:
            raise ValueError(f"Expected finite positive {key}")
    if type(cfg["memory"]) is not int or cfg["memory"] <= 0:
        raise ValueError("LBFGS memory must be a positive integer")
    scale, _ = load_calibration(Path(cfg["calibration"]))
    return RelayCase(cfg, topology, restart, xyz, active, qm, hydrogens, mapping, scale)


class CheckedPhysical:
    """Apply SCC and exactly-one-frame gates to each native evaluation."""

    def __init__(self, native: NativeSanderEvaluator, active: np.ndarray) -> None:
        """Retain evaluator and its exact active mask."""
        self.native, self.active = native, active

    def __call__(self, xyz: np.ndarray) -> tuple[float, np.ndarray]:
        """Evaluate and verify native force/geometry mapping at every call."""
        energy, gradient = self.native(xyz)
        tag = f"eval_{self.native.counter:05d}"
        directory = self.native.branch_dir
        mdout = (directory / f"{tag}.mdout").read_text()
        if re.search(r"(?:SCC|SCF)[^\n]*(?:failed|not converg|failure)", mdout, re.I):
            raise ValueError("Native SCC/SCF failure detected")
        frames = read_xyz_frames(
            directory / f"{tag}.path_force.xyz", int(self.active.sum())
        )
        if (
            len(re.findall("Energy calculation finished, energy:", mdout)) != 1
            or len(frames) != 1
        ):
            raise ValueError("Expected exactly one native energy and force frame")
        geometries = read_xyz_frames(
            directory / f"{tag}.path_active.xyz", int(self.active.sum())
        )
        if (
            len(geometries) != 1
            or np.max(np.abs(geometries[0] - xyz[self.active])) > 2e-6
        ):
            raise ValueError("Native active geometry mapping mismatch")
        with (directory / "native_evidence.jsonl").open("a") as stream:
            stream.write(json.dumps(self.native.evidence[-1], allow_nan=False) + "\n")
        return energy, gradient


def state_record(model: RelayCase, xyz: np.ndarray) -> dict[str, Any]:
    """Record relay coordination, end geometry and every real QM H owner."""
    ids = {key: value - 1 for key, value in model.atom_map.items()}
    values = {
        label: smooth_coordination(xyz, ids[key], model.hydrogens)[0]
        for label, key in (
            ("CN_N", "thr_n"),
            ("CN_OG1", "thr_og"),
            ("CN_water", "water_o"),
            ("CN_leave", "leaving_n"),
        )
    }
    values["q"] = values["CN_N"] - values["CN_OG1"]
    values["qCN_N-CN_OG1"] = values["q"]
    distances = {
        label: float(np.linalg.norm(xyz[ids[a]] - xyz[ids[b]]))
        for label, a, b in (
            ("rOG1C", "thr_og", "carbon"),
            ("rCN", "carbon", "leaving_n"),
            ("rCO", "carbon", "oxygen"),
            ("rNwater", "thr_n", "water_o"),
            ("rwaterLeave", "water_o", "leaving_n"),
        )
    }
    heavy = np.array(
        [i for i in model.qm if model.topology.atoms[i].atomic_number != 1]
    )
    owners = []
    for h in model.hydrogens:
        lengths = np.linalg.norm(xyz[heavy] - xyz[h], axis=1)
        nearest = np.argsort(lengths)[:2]
        first, second = (float(lengths[i]) for i in nearest)
        owners.append(
            {
                "hydrogen_1based": int(h + 1),
                "nearest_qm_heavy_1based": int(heavy[nearest[0]] + 1),
                "nearest_distance_A": first,
                "second_qm_heavy_1based": int(heavy[nearest[1]] + 1),
                "second_distance_A": second,
                "classification": "shared_ambiguous"
                if second <= 1.35
                else ("covalent" if first <= 1.25 else "unassigned"),
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
        **values,
        **distances,
        "distances_A": distances,
        "proton_owners": owners,
        "original_qm_bond_distances": bonds,
    }


def snapshot_stage_inputs(
    case_path: Path, start: Path, output: Path, model: RelayCase, lbfgs_source: Path
) -> None:
    paths = [
        case_path,
        start,
        Path(model.config["source_restart"]),
        Path(model.config["topology"]),
        Path(model.config["calibration"]),
        lbfgs_source,
    ]
    package = Path(__file__).parent
    paths += [
        package / name
        for name in (
            "relay_stage.py",
            "relay_case.py",
            "relay_stage_units.py",
            "reaction_bias.py",
            "native_sander_cli.py",
            "dlfind_lbfgs.py",
            "proton_coordination.py",
            "occupancy_bias.py",
        )
    ]
    paths.append(package.parents[1] / "scripts/prepare_dimer_continuation.py")
    amber = Path(model.config["amber_home"])
    binaries = [
        amber / name for name in ("bin/sander", "lib/libdlfind.so", "lib/libsander.so")
    ]
    snapshot = output / "source_snapshot"
    snapshot.mkdir()
    hashes = {}
    for index, path in enumerate(paths):
        digest = source_sha256(path)
        copied = snapshot / f"{index:02d}_{path.name}"
        shutil.copy2(path, copied)
        if source_sha256(copied) != digest:
            raise ValueError("Source snapshot hash mismatch")
        hashes[str(path)] = digest
    write_json(
        output / "source_manifest.json",
        {
            "source_sha256": hashes,
            "binaries": [
                {
                    "path": str(path),
                    "sha256": source_sha256(path),
                    "size_bytes": path.stat().st_size,
                    "retained_in_place": True,
                    "version_identity": "SHA256-pinned installed binary; not executed for version query",
                }
                for path in binaries
            ],
        },
    )
