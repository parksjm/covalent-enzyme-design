"""One-frame functional watershell control, independent of the main analysis."""

from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent


def _record(path: Path) -> dict[str, str | int]:
    stat = path.stat()
    return {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def main() -> None:
    """Compare water-centered periodic counts to independently computed counts."""
    started = time.time()
    prior_path = (
        PROJECT
        / "output_nyl12_pa6_es_audit_20261024/reuse_audit_attempt1/reuse_audit.json"
    )
    prior = json.loads(prior_path.read_text())
    config = json.loads((PROJECT / "configs/pa6_wt_es_reuse_20261024.json").read_text())
    seed = prior["selected_frames"][0]
    coordinates = Path(seed["solvated_coordinates"]["path"])
    coordinate_record = _record(coordinates)
    if coordinate_record["sha256"] != seed["solvated_coordinates"]["sha256"]:
        raise ValueError("Accepted seed coordinates changed")
    topology = Path(config["system"]) / "md_solvated.prmtop"
    topology_record = next(
        row for row in prior["sources"] if row["path"] == str(topology)
    )
    stat = topology.stat()
    if (stat.st_size, stat.st_mtime_ns) != (
        topology_record["bytes"],
        topology_record["mtime_ns"],
    ):
        raise ValueError("Original topology metadata changed")
    structure = AmberParm(str(topology), xyz=str(coordinates))
    if structure.get_coordinates().shape[0] != 1:
        raise ValueError("Expected one accepted coordinate frame")
    water = next(residue for residue in structure.residues if residue.name == "WAT")
    oxygen = [
        atom for atom in water.atoms if atom.name == "O" and atom.atomic_number == 8
    ]
    if len(oxygen) != 1:
        raise ValueError("Ambiguous positive-control water oxygen")
    oxygen = oxygen[0]
    mask = f":{water.idx + 1}@O"
    box = np.asarray(structure.box, dtype=float)
    if box.shape != (6,) or not np.isfinite(box).all() or np.any(box[:3] <= 0):
        raise ValueError("Independent NumPy control requires finite positive box")
    atoms = [
        atom.idx
        for atom in structure.atoms
        if atom.residue.name == "WAT" and atom.name == "O"
    ]
    xyz = np.asarray(structure.coordinates)
    difference = xyz[atoms] - xyz[oxygen.idx]
    alpha, beta, gamma = np.deg2rad(box[3:])
    a, b, c = box[:3]
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz2 = c * c - cx * cx - cy * cy
    if cz2 <= 0:
        raise ValueError("Invalid triclinic lattice")
    lattice = np.array([[a, 0, 0], [b * np.cos(gamma), b * np.sin(gamma), 0],
                        [cx, cy, np.sqrt(cz2)]])
    fractional = difference @ np.linalg.inv(lattice)
    reduced = (fractional - np.round(fractional)) @ lattice
    distances27 = np.full(len(atoms), np.inf)
    distances125 = np.full(len(atoms), np.inf)
    for shift in itertools.product(range(-2, 3), repeat=3):
        candidate = np.linalg.norm(reduced + np.asarray(shift) @ lattice, axis=1)
        distances125 = np.minimum(distances125, candidate)
        if max(abs(value) for value in shift) <= 1:
            distances27 = np.minimum(distances27, candidate)
    if not np.allclose(distances27, distances125, atol=1e-10, rtol=0):
        raise ValueError("27-image and 125-image nearest-image checks disagree")
    distances = distances125
    expected = {
        str(cut): int(np.count_nonzero(distances < cut)) for cut in (3.5, 4.0, 5.0)
    }
    if expected["3.5"] < 1:
        raise ValueError("Self-distance failed independent positive control")
    attempt = ROOT / "positive_control_attempt2"
    attempt.mkdir(exist_ok=False)
    script = attempt / "positive_control.cpptraj"
    script.write_text(
        f"parm {topology}\ntrajin {coordinates} 1 1\n"
        f"watershell {mask} out {attempt / 'lower35.dat'} lower 3.5 upper 5.0 :WAT@O\n"
        f"watershell {mask} out {attempt / 'lower4.dat'} lower 4.0 upper 5.0 :WAT@O\nrun\n"
    )
    amberhome = Path(config["amberhome"])
    executable = amberhome / "bin/cpptraj"
    env = os.environ.copy()
    env.update(
        AMBERHOME=str(amberhome),
        OMP_NUM_THREADS="1",
        LD_LIBRARY_PATH=str(amberhome / "lib"),
    )
    with (attempt / "cpptraj.log").open("w") as log:
        process = subprocess.run(
            [str(executable), "-i", str(script)],
            cwd=attempt,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    comparisons = {}
    for name, cutoff in (("lower35", 3.5), ("lower4", 4.0)):
        path = attempt / f"{name}.dat"
        lines = [
            line
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
        if len(lines) != 1 or len(lines[0].split()) != 3:
            raise ValueError("Expected one three-column control record")
        frame, lower, upper = map(int, lines[0].split())
        comparisons[name] = {
            "frame": frame,
            "lower": lower,
            "upper": upper,
            "independent_numpy_lower": expected[str(cutoff)],
            "independent_numpy_upper": expected["5.0"],
            "passed": frame == 1
            and lower == expected[str(cutoff)]
            and upper == expected["5.0"]
            and lower > 0,
        }
    accepted = process.returncode == 0 and all(
        row["passed"] for row in comparisons.values()
    )
    result = {
        "accepted": accepted,
        "process_exit_status": process.returncode,
        "elapsed_seconds": time.time() - started,
        "solute_mask": mask,
        "solvent_mask": ":WAT@O",
        "periodic_imaging": True, "box_lengths_angles": box.tolist(), "independent_image_search": "Triclinic lattice; 27 and125 nearest-image searches agree",
        "known_water_residue": water.idx + 1,
        "known_oxygen_index": oxygen.idx + 1,
        "water_oxygen_count": len(atoms),
        "comparison": comparisons,
        "sources": [
            _record(Path(__file__)),
            _record(executable),
            coordinate_record,
            {
                **topology_record,
                "integrity_policy": "Inherited SHA256; current size and mtime checked",
            },
        ],
        "scope": "One-frame positive functional control of water counting; self-distance included; no new dynamics or mechanism claim.",
    }
    (ROOT / "positive_control.json").write_text(json.dumps(result, indent=2) + "\n")
    if not accepted:
        raise ValueError("Positive functional control failed")


if __name__ == "__main__":
    main()
