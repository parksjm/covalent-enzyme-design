"""Audit final two retained calls for ordinary 10 Å cutoff crossings."""

from __future__ import annotations

import hashlib
import json
import logging
import sys
from itertools import product
from pathlib import Path

import numpy as np
from parmed import unit
from parmed.amber import Rst7
from parmed.geometry import box_lengths_and_angles_to_vectors
from scipy.spatial import cKDTree

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = ROOT / (
    "qmmm/pa6_flipped_model30/07_ash274_comparison/"
    "08_activation_attack_grid/occupancy_callback_pilot"
)
sys.path[:0] = [str(ROOT / "src"), str(ROOT), str(BASE)]
from runtime_support import load_model  # noqa: E402


def main() -> None:
    """Write final pair distances close to the declared 10 Å cutoff."""
    model = load_model()
    paths = [BASE / f"pilot/point_A/eval_{i:05d}.input.ncrst" for i in (283, 284)]
    restarts = [Rst7.open(str(path)) for path in paths]
    arrays = [np.asarray(rst.coordinates).reshape(-1, 3) for rst in restarts]
    cell = np.asarray(
        box_lengths_and_angles_to_vectors(*restarts[1].box).value_in_unit(
            unit.angstroms
        )
    )
    inverse = np.linalg.inv(cell)
    wrapped = ((arrays[1] @ inverse) % 1) @ cell
    shifts = np.asarray(list(product((-1, 0, 1), repeat=3))) @ cell
    tree = cKDTree(wrapped)
    active = np.flatnonzero(model.active)
    pairs = set()
    for start in range(0, len(active), 128):
        chunk = active[start : start + 128]
        images = (wrapped[chunk][None, :, :] + shifts[:, None, :]).reshape(-1, 3)
        sparse = cKDTree(images).sparse_distance_matrix(
            tree, 10.0000002, output_type="ndarray"
        )
        near = sparse[np.abs(sparse["v"] - 10.0) < 0.0000002]
        for item in near:
            a, b = int(chunk[item["i"] % len(chunk)]), int(item["j"])
            pairs.add(tuple(sorted((a, b))))
    rows = []
    for a, b in sorted(pairs):
        distances = []
        for xyz in arrays:
            delta = xyz[a] - xyz[b]
            delta -= np.rint(delta @ inverse) @ cell
            distances.append(
                float(np.linalg.norm(delta[None, :] + shifts, axis=1).min())
            )
        atoms = [model.topology.atoms[i] for i in (a, b)]
        rows.append(
            {
                "atoms_1based": [a + 1, b + 1],
                "names": [atom.name for atom in atoms],
                "residues": [
                    f"{atom.residue.name}{atom.residue.idx + 1}" for atom in atoms
                ],
                "charges_e": [atom.charge for atom in atoms],
                "qm": [bool(i in model.qm) for i in (a, b)],
                "distances_A": distances,
                "crosses_10_A": bool((distances[0] <= 10) != (distances[1] <= 10)),
            }
        )
    result = {
        "scope": "Final calls283/284 only; retained coordinates; no engine",
        "cutoff_shell_A": 0.0000002,
        "near_cutoff_pairs": rows,
        "sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths + [Path(__file__)]
        },
    }
    out = BASE / "line_search_diagnosis/cutoff_pairs.json"
    with out.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    logging.info(
        "Found %d near-cutoff pairs; %d crossings",
        len(rows),
        sum(row["crosses_10_A"] for row in rows),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
