"""Measure water1576 geometry in the converged corrected-shell endpoint."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import parmed

RUN = Path(__file__).resolve().parent
CASE = RUN.parent
TOPOLOGY = CASE.parent.parent / "00_system/ash274_solvated.prmtop"


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def coordinate(structure: Any, residue: int, atom_name: str) -> np.ndarray:
    """Return one exact atom coordinate."""
    matches = [
        atom
        for atom in structure.residues[residue - 1].atoms
        if atom.name.strip() == atom_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one atom {residue}:{atom_name}, found {len(matches)}"
        )
    return np.asarray(structure.coordinates[matches[0].idx], dtype=float)


def angle_degrees(first: np.ndarray, vertex: np.ndarray, last: np.ndarray) -> float:
    """Return the first-vertex-last angle in degrees."""
    left = first - vertex
    right = last - vertex
    cosine = np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def main() -> None:
    """Write hashed initial and converged endpoint water measurements."""
    output = RUN / "water_1576_geometry.json"
    if output.exists():
        raise FileExistsError(output)
    acceptance = json.loads((RUN / "acceptance.json").read_text())
    if acceptance["process"]["exit_status"] != 0 or not acceptance["force_convergence"]:
        raise ValueError(
            "Water endpoint measurement requires normal force-converged run"
        )
    topology = parmed.load_file(str(TOPOLOGY))
    initial = parmed.load_file(str(TOPOLOGY), xyz=str(RUN / "starting.ncrst"))
    final = parmed.load_file(str(TOPOLOGY), xyz=str(RUN / "optimized.ncrst"))
    names = {atom.name.strip() for atom in topology.residues[1575].atoms}
    oxygen = "O" if "O" in names else "OW"
    if not {oxygen, "H1", "H2"} <= names:
        raise ValueError(f"Unexpected WAT1576 atom inventory: {sorted(names)}")
    metrics: dict[str, Any] = {}
    for label, structure in (("initial", initial), ("converged_endpoint", final)):
        o = coordinate(structure, 1576, oxygen)
        h1 = coordinate(structure, 1576, "H1")
        h2 = coordinate(structure, 1576, "H2")
        og1 = coordinate(structure, 929, "OG1")
        metrics[label] = {
            "O_to_Thr929_OG1_A": float(np.linalg.norm(o - og1)),
            "O_H1_A": float(np.linalg.norm(o - h1)),
            "O_H2_A": float(np.linalg.norm(o - h2)),
            "H1_H2_A": float(np.linalg.norm(h1 - h2)),
            "H1_to_Thr929_OG1_A": float(np.linalg.norm(h1 - og1)),
            "H2_to_Thr929_OG1_A": float(np.linalg.norm(h2 - og1)),
            "O_H1_OG1_angle_deg": angle_degrees(o, h1, og1),
            "O_H2_OG1_angle_deg": angle_degrees(o, h2, og1),
        }
    metrics["converged_endpoint"]["O_displacement_from_start_A"] = float(
        np.linalg.norm(
            coordinate(final, 1576, oxygen) - coordinate(initial, 1576, oxygen)
        )
    )
    report = {
        "residue": 1576,
        "atom_names": sorted(names),
        "active_in_corrected_shell": True,
        "interpretation": "Final geometry is from the force-converged endpoint; water was not promoted to QM.",
        "measurements": metrics,
        "source_sha256": {
            str(path.resolve()): sha256(path)
            for path in (
                Path(__file__),
                TOPOLOGY,
                RUN / "starting.ncrst",
                RUN / "optimized.ncrst",
                RUN / "shell_audit.json",
                RUN / "acceptance.json",
            )
        },
    }
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
