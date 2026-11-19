"""Audit WAT1576 contacts in water-in-QM ASH274 start and endpoint."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import parmed

RUN = Path(__file__).resolve().parent
TOPOLOGY = RUN.parents[2] / "00_system" / "ash274_solvated.prmtop"
INPUTS = {
    "start": RUN / "starting.ncrst",
    "endpoint": RUN / "optimized.ncrst",
}
OUTPUT = RUN / "water_1576_full_hbond_audit_final.json"

# Residue IDs follow Amber topology numbering. The leaving-group nitrogen is
# the PA6 amide N at residue 1315, as used by the endpoint preparation audit.
TARGETS = {
    "Thr929:OG1": (929, "OG1"),
    "Thr929:HG1": (929, "HG1"),
    "Asp272:OD1": (968, "OD1"),
    "Asp272:OD2": (968, "OD2"),
    "Asp274:OD1": (970, "OD1"),
    "Asp274:OD2": (970, "OD2"),
    "Thr929:N": (929, "N"),
    "leaving_N": (1315, "N"),
}
PROTEIN_DONOR_H = {
    "Thr929:HG1": ("Thr929:OG1", 929, "HG1"),
    "Asp274:HD2": ("Asp274:OD2", 970, "HD2"),
    "Thr929:H1": ("Thr929:N", 929, "H1"),
    "Thr929:H2": ("Thr929:N", 929, "H2"),
}
ACCEPTORS = (
    "Thr929:OG1",
    "Asp272:OD1",
    "Asp272:OD2",
    "Asp274:OD1",
    "Asp274:OD2",
    "leaving_N",
)


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_coordinate(structure: Any, residue: int, name: str) -> np.ndarray:
    """Return a unique atom coordinate by Amber residue number and name."""
    matches = [
        atom
        for atom in structure.residues[residue - 1].atoms
        if atom.name.strip() == name
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected unique atom {residue}:{name}; found {len(matches)}")
    return np.asarray(structure.coordinates[matches[0].idx], dtype=float)


def distance(first: np.ndarray, second: np.ndarray) -> float:
    """Return Cartesian distance in angstroms."""
    return float(np.linalg.norm(first - second))


def dha_angle(donor: np.ndarray, hydrogen: np.ndarray, acceptor: np.ndarray) -> float:
    """Return donor-hydrogen-acceptor angle in degrees."""
    left = donor - hydrogen
    right = acceptor - hydrogen
    cosine = np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def measure(structure: Any) -> dict[str, Any]:
    """Measure target distances and both protein-to-water and water-to-target angles."""
    target_coords = {
        label: atom_coordinate(structure, residue, atom)
        for label, (residue, atom) in TARGETS.items()
    }
    water = {name: atom_coordinate(structure, 1576, name) for name in ("O", "H1", "H2")}
    # Include the protonated Asp274 donor hydrogen in the distance inventory.
    asp_hd2 = atom_coordinate(structure, 970, "HD2")
    distances = {
        target: {wat: distance(coord, wat_coord) for wat, wat_coord in water.items()}
        for target, coord in target_coords.items()
    }
    distances["Asp274:HD2"] = {
        wat: distance(asp_hd2, wat_coord) for wat, wat_coord in water.items()
    }

    angles: dict[str, float] = {}
    # Each listed protein donor hydrogen points toward water oxygen.
    for hydrogen_label, (donor_label, residue, atom) in PROTEIN_DONOR_H.items():
        hydrogen = (
            asp_hd2
            if hydrogen_label == "Asp274:HD2"
            else atom_coordinate(structure, residue, atom)
        )
        angles[f"{donor_label}-{hydrogen_label}...WAT1576:O"] = dha_angle(
            target_coords[donor_label], hydrogen, water["O"]
        )
    # Each water proton points toward every specified acceptor atom.
    for hname in ("H1", "H2"):
        for acceptor_label in ACCEPTORS:
            angles[f"WAT1576:O-{hname}...{acceptor_label}"] = dha_angle(
                water["O"], water[hname], target_coords[acceptor_label]
            )
    return {"distances_A": distances, "donor_h_acceptor_angles_deg": angles}


def main() -> None:
    """Write a provenance-hashed start/endpoint WAT1576 geometry report."""
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    acceptance = json.loads((RUN / "acceptance.json").read_text())
    if acceptance["process"]["exit_status"] != 0 or not acceptance["force_convergence"]:
        raise ValueError("Endpoint is not normally completed and force-converged")
    topology = parmed.load_file(str(TOPOLOGY))
    structures = {
        label: parmed.load_file(str(TOPOLOGY), xyz=str(path))
        for label, path in INPUTS.items()
    }
    inventory = {
        str(residue): sorted(
            atom.name.strip() for atom in topology.residues[residue - 1].atoms
        )
        for residue in (929, 968, 970, 1315, 1576)
    }
    if not {"O", "H1", "H2"}.issubset(inventory["1576"]):
        raise ValueError(f"Unexpected WAT1576 atom inventory: {inventory['1576']}")
    report = {
        "system": "force-converged water-in-QM ASH274 E-S reoptimization",
        "water_residue": 1576,
        "substrate_leaving_nitrogen": "Amber residue 1315 atom N",
        "residue_atom_inventory": inventory,
        "measurements": {
            label: measure(structure) for label, structure in structures.items()
        },
        "interpretation": (
            "One optimized endpoint snapshot can establish local geometric "
            "compatibility only; it does not establish water occupancy, proton "
            "transfer, or a catalytic relay. Angles are Cartesian, with no PBC "
            "minimum-image adjustment."
        ),
        "source_sha256": {
            str(path.resolve()): sha256(path)
            for path in (
                Path(__file__),
                TOPOLOGY,
                *INPUTS.values(),
                RUN / "acceptance.json",
                RUN / "preparation.json",
            )
        },
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
