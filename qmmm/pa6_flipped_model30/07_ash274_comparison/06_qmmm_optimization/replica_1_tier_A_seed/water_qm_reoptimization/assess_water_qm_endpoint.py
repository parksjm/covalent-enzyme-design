"""Assess proton retention and local water geometry after QM/MM relaxation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import parmed

RUN = Path(__file__).resolve().parent
TOPOLOGY = RUN.parents[2] / "00_system/ash274_solvated.prmtop"


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    """Write a non-overwritable water-specific scientific assessment."""
    output = RUN / "scientific_water_assessment.json"
    if output.exists():
        raise FileExistsError(output)
    process = json.loads((RUN / "process.json").read_text())
    acceptance = json.loads((RUN / "acceptance.json").read_text())
    audit = json.loads((RUN / "water_1576_full_hbond_audit_final.json").read_text())
    structure = parmed.load_file(str(TOPOLOGY), xyz=str(RUN / "optimized.ncrst"))

    def point(residue: int, name: str) -> np.ndarray:
        atoms = [
            atom for atom in structure.residues[residue - 1].atoms if atom.name == name
        ]
        if len(atoms) != 1:
            raise ValueError(f"Expected exactly one {residue}:{name}")
        return np.asarray(structure.coordinates[atoms[0].idx], dtype=float)

    def distance(first: tuple[int, str], second: tuple[int, str]) -> float:
        return float(np.linalg.norm(point(*first) - point(*second)))

    bond_lengths = {
        "water_O_H1_A": distance((1576, "O"), (1576, "H1")),
        "water_O_H2_A": distance((1576, "O"), (1576, "H2")),
        "Thr233_OG1_HG1_A": distance((929, "OG1"), (929, "HG1")),
        "ASH274_OD2_HD2_A": distance((970, "OD2"), (970, "HD2")),
    }
    endpoint = audit["measurements"]["endpoint"]
    start = audit["measurements"]["start"]
    distances = endpoint["distances_A"]
    angles = endpoint["donor_h_acceptor_angles_deg"]
    thr_to_water = (
        distances["Thr929:HG1"]["O"] <= 2.5
        and angles["Thr929:OG1-Thr929:HG1...WAT1576:O"] >= 150.0
    )
    water_to_asp272 = (
        distances["Asp272:OD1"]["H2"] <= 2.5
        and angles["WAT1576:O-H2...Asp272:OD1"] >= 150.0
    )
    gates = {
        "normal_force_converged_tier_a": process["returncode"] == 0
        and acceptance["status"] == "accepted"
        and acceptance["tier_a_b"]["tier"] == "A",
        "covalent_protons_retained": all(
            value <= 1.2 for value in bond_lengths.values()
        ),
        "thr_hg1_to_water_hbond_geometry": thr_to_water,
        "water_h2_to_asp272_hbond_geometry": water_to_asp272,
    }
    report = {
        "status": "local_geometry_accepted" if all(gates.values()) else "gate_failed",
        "gates": gates,
        "bond_lengths_A": bond_lengths,
        "start_thr_og1_to_water_O_A": start["distances_A"]["Thr929:OG1"]["O"],
        "endpoint_thr_og1_to_water_O_A": distances["Thr929:OG1"]["O"],
        "endpoint_thr_hg1_to_water_O_A": distances["Thr929:HG1"]["O"],
        "endpoint_water_H2_to_asp272_OD1_A": distances["Asp272:OD1"]["H2"],
        "endpoint_hbond_angles_deg": {
            "Thr_OG1_HG1_water_O": angles["Thr929:OG1-Thr929:HG1...WAT1576:O"],
            "water_O_H2_Asp272_OD1": angles["WAT1576:O-H2...Asp272:OD1"],
        },
        "leaving_N_to_water_O_A": distances["leaving_N"]["O"],
        "interpretation": "Local minimized geometry only; no occupancy, transfer, or activation-barrier inference.",
        "acceptance_erratum": (
            "acceptance.json water_1576.initial_og1_distance_A is inherited from the "
            "earlier pre-corrected-shell start; the actual water-QM start is reported here."
        ),
        "source_sha256": {
            str(path.resolve()): sha256(path)
            for path in (
                Path(__file__),
                TOPOLOGY,
                RUN / "process.json",
                RUN / "acceptance.json",
                RUN / "preparation.json",
                RUN / "starting.ncrst",
                RUN / "optimized.ncrst",
                RUN / "water_1576_full_hbond_audit_final.json",
            )
        },
    }
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
