"""Audit real QM hydrogen occupancy across completed activation-grid points."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import parmed

LOGGER = logging.getLogger("proton_occupancy_audit")
SCRIPT_DIR = Path(__file__).resolve().parent
HERE = (
    SCRIPT_DIR.parent
    if SCRIPT_DIR.name == "source_snapshot"
    or SCRIPT_DIR.name.startswith("source_snapshot_")
    else SCRIPT_DIR
)
GRID = HERE.parent
TOPOLOGY = GRID.parent / "00_system" / "ash274_solvated.prmtop"
BASELINE = -450493.0781
R_COV = 1.25
R_SHARED = 1.35
QM_ATOMS = [
    1621,
    1622,
    1623,
    1624,
    1625,
    1626,
    1627,
    1628,
    1629,
    1630,
    1631,
    1632,
    1633,
    1634,
    1635,
    1636,
    1637,
    1638,
    1639,
    1640,
    1641,
    1642,
    1643,
    1644,
    1645,
    1646,
    1647,
    1648,
    1649,
    1650,
    1651,
    1652,
    1653,
    1654,
    1655,
    2727,
    2728,
    2729,
    2730,
    2731,
    2732,
    2733,
    2734,
    *range(13905, 13932),
    *range(14555, 14561),
    *range(14574, 14581),
    *range(19918, 19922),
    *range(19929, 19936),
    *range(19941, 19946),
    20317,
    20318,
    20319,
]
INTEREST = {
    "Thr929_N": 13905,
    "Thr929_OG1": 13916,
    "water1576_O": 20317,
    "Asp968_OD1": 14559,
    "Asp968_OD2": 14560,
    "ASH970_OD1": 14578,
    "ASH970_OD2": 14579,
    "substrate_N": 19933,
    "substrate_O": 19921,
}


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def occupancy(hydrogens: np.ndarray, heavy: np.ndarray) -> np.ndarray:
    """Compute smooth coordination counts independent of hydrogen ordering."""
    distances = np.linalg.norm(hydrogens[:, None, :] - heavy[None, :, :], axis=2)
    return np.sum(1.0 / (1.0 + (distances / R_COV) ** 6), axis=0)


def check_permutation_invariance() -> None:
    """Assert hydrogen coordinate permutation leaves coordination unchanged."""
    hs = np.array([[0.0, 0, 0], [1.1, 0, 0], [0, 2.0, 0]])
    heavy = np.array([[0.1, 0, 0], [3.0, 0, 0]])
    if not np.allclose(
        occupancy(hs, heavy), occupancy(hs[[2, 0, 1]], heavy), atol=1e-14
    ):
        raise AssertionError("Occupancy is sensitive to hydrogen ordering")


def atom_element(atom: Any) -> str:
    """Get normalized element symbol from a ParmEd atom."""
    return str(atom.element_name).strip().upper()


def main() -> None:
    """Validate provenance, assign QM protons, and write tables and plot."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    check_permutation_invariance()
    topology = parmed.load_file(str(TOPOLOGY))
    grid = json.loads((GRID / "grid_points.json").read_text())
    protocol = json.loads((GRID / "protocol.json").read_text())
    if protocol["qm"]["real_atoms"] != 102 or len(QM_ATOMS) != 102:
        raise ValueError("Protocol real-atom count or explicit QM map mismatch")
    declared_map = [item["amber_atom_id"] for item in protocol["qm_atoms_ordered"]]
    if declared_map != QM_ATOMS:
        raise ValueError("Explicit QM atom map differs from protocol ordering")
    qm_indices = np.array(QM_ATOMS, dtype=int) - 1
    qm_h = np.array(
        [i for i in qm_indices if atom_element(topology.atoms[int(i)]) == "H"]
    )
    qm_heavy = np.array(
        [i for i in qm_indices if atom_element(topology.atoms[int(i)]) != "H"]
    )
    all_heavy = np.array(
        [i for i, atom in enumerate(topology.atoms) if atom_element(atom) != "H"]
    )
    records: list[dict[str, Any]] = []
    hydrogen_audit: list[dict[str, Any]] = []
    points = grid["points"]
    completed = [
        (key, value)
        for key, value in points.items()
        if value.get("status") == "completed"
    ]
    if len(completed) != 25:
        raise ValueError(f"Expected 25 completed grid points, found {len(completed)}")
    for key, point in completed:
        restart = GRID / point["restart"]
        if sha256(restart) != point["restart_sha256"]:
            raise ValueError(f"Restart digest mismatch: {key}")
        acceptance_path = GRID / point["acceptance_record"]
        acceptance = json.loads(acceptance_path.read_text())
        process_path = GRID / point["process_record"]
        process = json.loads(process_path.read_text()) if process_path.is_file() else {}
        if (
            not process
            or not acceptance
            or acceptance.get("restart_sha256") != point["restart_sha256"]
        ):
            raise ValueError(f"Missing or inconsistent source acceptance record: {key}")
        coords = np.asarray(parmed.load_file(str(restart)).coordinates[0]).reshape(
            (-1, 3)
        )
        if len(coords) != len(topology.atoms):
            raise ValueError(f"Restart/topology atom count mismatch at {key}")
        hxyz, qheavy, systemheavy = coords[qm_h], coords[qm_heavy], coords[all_heavy]
        nearest = np.linalg.norm(hxyz[:, None, :] - qheavy[None, :, :], axis=2)
        sorted_idx = np.argsort(nearest, axis=1)[:, :2]
        system_dist = np.linalg.norm(hxyz[:, None, :] - systemheavy[None, :, :], axis=2)
        for row, hydrogen_idx in enumerate(qm_h):
            d1, d2 = nearest[row, sorted_idx[row]]
            first, second = qm_heavy[sorted_idx[row]] + 1
            first_atom = topology.atoms[int(first) - 1]
            second_atom = topology.atoms[int(second) - 1]
            destination = (
                "shared_ambiguous"
                if d1 <= R_SHARED and d2 <= R_SHARED
                else ("covalent" if d1 <= R_COV else "unassigned")
            )
            if d1 <= R_COV and d2 <= R_SHARED:
                destination = "shared_ambiguous"
            full_arg = int(np.argmin(system_dist[row]))
            full_atom = topology.atoms[int(all_heavy[full_arg])]
            escaped = (
                int(all_heavy[full_arg]) not in set(qm_heavy.tolist())
                and system_dist[row, full_arg] + 1e-8 < d1
            )
            hydrogen_audit.append(
                {
                    "point": key,
                    "hydrogen_atom_id": int(hydrogen_idx + 1),
                    "nearest_qm_heavy_atom_id": int(first),
                    "nearest_qm_heavy_name": first_atom.name,
                    "nearest_qm_heavy_residue": f"{first_atom.residue.name}{first_atom.residue.idx + 1}",
                    "nearest_qm_distance_A": float(d1),
                    "second_qm_heavy_atom_id": int(second),
                    "second_qm_heavy_name": second_atom.name,
                    "second_qm_heavy_residue": f"{second_atom.residue.name}{second_atom.residue.idx + 1}",
                    "second_qm_distance_A": float(d2),
                    "classification": destination,
                    "nearest_full_system_heavy_atom_id": int(all_heavy[full_arg] + 1),
                    "nearest_full_system_heavy_name": full_atom.name,
                    "nearest_full_system_distance_A": float(system_dist[row, full_arg]),
                    "boundary_escape": bool(escaped),
                }
            )
        cn = occupancy(hxyz, coords[qm_heavy])
        by_id = {int(idx + 1): float(cn[n]) for n, idx in enumerate(qm_heavy)}
        hg1_audit = next(
            item
            for item in hydrogen_audit
            if item["point"] == key and item["hydrogen_atom_id"] == 13917
        )
        row: dict[str, Any] = {
            "point": key,
            "q_actual_A": point["q_actual_A"],
            "attack_actual_A": point["attack_actual_A"],
            "physical_energy_kcal_mol": point["physical_energy_kcal_mol"],
            "relative_energy_kcal_mol": point["physical_energy_kcal_mol"] - BASELINE,
            "restraint_energy_kcal_mol": point["restraint_energy_kcal_mol"],
            "accepted": point.get("accepted"),
            "acceptance_chemistry_class": acceptance.get("chemistry_class"),
            "process_accepted": acceptance.get("process_accepted"),
            "mapped_Thr_HG1_destination_atom_id": hg1_audit["nearest_qm_heavy_atom_id"],
            "mapped_Thr_HG1_destination": f"{hg1_audit['nearest_qm_heavy_name']}:{hg1_audit['nearest_qm_heavy_residue']}",
            "mapped_Thr_HG1_nearest_distance_A": hg1_audit["nearest_qm_distance_A"],
            "restart_sha256": point["restart_sha256"],
            "acceptance_record_sha256": sha256(acceptance_path),
            "process_record_sha256": sha256(process_path),
        }
        for label, atom_id in INTEREST.items():
            row[f"CN_{label}"] = by_id[atom_id]
        row["CN_N_minus_OG1"] = (
            by_id[INTEREST["Thr929_N"]] - by_id[INTEREST["Thr929_OG1"]]
        )
        records.append(row)
    fields = list(records[0])
    with (HERE / "occupancy.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    with (HERE / "hydrogen_destinations.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(hydrogen_audit[0]))
        writer.writeheader()
        writer.writerows(hydrogen_audit)
    metadata = {
        "thresholds_A": {
            "covalent": R_COV,
            "shared_or_ambiguous_second_neighbor": R_SHARED,
        },
        "classification": "shared_ambiguous if two nearest QM heavy atoms are <=1.35 A, or nearest <=1.25 A and second <=1.35 A; else covalent if nearest <=1.25 A; else unassigned",
        "QM_atom_ids_1based": QM_ATOMS,
        "QM_real_atom_count": 102,
        "QM_hydrogen_count": len(qm_h),
        "excluded_fictitious_link_atoms": 8,
        "mapped_catalytic_hydrogen": {
            "amber_atom_id_1based": 13917,
            "name": topology.atoms[13916].name,
            "residue_name": topology.atoms[13916].residue.name,
            "amber_residue_number_1based": topology.atoms[13916].residue.idx + 1,
            "parmed_residue_idx_0based": topology.atoms[13916].residue.idx,
        },
        "residue_label_convention": "Amber one-based residue numbering = ParmEd residue.idx + 1",
        "atom_map": {
            k: {
                "amber_atom_id_1based": v,
                "name": topology.atoms[v - 1].name,
                "residue_name": topology.atoms[v - 1].residue.name,
                "amber_residue_number_1based": topology.atoms[v - 1].residue.idx + 1,
                "parmed_residue_idx_0based": topology.atoms[v - 1].residue.idx,
            }
            for k, v in INTEREST.items()
        },
        "topology_sha256": sha256(TOPOLOGY),
        "protocol_sha256": sha256(GRID / "protocol.json"),
        "grid_points_sha256": sha256(GRID / "grid_points.json"),
        "energy_baseline_kcal_mol": BASELINE,
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "smooth_coordination_formula": "sum_H 1 / (1 + (r_H,heavy / 1.25 A)^6)",
        "energy_field": "physical_energy_kcal_mol (native ES minus restraint); relative to common baseline",
        "accepted_point_count": sum(bool(r["accepted"]) for r in records),
        "rejected_completed_points": [r["point"] for r in records if not r["accepted"]],
        "boundary_escape_count": sum(r["boundary_escape"] for r in hydrogen_audit),
        "hydrogen_permutation_test": "passed",
        "interpretation_limit": "Occupancy/counts are geometric descriptors; no pKa or proton charge inference.",
    }
    (HERE / "audit.json").write_text(
        json.dumps(
            {
                "metadata": metadata,
                "points": records,
                "hydrogen_destinations": hydrogen_audit,
            },
            indent=2,
        )
        + "\n"
    )
    fig, ax = plt.subplots(figsize=(9, 6))
    scatter = ax.scatter(
        [r["q_actual_A"] for r in records],
        [r["attack_actual_A"] for r in records],
        c=[r["relative_energy_kcal_mol"] for r in records],
        s=85,
        cmap="viridis",
        edgecolor="black",
    )
    for r in records:
        ax.annotate(
            f"{r['point']}\nHG1→{r['mapped_Thr_HG1_destination']}",
            (r["q_actual_A"], r["attack_actual_A"]),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=6,
        )
    ax.set(
        xlabel="q = d(OG1, HG1) - d(N, HG1) / Å",
        ylabel="Attack distance OG1–substrate C / Å",
        title="Activation/attack grid: energy and QM hydrogen occupancy",
    )
    ax.set_xlim(
        min(r["q_actual_A"] for r in records) - 0.2,
        max(r["q_actual_A"] for r in records) + 0.65,
    )
    fig.colorbar(
        scatter, ax=ax, label="Physical energy relative to baseline / kcal mol$^{-1}$"
    )
    fig.tight_layout()
    fig.savefig(HERE / "occupancy_energy.png", dpi=180)
    plt.close(fig)
    LOGGER.info(
        "Wrote occupancy audit for %d points and %d real QM hydrogens each",
        len(records),
        len(qm_h),
    )


if __name__ == "__main__":
    main()
