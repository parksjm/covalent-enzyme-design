"""Prepare a bounded independent-coordinate DFTB3/MM activation grid."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np
import parmed

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / "07_water_to_asp272_ti1_like"))
from prepare_water_to_asp272_path import (  # noqa: E402
    ATOM,
    BASE_MDIN,
    FORCE_CONSTANT,
    SOURCE_CASE,
    START_RESTART,
    TOPOLOGY,
    distance,
    driven_mdin,
    sha256,
    verify_source,
)

LOGGER = logging.getLogger(__name__)
Q_MAX = 1.05
Q_STEP = 0.15
ATTACK_STEP = 0.10
Q_COLUMNS = 5


def q_n(coords: np.ndarray) -> float:
    """Measure Thr N activation coordinate d(OG1-HG1)-d(N-HG1)."""
    return distance(coords, ATOM["thr_og1"], ATOM["thr_hg1"]) - distance(
        coords, ATOM["thr_n"], ATOM["thr_hg1"]
    )


def line(start: float, stop: float, max_step: float) -> list[float]:
    """Return endpoint-inclusive linear targets with bounded increments."""
    count = max(1, int(np.ceil(abs(stop - start) / max_step)))
    return [start + (stop - start) * i / count for i in range(1, count + 1)]


def rst(first: int, second: int, target: float) -> str:
    """Format one Amber harmonic distance restraint using the established syntax."""
    return (
        f"&rst iat={first},{second},0,0,0,0,0,0, rstwt=1,0,0,0, "
        f"r1={target - 2:.10f}, r2={target:.10f}, r3={target:.10f}, "
        f"r4={target + 2:.10f}, rk2={FORCE_CONSTANT:g}, rk3={FORCE_CONSTANT:g}, /\n"
    )


def write_immutable(path: Path, content: str) -> None:
    """Write a prepared text input, allowing identical retry after partial prep."""
    if path.exists():
        if path.read_text() != content:
            raise FileExistsError(f"Prepared artifact differs from retry: {path}")
        return
    path.write_text(content)


def q_rst(target: float) -> str:
    """Format signed q_N restraint from OG1-HG1 minus N-HG1 distances."""
    atoms = (ATOM["thr_og1"], ATOM["thr_hg1"], ATOM["thr_n"], ATOM["thr_hg1"])
    return (
        f"&rst iat={','.join(map(str, atoms))},0,0,0,0, rstwt=1,-1,0,0, "
        f"r1={target - 2:.10f}, r2={target:.10f}, r3={target:.10f}, "
        f"r4={target + 2:.10f}, rk2={FORCE_CONSTANT:g}, rk3={FORCE_CONSTANT:g}, /\n"
    )


def make_plan(seed_q: float, seed_attack: float) -> dict[str, Any]:
    """Build all q-only row preparations and attack-only grid point schedules."""
    q_targets = np.linspace(seed_q, Q_MAX, Q_COLUMNS).tolist()
    attack_targets = np.linspace(seed_attack, 1.5, Q_COLUMNS).tolist()
    rows = []
    for ri, qtarget in enumerate(q_targets, 1):
        stages: list[dict[str, Any]] = []
        q_targets_incremental = line(seed_q, qtarget, Q_STEP)
        for si, target in enumerate(q_targets_incremental, 1):
            stages.append(
                {
                    "stage": f"qprep_{si:03d}",
                    "kind": "q_prep",
                    "q_target_A": target,
                    "attack_target_A": seed_attack,
                    "grid_key": None,
                }
            )
        for ci, attack in enumerate(attack_targets, 1):
            previous = seed_attack if ci == 1 else attack_targets[ci - 2]
            for ai, target in enumerate(line(previous, attack, ATTACK_STEP), 1):
                key = f"q_{ri:03d}_attack_{ci:03d}" if target == attack else None
                stages.append(
                    {
                        "stage": f"attack_{ci:03d}_{ai:03d}",
                        "kind": "grid" if key else "attack_bridge",
                        "q_target_A": qtarget,
                        "attack_target_A": target,
                        "grid_key": key,
                    }
                )
        rows.append({"row": ri, "q_target_A": qtarget, "stages": stages})
    return {
        "rows": rows,
        "attack_targets": attack_targets,
        "stage_count": sum(len(row["stages"]) for row in rows),
    }


def prepare() -> dict[str, Any]:
    """Write immutable inputs and provenance; never launch Sander."""
    if (HERE / "protocol.json").exists():
        raise FileExistsError(f"Refusing overwrite in {HERE}")
    system, source, seed_geometry = verify_source()
    coords = np.asarray(system.coordinates, dtype=float)
    seed_q, seed_attack = q_n(coords), seed_geometry["og1_substrate_c_A"]
    plan = make_plan(seed_q, seed_attack)
    if plan["stage_count"] > 150:
        raise ValueError(f"Grid requires {plan['stage_count']} stages; maximum is 150")
    mdin = driven_mdin(BASE_MDIN.read_text())
    rows_out = []
    grid_points: dict[str, Any] = {}
    for row in plan["rows"]:
        folder = HERE / f"row_{row['row']:03d}"
        folder.mkdir(parents=True, exist_ok=True)
        if (folder / "starting.ncrst").exists():
            if sha256(folder / "starting.ncrst") != sha256(START_RESTART):
                raise ValueError(
                    f"Partial preparation restart differs from source: {folder}"
                )
        else:
            shutil.copy2(START_RESTART, folder / "starting.ncrst")
        stages = []
        for stage in row["stages"]:
            d = folder / stage["stage"]
            d.mkdir(exist_ok=True)
            write_immutable(d / "drive.mdin", mdin)
            rest = q_rst(stage["q_target_A"]) + rst(
                ATOM["thr_og1"], ATOM["substrate_c"], stage["attack_target_A"]
            )
            write_immutable(d / "restraints.RST", rest)
            item = stage | {
                "directory": str(d.relative_to(HERE)),
                "status": "prepared_not_run",
                "input_sha256": {
                    "drive.mdin": sha256(d / "drive.mdin"),
                    "restraints.RST": sha256(d / "restraints.RST"),
                },
            }
            stages.append(item)
            if stage["grid_key"]:
                grid_points[stage["grid_key"]] = {
                    "status": "prepared_not_run",
                    "q_target_A": stage["q_target_A"],
                    "attack_target_A": stage["attack_target_A"],
                    "q_actual_A": None,
                    "attack_actual_A": None,
                    "physical_energy_kcal_mol": None,
                    "bias_inclusive_energy_kcal_mol": None,
                    "restraint_energy_kcal_mol": None,
                    "restart": f"{item['directory']}/optimized.ncrst",
                    "restart_sha256": None,
                    "chemistry_class": "unclassified",
                    "process_record": f"{item['directory']}/process.json",
                    "acceptance_record": f"{item['directory']}/acceptance.json",
                }
        write_immutable(
            folder / "row.json",
            json.dumps({"row": row["row"], "stages": stages}, indent=2) + "\n",
        )
        rows_out.append(
            {"row": row["row"], "q_target_A": row["q_target_A"], "stages": len(stages)}
        )
    snapshot_dir = HERE / "source_snapshot_attempt_2"
    snapshot_dir.mkdir(exist_ok=True)
    prior_scripts = HERE.parent / "07_water_to_asp272_ti1_like"
    snapshots = [
        prior_scripts / name
        for name in (
            "prepare_water_to_asp272_path.py",
            "run_water_to_asp272_path.py",
            "export_water_path.py",
            "water_path_vmd_support.py",
        )
    ]
    for item in snapshots:
        target = snapshot_dir / item.name
        if target.exists() and sha256(target) != sha256(item):
            raise ValueError(f"Source snapshot differs from retry: {target}")
        if not target.exists():
            shutil.copy2(item, target)
    for source_path, target in (
        (Path(__file__), snapshot_dir / "prepare_grid.py"),
        (HERE / "run_row.py", snapshot_dir / "run_row.py"),
    ):
        if target.exists() and sha256(target) != sha256(source_path):
            raise ValueError(f"Source snapshot differs from retry: {target}")
        if not target.exists():
            shutil.copy2(source_path, target)
    sources = [
        TOPOLOGY,
        BASE_MDIN,
        START_RESTART,
        SOURCE_CASE / "preparation.json",
        SOURCE_CASE / "status.json",
        SOURCE_CASE / "process.json",
        *snapshots,
        Path(__file__),
        HERE / "run_row.py",
    ]
    qm_indices = (
        np.flatnonzero(
            parmed.amber.mask.AmberMask(system, source["qm_region"]["mask"]).Selection()
        )
        + 1
    )
    qm_atoms = [
        {"amber_atom_id": int(i), "element": system.atoms[int(i) - 1].element_name}
        for i in qm_indices
    ]
    qm_set = set(int(i) for i in qm_indices)
    pc_records = [
        {
            "amber_atom_id": atom.idx + 1,
            "charge_e": float(atom.charge),
            "xyz_A": [float(v) for v in coords[atom.idx]],
        }
        for atom in system.atoms
        if atom.idx + 1 not in qm_set
    ]
    pc_payload = {
        "description": "Full-system MM atom coordinates and topology charges from accepted source; raw Amber charges, not Sander adjust_q/PQEq embedding charges.",
        "source_restart_sha256": sha256(START_RESTART),
        "source_topology_sha256": sha256(TOPOLOGY),
        "qm_real_atom_ids": sorted(qm_set),
        "mm_atom_count": len(pc_records),
        "atoms": pc_records,
    }
    pc_path = HERE / "raw_mm_membership.json"
    write_immutable(pc_path, json.dumps(pc_payload, separators=(",", ":")) + "\n")
    protocol = {
        "schema_version": 1,
        "name": HERE.name,
        "status": "prepared_not_run",
        "method": "Amber Sander26 native DL-FIND / DFTB3-3OB-3-1 QM/MM",
        "coordinates": {
            "q_N": "d(Thr929 OG1,HG1)-d(Thr929 N,HG1)",
            "attack": "d(Thr929 OG1,substrate residue1314 C)",
            "seed_q_actual_A": seed_q,
            "seed_attack_actual_A": seed_attack,
            "q_target_actual_points_A": [r["q_target_A"] for r in plan["rows"]],
            "attack_target_points_A": plan["attack_targets"],
            "force_constant_kcal_mol_A2": FORCE_CONSTANT,
            "restraint_convention": "Amber signed distance difference (q_N) and separate OG1-C distance; harmonic rk2/rk3=100",
        },
        "qm": {"mask": source["qm_region"]["mask"], "real_atoms": 102, "charge": -1},
        "qm_atoms_ordered": qm_atoms,
        "link_pairs": source["qm_region"]["link_cuts"],
        "raw_mm_membership": {
            "path": pc_path.name,
            "sha256": sha256(pc_path),
            "interpretation": "reproducible raw MM source record; not adjusted native Sander embedding charges",
        },
        "active": {
            "mask": source["movable_region"]["mask"],
            "atoms": 4662,
            "radius_A": 12.0,
        },
        "microstate": {
            "ASP272": "ASP",
            "ASP274": "ASH(HD2-OD2)",
            "Thr929_N": "neutral NH2",
            "water1576": "QM",
        },
        "stage_count": plan["stage_count"],
        "maximum_stage_count": 150,
        "parallelism": "five independent row workers, one CPU each; do not launch until approved",
        "acceptance": {
            "process": "Sander exit 0, SCC converged, DL-FIND convergence marker",
            "chemistry": "record each covalent proton inventory; identity changes classified as alternate chemistry and not silently merged",
            "continuity": "retain all observed distances; flag jumps >0.35 A for review; do not use as automatic rejection criterion",
        },
        "source_sha256": {str(p): sha256(p) for p in sources},
        "source_snapshot": {
            p.name: {"path": str(p.relative_to(HERE)), "sha256": sha256(p)}
            for p in snapshot_dir.iterdir()
            if p.is_file()
        },
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "rows": rows_out,
    }
    write_immutable(HERE / "protocol.json", json.dumps(protocol, indent=2) + "\n")
    write_immutable(
        HERE / "grid_points.json",
        json.dumps({"schema_version": 1, "points": grid_points}, indent=2) + "\n",
    )
    with (HERE / "grid_targets.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["row", "q_target_A", "attack_target_A"]
        )
        writer.writeheader()
        writer.writerows(
            {"row": r["row"], "q_target_A": r["q_target_A"], "attack_target_A": a}
            for r in plan["rows"]
            for a in plan["attack_targets"]
        )
    return protocol


def main() -> None:
    """Prepare inputs only; no engine is invoked."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true", required=True)
    parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    result = prepare()
    LOGGER.info(
        "Prepared %d stages in %d rows; no calculation launched",
        result["stage_count"],
        len(result["rows"]),
    )


if __name__ == "__main__":
    main()
