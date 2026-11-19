"""Prepare a bounded water-to-Asp272 TI1-like diagnostic path.

This is an alternative protonation diagnostic, not the canonical NH3+ TI1
pathway. Preparation writes immutable per-window inputs; a separate runner
handles calculations after review.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
import re
import shutil
from typing import Any

import numpy as np
import parmed
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[4]
CASE = Path(__file__).resolve().parent
SOURCE_CASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/06_qmmm_optimization"
    / "replica_1_tier_A_seed/water_qm_reoptimization"
)
TOPOLOGY = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/00_system/ash274_solvated.prmtop"
)
BASE_MDIN = SOURCE_CASE / "minimize.mdin"
START_RESTART = SOURCE_CASE / "optimized.ncrst"
SANDER = Path("/data/jpc/code/ambertools26/venv/bin/sander")
JOURNAL = ROOT / "output_agent_actions/20261021.jsonl"
LOGGER = logging.getLogger(__name__)

# Amber one-based atom IDs resolved from the validated ASH274 topology.
ATOM = {
    "thr_n": 13905,
    "thr_h1": 13906,
    "thr_h2": 13907,
    "thr_og1": 13916,
    "thr_hg1": 13917,
    "water_o": 20317,
    "water_h1": 20318,
    "water_h2": 20319,
    "asp272_od1": 14559,
    "ash274_od2": 14579,
    "ash274_hd2": 14580,
    "phe113_n": 1638,
    "asn186_nd2": 2732,
    "substrate_c": 19920,
    "substrate_o": 19921,
    "leaving_n": 19933,
    "leaving_hn": 19941,
}
WINDOWS = 25
Q_FINAL = 1.50
ATTACK_FINAL = 1.65
FORCE_CONSTANT = 100.0
DISTANCE_KEYS = (
    "thr_og1_hg1_A",
    "water_o_hg1_A",
    "water_o_h1_A",
    "water_o_h2_A",
    "asp272_od1_h2_A",
    "q_thr_to_water_A",
    "q_water_to_asp272_A",
    "relay_Q_A",
    "og1_substrate_c_A",
    "substrate_c_o_A",
    "substrate_c_n_A",
    "thr_n_h1_A",
    "thr_n_h2_A",
    "leaving_n_hn_A",
    "leaving_n_hg1_A",
    "ash274_od2_hd2_A",
    "ash274_hd2_water_o_A",
    "thr_n_hg1_A",
    "thr_n_water_h1_A",
    "thr_n_water_h2_A",
    "leaving_n_water_h1_A",
    "leaving_n_water_h2_A",
    "substrate_o_phe113_n_A",
    "substrate_o_asn186_nd2_A",
    "attack_angle_deg",
)


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_xyz(coordinates: np.ndarray, atom_id: int) -> np.ndarray:
    """Return coordinates for one-based Amber atom ID."""
    return coordinates[atom_id - 1]


def distance(coordinates: np.ndarray, first: int, second: int) -> float:
    """Measure one Euclidean distance in angstrom."""
    return float(
        np.linalg.norm(atom_xyz(coordinates, first) - atom_xyz(coordinates, second))
    )


def geometry(coordinates: np.ndarray) -> dict[str, float]:
    """Measure relay, attack, and proton-inventory distances."""
    pairs = {
        "thr_og1_hg1_A": (ATOM["thr_og1"], ATOM["thr_hg1"]),
        "water_o_hg1_A": (ATOM["water_o"], ATOM["thr_hg1"]),
        "water_o_h1_A": (ATOM["water_o"], ATOM["water_h1"]),
        "water_o_h2_A": (ATOM["water_o"], ATOM["water_h2"]),
        "asp272_od1_h2_A": (ATOM["asp272_od1"], ATOM["water_h2"]),
        "og1_substrate_c_A": (ATOM["thr_og1"], ATOM["substrate_c"]),
        "substrate_c_o_A": (ATOM["substrate_c"], ATOM["substrate_o"]),
        "substrate_c_n_A": (ATOM["substrate_c"], ATOM["leaving_n"]),
        "thr_n_h1_A": (ATOM["thr_n"], ATOM["thr_h1"]),
        "thr_n_h2_A": (ATOM["thr_n"], ATOM["thr_h2"]),
        "leaving_n_hn_A": (ATOM["leaving_n"], ATOM["leaving_hn"]),
        "leaving_n_hg1_A": (ATOM["leaving_n"], ATOM["thr_hg1"]),
        "ash274_od2_hd2_A": (ATOM["ash274_od2"], ATOM["ash274_hd2"]),
        "ash274_hd2_water_o_A": (ATOM["ash274_hd2"], ATOM["water_o"]),
        "thr_n_hg1_A": (ATOM["thr_n"], ATOM["thr_hg1"]),
        "thr_n_water_h1_A": (ATOM["thr_n"], ATOM["water_h1"]),
        "thr_n_water_h2_A": (ATOM["thr_n"], ATOM["water_h2"]),
        "leaving_n_water_h1_A": (ATOM["leaving_n"], ATOM["water_h1"]),
        "leaving_n_water_h2_A": (ATOM["leaving_n"], ATOM["water_h2"]),
        "substrate_o_phe113_n_A": (ATOM["substrate_o"], ATOM["phe113_n"]),
        "substrate_o_asn186_nd2_A": (ATOM["substrate_o"], ATOM["asn186_nd2"]),
    }
    result = {name: distance(coordinates, *pair) for name, pair in pairs.items()}
    result["q_thr_to_water_A"] = result["thr_og1_hg1_A"] - result["water_o_hg1_A"]
    result["q_water_to_asp272_A"] = result["water_o_h2_A"] - result["asp272_od1_h2_A"]
    result["relay_Q_A"] = result["q_thr_to_water_A"] + result["q_water_to_asp272_A"]
    first = atom_xyz(coordinates, ATOM["thr_og1"]) - atom_xyz(
        coordinates, ATOM["substrate_c"]
    )
    second = atom_xyz(coordinates, ATOM["substrate_o"]) - atom_xyz(
        coordinates, ATOM["substrate_c"]
    )
    cosine = float(
        np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
    )
    result["attack_angle_deg"] = float(
        np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    )
    return result


def schedule(initial: dict[str, float]) -> list[dict[str, float | str]]:
    """Return 25 smooth, asynchronous composite-Q and attack targets."""
    return [
        {
            "window": f"window_{index:03d}",
            "relay_Q_target_A": initial["relay_Q_A"]
            + index / WINDOWS * (Q_FINAL - initial["relay_Q_A"]),
            "attack_target_A": initial["og1_substrate_c_A"]
            + index / WINDOWS * (ATTACK_FINAL - initial["og1_substrate_c_A"]),
        }
        for index in range(1, WINDOWS + 1)
    ]


def restraints(relay_target: float, attack_target: float) -> str:
    """Render validated Amber four-distance-Q and separate attack restraints."""
    q_atoms = (
        ATOM["thr_og1"],
        ATOM["thr_hg1"],
        ATOM["water_o"],
        ATOM["thr_hg1"],
        ATOM["water_o"],
        ATOM["water_h2"],
        ATOM["asp272_od1"],
        ATOM["water_h2"],
    )
    return (
        f"&rst iat={','.join(map(str, q_atoms))}, rstwt=1,-1,1,-1, "
        f"r1={relay_target - 2.0:.10f}, r2={relay_target:.10f}, "
        f"r3={relay_target:.10f}, r4={relay_target + 2.0:.10f}, "
        f"rk2={FORCE_CONSTANT:g}, rk3={FORCE_CONSTANT:g}, /\n"
        f"&rst iat={ATOM['thr_og1']},{ATOM['substrate_c']},0,0,0,0,0,0, "
        f"rstwt=1,0,0,0, r1={attack_target - 2.0:.10f}, "
        f"r2={attack_target:.10f}, r3={attack_target:.10f}, "
        f"r4={attack_target + 2.0:.10f}, rk2={FORCE_CONSTANT:g}, "
        f"rk3={FORCE_CONSTANT:g}, /\n"
    )


def driven_mdin(base: str) -> str:
    """Enable Amber distance restraints without changing QM/MM or active masks."""
    if "nmropt=0" not in base.replace(" ", "") or "ntr=0" not in base.replace(" ", ""):
        raise ValueError("Source mdin is not the accepted unbiased QM/MM input")
    if "&qmmm" not in base or "&dlfind" not in base:
        raise ValueError("Source mdin lacks its audited QM/MM or DL-FIND block")
    text, count = re.subn(r"\bnmropt\s*=\s*0", "nmropt=1", base, count=1)
    if count != 1:
        raise ValueError("Expected exactly one nmropt=0 control")
    control, dlfind = text.split("&dlfind", maxsplit=1)
    return (
        control.rstrip()
        + "\n&wt type='END' /\nDISANG=restraints.RST\n"
        + "DUMPAVE=restraint_values.dat\n&dlfind"
        + dlfind
    )


def verify_source() -> tuple[Any, dict[str, Any], np.ndarray]:
    """Validate the accepted source endpoint and return its system and geometry."""
    prep = json.loads((SOURCE_CASE / "preparation.json").read_text())
    status = json.loads((SOURCE_CASE / "status.json").read_text())
    process = json.loads((SOURCE_CASE / "process.json").read_text())
    if (
        status.get("exit_status") != 0
        or not status["scientific_acceptance"]["accepted"]
    ):
        raise ValueError("Water-QM source endpoint is not accepted")
    if prep["qm_region"]["real_atoms"] != 102 or prep["qm_region"]["charge"] != -1:
        raise ValueError("Water-QM source does not match 102-atom, charge -1 QM region")
    if prep["movable_region"]["atom_count"] != 4662:
        raise ValueError(
            "Water-QM source does not match audited 4,662-atom active region"
        )
    if process["output_sha256"]["optimized.ncrst"] != sha256(START_RESTART):
        raise ValueError("Accepted source restart hash differs from its process record")
    system = parmed.load_file(str(TOPOLOGY), xyz=str(START_RESTART))
    for name, atom_id in ATOM.items():
        atom = system.atoms[atom_id - 1]
        expected = {
            "thr_n": (929, "N"),
            "thr_h1": (929, "H1"),
            "thr_h2": (929, "H2"),
            "thr_og1": (929, "OG1"),
            "thr_hg1": (929, "HG1"),
            "water_o": (1576, "O"),
            "water_h1": (1576, "H1"),
            "water_h2": (1576, "H2"),
            "asp272_od1": (968, "OD1"),
            "ash274_od2": (970, "OD2"),
            "ash274_hd2": (970, "HD2"),
            "phe113_n": (113, "N"),
            "asn186_nd2": (186, "ND2"),
            "substrate_c": (1314, "C"),
            "substrate_o": (1314, "O"),
            "leaving_n": (1315, "N"),
            "leaving_hn": (1315, "HN"),
        }[name]
        if (atom.residue.idx + 1, atom.name.strip()) != expected:
            raise ValueError(f"Atom ID {atom_id} no longer maps to {name}={expected}")
    carbon = system.atoms[ATOM["substrate_c"] - 1]
    if system.atoms[ATOM["leaving_n"] - 1] not in carbon.bond_partners:
        raise ValueError("Substrate C1314-N1315 is not bonded in the input topology")
    qmask = prep["qm_region"]["mask"]
    active = prep["movable_region"]["mask"]
    qm_indices = {
        i for i, value in enumerate(AmberMask(system, qmask).Selection()) if value
    }
    active_indices = {
        i for i, value in enumerate(AmberMask(system, active).Selection()) if value
    }
    if (
        len(qm_indices) != 102
        or len(active_indices) != 4662
        or not qm_indices <= active_indices
    ):
        raise ValueError("Audited QM/active masks changed in count or containment")
    coords = np.asarray(system.coordinates, dtype=float)
    start_geometry = geometry(coords)
    if abs(start_geometry["relay_Q_A"] - (-1.6929592019)) > 0.08:
        raise ValueError("Starting relay Q differs from the verified E-S coordinate")
    return system, prep, start_geometry


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    """Write a deterministic CSV ledger."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def prepare() -> dict[str, Any]:
    """Write the 25 window inputs and explicit scientific acceptance gates."""
    if (CASE / "protocol.json").exists() or (CASE / "window_000").exists():
        raise FileExistsError("Path preparation already exists; refusing overwrite")
    _, prep, start = verify_source()
    base = driven_mdin(BASE_MDIN.read_text())
    windows = schedule(start)
    CASE.mkdir(parents=True, exist_ok=True)
    baseline = CASE / "window_000"
    baseline.mkdir()
    shutil.copy2(START_RESTART, baseline / "optimized.ncrst")
    baseline_record = {
        "window": "window_000",
        "status": "accepted_source_endpoint",
        "geometry": start,
        "restart_sha256": sha256(baseline / "optimized.ncrst"),
        "source_restart_sha256": sha256(START_RESTART),
    }
    (baseline / "acceptance.json").write_text(
        json.dumps(baseline_record, indent=2) + "\n"
    )
    manifest_windows = []
    for item in windows:
        folder = CASE / str(item["window"])
        folder.mkdir()
        (folder / "drive.mdin").write_text(base)
        (folder / "restraints.RST").write_text(
            restraints(float(item["relay_Q_target_A"]), float(item["attack_target_A"]))
        )
        record = dict(item) | {
            "start_from": "window_000/optimized.ncrst"
            if item["window"] == "window_001"
            else f"window_{int(str(item['window']).split('_')[1]) - 1:03d}/optimized.ncrst",
            "expected_outputs": [
                "minimized.mdout",
                "optimized.ncrst",
                "optimization.nc",
                "mdinfo",
                "restraint_values.dat",
            ],
            "status": "prepared_not_run",
        }
        (folder / "window.json").write_text(json.dumps(record, indent=2) + "\n")
        manifest_windows.append(record)
    gates = {
        "label": "alternative Asp272-protonated TI1-like diagnostic; not canonical NH3+ TI1",
        "source": "accepted 102-atom, charge -1 water-QM E-S endpoint",
        "force_convergence": "Sander exit 0 and DL-FIND Convergence reached / Optimisation converged in every window",
        "individual_relay_transfers": "Thr OG1-HG1 >= 1.30 A and water O-HG1 <= 1.20 A; water O-H2 >= 1.30 A and Asp272 OD1-H2 <= 1.20 A",
        "attack_and_tetrahedral_geometry": "OG1-C1314 <= 1.75 A and C1314-O >= 1.30 A",
        "scissile_bond_retention": "C1314-N1315 <= 1.55 A throughout driven stage",
        "no_unintended_n_protonation": "Thr929 N retains H1/H2 at <=1.25 A, has no H3; substrate N1315 retains HN <=1.25 A; Thr HG1 and either water H remain >1.30 A from both N atoms",
        "water_h1_retention": "water1576 O-H1 remains <=1.25 A",
        "ash274_proton_retention": "ASH274 OD2-HD2 remains <=1.25 A and HD2 stays >1.30 A from water O",
        "oxyanion_hole": "unrestrained Phe113:N and Asn186:ND2 each stay within 3.5 A of scissile O1314; distances are monitored without restraints",
        "continuity": "flag adjacent-window changes >0.35 A in any measured reactive distance; stop if the jump indicates a topology/geometry discontinuity",
        "energy": "record physical EPtot-RESTRAINT and bias-inclusive EPtot per window; flag energy jumps for review, with no arbitrary absolute barrier cutoff",
        "endpoint_release": "only after all endpoint gates pass, prepare a separate unrestrained DFTB3/MM relaxation from the driven endpoint; classify release separately",
        "canonical_ti1_status": "not accepted as native TI1; canonical TI1 requires Thr233 beta N (Thr929:N here) to be NH3+",
        "proton_inventory": "water1576 H1 stays on water; ASH274 HD2 stays on OD2; Thr929 N remains NH2; substrate N1315 retains HN and is not protonated by either relay hydrogen",
    }
    release_dir = CASE / "release_pending"
    release_dir.mkdir()
    shutil.copy2(BASE_MDIN, release_dir / "unrestrained_relax.mdin")
    protocol = {
        "name": CASE.name,
        "status": "prepared_not_run",
        "method": "Amber Sander26 native DL-FIND / DFTB3-3OB-3-1 QM/MM",
        "qm": {"mask": prep["qm_region"]["mask"], "real_atoms": 102, "charge": -1},
        "active": {
            "mask": prep["movable_region"]["mask"],
            "atoms": 4662,
            "radius_A": 12.0,
        },
        "coordinates": {
            "Q": "d(Thr929 OG1,HG1)-d(WAT1576 O,HG1)+d(WAT1576 O,H2)-d(Asp272 OD1,H2)",
            "attack": "d(Thr929 OG1,substrate residue1314 C)",
            "starting_geometry": start,
            "final_targets": {"Q_A": Q_FINAL, "attack_A": ATTACK_FINAL},
            "window_count": WINDOWS,
            "maximum_target_steps_A": {
                "Q": abs(Q_FINAL - start["relay_Q_A"]) / WINDOWS,
                "attack": abs(ATTACK_FINAL - start["og1_substrate_c_A"]) / WINDOWS,
            },
            "force_constant_kcal_mol_A2": FORCE_CONSTANT,
        },
        "atom_ids_amber_one_based": ATOM,
        "restraint_syntax_precedent": "eight-atom/four-distance rstwt syntax verified against direct_acyl_drive.py and executed asynchronous direct-acyl path windows",
        "acceptance_gates": gates,
        "outputs": {
            "energy_ledger": "energy_ledger.csv",
            "distance_ledger": "distance_ledger.csv",
            "vmd_progress": "vmd_progress/",
        },
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sources_sha256": {
            str(path): sha256(path)
            for path in (
                TOPOLOGY,
                BASE_MDIN,
                START_RESTART,
                SOURCE_CASE / "preparation.json",
                SOURCE_CASE / "status.json",
                SOURCE_CASE / "process.json",
                Path(__file__),
                CASE / "run_water_to_asp272_path.py",
                CASE / "export_water_path.py",
                CASE / "water_path_vmd_support.py",
            )
        },
        "release_input": "release_pending/unrestrained_relax.mdin; separate release only after endpoint gates pass and parent review",
        "windows": manifest_windows,
    }
    (CASE / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    (CASE / "acceptance_criteria.json").write_text(json.dumps(gates, indent=2) + "\n")
    return protocol


def main() -> None:
    """Prepare inputs without launching any calculation."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true", required=True)
    parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    prepare()
    LOGGER.info("Prepared path inputs only; no calculation launched")


if __name__ == "__main__":
    main()
