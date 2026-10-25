"""Strict source and representative-force proof gate for ES relaxation."""

from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
from nyl12_covalent.es_dftb_force_evidence import check_force_evidence
from nyl12_covalent.es_dftb_analytic_pilot import check_analytic_pilot


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _compact_mask(indices: list[int]) -> str:
    """Express identical atom indices as ranges within Amber's mask limit."""
    if not indices or any(type(value) is not int or value < 1 for value in indices):
        raise ValueError("Positive integer atom indices required")
    if len(set(indices)) != len(indices):
        raise ValueError("Atom indices must be unique")
    ordered = sorted(indices)
    first = previous = ordered[0]
    intervals = []
    for value in ordered[1:] + [ordered[-1] + 2]:
        if value != previous + 1:
            intervals.append(str(first) if first == previous else f"{first}-{previous}")
            first = value
        previous = value
    mask = "@" + ",".join(intervals)
    if len(mask) + 3 > 255:
        raise ValueError(
            "Compressed restraint mask exceeds Amber's 256-character field"
        )
    return mask


def _force_gate(
    engine_path: Path,
    force_path: Path,
    audit: dict,
    inventory_source: Path,
    pilot_decision: Path | None = None,
) -> dict:
    engine = json.loads(engine_path.read_text())
    check = json.loads(force_path.read_text())
    approved = engine["approved_inventory_source"]
    if Path(approved["path"]).resolve() != inventory_source.resolve() or approved[
        "sha256"
    ] != _sha(inventory_source):
        raise ValueError(
            "Force proof does not pin the independently audited inventory config"
        )
    if pilot_decision is None:
        if check["accepted"] is not True or check["sources_unchanged"] is not True:
            raise ValueError(
                "Accepted and source-unchanged enzyme force check required"
            )
        proof_paths = check_force_evidence(engine, check, force_path.parent)
        baseline, decision = check["baseline"], None
    else:
        baseline, proof_paths, decision = check_analytic_pilot(
            engine, force_path, pilot_decision, inventory_source
        )
    records = json.loads((force_path.parent / "source_manifest.json").read_text())
    records += baseline["api_sources"]
    paths = {}
    for record in records:
        path = Path(record["path"]).resolve()
        if (
            not path.is_file()
            or _sha(path) != record["sha256"]
            or path.stat().st_size != record["bytes"]
        ):
            raise ValueError("Force-check source manifest mismatch")
        paths[path] = record["sha256"]
    required = [
        engine_path,
        Path(engine["topology"]),
        Path(engine["coordinates"]),
        Path(engine["python_executable"]),
        *[Path(p) for p in engine["parameter_files"]],
    ]
    if any(path.resolve() not in paths for path in required):
        raise ValueError("Required force-check sources are not pinned")
    snapshot = force_path.parent / "source/declared_config.json"
    if _sha(snapshot) != _sha(engine_path):
        raise ValueError("Force-check config snapshot differs from live engine config")
    mm, qm = engine["mm_options"], engine["qm_options"]
    expected_mm = {
        "ntb": 1,
        "igb": 0,
        "ntc": 1,
        "ntf": 1,
        "ntr": 0,
        "ifqnt": 1,
        "cut": 10.0,
        "jfastw": 4,
    }
    expected_qm = {
        "qm_theory": "DFTB3",
        "qmcharge": -2,
        "qmmm_int": 1,
        "adjust_q": 2,
        "qm_ewald": 1,
        "qm_pme": 1,
        "qmcut": 10.0,
        "scfconv": 1e-10,
        "lnk_dis": 1.09,
        "lnk_atomic_no": 1,
        "lnk_method": 1,
        "qmmask": audit["qmmask"],
    }
    if any(mm.get(key) != value for key, value in expected_mm.items()) or any(
        qm.get(key) != value for key, value in expected_qm.items()
    ):
        raise ValueError(
            "Force-check Hamiltonian differs from approved relaxation physics"
        )
    wanted = {
        "real_atoms": 98,
        "link_atoms": 8,
        "electron_count": 368,
        "boundary_pairs": sorted(
            [[row["qm_index"], row["mm_index"]] for row in audit["boundaries"]]
        ),
    }
    for observed in [engine["expected_inventory"], baseline["inventory"]]:
        for key, value in wanted.items():
            actual = sorted(observed[key]) if key == "boundary_pairs" else observed[key]
            if actual != value:
                raise ValueError(
                    "Force-check inventory does not match independent region audit"
                )
    for key, value in expected_mm.items():
        if baseline["realized_options"]["mm"][key] != value:
            raise ValueError("Realized force-check MM physics mismatch")
    for key, value in expected_qm.items():
        if baseline["realized_options"]["qm"][key] != value:
            raise ValueError("Realized force-check QM physics mismatch")
    if engine["environment"]["OMP_NUM_THREADS"] != "1":
        raise ValueError("Explicit serial Amber environment required")
    engine["_force_source_paths"] = [str(path) for path in paths] + [
        str(path) for path in proof_paths
    ]
    engine["_baseline_box"] = baseline["box"]
    engine["_baseline_energy_kcal_per_mol"] = baseline["energy_kcal_per_mol"]
    force_array_path = (
        (force_path.parent / "baseline/analytic_forces.npy")
        if decision is None
        else Path(decision["baseline_directory"]) / "analytic_forces.npy"
    )
    engine["_baseline_force_atom_count"] = len(
        np.load(force_array_path, allow_pickle=False)
    )
    engine["_pilot_decision"] = decision
    return engine
