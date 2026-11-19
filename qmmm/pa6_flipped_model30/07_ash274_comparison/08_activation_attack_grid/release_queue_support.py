"""Read-only selection and gate checks for candidate release preparation."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import parmed

GRID = Path(__file__).resolve().parent
ROOT = GRID.parents[3]
PROTOCOL = GRID / "protocol.json"
POINTS = GRID / "grid_points.json"
QUEUE_GATE = GRID / "reference_checks/queue_status.json"
STATUS_PATH = GRID / "release_queue_status.json"
TIMEOUT_SECONDS = 48 * 60 * 60
POLL_SECONDS = 120
TOPOLOGY = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/00_system/ash274_solvated.prmtop"
)
ATOM = {
    "thr_n": 13905,
    "thr_h1": 13906,
    "thr_h2": 13907,
    "thr_og1": 13916,
    "thr_hg1": 13917,
    "water_o": 20317,
    "water_h1": 20318,
    "ash274_od2": 14579,
    "ash274_hd2": 14580,
    "substrate_c": 19920,
    "substrate_o": 19921,
    "leaving_n": 19933,
    "leaving_hn": 19941,
}


def digest(path: Path) -> str:
    """Return SHA-256 for one file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    """Load one JSON object and reject non-object documents."""
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def geometry(restart: Path) -> dict[str, float]:
    """Measure reaction distances and retained proton inventory."""
    system = parmed.load_file(str(TOPOLOGY), xyz=str(restart))
    xyz = np.asarray(system.coordinates, dtype=float)

    def dist(first: str, second: str) -> float:
        return float(np.linalg.norm(xyz[ATOM[first] - 1] - xyz[ATOM[second] - 1]))

    return {
        "q_N_A": dist("thr_og1", "thr_hg1") - dist("thr_n", "thr_hg1"),
        "attack_A": dist("thr_og1", "substrate_c"),
        "thr_n_hg1_A": dist("thr_n", "thr_hg1"),
        "thr_n_h1_A": dist("thr_n", "thr_h1"),
        "thr_n_h2_A": dist("thr_n", "thr_h2"),
        "ash274_od2_hd2_A": dist("ash274_od2", "ash274_hd2"),
        "ash274_hd2_water_o_A": dist("ash274_hd2", "water_o"),
        "water_o_h1_A": dist("water_o", "water_h1"),
        "leaving_n_hn_A": dist("leaving_n", "leaving_hn"),
        "substrate_c_o_A": dist("substrate_c", "substrate_o"),
        "substrate_c_n_A": dist("substrate_c", "leaving_n"),
    }


def inventory_retained(values: dict[str, float]) -> bool:
    """Require Thr N H1/H2, ASH274 HD2 and other declared protons retained."""
    return (
        values["thr_n_h1_A"] <= 1.25
        and values["thr_n_h2_A"] <= 1.25
        and values["ash274_od2_hd2_A"] <= 1.25
        and values["ash274_hd2_water_o_A"] >= 1.30
        and values["water_o_h1_A"] <= 1.25
        and values["leaving_n_hn_A"] <= 1.25
    )


def candidate_classes(values: dict[str, float]) -> list[str]:
    """Return declared activated and TI-like classes for one geometry."""
    classes = []
    if values["q_N_A"] > 0.5 and values["thr_n_hg1_A"] < 1.25:
        classes.append("activated")
    if (
        values["thr_n_hg1_A"] < 1.25
        and values["attack_A"] <= 1.75
        and values["substrate_c_o_A"] >= 1.30
        and values["substrate_c_n_A"] <= 1.55
    ):
        classes.append("ti_like")
    return classes


def rows_terminal() -> tuple[bool, dict[str, Any]]:
    """Read all five row statuses and identify terminal completion or failure."""
    rows: dict[str, Any] = {}
    for row in range(1, 6):
        path = GRID / f"row_{row:03d}/row_status.json"
        if not path.is_file():
            return False, {"missing_row_status": str(path)}
        record = load_json(path)
        rows[str(row)] = record
        if record.get("status") not in {"completed", "failed"}:
            return False, rows
    return True, rows


def validate_reference_gate() -> tuple[bool, str, dict[str, Any] | None]:
    """Accept only the exact completed, root-reviewed qualitative pass gate."""
    if not QUEUE_GATE.is_file():
        return False, "reference_queue_missing", None
    gate = load_json(QUEUE_GATE)
    method = gate.get("method_comparison", {})
    states = gate.get("required_state_count", {})
    if gate.get("schema_version") != 1 or gate.get("status") != "complete":
        return False, "reference_queue_not_complete", gate
    if (
        method.get("status") in {"fail", "inconclusive"}
        or method.get("qualitative_agreement") is False
    ):
        return False, "method_comparison_failed_or_inconclusive", gate
    if method.get("status") != "pass" or method.get("reviewed_by") != "root":
        return False, "method_comparison_not_root_reviewed_pass", gate
    if method.get("qualitative_agreement") is not True:
        return False, "qualitative_agreement_not_pass", gate
    if not isinstance(states.get("expected"), int) or not 1 <= states["expected"] <= 4:
        return False, "invalid_reference_state_count", gate
    if (
        not isinstance(states.get("completed"), int)
        or states["completed"] < states["expected"]
    ):
        return False, "reference_states_incomplete", gate
    if gate.get("source_grid_protocol_sha256") != digest(PROTOCOL):
        return False, "grid_protocol_hash_mismatch", gate
    pc_manifest = gate.get("fixed_pc_manifest_path")
    if not isinstance(pc_manifest, str) or not pc_manifest:
        return False, "adjusted_fixed_pc_manifest_path_missing", gate
    pc_path = Path(pc_manifest)
    if not pc_path.is_absolute():
        pc_path = GRID / pc_path
    if not pc_path.is_file() or gate.get("fixed_pc_manifest_sha256") != digest(pc_path):
        return False, "fixed_pc_manifest_hash_mismatch", gate
    raw = Path(method.get("raw_results", ""))
    if not raw.is_absolute():
        raw = GRID / raw
    if not raw.is_file():
        return False, "reference_raw_results_missing", gate
    records = gate.get("process_records")
    if not isinstance(records, list) or not records:
        return False, "reference_process_records_missing", gate
    for item in records:
        path = Path(item)
        if not path.is_absolute():
            path = GRID / path
        if not path.is_file():
            return False, f"reference_process_record_missing:{item}", gate
    return True, "pass", gate


def select_candidates() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Select lowest-energy unique restart in each retained chemistry class."""
    points = load_json(POINTS).get("points")
    if not isinstance(points, dict) or len(points) != 25:
        raise ValueError("Expected exactly 25 grid points")
    eligible = []
    for key, point in points.items():
        if point.get("status") != "completed":
            continue
        acceptance_path = GRID / point["acceptance_record"]
        restart = GRID / point["restart"]
        acceptance = load_json(acceptance_path)
        restart_hash = point.get("restart_sha256")
        process_path = GRID / point["process_record"]
        process = load_json(process_path)
        if (
            acceptance.get("accepted") is not True
            or acceptance.get("force_converged") is not True
            or acceptance.get("scc_converged") is not True
            or process.get("status") != "finished"
            or process.get("exit_status") != 0
            or process.get("force_converged") is not True
            or process.get("scc_converged") is not True
            or not restart.is_file()
            or digest(restart) != restart_hash
            or acceptance.get("restart_sha256") != restart_hash
        ):
            continue
        values = acceptance.get("geometry", {})
        if not inventory_retained(values):
            continue
        energy = acceptance.get("energies", {}).get("physical_energy_kcal_mol")
        if not isinstance(energy, (int, float)) or not math.isfinite(energy):
            continue
        eligible.append(
            {
                "key": key,
                "point": point,
                "restart_path": restart,
                "acceptance_path": acceptance_path,
                "geometry": values,
                "energy": float(energy),
                "classes": candidate_classes(values),
            }
        )
    candidates = []
    seen = set()
    decisions = {"classes": {}}
    for name in ("activated", "ti_like"):
        ranked = sorted(
            (x for x in eligible if name in x["classes"]),
            key=lambda x: (x["energy"], x["key"]),
        )
        if not ranked:
            decisions["classes"][name] = {
                "status": "omitted",
                "reason": "no accepted retained-chemistry geometry",
            }
            continue
        candidate = ranked[0]
        restart_hash = candidate["point"]["restart_sha256"]
        if restart_hash in seen:
            decisions["classes"][name] = {
                "status": "omitted",
                "reason": "lowest-energy candidate duplicates selected restart",
                "point": candidate["key"],
                "restart_sha256": restart_hash,
            }
            continue
        seen.add(restart_hash)
        candidate["class_name"] = name
        candidates.append(candidate)
        decisions["classes"][name] = {
            "status": "selected",
            "point": candidate["key"],
            "physical_energy_kcal_mol": candidate["energy"],
            "restart_sha256": restart_hash,
            "class_match_count": len(ranked),
        }
    decisions["selected_count"] = len(candidates)
    return candidates, decisions


def wait_until_ready() -> tuple[bool, str, dict[str, Any]]:
    """Wait at 120-second intervals for grid and reference gate, up to 48 hours."""
    deadline = time.monotonic() + TIMEOUT_SECONDS
    rows: dict[str, Any] = {}
    while True:
        terminal, rows = rows_terminal()
        if terminal:
            points = load_json(POINTS).get("points", {})
            states = {p.get("status") for p in points.values()}
            if len(points) != 25 or states != {"completed"}:
                return (
                    False,
                    "grid_rows_terminal_but_grid_points_incomplete",
                    {"rows": rows, "point_statuses": sorted(states)},
                )
            ready, reason, gate = validate_reference_gate()
            if ready:
                return True, "pass", {"rows": rows, "gate": gate}
            method = gate.get("method_comparison", {}) if gate else {}
            if gate is not None and (
                gate.get("status") == "failed"
                or method.get("status") in {"fail", "inconclusive"}
                or method.get("qualitative_agreement") is False
            ):
                return False, reason, {"rows": rows, "gate": gate}
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return False, "bounded_wait_timeout", {"rows": rows}
        time.sleep(min(POLL_SECONDS, remaining))
