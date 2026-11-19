#!/usr/bin/env python3
"""Wait a bounded time for the activation grid and select sparse references."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import time
from typing import Any

HERE = Path(__file__).resolve().parent
GRID = HERE.parent
QUEUE_PATH = HERE / "queue.json"
LOGGER = logging.getLogger(__name__)
POLL_SECONDS = 300
MAX_WAIT_SECONDS = 48 * 60 * 60


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object from disk."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return data


def row_is_terminal(row_number: int) -> bool:
    """Return whether a prepared grid row has finished or failed explicitly."""
    row_path = GRID / f"row_{row_number:03d}" / "row.json"
    status_path = GRID / f"row_{row_number:03d}" / "row_status.json"
    if not row_path.is_file():
        return False
    if not status_path.is_file():
        return False
    status = read_json(status_path).get("status")
    return status in {"completed", "failed"}


def all_rows_terminal() -> bool:
    """Return whether all five grid row workers have terminal process records."""
    return all(row_is_terminal(number) for number in range(1, 6))


def resolved_grid_path(value: str) -> Path:
    """Resolve a grid artifact path and ensure it remains under the grid root."""
    path = (GRID / value).resolve(strict=True)
    if not path.is_relative_to(GRID.resolve()):
        raise ValueError(f"Grid artifact escapes grid directory: {value}")
    return path


def validate_grid_identity(queue: dict[str, Any]) -> dict[str, Any]:
    """Require matching QM identity before consuming grid points."""
    protocol_path = GRID / "protocol.json"
    protocol = read_json(protocol_path)
    qm = protocol.get("qm", {})
    if qm.get("real_atoms") != queue["qm_real_atoms"]:
        raise ValueError("Grid QM atom count does not match reference queue")
    if qm.get("charge") != queue["qm_charge"]:
        raise ValueError("Grid QM charge does not match reference queue")
    if qm.get("mask") != queue["qm_mask"]:
        raise ValueError("Grid QM mask does not match reference queue")
    if protocol.get("method") != "Amber Sander26 native DL-FIND / DFTB3-3OB-3-1 QM/MM":
        raise ValueError("Unexpected grid method/protocol")
    return protocol


def accepted_points(queue: dict[str, Any]) -> list[dict[str, Any]]:
    """Return only completed, accepted points with verified process and restart."""
    validate_grid_identity(queue)
    point_file = GRID / "grid_points.json"
    points = read_json(point_file).get("points")
    if not isinstance(points, dict):
        raise ValueError("grid_points.json has no point map")
    accepted: list[dict[str, Any]] = []
    for key, point in points.items():
        if point.get("status") != "completed":
            continue
        process = read_json(resolved_grid_path(point["process_record"]))
        acceptance = read_json(resolved_grid_path(point["acceptance_record"]))
        if process.get("status") != "finished" or process.get("exit_status") != 0:
            continue
        if (
            acceptance.get("accepted") is not True
            or acceptance.get("force_converged") is not True
        ):
            continue
        if point.get("chemistry_class") != "retained_candidate":
            continue
        restart = resolved_grid_path(point["restart"])
        digest = sha256(restart)
        if digest != point.get("restart_sha256"):
            raise ValueError(f"Restart hash mismatch for {key}")
        if digest != acceptance.get("restart_sha256"):
            raise ValueError(f"Acceptance/restart hash mismatch for {key}")
        geometry = acceptance.get("geometry", {})
        q_actual = point.get("q_actual_A")
        attack_actual = point.get("attack_actual_A")
        if q_actual is None or attack_actual is None:
            raise ValueError(f"Missing measured reaction coordinates for {key}")
        accepted.append(
            {
                "key": key,
                "restart": str(restart),
                "restart_sha256": digest,
                "q_actual_A": float(q_actual),
                "attack_actual_A": float(attack_actual),
                "carbonyl_co_A": float(geometry["substrate_c_o_A"]),
                "scissile_cn_A": float(geometry["substrate_c_n_A"]),
                "thr_n_hg1_A": float(geometry["thr_n_hg1_A"]),
                "chemistry_class": point.get("chemistry_class"),
                "process_record": str(resolved_grid_path(point["process_record"])),
                "acceptance_record": str(
                    resolved_grid_path(point["acceptance_record"])
                ),
            }
        )
    return accepted


def select_points(points: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Select at most one activation, attack, and tetrahedral-like geometry."""
    selected: list[dict[str, Any]] = []
    activated = [
        point
        for point in points
        if point["q_actual_A"] > 0.0
        and point["attack_actual_A"] >= 2.8
        and point["thr_n_hg1_A"] <= 1.25
    ]
    if activated:
        selected.append(
            max(activated, key=lambda point: point["q_actual_A"])
            | {"class": "activated"}
        )

    attacks = [
        point
        for point in points
        if point["key"] not in {row["key"] for row in selected}
        and point["q_actual_A"] < 0.0
        and point["attack_actual_A"] <= 2.2
        and point["carbonyl_co_A"] < 1.30
        and point["scissile_cn_A"] <= 1.55
    ]
    if attacks:
        selected.append(
            min(attacks, key=lambda point: point["attack_actual_A"])
            | {"class": "attack"}
        )

    tetrahedral = [
        point
        for point in points
        if point["key"] not in {row["key"] for row in selected}
        and point["attack_actual_A"] <= 1.75
        and point["carbonyl_co_A"] >= 1.30
        and point["scissile_cn_A"] <= 1.55
        and point["thr_n_hg1_A"] <= 1.25
    ]
    if tetrahedral:
        selected.append(
            max(tetrahedral, key=lambda point: point["carbonyl_co_A"])
            | {"class": "tetrahedral_like"}
        )
    return selected[:3]


def wait_and_select() -> dict[str, Any]:
    """Wait up to 48 hours for grid rows, then write the selection manifest."""
    queue = read_json(QUEUE_PATH)
    started = time.monotonic()
    while not all_rows_terminal():
        elapsed = time.monotonic() - started
        if elapsed >= MAX_WAIT_SECONDS:
            record = {
                "status": "timed_out",
                "maximum_wait_hours": MAX_WAIT_SECONDS / 3600,
                "selected_grid_states": [],
                "reason": "Not all five grid rows reached a terminal process state within 48 hours.",
            }
            (HERE / "selection.json").write_text(json.dumps(record, indent=2) + "\n")
            return record
        LOGGER.info(
            "Grid rows still running or not prepared; elapsed %.1f h", elapsed / 3600
        )
        time.sleep(min(POLL_SECONDS, MAX_WAIT_SECONDS - elapsed))

    points = accepted_points(queue)
    selected = select_points(points)
    record = {
        "status": "grid_terminal",
        "maximum_wait_hours": MAX_WAIT_SECONDS / 3600,
        "grid_protocol_sha256": sha256(GRID / "protocol.json"),
        "grid_points_sha256": sha256(GRID / "grid_points.json"),
        "accepted_point_count": len(points),
        "selected_grid_states": selected,
        "maximum_selected_grid_states": 3,
        "reference_queue_status": "selected_pending_PC_membership_and_calculator_preparation",
    }
    (HERE / "selection.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    """Run the bounded grid wait and print a compact status message."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    result = wait_and_select()
    LOGGER.info(
        "Queue status=%s selected=%d",
        result["status"],
        len(result["selected_grid_states"]),
    )


if __name__ == "__main__":
    main()
