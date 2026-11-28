"""Nonblocking orchestration of optional native-prediction confidence scoring."""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from .ipsae2_scoring import (
    file_sha256,
    require_qualification,
    score_model_raw,
    verify_runtime,
)

LOGGER = logging.getLogger(__name__)


def read_json(path: str | Path) -> dict[str, Any]:
    """Read an explicit JSON object."""
    result = json.loads(Path(path).read_text())
    if not isinstance(result, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return result


def _journal(config: dict[str, Any], phase: str, action: str) -> None:
    subprocess.run(
        [
            *config["journal_command"],
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Optional ipSAE2 sidecar; prediction and Amber jobs remain independent",
            "--rollback",
            "Retain outputs; stop only this sidecar to disable optional scoring",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def poll_cases(config: dict[str, Any], *, plan: bool = False) -> dict[str, Any]:
    """Poll once, retaining pending cases and isolating each scoring failure.

    A journal failure aborts dependent mutations. Scientific or scoring failures
    are recorded per case and cannot signal or stop prediction/Amber processes.
    """
    for path, expected in config.get("input_sha256", {}).items():
        if file_sha256(path) != expected:
            raise ValueError(f"Sidecar input hash mismatch: {path}")
    runtime = read_json(config["runtime"])
    verify_runtime(runtime)
    qualification_path = Path(config["qualification"])
    if not qualification_path.exists():
        return {"status": "pending_qualification", "cases": {}}
    require_qualification(read_json(qualification_path), runtime)
    states: dict[str, Any] = {}
    for case in config["cases"]:
        case_id = case["case_id"]
        if Path(case_id).name != case_id or case_id in {".", ".."}:
            raise ValueError("case_id must be a single descriptive directory name")
        output = Path(config["scoring_root"]) / case_id
        failure_path = (
            Path(config["scoring_root"]) / "case_failures" / f"{case_id}.json"
        )
        if failure_path.exists():
            states[case_id] = read_json(failure_path)
            continue
        completed = output / "scores.json"
        if completed.exists():
            states[case_id] = {"status": read_json(completed)["status"]}
            continue
        native_path = Path(case["native_receipt"])
        if not native_path.exists():
            states[case_id] = {"status": "pending_native_prediction"}
            continue
        native = read_json(native_path)
        native_status = native.get("status", "")
        if "failed" in native_status or native.get("exit_status", 0) != 0:
            states[case_id] = {
                "status": "failed",
                "error": f"Native prediction failed: {native_status}",
            }
            continue
        if native_status != "native_prediction_completed":
            states[case_id] = {
                "status": "pending_native_prediction",
                "native_status": native.get("status"),
            }
            continue
        mapping_path = Path(case["mapping_receipt"])
        if plan:
            states[case_id] = {
                "status": "ready" if mapping_path.exists() else "pending_mapping",
                "output": str(output),
            }
            continue
        if not mapping_path.exists() and "mapping_command" not in case:
            states[case_id] = {"status": "pending_mapping"}
            continue
        # The journal is deliberately outside the scientific-failure handler.
        _journal(config, "intent", f"Map and score ipSAE2 case {case_id}")
        try:
            if (
                "native_receipt_sha256" in case
                and file_sha256(native_path) != case["native_receipt_sha256"]
            ):
                raise ValueError("Native completion receipt hash mismatch")
            model = read_json(case["model_manifest"])
            if model["case_id"] != case_id or native["case"] != case_id:
                raise ValueError("Native receipt/model case identity mismatch")
            for artifact in model["artifacts"].values():
                if file_sha256(artifact["path"]) != artifact["sha256"]:
                    raise ValueError(
                        f"Native artifact hash mismatch: {artifact['path']}"
                    )
            if not mapping_path.exists():
                process = subprocess.run(
                    case["mapping_command"], capture_output=True, text=True, check=False
                )
                # Mapping tools retain their own receipt; retain command logs separately.
                log_root = Path(config["scoring_root"]) / "mapping_processes"
                log_root.mkdir(parents=True, exist_ok=True)
                (log_root / f"{case_id}.json").write_text(
                    json.dumps(
                        {
                            "command": case["mapping_command"],
                            "returncode": process.returncode,
                            "stdout": process.stdout,
                            "stderr": process.stderr,
                        },
                        indent=2,
                    )
                    + "\n"
                )
                if process.returncode:
                    raise RuntimeError(f"Mapping command exited {process.returncode}")
            mapping = read_json(mapping_path)
            if "physical_pairing_receipt" in case:
                if (
                    file_sha256(case["physical_pairing_receipt"])
                    != case["physical_pairing_sha256"]
                ):
                    raise ValueError("Physical pairing receipt hash mismatch")
                mapping["physical_pairing_evidence_sha256"] = case[
                    "physical_pairing_sha256"
                ]
            if "protomer_groups" in case:
                mapping.update(
                    protomer_groups=case["protomer_groups"],
                    physical_pairs_verified=case["physical_pairs_verified"],
                )
            result = score_model_raw(model, mapping, runtime, output)
            states[case_id] = {"status": result["status"], "output": str(output)}
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            states[case_id] = {"status": "failed", "error": str(exc)}
            LOGGER.exception("Optional ipSAE2 case failed: %s", case_id)
            failure_path.parent.mkdir(parents=True, exist_ok=True)
            failure_path.write_text(json.dumps(states[case_id], indent=2) + "\n")
        _journal(config, "outcome", f"ipSAE2 case {case_id}: {states[case_id]}")
    return {"status": "polled", "cases": states}
