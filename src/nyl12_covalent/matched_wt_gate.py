"""Strict process-completion evidence for the twelve matched WT replicas."""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any


class ReceiptReader:
    """Read future receipts and bound retries for status files being written."""

    def __init__(self, partial_limit: int = 3) -> None:
        self.partial_limit = partial_limit
        self.partial_reads: dict[str, int] = {}

    def read(self, name: str, *, status: bool = False) -> dict[str, Any] | None:
        """Read an object; missing files and bounded partial status writes wait."""
        path = Path(name)
        if not path.exists():
            return None
        raw = path.read_bytes()
        if status and not raw.endswith(b"\n"):
            self.partial_reads[name] = self.partial_reads.get(name, 0) + 1
            if self.partial_reads[name] <= self.partial_limit:
                return None
            raise ValueError(f"Status remained unterminated: {name}")
        self.partial_reads.pop(name, None)
        record = json.loads(raw)
        if not isinstance(record, dict):
            raise ValueError(f"Receipt must be an object: {name}")
        return record


def digest(path: str) -> str:
    """Return a source or topology SHA256 without loading it all in memory."""
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify_sources(config: dict[str, Any]) -> None:
    """Verify every explicitly frozen source and the analysis manifest."""
    sources = config["source_sha256"]
    if not isinstance(sources, dict) or not sources:
        raise ValueError("A nonempty source_sha256 mapping is required")
    expected = {**sources, config["source_manifest"]: config["source_manifest_sha256"]}
    for name, sha256 in expected.items():
        if digest(name) != sha256:
            raise ValueError(f"Frozen source SHA256 mismatch: {name}")


def _terminal_job(record: dict[str, Any], name: str) -> bool:
    if "ended_utc" not in record:
        if "exit_status" in record and record["exit_status"] is not None:
            raise ValueError(f"Job exit lacks ended_utc: {name}")
        return False
    if record["ended_utc"] is None:
        return False
    ended = datetime.fromisoformat(record["ended_utc"].replace("Z", "+00:00"))
    if ended.utcoffset() is None or ended.utcoffset().total_seconds() != 0:
        raise ValueError(f"Job ended_utc must be UTC: {name}")
    if type(record["exit_status"]) is not int or record["exit_status"] != 0:
        raise RuntimeError(
            f"Upstream replica job failed: {name}: {record['exit_status']}"
        )
    return True


def _stage(record: dict[str, Any], stage: int, name: str) -> None:
    steps, time_ps = (5_000_000, 10500) if stage == 5 else (45_000_000, 100500)
    if type(record["completed_steps"]) is not int or record["completed_steps"] != steps:
        raise ValueError(f"Wrong completed_steps: {name}")
    value = record["final_time_ps"]
    if isinstance(value, bool) or not math.isclose(
        value, time_ps, rel_tol=0, abs_tol=1e-6
    ):
        raise ValueError(f"Wrong final_time_ps: {name}")
    temperature = record["final_temperature_K"]
    if isinstance(temperature, bool) or not 250 <= temperature <= 350:
        raise ValueError(f"Temperature outside gate: {name}")
    if record["temperature_gate_pass"] is not True:
        raise ValueError(f"Temperature gate failed: {name}")
    if stage == 6 and (
        isinstance(record["total_production_ns"], bool)
        or record["total_production_ns"] != 100
    ):
        raise ValueError(f"Wrong total_production_ns: {name}")


def inspect_cases(config: dict[str, Any], reader: ReceiptReader) -> list[str]:
    """Return accepted replica IDs; reject invalid existing terminal evidence."""
    for name in config["worker_statuses"]:
        record = reader.read(name, status=True)
        if record is not None and record["status"] == "failed":
            raise RuntimeError(f"Upstream queue worker failed: {name}")
    accepted = []
    for case in config["cases"]:
        job = reader.read(case["replica_job_status"], status=True)
        terminal = job is not None and _terminal_job(job, case["id"])
        stages = [reader.read(case[f"acceptance_0{stage}"]) for stage in (5, 6)]
        for stage, record in zip((5, 6), stages, strict=True):
            if record is not None:
                _stage(record, stage, case["id"])
        complete = reader.read(case["replica_complete"])
        if complete is not None:
            if "completed" not in complete or type(complete["completed"]) is not bool:
                raise ValueError(f"Invalid completed flag: {case['id']}")
            if terminal and complete["completed"] is not True:
                raise ValueError(f"Terminal replica incomplete: {case['id']}")
        if terminal and all(record is not None for record in stages) and complete:
            if complete["completed"] is True:
                accepted.append(case["id"])
    return accepted


def verify_topologies(config: dict[str, Any]) -> None:
    """Require every production topology to match its frozen canonical hash."""
    for case in config["cases"]:
        if digest(case["production_topology"]) != case["expected_topology_sha256"]:
            raise ValueError(f"Production topology SHA256 mismatch: {case['id']}")
