#!/usr/bin/env python3
"""Append or validate a small, immutable JSONL agent action journal."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import logging
import os
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)
CORE_FIELDS = (
    "timestamp_utc",
    "agent",
    "phase",
    "action",
    "permission_mode",
    "reason",
    "rollback",
)
PHASES = {"intent", "outcome", "note"}
PERMISSION_MODES = {"sandbox", "escalated", "tool"}
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
LEGACY_COMMIT_FILE_VALUES = {
    "a3346f0a-c876-4989-88bf-6ca3193e7ca3": "7845d12",
    "01b88732-2acc-4fda-a955-8f91f4f11af2": "b7fdc7a",
    "6f139086-963d-4d90-9a85-62f7a95c656c": "eccf714",
    "dd1493d9-e027-4c4b-afcc-c9344a5d08f8": "eccf714",
    "b0ec56c1-59a2-4fbc-97c0-59522972751a": "4998ae7",
    "89f5cde5-0000-4139-a2c3-8573cde3ced0": "8b34c06",
}


def file_hashes(paths: list[str]) -> dict[str, str]:
    """Hash explicitly named files, recording missing paths as absent."""
    result = {}
    for name in paths:
        path = Path(name)
        if name in result:
            raise ValueError(f"Duplicate --file path: {name}")
        if not path.exists():
            result[name] = "absent"
        elif not path.is_file():
            raise ValueError(f"--file path is not a regular file: {name}")
        else:
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            result[name] = digest.hexdigest()
    return result


def validate_record(record: Any, line_number: int) -> None:
    """Validate new schema strictly and tolerate sparse historical messages."""
    if not isinstance(record, dict):
        raise ValueError(f"Journal line {line_number} is not a JSON object")
    if record.get("schema_version") != 1:
        for field in ("agent", "action"):
            if field in record and (
                not isinstance(record[field], str) or not record[field].strip()
            ):
                raise ValueError(f"Journal line {line_number} has invalid {field}")
        return
    if set(CORE_FIELDS) - record.keys():
        missing = sorted(set(CORE_FIELDS) - record.keys())
        raise ValueError(f"Journal line {line_number} lacks fields: {missing}")
    required = CORE_FIELDS
    for field in required:
        if not isinstance(record.get(field), str) or not record[field]:
            raise ValueError(f"Journal line {line_number} has invalid {field}")
    for field in set(CORE_FIELDS) - set(required):
        if field in record and not isinstance(record[field], str):
            raise ValueError(f"Journal line {line_number} has invalid {field}")
    if record["phase"] not in PHASES:
        raise ValueError(f"Journal line {line_number} has invalid phase")
    if record.get("permission_mode") not in (None, *PERMISSION_MODES):
        raise ValueError(f"Journal line {line_number} has invalid permission_mode")
    if not isinstance(record.get("event_uuid"), str):
        raise ValueError(f"Journal line {line_number} has invalid event_uuid")
    try:
        uuid.UUID(record["event_uuid"])
    except (ValueError, TypeError, AttributeError) as error:
        raise ValueError(
            f"Journal line {line_number} has invalid event_uuid"
        ) from error
    files = record.get("files", {})
    if not isinstance(files, dict):
        raise ValueError(f"Journal line {line_number} has invalid files map")
    for name, digest in files.items():
        if not isinstance(name, str) or not name:
            raise ValueError(f"Journal line {line_number} has invalid file path")
        if (
            name == "commit"
            and LEGACY_COMMIT_FILE_VALUES.get(record["event_uuid"]) == digest
        ):
            continue
        if digest != "absent" and (
            not isinstance(digest, str) or not SHA256_PATTERN.fullmatch(digest)
        ):
            raise ValueError(f"Journal line {line_number} has invalid file digest")


def verify_journal(log_path: Path) -> int:
    """Read-only validate JSONL records and return their count."""
    count = 0
    with log_path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                raise ValueError(f"Journal contains a blank line at {line_number}")
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Malformed JSON at journal line {line_number}"
                ) from error
            validate_record(record, line_number)
            count += 1
    return count


def append_action(
    log_path: Path,
    agent: str,
    phase: str,
    action: str,
    permission_mode: str,
    reason: str,
    rollback: str,
    result: str | None = None,
    files: list[str] | None = None,
) -> dict[str, Any]:
    """Append one fsynced action record while holding an exclusive file lock."""
    if (
        not agent.strip()
        or not action.strip()
        or not reason.strip()
        or not rollback.strip()
    ):
        raise ValueError("Agent, action, reason, and rollback must be non-empty")
    if phase not in PHASES or permission_mode not in PERMISSION_MODES:
        raise ValueError("Invalid phase or permission mode")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a+", encoding="utf-8") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        stream.seek(0)
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                raise ValueError(f"Journal contains a blank line at {line_number}")
            try:
                validate_record(json.loads(line), line_number)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Malformed JSON at journal line {line_number}"
                ) from error
        record = {
            "schema_version": 1,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "event_uuid": str(uuid.uuid4()),
            "agent": agent,
            "phase": phase,
            "action": action,
            "permission_mode": permission_mode,
            "reason": reason,
            "rollback": rollback,
            "result": result,
            "files": file_hashes(files or []),
        }
        stream.seek(0, 2)
        stream.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
        fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    return record


def main() -> None:
    """Append one action record or validate a journal without modifying it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--agent")
    parser.add_argument("--phase", choices=sorted(PHASES))
    parser.add_argument("--action")
    parser.add_argument("--permission-mode", choices=sorted(PERMISSION_MODES))
    parser.add_argument("--reason")
    parser.add_argument("--rollback")
    parser.add_argument("--result")
    parser.add_argument("--file", action="append", default=[])
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.verify:
        count = verify_journal(args.log)
        LOGGER.info("Validated %d journal records in %s", count, args.log)
        return
    required = ("agent", "phase", "action", "permission_mode", "reason", "rollback")
    missing = [name for name in required if not getattr(args, name)]
    if missing:
        parser.error(
            "append mode requires: "
            + ", ".join("--" + name.replace("_", "-") for name in missing)
        )
    record = append_action(
        args.log,
        args.agent,
        args.phase,
        args.action,
        args.permission_mode,
        args.reason,
        args.rollback,
        args.result,
        args.file,
    )
    LOGGER.info("Appended %s event %s", record["phase"], record["event_uuid"])


if __name__ == "__main__":
    main()
