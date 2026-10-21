"""Concurrency, schema, immutability, and hash tests for the action journal."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from scripts.record_agent_action import append_action, file_hashes, verify_journal
from scripts import record_agent_action as journal_cli

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/record_agent_action.py"


def _append_cli(log_path: Path, sequence: int) -> None:
    """Append a record in a child process for the lock regression test."""
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--log",
            str(log_path),
            "--agent",
            "concurrency-test",
            "--phase",
            "note",
            "--action",
            f"append {sequence}",
            "--permission-mode",
            "tool",
            "--reason",
            "concurrent append test",
            "--rollback",
            "no repository change",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def test_file_hashes_record_sha256_and_absent_paths(tmp_path: Path) -> None:
    present = tmp_path / "input.dat"
    present.write_bytes(b"known input")
    absent = tmp_path / "missing.dat"

    result = file_hashes([str(present), str(absent)])

    assert result[str(present)] == hashlib.sha256(b"known input").hexdigest()
    assert result[str(absent)] == "absent"


def test_duplicate_and_directory_file_arguments_fail(tmp_path: Path) -> None:
    directory = tmp_path / "folder"
    directory.mkdir()
    with pytest.raises(ValueError, match="Duplicate"):
        file_hashes(["same", "same"])
    with pytest.raises(ValueError, match="regular file"):
        file_hashes([str(directory)])


def test_append_and_verify_preserve_sparse_legacy_records(tmp_path: Path) -> None:
    log_path = tmp_path / "journal.jsonl"
    log_path.write_text(
        json.dumps({"agent": "bootstrap", "action": "older sparse event"}) + "\n",
        encoding="utf-8",
    )
    source = tmp_path / "source.py"
    source.write_text("content", encoding="utf-8")

    record = append_action(
        log_path,
        "test-agent",
        "outcome",
        "write record",
        "sandbox",
        "test action",
        "no mutation to undo",
        result="done",
        files=[str(source)],
    )

    assert record["schema_version"] == 1
    assert len(record["event_uuid"]) == 36
    assert record["files"][str(source)] == file_hashes([str(source)])[str(source)]
    assert verify_journal(log_path) == 2


def test_malformed_existing_log_fails_without_appending(tmp_path: Path) -> None:
    log_path = tmp_path / "journal.jsonl"
    log_path.write_text("{truncated\n", encoding="utf-8")
    before = log_path.read_bytes()

    with pytest.raises(ValueError, match="Malformed JSON"):
        append_action(
            log_path,
            "test-agent",
            "note",
            "must not append",
            "sandbox",
            "test malformed journal handling",
            "nothing to undo",
        )

    assert log_path.read_bytes() == before


def test_concurrent_cli_appends_preserve_every_record(tmp_path: Path) -> None:
    log_path = tmp_path / "journal.jsonl"
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(lambda number: _append_cli(log_path, number), range(24)))

    assert verify_journal(log_path) == 24
    records = [json.loads(line) for line in log_path.read_text().splitlines()]
    assert {record["action"] for record in records} == {
        f"append {number}" for number in range(24)
    }
    assert len({record["event_uuid"] for record in records}) == 24


def test_verify_rejects_invalid_digest_without_modifying_log(tmp_path: Path) -> None:
    log_path = tmp_path / "journal.jsonl"
    record = {
        "schema_version": 1,
        "timestamp_utc": "2026-10-21T00:00:00+00:00",
        "event_uuid": "55b28490-410a-419f-a88d-9b928558d353",
        "agent": "test-agent",
        "phase": "note",
        "action": "bad digest",
        "permission_mode": "tool",
        "reason": "negative test",
        "rollback": "none",
        "files": {"input": "not-a-sha"},
    }
    log_path.write_text(json.dumps(record) + "\n", encoding="utf-8")
    before = log_path.read_bytes()

    with pytest.raises(ValueError, match="digest"):
        verify_journal(log_path)

    assert log_path.read_bytes() == before


def test_cli_append_verify_and_missing_options(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise the real CLI without losing child-process coverage."""
    log = tmp_path / "actions.jsonl"
    argv = ["journal", "--log", str(log)]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit) as failure:
        journal_cli.main()
    assert failure.value.code == 2
    assert not log.exists()
    monkeypatch.setattr(
        sys,
        "argv",
        argv
        + [
            "--agent",
            "cli-test",
            "--phase",
            "intent",
            "--action",
            "read a file",
            "--permission-mode",
            "sandbox",
            "--reason",
            "inspect data",
            "--rollback",
            "read-only",
            "--file",
            str(tmp_path / "absent"),
        ],
    )
    journal_cli.main()
    before = log.read_bytes()
    assert json.loads(before)["files"] == {str(tmp_path / "absent"): "absent"}
    monkeypatch.setattr(sys, "argv", argv + ["--verify"])
    journal_cli.main()
    assert log.read_bytes() == before


@pytest.mark.parametrize(
    "field,value",
    [
        ("event_uuid", "invalid"),
        ("phase", "unknown"),
        ("permission_mode", "unknown"),
        ("agent", ""),
    ],
)
def test_new_schema_rejects_invalid_fields(
    tmp_path: Path, field: str, value: str
) -> None:
    """Legacy tolerance must not weaken validation of new action records."""
    log = tmp_path / "actions.jsonl"
    record = append_action(log, "test", "note", "check", "tool", "test", "none")
    record[field] = value
    with pytest.raises(ValueError):
        journal_cli.validate_record(record, 1)
