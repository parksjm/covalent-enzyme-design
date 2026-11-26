"""Exercise the actual flat native receipt schema and rejection gates."""

import hashlib
import json
from pathlib import Path

import pytest

from nyl12_covalent.hydration_receipt_adapter import accepted_extensions


def fixture_block(base: Path, index: int = 1) -> Path:
    """Create a small block matching the production extension writer."""
    segment = base / f"extend_{index:02d}"
    segment.mkdir()
    (segment / "restart.ncrst").write_bytes(b"accepted restart")
    (segment / "trajectory.nc").write_bytes(b"accepted trajectory")
    record = {
        "numeric": {
            "completed_steps": 10000,
            "final_temperature_K": 298.0,
            "electronic_history": {"verified_scc_history": True},
        },
        "boundary": {"box_unchanged": True, "frozen_max_displacement_A": 0.0},
        "C_OG_A": 1.373,
        "scc_verified": True,
        "acyl_retained": True,
        "finite_temperature_passed": True,
        "native_process_passed": True,
    }
    (segment / "scientific_acceptance.json").write_text(json.dumps(record))
    (segment / "process.json").write_text(
        json.dumps({"status": "process_completed", "exit_code": 0})
    )
    (segment / "status.json").write_text(
        json.dumps(
            {
                "status": "accepted",
                "restart_sha256": hashlib.sha256(b"accepted restart").hexdigest(),
            }
        )
    )
    return segment


def test_flat_schema_and_pending_block(tmp_path: Path) -> None:
    """Accept real schema and exclude a running successor."""
    first = fixture_block(tmp_path)
    second = tmp_path / "extend_02"
    second.mkdir()
    (second / "status.json").write_text('{"status":"running"}')
    assert accepted_extensions(tmp_path) == [first]


@pytest.mark.parametrize(
    "field,value",
    [
        ("C_OG_A", 1.9),
        ("C_OG_A", float("nan")),
        ("acyl_retained", False),
        ("scc_verified", False),
        ("finite_temperature_passed", False),
        ("native_process_passed", False),
    ],
)
def test_bad_native_gates(tmp_path: Path, field: str, value: object) -> None:
    """Never relax native acceptance to fix a processing mismatch."""
    segment = fixture_block(tmp_path)
    path = segment / "scientific_acceptance.json"
    record = json.loads(path.read_text())
    record[field] = value
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="acceptance failed"):
        accepted_extensions(tmp_path)


@pytest.mark.parametrize("failure", ["restart", "trajectory", "gap", "process"])
def test_identity_and_completeness(tmp_path: Path, failure: str) -> None:
    """Reject changed native artifacts and orphaned accepted blocks."""
    segment = fixture_block(tmp_path, 2 if failure == "gap" else 1)
    if failure == "restart":
        (segment / "restart.ncrst").write_bytes(b"changed")
    elif failure == "trajectory":
        (segment / "trajectory.nc").write_bytes(b"")
    elif failure == "process":
        (segment / "process.json").write_text(
            '{"status":"process_failed","exit_code":1}'
        )
    with pytest.raises(ValueError):
        accepted_extensions(tmp_path)
