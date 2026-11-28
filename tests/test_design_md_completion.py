"""Check replica routing, process gates and diagnostic-only selection."""

import json
from pathlib import Path

import pytest

from nyl12_covalent.design_md_completion import (
    add_selection_eligibility,
    replica_path,
    sha,
    validate_completion,
)


def make_replica(root: Path, dataset: str, number: str) -> tuple[Path, str, str]:
    """Create a miniature but structurally faithful completion record."""
    replica = replica_path(root, dataset, number)
    dynamics = replica / "dynamics"
    dynamics.mkdir(parents=True)
    (replica / "100ns_process_complete.json").write_text('{"completed": true}')
    for name in (
        "system.prmtop",
        "04_equil_npt.rst7",
        "05_production_10ns.in",
        "05_production_10ns.rst7",
        "06_production_90ns.in",
        "06_production_90ns.rst7",
    ):
        (dynamics / name).write_text(name)
    top_hash = sha(dynamics / "system.prmtop")
    start_hash = "minimized-source-hash"

    def inputs(names: tuple[str, ...]) -> dict[str, str]:
        return {str(dynamics / name): sha(dynamics / name) for name in names}

    first = {
        "returncode": 0,
        "reference_sha256": start_hash,
        "input_sha256": inputs(
            ("system.prmtop", "04_equil_npt.rst7", "05_production_10ns.in")
        ),
    }
    continuation = {
        "returncode": 0,
        "total_production_ns": 100,
        "input_sha256": inputs(
            ("system.prmtop", "05_production_10ns.rst7", "06_production_90ns.in")
        ),
    }
    (dynamics / "05_production_10ns.process.json").write_text(json.dumps(first))
    (dynamics / "extension_100ns.json").write_text(json.dumps(continuation))
    for stem, steps, time in (
        ("05_production_10ns", 5000000, 10500.0),
        ("06_production_90ns", 45000000, 100500.0),
    ):
        gate = {
            "completed_steps": steps,
            "final_time_ps": time,
            "temperature_gate_pass": True,
            "restart_sha256": sha(dynamics / f"{stem}.rst7"),
        }
        (dynamics / f"{stem}.acceptance.json").write_text(json.dumps(gate))
    return replica, top_hash, start_hash


def test_replica_routing_and_native_diagnostic(tmp_path: Path) -> None:
    assert replica_path(tmp_path, "V111C", "02") == tmp_path / "md/V111C/replica_02"
    assert replica_path(tmp_path, "WT", "01") == tmp_path / "md/WT/replica_01"
    with pytest.raises(ValueError):
        replica_path(tmp_path, "V111C", "03")
    row = add_selection_eligibility(
        {"tier": "B", "native_pass": False, "joint_pass": False}
    )
    assert row["selection_eligible"] is True and row["native_pass"] is False
    assert add_selection_eligibility({"tier": "reject"})["selection_eligible"] is False


def test_accepted_and_failed_receipts(tmp_path: Path) -> None:
    replica, top, start = make_replica(tmp_path, "V111C", "02")
    assert len(validate_completion(replica, top, start)) == 2
    receipt = replica / "100ns_process_complete.json"
    receipt.write_text('{"completed": false}')
    with pytest.raises(ValueError, match="completion"):
        validate_completion(replica, top, start)
    receipt.write_text('{"completed": true}')
    gate = replica / "dynamics/06_production_90ns.acceptance.json"
    original = json.loads(gate.read_text())
    bad = {**original, "temperature_gate_pass": False}
    gate.write_text(json.dumps(bad))
    with pytest.raises(ValueError, match="scientific"):
        validate_completion(replica, top, start)
    gate.write_text(json.dumps(original))
    (replica / "dynamics/05_production_10ns.rst7").write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        validate_completion(replica, top, start)
