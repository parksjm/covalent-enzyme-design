"""Guard costly minimum verification against incomplete numerical qualification."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from scripts import run_relay_minimum_hessian as driver
from scripts.run_relay_minimum_hessian import qualify


def _report() -> dict:
    return {
        "dimension": 12495,
        "engineering_relative_column_gate": 1e-4,
        "two_adjacent_fine_pair_passes": True,
        "HB_yz_reciprocity_improved_at_least_tenfold": True,
        "adjacent_sampled_step_comparisons": [
            {"coarse_step_A": 0.000125, "fine_step_A": 0.0000625, "all_nine_pass": True}
        ],
    }


def test_accepts_unchanged_precision_plateau() -> None:
    qualify(_report())


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("dimension", 300),
        ("engineering_relative_column_gate", 1e-3),
        ("two_adjacent_fine_pair_passes", False),
        ("HB_yz_reciprocity_improved_at_least_tenfold", False),
    ],
)
def test_rejects_unqualified_or_relaxed_precision(key: str, value: object) -> None:
    report = _report()
    report[key] = value
    with pytest.raises(ValueError):
        qualify(report)


def test_rejects_selected_step_without_matching_fine_pair() -> None:
    report = deepcopy(_report())
    report["adjacent_sampled_step_comparisons"][0]["all_nine_pass"] = False
    with pytest.raises(ValueError):
        qualify(report)


def _prepare(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    stage = tmp_path / "05_fine_step_diagnostic"
    stage.mkdir()
    report = _report()
    report["source_column_sha256"] = {}
    (stage / "diagnostic_acceptance.json").write_text(json.dumps(report))
    monkeypatch.setattr("sys.argv", ["driver", "--root", str(tmp_path)])
    monkeypatch.setattr(driver, "verify_frozen", lambda path: None)
    monkeypatch.setattr(driver, "resource_gate", lambda *args: {"passed": True})
    monkeypatch.setattr(driver, "journal", lambda *args: None)
    return stage / "diagnostic_acceptance.json"


@pytest.mark.parametrize("passed", [True, False])
def test_full_result_controls_minimum_claim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, passed: bool
) -> None:
    _prepare(tmp_path, monkeypatch)
    calls = []

    def hessian(*args: object, **kwargs: object) -> dict:
        calls.append((args, kwargs))
        return {"passed": passed}

    monkeypatch.setattr(driver, "full_hessian", hessian)
    driver.main()
    summary = json.loads((tmp_path / "verification_summary.json").read_text())
    assert summary["numerically_verified"] is passed
    assert calls[0][1]["expected_negative"] == 0
    assert calls[0][0][4] == 0.000125


def test_preflight_only_does_not_launch_matrix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _prepare(tmp_path, monkeypatch)
    monkeypatch.setattr(
        "sys.argv", ["driver", "--root", str(tmp_path), "--preflight-only"]
    )
    driver.main()
    assert not (tmp_path / "verification_summary.json").exists()


def test_matrix_exception_retains_incomplete_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _prepare(tmp_path, monkeypatch)

    def fail(*args: object, **kwargs: object) -> dict:
        raise RuntimeError("Native gradient failed")

    monkeypatch.setattr(driver, "full_hessian", fail)
    with pytest.raises(RuntimeError, match="gradient failed"):
        driver.main()
    summary = json.loads((tmp_path / "verification_summary.json").read_text())
    assert summary["status"] == "incomplete"
    assert not summary["numerically_verified"]


def test_changed_column_evidence_blocks_matrix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _prepare(tmp_path, monkeypatch)
    column = tmp_path / "column.npy"
    column.write_bytes(b"altered")
    report = json.loads(source.read_text())
    report["source_column_sha256"] = {"column.npy": "incorrect"}
    source.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="evidence changed"):
        driver.main()


def test_no_resource_headroom_blocks_matrix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _prepare(tmp_path, monkeypatch)
    monkeypatch.setattr(driver, "resource_gate", lambda *args: {"passed": False})
    with pytest.raises(RuntimeError, match="headroom"):
        driver.main()


def test_rejects_excess_workers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "sys.argv", ["driver", "--root", str(tmp_path), "--workers", "9"]
    )
    with pytest.raises(ValueError, match="eight"):
        driver.main()
