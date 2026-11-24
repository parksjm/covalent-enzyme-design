"""Scientific engineering gates for the bounded fine-step diagnostic."""

from __future__ import annotations

import numpy as np
import pytest

from scripts.diagnose_relay_minimum_steps import COLUMNS, FINE_STEPS, analyze


def blocks() -> list[np.ndarray]:
    result = []
    for i in range(6):
        block = np.ones((12495, len(COLUMNS)))
        block[6587, COLUMNS.index(6586)] += 1e-3 if i == 0 else 1e-7
        result.append(block)
    return result


def test_plateau_and_reciprocity_are_separate_from_minimum() -> None:
    result = analyze(blocks(), [0.001, 0.0005, *FINE_STEPS])
    assert result["two_adjacent_fine_pair_passes"]
    assert result["HB_yz_reciprocity_improved_at_least_tenfold"]
    assert result["status"] == "ready_for_lead_review"
    assert not result["numerically_verified_minimum"]
    assert result["engineering_relative_column_gate"] == 1e-4


def test_unstable_columns_prevent_plateau_even_if_hb_reciprocity_improves() -> None:
    data = blocks()
    for i, block in enumerate(data):
        block[0, 0] += i * 0.1
    result = analyze(data, [0.001, 0.0005, *FINE_STEPS])
    assert not result["two_adjacent_fine_pair_passes"]
    assert result["status"] == "not_stable_or_reciprocity_unresolved"


def test_stable_asymmetry_stays_unresolved() -> None:
    data = blocks()
    for block in data:
        block[6587, COLUMNS.index(6586)] = 1.01
    result = analyze(data, [0.001, 0.0005, *FINE_STEPS])
    assert result["two_adjacent_fine_pair_passes"]
    assert not result["HB_yz_reciprocity_improved_at_least_tenfold"]
    assert result["reciprocity_by_step"][-1][
        "HB_yz_absolute_difference"
    ] == pytest.approx(0.01)


@pytest.mark.parametrize(
    "data,steps",
    [
        ([], []),
        ([np.ones((2, 2))] * 3, [0.1, 0.05, 0.025]),
        ([np.ones((12495, 9))] * 3, [0.1, 0.1, 0.025]),
        ([np.ones((12495, 9))] * 3, [0.1, 0.05, np.nan]),
    ],
)
def test_invalid_or_changed_column_specification_rejected(
    data: list[np.ndarray], steps: list[float]
) -> None:
    with pytest.raises(ValueError):
        analyze(data, steps)


def test_complete_bounded_run_preserves_original_columns(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    import json
    from concurrent.futures import ThreadPoolExecutor
    from scripts import diagnose_relay_minimum_steps as job

    stage = tmp_path / "05_fine_step_diagnostic"
    stage.mkdir()
    frozen = tmp_path / "01_frozen"
    for parent in ("inputs", "library_build"):
        (frozen / parent).mkdir(parents=True)
    (frozen / "inputs/optimized.ncrst").write_bytes(b"accepted")
    (frozen / "library_build/libdlfind.so").write_bytes(b"unchanged")
    original = tmp_path / "03_step_probes"
    original.mkdir()
    (original / "scientific_acceptance.json").write_text(
        json.dumps({"columns": {str(c): {} for c in COLUMNS}})
    )
    for label in ("h_0p001A", "h_0p0005A"):
        for c in COLUMNS:
            directory = original / label / "columns" / f"column_{c:05d}"
            directory.mkdir(parents=True)
            np.save(directory / "column.npy", np.ones(12495))
    expected = (original / "h_0p001A/columns/column_00000/column.npy").read_bytes()
    monkeypatch.setattr(job, "verify_frozen", lambda *args: None)
    monkeypatch.setattr(job, "verify_prerequisites", lambda *args: None)
    monkeypatch.setattr(job, "resource_gate", lambda *args: {"passed": True})
    monkeypatch.setattr(job, "journal", lambda *args: None)
    monkeypatch.setattr(job, "ProcessPoolExecutor", ThreadPoolExecutor)
    calls = []

    def worker(task: tuple) -> str:
        root, step_stage, restart, c, h = task
        directory = step_stage / "columns" / f"column_{c:05d}"
        directory.mkdir(parents=True)
        np.save(directory / "column.npy", np.ones(12495))
        calls.append((c, h))
        return str(directory / "column.npy")

    monkeypatch.setattr(job, "worker_column", worker)
    result = job.run(tmp_path, 8)
    assert len(calls) == 36
    assert set(h for c, h in calls) == set(FINE_STEPS)
    assert not result["numerically_verified_minimum"]
    assert len(result["source_column_sha256"]) == 54
    assert (
        expected == (original / "h_0p001A/columns/column_00000/column.npy").read_bytes()
    )
    with pytest.raises(ValueError):
        job.run(tmp_path, 9)
