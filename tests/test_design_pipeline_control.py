"""Verify bounded control orchestration without launching scientific engines."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from nyl12_covalent import design_pipeline_control as control


def _write(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def _fixture(root: Path) -> dict[str, Any]:
    pinned = root / "pinned"
    pinned.write_text("original")
    plan = {
        "candidates": [
            dict(id=name, alpha_sequence="A", beta_sequence="T")
            for name in ("WT", "variant_01")
        ],
        "candidate_count": 2,
        "physical_gpu": 0,
        "new_sequence_samples": 0,
        "input_sha256": {str(pinned): hashlib.sha256(b"original").hexdigest()},
        "matched_seed": 9,
        "scientific_scope": "technical control",
    }
    _write(root / "plan.json", plan)
    _write(root / "prediction_job/status.json", {"exit_status": 0})
    _write(root / "predictions/WORKER_COMPLETE.json", {"jobs": 2, "complete": 2})
    _write(root / "relaxation.json", {"mol2": "reference.mol2"})
    _write(
        root / "execution.json",
        {
            "relaxation_config": str(root / "relaxation.json"),
            "openmm_python": "/explicit/python",
            "dependency_sha256": {},
            "environment": {"CUDA_VISIBLE_DEVICES": "0"},
            "reference_cif": str(root / "reference.cif"),
        },
    )
    return plan


def test_frozen_inputs_and_budget_are_enforced(tmp_path: Path) -> None:
    """Reject expanded candidate/GPU budgets and changed scientific inputs."""
    plan = _fixture(tmp_path)
    control.verify_inputs(plan)
    plan["candidate_count"] = 3
    with pytest.raises(ValueError, match="exactly"):
        control.verify_inputs(plan)
    plan["candidate_count"] = 2
    plan["physical_gpu"] = 1
    with pytest.raises(ValueError, match="GPU0"):
        control.verify_inputs(plan)
    plan["physical_gpu"] = 0
    (tmp_path / "pinned").write_text("changed")
    with pytest.raises(ValueError, match="hash changed"):
        control.verify_inputs(plan)


def test_prediction_exit_completeness_and_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject failed/incomplete predictions and stop after a bounded wait."""
    _fixture(tmp_path)
    _write(tmp_path / "prediction_job/status.json", {"exit_status": 1})
    with pytest.raises(RuntimeError, match="exited 1"):
        control.wait_for_predictions(tmp_path)
    _write(tmp_path / "prediction_job/status.json", {"exit_status": 0})
    _write(tmp_path / "predictions/WORKER_COMPLETE.json", {"jobs": 2, "complete": 1})
    with pytest.raises(ValueError, match="completeness"):
        control.wait_for_predictions(tmp_path)
    _write(tmp_path / "prediction_job/status.json", {"exit_status": None})
    times = iter([0, 1, 3])
    monkeypatch.setattr(control.time, "monotonic", lambda: next(times))
    monkeypatch.setattr(control.time, "sleep", lambda _: None)
    with pytest.raises(TimeoutError):
        control.wait_for_predictions(tmp_path, timeout_seconds=2)


def test_control_connects_two_predictions_to_openmm_and_geometry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Preserve matched candidates, explicit GPU and separate convergence status."""
    _fixture(tmp_path)
    calls = []

    def normalize(_source: Path, _component: Path, _mol2: Path, output: Path) -> None:
        output.mkdir()
        (output / "normalized.pdb").write_text("canonical")

    def prepare(_raw: Path, output: Path, _config: dict, **kwargs: Any) -> dict:
        assert kwargs == {"preparation_only": True}
        output.mkdir()
        return dict(
            topology="topology", coordinates="coordinates", conversion="mapping"
        )

    def run(argv: list[str], **kwargs: Any) -> None:
        assert kwargs["env"]["CUDA_VISIBLE_DEVICES"] == "0"
        calls.append(argv)
        output = Path(argv[-1])
        output.mkdir()
        _write(output / "acceptance.json", {"force_converged": False})
        np.save(output / "positions_nm.npy", np.zeros((2, 3)))

    class Endpoint:
        coordinates: Any

        def strip(self, mask: str) -> None:
            assert mask == ":WAT,Na+,Cl-"

        def save(self, path: str, **_kwargs: Any) -> None:
            Path(path).write_text("endpoint")

    monkeypatch.setattr(control, "normalize_model", normalize)
    monkeypatch.setattr(
        control,
        "_pdb_sequences",
        lambda _: {**dict.fromkeys("ABCD", "A"), **dict.fromkeys("EFGH", "T")},
    )
    monkeypatch.setattr(control, "relax_model", prepare)
    monkeypatch.setattr(control.subprocess, "run", run)
    monkeypatch.setattr(control.parmed, "load_file", lambda *_a, **_k: Endpoint())
    monkeypatch.setattr(
        control, "export_canonical_heavy", lambda *args: args[-1].write_text("pdb")
    )
    observed = []
    monkeypatch.setattr(
        control,
        "compare_models",
        lambda entries, *_: observed.extend(entries) or {"models": 4},
    )
    result = control.run_control(tmp_path, tmp_path / "execution.json")
    assert len(calls) == 2
    assert {row["variant"] for row in observed} == {"WT", "variant_01"}
    assert {row["seed"] for row in observed} == {9}
    assert result["process_completed"] is True
    assert result["claim_of_improved_activity"] is False
    assert result["cases"][0]["openmm_acceptance"]["force_converged"] is False
    assert (tmp_path / "COMPLETED.json").is_file()


def test_control_rejects_changed_predicted_sequence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A prediction with changed protein identity never reaches topology prep."""
    _fixture(tmp_path)
    monkeypatch.setattr(control, "normalize_model", lambda *_: None)
    monkeypatch.setattr(
        control, "_pdb_sequences", lambda _: dict.fromkeys("ABCDEFGH", "X")
    )
    with pytest.raises(ValueError, match="sequence mismatch"):
        control.run_control(tmp_path, tmp_path / "execution.json")
