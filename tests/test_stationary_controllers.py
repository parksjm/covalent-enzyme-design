"""Controller gating tests with analytic derivatives and explicit stage receipts."""

import gzip
import json
import sys

import numpy as np
import pytest

from scripts import run_stationary_accuracy_diagnostic as accuracy
from scripts import run_stationary_validation as loop


@pytest.mark.parametrize("stable", [True, False])
def test_accuracy_controller_gates_full_loop(tmp_path, monkeypatch, stable):
    monkeypatch.setattr(sys, "argv", ["accuracy", "--root", str(tmp_path)])
    monkeypatch.setattr(accuracy, "NativeState", lambda root: object())
    calls = []
    monkeypatch.setattr(
        accuracy.subprocess, "run", lambda *args, **kwargs: calls.append(args)
    )

    def derivative(state, column, step, output):
        output.mkdir(parents=True)
        with gzip.open(output / "native.mdout.gz", "wt") as stream:
            stream.write("QMMM SCC-DFTB: 42 -1.0 0.0 0.0 1 2 H 0.0\n")
        return np.ones(3) + (step**2 if stable else step * 100)

    monkeypatch.setattr(accuracy, "column_difference", derivative)
    if stable:
        accuracy.main()
        assert len(calls) == 1
    else:
        with pytest.raises(RuntimeError, match="stability"):
            accuracy.main()
        assert not calls
    receipt = json.loads(
        (tmp_path / "accuracy_diagnostic/scientific_acceptance.json").read_text()
    )
    assert receipt["precision_screen_passed"] is stable
    assert receipt["scfconv"] == 1e-12
    assert receipt["electronic_records"][0]["last_verbose_scc_row"][0] == "42"


@pytest.mark.parametrize("precision_passed", [True, False])
def test_validation_controller_waits_for_precision(
    tmp_path, monkeypatch, precision_passed
):
    monkeypatch.setattr(
        sys, "argv", ["loop", "--root", str(tmp_path), "--precision-only"]
    )
    monkeypatch.setattr(
        loop, "precision_gate", lambda *args: {"passed": precision_passed}
    )
    if precision_passed:
        loop.main()
    else:
        with pytest.raises(RuntimeError, match="Precision"):
            loop.main()


def test_validation_controller_requires_index_one_before_downhill(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(sys, "argv", ["loop", "--root", str(tmp_path)])
    (tmp_path / "precision_acceptance.json").write_text(json.dumps({"passed": True}))
    monkeypatch.setattr(loop, "full_hessian", lambda *args: {"passed": False})
    with pytest.raises(RuntimeError, match="index one"):
        loop.main()


def test_validation_controller_completes_only_with_both_minimum_hessians(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(sys, "argv", ["loop", "--root", str(tmp_path)])
    (tmp_path / "precision_acceptance.json").write_text(json.dumps({"passed": True}))
    (tmp_path / "precision/precise").mkdir(parents=True)
    (tmp_path / "precision/precise/receipt.json").write_text(
        json.dumps({"energy_hartree": -1.0})
    )
    restarts = [tmp_path / "minus/restart", tmp_path / "plus/restart"]
    monkeypatch.setattr(loop, "endpoints", lambda *args: restarts)
    calls = []

    def full(*args):
        calls.append(args)
        return {"passed": True}

    monkeypatch.setattr(loop, "full_hessian", full)
    loop.main()
    assert len(calls) == 3
    assert [args[-1] for args in calls] == [1, 0, 0]
    receipt = json.loads((tmp_path / "scientific_acceptance.json").read_text())
    assert receipt["passed"] is True
