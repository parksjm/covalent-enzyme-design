"""Bounded support-grid preparation, immutable hashes, and cross-bias contracts."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from test_relay_targeted import plan as plan
from nyl12_covalent.relay_activation_support import (
    prepare_support_targets,
    validate_support_grid,
)
from nyl12_covalent.relay_targeted import verify_target
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies


def support_plan(plan: Path) -> tuple[Path, dict]:
    """Adapt the existing validated fixture to the exact approved33nodes."""
    data = json.loads(plan.read_text())
    base = data["targets"][0]
    data["output"] = str(plan.parent / "activation_support_20261127")
    data["targets"] = [
        {
            **base,
            "label": f"support_{row}_{col:02d}",
            "feature_names": ["q_acyl", "q_act"],
            "centers_A": [x, round(-1 + 0.2 * col, 8)],
            "initial_ig": 10000 + 10 * (11 * row + col),
        }
        for row, x in enumerate((-1.25, -1.5, -1.75))
        for col in range(11)
    ]
    plan.write_text(json.dumps(data))
    return plan, data


def test_exact_grid_and_immutable_preparation(plan: Path) -> None:
    """All33nativeinputsreuseauditedmasksandcanbeverifiedbytheoldruntime."""
    path, data = support_plan(plan)
    prepared = prepare_support_targets(path)
    digest = prepared.with_suffix(".sha256").read_text().strip()
    cfg = verify_target(prepared, digest)
    assert len(cfg["windows"]) == 33
    assert cfg["physical_hamiltonian_sha256"] == "parent Hamiltonian"
    assert cfg["system"]["active_mask"] == ":1-20,25"
    assert cfg["production_horizon_ps"] == 30
    stages = cfg["windows"][0]["stages"]
    first = (prepared.parent / stages[0]["input"]).read_text()
    final = (prepared.parent / stages[-1]["input"]).read_text()
    assert "irest=0, ntx=1, tempi=300.0" in first
    assert "irest=1, ntx=5" in final
    assert "qmcharge=-2" in final and "DFTB3" in final
    assert [s["steps"] for s in stages] == [2000, 8000, *[10000] * 6]
    with pytest.raises(FileExistsError):
        prepare_support_targets(path)
    seed = prepared.parent / cfg["windows"][0]["seed"]
    seed.write_text("changed")
    with pytest.raises(ValueError, match="Pinned file"):
        verify_target(prepared, digest)


@pytest.mark.parametrize("change", ["missing", "third_cv", "force"])
def test_reject_unapproved_grid(plan: Path, change: str) -> None:
    """Rejectmissingstates,addedq_delbias,andchangedforceconstants."""
    _, data = support_plan(plan)
    targets = data["targets"]
    if change == "missing":
        targets.pop()
    elif change == "third_cv":
        targets[0]["feature_names"] = ["q_acyl", "q_act", "q_del"]
    else:
        targets[0]["k"] = 50
    with pytest.raises(ValueError):
        validate_support_grid(targets)


def test_native_restraint_and_joint_cross_bias_contract(plan: Path) -> None:
    """The33newstatesandoldtwobiasstatecanbeevaluatedonthefullsameframepool."""
    path, _ = support_plan(plan)
    cfg = json.loads(prepare_support_targets(path).read_text())
    biases = [RelayBiasWindow(**row["bias"]) for row in cfg["windows"]]
    biases.append(RelayBiasWindow(("q_acyl", "q_del"), (-0.5, -0.7), (25.0, 25.0)))
    raw = {"q_acyl": [-1.5, -1.7], "q_act": [-0.4, 0.1], "q_del": [-2.2, -1.8]}
    table = cross_bias_energies(raw, biases)
    assert table.shape == (34, 2) and np.isfinite(table).all()
    for row, bias in enumerate(biases):
        expected = [
            sum(
                k * (raw[name][i] - center) ** 2
                for name, center, k in zip(
                    bias.feature_names,
                    bias.centers_A,
                    bias.force_constants,
                    strict=True,
                )
            )
            for i in range(2)
        ]
        np.testing.assert_allclose(table[row], expected)
    rst = biases[0].restraint_text()
    assert rst.count("&rst") == 2 and "rk2=25" in rst
    assert "13916,13917,13905,13917" in rst


def test_runtime_reaches30ps_after10ps_diagnostic_without_restart(
    plan: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Allsixproductionblocksruninonecallwithuncensoreddiagnostics."""
    import test_relay_targeted as existing_tests
    from nyl12_covalent import relay_activation_support_runtime as runtime
    from nyl12_covalent.neb_umbrella_campaign import digest

    path, _ = support_plan(plan)
    prepared = prepare_support_targets(path)
    monkeypatch.setattr(existing_tests, "runtime", runtime)
    checkpoints, diagnostics = existing_tests._fake_runtime(monkeypatch)
    label = "support_0_00"
    runtime.run_support_target(prepared, label, digest(prepared))
    folder = prepared.parent / label
    status = json.loads((folder / "window_status.json").read_text())
    assert status["production_ps"] == 30 and status["status"] == "completed"
    assert checkpoints.call_count == diagnostics.call_count == 8
    checkpoint = json.loads((folder / "diagnostic_checkpoint_10ps.json").read_text())
    assert checkpoint["production_ps"] == 10
    assert checkpoint["automatic_continuation_to_ps"] == 30
    assert checkpoint["checkpoint"].endswith("04_production_5ps/restart.ncrst")
    with pytest.raises(FileExistsError, match="already attempted"):
        runtime.run_support_target(prepared, label, digest(prepared))
