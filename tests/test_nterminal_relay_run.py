"""Branch pipeline tests using mock native execution and publisher."""

import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from nyl12_covalent import nterminal_relay_run as runner


def setup_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Mock, Mock, Mock]:
    """Create a campaign and deterministic helper substitutes."""
    campaign = tmp_path / "campaign"
    model = campaign / "model/inputs"
    model.mkdir(parents=True)
    (model / "optimized.ncrst").write_text("start")
    (campaign / "config.json").write_text(
        json.dumps(
            {
                "model": "model",
                "branches": {
                    "reverse_activation": "reverse",
                    "forward_acylation": "forward",
                },
                "reference_energy_hartree": -1,
                "path_final_A": 1.2,
                "path_maxincrement_A": 0.1,
            }
        )
    )
    monkeypatch.setattr(
        runner,
        "NativeState",
        lambda path: SimpleNamespace(xyz=np.zeros((2, 3)), mdin_source="source"),
    )
    monkeypatch.setattr(runner, "geometry", lambda xyz: {"state": "initial"})
    monkeypatch.setattr(runner, "coordinate", lambda xyz, branch, proton: 1.0)
    monkeypatch.setattr(runner, "_restart_xyz", lambda path: np.zeros((2, 3)))
    monkeypatch.setattr(runner, "_driven_mdin", lambda source, cycles: "driven")
    monkeypatch.setattr(
        runner, "minimization_input", lambda source, maxcycles: "unbiased"
    )
    monkeypatch.setattr(
        runner, "amber_restraint", lambda target, branch, proton: "restraint"
    )
    pub = Mock()
    monkeypatch.setattr(runner, "publish", pub)
    ts = Mock(return_value={"status": "ok"})

    def native(
        campaign: Path,
        folder: Path,
        start: Path,
        mdin: str,
        restraints: str | None = None,
    ) -> dict:
        folder.mkdir(parents=True, exist_ok=True)
        restart = folder / "optimized.ncrst"
        restart.write_text("restart")
        return {
            "restart": str(restart),
            "accepted": True,
            "native_converged": True,
            "geometry": {"state": "candidate"},
            "total_kcal_mol": 1.0,
            "bias_kcal_mol": 0.0,
            "physical_kcal_mol": 1.0,
        }

    native_mock = Mock(side_effect=native)
    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.nterminal_relay_native",
        SimpleNamespace(run_native=native_mock),
    )
    monkeypatch.setitem(
        sys.modules, "nyl12_covalent.nterminal_relay_ts", SimpleNamespace(run_ts=ts)
    )
    monkeypatch.setattr(
        runner, "candidate", lambda geometry, branch: geometry["state"] == "candidate"
    )
    return campaign, native_mock, pub, ts


def test_success_and_release(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign, native, pub, ts = setup_branch(tmp_path, monkeypatch)
    result = runner.run_branch(campaign, "reverse_activation")
    assert result["status"] == "completed"
    assert native.call_count == 2 and pub.call_count == 2 and ts.call_count == 1
    assert native.call_args_list[0].kwargs["restraints"] == "restraint"
    assert native.call_args_list[1].kwargs["restraints"] is None
    assert (
        json.loads((campaign / "reverse/branch_status.json").read_text())["status"]
        == "completed"
    )


def test_rejection_and_exhaustion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign, native, pub, ts = setup_branch(tmp_path, monkeypatch)
    original = native.side_effect

    def reject(*args: object, **kwargs: object) -> dict:
        result = original(*args, **kwargs)
        result["accepted"] = False
        return result

    native.side_effect = reject
    assert (
        runner.run_branch(campaign, "forward_acylation")["reason"]
        == "native_window_rejected"
    )
    assert pub.call_count == 1 and ts.call_count == 0
    native.side_effect = original
    monkeypatch.setattr(runner, "candidate", lambda geometry, branch: False)
    assert (
        runner.run_branch(campaign, "forward_acylation")["reason"]
        == "schedule_exhausted_without_candidate"
    )


def test_failed_release_and_exception(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign, native, _, ts = setup_branch(tmp_path, monkeypatch)
    original = native.side_effect

    def release_fails(*args: object, **kwargs: object) -> dict:
        result = original(*args, **kwargs)
        if kwargs["restraints"] is None:
            result["accepted"] = False
        return result

    native.side_effect = release_fails
    assert (
        runner.run_branch(campaign, "reverse_activation")["reason"]
        == "unrestrained_endpoint_not_retained"
    )
    assert ts.call_count == 0
    native.side_effect = RuntimeError("boom")
    with pytest.raises(RuntimeError):
        runner.run_branch(campaign, "reverse_activation")
    assert (
        json.loads((campaign / "reverse/branch_status.json").read_text())["reason"]
        == "exception"
    )
    with pytest.raises(ValueError):
        runner.run_branch(campaign, "unknown")


def test_driven_tolerance(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        runner, "minimization_input", lambda source, maxcycles: "tol=1.d-5, tole=1.d-7"
    )
    monkeypatch.setattr(runner, "driven_input", lambda base, max_cycles: base)
    assert "tol=1.d-4" in runner._driven_mdin("source", 4000)
    monkeypatch.setattr(
        runner, "minimization_input", lambda source, maxcycles: "invalid"
    )
    with pytest.raises(ValueError):
        runner._driven_mdin("source", 4000)
