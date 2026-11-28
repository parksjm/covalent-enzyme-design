"""Exact six-check target contracts and shared preparation/runtime entrypoints."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from test_relay_targeted import plan as plan
from nyl12_covalent import relay_mixing_checks as mixing
from nyl12_covalent.neb_umbrella_campaign import digest


def mixing_plan(plan: Path) -> Path:
    """Adapt the validated fixture to six exact parent main-window targets."""
    data = json.loads(plan.read_text())
    base = data["targets"][0]
    parent_path = Path(data["parent_prepared"])
    parent = json.loads(parent_path.read_text())
    parent["windows"] = [
        dict(
            label=label,
            bias=dict(
                feature_names=base["feature_names"],
                centers_A=base["centers_A"],
                force_constants=[25.0, 25.0],
            ),
        )
        for label in set(mixing.TARGETS.values())
    ]
    parent_path.write_text(json.dumps(parent))
    parent_path.with_suffix(".sha256").write_text(digest(parent_path))
    data["output"] = str(plan.parent / "mixing_checks_20261127")
    data["targets"] = [
        {**base, "label": label, "initial_ig": 1000 + 10 * i}
        for i, label in enumerate(mixing.TARGETS)
    ]
    plan.write_text(json.dumps(data))
    return plan


def test_preparation_and_runtime_delegate(
    plan: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Six inputs retain exact target bias, fresh velocities and eight stages."""
    path = mixing_plan(plan)
    prepared = mixing.prepare_mixing_checks(path)
    cfg = json.loads(prepared.read_text())
    assert len(cfg["windows"]) == 6
    assert all(len(r["stages"]) == 8 for r in cfg["windows"])
    first = prepared.parent / cfg["windows"][0]["stages"][0]["input"]
    assert "irest=0, ntx=1, tempi=300.0" in first.read_text()
    execute = Mock()
    monkeypatch.setattr(mixing, "run_reviewed_target", execute)
    mixing.run_mixing_check(prepared, "A026_original", digest(prepared))
    execute.assert_called_once_with(prepared, "A026_original", digest(prepared), 6)
    Path(cfg["parent_prepared"]).write_text("{}")
    with pytest.raises(ValueError, match="Parent hash"):
        mixing.prepare_mixing_checks(path)


@pytest.mark.parametrize("change", ["missing", "label", "center", "feature", "force"])
def test_reject_changed_native_contract(plan: Path, change: str) -> None:
    """Missing targets or changed native Hamiltonian fail explicitly."""
    path = mixing_plan(plan)
    data = json.loads(path.read_text())
    if change == "missing":
        data["targets"].pop()
    elif change == "label":
        data["targets"][0]["label"] = "unexpected"
    elif change == "center":
        data["targets"][0]["centers_A"] = [0, 0]
    elif change == "feature":
        data["targets"][0]["feature_names"] = ["q_act", "q_del"]
    else:
        data["targets"][0]["k"] = 50
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        mixing.prepare_mixing_checks(path)


@pytest.mark.parametrize(
    "args,operation",
    [
        (["--plan", "p"], "prepare"),
        (["--prepared", "p", "--label", "A026_late", "--sha256", "h"], "run"),
        ([], "error"),
    ],
)
def test_cli(monkeypatch: pytest.MonkeyPatch, args: list[str], operation: str) -> None:
    """CLI requires one complete mode and delegates exact requested arguments."""
    import sys

    prepare = Mock()
    run = Mock()
    monkeypatch.setattr(mixing, "prepare_mixing_checks", prepare)
    monkeypatch.setattr(mixing, "run_mixing_check", run)
    monkeypatch.setattr(sys, "argv", ["mixing", *args])
    if operation == "error":
        with pytest.raises(SystemExit):
            mixing.main()
    else:
        mixing.main()
        assert (prepare if operation == "prepare" else run).call_count == 1
