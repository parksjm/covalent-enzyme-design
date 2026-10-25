"""Check continuation preserves old windows and stops for orientation review."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _module() -> object:
    path = ROOT / "provenance/asp272_relay_20261025/run_relay_warning_only.py"
    spec = importlib.util.spec_from_file_location("warning_relay", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _prepare(path: Path) -> None:
    original = {
        "state": "stopped",
        "completed_windows": [f"stage1_{i:02d}" for i in range(13)],
        "reason": "stage1_13: convergence=True, coordinates=True, target_tracking=False",
    }
    (path / "status.json").write_text(json.dumps(original))
    inputs = {
        "sha256": {},
        "initial_distances": {"attack": 3.14, "OG_H": 1.0, "Asp_H": 2.0},
    }
    (path / "inputs.json").write_text(json.dumps(inputs))
    folder = path / "stage1_13"
    folder.mkdir()
    (folder / "optimized.ncrst").write_text("preserved restart")
    (folder / "process.json").write_text(
        json.dumps(
            {
                "exit_status": 0,
                "attack_target_A": 1.85,
                "proton_difference_target_A": 0.46,
            }
        )
    )
    (folder / "acceptance.json").write_text(
        json.dumps(
            {
                "native_converged": True,
                "finite_coordinates_and_fixed_box": True,
                "distances_A": {"attack": 2.06},
                "actual_proton_difference_A": 0.41,
                "target_tracking_within_0.20_A": False,
            }
        )
    )


def test_warning_only_continues_phase1_and_holds_phase2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _module()
    _prepare(tmp_path)
    prior = (tmp_path / "stage1_13/acceptance.json").read_bytes()
    calls = []

    def window(
        output: Path,
        label: str,
        attack: float,
        proton: float,
        stage: int,
        restart: Path,
        inputs: dict,
    ) -> Path:
        calls.append((label, stage))
        return output / label / "optimized.ncrst"

    monkeypatch.setattr(module, "window", window)
    monkeypatch.setattr(module, "publish_inspection", lambda *args: None)
    gates = []
    monkeypatch.setattr(
        module, "gate", lambda output, label, stage: gates.append((label, stage))
    )
    module.run(tmp_path)
    assert calls == [("stage1_14", 1), ("stage1_15", 1), ("stage1_16", 1)]
    assert gates == [("stage1_16", 1)]
    state = json.loads((tmp_path / "status.json").read_text())
    assert state["phase1_complete"] is True
    assert state["tracking_policy"]["stop_A"] is None
    assert state["state"] == "stopped"
    assert "orientation review" in state["reason"]
    assert (tmp_path / "stage1_13/acceptance.json").read_bytes() == prior


@pytest.mark.parametrize(
    "flag", ["native_converged", "finite_coordinates_and_fixed_box"]
)
def test_relaxed_tracking_still_requires_valid_convergence(
    tmp_path: Path, flag: str
) -> None:
    module = _module()
    _prepare(tmp_path)
    p = tmp_path / "stage1_13/acceptance.json"
    a = json.loads(p.read_text())
    a[flag] = False
    p.write_text(json.dumps(a))
    with pytest.raises(ValueError, match="does not qualify"):
        module.run(tmp_path)
    assert not (tmp_path / "warning_only_tracking_policy.json").exists()
