"""Tests for the multi-block PyMBAR analysis entry point."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType

import numpy as np


def _load_script() -> ModuleType:
    """Load the MBAR analysis entry point as a test module."""
    path = Path(__file__).parents[1] / "scripts/run_mbar_analysis.py"
    spec = importlib.util.spec_from_file_location("test_run_mbar_analysis", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _campaign(root: Path, offset: float) -> Path:
    """Create one complete synthetic three-window production block."""
    campaign = root / "campaign"
    windows = []
    for index, target in enumerate((-0.5, 0.0, 0.5)):
        label = f"window_{index:02d}"
        windows.append({"label": label, "target_A": target})
        stage = campaign / label / "production"
        stage.mkdir(parents=True)
        samples = np.linspace(target - 0.35, target + 0.35, 80) + offset
        (stage / "umbrella_coordinate.dat").write_text(
            "".join(
                f"{step * 10:8d} {value:10.5f}\n" for step, value in enumerate(samples)
            )
        )
        (stage / "process.json").write_text(json.dumps({"exit_status": 0}))
    (campaign / "inputs.json").write_text(
        json.dumps(
            {
                "temperature_K": 300.0,
                "force_constant_kcal_mol_A2": 10.0,
                "windows": windows,
            }
        )
    )
    return campaign


def test_profile_change_reports_finite_differences() -> None:
    """Compare profiles on a common interpolated grid."""
    module = _load_script()
    previous = [
        {"coordinate_A": value, "free_energy_kcal_mol": value * value}
        for value in np.linspace(-1.0, 1.0, 10)
    ]
    current = [
        {"coordinate_A": value, "free_energy_kcal_mol": 1.1 * value * value}
        for value in np.linspace(-0.9, 0.9, 12)
    ]
    result = module._profile_change(previous, current)
    assert result["common_point_count"] >= 4
    assert result["rms_change_kcal_mol"] > 0.0


def test_analyze_writes_hashed_multiblock_outputs(
    tmp_path: Path, monkeypatch: object
) -> None:
    """Run cumulative MBAR and require the declared result bundle."""
    module = _load_script()
    first = _campaign(tmp_path / "first", 0.0)
    second = _campaign(tmp_path / "second", 0.01)

    def fake_run(command: list[str], **_: object) -> None:
        Path(command[3]).write_bytes(b"plot")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    output = tmp_path / "analysis"
    module.analyze(
        [(first, "production"), (second, "production")],
        output,
        24,
        Path("python"),
        "test PMF",
    )
    summary = json.loads((output / "mbar_summary.json").read_text())
    assert summary["raw_sample_count"] == 480
    assert len(summary["cumulative_convergence"]) == 2
    assert "change_from_previous" in summary["cumulative_convergence"][1]
    assert (output / "pmf.csv").is_file()
    assert (output / "mbar_overlap_matrix.csv").is_file()
    assert (output / "MANIFEST.sha256").is_file()
