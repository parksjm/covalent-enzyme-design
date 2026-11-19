from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.umbrella_mbar import (
    load_umbrella_sample_blocks,
    load_umbrella_samples,
    run_mbar_smoke,
)


def _campaign(tmp_path: Path, *, complete: bool) -> Path:
    campaign = tmp_path / "campaign"
    windows = []
    for index, target in enumerate((-0.5, 0.0, 0.5)):
        label = f"window_{index:02d}"
        windows.append({"label": label, "target_A": target})
        stage = campaign / label / "production"
        stage.mkdir(parents=True)
        samples = np.linspace(target - 0.35, target + 0.35, 80)
        (stage / "umbrella_coordinate.dat").write_text(
            "".join(
                f"{step * 10:8d} {value:10.5f}\n" for step, value in enumerate(samples)
            )
        )
        if complete:
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


def test_load_requires_completed_process_by_default(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path, complete=False)
    with pytest.raises(ValueError, match="Missing process record"):
        load_umbrella_samples(
            campaign,
            "production",
            require_successful_process=True,
        )
    inputs, windows = load_umbrella_samples(
        campaign,
        "production",
        require_successful_process=False,
    )
    assert inputs["temperature_K"] == 300.0
    assert len(windows) == 3
    assert windows[0]["samples_A"].size == 80


def test_mbar_smoke_returns_finite_diagnostics(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path, complete=True)
    _, windows = load_umbrella_samples(
        campaign,
        "production",
        require_successful_process=True,
    )
    result = run_mbar_smoke(
        windows,
        temperature_k=300.0,
        force_constant_kcal_mol_a2=10.0,
        bins=24,
    )
    assert result["window_count"] == 3
    assert result["raw_sample_count"] == 240
    assert 0 < result["decorrelated_sample_count"] <= 240
    assert np.asarray(result["mbar_overlap_matrix"]).shape == (3, 3)
    assert np.isfinite(result["mbar_overlap_scalar"])
    assert min(row["free_energy_kcal_mol"] for row in result["provisional_pmf"]) == 0.0
    assert all(
        np.isfinite(row["uncertainty_kcal_mol"]) for row in result["provisional_pmf"]
    )


def test_load_and_decorrelate_multiple_blocks(tmp_path: Path) -> None:
    """Combine matching stages while decorrelating each block independently."""
    first = _campaign(tmp_path / "first", complete=True)
    second = _campaign(tmp_path / "second", complete=True)
    inputs, windows = load_umbrella_sample_blocks(
        [(first, "production"), (second, "production")],
        require_successful_process=True,
    )
    assert inputs["temperature_K"] == 300.0
    assert len(windows) == 3
    assert windows[0]["samples_A"].size == 160
    assert len(windows[0]["sample_blocks_A"]) == 2

    result = run_mbar_smoke(
        windows,
        temperature_k=300.0,
        force_constant_kcal_mol_a2=10.0,
        bins=24,
    )
    assert result["raw_sample_count"] == 480
    assert all(len(window["blocks"]) == 2 for window in result["windows"])


def test_combined_loader_rejects_different_grid(tmp_path: Path) -> None:
    """Reject blocks that do not share exactly the same target grid."""
    first = _campaign(tmp_path / "first", complete=True)
    second = _campaign(tmp_path / "second", complete=True)
    inputs_path = second / "inputs.json"
    inputs = json.loads(inputs_path.read_text())
    inputs["windows"][0]["target_A"] = -0.4
    inputs_path.write_text(json.dumps(inputs))
    with pytest.raises(ValueError, match="window grids differ"):
        load_umbrella_sample_blocks(
            [(first, "production"), (second, "production")],
            require_successful_process=True,
        )


@pytest.mark.parametrize(
    ("temperature", "force_constant", "bins", "message"),
    [
        (0.0, 10.0, 24, "positive"),
        (300.0, 0.0, 24, "positive"),
        (300.0, 10.0, 3, "four"),
    ],
)
def test_mbar_smoke_rejects_invalid_controls(
    tmp_path: Path,
    temperature: float,
    force_constant: float,
    bins: int,
    message: str,
) -> None:
    campaign = _campaign(tmp_path, complete=True)
    _, windows = load_umbrella_samples(
        campaign,
        "production",
        require_successful_process=True,
    )
    with pytest.raises(ValueError, match=message):
        run_mbar_smoke(
            windows,
            temperature_k=temperature,
            force_constant_kcal_mol_a2=force_constant,
            bins=bins,
        )
