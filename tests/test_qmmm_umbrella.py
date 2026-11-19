"""Tests for proton-transfer DFTB3/MM umbrella inputs."""

import json
from pathlib import Path

import pytest

from nyl12_covalent.qmmm_umbrella import (
    asynchronous_umbrella_restraint,
    md_input,
    umbrella_restraint,
)
from scripts.run_qmmm_umbrella_campaign import stage_plan, window_plan

AUDITED = """accepted
&qmmm
 qmmask='@1-3', qmcharge=-2, qm_theory='DFTB3',
/
"""


def test_umbrella_restraint_is_harmonic_difference_coordinate() -> None:
    text = umbrella_restraint(0.25)
    assert "iat=13916,13917,13905,13917" in text
    assert "rstwt=1,-1,0,0" in text
    assert "r2=0.25000000, r3=0.25000000" in text
    assert "rk2=50.0, rk3=50.0" in text


def test_asynchronous_umbrella_uses_one_three_distance_coordinate() -> None:
    text = asynchronous_umbrella_restraint(-2.5)
    assert "iat=13916,13917,13905,13917,13916,19916,0,0" in text
    assert "rstwt=1,-1,-1,0" in text
    assert "r2=-2.50000000, r3=-2.50000000" in text
    assert "rk2=50.0, rk3=50.0" in text


def test_heating_input_is_nvt_and_records_coordinate_frequently() -> None:
    text = md_input(AUDITED, "heating", 2000)
    assert "irest=0, ntx=1, tempi=25.0" in text
    assert "value1=25.0, value2=300.0" in text
    assert "type='DUMPFREQ', istep1=10" in text
    assert "ntr=0, nmropt=1" in text


def test_production_restarts_and_has_no_temperature_ramp() -> None:
    text = md_input(AUDITED, "production", 4000)
    assert "irest=1, ntx=5, tempi=300.0" in text
    assert "TEMP0" not in text
    assert "ntwx=100" in text
    assert "restraint_wt" not in text


def test_asynchronous_window_plan_includes_es_and_driven_seeds(
    tmp_path: Path,
) -> None:
    path_run = tmp_path / "path"
    drive = path_run / "01_coupled_coordinate_drive"
    drive.mkdir(parents=True)
    es_start = tmp_path / "es.ncrst"
    es_start.write_text("ES")
    schedule = [
        {"label": "asynchronous_01", "composite_target_A": -2.0},
        {"label": "asynchronous_02", "composite_target_A": -1.0},
    ]
    for row in schedule:
        folder = drive / row["label"]
        folder.mkdir()
        (folder / "optimized.ncrst").write_text(row["label"])
    (path_run / "status.json").write_text(
        json.dumps({"completed_windows": [row["label"] for row in schedule]})
    )
    (path_run / "inputs.json").write_text(
        json.dumps(
            {
                "initial_geometry_A": {
                    "proton_difference_A": -1.0,
                    "OG1_carbonyl_C_A": 3.0,
                },
                "schedule": schedule,
            }
        )
    )
    rows = window_plan(path_run, es_start, "asynchronous")
    assert [row["label"] for row in rows] == [
        "window_00",
        "window_01",
        "window_02",
    ]
    assert [row["target_A"] for row in rows] == [-4.0, -2.0, -1.0]
    assert all(row["start"].is_file() for row in rows)


def test_production30_stage_plan_is_segmented() -> None:
    """Full production uses heating, equilibration, and six restart blocks."""
    stages = stage_plan("production30")
    assert stages[0] == ("01_heating_1ps", "heating", 2_000)
    assert stages[1] == ("02_equilibration_5ps", "equilibration", 10_000)
    assert stages[-1] == ("03_production_block_06_5ps", "production", 10_000)
    assert len(stages) == 8


def test_asynchronous_extension_window_plan_uses_only_new_seeds(
    tmp_path: Path,
) -> None:
    path_run = tmp_path / "path"
    drive = path_run / "03_coordinate_drive_extension_to_ti1"
    drive.mkdir(parents=True)
    schedule = [
        {"label": "asynchronous_25", "composite_target_A": -0.65},
        {"label": "asynchronous_26", "composite_target_A": -0.50},
    ]
    for row in schedule:
        folder = drive / row["label"]
        folder.mkdir()
        (folder / "optimized.ncrst").write_text(row["label"])
    (path_run / "extension_inputs.json").write_text(json.dumps({"schedule": schedule}))
    (path_run / "extension_status.json").write_text(
        json.dumps(
            {
                "state": "completed",
                "completed_windows": [row["label"] for row in schedule],
            }
        )
    )
    rows = window_plan(path_run, tmp_path / "unused.ncrst", "asynchronous_extension")
    assert [row["label"] for row in rows] == ["window_25", "window_26"]
    assert [row["target_A"] for row in rows] == [-0.65, -0.50]
    assert all(row["start"].is_file() for row in rows)


def test_stage_plan_rejects_unknown_protocol() -> None:
    """Unknown schedules fail rather than silently falling back."""
    with pytest.raises(ValueError, match="Unsupported umbrella protocol"):
        stage_plan("unknown")
