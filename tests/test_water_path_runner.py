"""Focused tests for water-path continuity and frozen-window resume checks."""

from __future__ import annotations

import importlib
import json
from pathlib import Path
import sys
from typing import Any

import pytest

RUNNER_DIR = (
    Path(__file__).resolve().parents[1]
    / "qmmm/pa6_flipped_model30/07_ash274_comparison"
    / "07_water_to_asp272_ti1_like"
)
sys.path.insert(0, str(RUNNER_DIR))
runner = importlib.import_module("run_water_to_asp272_path")


def write_json(path: Path, value: dict[str, Any]) -> None:
    """Write a compact JSON fixture file."""
    path.write_text(json.dumps(value))


def test_continuity_ignores_spectator_but_retains_complete_delta_map() -> None:
    """A spectator movement does not fail the reactive-distance gate."""
    changes = {key: 0.1 for key in runner.DISTANCE_KEYS}
    changes["ash274_hd2_water_o_A"] = 0.4448
    assert not runner.continuity_passes({**changes, "og1_substrate_c_A": 0.35001})
    assert runner.continuity_passes(changes)
    assert len(changes) == len(runner.DISTANCE_KEYS)


def test_resume_history_validates_frozen_outputs_and_start_restart(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Resume accepts exact provenance and rejects an altered window output."""
    monkeypatch.setattr(runner, "CASE", tmp_path)
    for index in range(0, 6):
        window = tmp_path / f"window_{index:03d}"
        window.mkdir()
        (window / "optimized.ncrst").write_bytes(f"restart-{index}".encode())
    initial_geometry = {key: 1.0 for key in runner.DISTANCE_KEYS}
    protocol: dict[str, Any] = {
        "status": "path_running",
        "last_accepted_window": "window_004",
        "coordinates": {"starting_geometry": initial_geometry},
    }
    output_names = ("minimized.mdout", "optimized.ncrst", "optimization.nc")
    for index in range(1, 6):
        window = tmp_path / f"window_{index:03d}"
        for name in output_names:
            content = (
                b"Energy calculation finished\nConvergence reached\nOptimisation converged\n"
                if name == "minimized.mdout"
                else f"{index}:{name}".encode()
            )
            (window / name).write_bytes(content)
        previous_restart = tmp_path / f"window_{index - 1:03d}" / "optimized.ncrst"
        hashes = {name: runner.sha256(window / name) for name in output_names}
        write_json(
            window / "process.json",
            {
                "status": "finished",
                "exit_status": 0,
                "starting_restart_sha256": runner.sha256(previous_restart),
                "drive_mdin_sha256": runner.sha256(window / "drive.mdin")
                if (window / "drive.mdin").exists()
                else "unused",
                "restraints_sha256": runner.sha256(window / "restraints.RST")
                if (window / "restraints.RST").exists()
                else "unused",
            },
        )
        if index == 5:
            (window / "starting.ncrst").write_bytes(previous_restart.read_bytes())
            (window / "drive.mdin").write_bytes(b"drive")
            (window / "restraints.RST").write_bytes(b"restraints")
            process = json.loads((window / "process.json").read_text())
            process["drive_mdin_sha256"] = runner.sha256(window / "drive.mdin")
            process["restraints_sha256"] = runner.sha256(window / "restraints.RST")
            write_json(window / "process.json", process)
            (window / "acceptance.json").write_text(
                json.dumps(
                    {
                        "accepted": False,
                        "restart_sha256": runner.sha256(window / "optimized.ncrst"),
                        "output_sha256": hashes,
                    }
                )
            )
            (window / "window.json").write_text("{}")
        else:
            (window / "acceptance.json").write_text(
                json.dumps(
                    {
                        "window": f"window_{index:03d}",
                        "accepted": True,
                        "geometry": initial_geometry,
                        "energy": initial_energy_fixture(),
                        "restart_sha256": runner.sha256(window / "optimized.ncrst"),
                        "output_sha256": hashes,
                    }
                )
            )
    history = runner.resume_history(protocol, initial_energy_fixture())
    assert [row["window"] for row in history] == [
        "window_000",
        "window_001",
        "window_002",
        "window_003",
        "window_004",
    ]

    window_005 = tmp_path / "window_005"
    original = (window_005 / "acceptance.json").read_bytes()
    (window_005 / "acceptance_revision_1.json").write_bytes(original)
    write_json(
        window_005 / "acceptance.json",
        {
            "window": "window_005",
            "accepted": True,
            "re_evaluated_under_protocol_revision": 5,
            "restart_sha256": runner.sha256(window_005 / "optimized.ncrst"),
            "output_sha256": {
                name: runner.sha256(window_005 / name) for name in output_names
            },
        },
    )
    protocol["last_accepted_window"] = "window_005"
    resumed_history = runner.resume_history(protocol, initial_energy_fixture())
    assert resumed_history[-1]["window"] == "window_005"

    (tmp_path / "window_002" / "minimized.mdout").write_bytes(b"changed")
    with pytest.raises(ValueError, match="output hash mismatch"):
        runner.resume_history(protocol, initial_energy_fixture())


def initial_energy_fixture() -> dict[str, float]:
    """Return the minimal energy fields needed to rebuild resume history."""
    return {
        "physical_kcal_mol": 0.0,
        "EPtot_bias_inclusive_kcal_mol": 0.0,
    }
