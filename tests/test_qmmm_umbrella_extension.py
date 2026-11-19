"""Tests for segmented QM/MM umbrella production extensions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from scripts.run_qmmm_umbrella_extension import (
    extension_stages,
    extension_windows,
    run,
)


def _pilot(tmp_path: Path, coordinate: str = "direct_pt") -> Path:
    pilot = tmp_path / "pilot"
    restart = pilot / "window_00/03_production_2ps/restart.ncrst"
    restart.parent.mkdir(parents=True)
    restart.write_text("restart\n")
    reaction_coordinate = {
        "direct_pt": "d(OG1,HG1)-d(N,HG1)",
        "asynchronous": "d(OG1,HG1)-d(N,HG1)-d(OG1,carbonyl_C)",
    }[coordinate]
    (pilot / "inputs.json").write_text(
        json.dumps(
            {
                "reaction_coordinate": reaction_coordinate,
                "windows": [{"label": "window_00", "target_A": -1.0}],
            }
        )
    )
    (pilot / "status.json").write_text(
        json.dumps(
            {
                "state": "completed",
                "failures": [],
                "windows": {"window_00": {"state": "completed"}},
            }
        )
    )
    return pilot


def test_extension_stages_are_five_ps_blocks() -> None:
    """The default schedule contains 5 ps equilibration plus 30 ps production."""
    stages = extension_stages(5, 30, 5)
    assert stages[0] == ("01_equilibration_5ps", "equilibration", 10_000)
    assert len(stages) == 7
    assert stages[-1] == ("02_production_block_06_5ps", "production", 10_000)


def test_extension_stages_reject_nondivisible_production() -> None:
    """Partial blocks are rejected rather than silently shortened."""
    with pytest.raises(ValueError, match="divisible"):
        extension_stages(5, 31, 5)


def test_extension_windows_requires_completed_matching_pilot(tmp_path: Path) -> None:
    """Only completed pilots with the declared coordinate can seed production."""
    pilot = _pilot(tmp_path)
    rows = extension_windows(pilot, "direct_pt")
    assert rows[0]["start"].name == "restart.ncrst"
    with pytest.raises(ValueError, match="identity"):
        extension_windows(pilot, "asynchronous")
    status_path = pilot / "status.json"
    status = json.loads(status_path.read_text())
    status["state"] = "running"
    status_path.write_text(json.dumps(status))
    with pytest.raises(ValueError, match="not completed"):
        extension_windows(pilot, "direct_pt")


def test_run_writes_segmented_campaign_with_fake_sander(tmp_path: Path) -> None:
    """The controller executes each stage and records successful provenance."""
    pilot = _pilot(tmp_path)
    topology = tmp_path / "system.parm7"
    topology.write_text("topology\n")
    audited = tmp_path / "audited.mdin"
    audited.write_text("&qmmm\n qm_theory='DFTB3', qmcharge=-2,\n/\n")
    sander = tmp_path / "sander"
    sander.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "for flag, text in [('-r', 'restart\\n'), ('-x', 'trajectory\\n'), "
        "('-o', 'output\\n'), ('-inf', 'info\\n')]:\n"
        "    pathlib.Path(args[args.index(flag) + 1]).write_text(text)\n"
    )
    sander.chmod(0o750)
    output = tmp_path / "extension"
    run(
        argparse.Namespace(
            pilot=pilot,
            topology=topology,
            audited_mdin=audited,
            sander=sander,
            output=output,
            workers=1,
            equilibration_ps=1,
            production_ps=1,
            block_ps=1,
            coordinate="direct_pt",
        )
    )
    status = json.loads((output / "status.json").read_text())
    assert status["state"] == "completed"
    window = output / "window_00"
    assert (window / "01_equilibration_1ps/restart.ncrst").is_file()
    assert (window / "02_production_block_01_1ps/restart.ncrst").is_file()
    process = json.loads(
        (window / "02_production_block_01_1ps/process.json").read_text()
    )
    assert process["exit_status"] == 0
    assert "restart_sha256" in process
