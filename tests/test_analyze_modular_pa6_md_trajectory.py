"""Tests for modular PA6 trajectory geometry analysis."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts/analyze_modular_pa6_md_trajectory.py"
SPEC = importlib.util.spec_from_file_location(
    "analyze_modular_pa6_md_trajectory", SCRIPT
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_cpptraj_uses_modular_scissile_mapping(tmp_path: Path) -> None:
    text = MODULE.cpptraj_text(Path("top.prmtop"), Path("prod.nc"), tmp_path, 700)
    assert "trajin" in text and " 1 700 1" in text
    assert "distance d1 :929@OG1 :1314@C" in text
    assert "distance d2 :113@N :1314@O" in text
    assert ":1315@N" in text and ":1314@C5" in text
    assert "rms ligand :1313-1316&!@H= first nofit" in text


def test_summary_applies_all_five_metrics_and_detects_drift() -> None:
    favorable = {
        "frame": 1.0,
        "time_ns": 0.0,
        "d1": 3.1,
        "d2": 3.7,
        "d3": 3.7,
        "attack": 100.0,
        "projected_phi": 90.0,
        "ligand_rmsd": 1.0,
    }
    unfavorable = {
        **favorable,
        "frame": 2.0,
        "time_ns": 0.05,
        "d1": 6.0,
        "ligand_rmsd": 8.0,
    }
    assert MODULE.passes(favorable, MODULE.TIERS["tier_a"])
    assert not MODULE.passes(unfavorable, MODULE.TIERS["tier_c"])
    summary = MODULE.summarize([favorable, unfavorable])
    assert summary["tier_c_fraction"] == 0.5
    assert summary["tier_c_first_half"] == 1.0
    assert summary["tier_c_second_half"] == 0.0
    assert not summary["interim_persistence_rule_passed"]
