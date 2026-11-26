"""Preparation helpers for the relay umbrella pilot."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.relay_umbrella_prepare import (
    geometry_features,
    pilot_stage_input,
    pilot_window_specs,
    seed_candidates,
    select_seed,
    write_json,
)


OLD_INPUT = """validated source
&cntrl
 imin=1,
/
&qmmm
 qm_theory='DFTB3', qmcharge=-2, dftb_telec=0.0,
 qmmask='@1-98', scfconv=1.d-10, qmmm_int=1,
/
&dlfind
 active=':1-20,25', neb='FROZEN',
/
"""


def _acceptance(path: Path) -> None:
    restart = path.parent / "optimized.ncrst"
    restart.write_text(f"seed {path.parent.name}\n")
    payload = {
        "restart": str(restart.resolve()),
        "native_converged": True,
        "electronic_history": {"passed": True},
        "geometry": {
            "distances_A": {
                "C_NL": 1.2,
                "C_OG": 1.4,
                "HG_OG": 2.0,
                "HG_N": 1.0,
                "HG_NL": 2.5,
            }
        },
    }
    path.write_text(json.dumps(payload))


def test_geometry_features_and_rejected_seed_receipts() -> None:
    receipt = {
        "native_converged": True,
        "electronic_history": {"passed": True},
        "geometry": {
            "distances_A": {
                "C_NL": 1.2,
                "C_OG": 1.4,
                "HG_OG": 2.0,
                "HG_N": 1.0,
                "HG_NL": 2.5,
            }
        },
    }
    assert geometry_features(receipt) == {
        "q_acyl": pytest.approx(-0.2),
        "q_act": pytest.approx(1.0),
        "q_del": pytest.approx(-1.5),
    }
    with pytest.raises(ValueError, match="native convergence"):
        geometry_features({**receipt, "native_converged": False})
    with pytest.raises(ValueError, match="electronic-history"):
        geometry_features({**receipt, "electronic_history": {"passed": False}})
    with pytest.raises(ValueError, match="lacks required"):
        geometry_features({**receipt, "geometry": {"distances_A": {}}})


def test_seed_candidate_hashes_and_nearest_cv_selection(tmp_path: Path) -> None:
    branch = tmp_path / "branch"
    evidence_dir = branch / "01_coordinate_drive" / "window_001"
    evidence_dir.mkdir(parents=True)
    evidence = evidence_dir / "acceptance.json"
    _acceptance(evidence)
    candidates = seed_candidates(branch)
    selected, distance = select_seed(candidates, ("q_acyl", "q_act"), (-0.2, 1.1))
    assert selected["window"] == "window_001"
    assert selected["sha256"]
    assert distance == pytest.approx(0.1)


def test_empty_or_malformed_seed_candidate_sets_fail(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="No accepted restart"):
        seed_candidates(tmp_path / "empty")
    with pytest.raises(ValueError, match="two CVs"):
        select_seed([], ("q_acyl",), (0.0,))


def test_pilot_grid_uses_acyl_axis_on_both_charts() -> None:
    specs = pilot_window_specs()
    assert len(specs) == 18
    assert sum(name.startswith("activation_") for name, _, _ in specs) == 9
    assert sum(name.startswith("delivery_") for name, _, _ in specs) == 9
    assert {features for _, features, _ in specs} == {
        ("q_acyl", "q_act"),
        ("q_acyl", "q_del"),
    }
    assert ("activation_05", ("q_acyl", "q_act"), (-0.25, 1.10)) in specs
    assert ("delivery_05", ("q_acyl", "q_del"), (-0.25, -0.90)) in specs


@pytest.mark.parametrize(
    "stage,expected_steps",
    [("heat", 2000), ("equil", 2000), ("production", 4000)],
)
def test_pilot_stage_input_reuses_native_builder_with_bounded_steps(
    stage: str, expected_steps: int
) -> None:
    text = pilot_stage_input(
        OLD_INPUT, stage, 12345, Path("umbrella.RST"), active_mask=":1-20,25"
    )
    assert f"nstlim={expected_steps}," in text
    assert "dt=0.0005000" in text
    assert "qmshake=0" in text and "ibelly=1" in text
    assert "qmmask='@1-98'" in text
    assert "irest=0, ntx=1" in text if stage == "heat" else "irest=1, ntx=5" in text


def test_unsupported_pilot_stage_is_explicit() -> None:
    with pytest.raises(ValueError, match="Unsupported"):
        pilot_stage_input(OLD_INPUT, "production30", 1, Path("u.RST"))


def test_write_json_is_finite_and_sorted(tmp_path: Path) -> None:
    path = tmp_path / "record.json"
    write_json(path, {"z": 2, "a": 1})
    assert path.read_text().startswith('{\n  "a": 1,')
    with pytest.raises(ValueError):
        write_json(path, {"invalid": float("nan")})
