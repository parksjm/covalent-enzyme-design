"""Reject stale prefixes and unverifiable source claims in supplemental movies."""

from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

SOURCE = (
    Path(__file__).resolve().parents[1]
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/supplemental_movies_post5ps_prepare_20261127/validate_source_selection.py"
)
spec = importlib.util.spec_from_file_location("movie_source_selection", SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fixture() -> tuple[dict, dict, dict, dict]:
    source = "/campaign/A/04_production_5ps/trajectory.nc"
    ref = {
        "original_atom_ids_1based": list(range(1, 19968)),
        "qm_atom_ids_1based": list(range(1, 99)),
        "roles_1based": {"carbon": 1},
        "parameter_flag_sha256": "parameters",
        "source_topology_path": "/full/source.parm7",
        "source_topology_sha256": "fulltopology",
    }
    part = {
        "campaign": "/campaign",
        "expected_frames": 100,
        "full_horizon_ps": 10,
        "selected_blocks": [{"window": "A", "blocks": ["04_production_5ps"]}],
        "source_observables": [
            {
                "window": "A",
                "block": "04_production_5ps",
                "trajectory_sha256": "accepted-source",
            }
        ],
    }
    rows = [
        {
            "window": "A",
            "block": "04_production_5ps",
            "source_frame1based": i,
            "output_frame0based": i - 1,
            "source": source,
            "window_time_ps": 50 + i * 0.05,
        }
        for i in range(1, 101)
    ]
    manifest = {
        **ref,
        "frames": 100,
        "atoms": 19967,
        "frame_plan": rows,
        "source_sha256": {
            source: "accepted-source",
            "/full/source.parm7": "fulltopology",
        },
    }
    return manifest, part, ref, {source: [50 + i * 0.05 for i in range(1, 101)]}


def test_all_frames_postdiscard() -> None:
    records = m.verify_selection(*fixture())
    assert len(records) == 100
    assert records[0]["MBAR_production_time_ps"] == 5.05
    assert records[-1]["MBAR_production_time_ps"] == 10


@pytest.mark.parametrize(
    "change", ["prefix", "count", "map", "hash", "time", "order", "future"]
)
def test_reject_invalid_actual_output(change: str) -> None:
    manifest, part, ref, times = deepcopy(fixture())
    if change == "prefix":
        manifest["frame_plan"][0]["block"] = "03_production_5ps"
    elif change == "count":
        manifest["frames"] = 99
    elif change == "map":
        manifest["qm_atom_ids_1based"] = [2]
    elif change == "hash":
        manifest["source_sha256"] = {}
    elif change == "time":
        manifest["frame_plan"][0]["window_time_ps"] = 0
    elif change == "order":
        manifest["frame_plan"][0]["output_frame0based"] = 1
    elif change == "future":
        part["full_horizon_ps"] = 9
    with pytest.raises(ValueError):
        m.verify_selection(manifest, part, ref, times)


def test_truncated_frame_plan() -> None:
    manifest, part, ref, times = fixture()
    manifest["frame_plan"].pop()
    with pytest.raises(ValueError):
        m.verify_selection(manifest, part, ref, times)


def test_style_policy() -> None:
    m.verify_style(
        "Licorice 0.2 Q orange display depthcue off mol showrep $m [expr {[molinfo $m get numreps] - 1}] off"
    )
    with pytest.raises(ValueError):
        m.verify_style("incorrect")


def test_Amber_flag_metadata_only() -> None:
    reference = {"TITLE": ["original"], "CHARGE": [1.0], "MASS": [12.0]}
    actual = {**reference, "TITLE": ["display title"]}
    result = m.verify_parameter_flags(reference, actual)
    assert result["all_Amber_parameter_flags_equal"]
    assert list(result["excluded_metadata_differences"]) == ["TITLE"]
    assert (
        m.verify_parameter_flags(reference, reference)["parameter_flag_sha256"]
        == result["parameter_flag_sha256"]
    )


@pytest.mark.parametrize(
    "flag",
    [
        "CHARGE",
        "MASS",
        "AMBER_ATOM_TYPE",
        "BONDS_INC_HYDROGEN",
        "ANGLE_FORCE_CONSTANT",
        "POINTERS",
    ],
)
def test_any_real_parameter_change_rejected(flag: str) -> None:
    with pytest.raises(ValueError, match="parameter data"):
        m.verify_parameter_flags({flag: [1]}, {flag: [2]})


def test_parameter_flag_removal_rejected() -> None:
    with pytest.raises(ValueError, match="flag set"):
        m.verify_parameter_flags({"CHARGE": [1]}, {})
