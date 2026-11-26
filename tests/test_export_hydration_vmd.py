"""Focused tests for the portable deacylation hydration viewer export."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from parmed.amber import AmberParm

from nyl12_covalent.hydration_vmd_export_io import (
    accepted_extensions,
    cpptraj_input,
    dcd_dimensions,
    selected_indices,
    topology_cpptraj_input,
    trajectory_times,
    water_residues,
)
from nyl12_covalent.hydration_vmd_exporter import export
from nyl12_covalent.hydration_vmd_viewer import viewer_text

SOURCE = Path("qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance")


def test_water_residues_validate_fixed_complete_selection(tmp_path: Path) -> None:
    parm = AmberParm(str(SOURCE / "inputs/system.prmtop"))
    summary = tmp_path / "waters.json"
    summary.write_text(
        json.dumps(
            {
                "waters": [
                    {"water_residue_index_1based": index}
                    for index in (23286, 23285, 23287, 23284, 23288)
                ]
            }
        ),
        encoding="utf-8",
    )
    assert water_residues(summary, parm) == [23286, 23285, 23287, 23284, 23288]


def test_water_residues_reject_duplicate_or_nonwater_ids(tmp_path: Path) -> None:
    parm = AmberParm(str(SOURCE / "inputs/system.prmtop"))
    summary = tmp_path / "waters.json"
    summary.write_text(
        json.dumps({"waters": [{"water_residue_index_1based": 1}] * 5}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unique"):
        water_residues(summary, parm)


def test_selection_retains_qm_atoms_and_exact_reactive_atoms() -> None:
    parm = AmberParm(str(SOURCE / "inputs/system.prmtop"))
    kept = selected_indices(parm, [23286, 23285, 23287, 23284, 23288])
    assert len(kept) == 19967 + 15  # full protein/substrate plus five complete waters
    assert {1620, 13904, 13915, 19915, 19928} <= set(kept)


def test_cpptraj_input_orders_all_segments_and_images_before_strip() -> None:
    text = cpptraj_input(
        Path("system.prmtop"),
        Path("endpoint.ncrst"),
        Path("heat.nc"),
        Path("equil.nc"),
        [23286],
        [Path("extend_01"), Path("extend_02")],
        Path("movie.dcd"),
    )
    assert text.index("trajin endpoint.ncrst") < text.index("trajin heat.nc")
    assert text.index("trajin equil.nc") < text.index("trajin extend_01/trajectory.nc")
    assert text.index("trajin extend_01/trajectory.nc") < text.index(
        "trajin extend_02/trajectory.nc"
    )
    assert text.index("autoimage") < text.index("strip !:")
    topology_text = topology_cpptraj_input(
        Path("system.prmtop"), [23286], Path("stripped.parm7")
    )
    assert "parmstrip !:1-1316,23286" in topology_text
    assert "parmwrite out stripped.parm7 parm7" in topology_text


def test_dcd_dimensions_readback(tmp_path: Path) -> None:
    path = tmp_path / "minimal.dcd"

    def record(data: bytes) -> bytes:
        marker = len(data).to_bytes(4, "little", signed=True)
        return marker + data + marker

    controls = [0] * 20
    controls[0] = 7
    header = b"CORD" + b"".join(
        value.to_bytes(4, "little", signed=True) for value in controls
    )
    path.write_bytes(
        record(header)
        + record(b"title")
        + record((12).to_bytes(4, "little", signed=True))
    )
    assert dcd_dimensions(path) == (7, 12)


def test_native_heat_and_equil_times_are_contiguous_and_chronological() -> None:
    heat = trajectory_times(SOURCE / "heat/trajectory.nc")
    equil = trajectory_times(SOURCE / "equil/trajectory.nc")
    assert (len(heat), heat[0], heat[-1]) == pytest.approx((20, 0.05, 1.0))
    assert (len(equil), equil[0], equil[-1]) == pytest.approx((100, 1.05, 6.0))
    assert heat[-1] < equil[0]


def test_export_bundle_refreshes_atomically_and_retains_prior_manifest(
    tmp_path: Path,
) -> None:
    pytest.importorskip("parmed")
    output = tmp_path / "vmd_progress"
    cpptraj = Path("/data/jpc/code/ambertools26/venv/bin/cpptraj")
    summary = SOURCE / "analysis/per_water_summary.json"
    assert cpptraj.is_file()
    first = export(SOURCE, output, cpptraj, summary)
    prior_manifest = (output / "MANIFEST.sha256").read_text(encoding="utf-8")
    second = export(SOURCE, output, cpptraj, summary)
    assert first == second
    assert first["atoms"] == 19982
    expected_waters = [
        item["water_residue_index_1based"]
        for item in json.loads(summary.read_text(encoding="utf-8"))["waters"][:5]
    ]
    assert first["waters"] == expected_waters
    frame_manifest = json.loads((output / "frame_manifest.json").read_text())
    assert first["frames"] == frame_manifest["frame_count"]
    assert first["frames"] == len(frame_manifest["frame_times_ps"])
    assert frame_manifest["source_time_ps"]["all_frames_chronological"]
    history = json.loads((output / "manifest_history.json").read_text(encoding="utf-8"))
    assert history[-1]["previous_manifest"] == prior_manifest
    assert (tmp_path / "clutter").is_dir()
    assert (output / "frame_manifest.json").is_file()


def test_only_scientifically_accepted_extension_segments_are_selected(
    tmp_path: Path,
) -> None:
    accepted = tmp_path / "extend_01"
    accepted.mkdir()
    (accepted / "trajectory.nc").write_bytes(b"trajectory")
    (accepted / "status.json").write_text(
        json.dumps({"status": "accepted", "exit_code": 0}), encoding="utf-8"
    )
    acceptance = {
        "numeric": {"electronic_history": {"verified_scc_history": True}},
        "boundary": {"box_unchanged": True, "frozen_max_displacement_A": 0.0},
        "acyl_retained": True,
        "scc_verified": True,
        "sampling_converged": False,
    }
    (accepted / "scientific_acceptance.json").write_text(
        json.dumps(acceptance), encoding="utf-8"
    )
    pending = tmp_path / "extend_02"
    pending.mkdir()
    (pending / "status.json").write_text(
        json.dumps({"status": "running"}), encoding="utf-8"
    )
    assert accepted_extensions(tmp_path) == [accepted]
    acceptance["acyl_retained"] = False
    (accepted / "scientific_acceptance.json").write_text(
        json.dumps(acceptance), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="failed validation"):
        accepted_extensions(tmp_path)


def test_viewer_suppresses_labile_edges_and_labels_hydration_geometry() -> None:
    pairs = {
        "OG1_C199": (1, 2),
        "C199_Nleave": (2, 3),
        "OG1_HG1": (1, 4),
        "Nalpha_HG1": (5, 4),
        "Nleave_HG1": (3, 4),
        "O199_Phe113N": (6, 7),
        "O199_Asn186ND2": (6, 8),
        "candidate_O_C199": (9, 2),
        "candidate_O_Nalpha": (9, 5),
    }
    text = viewer_text(
        7,
        12,
        [1, 2, 3],
        [9, 10, 11],
        4,
        pairs,
        (1, 2, 6),
        [1, 2, 9],
        substrate_indices=[2, 3, 6],
        cc_bonds=[(2, 6)],
        chain_groups={},
    )
    assert "::nyl12_movie_policy::apply $top {1 2 3} {2 3 6} {{2 6}}" in text
    assert "VDW 0.28 12" in text
    assert "DynamicBonds $cutoff 0.1 12" in text
    assert "removebond $top 1 4" in text
    assert "removebond $top 2 3" in text
    assert "label add Bonds $top/9 $top/2" in text
    assert "label add Bonds $top/9 $top/5" in text
    assert "label add Angles $top/1 $top/2 $top/6" in text
    assert "source " not in text
    assert "Unexpected atom count" in text
    assert "center_matrix" in text
