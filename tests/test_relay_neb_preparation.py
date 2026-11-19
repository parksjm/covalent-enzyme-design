"""Trace order, atom-labeled proton gates and native NEB initialization."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.relay_neb_native import native_initial_band, native_neb_input
from nyl12_covalent.relay_neb_trace import (
    proton_audit,
    resample_trace,
    stage_trace_paths,
)


def test_qm_arc_ignores_dominating_mm_motion_and_retains_curve() -> None:
    frames = np.zeros((3, 2, 3))
    frames[:, 0, :2] = [[0, 0], [1, 0], [1, 1]]
    frames[:, 1, 0] = [0, 1000, 1001]
    controls, brackets, arcs = resample_trace(frames, np.array([0]), 5)
    np.testing.assert_allclose(
        controls[:, 0, :2], [[0, 0], [0.5, 0], [1, 0], [1, 0.5], [1, 1]]
    )
    assert arcs["qm_arc_A"] == 2
    assert arcs["active_arc_A"] > 1000
    assert brackets[1]["fraction"] == 0.5
    np.testing.assert_array_equal(controls[[0, -1]], frames[[0, -1]])


@pytest.mark.parametrize("bad", ["count", "shape", "nan", "indices", "zero"])
def test_resampling_rejects_invalid_trace(bad: str) -> None:
    frames = np.zeros((3, 2, 3))
    frames[:, 0, 0] = [0, 1, 2]
    count, ids = 3, np.array([0])
    if bad == "count":
        count = 1
    elif bad == "shape":
        frames = frames[:, :, 0]
    elif bad == "nan":
        frames[0, 0, 0] = np.nan
    elif bad == "indices":
        ids = np.array([2])
    elif bad == "zero":
        frames[:] = 0
    with pytest.raises(ValueError):
        resample_trace(frames, ids, count)


def controls() -> np.ndarray:
    """Return 13 non-rigid 3D controls with a removable global translation."""
    base = np.array([[0.0, 0, 0], [10, 0, 0], [0, 10, 0], [0, 0, 10], [1, 2, 3]])
    frames = np.repeat(base[None], 13, axis=0)
    for i in range(13):
        frames[i] += [i * 0.1, i * 0.2, -0.1 * i]
        frames[i, 4, 0] += i * 0.05
    return frames


def test_native_alignment_slots_endpoints_and_regular_count() -> None:
    frames = controls()
    band, report = native_initial_band(frames, np.arange(5))
    assert band.shape == (24, 5, 3)
    slots = report["source_controls_to_regular_slots_1based"]
    assert slots[0] == 1 and slots[-1] == 24
    assert slots == sorted(set(slots))
    assert report["dropped_controls_zero_based"] == []
    assert report["max_active_alignment_shift_A"] > 1
    np.testing.assert_array_equal(band[0], frames[0])
    x = frames[-1] - frames[-1].mean(axis=0)
    y = frames[0] - frames[0].mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    if np.linalg.det(u @ vt) < 0:
        u[:, -1] *= -1
    expected_aligned_last = x @ (u @ vt) + frames[0].mean(axis=0)
    np.testing.assert_allclose(band[-1], expected_aligned_last, atol=1e-12)
    assert report["native_slot24_exact_aligned_last_control"] is True

    for i, slot in enumerate(slots):
        np.testing.assert_allclose(
            np.linalg.norm(band[slot - 1, 0] - band[slot - 1, 4]),
            np.linalg.norm(frames[i, 0] - frames[i, 4]),
            atol=1e-12,
        )


@pytest.mark.parametrize("bad", ["count", "nan", "indices", "zero"])
def test_native_initializer_rejects_invalid_controls(bad: str) -> None:
    frames, ids = controls(), np.arange(5)
    if bad == "count":
        frames = frames[:-1]
    elif bad == "nan":
        frames[0, 0, 0] = np.nan
    elif bad == "indices":
        ids = np.array([20])
    elif bad == "zero":
        frames[:] = frames[0]
    with pytest.raises(ValueError):
        native_initial_band(frames, ids)


BASE = "&cntrl nmropt=0,ntr=0,maxcyc=1 /\n&qmmm qmmask='@1-102',qmcharge=-1,dftb_telec=0.0 /\n&dlfind active='@1-4662',maxcycle=1 /\n"


def test_native_input_preserves_hamiltonian_and_all_active_forces() -> None:
    images = [Path(f"/new/control_{i:02d}.ncrst") for i in range(13)]
    text = native_neb_input(BASE, images)
    assert "wtmask='@1-102'" in text and "active='@1-4662'" in text
    assert "qmcharge=-1,dftb_telec=0.0" in text
    assert text.count("coords2(") == 12
    assert "neb_nimage=25" in text and "neb='FROZEN'" in text
    assert "maxcycle=16000" in text and "maxcyc=100000" in text
    assert "maxeneval" not in text and "DISANG" not in text


@pytest.mark.parametrize("bad", ["count", "restraint", "nmropt", "mask", "cycles"])
def test_native_input_rejects_bad_sources(bad: str) -> None:
    images = [Path(f"/i{i}") for i in range(13)]
    text = BASE
    if bad == "count":
        images = images[:-1]
    elif bad == "restraint":
        text += "DISANG=anything"
    elif bad == "nmropt":
        text = text.replace("nmropt=0", "nmropt=1")
    elif bad == "mask":
        text = text.replace("active='@1-4662'", "")
    elif bad == "cycles":
        text = text.replace("maxcyc=1", "")
    with pytest.raises(ValueError):
        native_neb_input(text, images)


@pytest.mark.parametrize("bad", [None, "far_H", "compressed_HH", "acid"])
def test_proton_audit_retains_labels_and_rejects_bad_geometry(bad: str | None) -> None:
    xyz = np.array([[0.0, 0, 0], [4, 0, 0], [1, 0, 0], [3, 0, 0]])
    atoms = [SimpleNamespace(atomic_number=z) for z in (8, 8, 1, 1)]
    model = SimpleNamespace(
        qm=np.arange(4),
        hydrogens=np.array([2, 3]),
        topology=SimpleNamespace(atoms=atoms),
        atom_map={"ash_proton": 3, "ash_oxygen": 1},
    )
    if bad == "far_H":
        xyz[3] = [10, 0, 0]
    elif bad == "compressed_HH":
        xyz[3] = [1.1, 0, 0]
    elif bad == "acid":
        model.atom_map["ash_oxygen"] = 2
    if bad:
        with pytest.raises(ValueError):
            proton_audit(model, xyz)
    else:
        row = proton_audit(model, xyz)
        assert row["hydrogen_parents_1based"] == {"3": 1, "4": 2}


def test_source_order_reverses_all_accepted_transition_checkpoints(
    tmp_path: Path,
) -> None:
    stages = [tmp_path / "02_acyl/02_unbiased_release"]
    stages += [tmp_path / f"03_reverse_path/window_{i:03d}" for i in range(1, 37)]
    stages += [tmp_path / "03_reverse_path/unbiased_ES_release"]
    for stage in stages:
        stage.mkdir(parents=True)
        (stage / "process.json").write_text(json.dumps({"completed": True}))
        (stage / "scientific_acceptance.json").write_text(
            json.dumps({"force_converged": True})
        )
        for name in ("accepted_00000.ncrst", "accepted_00025.ncrst", "optimized.ncrst"):
            (stage / name).touch()
    paths = stage_trace_paths(tmp_path)
    assert paths[0] == stages[-1] / "optimized.ncrst"
    assert paths[-1] == stages[0] / "optimized.ncrst"
    transition = [p.name for p in paths if p.parent.name == "window_025"]
    assert transition == [
        "optimized.ncrst",
        "accepted_00025.ncrst",
        "accepted_00000.ncrst",
    ]
    (stages[4] / "scientific_acceptance.json").write_text(
        json.dumps({"force_converged": False})
    )
    with pytest.raises(ValueError):
        stage_trace_paths(tmp_path)
