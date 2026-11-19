"""Frozen-band spacing controls and variable native-image mapping tests."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import neb_spacing as spacing
from nyl12_covalent.relay_neb_native import native_initial_band
from test_relay_neb_preparation import controls


def source_band() -> np.ndarray:
    """Return24 curved active images with large unrelated MM displacements."""
    xyz = np.zeros((24, 3, 3))
    t = np.linspace(0, 1, 24)
    xyz[:, 0, 0] = t
    xyz[:, 0, 1] = t * t
    xyz[:, 1, 2] = 1000 * t
    return xyz


def test_uniform_and_adaptive_preserve_endpoint_arrays_and_focus_segment() -> None:
    source = source_band()
    original = source.copy()
    uniform, u = spacing.spacing_controls(source, np.array([0]))
    adaptive, a = spacing.spacing_controls(source, np.array([0]), focus_weight=4.0)
    for band in (uniform, adaptive):
        np.testing.assert_array_equal(band[[0, -1]], source[[0, -1]])
    np.testing.assert_array_equal(source, original)
    assert u["active_arc_A"] > 100 * u["qm_arc_A"]
    assert sum(
        row["source_images_1based"] == [9, 10] for row in a["control_brackets"]
    ) > sum(row["source_images_1based"] == [9, 10] for row in u["control_brackets"])
    assert not np.array_equal(uniform[1:-1], adaptive[1:-1])


@pytest.mark.parametrize("bad", ["shape", "metric", "density", "zero"])
def test_sampling_rejects_invalid_inputs(bad: str) -> None:
    xyz, qm, density = source_band(), np.array([0]), 1.0
    if bad == "shape":
        xyz = xyz[:-1]
    elif bad == "metric":
        qm = np.array([10])
    elif bad == "density":
        density = 0
    elif bad == "zero":
        xyz[:] = 0
    with pytest.raises(ValueError):
        spacing.spacing_controls(xyz, qm, focus_weight=density)


def test_generic24_predictor_matches_validated_native_predictor() -> None:
    xyz = controls()
    qm = np.arange(5)
    new, report = spacing.predict_spacing_band(xyz, qm, 24)
    old, old_report = native_initial_band(xyz, qm)
    np.testing.assert_array_equal(new, old)
    assert (
        report["slots_1based"] == old_report["source_controls_to_regular_slots_1based"]
    )
    larger, evidence = spacing.predict_spacing_band(xyz, qm, 48)
    assert larger.shape == (48, 5, 3) and evidence["slots_1based"][-1] == 48
    np.testing.assert_array_equal(larger[[0, -1]], new[[0, -1]])
    assert evidence["endpoint_within_1e_minus6_A"] is False


@pytest.mark.parametrize("bad", ["shape", "metric", "count", "zero"])
def test_native_prediction_rejects_invalid_inputs(bad: str) -> None:
    xyz, qm, count = controls(), np.arange(5), 24
    if bad == "shape":
        xyz = xyz[:-1]
    elif bad == "metric":
        qm = np.array([-1])
    elif bad == "count":
        count = 25
    elif bad == "zero":
        xyz[:] = xyz[0]
    with pytest.raises(ValueError):
        spacing.predict_spacing_band(xyz, qm, count)


@pytest.mark.parametrize("bad", [None, "hash", "topology", "restraint", "active"])
def test_source_loader_verifies_immutable_payload_and_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bad: str | None
) -> None:
    config = tmp_path / "config.json"
    reference = tmp_path / "reference.ncrst"
    topology = tmp_path / "top.prmtop"
    reference.write_text("reference")
    topology.write_text("topology")

    def sha(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    cfg = dict(
        reactant=str(reference),
        topology=str(topology),
        active_mask="@1-4662",
        case_path=str(tmp_path / "case.json"),
    )
    config.write_text(json.dumps(cfg))
    payload = b"fixed XYZ payload\n"
    (tmp_path / "source_nebpath.xyz.gz").write_bytes(gzip.compress(payload))
    meta = dict(
        reference={"path": str(reference), "sha256": sha(reference)},
        topology={"path": str(topology), "sha256": sha(topology)},
        restraints=False,
        active_mask="@1-4662",
        source_xyz_sha256=hashlib.sha256(payload).hexdigest(),
        active_indices_0based=list(range(4662)),
        report={"cycle": 268},
    )
    if bad == "hash":
        meta["source_xyz_sha256"] = "incorrect"
    elif bad == "topology":
        meta["topology"]["sha256"] = "incorrect"
    elif bad == "restraint":
        meta["restraints"] = True
    elif bad == "active":
        meta["active_indices_0based"] = [1]
    (tmp_path / "band_snapshot.json").write_text(json.dumps(meta))
    bonds = [
        SimpleNamespace(atom1=SimpleNamespace(idx=0), atom2=SimpleNamespace(idx=200))
        for _ in range(8)
    ]
    model = SimpleNamespace(
        qm=np.arange(102),
        hydrogens=np.arange(54),
        active=np.ones(4662, dtype=bool),
        topology=SimpleNamespace(bonds=bonds),
        config={"qm_charge": -1},
    )
    monkeypatch.setattr(spacing, "load_case", lambda _: model)
    monkeypatch.setattr(spacing.Rst7, "open", lambda _: SimpleNamespace())
    monkeypatch.setattr(
        spacing,
        "reconstruct_xyz_band",
        lambda *args: ([np.zeros((4662, 3)) for _ in range(24)], list(range(4662))),
    )
    if bad:
        with pytest.raises(ValueError):
            spacing.load_spacing_source(config, tmp_path)
    else:
        _, _, frames, evidence = spacing.load_spacing_source(config, tmp_path)
        assert frames.shape == (24, 4662, 3)
        assert evidence["numbered_native_restarts_used"] is False
        assert evidence["cycle"] == 268
