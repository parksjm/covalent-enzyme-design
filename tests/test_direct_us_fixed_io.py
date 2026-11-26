"""Synthetic validation tests for the direct-US fixed-analysis adapter."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import direct_us_fixed_io as direct_io
from nyl12_covalent.neb_umbrella_bias import HarmonicWindow, amber_restraints
from nyl12_covalent.relay_umbrella_analysis import KB
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_direct_reference_hashes_and_restores_seed_image_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = []
    expected = {}
    for image in reversed(range(1, 25)):
        path = tmp_path / f"image_{image:03}.rst7"
        path.write_text(str(image))
        sha = _digest(path)
        entries.append(
            dict(seed_image=image, source_restart=str(path), source_restart_sha256=sha)
        )
        expected[str(path)] = sha

    class FakeRst7:
        @staticmethod
        def open(filename: str) -> SimpleNamespace:
            return SimpleNamespace(
                coordinates=np.array([int(Path(filename).stem[-3:]), 0.0, 0.0]),
                box=None,
            )

    monkeypatch.setattr(direct_io, "Rst7", FakeRst7)
    monkeypatch.setattr(
        direct_io,
        "raw_features",
        lambda coordinates, box: dict(
            q_acyl=float(coordinates[0, 0]),
            q_act=float(coordinates[0, 0] + 100),
            q_del=float(coordinates[0, 0] + 200),
        ),
    )
    points, hashes = direct_io.direct_reference(dict(windows=entries))
    np.testing.assert_array_equal(points[:, 0], np.arange(1, 25))
    np.testing.assert_array_equal(points[:, 1], np.arange(101, 125))
    np.testing.assert_array_equal(points[:, 2], np.arange(201, 225))
    assert hashes == expected


def test_direct_reference_rejects_wrong_image_set_and_hash_mismatch(
    tmp_path: Path,
) -> None:
    rows = [
        dict(seed_image=i, source_restart=f"image_{i}", source_restart_sha256="sha")
        for i in range(1, 24)
    ]
    with pytest.raises(ValueError, match="exactly the original 24"):
        direct_io.direct_reference(dict(windows=rows))

    path = tmp_path / "changed.rst7"
    path.write_text("changed")
    # Reuse the fixture path while keeping all 24 declared image IDs distinct.
    rows = [
        dict(
            seed_image=i,
            source_restart=str(path),
            source_restart_sha256="bad",
        )
        for i in range(1, 25)
    ]
    with pytest.raises(ValueError, match="restart changed"):
        direct_io.direct_reference(dict(windows=rows))


def test_direct_reference_rejects_conflicting_duplicate_image() -> None:
    rows = [
        dict(seed_image=i, source_restart=f"image_{i}", source_restart_sha256="sha")
        for i in range(1, 25)
    ]
    rows.append(dict(seed_image=1, source_restart="other", source_restart_sha256="sha"))
    with pytest.raises(ValueError, match="conflicting restart identities"):
        direct_io.direct_reference(dict(windows=rows))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _dataset_fixture(root: Path) -> tuple[Path, Path, Path, dict]:
    campaign = root / "campaign"
    augmented = root / "augmented"
    campaign.mkdir(parents=True)
    (campaign / "inputs").mkdir()
    topology = campaign / "inputs/system.prmtop"
    topology.write_text("synthetic topology identity")
    windows = []
    frame = dict(
        q_acyl=0.1,
        q_act=0.2,
        q_del=-0.05,
        **{name: float(i) for i, name in enumerate(direct_io.DIAGNOSTICS)},
        time_ps="5.0",
    )
    for index in range(63):
        label = f"window_{index + 1:03}"
        bias = dict(
            q_acyl_center=-0.5 + index / 62,
            q_h_center=0.1,
            k_acyl=100.0,
            k_h=100.0,
            outer_extent=10.0,
        )
        windows.append(
            dict(
                label=label,
                bias=bias,
                seed_image=(index % 24) + 1,
                source_restart=f"restart_{(index % 24) + 1}",
                source_restart_sha256=f"restart_sha_{(index % 24) + 1}",
                stages=[dict(kind="production", label="prod")],
            )
        )
        native = campaign / label / "prod"
        for filename in (
            "process.json",
            "acceptance.json",
            "observables/observables.json",
            "observables/observables.csv",
        ):
            path = native / filename
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("synthetic compact source")
        expected_sources = {
            str(native / name): _digest(native / name)
            for name in (
                "process.json",
                "acceptance.json",
                "observables/observables.json",
                "observables/observables.csv",
            )
        }
        (native / "umbrella.RST").write_text(
            amber_restraints(
                HarmonicWindow(**bias),
                dict(carbon=1, leaving_n=2, thr_og=3, proton=4),
            )
        )
        row = dict(frame)
        augmented_stage = augmented / label / "prod"
        csv_path = augmented_stage / "observables.csv"
        _write_csv(csv_path, [row])
        receipt = dict(csv_sha256=_digest(csv_path), source_sha256=expected_sources)
        (augmented_stage / "augmentation.json").write_text(json.dumps(receipt))

    native_biases = [
        RelayBiasWindow(
            feature_names=("q_acyl", "q_h"),
            centers_A=(window["bias"]["q_acyl_center"], window["bias"]["q_h_center"]),
            force_constants=(100.0, 100.0),
        )
        for window in windows
    ]
    features = dict(
        q_acyl=np.array([float(frame["q_acyl"])]),
        q_act=np.array([float(frame["q_act"])]),
        q_del=np.array([float(frame["q_del"])]),
        q_h=np.array([float(frame["q_act"]) + float(frame["q_del"])]),
    )
    cross = cross_bias_energies(features, native_biases)
    for window in windows:
        native = campaign / window["label"] / "prod"
        np.savez(
            native / "samples_and_cross_bias.npz",
            cv_values=np.array([[frame["q_acyl"], features["q_h"][0]]]),
            bias_kcal_mol=cross,
            u_kn_bias_only=cross / (KB * 300.0),
        )
    prepared_path = campaign / "prepared.json"
    prepared_path.write_text(
        json.dumps(
            dict(
                temperature_K=300.0,
                windows=windows,
                roles_1based=dict(carbon=1, leaving_n=2, thr_og=3, proton=4),
            )
        )
    )
    (campaign / "prepared.sha256").write_text(_digest(prepared_path))
    relay_request = root / "relay_request.json"
    relay_request.write_text(
        json.dumps(
            dict(
                temperature_K=300.0,
                s_bin_edges=[-0.1, 0.0, 0.1],
                es_search_bounds=[0.0, 0.1],
                acyl_search_bounds=[0.9, 1.0],
            )
        )
    )
    return campaign, augmented, relay_request, frame


def test_load_direct_dataset_checks_native_bias_parity_and_preserves_all_windows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign, augmented, request, _ = _dataset_fixture(tmp_path)
    reference = np.column_stack((np.linspace(-1, 1, 24), np.zeros((24, 2))))
    monkeypatch.setattr(
        direct_io, "direct_reference", lambda prepared: (reference, {"ref.rst7": "abc"})
    )
    config, blocks, sources = direct_io.load_direct_dataset(
        campaign, augmented, request
    )
    assert len(config["windows"]) == len(blocks) == 63
    np.testing.assert_array_equal(config["path_3cv_A"], reference)
    assert config["fixed_reference_s"] == 0.0125
    assert config["bin_extension"]["frames_censored"] == 0
    assert config["bin_extension"]["added_right_bins"] > 0
    assert config["s_bin_edges"][0] == pytest.approx(-0.1)
    assert blocks[0]["values"].shape == (1, 13)
    assert blocks[0]["values"][0, 1] + blocks[0]["values"][0, 2] == pytest.approx(0.15)
    assert blocks[0]["native_bias_kcal_mol"].shape == (1,)
    assert sources["ref.rst7"] == "abc"


@pytest.mark.parametrize(
    ("corruption", "message"),
    [
        ("temperature", "Physical temperatures differ"),
        ("csv_hash", "Augmented direct CSV changed"),
        ("source_hash", "Compact source changed"),
        ("restraint", "native restraint differs"),
        ("cv", "CV archive differs"),
        ("bias", "cross-bias parity failed"),
        ("reduced", "units or temperature differ"),
        ("window_count", "Expected 63 unique"),
    ],
)
def test_load_direct_dataset_rejects_compact_or_physical_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    corruption: str,
    message: str,
) -> None:
    campaign, augmented, request, _ = _dataset_fixture(tmp_path)
    prepared_path = campaign / "prepared.json"
    prepared = json.loads(prepared_path.read_text())
    relay = json.loads(request.read_text())
    native = campaign / "window_001/prod"
    derived = augmented / "window_001/prod"
    if corruption == "temperature":
        relay["temperature_K"] = 301.0
        request.write_text(json.dumps(relay))
    elif corruption == "csv_hash":
        with (derived / "observables.csv").open("a") as stream:
            stream.write("changed\n")
    elif corruption == "source_hash":
        (native / "process.json").write_text("changed")
    elif corruption == "restraint":
        (native / "umbrella.RST").write_text("changed restraint")
    elif corruption == "cv":
        with np.load(native / "samples_and_cross_bias.npz") as old:
            contents = {key: old[key] for key in old.files}
        contents["cv_values"][0, 1] += 1.0
        np.savez(native / "samples_and_cross_bias.npz", **contents)
    elif corruption == "bias":
        with np.load(native / "samples_and_cross_bias.npz") as old:
            contents = {key: old[key] for key in old.files}
        contents["bias_kcal_mol"][0, 0] += 1.0
        np.savez(native / "samples_and_cross_bias.npz", **contents)
    elif corruption == "reduced":
        with np.load(native / "samples_and_cross_bias.npz") as old:
            contents = {key: old[key] for key in old.files}
        contents["u_kn_bias_only"][0, 0] += 1.0
        np.savez(native / "samples_and_cross_bias.npz", **contents)
    elif corruption == "window_count":
        prepared["windows"].pop()
        prepared_path.write_text(json.dumps(prepared))
        (campaign / "prepared.sha256").write_text(_digest(prepared_path))
    monkeypatch.setattr(
        direct_io,
        "direct_reference",
        lambda value: (np.zeros((24, 3)), {}),
    )
    with pytest.raises(ValueError, match=message):
        direct_io.load_direct_dataset(campaign, augmented, request)
