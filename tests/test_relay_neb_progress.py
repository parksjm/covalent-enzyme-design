"""Verify coherent relay NEB publication and terminal scientific acceptance."""

from pathlib import Path
import json
import shlex
import shutil
import subprocess
from types import SimpleNamespace
from typing import Any

import numpy as np
from parmed.amber import Rst7
import pytest
from test_neb_progress import native_report

from nyl12_covalent import relay_neb_progress as helper
from nyl12_covalent.neb_progress import parse_last_complete_neb_report


def test_publish_preserves_native_proof_and_small_image_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    band = tmp_path / "band"
    band.mkdir()
    output = tmp_path / "04_neb/progress"
    output.mkdir(parents=True)
    monkeypatch.setattr(
        helper,
        "image_geometry",
        lambda *_: {
            "rWaterH1_OG1": 1.0,
            "rWaterH1_Nleave": 2.0,
            "rO_PheN": 3.0,
            "rO_AsnND2": 4.0,
            "BD_angle_deg": 110.0,
        },
    )
    for i in (1, 2):
        restart = Rst7(natom=1)
        restart.coordinates = np.array([[float(i), 0, 0]])
        restart.write(str(band / f"optimized.{i:04d}.ncrst"), netcdf=True)
    payload = native_report([-10.0, -9.9]) + b"Testing convergence in cycle 3\n"
    report = parse_last_complete_neb_report(payload, 2)
    (band / "source_nebpath.xyz").write_bytes(b"original native xyz\n")
    (band / "source_output_snapshot.mdout").write_bytes(payload)
    (band / "band_snapshot.json").write_text(json.dumps({"source_xyz_sha256": "abc"}))
    monkeypatch.setattr(
        helper,
        "load_case",
        lambda _: SimpleNamespace(
            config={},
            atom_map={
                key: 1
                for key in (
                    "water_h1",
                    "proton",
                    "thr_og",
                    "leaving_n",
                    "carbon",
                    "oxygen",
                    "phe_n",
                    "asn_nd2",
                )
            },
        ),
    )
    monkeypatch.setattr(
        helper,
        "state_record",
        lambda model, xyz: {
            key: float(xyz[0, 0])
            for key in (
                "rOG1C",
                "rCN",
                "rCO",
                "rNwater",
                "rwaterLeave",
                "CN_N",
                "CN_OG1",
                "CN_water",
                "CN_leave",
            )
        },
    )

    def movie(folder: Path, frames: list, case: dict) -> None:
        assert len(frames) == 2
        bundle = folder / "vmd"
        bundle.mkdir()
        (bundle / "manifest.json").write_text("{}")
        (bundle / "view_path.vmd").write_text("movie")
        (bundle / "MANIFEST.sha256").write_text(
            "old  manifest.json\nold  view_path.vmd\n"
        )

    monkeypatch.setattr(helper, "export_movie", movie)
    result = helper.publish_band(
        {"case_path": "case"},
        band,
        report,
        output,
        {"id": 1},
        False,
        tmp_path / "sources",
    )
    latest = output / "latest"
    assert latest.is_symlink()
    rows = json.loads((latest / "images.json").read_text())
    assert len(rows) == 2 and rows[1]["state"]["CN_leave"] == 2
    assert rows[1]["physical_energy_kcal_mol"] < 0
    assert "not a verified TS" in (latest / "energy.svg").read_text()
    for key in ("rWaterH1_OG1", "rWaterH1_Nleave", "rO_PheN", "rO_AsnND2"):
        assert key in (latest / "reaction_distances.svg").read_text()
        assert key in rows[0]["state"]
    assert "BD_angle_deg" in (latest / "attack_angle.svg").read_text()
    assert result["geometry_atom_indices_0based"]["water_h1"] == 0
    portable_names = [
        line.split("  ", 1)[1]
        for line in (latest / "MANIFEST.sha256").read_text().splitlines()
    ]
    assert "native_report.mdout" not in portable_names
    assert "source_nebpath.xyz.gz" not in portable_names
    for name, digest in result["local_native_evidence_sha256"].items():
        assert helper.digest(latest / name) == digest
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync required for actual portable-manifest transfer regression")
    text = (helper.ROOT / "scripts/pull_nyl12_inspection.sh").read_text()
    block = text.split("# Relay acylation campaign:", 1)[1].split(
        "# Final embedded-method", 1
    )[0]
    options = shlex.split(
        block.split('rsync "${vmd_rsync_args[@]}"', 1)[1]
        .split('"$remote:$relay_campaign_source/"', 1)[0]
        .replace("\\\n", "")
    )
    target = tmp_path / "pulled"
    subprocess.run(
        [rsync, "-rt", *options, str(tmp_path) + "/", str(target) + "/"], check=True
    )
    copied = target / "04_neb/progress/latest"
    for manifest in (copied / "MANIFEST.sha256", copied / "vmd/MANIFEST.sha256"):
        assert manifest.is_file()
        for line in manifest.read_text().splitlines():
            expected, name = line.split("  ", 1)
            assert (manifest.parent / name).is_file(), name
            assert helper.digest(manifest.parent / name) == expected, name
    assert not (copied / "source_nebpath.xyz.gz").exists()
    assert not (copied / "native_report.mdout").exists()
    assert "numbered native restart" in result["coordinate_source"]
    assert (latest / "source_nebpath.xyz.gz").is_file()
    for line in (latest / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        assert helper.digest(latest / name) == expected
    assert (
        "path candidate"
        in json.loads((latest / "vmd/manifest.json").read_text())["state"]
    )


def test_terminal_acceptance_keeps_chemistry_and_curvature_separate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Terminal publication must retain native convergence/SCC and provisional chemistry."""
    from scripts import run_relay_acyl_campaign

    band = tmp_path / "band"
    band.mkdir()
    stage = tmp_path / "04_neb"
    output = stage / "progress"
    output.mkdir(parents=True)
    (stage / "job").mkdir()
    (stage / "job/status.json").write_text(
        json.dumps({"ended_utc": "done", "exit_status": 0})
    )
    monkeypatch.setattr(
        helper,
        "image_geometry",
        lambda *_: {
            "rWaterH1_OG1": 1.0,
            "rWaterH1_Nleave": 2.0,
            "rO_PheN": 3.0,
            "rO_AsnND2": 4.0,
            "BD_angle_deg": 110.0,
        },
    )
    for i in (1, 2):
        restart = Rst7(natom=1)
        restart.coordinates = np.array([[float(i), 0, 0]])
        restart.write(str(band / f"optimized.{i:04d}.ncrst"), netcdf=True)
    payload = (
        native_report([-10.0, -9.9])
        + b"Testing convergence in cycle 3\nOptimisation converged\n"
    )
    report = parse_last_complete_neb_report(payload, 2)
    (band / "source_nebpath.xyz").write_bytes(b"nativexyz\n")
    (band / "source_output_snapshot.mdout").write_bytes(payload)
    (band / "band_snapshot.json").write_text(json.dumps({"source_xyz_sha256": "abc"}))
    monkeypatch.setattr(
        helper,
        "load_case",
        lambda _: SimpleNamespace(
            config={},
            atom_map={
                key: 1
                for key in (
                    "water_h1",
                    "proton",
                    "thr_og",
                    "leaving_n",
                    "carbon",
                    "oxygen",
                    "phe_n",
                    "asn_nd2",
                )
            },
        ),
    )
    monkeypatch.setattr(
        helper,
        "state_record",
        lambda model, xyz: {
            key: 1.0
            for key in (
                "rOG1C",
                "rCN",
                "rCO",
                "rNwater",
                "rwaterLeave",
                "CN_N",
                "CN_OG1",
                "CN_water",
                "CN_leave",
            )
        },
    )
    monkeypatch.setattr(
        run_relay_acyl_campaign,
        "classify",
        lambda *_: {"state": {}, "product_geometry_retained": True},
    )

    def movie(folder: Path, *_: Any) -> None:
        bundle = folder / "vmd"
        bundle.mkdir()
        (bundle / "manifest.json").write_text("{}")
        (bundle / "MANIFEST.sha256").write_text("old  manifest.json\n")

    monkeypatch.setattr(helper, "export_movie", movie)
    helper.publish_band(
        {"case_path": "case"}, band, report, output, {}, True, tmp_path / "sources"
    )
    record = json.loads((stage / "scientific_acceptance.json").read_text())
    assert record["native_optimisation_converged"] and record["scc_error_count"] == 0
    assert record["provisional_endpoint_chemistry"]["product_geometry_retained"]
    assert (
        record["curvature_verified"] is False and record["automatic_ts_launch"] is False
    )
    assert record["profile_maximum_is_ts_barrier"] is False
    cfg = stage / "config.json"
    cfg.write_text("{}")
    native = stage / "calculation"
    native.mkdir()
    (native / "minimized.mdout").write_text("SCC failed to converge\n")
    assert (
        helper.record_terminal_failure(cfg, ValueError("invalid final band"))
        == stage / "scientific_acceptance.json"
    )
    failed = json.loads((stage / "scientific_acceptance.json").read_text())
    assert (
        failed["scc_error_count"] == 1
        and failed["provisional_endpoint_chemistry"] is None
    )


def test_image_geometry_uses_actual_transfer_h_and_explicit_case_mapping() -> None:
    """Two H atoms at different coordinates must never be silently interchanged."""
    atom_map = {
        "water_h1": 8,
        "proton": 2,
        "thr_og": 5,
        "leaving_n": 3,
        "carbon": 7,
        "oxygen": 1,
        "phe_n": 6,
        "asn_nd2": 4,
    }
    xyz = np.zeros((8, 3))
    xyz[7] = [1.0, 0.0, 0.0]  # Actual transferring water H1.
    xyz[1] = [9.0, 9.0, 9.0]  # Original HG1, deliberately far away.
    xyz[4] = [0.0, 0.0, 0.0]
    xyz[2] = [3.0, 0.0, 0.0]
    xyz[6] = [0.0, 1.0, 0.0]
    xyz[0] = [1.0, 1.0, 0.0]
    xyz[5] = [1.0, 4.0, 0.0]
    xyz[3] = [1.0, 1.0, 4.0]
    values = helper.image_geometry(atom_map, xyz)
    assert values == pytest.approx(
        {
            "rWaterH1_OG1": 1.0,
            "rWaterH1_Nleave": 2.0,
            "rO_PheN": 3.0,
            "rO_AsnND2": 4.0,
            "BD_angle_deg": 90.0,
        }
    )
    xyz[7] = [2.0, 0.0, 0.0]
    moved = helper.image_geometry(atom_map, xyz)
    assert moved["rWaterH1_OG1"] == pytest.approx(2.0)
    assert moved["rWaterH1_Nleave"] == pytest.approx(1.0)
    assert moved["rO_PheN"] == values["rO_PheN"]
