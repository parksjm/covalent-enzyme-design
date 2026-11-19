"""Accepted-frame, atomic progress and bounded laptop-transfer contracts."""

import json
from pathlib import Path
import shlex
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import relay_progress as helper
from scripts import export_relay_campaign_progress as monitor


def record(iteration: int = 0) -> dict:
    """Return a complete synthetic accepted state with no scientific engine call."""
    return {
        "iteration": iteration,
        "physical_energy_kcal_mol": -100.0,
        "total_energy_kcal_mol": -98.0,
        "bias_energy_kcal_mol": 2.0,
        "state": {
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
    }


def stage(campaign: Path, name: str, iterations: list[int]) -> Path:
    """Create synthetic accepted-state records and periodic snapshot placeholders."""
    folder = campaign / name
    folder.mkdir(parents=True)
    (folder / "accepted_iterations.jsonl").write_text(
        "".join(json.dumps(record(i)) + "\n" for i in iterations)
    )
    for i in iterations:
        (folder / f"accepted_{i:05d}.ncrst").write_text(str(i))
    return folder


def test_accepted_rows_frames_partial_and_invalid(tmp_path: Path) -> None:
    folder = stage(tmp_path, "01_test", [0, 1, 25, 27])
    with (folder / "accepted_iterations.jsonl").open("a") as stream:
        stream.write('{"incomplete":')
    rows = helper.accepted_rows(tmp_path)
    assert [r["iteration"] for r in rows] == [0, 1, 25, 27]
    assert [p.name for p, _ in helper.snapshot_frames(tmp_path, rows)] == [
        "accepted_00000.ncrst",
        "accepted_00025.ncrst",
    ]
    (folder / "process.json").write_text("{}")
    (folder / "optimized.ncrst").write_text("final")
    assert helper.snapshot_frames(tmp_path, rows)[-1][1]["iteration"] == 27
    copy = tmp_path / "progress/update_000001"
    copy.mkdir(parents=True)
    (copy / "accepted_iterations.jsonl").write_text("invalid\n")
    assert len(helper.accepted_rows(tmp_path)) == 4
    for entries in [
        [record(0), record(0)],
        [{**record(0), "physical_energy_kcal_mol": np.nan}],
    ]:
        (folder / "accepted_iterations.jsonl").write_text(
            "".join(json.dumps(r) + "\n" for r in entries)
        )
        with pytest.raises(ValueError):
            helper.accepted_rows(tmp_path)


def test_plots_and_tables(tmp_path: Path) -> None:
    rows = [{**record(i), "stage": f"stage_{i}", "nested": [1]} for i in (0, 1)]
    helper.write_table(tmp_path / "records.csv", rows)
    helper.plot_progress(tmp_path, rows)
    assert "stage_0:0" in (tmp_path / "energy.svg").read_text()
    assert "CN_leave" in (tmp_path / "coordination.svg").read_text()
    assert "nested" not in (tmp_path / "records.csv").read_text()


def test_atomic_refresh_stability_and_terminal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "case.json").write_text("{}")
    stage(campaign, "01_test", [0])
    output = campaign / "progress"
    seen = []
    monkeypatch.setattr(
        monitor, "export_movie", lambda out, frames, case: seen.append(frames)
    )
    assert monitor.refresh(campaign, output) is False
    assert seen[-1] == []
    first = (output / "latest").resolve()
    assert monitor.refresh(campaign, output) is False
    assert len(seen[-1]) == 1
    second = (output / "latest").resolve()
    assert first != second and first.is_dir()
    assert monitor.refresh(campaign, output) is False
    assert (output / "latest").resolve() == second
    (campaign / "campaign_status.json").write_text('{"finished": true}')
    assert monitor.refresh(campaign, output) is True
    latest = (output / "latest").resolve()
    assert json.loads((latest / "progress.json").read_text())["terminal"] is True
    assert (latest / "source/relay_progress_vmd.py").is_file()
    checksum = (latest / "MANIFEST.sha256").read_text()
    assert "  accepted_iterations.json\n" not in checksum
    assert (latest / "accepted_iterations.json").is_file()
    report = json.loads((latest / "progress.json").read_text())
    assert report["remote_evidence_sha256"][
        "accepted_iterations.json"
    ] == monitor.digest(latest / "accepted_iterations.json")
    assert not list(output.glob(".building_*"))


def test_explicit_terminal_status(tmp_path: Path) -> None:
    assert not monitor.finished(tmp_path)
    for status, expected in [("running", False), ("failed", True), ("completed", True)]:
        (tmp_path / "campaign_status.json").write_text(json.dumps({"status": status}))
        assert monitor.finished(tmp_path) is expected


def test_relay_pull_filter(tmp_path: Path) -> None:
    rsync = shutil.which("rsync")
    if rsync is None:
        pytest.skip("rsync unavailable")
    text = (monitor.ROOT / "scripts/pull_nyl12_inspection.sh").read_text()
    block = text.split("# Relay acylation campaign:", 1)[1].split(
        "# Final embedded-method", 1
    )[0]
    arguments = block.split('rsync "${vmd_rsync_args[@]}"', 1)[1].split('"$remote:', 1)[
        0
    ]
    options = shlex.split(arguments.replace("\\\n", ""))
    source, target = tmp_path / "source", tmp_path / "target"
    keep = [
        "case.json",
        "00_screen/report.json",
        "01_test/accepted_iterations.jsonl",
        "01_test/source_snapshot/code.py",
        "progress/update_000001/vmd/path.dcd",
        "progress/update_000001/vmd/path.parm7",
        "progress/update_000001/MANIFEST.sha256",
    ]
    reject = [
        "01_test/optimized.ncrst",
        "01_test/evaluations/report.json",
        "01_test/source_snapshot/native.so",
        "03_md/production.nc",
        "progress/update_000001/accepted_iterations.json",
        "progress/update_000001/accepted_iterations.jsonl",
    ]
    for name in keep + reject:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    (source / "progress/latest").symlink_to("update_000001")
    subprocess.run(
        [rsync, "-rt", *options, str(source) + "/", str(target) + "/"], check=True
    )
    assert (target / "progress/latest/vmd/path.dcd").is_file()
    assert not (target / "progress/latest").is_symlink()
    assert not (target / "progress/update_000001").exists()
    assert all((target / name).is_file() for name in (keep[0], keep[1], keep[3]))
    assert not (target / keep[2]).exists()  # Raw accepted ledgers are excluded.
    assert not any((target / name).exists() for name in reject)
    assert not (target / "progress/latest/accepted_iterations.json").exists()
    assert not (target / "progress/latest/accepted_iterations.jsonl").exists()


def test_mapping_checks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    qm = [*range(99), 20316, 20317, 20318]
    atoms = [SimpleNamespace(atomic_number=6) for _ in range(20319)]
    for i in [*range(52), 20317, 20318]:
        atoms[i].atomic_number = 1
    residues = [SimpleNamespace(atoms=[]) for _ in range(1576)]
    residues[1575].atoms = [SimpleNamespace(idx=i) for i in (20316, 20317, 20318)]
    topology = SimpleNamespace(atoms=atoms, residues=residues)
    monkeypatch.setattr(helper.parmed, "load_file", lambda path: topology)
    monkeypatch.setattr(
        helper, "AmberMask", lambda top, mask: SimpleNamespace(Selected=lambda: qm)
    )
    monkeypatch.setattr(helper, "COMPARISON", tmp_path)
    prep = (
        tmp_path
        / "06_qmmm_optimization/replica_1_tier_A_seed/water_qm_reoptimization/preparation.json"
    )
    prep.parent.mkdir(parents=True)
    prep.write_text('{"qm_region":{"mask":"explicit"}}')
    case = {"topology": "fixture", "qm_atom_ids": [i + 1 for i in qm]}
    assert helper.verify_mapping(case) is topology
    with pytest.raises(ValueError, match="QM102"):
        helper.verify_mapping(
            {**case, "qm_atom_ids": list(reversed(case["qm_atom_ids"]))}
        )
    residues[1575].atoms = [SimpleNamespace(idx=0)]
    with pytest.raises(ValueError, match="water1576"):
        helper.verify_mapping(case)


def test_movie_reuses_unique_geometry_and_original_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    helper.export_movie(tmp_path, [], {})
    paths = [tmp_path / f"{i}.ncrst" for i in range(2)]
    for path in paths:
        path.write_text("same accepted geometry")
    seen = []

    def export(
        topology: Path, links: list, energy: Path, distances: Path, output: Path
    ) -> None:
        seen.extend(links)
        assert energy.is_file() and distances.is_file()
        output.mkdir()
        (output / "view_path.vmd").write_text("old view")
        (output / "manifest.json").write_text("{}")
        (output / "MANIFEST.sha256").write_text(
            "old  view_path.vmd\nold  manifest.json\n"
        )

    fake = SimpleNamespace(sha256=monitor.digest, export=export)
    monkeypatch.setattr(helper, "legacy_exporter", lambda: fake)
    monkeypatch.setattr(helper, "verify_mapping", lambda case: "verified-topology")
    monkeypatch.setattr(
        helper, "vmd_script", lambda manifest, case, topology: "new view"
    )
    frames = [(p, {**record(i), "stage": f"{i:02d}_test"}) for i, p in enumerate(paths)]
    helper.export_movie(tmp_path, frames, {"topology": "fixture"})
    assert len(seen) == 1
    manifest = json.loads((tmp_path / "vmd/manifest.json").read_text())
    assert [
        r["movie_frame_0based"] for r in manifest["accepted_stage_frame_records"]
    ] == [0, 0]
    assert [r["path"] for r in manifest["accepted_restarts"]] == [str(p) for p in paths]
    assert (
        monitor.digest(tmp_path / "vmd/view_path.vmd")
        in (tmp_path / "vmd/MANIFEST.sha256").read_text()
    )


def test_main_and_journal(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    commands = []
    monkeypatch.setattr(
        monitor.subprocess, "run", lambda args, check: commands.append(args)
    )
    monitor.journal("intent", "fixture read-only check")
    assert "--phase" in commands[0] and "intent" in commands[0]
    monkeypatch.setattr(monitor, "journal", lambda *args: None)
    monkeypatch.setattr(monitor, "refresh", lambda *args: True)
    monkeypatch.setattr(monitor.sys, "argv", ["export", "--campaign", str(tmp_path)])
    monitor.main()
    monkeypatch.setattr(monitor.sys, "argv", ["export", "--interval", "0"])
    with pytest.raises(SystemExit):
        monitor.main()
    monkeypatch.setattr(monitor.sys, "argv", ["export", "--watch"])
    monkeypatch.setattr(monitor, "refresh", lambda *args: False)
    ticks = iter([0.0, 48.0 * 3600 + 1.0])
    monkeypatch.setattr(monitor.time, "monotonic", lambda: next(ticks))
    with pytest.raises(TimeoutError):
        monitor.main()

    def failing_refresh(*args: object) -> bool:
        raise ValueError("explicit failure")

    monkeypatch.setattr(monitor, "refresh", failing_refresh)
    monkeypatch.setattr(monitor.time, "monotonic", lambda: 0.0)
    with pytest.raises(ValueError, match="explicit failure"):
        monitor.main()


def test_view_only_update_preserves_verified_base(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from scripts import refresh_relay_campaign_view as view

    campaign = tmp_path / "campaign"
    base = campaign / "progress/update_000001"
    bundle = base / "vmd"
    bundle.mkdir(parents=True)
    case = {
        "atom_map": {
            key: i + 1
            for i, key in enumerate(
                (
                    "carbon",
                    "leaving_n",
                    "thr_og",
                    "water_h1",
                    "oxygen",
                    "phe_n",
                    "asn_nd2",
                )
            )
        }
    }
    (campaign / "case.json").write_text(json.dumps(case))
    names = [
        "path.dcd",
        "path.parm7",
        "path.csv",
        "view_path.vmd",
        "energy.svg",
        "distances.svg",
        "angles.svg",
        "export.cpptraj",
        "cpptraj.log",
    ]
    for name in names:
        (bundle / name).write_text("old " + name)
    (bundle / "manifest.json").write_text("{}")
    names.append("manifest.json")
    (bundle / "MANIFEST.sha256").write_text(
        "".join(f"{monitor.digest(bundle / name)}  {name}\n" for name in names)
    )
    (base / "accepted_iterations.json").write_text('[{"authoritative":"raw"}]')
    (base / "progress.json").write_text(
        json.dumps(
            {
                "identity": {
                    "source_sha256": {
                        str(Path(monitor.__file__)): monitor.digest(
                            Path(monitor.__file__)
                        )
                    }
                },
                "latest_update": base.name,
            }
        )
    )
    payload = [p for p in base.rglob("*") if p.is_file()]
    (base / "MANIFEST.sha256").write_text(
        "".join(f"{monitor.digest(p)}  {p.relative_to(base)}\n" for p in payload)
    )
    original = {str(p): p.read_bytes() for p in payload}
    monkeypatch.setattr(view, "verify_mapping", lambda case: object())
    monkeypatch.setattr(view, "substrate_indices", lambda *args: {1})
    monkeypatch.setattr(view, "vmd_script", lambda *args: "corrected view")
    result = view.refresh_view(campaign, base)
    assert result.name == "update_000002"
    assert (campaign / "progress/latest").resolve() == result
    assert (result / "vmd/view_path.vmd").read_text() == "corrected view"
    assert (result / "vmd/path.dcd").stat().st_ino == (
        bundle / "path.dcd"
    ).stat().st_ino
    assert all(Path(path).read_bytes() == content for path, content in original.items())
    assert (
        "  accepted_iterations.json\n" not in (result / "MANIFEST.sha256").read_text()
    )
    report = json.loads((result / "progress.json").read_text())
    assert report["remote_evidence_sha256"][
        "accepted_iterations.json"
    ] == monitor.digest(base / "accepted_iterations.json")
    assert report["view_correction"]["coordinate_payload_reused"] is True
    view.verify(base)
    view.verify(result)
    view.verify(result / "vmd")
    (bundle / "path.dcd").write_text("tampered")
    with pytest.raises(ValueError, match="checksum mismatch"):
        view.verify(base)


def test_view_only_cli_journals_success_and_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from scripts import refresh_relay_campaign_view as view

    calls = []
    monkeypatch.setattr(view, "journal", lambda *args: calls.append(args))
    monkeypatch.setattr(
        view.sys,
        "argv",
        ["view", "--campaign", str(tmp_path), "--base", str(tmp_path / "base")],
    )
    monkeypatch.setattr(view, "refresh_view", lambda *args: tmp_path / "result")
    view.main()
    assert calls[-1][0] == "outcome"

    def fail(*args: object) -> Path:
        raise ValueError("bad base")

    monkeypatch.setattr(view, "refresh_view", fail)
    with pytest.raises(ValueError, match="bad base"):
        view.main()
    assert "failed" in calls[-1][1]
