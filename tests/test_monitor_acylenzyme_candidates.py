"""Mock all process launches while testing coherent-band candidate monitoring."""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from test_neb_progress import native_report

from nyl12_covalent.neb_progress import parse_last_complete_neb_report
from scripts import acylenzyme_candidate_tools as tools
from scripts import monitor_acylenzyme_candidates as monitor


def campaign(tmp_path: Path) -> Path:
    """Create fixed-reference campaign configuration."""
    config = tmp_path / "config.json"
    cfg = {
        "sander": "/amber/bin/sander",
        "topology": str(tmp_path / "top.parm7"),
        "base_input": str(tmp_path / "base.mdin"),
        "reactant": str(tmp_path / "ref.ncrst"),
    }
    for key in ("topology", "reactant"):
        Path(cfg[key]).write_bytes(b"input")
    Path(cfg["base_input"]).write_text("&cntrl nmropt=0,/\n&dlfind maxcycle=20,/")
    config.write_text(json.dumps(cfg))
    return config


def mock_snapshot(destination: Path, report_index: int, cycle: int) -> object:
    """Create a hash-consistent immutable fake band."""
    import hashlib

    destination.mkdir()
    rows = []
    for image in range(1, 6):
        path = destination / f"optimized.{image:04d}.ncrst"
        path.write_bytes(f"image {image}".encode())
        rows.append(
            {
                "image": image,
                "name": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    (destination / "band_snapshot.json").write_text(json.dumps({"restarts": rows}))
    report = parse_last_complete_neb_report(native_report([-10, -8, -9, -7, -11]), 5)
    return replace(report, report_index=report_index, cycle=cycle)


def test_observation_policy_is_distinct_bounded_and_provisional() -> None:
    state = {"reserved_images": []}
    assert monitor.observe_candidates(state, 1, 4, (2, 3, 4, 5)) == []
    assert monitor.observe_candidates(state, 1, 5, (2, 3)) == []
    assert monitor.observe_candidates(state, 2, 5, (2, 3, 4, 5)) == [2, 3, 4]
    state["reserved_images"] = [2, 3]
    assert monitor.observe_candidates(state, 3, 6, (2, 3, 4, 5)) == [4]
    state["reserved_images"].append(4)
    assert monitor.observe_candidates(state, 4, 7, (5,)) == []


def test_proc_selection_requires_unique_exact_cwd(tmp_path: Path) -> None:
    proc = tmp_path / "proc"
    proc.mkdir()
    native = tmp_path / "native"
    native.mkdir()
    executable = tmp_path / "sander"
    executable.touch()
    (proc / "self").mkdir()
    (proc / "99").mkdir()
    assert monitor.find_sander_pid(native, proc) is None
    (proc / "99/exe").symlink_to(executable)
    (proc / "99/cwd").symlink_to(native)
    assert monitor.find_sander_pid(native, proc) == 99
    (proc / "100").mkdir()
    (proc / "100/exe").symlink_to(executable)
    (proc / "100/cwd").symlink_to(native)
    with pytest.raises(RuntimeError, match="Multiple"):
        monitor.find_sander_pid(native, proc)


def test_monitor_reserves_only_matched_persistent_candidate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = campaign(tmp_path)
    monkeypatch.setattr(monitor, "audit", lambda *args: None)
    monkeypatch.setattr(monitor, "find_sander_pid", lambda path: 123)
    observed = iter([(1, 5), (2, 10), (2, 10)])

    def snapshot(
        mdout: Path, base: Path, dest: Path, cfg: dict, count: int, **kwargs: object
    ) -> object:
        index, cycle = next(observed)
        return mock_snapshot(dest, index, cycle)

    monkeypatch.setattr(monitor, "snapshot_complete_band_from_xyz", snapshot)
    exports = []
    monkeypatch.setattr(
        monitor,
        "export_band",
        lambda cfg, snap, dest, cycle: exports.append(cycle) or dest,
    )
    launches = []
    monkeypatch.setattr(
        monitor, "launch_candidate", lambda cfg, cand, journal: launches.append(cand)
    )
    for _ in range(3):
        monitor.monitor_once(config, tmp_path / "journal")
    assert len(launches) == 1
    assert launches[0].name == "image_03"
    assert exports == [5]
    metadata = json.loads((launches[0] / "candidate.json").read_text())
    assert metadata["candidate_image"] == 3
    assert (launches[0] / "candidate.ncrst").read_bytes() == b"image 3"
    state = json.loads((tmp_path / "candidate_monitor_state.json").read_text())
    assert state["reserved_images"] == [3]


def test_snapshot_waits_without_launch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = campaign(tmp_path)
    monkeypatch.setattr(monitor, "find_sander_pid", lambda path: None)
    assert "waiting" in monitor.monitor_once(config, tmp_path / "journal")
    monkeypatch.setattr(monitor, "find_sander_pid", lambda path: 99)
    monkeypatch.setattr(monitor, "audit", lambda *args: None)

    def unavailable(*args: object, **kwargs: object) -> None:
        raise ValueError("unfinished save loop")

    monkeypatch.setattr(monitor, "snapshot_complete_band_from_xyz", unavailable)
    assert "deferred" in monitor.monitor_once(config, tmp_path / "journal")
    assert not (tmp_path / "candidate_minima").exists()


def test_export_has_dynamic_leaving_bonds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = campaign(tmp_path)
    snapshot = tmp_path / "snapshot"
    mock_snapshot(snapshot, 1, 10)
    monkeypatch.setattr(tools, "run_cpptraj", lambda *args: None)
    monkeypatch.setattr(tools, "dcd_counts", lambda path: (5, 19967))
    generation = tools.export_band(
        json.loads(config.read_text()), snapshot, tmp_path / "vmd_progress", 10
    )
    script = (generation / "reaction_path.vmd").read_text()
    for row in json.loads((snapshot / "band_snapshot.json").read_text())["restarts"]:
        retained = generation / "source_restarts" / row["name"]
        assert retained.read_bytes() == (snapshot / row["name"]).read_bytes()
        assert (
            str(retained.resolve()) in (generation / "trajectory.cpptraj").read_text()
        )
    assert "remove_reactive_bond $molid 19915 19928" in script
    assert 'mol selection "index 13916 19928"' in script
    assert "label add Bonds $molid/19915 $molid/19928" in script
    assert (tmp_path / "vmd_progress/latest").resolve() == generation
    monkeypatch.setattr(tools, "dcd_counts", lambda path: (4, 19967))
    with pytest.raises(ValueError, match="frame/atom"):
        tools.export_band(
            json.loads(config.read_text()), snapshot, tmp_path / "bad_vmd", 10
        )


def test_candidate_launch_command_and_unbiased_release(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = campaign(tmp_path)
    snapshot = tmp_path / "snapshot"
    mock_snapshot(snapshot, 1, 10)
    candidate = tools.prepare_candidate(snapshot, tmp_path, 3)
    commands = []
    monkeypatch.setattr(
        tools.subprocess, "run", lambda cmd, **kwargs: commands.append(cmd)
    )
    tools.launch_candidate(config, candidate, tmp_path / "journal")
    assert "--journal" in commands[0]
    assert "--release" in commands[0]
    assert "launch" in commands[0]

    def engine(cfg: dict, folder: Path, start: Path, text: str) -> Path:
        assert "DISANG" not in text and "nmropt=0" in text
        assert f"refcrd='{tmp_path / 'ref.ncrst'}'" in text
        folder.mkdir()
        result = folder / "optimized.ncrst"
        result.write_bytes(b"released")
        return result

    monkeypatch.setattr(tools, "run_engine", engine)
    monkeypatch.setattr(tools, "geometry", lambda path: {"OG1_C": 4.0})
    monkeypatch.setattr(tools, "native_converged", lambda path: True)
    tools.release_candidate(config, candidate)
    accepted = json.loads((candidate / "trial_acceptance.json").read_text())
    assert accepted["native_converged"]
    assert not accepted["acylenzyme_gate_applied"]
    assert not accepted["stationary_point_claim"]
    (candidate / "candidate.ncrst").write_bytes(b"changed")
    with pytest.raises(ValueError, match="hash changed"):
        tools.release_candidate(config, candidate)
    (snapshot / "optimized.0002.ncrst").write_bytes(b"changed")
    with pytest.raises(ValueError, match="snapshot hash"):
        tools.prepare_candidate(snapshot, tmp_path, 2)


def test_main_once_and_release(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = campaign(tmp_path)
    monkeypatch.setattr(monitor, "audit", lambda *args: None)
    monkeypatch.setattr(monitor, "monitor_once", lambda *args: "waiting")
    assert monitor.main(["--config", str(config), "--once"]) == 0
    calls = []
    monkeypatch.setattr(monitor, "release_candidate", lambda *args: calls.append(args))
    assert monitor.main(["--config", str(config), "--release", str(tmp_path)]) == 0
    assert calls
    with pytest.raises(SystemExit):
        monitor.main(["--config", str(config), "--interval-seconds", "0"])


def test_duplicate_monitor_cannot_launch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import fcntl

    config = campaign(tmp_path)
    monkeypatch.setattr(monitor, "audit", lambda *args: None)
    with (tmp_path / "candidate_monitor.lock").open("a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RuntimeError, match="campaign lock"):
            monitor.main(["--config", str(config), "--once"])


@pytest.mark.parametrize("phase", ["finished", "stopped"])
def test_monitor_stops_after_terminal_campaign(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, phase: str
) -> None:
    config = campaign(tmp_path)
    (tmp_path / "status.json").write_text(json.dumps({"phase": phase}))
    monkeypatch.setattr(monitor, "find_sander_pid", lambda path: None)
    monkeypatch.setattr(monitor, "audit", lambda *args: None)
    assert (
        monitor.monitor_once(config, tmp_path / "journal")
        == f"campaign terminal: {phase}"
    )
    assert monitor.main(["--config", str(config)]) == 0
