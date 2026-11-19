"""Tests for the corrected post-proton-transfer attack runner."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import run_corrected_attack_to_ti1 as runner
from scripts.run_corrected_attack_to_ti1 import _proton_identity_retained


def test_proton_identity_gate_accepts_three_terminal_n_h_bonds() -> None:
    values = {
        "N_HG1_A": 1.02,
        "N_H1_A": 1.01,
        "N_H2_A": 1.03,
        "OG1_HG1_A": 2.20,
        "OG1_H1_A": 2.60,
        "OG1_H2_A": 2.90,
    }
    assert _proton_identity_retained(values)


def test_proton_identity_gate_rejects_any_swapped_proton() -> None:
    values = {
        "N_HG1_A": 1.02,
        "N_H1_A": 1.52,
        "N_H2_A": 1.03,
        "OG1_HG1_A": 2.20,
        "OG1_H1_A": 1.01,
        "OG1_H2_A": 2.90,
    }
    assert not _proton_identity_retained(values)


@pytest.mark.parametrize("attack_only", [False, True])
def test_run_publishes_accepted_sequential_attack_windows(
    tmp_path: Path, monkeypatch, attack_only: bool
) -> None:
    sources = {
        name: tmp_path / name for name in ("topology", "start", "sander", "cpptraj")
    }
    for path in sources.values():
        path.write_text(path.name)
    base_mdin = tmp_path / "base.mdin"
    base_mdin.write_text(
        "&cntrl\n nmropt=0, ntr=0, qm_theory='DFTB3', maxcyc=10,\n/\n"
        "&dlfind active='@*', maxcycle=10, /\n"
    )
    output = tmp_path / "corrected_attack"
    initial = {
        "OG1_carbonyl_C_A": 2.1036,
        "carbonyl_C_O_A": 1.24,
        "carbonyl_C_N_A": 1.39,
        "N_HG1_A": 1.02,
        "N_H1_A": 1.01,
        "N_H2_A": 1.03,
        "OG1_HG1_A": 2.20,
        "OG1_H1_A": 2.60,
        "OG1_H2_A": 2.90,
    }
    schedule = runner.attack_only_schedule(initial["OG1_carbonyl_C_A"])
    targets = {row["label"]: row["attack_target_A"] for row in schedule}

    def fake_geometry(_topology: Path, coordinates: Path) -> dict[str, float]:
        values = dict(initial)
        if coordinates.name == "optimized.ncrst":
            values["OG1_carbonyl_C_A"] = (
                1.50
                if coordinates.parent.name == "02_unrestrained_endpoint"
                else targets[coordinates.parent.name]
            )
            values["carbonyl_C_O_A"] = 1.31
        return values

    def fake_engine(
        _command: list[str], folder: Path, _environment: dict[str, str]
    ) -> None:
        (folder / "optimized.ncrst").write_text(folder.name)
        (folder / "minimized.mdout").write_text("converged")

    monkeypatch.setattr(runner, "coupled_geometry", fake_geometry)
    monkeypatch.setattr(runner, "_run_engine", fake_engine)
    monkeypatch.setattr(runner, "native_converged", lambda _path: True)
    monkeypatch.setattr(
        runner, "export", lambda *_args: output / "vmd_progress" / "latest"
    )
    runner.run(
        argparse.Namespace(
            topology=sources["topology"],
            start=sources["start"],
            base_mdin=base_mdin,
            sander=sources["sander"],
            cpptraj=sources["cpptraj"],
            output=output,
            final_attack=1.50,
            maximum_increment=0.075,
            max_cycles=2000,
            attack_only=attack_only,
            release_endpoint=attack_only,
        )
    )

    status = json.loads((output / "status.json").read_text())
    assert status["state"] == "completed"
    assert status["completed_windows"] == [row["label"] for row in schedule]
    assert status["endpoint_is_stationary_point"] is False
    for row in schedule:
        acceptance = json.loads(
            (output / "windows" / row["label"] / "acceptance.json").read_text()
        )
        assert acceptance["accepted"] is True
        assert acceptance["proton_identity_retained"] is True
        restraints = (output / "windows" / row["label"] / "restraints.RST").read_text()
        assert restraints.count("&rst") == (1 if attack_only else 4)
    if attack_only:
        release = output / "02_unrestrained_endpoint"
        assert "nmropt=0" in (release / "minimize.mdin").read_text()
        assert not (release / "restraints.RST").exists()
        assert json.loads((release / "acceptance.json").read_text())["accepted"]


@pytest.mark.parametrize("exit_status", [0, 2])
def test_run_engine_records_process_outcome(
    tmp_path: Path, monkeypatch, exit_status: int
) -> None:
    folder = tmp_path / "window"
    folder.mkdir()
    monkeypatch.setattr(
        runner.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=exit_status),
    )
    if exit_status:
        with pytest.raises(RuntimeError, match="Sander failed"):
            runner._run_engine(["sander"], folder, {"OMP_NUM_THREADS": "1"})
    else:
        runner._run_engine(["sander"], folder, {"OMP_NUM_THREADS": "1"})
    record = json.loads((folder / "process.json").read_text())
    assert record["exit_status"] == exit_status
    assert record["command"] == ["sander"]


def test_resume_promotes_only_converged_identity_safe_tracking_stop(
    tmp_path: Path,
) -> None:
    output = tmp_path / "calculation"
    folder = output / "windows" / "attack_01"
    folder.mkdir(parents=True)
    endpoint = folder / "optimized.ncrst"
    endpoint.write_text("endpoint")
    (folder / "acceptance.json").write_text(
        json.dumps(
            {
                "accepted": False,
                "native_converged": True,
                "proton_identity_retained": True,
                "attack_tracking_error_A": 0.115,
            }
        )
    )
    (output / "status.json").write_text(
        json.dumps(
            {
                "state": "stopped",
                "completed_windows": [],
                "reason": "Scientific gate failed in attack_01",
            }
        )
    )
    rows = [{"label": "attack_01", "attack_target_A": 1.9}]

    status, restart, next_index = runner._resume_state(output, rows)

    assert restart == endpoint
    assert next_index == 1
    assert status["state"] == "running"
    assert status["completed_windows"] == ["attack_01"]
    assert "reason" not in status
    revision = json.loads((folder / "acceptance_revision.json").read_text())
    assert revision["accepted_for_continuation"] is True
    assert revision["attack_tracking_error_A"] == 0.115


@pytest.mark.parametrize(
    ("status", "rows", "make_endpoint", "message"),
    [
        (None, [{"label": "attack_01"}], False, "status.json"),
        (
            {"state": "running", "completed_windows": []},
            [{"label": "attack_01"}],
            False,
            "stopped corrected-attack",
        ),
        (
            {"state": "stopped", "completed_windows": ["wrong"]},
            [{"label": "attack_01"}],
            False,
            "disagrees with the schedule",
        ),
        (
            {"state": "stopped", "completed_windows": ["attack_01"]},
            [{"label": "attack_01"}],
            False,
            "No corrected-attack windows remain",
        ),
        (
            {"state": "stopped", "completed_windows": []},
            [{"label": "attack_01"}],
            True,
            "not safe",
        ),
    ],
)
def test_resume_rejects_invalid_state(
    tmp_path: Path,
    status: dict[str, object] | None,
    rows: list[dict[str, str]],
    make_endpoint: bool,
    message: str,
) -> None:
    output = tmp_path / "calculation"
    output.mkdir()
    if status is not None:
        (output / "status.json").write_text(json.dumps(status))
    if make_endpoint:
        folder = output / "windows" / "attack_01"
        folder.mkdir(parents=True)
        (folder / "optimized.ncrst").write_text("endpoint")
        (folder / "acceptance.json").write_text(
            json.dumps(
                {
                    "native_converged": False,
                    "proton_identity_retained": True,
                    "attack_tracking_error_A": 0.1,
                }
            )
        )
    with pytest.raises((FileNotFoundError, ValueError), match=message):
        runner._resume_state(output, rows)


@pytest.mark.parametrize("plots", [False, True])
def test_publish_calls_only_one_export_route(
    tmp_path: Path, monkeypatch, plots: bool
) -> None:
    calls = []
    bundle = tmp_path / "bundle"
    monkeypatch.setattr(
        runner, "export", lambda *_args: calls.append("export") or bundle
    )
    monkeypatch.setattr(
        runner, "publish_progress", lambda *_args: calls.append("plots") or bundle
    )
    args = argparse.Namespace(
        output=tmp_path,
        topology=tmp_path,
        cpptraj=tmp_path,
        plot_python=tmp_path if plots else None,
    )
    assert runner._publish(args) == bundle
    assert calls == ["plots" if plots else "export"]


@pytest.mark.parametrize("accepted", [False, True])
def test_resume_after_publication_failure_uses_only_accepted_checkpoint(
    tmp_path: Path, accepted: bool
) -> None:
    folder = tmp_path / "windows/attack_01"
    folder.mkdir(parents=True)
    endpoint = folder / "optimized.ncrst"
    endpoint.write_text("checkpoint")
    (folder / "acceptance.json").write_text(
        json.dumps(
            {
                "accepted": accepted,
                "native_converged": True,
                "proton_identity_retained": True,
                "scissile_connectivity_retained": True,
            }
        )
    )
    (tmp_path / "status.json").write_text(
        json.dumps(
            {
                "state": "stopped",
                "completed_windows": ["attack_01"],
                "reason": "duplicate publication",
            }
        )
    )
    rows = [{"label": "attack_01"}, {"label": "attack_02"}]
    if not accepted:
        with pytest.raises(ValueError, match="not safe"):
            runner._resume_state(tmp_path, rows)
        return
    status, restart, index = runner._resume_state(tmp_path, rows)
    assert (restart, index) == (endpoint, 1)
    assert status["completed_windows"] == ["attack_01"]
    assert status["resume_history"][0]["accepted_checkpoint"] == "attack_01"
    assert not (folder / "acceptance_revision.json").exists()
