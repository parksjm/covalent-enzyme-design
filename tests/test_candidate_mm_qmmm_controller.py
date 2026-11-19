"""Tests for the candidate MM continuation and conditional QM/MM handoff."""

import json
from pathlib import Path

import pytest

from nyl12_covalent.candidate_mm_qmmm_controller import (
    CandidateControllerError,
    run_candidate,
)
from scripts.run_modular_es_candidate import launch_qmmm_shells, validate_active_radii


def _amber_files(tmp_path: Path) -> tuple[Path, Path, Path]:
    paths = tuple(
        tmp_path / name for name in ("system.prmtop", "start.ncrst", "pmemd.cuda")
    )
    for path in paths:
        path.write_text("fixture\n", encoding="utf-8")
    return paths


def _fake_runner_factory(rms_values: list[float], commands: list[list[str]]):
    remaining = iter(rms_values)

    def run(command: list[str]) -> int:
        commands.append(list(command))

        def get_value(flag: str) -> Path:
            return Path(command[command.index(flag) + 1])

        mdout = get_value("-o")
        restart = get_value("-r")
        mdout.write_text(
            f"FINAL RESULTS\n NSTEP ENERGY RMS GMAX\n 10 -1.0 {next(remaining)} 0.05\n",
            encoding="utf-8",
        )
        restart.write_text("restart\n", encoding="utf-8")
        return 0

    return run


def _run(
    tmp_path: Path,
    rms_values: list[float],
    gate_accepted: bool = True,
) -> tuple[dict, list[list[str]], list[tuple[Path, Path]]]:
    topology, coordinates, engine = _amber_files(tmp_path)
    commands: list[list[str]] = []
    qmmm_calls: list[tuple[Path, Path]] = []
    report = run_candidate(
        topology,
        coordinates,
        engine,
        tmp_path / "selected_candidate_run",
        "selected WT frame",
        segment_cycles=10,
        total_cycle_ceiling=20,
        endpoint_gate=lambda _endpoint: {"accepted": gate_accepted},
        qmmm_launcher=lambda mm, output: (
            qmmm_calls.append((mm, output)) or {"exit_status": 0}
        ),
        engine_runner=_fake_runner_factory(rms_values, commands),
    )
    return report, commands, qmmm_calls


def test_converged_first_segment_runs_gate_then_qmmm(tmp_path: Path) -> None:
    report, commands, qmmm_calls = _run(tmp_path, [0.009])

    assert report["status"] == "complete"
    assert len(report["segments"]) == 1
    assert len(commands) == 1
    assert report["qmmm_launched"] is True
    assert len(qmmm_calls) == 1
    assert qmmm_calls[0][0].name == "segment_001_cycles_00010"
    mdin = qmmm_calls[0][0] / "minimize.mdin"
    assert "ntr=0" in mdin.read_text(encoding="utf-8")
    assert "nmropt=0" in mdin.read_text(encoding="utf-8")


def test_unconverged_segment_continues_from_immutable_restart(tmp_path: Path) -> None:
    report, commands, qmmm_calls = _run(tmp_path, [0.2, 0.005])

    assert [row["cycle_budget"] for row in report["segments"]] == [10, 10]
    assert len(commands) == 2
    first_restart = Path(commands[0][commands[0].index("-r") + 1]).resolve()
    second_input = Path(commands[1][commands[1].index("-c") + 1]).resolve()
    assert second_input == first_restart
    assert report["qmmm_launched"] is True
    assert len(qmmm_calls) == 1


def test_total_cycle_ceiling_fails_without_qmmm(tmp_path: Path) -> None:
    topology, coordinates, engine = _amber_files(tmp_path)
    commands: list[list[str]] = []
    qmmm_calls: list[tuple[Path, Path]] = []
    output = tmp_path / "selected_candidate_run"

    with pytest.raises(CandidateControllerError, match="cycle ceiling"):
        run_candidate(
            topology,
            coordinates,
            engine,
            output,
            "selected WT frame",
            segment_cycles=10,
            total_cycle_ceiling=15,
            endpoint_gate=lambda _endpoint: {"accepted": True},
            qmmm_launcher=lambda mm, out: qmmm_calls.append((mm, out) or {}),
            engine_runner=_fake_runner_factory([0.2, 0.1], commands),
        )

    report = json.loads((output / "controller_report.json").read_text())
    assert report["status"] == "mm_cycle_ceiling_failed"
    assert [row["cycle_budget"] for row in report["segments"]] == [10, 5]
    assert qmmm_calls == []


def test_failed_scientific_gate_blocks_qmmm(tmp_path: Path) -> None:
    topology, coordinates, engine = _amber_files(tmp_path)
    commands: list[list[str]] = []
    qmmm_calls: list[tuple[Path, Path]] = []
    output = tmp_path / "selected_candidate_run"

    with pytest.raises(CandidateControllerError, match="gate"):
        run_candidate(
            topology,
            coordinates,
            engine,
            output,
            "selected WT frame",
            segment_cycles=10,
            total_cycle_ceiling=20,
            endpoint_gate=lambda _endpoint: {"accepted": False, "tier": "none"},
            qmmm_launcher=lambda mm, out: qmmm_calls.append((mm, out) or {}),
            engine_runner=_fake_runner_factory([0.005], commands),
        )

    report = json.loads((output / "controller_report.json").read_text())
    assert report["status"] == "scientific_gate_failed"
    assert report["qmmm_launched"] is False
    assert qmmm_calls == []


def test_one_accepted_mm_endpoint_launches_both_shells(tmp_path: Path) -> None:
    topology, coordinates, engine = _amber_files(tmp_path)
    mm_commands: list[list[str]] = []
    queue_commands: list[list[str]] = []
    results = iter([0, 0])

    def run_queue(command: list[str], _cwd: Path) -> int:
        queue_commands.append(command)
        return next(results)

    def launch(mm_directory: Path, output: Path) -> dict:
        return launch_qmmm_shells(
            queue_script=Path("queue.py"),
            common_arguments=["--mm-directory", str(mm_directory)],
            active_radii=[8.0, 12.0],
            output=output,
            cwd=tmp_path,
            runner=run_queue,
        )

    report = run_candidate(
        topology=topology,
        initial_coordinates=coordinates,
        engine=engine,
        output=tmp_path / "selected_candidate_run",
        case_label="selected WT frame",
        segment_cycles=10,
        total_cycle_ceiling=20,
        endpoint_gate=lambda _endpoint: {"accepted": True},
        qmmm_launcher=launch,
        engine_runner=_fake_runner_factory([0.009], mm_commands),
    )

    assert report["status"] == "complete"
    assert len(report["segments"]) == 1
    assert len(mm_commands) == 1
    assert len(queue_commands) == 2
    endpoint = Path(report["segments"][0]["endpoint_coordinates"])
    assert all(
        Path(command[command.index("--mm-directory") + 1]) / "minimized.ncrst"
        == endpoint
        for command in queue_commands
    )
    assert [
        command[command.index("--active-radius") + 1] for command in queue_commands
    ] == [
        "8.0",
        "12.0",
    ]
    outputs = [command[command.index("--output") + 1] for command in queue_commands]
    assert outputs[0].endswith("qmmm_dftb3_shell_8A")
    assert outputs[1].endswith("qmmm_dftb3_shell_12A")
    assert outputs[0] != outputs[1]
    assert [shell["exit_status"] for shell in report["qmmm"]["shells"]] == [0, 0]
    persisted = json.loads(
        (
            tmp_path
            / "selected_candidate_run"
            / "qmmm_shells"
            / "qmmm_shell_launch_report.json"
        ).read_text()
    )
    assert [shell["exit_status"] for shell in persisted["shells"]] == [0, 0]


def test_shell_launch_stops_and_records_nonzero_queue(tmp_path: Path) -> None:
    commands: list[list[str]] = []

    def run_queue(command: list[str], _cwd: Path) -> int:
        commands.append(command)
        return 7

    with pytest.raises(RuntimeError, match="8 A shell exited with status 7"):
        launch_qmmm_shells(
            queue_script=Path("queue.py"),
            common_arguments=[],
            active_radii=[8.0, 12.0],
            output=tmp_path / "qmmm_shells",
            cwd=tmp_path,
            runner=run_queue,
        )

    assert len(commands) == 1
    persisted = json.loads(
        (tmp_path / "qmmm_shells" / "qmmm_shell_launch_report.json").read_text()
    )
    assert persisted["shells"][0]["exit_status"] == 7


@pytest.mark.parametrize("radii", [[], [0.0], [-1.0], [8.0, 8.0], [float("nan")]])
def test_active_radii_must_be_unique_positive_finite_values(
    radii: list[float],
) -> None:
    with pytest.raises(ValueError):
        validate_active_radii(radii)
