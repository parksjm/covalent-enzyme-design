"""Focused no-calculation tests for native-Sander path orchestration."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from parmed.amber import Rst7

from nyl12_covalent import native_irc_driver as driver
from nyl12_covalent import native_sander_cli
from nyl12_covalent.native_sander_cli import (
    FORCE_SOURCE_FILES,
    PATH_FORCE_TO_GRADIENT_SCALE,
    NativeSanderEvaluator,
    load_calibration,
    prepare_native_input,
    source_sha256,
)

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / (
    "qmmm/pa6_flipped_model30/05_neutral_reaction_path/"
    "07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113"
)
SANDER = Path("/data/jpc/code/ambertools26/venv/bin/sander")


def test_native_input_keeps_qmmm_and_disables_temperature_and_dimer() -> None:
    source = (BASE / "calculation/dimer.mdin").read_text()
    prepared = prepare_native_input(source)

    assert "dftb_telec=0.0" in prepared
    assert "dimer=.false." in prepared
    assert "maxcycle=1" in prepared
    for key in (
        "qmmask",
        "qmcharge",
        "qm_theory",
        "qm_ewald",
        "qm_pme",
        "qmcut",
        "active",
    ):
        assert key in prepared
    assert "&dlfind" in prepared
    assert "coords2" not in prepared


def test_calibration_fails_closed_until_validated(tmp_path: Path) -> None:
    path = tmp_path / "calibration.json"
    path.write_text(json.dumps({"status": "pending", "passed": False}))

    with pytest.raises(ValueError, match="not validated"):
        load_calibration(path)


def test_native_cli_evaluator_parses_energy_and_calibrates_mock_force(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        native_sander_cli,
        "write_restart",
        lambda path, _coordinates, _restart: path.write_bytes(b"restart"),
    )

    def mock_sander(
        command: list[str], cwd: Path, stdout: object, stderr: object, check: bool
    ) -> None:
        assert check is True
        assert "-O" in command
        output = Path(command[command.index("-o") + 1])
        output.write_text("Energy calculation finished, energy: -1.25D+00\n")
        for name in ("path_force.xyz", "path.xyz", "path_active.xyz"):
            (cwd / name).write_text("2\nraw force\nH 1 2 3\nH -4 -5 -6\n")

    monkeypatch.setattr(native_sander_cli.subprocess, "run", mock_sander)
    active = np.array([True, False, True])
    evaluator = NativeSanderEvaluator(
        tmp_path / "system.prmtop",
        None,
        active,
        "one-cycle input",
        tmp_path / "sander",
        tmp_path,
        PATH_FORCE_TO_GRADIENT_SCALE,
    )
    energy, gradient = evaluator(np.zeros((3, 3)))

    assert energy == pytest.approx(-1.25)
    assert gradient[0].tolist() == pytest.approx(
        (PATH_FORCE_TO_GRADIENT_SCALE * np.array([1, 2, 3])).tolist()
    )
    assert gradient[1].tolist() == [0.0, 0.0, 0.0]
    assert gradient[2].tolist() == pytest.approx(
        (PATH_FORCE_TO_GRADIENT_SCALE * np.array([-4, -5, -6])).tolist()
    )
    for name in ("path_force.xyz", "path.xyz", "path_active.xyz"):
        assert (tmp_path / f"eval_00001.{name}").is_file()
    assert "eval_00001.path.xyz" in evaluator.evidence[0]["artifacts_sha256"]
    resumed = NativeSanderEvaluator(
        tmp_path / "system.prmtop",
        None,
        active,
        "one-cycle input",
        tmp_path / "sander",
        tmp_path,
        PATH_FORCE_TO_GRADIENT_SCALE,
    )
    resumed(np.zeros((3, 3)))
    assert (tmp_path / "eval_00002.mdout").is_file()


def test_two_branch_runner_uses_mock_engine_and_resumes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    class MockEngine:
        def __init__(
            self,
            _topology: Path,
            source_restart: Rst7,
            active: np.ndarray,
            _mdin: str,
            _exe: Path,
            _dir: Path,
            _scale: float,
        ) -> None:
            self.reference = (
                np.asarray(source_restart.coordinates).reshape((-1, 3)).copy()
            )
            self.active = active
            self.evidence: list[dict[str, object]] = []

        def __call__(self, coordinates: np.ndarray) -> tuple[float, np.ndarray]:
            displacement = coordinates - self.reference
            energy = -0.5 * float(np.sum(displacement[self.active] ** 2))
            gradient = np.zeros_like(coordinates)
            gradient[self.active] = -displacement[self.active]
            self.evidence.append({"energy_hartree": energy, "mock": True})
            return energy, gradient

    monkeypatch.setattr(driver, "NativeSanderEvaluator", MockEngine)
    calibration = tmp_path / "calibration.json"
    calibration.write_text(
        json.dumps(
            {
                "status": "validated",
                "passed": True,
                "active_gradient": {
                    "observed_force_to_gradient_scale": PATH_FORCE_TO_GRADIENT_SCALE,
                    "writer_source": {
                        "source_hashes": [
                            {
                                "path": str(path.resolve()),
                                "sha256": source_sha256(path),
                            }
                            for path in FORCE_SOURCE_FILES
                        ]
                    },
                },
                "validation_sha256": "mocked-evidence",
            }
        )
    )
    distances = tmp_path / "distances.json"
    distances.write_text(
        json.dumps([{"label": "sample", "atom_i_1based": 1, "atom_j_1based": 2}])
    )
    output = tmp_path / "path"

    common = {
        "topology_path": BASE / "system.prmtop",
        "restart_path": BASE / "calculation/optimized.ncrst",
        "mdin_path": BASE / "calculation/dimer.mdin",
        "path_xyz": BASE / "calculation/path.xyz",
        "mode_path": BASE / "calculation/path_tsmode.xyz",
        "sander": SANDER,
        "output_dir": output,
        "calibration_path": calibration,
        "distance_pairs_path": distances,
    }
    driver.run_two_branches(**common, accepted_steps=1)
    for branch in ("mode_plus", "mode_minus"):
        rows = [
            json.loads(line)
            for line in (output / branch / "steps.jsonl").read_text().splitlines()
        ]
        assert len(rows) == 1
        assert rows[0]["accepted"] is True
        assert rows[0]["energy_hartree"] < 0.0
        assert "sample" in rows[0]["distances_angstrom"]
        end = json.loads((output / branch / "endpoint_candidate.json").read_text())
        assert end["terminal_status"] == "maximum_accepted_steps"
        assert "not a formal IRC" in end["interpretation"]
        checkpoint = Rst7.open(str(output / branch / "checkpoint.ncrst"))
        assert checkpoint.coordinates.size > 0

    driver.run_two_branches(**common, accepted_steps=2)
    for branch in ("mode_plus", "mode_minus"):
        rows = [
            json.loads(line)
            for line in (output / branch / "steps.jsonl").read_text().splitlines()
        ]
        assert sum(row["accepted"] is True for row in rows) == 2
