"""Tests for single-model Amber orchestration without external execution."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent import design_model_relaxation as runner

RELEASE = Path(
    "/data/jpc/polyamide-amber-ff/parameters/releases/"
    "pa6_l4_modular_resp_gaff2_v1_20261031"
)


class FakeAmber:
    """Small ParmEd stand-in for box and endpoint handling."""

    box = (10.0, 20.0, 30.0, 90.0, 90.0, 90.0)

    def strip(self, _mask: str) -> None:
        self.stripped = True

    def save(self, filename: str, **_kwargs: Any) -> None:
        Path(filename).write_text("synthetic artifact\n")


def _config(tmp_path: Path) -> dict[str, str]:
    """Create explicit mock Amber paths and immutable parameter inputs."""
    amberhome = tmp_path / "amber"
    amberhome.mkdir()
    mol2 = tmp_path / "source.mol2"
    mapping = tmp_path / "mapping.json"
    tleap = tmp_path / "tleap"
    engine = tmp_path / "bin" / "pmemd.cuda"
    engine.parent.mkdir()
    for path in (mol2, mapping, tleap, engine):
        path.write_text("test input\n")
    return {
        "release": str(RELEASE),
        "mol2": str(mol2),
        "mapping": str(mapping),
        "tleap": str(tleap),
        "engine": str(engine),
        "amberhome": str(amberhome),
        "gpu": "2",
    }


def _mock_mapping(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace only scientific mapping/export with deterministic artifacts."""

    def convert(
        source: Path, _mol2: Path, _mapping: Path, _release: Path, system: Path
    ) -> dict[str, Any]:
        system.mkdir()
        (system / "conversion.json").write_text(json.dumps({"atom_mapping": []}))
        (system / "complex_heavy_modular.pdb").write_text(source.read_text())
        return {"status": "accepted"}

    def export(
        _top: Path, _coords: Path, _mapping: Path, output: Path
    ) -> dict[str, str]:
        output.write_text("ATOM\n")
        return {"output": str(output)}

    monkeypatch.setattr(runner, "convert_four_ligands", convert)
    monkeypatch.setattr(runner, "export_canonical_heavy", export)


def _mock_subprocess(
    monkeypatch: pytest.MonkeyPatch, *, leap_error: bool = False, rms: str = "0.03"
) -> list[dict[str, Any]]:
    """Emulate LEaP and pmemd file outputs while retaining exact argv/env."""
    calls: list[dict[str, Any]] = []

    def run(argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[bytes]:
        calls.append({"argv": argv, "env": kwargs["env"], "cwd": kwargs["cwd"]})
        if Path(argv[0]).name == "tleap":
            stem = "noions" if "count.in" in argv else "complex"
            directory = Path(kwargs["cwd"])
            for suffix in ("prmtop", "rst7"):
                (directory / f"{stem}.{suffix}").write_text("fixture\n")
            (directory / "dry.prmtop").write_text("fixture\n")
            (directory / "dry.rst7").write_text("fixture\n")
            (directory / "noions.prmtop").write_text("fixture\n")
            (directory / "noions.rst7").write_text("fixture\n")
            kwargs["stdout"].write(
                "Errors = 1\n" if leap_error and "count.in" in argv else "Errors = 0\n"
            )
        else:
            output = Path(argv[argv.index("-o") + 1])
            restart = Path(argv[argv.index("-r") + 1])
            output.write_text(f"FINAL RESULTS\n 100 -123.0 {rms} 0.4 N 1\n")
            restart.write_text("restart\n")
            kwargs["stdout"].write("pmemd complete\n")
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(runner.subprocess, "run", run)
    return calls


def _patch_science(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mock topology acceptance and return a deterministic atom inventory."""
    monkeypatch.setattr(
        runner,
        "validate_four_ligand_topology",
        lambda _path, _release, solvated: {
            "release_id": "pa6_l4_modular_resp_gaff2_v1_20261031",
            "atom_count": 9000 if solvated else 8000,
            "solvated": solvated,
        },
    )
    monkeypatch.setattr(runner, "_volume", lambda _top, _coords: 100000.0)
    monkeypatch.setattr(
        runner.parmed, "load_file", lambda *_args, **_kwargs: FakeAmber()
    )


def test_relax_model_records_setup_commands_and_nonconverged_endpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Run mocked orchestration through endpoint export without claiming convergence."""
    source = tmp_path / "source.pdb"
    source.write_text("source\n")
    config = _config(tmp_path)
    _mock_mapping(monkeypatch)
    _patch_science(monkeypatch)
    calls = _mock_subprocess(monkeypatch, rms="0.08")
    output = tmp_path / "run"
    result = runner.relax_model(source, output, config)
    assert len(calls) == 3
    assert calls[0]["env"]["CUDA_VISIBLE_DEVICES"] == "2"
    assert calls[0]["env"]["OMP_NUM_THREADS"] == "1"
    assert result["force_converged"] is False
    assert result["salt"]["sodium_count"] == result["salt"]["chloride_count"] + 8
    mdin = (output / "minimization/min.in").read_text()
    assert "maxcyc=10000" in mdin and "drms=0.05" in mdin
    assert "ntr=0, nmropt=0, ntc=1, ntf=1" in mdin
    assert "ntwprt=9000" in mdin
    assert len((output / "command_log.jsonl").read_text().splitlines()) == 3
    assert (
        json.loads((output / "process_acceptance.json").read_text())["status"]
        == "completed"
    )
    assert (output / "acceptance.json").is_file()
    assert (output / "relaxed.pdb").is_file()


def test_relax_model_records_and_raises_on_leap_errors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Do not proceed when LEaP exits zero but reports errors."""
    source = tmp_path / "source.pdb"
    source.write_text("source\n")
    _mock_mapping(monkeypatch)
    _patch_science(monkeypatch)
    _mock_subprocess(monkeypatch, leap_error=True)
    output = tmp_path / "failed"
    with pytest.raises(ValueError, match="LEaP log is not clean"):
        runner.relax_model(source, output, _config(tmp_path))
    process = json.loads((output / "process_acceptance.json").read_text())
    assert process["status"] == "failed"
    assert process["stage"] == "dry_leap"
    assert not (output / "acceptance.json").exists()


def test_box_volume_and_amber_final_results_validation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Validate triclinic volume and fail on missing or non-finite RMS data."""
    assert round(runner.ION_PAIR_FACTOR * 100000.0) == 9
    monkeypatch.setattr(
        runner.parmed, "load_file", lambda *_args, **_kwargs: FakeAmber()
    )
    assert runner._volume(Path("top"), Path("coords")) == pytest.approx(6000.0)
    energy, rms = runner._final_results("FINAL RESULTS\n 10 -12.3 0.05 0.1 N 2\n")
    assert energy == pytest.approx(-12.3)
    assert rms == pytest.approx(0.05)
    with pytest.raises(ValueError, match="FINAL RESULTS"):
        runner._final_results("no result")
    with pytest.raises(ValueError, match="non-finite"):
        runner._final_results("FINAL RESULTS\n 10 1e999 0.1 0.2 N 2\n")


def test_existing_output_and_bad_gpu_fail_before_execution(tmp_path: Path) -> None:
    """Reject resumed output directories and missing explicit GPU selection."""
    source = tmp_path / "source.pdb"
    source.write_text("source\n")
    output = tmp_path / "exists"
    output.mkdir()
    with pytest.raises(FileExistsError):
        runner.relax_model(source, output, {})
    config = _config(tmp_path)
    config["gpu"] = " "
    with pytest.raises(ValueError, match="gpu"):
        runner.relax_model(source, tmp_path / "bad", config)
