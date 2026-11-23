"""Exercise controller plumbing with explicit fake engines, never chemistry."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import pytest

SOURCE = Path(__file__).resolve().parents[1] / (
    "analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/"
    "repaired_comparison_20261120/source_revision03/run_comparison.py"
)
spec = importlib.util.spec_from_file_location("fake_native_orchestration", SOURCE)
assert spec is not None and spec.loader is not None
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def executable(path: Path, body: str) -> Path:
    """Create an explicitly fake executable for subprocess-boundary tests."""
    path.write_text(f"#!{sys.executable}\n" + body)
    path.chmod(0o755)
    return path


@pytest.fixture
def native_case(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple:
    """Build accepted synthetic fixtures and fake native executable artifacts."""
    software = tmp_path / "software"
    bin_dir = software / "bin"
    bin_dir.mkdir(parents=True)
    executable(bin_dir / "delphi", "pass\n")
    executable(
        bin_dir / "mcce",
        """from pathlib import Path
phase=Path('run.prm').read_text()
if phase=='step1': Path('step1_out.pdb').write_text('synthetic plumbing fixture')
if phase=='step2': Path('step2_out.pdb').write_text('synthetic plumbing fixture')
if phase=='step4':
 Path('fort.38').write_text('ph '+' '.join(map(str,range(15)))+'\\nASP01A0001_001 '+' '.join(['1']*15)+'\\n')
 Path('pK.out').write_text('synthetic fake-native pKa output')
 Path('run.prm.record').write_text('synthetic fake-native parameter history')
""",
    )
    executable(
        bin_dir / "step3.py",
        """from pathlib import Path
Path('head3.lst').write_text('header\\n1 ASP01A0001_001 '+' '.join(['0']*14)+' t\\n')
Path('energies').mkdir()
Path('energies/ASP01A0001_001.opp').write_text('fake pair energies')
""",
    )
    heavy = tmp_path / "heavy.pdb"
    heavy.write_text(
        "ATOM      1  OXT GLY A   1       1.000   2.000   3.000  1.00  0.00           O\n"
    )
    original = tmp_path / "original.pdb"
    original.write_text("SEQRES   1 A    1  GLY\n")
    termini = {"A": {"n": 1, "c": 1}}
    receipt = tmp_path / "receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "state": "accepted",
                "geometry": {"failures": []},
                "declared_termini": termini,
                "output_sha256": m.sha(heavy),
                "source_hashes": {},
            }
        )
    )
    upstream = tmp_path / "status.json"
    upstream.write_text(json.dumps({"exit_status": 0}))
    profile = {
        "root": str(tmp_path),
        "output": str(tmp_path / "result"),
        "source_sha256": {},
        "upstream_status": str(upstream),
        "upstream_exit_key": "exit_status",
        "receipt": str(receipt),
        "heavy": str(heavy),
        "original": str(original),
        "expected_termini": termini,
        "adoption_contract": "explicit_refinement_termini",
        "software": str(software),
        "python": sys.executable,
        "seed": 12345,
        "parameter_model": "unchanged_native_FTPL_diagnostic_comparison",
    }
    profile_path = tmp_path / "profile.json"
    profile_path.write_text(json.dumps(profile))
    monkeypatch.setattr(m.Controller, "journal", lambda *args: None)
    monkeypatch.chdir(tmp_path)

    def prepare(source: Path, folder: Path, digest: str, expected: dict) -> None:
        assert m.sha(source) == digest
        assert expected == termini
        folder.mkdir()
        shutil.copy2(source, folder / "protein.pdb")

    import nyl12_covalent.mcce_preparation as preparation
    import nyl12_covalent.mcce_validation as validation

    monkeypatch.setattr(preparation, "prepare_input", prepare)
    monkeypatch.setattr(
        validation, "audit_preparation", lambda *args: {"accepted": True}
    )

    def write_params(phase: str, args: SimpleNamespace) -> None:
        Path("run.prm").write_text(phase)
        Path(f"{phase}_args.json").write_text(json.dumps(vars(args)))

    monkeypatch.setitem(
        sys.modules,
        "step1",
        SimpleNamespace(
            write_runprm=lambda args: write_params("step1", args),
            identify_ligands=lambda lines: lines,
            fix_format=lambda *args, **kwargs: ["fixture"],
        ),
    )
    monkeypatch.setitem(
        sys.modules,
        "step2",
        SimpleNamespace(write_runprm=lambda args: write_params("step2", args)),
    )
    monkeypatch.setitem(
        sys.modules,
        "step4",
        SimpleNamespace(write_runprm=lambda args: write_params("step4", args)),
    )
    return m.Controller(profile_path), profile, tmp_path


def test_full_fake_native_pipeline_and_parameter_identity(native_case: tuple) -> None:
    """Exercise real timeout/subprocess boundaries and exact native parameter calls."""
    controller, profile, root = native_case
    controller.run()
    result = root / "result"
    assert (
        json.loads((result / "controller_status.json").read_text())["state"]
        == "native_completed_pka_unqualified"
    )
    completion = json.loads((result / "native_completion.json").read_text())
    assert completion["assigned_states"] is False
    assert completion["MS_OUT"] is False
    args = json.loads((result / "prepared/step4_args.json").read_text())
    assert args == {
        "t": "ph",
        "i": "0.0",
        "d": "1.0",
        "n": "15",
        "xts": True,
        "ms": False,
        "u": "MONTE_SEED=12345,MONTE_T=298.15",
        "load_runprm": "",
    }
    pbe = json.loads((result / "step3_process_intent.json").read_text())["command"]
    assert pbe[-12:] == [
        "-s",
        "delphi",
        "-d",
        "4",
        "-do",
        "80",
        "-salt",
        "0.15",
        "-p",
        "2",
        "-t",
        "pbe_tmp",
    ]
    assert json.loads((result / "prepared/step2_args.json").read_text())["l"] == 1
    assert set(completion["artifact_sha256"]) == {
        "head3.lst",
        "fort.38",
        "pK.out",
        "run.prm.record",
    }


def test_optional_refinement_uses_declared_inputs(native_case: tuple) -> None:
    """Run a fake refinement child and retain the precise argument vector."""
    c, p, root = native_case
    fake = executable(root / "refine.py", "pass\n")
    c.profile["refinement"] = {
        "python": sys.executable,
        "script": str(fake),
        "candidate": "declared_candidate",
        "original": "declared_original",
        "output": "fresh_refine",
    }
    c.profile["termini_path"] = "declared_termini.json"
    c.run()
    command = json.loads((root / "result/refinement_process_intent.json").read_text())[
        "command"
    ]
    assert command[-8:] == [
        "--candidate",
        "declared_candidate",
        "--original",
        "declared_original",
        "--outputdir",
        "fresh_refine",
        "--termini",
        "declared_termini.json",
    ]


@pytest.mark.parametrize(
    "body,timeout,expected",
    [
        ("import sys; sys.exit(7)", 5, 7),
        ("import time; time.sleep(30)", 1, 124),
        ('print("ERROR: synthetic native failure")', 5, 0),
    ],
)
def test_failed_native_child_has_actual_exit_record(
    native_case: tuple, body: str, timeout: int, expected: int
) -> None:
    """Failures preserve actual native exits, including GNU timeout and error logs."""
    c, p, root = native_case
    fake = executable(root / "failing.py", body)
    with pytest.raises(RuntimeError, match="failed"):
        c.stage("fake_failure", [str(fake)], timeout)
    outcome = json.loads(
        (root / "result/fake_failure_process_outcome.json").read_text()
    )
    assert outcome["exit_status"] == expected
    assert Path("fake_failure.exit").read_text().strip() == str(expected)


@pytest.mark.parametrize(
    "fault,match",
    [
        ("source", "Source changed"),
        ("receipt", "scientific acceptance"),
        ("digest", "digest mismatch"),
        ("prep", "scientific gate"),
        ("step2", "step2_out"),
        ("pbe_failure", "PBE conformers"),
        ("table", "conformer table"),
        ("pair", "pair-energy"),
    ],
)
def test_pipeline_rejects_failed_output_gates(
    native_case: tuple, monkeypatch: pytest.MonkeyPatch, fault: str, match: str
) -> None:
    """Stop each failed artifact boundary rather than propagating invalid output."""
    c, p, root = native_case
    if fault == "source":
        c.profile["source_sha256"] = {p["heavy"]: "wrong"}
    if fault in {"receipt", "digest"}:
        receipt = json.loads(Path(p["receipt"]).read_text())
        receipt["state" if fault == "receipt" else "output_sha256"] = "invalid"
        Path(p["receipt"]).write_text(json.dumps(receipt))
    if fault == "prep":
        import nyl12_covalent.mcce_validation as validation

        monkeypatch.setattr(
            validation,
            "audit_preparation",
            lambda *args: {"accepted": False, "errors": ["synthetic"]},
        )
    original_stage = c.stage

    def stage(label: str, argv: list[str], seconds: int) -> None:
        original_stage(label, argv, seconds)
        if label == "step2" and fault == "step2":
            Path("step2_out.pdb").unlink()
        if label == "step3":
            if fault == "pbe_failure":
                Path("failed_conformers.txt").write_text("synthetic failed conformer")
            if fault == "table":
                Path("head3.lst").write_text("header\n1 truncated\n")
            if fault == "pair":
                Path("energies/ASP01A0001_001.opp").unlink()

    monkeypatch.setattr(c, "stage", stage)
    with pytest.raises((ValueError, FileNotFoundError), match=match):
        c.run()
