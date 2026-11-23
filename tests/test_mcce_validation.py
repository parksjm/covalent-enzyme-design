"""Tests for MCCE step1 process and scientific acceptance auditing."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.mcce_validation import audit_preparation


def _pdb_line(
    record: str,
    serial: int,
    atom: str,
    residue: str,
    chain: str,
    residue_number: int,
    xyz: tuple[float, float, float],
    conformer: int = 0,
    native: bool = True,
) -> str:
    """Format one standard or MCCE native fixed-column atom record."""
    prefix = f"{record:<6}{serial:5d} {atom:4s} {residue:3s} {chain:1s}"
    if native:
        prefix += f"{residue_number:04d}_{conformer:03d}"
    else:
        prefix += f"{residue_number:4d}    "
    return f"{prefix}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}\n"


def _fixture(
    tmp_path: Path,
    native_records: list[str] | None = None,
    log: str = "step1 complete\n",
    exit_code: str = "0\n",
) -> tuple[Path, Path]:
    """Create a two-chain input and corresponding minimal MCCE artifacts."""
    input_path = tmp_path / "input.pdb"
    originals = [
        ("N", "ALA", "A", 1, (0.0, 0.0, 0.0)),
        ("CA", "ALA", "A", 1, (1.0, 0.0, 0.0)),
        ("OXT", "GLY", "B", 9, (0.0, 1.0, 0.0)),
    ]
    input_path.write_text(
        "".join(
            _pdb_line("ATOM", index, atom, residue, chain, number, xyz, native=False)
            for index, (atom, residue, chain, number, xyz) in enumerate(originals, 1)
        ),
        encoding="ascii",
    )
    folder = tmp_path / "run"
    folder.mkdir()
    if native_records is None:
        native_records = [
            _pdb_line("ATOM", 1, "N", "ALA", "A", 1, (0, 0, 0), 0),
            _pdb_line("ATOM", 2, "N", "ALA", "A", 1, (0, 0, 0), 1),
            _pdb_line("ATOM", 3, "CA", "ALA", "A", 1, (1, 0, 0), 0),
            _pdb_line("ATOM", 4, "OXT", "GLY", "B", 9, (0, 1, 0), 0),
            _pdb_line("ATOM", 5, "H", "ALA", "A", 1, (0, 0, 1), 0),
            _pdb_line("ATOM", 6, "N", "NTR", "A", 1, (0, 0, 0), 0),
            _pdb_line("ATOM", 7, "OXT", "NTR", "B", 9, (0, 1, 0), 0),
            _pdb_line("ATOM", 8, "N", "CTR", "A", 1, (0, 0, 0), 0),
            _pdb_line("ATOM", 9, "OXT", "CTR", "B", 9, (0, 1, 0), 0),
        ]
    (folder / "step1_out.pdb").write_text("".join(native_records), encoding="ascii")
    (folder / "step1.log").write_text(log, encoding="utf-8")
    (folder / "step1.exit").write_text(exit_code, encoding="ascii")
    return folder, input_path


@pytest.fixture
def termini() -> dict[str, dict[str, int]]:
    """Return explicit endpoints for the two test chains."""
    return {"A": {"n": 1, "c": 1}, "B": {"n": 9, "c": 9}}


def test_accepts_exact_termini_and_duplicate_conformer_coordinates(
    tmp_path: Path, termini: dict[str, dict[str, int]]
) -> None:
    """Conformer duplicates are accepted when an original position survives."""
    folder, input_path = _fixture(tmp_path)
    report = audit_preparation(folder, input_path, termini)
    assert report["accepted"] is True, report
    assert report["criteria"]["max_minimum_coordinate_displacement_angstrom"] == 0
    assert report["counts"]["native_heavy_atom_identities"] == 3
    assert len(report["observed_NTR"]) == 2
    assert len(report["observed_CTR"]) == 2
    assert len(report["artifact_sha256"]) == 4


@pytest.mark.parametrize(
    "terminal,records",
    [("NTR", []), ("NTR", [_pdb_line("ATOM", 8, "N", "NTR", "Z", 12, (0, 0, 0))])],
)
def test_missing_or_spurious_n_terminus_fails(
    tmp_path: Path,
    termini: dict[str, dict[str, int]],
    terminal: str,
    records: list[str],
) -> None:
    """Missing and extra N-terminal pseudo-residues fail exact matching."""
    folder, input_path = _fixture(tmp_path)
    output = folder / "step1_out.pdb"
    lines = output.read_text(encoding="ascii").splitlines(keepends=True)
    if not records:
        lines = [line for line in lines if line[17:20].strip() != terminal]
    else:
        lines.extend(records)
    output.write_text("".join(lines), encoding="ascii")
    report = audit_preparation(folder, input_path, termini)
    assert report["accepted"] is False
    assert "n_termini_match" in report["errors"]


def test_error_log_fails_even_with_zero_exit(
    tmp_path: Path, termini: dict[str, dict[str, int]]
) -> None:
    """Diagnostic errors in the log veto a zero process exit code."""
    folder, input_path = _fixture(tmp_path, log="Error! bad parameter\n")
    report = audit_preparation(folder, input_path, termini)
    assert report["criteria"]["exit_code_zero"] is True
    assert report["accepted"] is False
    assert report["criteria"]["error_log_lines"] == ["Error! bad parameter"]


@pytest.mark.parametrize("mutation", ["lost_oxt", "moved", "extra"])
def test_heavy_atom_loss_displacement_and_addition_fail(
    tmp_path: Path, termini: dict[str, dict[str, int]], mutation: str
) -> None:
    """Original atom identity/position and exact native atom inventories gate."""
    folder, input_path = _fixture(tmp_path)
    path = folder / "step1_out.pdb"
    lines = path.read_text(encoding="ascii").splitlines(keepends=True)
    if mutation == "lost_oxt":
        lines = [line for line in lines if line[12:16].strip() != "OXT"]
    elif mutation == "moved":
        lines = [
            _pdb_line("ATOM", 3, "CA", "ALA", "A", 1, (1.01, 0, 0))
            if line[12:16].strip() == "CA"
            else line
            for line in lines
        ]
    else:
        lines.append(_pdb_line("ATOM", 99, "C2", "ALA", "A", 1, (3, 3, 3)))
    path.write_text("".join(lines), encoding="ascii")
    report = audit_preparation(folder, input_path, termini)
    assert report["accepted"] is False
    if mutation == "extra":
        assert report["counts"]["extra_native_heavy_atom_identities"] == 1
    else:
        assert report["counts"]["missing_or_moved_heavy_atoms"] == 1


def test_nonzero_exit_fails(tmp_path: Path, termini: dict[str, dict[str, int]]) -> None:
    """A nonzero step1 exit code vetoes otherwise valid scientific data."""
    folder, input_path = _fixture(tmp_path, exit_code="2\n")
    report = audit_preparation(folder, input_path, termini)
    assert report["accepted"] is False
    assert "exit_code_zero" in report["errors"]


def test_malformed_exit_fails_acceptance_without_raising(
    tmp_path: Path, termini: dict[str, dict[str, int]]
) -> None:
    """A malformed process marker is a failed gate, not a parser exception."""
    folder, input_path = _fixture(tmp_path, exit_code="finished\n")
    report = audit_preparation(folder, input_path, termini)
    assert report["accepted"] is False
    assert "exit_code_is_integer" in report["errors"]


def test_native_extra_hydrogen_chain_fails(
    tmp_path: Path, termini: dict[str, dict[str, int]]
) -> None:
    """Chain counts include hydrogen-only native ATOM records."""
    folder, input_path = _fixture(tmp_path)
    output = folder / "step1_out.pdb"
    with output.open("a", encoding="ascii") as stream:
        stream.write(_pdb_line("ATOM", 10, "H", "ALA", "Z", 4, (0, 0, 0)))
    report = audit_preparation(folder, input_path, termini)
    assert report["accepted"] is False
    assert "native_chain_count_matches" in report["errors"]


def test_malformed_pdb_raises(
    tmp_path: Path, termini: dict[str, dict[str, int]]
) -> None:
    """Malformed fixed-column input is a parser error, not a failed gate."""
    folder, input_path = _fixture(tmp_path)
    input_path.write_text("ATOM\n", encoding="ascii")
    with pytest.raises(ValueError, match="ATOM record"):
        audit_preparation(folder, input_path, termini)


def test_cli_writes_exclusive_json_report(
    tmp_path: Path, termini: dict[str, dict[str, int]]
) -> None:
    """The command interface persists a JSON report at a fresh path."""
    from scripts.audit_mcce_preparation import main

    folder, input_path = _fixture(tmp_path)
    termini_path = tmp_path / "termini.json"
    termini_path.write_text(json.dumps(termini), encoding="utf-8")
    output = tmp_path / "report.json"
    result = main(
        [
            "--folder",
            str(folder),
            "--input",
            str(input_path),
            "--termini",
            str(termini_path),
            "--output",
            str(output),
        ]
    )
    assert result == 0
    assert json.loads(output.read_text(encoding="utf-8"))["accepted"] is True
    with pytest.raises(FileExistsError):
        main(
            [
                "--folder",
                str(folder),
                "--input",
                str(input_path),
                "--termini",
                str(termini_path),
                "--output",
                str(output),
            ]
        )
