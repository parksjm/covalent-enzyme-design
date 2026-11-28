"""Focused tests for intact precursor PDBFixer input preparation."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_nylc_precursor_pdbfixer import parse_seqres, prepare_original  # noqa: E402


def _seqres(chain: str, sequence: list[str]) -> list[str]:
    return [
        f"SEQRES {offset // 13 + 1:3d} {chain} {len(sequence):4d}  "
        + " ".join(sequence[offset : offset + 13])
        for offset in range(0, len(sequence), 13)
    ]


def _atom(
    serial: int,
    atom: str,
    altloc: str,
    chain: str,
    residue: int,
    occupancy: float,
    element: str = "C",
) -> str:
    return (
        f"ATOM  {serial:5d} {atom:4s}{altloc:1s}ALA {chain:1s}{residue:4d}    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{occupancy:6.2f}{20.0:6.2f}"
        f"          {element:>2s}  "
    )


def _source_lines() -> list[str]:
    sequence = ["ALA"] * 355
    return [
        *_seqres("A", sequence),
        *_seqres("B", sequence),
        _atom(1, " CA ", " ", "A", 1, 1.0),
        _atom(2, " CB ", " ", "A", 267, 1.0),
        _atom(3, " CB ", "A", "A", 267, 0.4),
        _atom(4, " CB ", "B", "A", 267, 0.6),
        _atom(5, " CA ", " ", "B", 355, 1.0),
        _atom(6, " H  ", " ", "B", 355, 1.0, "H"),
    ]


def test_parse_seqres_keeps_both_intact_precursor_sequences() -> None:
    """Requires both complete, 355-residue precursor chains."""
    sequences = parse_seqres(_source_lines())
    assert list(sequences) == ["A", "B"]
    assert len(sequences["A"]) == len(sequences["B"]) == 355


def test_prepare_original_prefers_altloc_a_and_retains_internal_267(
    tmp_path: Path,
) -> None:
    """Writes both unsplit chains, records altloc choice and excludes hydrogen."""
    source = tmp_path / "source.pdb"
    destination = tmp_path / "original_protein.pdb"
    source.write_text("\n".join(_source_lines()) + "\n")
    sequences, choices = prepare_original(source, destination)
    rows = destination.read_text().splitlines()
    assert len(sequences["A"]) == len(sequences["B"]) == 355
    assert {line[21:22] for line in rows if line.startswith("ATOM  ")} == {"A", "B"}
    assert choices["A:267:CB"] == "A"
    cb = next(
        line
        for line in rows
        if line.startswith("ATOM  ")
        and line[21:22] == "A"
        and line[22:26] == " 267"
        and line[12:16].strip() == "CB"
    )
    assert cb[54:60].strip() == "0.40"
    assert cb[16:17] == " "
    assert not any(
        line[76:78].strip() == "H" for line in rows if line.startswith("ATOM  ")
    )
    assert sum(line.startswith("TER") for line in rows) == 2
    assert any(line.startswith("SEQRES") and line[11:12] == "A" for line in rows)
