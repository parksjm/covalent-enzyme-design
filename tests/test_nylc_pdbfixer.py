"""Focused tests for the 5XYG preprocessing policy."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from repair_nylc_5xyg_pdbfixer import (
    parse_seqres,
    select_protein_atoms,
    write_split_input,
)  # noqa: E402


def _seqres_lines(chain_id: str, sequence: list[str]) -> list[str]:
    """Build fixed-column SEQRES lines for a synthetic deposited sequence."""
    result = []
    for offset in range(0, len(sequence), 13):
        serial = offset // 13 + 1
        result.append(
            f"SEQRES {serial:3d} {chain_id} {len(sequence):4d}  "
            + " ".join(sequence[offset : offset + 13])
        )
    return result


def _atom_line(
    serial: int,
    atom_name: str,
    altloc: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    occupancy: float,
    element: str = "C",
) -> str:
    """Build one synthetic ATOM line with PDB fixed columns."""
    return (
        f"ATOM  {serial:5d} {atom_name:4s}{altloc:1s}{residue_name:3s} "
        f"{chain_id:1s}{residue_number:4d}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"{occupancy:6.2f}{20.0:6.2f}          {element:>2s}  "
    )


def _synthetic_lines() -> list[str]:
    sequence = ["ALA"] * 355
    lines = _seqres_lines("A", sequence) + _seqres_lines("B", sequence)
    lines.extend(
        [
            _atom_line(1, " CA ", " ", "ALA", "A", 1, 1.0),
            _atom_line(2, " CB ", "A", "ALA", "A", 266, 0.40),
            _atom_line(3, " CB ", "B", "ALA", "A", 266, 0.60),
            _atom_line(4, " CA ", " ", "ALA", "A", 267, 1.0),
            _atom_line(5, " CA ", " ", "ALA", "B", 1, 1.0),
            _atom_line(6, " CA ", " ", "ALA", "B", 355, 1.0),
        ]
    )
    return lines


def test_parse_seqres_requires_matching_355_residue_chains() -> None:
    """Rejects mismatched or incomplete deposited sequence metadata."""
    sequences = parse_seqres(_synthetic_lines())
    assert len(sequences["A"]) == 355
    bad = _synthetic_lines()
    bad[0] = bad[0].replace("ALA", "GLY", 1)
    with pytest.raises(ValueError, match="identical SEQRES"):
        parse_seqres(bad)


def test_select_protein_atoms_splits_cleavage_and_prefers_altloc_a() -> None:
    """Splits 266/267 and deterministically selects conformer A."""
    selected, choices = select_protein_atoms(_synthetic_lines())
    chain_ids = {line[21:22] for lines in selected.values() for line in lines}
    assert chain_ids == {"A", "B", "C", "D"}
    assert {int(line[22:26]) for line in selected["A"]} == {1, 266}
    assert {int(line[22:26]) for line in selected["C"]} == {267}
    chosen_cb = next(line for line in selected["A"] if line[12:16].strip() == "CB")
    assert chosen_cb[54:60].strip() == "0.40"
    assert chosen_cb[16:17] == " "
    assert choices["A:266:CB"].startswith("altloc=A")


def test_write_split_input_emits_four_complete_split_sequences(tmp_path: Path) -> None:
    """Writes the mature chain plan and preserves original author numbering."""
    source = tmp_path / "source.pdb"
    destination = tmp_path / "split.pdb"
    source.write_text("\n".join(_synthetic_lines()) + "\n")
    metadata = write_split_input(source, destination)
    seqres: dict[str, list[str]] = {chain: [] for chain in "ABCD"}
    atoms: dict[str, set[int]] = {chain: set() for chain in "ABCD"}
    for line in destination.read_text().splitlines():
        if line.startswith("SEQRES"):
            seqres[line[11:12]].extend(line[19:70].split())
        elif line.startswith("ATOM  "):
            atoms[line[21:22]].add(int(line[22:26]))
    assert {chain: len(seq) for chain, seq in seqres.items()} == {
        "A": 266,
        "B": 266,
        "C": 89,
        "D": 89,
    }
    assert atoms["A"] == {1, 266}
    assert atoms["C"] == {267}
    assert metadata["input_atom_counts"] == {"A": 2, "B": 1, "C": 1, "D": 1}
    assert "A:266:CB" in metadata["selected_altlocs"]
