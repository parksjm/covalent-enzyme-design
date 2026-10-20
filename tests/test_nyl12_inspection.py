"""Tests for canonical full-system inspection PDB export."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.inspection import write_inspection_pdb


def _topology() -> SimpleNamespace:
    residues = [SimpleNamespace(name="ALA", atoms=[]) for _ in range(1314)]
    residues[928].name = "THR"
    residues[1312].name = "LIG"
    residues[1313].name = "WAT"
    atom = lambda name, residue, element: SimpleNamespace(  # noqa: E731
        name=name, residue=residue, element_name=element
    )
    og1 = atom("OG1", residues[928], "O")
    carbon = atom("C1", residues[1312], "C")
    hydrogen = atom("H1", residues[1312], "H")
    water_o = atom("O", residues[1313], "O")
    residues[928].atoms = [og1]
    residues[1312].atoms = [carbon, hydrogen]
    residues[1313].atoms = [water_o]
    for residue in residues:
        for item in residue.atoms:
            item.residue = residue
    return SimpleNamespace(
        residues=residues,
        atoms=[og1, carbon, hydrogen, water_o],
        bonds=[
            SimpleNamespace(atom1=og1, atom2=carbon),
            SimpleNamespace(atom1=carbon, atom2=hydrogen),
        ],
    )


def test_writes_explicit_atoms_mapping_chain_ters_and_covalent_bonds(tmp_path: Path) -> None:
    topology = _topology()
    destination = tmp_path / "nested" / "inspection.pdb"
    mapping = write_inspection_pdb(topology, np.arange(12).reshape(4, 3), destination)
    text = destination.read_text(encoding="ascii")

    assert mapping[1] == ("A", 1)
    assert mapping[232] == ("A", 232)
    assert mapping[233] == ("B", 1)
    assert mapping[929] == ("E", 1)
    assert mapping[1313] == ("I", 1)
    assert mapping[1314] == ("J", 1)
    assert " OG1 THR E   1" in text
    assert "  C1 LIG I   1" in text
    assert "  H1 LIG I   1" in text
    assert "TER" in text
    assert "CONECT    1    2" in text
    assert "CONECT    2    3" in text


def test_rejects_unexpected_residue_count_and_coordinate_shape(tmp_path: Path) -> None:
    topology = _topology()
    with pytest.raises(ValueError, match="1314 residues"):
        write_inspection_pdb(SimpleNamespace(residues=topology.residues[:-1], atoms=topology.atoms), np.zeros((4, 3)), tmp_path / "x.pdb")
    with pytest.raises(ValueError, match="coordinates must have shape"):
        write_inspection_pdb(topology, np.zeros((3, 3)), tmp_path / "x.pdb")
