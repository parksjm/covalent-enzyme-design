from pathlib import Path

import pytest

from nyl12_covalent.es_system_preparation import (
    _ion_counts,
    off_atom_names,
    write_leap_pdb,
)


def test_off_atom_names_reads_lig_table(tmp_path: Path) -> None:
    library = tmp_path / "lig.lib"
    library.write_text(
        '!entry.LIG.unit.atoms table str name str type\n"C1" "c3"\n'
        '"H1" "hc"\n!entry.LIG.unit.positions table\n'
    )
    assert off_atom_names(library) == {"C1", "H1"}


def test_ion_counts_neutralizes_and_adds_salt() -> None:
    sodium, chloride, pairs = _ion_counts(-8.0, 1_000_000.0)
    assert sodium == pairs + 8
    assert chloride == pairs


def test_ion_counts_rejects_nonintegral_charge() -> None:
    with pytest.raises(ValueError, match="Nonintegral"):
        _ion_counts(-8.2, 1_000_000.0)


def test_write_leap_pdb_rejects_missing_chain_terminators(tmp_path: Path) -> None:
    source = tmp_path / "source.pdb"
    source.write_text(
        "HETATM    1  C1  LIG A1313       0.000   0.000   0.000  1.00  0.00           C  \n"
    )
    with pytest.raises(ValueError, match="chain terminators"):
        write_leap_pdb(source, tmp_path / "out.pdb", {"C1"})
