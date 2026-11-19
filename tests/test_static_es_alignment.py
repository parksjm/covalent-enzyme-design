"""Tests for canonical static E-S catalytic alignment."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent import static_es_alignment as static_alignment


def _pdb_atom(serial: int, atom: str, residue: int, xyz: tuple[float, ...]) -> str:
    element = atom[0]
    return (
        f"ATOM  {serial:5d} {atom:>4s} THR A{residue:4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}"
        f"          {element:>2s}\n"
    )


def test_read_canonical_pdb_site_selects_heavy_atoms(tmp_path: Path) -> None:
    pdb = tmp_path / "model.pdb"
    pdb.write_text(
        _pdb_atom(1, "N", 233, (1.0, 2.0, 3.0))
        + _pdb_atom(2, "H", 233, (1.0, 2.0, 4.0))
        + _pdb_atom(3, "OG1", 233, (2.0, 2.0, 3.0)),
        encoding="utf-8",
    )

    atoms = static_alignment.read_canonical_pdb_site(pdb, "A", {233})

    assert set(atoms) == {(233, "N"), (233, "OG1")}
    assert atoms[(233, "OG1")][1] == pytest.approx([2.0, 2.0, 3.0])


def test_alignment_writes_report_and_refuses_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = {
        (113, "N"): ("N", np.array([0.0, 0.0, 0.0])),
        (113, "CA"): ("C", np.array([1.0, 0.0, 0.0])),
        (113, "C"): ("C", np.array([1.0, 1.0, 0.0])),
        (113, "O"): ("O", np.array([1.0, 1.0, 1.0])),
        (186, "N"): ("N", np.array([0.0, 3.0, 0.0])),
        (186, "CA"): ("C", np.array([1.0, 3.0, 0.0])),
        (186, "C"): ("C", np.array([1.0, 4.0, 0.0])),
        (186, "O"): ("O", np.array([1.0, 4.0, 1.0])),
        (186, "ND2"): ("N", np.array([2.0, 3.0, 0.0])),
        (233, "N"): ("N", np.array([4.0, 0.0, 0.0])),
        (233, "CA"): ("C", np.array([5.0, 0.0, 0.0])),
        (233, "C"): ("C", np.array([5.0, 1.0, 0.0])),
        (233, "O"): ("O", np.array([5.0, 2.0, 0.0])),
        (233, "CB"): ("C", np.array([5.0, 0.0, 1.0])),
        (233, "OG1"): ("O", np.array([5.0, 1.0, 1.0])),
        (233, "CG2"): ("C", np.array([6.0, 0.0, 1.0])),
        (234, "N"): ("N", np.array([5.0, 1.0, 1.0])),
        (272, "N"): ("N", np.array([0.0, 6.0, 0.0])),
        (272, "CA"): ("C", np.array([1.0, 6.0, 0.0])),
        (272, "C"): ("C", np.array([1.0, 7.0, 0.0])),
        (272, "O"): ("O", np.array([1.0, 7.0, 1.0])),
        (272, "OD1"): ("O", np.array([2.0, 6.0, 0.0])),
        (272, "OD2"): ("O", np.array([2.0, 7.0, 0.0])),
        (274, "N"): ("N", np.array([0.0, 9.0, 0.0])),
        (274, "CA"): ("C", np.array([1.0, 9.0, 0.0])),
        (274, "C"): ("C", np.array([1.0, 10.0, 0.0])),
        (274, "O"): ("O", np.array([1.0, 10.0, 1.0])),
        (274, "OD1"): ("O", np.array([2.0, 9.0, 0.0])),
        (274, "OD2"): ("O", np.array([2.0, 10.0, 0.0])),
    }
    model = {
        key: (element, coordinate + np.array([8.0, -3.0, 2.0]))
        for key, (element, coordinate) in reference.items()
    }
    monkeypatch.setattr(static_alignment, "_cif_atoms", lambda *args: reference)
    monkeypatch.setattr(
        static_alignment,
        "read_canonical_pdb_site",
        lambda *args: model,
    )
    output = tmp_path / "output"

    report = static_alignment.align_static_es_model(
        tmp_path / "reference.cif", tmp_path / "model.pdb", output
    )

    assert report["fit"]["rmsd_angstrom"] == pytest.approx(0.0)
    assert report["thr233_held_out"]["rmsd_angstrom"] == pytest.approx(0.0)
    assert (output / "alignment_report.json").is_file()
    assert (output / "catalytic_site_overlay.pdb").is_file()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        static_alignment.align_static_es_model(
            tmp_path / "reference.cif", tmp_path / "model.pdb", output
        )
