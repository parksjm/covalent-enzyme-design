"""Focused tests for catalytic scaffold alignment reporting."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from nyl12_covalent import catalytic_alignment as alignment
from nyl12_covalent import catalytic_alignment_geometry as geometry


def _reference_site() -> geometry.AtomMap:
    coordinates = {
        (113, "N"): (0.0, 0.0, 0.0),
        (113, "CA"): (1.0, 0.0, 0.0),
        (113, "C"): (1.0, 1.0, 0.0),
        (113, "O"): (1.0, 1.0, 1.0),
        (186, "N"): (0.0, 3.0, 0.0),
        (186, "CA"): (1.0, 3.0, 0.0),
        (186, "C"): (1.0, 4.0, 0.0),
        (186, "O"): (1.0, 4.0, 1.0),
        (186, "ND2"): (2.0, 3.0, 0.0),
        (233, "N"): (4.0, 0.0, 0.0),
        (233, "CA"): (5.0, 0.0, 0.0),
        (233, "C"): (5.0, 1.0, 0.0),
        (233, "O"): (5.0, 2.0, 0.0),
        (233, "CB"): (5.0, 0.0, 1.0),
        (233, "OG1"): (5.0, 1.0, 1.0),
        (233, "CG2"): (6.0, 0.0, 1.0),
        (234, "N"): (5.0, 1.0, 1.0),
        (272, "N"): (0.0, 6.0, 0.0),
        (272, "CA"): (1.0, 6.0, 0.0),
        (272, "C"): (1.0, 7.0, 0.0),
        (272, "O"): (1.0, 7.0, 1.0),
        (272, "OD1"): (2.0, 6.0, 0.0),
        (272, "OD2"): (2.0, 7.0, 0.0),
        (274, "N"): (0.0, 9.0, 0.0),
        (274, "CA"): (1.0, 9.0, 0.0),
        (274, "C"): (1.0, 10.0, 0.0),
        (274, "O"): (1.0, 10.0, 1.0),
        (274, "OD1"): (2.0, 9.0, 0.0),
        (274, "OD2"): (2.0, 10.0, 0.0),
    }
    return {
        key: (
            "N" if atom.startswith("N") else "O" if atom.startswith("O") else "C",
            np.array(value),
        )
        for key, value in coordinates.items()
        for atom in (key[1],)
    }


def test_kabsch_fit_removes_rigid_transform() -> None:
    reference = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    mobile = reference @ rotation + np.array([4.0, -2.0, 8.0])
    found_rotation, mobile_center, reference_center, rmsd = geometry.fit_coordinates(
        mobile, reference
    )
    transformed = (mobile - mobile_center) @ found_rotation + reference_center
    assert rmsd == pytest.approx(0.0, abs=1e-12)
    assert transformed == pytest.approx(reference)


def test_dihedral_delta_wraps_across_boundary() -> None:
    assert geometry.wrapped_delta(-179.0, 179.0) == pytest.approx(2.0)


def test_thr233_psi_uses_next_backbone_nitrogen(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(geometry, "CONTACTS", ())
    reference_coordinates = {
        (233, "N"): np.array([0.0, 0.0, 0.0]),
        (233, "CA"): np.array([1.0, 0.0, 0.0]),
        (233, "C"): np.array([1.0, 1.0, 0.0]),
        (233, "CB"): np.array([1.0, 0.0, 1.0]),
        (233, "OG1"): np.array([1.0, 1.0, 1.0]),
        (234, "N"): np.array([1.0, 1.0, 1.0]),
    }
    reference = {
        key: ("N" if key[1] == "N" else "C", value)
        for key, value in reference_coordinates.items()
    }
    model = {
        key: (element, value.copy()) for key, (element, value) in reference.items()
    }
    model[(234, "N")] = ("N", np.array([2.0, 1.0, 0.0]))

    metrics = geometry.state_metrics(reference, model)["thr233_held_out"]

    assert metrics["psi_delta_es_minus_reference_degrees"] != pytest.approx(0.0)


def test_topology_mapping_rejects_wrong_residue() -> None:
    residues = [SimpleNamespace(name="ALA", atoms=[]) for _ in range(113)]
    residues[112] = SimpleNamespace(name="GLY", atoms=[])
    structure = SimpleNamespace(
        residues=residues, atoms=[], coordinates=np.empty((0, 3))
    )
    with pytest.raises(ValueError, match="expected PHE"):
        alignment._topology_atoms(structure, ((113, 113, "PHE"),))


def test_cif_and_topology_atom_readers(tmp_path: Path) -> None:
    cif = tmp_path / "reference.cif"
    cif.write_text(
        """data_reference
loop_
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_atom_id
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
A 233 THR OG1 O 1.0 2.0 3.0 1
#
""",
        encoding="utf-8",
    )
    selected = alignment._cif_atoms(cif, "A", {233})
    assert selected[(233, "OG1")][0] == "O"
    assert selected[(233, "OG1")][1] == pytest.approx([1.0, 2.0, 3.0])

    atom = SimpleNamespace(name="OG1", atomic_number=8, element=8, idx=0)
    structure = SimpleNamespace(
        residues=[SimpleNamespace(name="THR", atoms=[atom])],
        atoms=[atom],
        coordinates=np.array([[1.0, 2.0, 3.0]]),
    )
    mapped = alignment._topology_atoms(structure, ((1, 233, "THR"),))
    assert mapped[(233, "OG1")][0] == "O"
    assert mapped[(233, "OG1")][1] == pytest.approx([1.0, 2.0, 3.0])


def test_existing_output_is_not_overwritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = tmp_path / "reference.cif"
    reference.write_text("data_missing_loop\n", encoding="utf-8")
    output = tmp_path / "out"
    output.mkdir()
    existing = output / "alignment_report.json"
    existing.write_text("preserve", encoding="utf-8")
    with pytest.raises(ValueError, match="No atom_site loop"):
        alignment.align_catalytic_site(
            reference,
            tmp_path / "top.prmtop",
            tmp_path / "start.ncrst",
            tmp_path / "end.ncrst",
            output,
        )
    assert existing.read_text() == "preserve"


def test_alignment_writes_start_and_endpoint_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = _reference_site()
    translation = np.array([10.0, -2.0, 3.0])
    mobile = {
        key: (element, coordinate + translation)
        for key, (element, coordinate) in reference.items()
    }

    def fake_load_file(topology: str, coordinates: str) -> Any:
        return SimpleNamespace(source=coordinates)

    def fake_topology_atoms(
        structure: Any, residues: tuple[tuple[int, int, str], ...]
    ) -> geometry.AtomMap:
        del structure, residues
        return {
            key: (element, coordinate.copy())
            for key, (element, coordinate) in mobile.items()
        }

    monkeypatch.setattr(alignment, "_cif_atoms", lambda *args: reference)
    monkeypatch.setattr(alignment.parmed, "load_file", fake_load_file)
    monkeypatch.setattr(alignment, "_topology_atoms", fake_topology_atoms)
    output = tmp_path / "alignment"

    report = alignment.align_catalytic_site(
        tmp_path / "reference.cif",
        tmp_path / "top.prmtop",
        tmp_path / "start.ncrst",
        tmp_path / "end.ncrst",
        output,
    )

    assert report["states"]["starting"]["fit_rmsd_angstrom"] == pytest.approx(0.0)
    assert report["states"]["continuation_endpoint"][
        "fit_rmsd_angstrom"
    ] == pytest.approx(0.0)
    assert (output / "alignment_report.json").is_file()
    overlay = (output / "catalytic_site_overlay.pdb").read_text(encoding="utf-8")
    assert "Chain A is n12c.cif reference" in overlay
    assert " B 233" in overlay
    assert " C 233" in overlay
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        alignment.align_catalytic_site(
            tmp_path / "reference.cif",
            tmp_path / "top.prmtop",
            tmp_path / "start.ncrst",
            tmp_path / "end.ncrst",
            output,
        )
