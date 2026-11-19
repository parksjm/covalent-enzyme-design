"""Tests for Nyl12 active-site geometry and Amber atom selection."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.active_site_geometry import (
    coordinates_from_amber,
    distance,
    measure_active_site,
    projected_herschlag_phi,
    source_sha256,
)


def _coordinates() -> dict[str, tuple[float, float, float]]:
    return {
        "Thr929.OG1": (1.0, 0.0, 1.0),
        "Thr929.CB": (1.0, 0.0, 2.0),
        "Phe113.N": (0.0, 1.0, 1.0),
        "Asn186.ND2": (1.0, 1.0, 0.0),
        "substrate.C54": (0.0, 0.0, 1.0),
        "substrate.O49": (0.0, 1.0, 0.0),
        "substrate.N59": (0.0, 2.0, 0.0),
        "substrate.C62": (1.0, 2.0, 0.0),
        "substrate.C033": (0.0, 0.0, 1.0),
        "substrate.O034": (0.0, 1.0, 0.0),
        "substrate.N035": (0.0, 2.0, 0.0),
        "substrate.C032": (1.0, 2.0, 0.0),
    }


def test_distance_and_measure_report_metrics_and_provenance(tmp_path: Path) -> None:
    source = tmp_path / "reference.py"
    source.write_bytes(b"reference math\n")
    report = measure_active_site(_coordinates(), "pa6", "ti1", source)
    assert distance((0, 0, 0), (1, 0, 0)) == pytest.approx(1.0)
    metrics = report["metrics"]
    assert metrics["d1_thr929_og1_scissile_c"] == pytest.approx(1.0)
    assert metrics["d2_phe113_n_scissile_o"] == pytest.approx(1.0)
    assert metrics["d3_asn186_nd2_scissile_o"] == pytest.approx(1.0)
    assert metrics["a_BD_og1_c_o"] == pytest.approx(90.0)
    assert metrics["theta_cb_og1_c"] == pytest.approx(90.0)
    assert metrics["herschlag_phi_projected"] > 0.0
    assert (
        report["geometry_source_sha256"]
        == hashlib.sha256(source.read_bytes()).hexdigest()
    )
    assert "do not classify" in report["state_interpretation"]


def test_poly_uses_its_own_named_atoms_and_es_state(tmp_path: Path) -> None:
    source = tmp_path / "reference.py"
    source.write_text("math", encoding="utf-8")
    report = measure_active_site(_coordinates(), "poly", "es", source)
    assert report["state"] == "es"
    assert report["atom_names"]["amide_n"] == "N035"
    assert report["metrics"]["d1_thr929_og1_scissile_c"] == pytest.approx(1.0)


@pytest.mark.parametrize("state", ["", "transition", "ES"])
def test_state_must_be_explicit_and_supported(state: str, tmp_path: Path) -> None:
    source = tmp_path / "reference.py"
    source.write_text("math", encoding="utf-8")
    with pytest.raises(ValueError, match="state"):
        measure_active_site(_coordinates(), "pa6", state, source)


def test_unknown_substrate_and_missing_coordinates_fail(tmp_path: Path) -> None:
    source = tmp_path / "reference.py"
    source.write_text("math", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported substrate"):
        measure_active_site(_coordinates(), "other", "es", source)
    coordinates = _coordinates()
    del coordinates["substrate.O49"]
    with pytest.raises(ValueError, match="Missing required"):
        measure_active_site(coordinates, "pa6", "es", source)


@pytest.mark.parametrize(
    "bad",
    [(0.0, 0.0, float("nan")), (1.0, 2.0)],
)
def test_invalid_atom_coordinates_fail(bad: tuple[float, ...], tmp_path: Path) -> None:
    source = tmp_path / "reference.py"
    source.write_text("math", encoding="utf-8")
    coordinates = _coordinates()
    coordinates["Thr929.OG1"] = bad  # type: ignore[assignment]
    with pytest.raises(ValueError, match="finite|three"):
        measure_active_site(coordinates, "pa6", "es", source)


def test_projected_phi_has_signed_convention() -> None:
    positive = projected_herschlag_phi(
        (1, 0, 1), (0, 1, 0), (0, 0, 1), (0, 2, 0), (1, 2, 0)
    )
    negative = projected_herschlag_phi(
        (1, 0, -1), (0, 1, 0), (0, 0, 1), (0, 2, 0), (1, 2, 0)
    )
    assert positive > 0
    assert negative < 0


@pytest.mark.parametrize(
    "points",
    [
        ((1, 0, 0), (0, 0, 0), (0, 0, 1), (1, 0, 0), (2, 0, 0)),
        ((0, 2, 0), (0, 1, 0), (0, 0, 1), (0, 2, 0), (1, 2, 0)),
    ],
)
def test_projected_phi_rejects_degenerate_plane_or_dihedral(points: tuple) -> None:
    with pytest.raises(ValueError, match="Degenerate"):
        projected_herschlag_phi(*points)


def test_source_sha256_reads_exact_bytes(tmp_path: Path) -> None:
    source = tmp_path / "source.py"
    source.write_bytes(b"provenance")
    assert source_sha256(source) == hashlib.sha256(b"provenance").hexdigest()


@dataclass
class _Atom:
    name: str
    idx: int


@dataclass
class _Residue:
    name: str
    atoms: list[_Atom]


def _amber_model() -> SimpleNamespace:
    residues = [_Residue("ALA", []) for _ in range(1314)]
    xyz: list[list[float]] = []
    index = 0
    for residue_index, residue_name, names in (
        (112, "PHE", ["N"]),
        (185, "ASN", ["ND2"]),
        (928, "THR", ["OG1", "CB"]),
        (1312, "LIG", ["C54", "O49", "N59", "C62", "C033", "O034", "N035", "C032"]),
    ):
        residue_atoms = []
        for atom_name in names:
            residue_atoms.append(_Atom(atom_name, index))
            xyz.append([float(index), 0.0, 0.0])
            index += 1
        residues[residue_index] = _Residue(residue_name, residue_atoms)
    atoms = [atom for residue in residues for atom in residue.atoms]
    return SimpleNamespace(residues=residues, atoms=atoms, coordinates=np.asarray(xyz))


def test_coordinates_from_amber_maps_required_residues_and_target_atoms() -> None:
    coordinates = coordinates_from_amber(_amber_model(), "pa6")
    assert set(coordinates) == {
        "Thr929.OG1",
        "Thr929.CB",
        "Phe113.N",
        "Asn186.ND2",
        "substrate.C54",
        "substrate.O49",
        "substrate.N59",
        "substrate.C62",
    }


def test_coordinates_from_amber_accepts_appended_solvent_residues() -> None:
    model = _amber_model()
    water_oxygen = _Atom("O", len(model.atoms))
    model.residues.append(_Residue("WAT", [water_oxygen]))
    model.atoms.append(water_oxygen)
    model.coordinates = np.vstack((model.coordinates, [[99.0, 99.0, 99.0]]))
    coordinates = coordinates_from_amber(model, "pa6")
    assert coordinates["substrate.C54"] == (4.0, 0.0, 0.0)


def test_coordinates_from_amber_rejects_identity_missing_duplicates_and_shape() -> None:
    model = _amber_model()
    model.residues[928].name = "SER"
    with pytest.raises(ValueError, match="must be THR"):
        coordinates_from_amber(model, "pa6")
    model = _amber_model()
    model.residues[928].atoms.append(_Atom("OG1", 0))
    with pytest.raises(ValueError, match="exactly one OG1"):
        coordinates_from_amber(model, "pa6")
    model = _amber_model()
    model.coordinates = model.coordinates[:-1]
    with pytest.raises(ValueError, match="coordinates must have shape"):
        coordinates_from_amber(model, "pa6")
    model = _amber_model()
    model.coordinates[0, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        coordinates_from_amber(model, "pa6")


def test_coordinates_from_amber_rejects_residue_count_and_target() -> None:
    model = _amber_model()
    model.residues.pop()
    with pytest.raises(ValueError, match="at least 1314 Amber residues"):
        coordinates_from_amber(model, "pa6")
    with pytest.raises(ValueError, match="Unsupported substrate"):
        coordinates_from_amber(_amber_model(), "unknown")
