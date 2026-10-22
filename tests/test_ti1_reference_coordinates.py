"""Focused tests for crystal-derived TI1 torsion references."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import ti1_reference_coordinates as reference


COORDINATES = {
    "N": (0.0, 0.0, 0.0),
    "CA": (0.1, 1.45, 0.0),
    "C": (1.45, 1.4, 0.4),
    "O": (1.6, 1.1, 1.6),
    "CB": (-0.6, 2.0, 0.2),
    "OG1": (-1.6, 1.9, 0.1),
    "CG2": (-0.5, 3.3, 0.3),
}


def _atom(name: str, xyz: tuple[float, float, float]) -> SimpleNamespace:
    return SimpleNamespace(
        name=name, xx=xyz[0], xy=xyz[1], xz=xyz[2], other_locations={}
    )


def _residue(
    number: int, name: str, atoms: dict[str, tuple[float, float, float]]
) -> SimpleNamespace:
    return SimpleNamespace(
        chain="A",
        number=number,
        name=name,
        atoms=[_atom(atom_name, xyz) for atom_name, xyz in atoms.items()],
    )


def _structure(
    *,
    current: dict[str, tuple[float, float, float]] | None = None,
    residues: list[SimpleNamespace] | None = None,
    model_count: int = 1,
) -> SimpleNamespace:
    if residues is None:
        current = COORDINATES if current is None else current
        previous = _residue(232, "GLY", {"C": (-1.2, 0.0, 0.0)})
        target = _residue(233, "THR", current)
        following = _residue(234, "GLY", {"N": (2.75, 1.6, 0.4)})
        residues = [previous, target, following]
    atom_count = sum(len(residue.atoms) for residue in residues)
    return SimpleNamespace(
        residues=residues,
        get_coordinates=lambda: np.zeros((model_count, atom_count, 3), dtype=float),
    )


def _parse_mock(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, structure: SimpleNamespace
) -> Path:
    source = tmp_path / "reference.cif"
    source.write_text("synthetic mmCIF source", encoding="utf-8")
    monkeypatch.setattr(reference.CIFFile, "parse", lambda *_args, **_kwargs: structure)
    return source


def test_extracts_threonine_angles_and_source_identity(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    source = _parse_mock(monkeypatch, tmp_path, _structure())
    result = reference.crystal_backbone_reference(source, "A", 233)

    assert result["chain"] == "A"
    assert result["residue_number"] == 233
    assert len(result["sha256"]) == 64
    assert result["phi_degrees"] is not None
    assert np.isfinite(
        [result["phi_degrees"], result["psi_degrees"], result["chi1_degrees"]]
    ).all()
    assert result["previous_C_to_N_A"] == pytest.approx(1.2)
    assert result["C_to_next_N_A"] == pytest.approx(
        np.linalg.norm(np.subtract((2.75, 1.6, 0.4), COORDINATES["C"]))
    )
    assert result["thr_coordinates_A"]["OG1"] == pytest.approx(COORDINATES["OG1"])


def test_missing_previous_residue_leaves_phi_unset(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    source = _parse_mock(
        monkeypatch,
        tmp_path,
        _structure(
            residues=[
                _residue(233, "THR", COORDINATES),
                _residue(234, "GLY", {"N": (2.75, 1.6, 0.4)}),
            ]
        ),
    )
    result = reference.crystal_backbone_reference(source, "A", 233)
    assert result["phi_degrees"] is None


@pytest.mark.parametrize(
    "mutate,chain,number",
    [
        (lambda residues: residues, "Z", 233),
        (lambda residues: residues, "A", 999),
        (lambda residues: [*residues, _residue(233, "THR", COORDINATES)], "A", 233),
        (
            lambda residues: [
                residues[0],
                _residue(233, "SER", COORDINATES),
                residues[2],
            ],
            "A",
            233,
        ),
        (
            lambda residues: [
                residues[0],
                _residue(
                    233, "THR", {k: v for k, v in COORDINATES.items() if k != "OG1"}
                ),
                residues[2],
            ],
            "A",
            233,
        ),
    ],
)
def test_rejects_invalid_selector_or_threonine(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    mutate: object,
    chain: str,
    number: int,
) -> None:
    residues = _structure().residues
    changed = mutate(residues)  # type: ignore[operator]
    source = _parse_mock(monkeypatch, tmp_path, _structure(residues=changed))
    with pytest.raises((ValueError, LookupError)):
        reference.crystal_backbone_reference(source, chain, number)


def test_rejects_duplicate_atom_and_altloc(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    residues = _structure().residues
    duplicate_atoms = [*residues[1].atoms, _atom("CA", COORDINATES["CA"])]
    duplicate = SimpleNamespace(
        chain="A", number=233, name="THR", atoms=duplicate_atoms
    )
    source = _parse_mock(
        monkeypatch,
        tmp_path,
        _structure(residues=[residues[0], duplicate, residues[2]]),
    )
    with pytest.raises(ValueError):
        reference.crystal_backbone_reference(source, "A", 233)

    residues = _structure().residues
    residues[1].atoms[0].other_locations = {"A": object()}
    source = _parse_mock(monkeypatch, tmp_path, _structure(residues=residues))
    with pytest.raises(ValueError):
        reference.crystal_backbone_reference(source, "A", 233)


def test_rejects_multimodel_and_nonfinite_coordinates(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    source = _parse_mock(monkeypatch, tmp_path, _structure(model_count=2))
    with pytest.raises(ValueError):
        reference.crystal_backbone_reference(source, "A", 233)

    nonfinite = dict(COORDINATES)
    nonfinite["CA"] = (float("nan"), 1.45, 0.0)
    source = _parse_mock(monkeypatch, tmp_path, _structure(current=nonfinite))
    with pytest.raises(ValueError):
        reference.crystal_backbone_reference(source, "A", 233)
