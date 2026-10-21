"""Tests for the explicit-input relay-water repair CLI."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import parmed as pmd
import pytest

from scripts import repair_nyl12_relay_water as cli


def _model() -> pmd.Structure:
    """Build a real ParmEd topology containing the indexed repair residues."""
    structure = pmd.Structure()
    residues: dict[int, list[pmd.Atom]] = {}
    positions: list[tuple[float, float, float]] = []
    for residue_index in range(1314):
        names = {
            185: ("CA", "CB", "CG"),
            928: ("N", "H1", "H2", "H3"),
            1313: ("O", "H1", "H2"),
        }.get(residue_index, ("CA",))
        residue_atoms: list[pmd.Atom] = []
        for name in names:
            element = 7 if name == "N" else 8 if name == "O" else 1 if name.startswith("H") else 6
            atom = pmd.Atom(name=name, atomic_number=element)
            structure.add_atom(atom, "ASN" if residue_index == 185 else "THR" if residue_index == 928 else "WAT" if residue_index == 1313 else "ALA", residue_index + 1)
            residue_atoms.append(atom)
            positions.append((float(residue_index * 10 + len(residue_atoms)), 0.0, 0.0))
        residues[residue_index] = residue_atoms

    def bond(first: pmd.Atom, second: pmd.Atom) -> None:
        structure.bonds.append(pmd.Bond(first, second))

    asn_ca, asn_cb, asn_cg = residues[185]
    bond(asn_ca, asn_cb)
    bond(asn_cb, asn_cg)
    thr_n, *thr_hydrogens = residues[928]
    for hydrogen in thr_hydrogens:
        bond(thr_n, hydrogen)
    water_o, water_h1, water_h2 = residues[1313]
    bond(water_o, water_h1)
    bond(water_o, water_h2)
    structure.coordinates = np.asarray(positions, dtype=float)
    return structure


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    """Write explicit placeholder input files used by the mocked pipeline."""
    topology = tmp_path / "model.parm7"
    coordinates = tmp_path / "model.rst7"
    topology.write_bytes(b"synthetic topology")
    coordinates.write_bytes(b"synthetic coordinates")
    return topology, coordinates


def _patch_pipeline(
    monkeypatch: pytest.MonkeyPatch,
    model: pmd.Structure,
    candidate: np.ndarray | None,
    *,
    preflights: list[dict[str, Any]] | None = None,
) -> None:
    """Mock chemistry-dependent IO while retaining real topology invariance."""
    monkeypatch.setattr(cli, "AmberParm", lambda *args, **kwargs: model)
    monkeypatch.setattr(cli, "check_geometry", lambda *args: (preflights or [{"accepted": True}]).pop(0))
    monkeypatch.setattr(cli, "coordinates_from_amber", lambda loaded, kind: {"coordinates": loaded.coordinates})
    monkeypatch.setattr(
        cli,
        "measure_active_site",
        lambda coordinates, kind, state: {"target": kind, "state": state},
    )
    monkeypatch.setattr(cli, "search_repair", lambda *args, **kwargs: ({"evaluated": 1}, candidate))

    def save(path: str, **kwargs: Any) -> None:
        Path(path).write_bytes(b"synthetic restart")

    monkeypatch.setattr(model, "save", save, raising=False)
    monkeypatch.setattr(
        cli,
        "write_inspection_pdb",
        lambda loaded, coordinates, path: (Path(path).write_text("pdb\n"), {"WAT": 1313})[1],
    )


def test_check_invariance_accepts_water_motion_and_reports_indices() -> None:
    """Rigid water motion is allowed and fixed atoms remain unchanged."""
    model = _model()
    before = np.asarray(model.coordinates).copy()
    after = before.copy()
    water_indices = [atom.idx for atom in model.residues[1313].atoms]
    after[water_indices] += np.array([1.0, 0.0, 0.0])

    result = cli.check_invariance(model, before, after)

    assert result["accepted"] is True
    assert result["changed_atom_indices_zero_based"] == water_indices
    assert set(water_indices).issubset(result["allowed_atom_indices_zero_based"])
    assert result["fixed_atom_maximum_coordinate_change_A"] == 0.0
    assert result["nonwater_maximum_bond_length_change_A"] == 0.0


def test_check_invariance_rejects_fixed_atom_and_asn_bond_changes() -> None:
    """Fixed-coordinate edits and stretched permitted Asn bonds both fail."""
    model = _model()
    before = np.asarray(model.coordinates).copy()
    fixed_atom = model.residues[0].atoms[0].idx
    changed = before.copy()
    changed[fixed_atom, 0] += 0.01
    assert cli.check_invariance(model, before, changed)["accepted"] is False

    changed = before.copy()
    cb = model.residues[185].atoms[1].idx
    changed[cb, 0] += 0.2
    result = cli.check_invariance(model, before, changed)
    assert result["fixed_atom_maximum_coordinate_change_A"] == 0.0
    assert result["nonwater_maximum_bond_length_change_A"] > 1e-8
    assert result["accepted"] is False


def test_repair_files_exports_accepted_seed_and_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Accepted preflight and invariance produce restart and inspection exports."""
    topology, coordinates = _inputs(tmp_path)
    model = _model()
    original = np.asarray(model.coordinates).copy()
    candidate = original.copy()
    water_indices = [atom.idx for atom in model.residues[1313].atoms]
    candidate[water_indices] += [1.0, 0.0, 0.0]
    _patch_pipeline(monkeypatch, model, candidate, preflights=[{"accepted": True}, {"accepted": True}])
    output = tmp_path / "accepted"

    result = cli.repair_files(topology, coordinates, "pa6_l4", output, 15)

    assert result["accepted"] is True
    assert (output / "candidate.rst7").read_bytes() == b"synthetic restart"
    assert (output / "pa6_l4_ti1_starting.pdb").read_text() == "pdb\n"
    assert json.loads((output / "search.json").read_text()) == {"evaluated": 1}
    provenance = json.loads((output / "pa6_l4_ti1_provenance.json").read_text())
    assert provenance["inputs"]["topology"]["sha256"] == hashlib.sha256(topology.read_bytes()).hexdigest()
    assert provenance["inputs"]["coordinates"]["sha256"] == hashlib.sha256(coordinates.read_bytes()).hexdigest()
    assert provenance["archival_target_alias"] == "pa6_l2"
    assert json.loads((output / "pa6_l4_residue_map.json").read_text()) == {"WAT": 1313}
    assert json.loads((output / "acceptance.json").read_text())["accepted"] is True


def test_repair_files_reports_no_candidate_without_export(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A completed search without a seed is recorded as an unaccepted result."""
    topology, coordinates = _inputs(tmp_path)
    _patch_pipeline(monkeypatch, _model(), None)
    output = tmp_path / "no-candidate"

    result = cli.repair_files(topology, coordinates, "maxa_maxa", output, 30)

    assert result == {"accepted": False, "candidate_found": False}
    assert not (output / "candidate.rst7").exists()
    assert json.loads((output / "acceptance.json").read_text()) == result


def test_repair_files_rejects_invariance_violation_and_keeps_acceptance_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A candidate that edits a fixed atom is rejected before restart export."""
    topology, coordinates = _inputs(tmp_path)
    model = _model()
    candidate = np.asarray(model.coordinates).copy()
    candidate[model.residues[0].atoms[0].idx, 0] += 1.0
    _patch_pipeline(monkeypatch, model, candidate)
    output = tmp_path / "invalid-invariance"

    with pytest.raises(ValueError, match="altered fixed atoms"):
        cli.repair_files(topology, coordinates, "pa6_l4", output, 15)

    assert json.loads((output / "acceptance.json").read_text())["accepted"] is False
    assert "invariance" in json.loads((output / "acceptance.json").read_text())
    assert not (output / "candidate.rst7").exists()


def test_repair_files_does_not_export_after_preflight_rejection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A geometrically rejected saved seed never gets the inspection PDB."""
    topology, coordinates = _inputs(tmp_path)
    model = _model()
    candidate = np.asarray(model.coordinates).copy()
    water_indices = [atom.idx for atom in model.residues[1313].atoms]
    candidate[water_indices] += [1.0, 0.0, 0.0]
    _patch_pipeline(
        monkeypatch,
        model,
        candidate,
        preflights=[{"accepted": True}, {"accepted": False}],
    )
    output = tmp_path / "preflight-rejected"

    result = cli.repair_files(topology, coordinates, "pa6_l4", output, 15)

    assert result["accepted"] is False
    assert (output / "candidate.rst7").exists()
    assert not list(output.glob("*_starting.pdb"))


def test_repair_files_refuses_existing_output_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An existing result directory is protected before input loading."""
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("preserve")

    with pytest.raises(FileExistsError):
        cli.repair_files(Path("missing-topology"), Path("missing-coordinates"), "pa6_l4", output, 15)

    assert marker.read_text() == "preserve"


def test_write_json_is_exclusive_and_rejects_nonfinite_values(tmp_path: Path) -> None:
    """JSON artifacts cannot overwrite and reject nonstandard NaN values."""
    path = tmp_path / "report.json"
    cli.write_json(path, {"ok": True})
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        cli.write_json(path, {"ok": False})
    assert path.read_bytes() == original
    with pytest.raises(ValueError, match="Out of range float values"):
        cli.write_json(tmp_path / "nan.json", {"metric": float("nan")})
