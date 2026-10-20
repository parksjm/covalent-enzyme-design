"""Focused tests for provenance-safe relay-water contact reports."""

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import diagnose_nyl12_relay_water_contacts as diagnostic  # noqa: E402
from diagnose_nyl12_relay_water_contacts import (  # noqa: E402
    _atom_label,
    _verify_search_inputs,
    _write_new_report,
    diagnose,
)


def _atom(topology_residue: int) -> SimpleNamespace:
    residue = SimpleNamespace(
        idx=topology_residue - 1,
        number=topology_residue - 1,
        name="ALA",
        chain="",
    )
    return SimpleNamespace(
        idx=0,
        name="CA",
        atomic_number=6,
        residue=residue,
        bond_partners=[],
    )


@pytest.mark.parametrize(
    ("topology_residue", "chain", "local", "canonical"),
    [(186, "A", 186, 186), (929, "E", 1, 233)],
)
def test_atom_labels_include_canonical_and_chain_local_numbers(
    topology_residue: int, chain: str, local: int, canonical: int
) -> None:
    label = _atom_label(_atom(topology_residue))
    assert label["canonical_chain"] == chain
    assert label["chain_residue_number"] == local
    assert label["canonical_residue_number"] == canonical


def test_search_input_hashes_must_match(tmp_path: Path) -> None:
    topology = tmp_path / "top.prmtop"
    coordinates = tmp_path / "coords.inpcrd"
    topology.write_text("topology")
    coordinates.write_text("coordinates")
    search = {
        "input_sha256": {
            str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (topology, coordinates)
        }
    }
    _verify_search_inputs(search, topology, coordinates)
    coordinates.write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        _verify_search_inputs(search, topology, coordinates)


def test_report_writer_refuses_overwrite(tmp_path: Path) -> None:
    destination = tmp_path / "report.json"
    _write_new_report(destination, {"status": "new"})
    assert json.loads(destination.read_text()) == {"status": "new"}
    with pytest.raises(FileExistsError):
        _write_new_report(destination, {"status": "replacement"})


def _synthetic_topology() -> SimpleNamespace:
    residues = [SimpleNamespace(idx=i, number=i, name="ALA", chain="") for i in range(1314)]
    atoms = []

    def add_atom(name: str, residue_idx: int, element: int, xyz: tuple[float, float, float], partners=()) -> None:
        atom = SimpleNamespace(
            idx=len(atoms), name=name, residue=residues[residue_idx],
            atomic_number=element, bond_partners=list(partners),
        )
        atoms.append(atom)
        residues[residue_idx].atoms = getattr(residues[residue_idx], "atoms", []) + [atom]
        points.append(xyz)

    points: list[tuple[float, float, float]] = []
    add_atom("N", 928, 7, (10.0, 0.0, 0.0))
    add_atom("OG1", 928, 8, (12.0, 0.0, 0.0))
    carbon = SimpleNamespace(name="N", atomic_number=7)
    add_atom("H1", 185, 1, (10.0, 0.0, 0.0), [carbon])
    add_atom("CA", 185, 6, (10.5, 0.0, 0.0))
    return SimpleNamespace(residues=residues, atoms=atoms, coordinates=np.asarray(points))


def test_diagnose_reports_sorted_threshold_contacts_and_nh_parent() -> None:
    topology = _synthetic_topology()
    candidates = {
        "method": "synthetic",
        "water_atom_indices_zero_based": [4, 5, 6],
        "candidates": [{
            "accepted": False,
            "maximum_contact_shortfall_A": 0.75,
            "water_xyz_A": [[10.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
        }],
    }
    report = diagnose(topology, candidates, {"input": "sha256"})
    candidate = report["candidates"][0]
    assert report["local_environment_atom_count"] == 4
    assert candidate["offending_contact_count"] == 10
    assert candidate["max_shortfall_A"] == pytest.approx(2.5)
    assert candidate["n_h_contacts"] == 3
    assert any(
        contact["canonical_chain"] == "A"
        and contact["canonical_residue_number"] == 186
        and contact["hydrogen_parents"] == [{"name": "N", "element": 7}]
        for contact in candidate["contacts"]
    )


def test_cli_writes_contact_report_and_provenance(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    topology_path = tmp_path / "system.prmtop"
    coordinate_path = tmp_path / "system.inpcrd"
    search_path = tmp_path / "search.json"
    output_path = tmp_path / "report.json"
    topology_path.write_text("topology")
    coordinate_path.write_text("coordinates")
    search = {
        "input_sha256": {
            str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (topology_path, coordinate_path)
        },
        "method": "cli synthetic",
        "water_atom_indices_zero_based": [4, 5, 6],
        "candidates": [],
    }
    search_path.write_text(json.dumps(search))
    monkeypatch.setattr(diagnostic, "AmberParm", lambda *args, **kwargs: _synthetic_topology())
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "diagnostic", "--topology", str(topology_path), "--coordinates", str(coordinate_path),
            "--search-json", str(search_path), "--output", str(output_path),
        ],
    )

    diagnostic.main()
    report = json.loads(output_path.read_text())
    assert report["source_method"] == "cli synthetic"
    assert report["candidates"] == []
    assert report["input_sha256"][str(search_path.resolve())]
