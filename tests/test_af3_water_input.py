"""Tests for the standard-library AF3 water-test input builder."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from scripts.launch_nyl12_af3_water_test import prepare_input


def _write_fixture(
    tmp_path: Path, mutation: str | None = None
) -> tuple[Path, Path, dict[str, object]]:
    """Write a minimal valid original JSON payload and CCD file."""
    ccd = tmp_path / "substrate.cif"
    ccd.write_text("data_substrate\n# original CCD\n", encoding="utf-8")
    proteins = [
        {
            "protein": {
                "id": ["A", "B", "C", "D"],
                "sequence": "A" * 232,
                "unpairedMsa": "alpha unpaired",
                "pairedMsa": "alpha paired",
                "templates": [],
            }
        },
        {
            "protein": {
                "id": ["E", "F", "G", "H"],
                "sequence": "T" * 96,
                "unpairedMsa": "beta unpaired",
                "pairedMsa": "beta paired",
                "templates": [],
            }
        },
    ]
    if mutation == "chain":
        proteins[1]["protein"]["id"] = ["X", "F", "G", "H"]
    if mutation == "missing_msa":
        del proteins[1]["protein"]["pairedMsa"]
    payload: dict[str, object] = {
        "name": "original",
        "sequences": [*proteins, {"ligand": {"id": "I", "ccdCodes": ["L4"]}}],
        "bondedAtomPairs": [[["E", 1, "OG1"], ["I", 1, "C54"]]],
        "userCCDPath": str(ccd),
        "modelSeeds": [7],
    }
    if mutation == "graph":
        payload["bondedAtomPairs"] = []
    original = tmp_path / "original.json"
    original.write_text(json.dumps(payload), encoding="utf-8")
    return original, ccd, payload


def test_prepare_input_keeps_required_entities_and_contact_constraints(
    tmp_path: Path,
) -> None:
    """The prepared input has A/E chains, intact covalency, and water contacts."""
    original, ccd, original_payload = _write_fixture(tmp_path)
    output = tmp_path / "attempt"
    output.mkdir()

    input_path = prepare_input(original, output)

    prepared = json.loads(input_path.read_text(encoding="utf-8"))
    proteins = [item["protein"] for item in prepared["sequences"] if "protein" in item]
    assert [protein["id"] for protein in proteins] == ["A", "E"]
    assert [len(protein["sequence"]) for protein in proteins] == [232, 96]
    assert prepared["bondedAtomPairs"] == original_payload["bondedAtomPairs"]
    assert Path(prepared["userCCDPath"]) == output / "tetrahedral_substrate.cif"
    assert (output / "tetrahedral_substrate.cif").read_bytes() == ccd.read_bytes()
    water = [
        item["ligand"]
        for item in prepared["sequences"]
        if item.get("ligand", {}).get("id") == "J"
    ]
    assert water == [{"id": "J", "ccdCodes": ["HOH"]}]
    contacts = prepared["contactConstraints"]
    assert [contact["token1"] for contact in contacts] == [
        ["E", 1, "N"],
        ["E", 1, "OG1"],
    ]
    for contact in contacts:
        assert contact["token2"] == ["J", 1, "O"]
        assert contact["minDistance"] == 2.6
        assert contact["maxDistance"] == 3.2
        assert contact["mobileEndpoint"] == "token2"
    untouched = json.loads(original.read_text(encoding="utf-8"))
    assert untouched["sequences"][0]["protein"]["id"] == ["A", "B", "C", "D"]
    assert untouched["name"] == "original"


@pytest.mark.parametrize("mutation", ["graph", "chain", "missing_msa"])
def test_prepare_input_rejects_invalid_graph_chain_or_msa(
    tmp_path: Path, mutation: str
) -> None:
    """Invalid covalency, chain IDs, or precomputed sequence inputs fail."""
    original, _, _ = _write_fixture(tmp_path, mutation)
    output = tmp_path / "attempt"
    output.mkdir()
    expected_message = {
        "graph": "audited PA6 L4 TI1 covalent graph",
        "chain": "chain identity or length",
        "missing_msa": "Missing precomputed pairedMsa",
    }[mutation]
    with pytest.raises(ValueError, match=expected_message):
        prepare_input(original, output)
    assert not (output / "input.json").exists()


def test_prepare_input_rejects_missing_precomputed_fields(
    tmp_path: Path, field: str = "unpairedMsa"
) -> None:
    """Every required precomputed MSA/template field is checked explicitly."""
    original, _, payload = _write_fixture(tmp_path)
    changed = copy.deepcopy(payload)
    del changed["sequences"][0]["protein"][field]
    original.write_text(json.dumps(changed), encoding="utf-8")
    output = tmp_path / "attempt"
    output.mkdir()
    with pytest.raises(ValueError, match=f"Missing precomputed {field}"):
        prepare_input(original, output)


def _write_exclusion_audit(
    path: Path, rows: list[dict[str, object]]
) -> Path:
    """Write a minimal reference-water exclusion audit."""
    path.write_text(
        json.dumps({"water_neighbors_within_5A": rows}), encoding="utf-8"
    )
    return path


def test_prepare_input_adds_local_exclusions_and_copies_reference_audit(
    tmp_path: Path,
) -> None:
    """Valid reference neighbors become mobile-water restraints with provenance."""
    original, _, _ = _write_fixture(tmp_path)
    output = tmp_path / "attempt"
    output.mkdir()
    audit = _write_exclusion_audit(
        tmp_path / "reference.json",
        [
            {"atom": "E:1:THR:N", "distance_A": 1.2, "minimum_A": 2.5},
            {"atom": "E:1:THR:OG1", "distance_A": 1.4, "minimum_A": 2.75},
            {"atom": "A:40:ALA:O", "distance_A": 3.1, "minimum_A": 2.5},
            {"atom": "I:1:LIG:C54", "distance_A": 3.4, "minimum_A": 2.75},
        ],
    )

    prepared = json.loads(prepare_input(original, output, audit).read_text())

    assert len(prepared["contactConstraints"]) == 4
    assert [item["token1"] for item in prepared["contactConstraints"]] == [
        ["E", 1, "N"],
        ["E", 1, "OG1"],
        ["A", 40, "O"],
        ["I", 1, "C54"],
    ]
    for item in prepared["contactConstraints"][2:]:
        assert item["token2"] == ["J", 1, "O"]
        assert item["maxDistance"] == 1000.0
        assert item["mobileEndpoint"] == "token2"
    assert [item["minDistance"] for item in prepared["contactConstraints"][2:]] == [
        2.5,
        2.75,
    ]
    assert (output / "exclusion_reference_audit.json").read_bytes() == audit.read_bytes()


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ({"atom": "Q:10:ALA:O", "distance_A": 3.0, "minimum_A": 2.5}, "Invalid atom or distance"),
        ({"atom": "A:10:ALA:O", "distance_A": 5.1, "minimum_A": 2.5}, "Invalid atom or distance"),
        ({"atom": "A:10:ALA:O", "distance_A": 3.0, "minimum_A": 3.0}, "Unexpected heavy-atom exclusion threshold"),
    ],
)
def test_prepare_input_rejects_invalid_exclusion_rows(
    tmp_path: Path, row: dict[str, object], message: str
) -> None:
    """Out-of-scope atoms, distances, and thresholds fail explicitly."""
    original, _, _ = _write_fixture(tmp_path)
    output = tmp_path / "attempt"
    output.mkdir()
    audit = _write_exclusion_audit(tmp_path / "reference.json", [row])

    with pytest.raises(ValueError, match=message):
        prepare_input(original, output, audit)

    assert not (output / "input.json").exists()


def test_prepare_input_rejects_exclusion_set_over_capacity(
    tmp_path: Path,
) -> None:
    """More than 64 total restraints fail before writing an inference input."""
    original, _, _ = _write_fixture(tmp_path)
    output = tmp_path / "attempt"
    output.mkdir()
    rows = [
        {"atom": f"A:{residue}:ALA:O", "distance_A": 3.0, "minimum_A": 2.5}
        for residue in range(1, 64)
    ]
    audit = _write_exclusion_audit(tmp_path / "reference.json", rows)

    with pytest.raises(ValueError, match="64-restraint capacity"):
        prepare_input(original, output, audit)

    assert not (output / "input.json").exists()
