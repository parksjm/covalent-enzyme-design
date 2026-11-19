"""Tests for template-conditioned Boltz pilot input preparation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from nyl12_covalent.template_boltz import (
    QUERY_CHAINS,
    TEMPLATE_CHAINS,
    build_template_yaml,
    prepare_pilot,
    validate_source_yaml,
)


def _source(path: Path) -> Path:
    path.write_text(
        "version: 1\n"
        "sequences:\n"
        "- protein:\n"
        "    id: [A, B, C, D]\n"
        "    sequence: AAA\n"
        "- ligand:\n"
        "    id: I\n"
        "    smiles: '[NH3+]CCCCCC'\n"
        "constraints:\n"
        "- contact:\n"
        "    token1: [E, 1]\n"
        "- contact:\n"
        "    token1: [A, 186]\n",
        encoding="utf-8",
    )
    return path


def test_build_template_yaml_preserves_contacts_and_adds_mapping(
    tmp_path: Path,
) -> None:
    source = _source(tmp_path / "source.yaml")
    template = tmp_path / "n12c.cif"
    template.write_text("data_test\n", encoding="utf-8")
    text = build_template_yaml(source, template)
    assert text.count("- ligand:") == 1
    assert text.count("- contact:") == 2
    assert f"chain_id: [{', '.join(QUERY_CHAINS)}]" in text
    assert f"template_id: [{', '.join(TEMPLATE_CHAINS)}]" in text
    assert str(template.resolve()) in text


@pytest.mark.parametrize(
    "text, message",
    [
        ("version: 1\n", "one ligand"),
        ("- ligand:\n  smiles: '[NH3+]CCCCCC'\n", "two contacts"),
        (
            "- ligand:\n  smiles: '[NH3+]CCCCCC'\n- contact:\n- contact:\ntemplates:\n",
            "already contains",
        ),
    ],
)
def test_validate_source_yaml_fails_explicitly(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        validate_source_yaml(text)


def test_prepare_pilot_writes_twelve_seeds_and_provenance(tmp_path: Path) -> None:
    source = _source(tmp_path / "source.yaml")
    template = tmp_path / "n12c.cif"
    template.write_text("data_test\n", encoding="utf-8")
    output = tmp_path / "pilot"
    report = prepare_pilot(source, template, output, 100)
    with (output / "seeds.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 12
    assert [int(row["seed"]) for row in rows] == list(range(100, 112))
    assert {row["campaign"] for row in rows} == {"template_conditioned_es"}
    assert {row["sequence_name"] for row in rows} == {"nyl12_pa6_l4_to_l2_template"}
    assert report["calculation_launched"] is False
    saved = json.loads((output / "provenance.json").read_text(encoding="utf-8"))
    assert saved == report
    assert len(str(report["generated_yaml_sha256"])) == 64
