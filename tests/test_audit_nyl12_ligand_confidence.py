"""Focused tests for the Nyl12 AF3 ligand confidence audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent.chemistry import build_poly_xama_maxa
from scripts.audit_nyl12_ligand_confidence import (
    _distribution,
    _single_match,
    audit_model,
)


def _write_model(directory: Path, target: str = "poly_xama_maxa") -> tuple[Path, list[str]]:
    """Write a compact, schema-faithful AF3 model fixture."""
    directory.mkdir()
    if target == "poly_xama_maxa":
        molecule = build_poly_xama_maxa().mol
        names = [atom.GetProp("name") for atom in molecule.GetAtoms()
                 if atom.GetAtomicNum() > 1]
        substrate = build_poly_xama_maxa()
        reactive = {substrate.carbonyl_c: 83.84, substrate.carbonyl_o: 54.43,
                    substrate.leaving_n: 60.16}
        atom_values: dict[str, float] = {}
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            name = atom.GetProp("name")
            position = atom.GetIntProp("monad_index")
            atom_values[name] = reactive.get(name, 20.0 if position in (1, 8) else 65.0)
    else:
        names = ["C54", "O49", "N59"]
        atom_values = {"C54": 71.0, "O49": 62.0, "N59": 58.0}
    chains = ["A", *(["I"] * len(names))]
    values = [45.0, *(atom_values[name] for name in names)]
    atom_headers = ["_atom_site.label_asym_id", "_atom_site.label_atom_id",
                    "_atom_site.B_iso_or_equiv"]
    cif_rows = ["loop_", *atom_headers,
                "A PRO 45.00",
                *(f"I {name} {value:.2f}" for name, value in zip(names, values[1:], strict=True)), "#"]
    (directory / "fixture_model.cif").write_text("\n".join(cif_rows) + "\n")
    confidence = {
        "atom_chain_ids": chains,
        "atom_plddts": values,
        "token_chain_ids": ["A", "I"],
        "pae": [[0.8, 4.0], [7.0, 0.8]],
    }
    (directory / "fixture_confidences.json").write_text(json.dumps(confidence))
    summary: dict[str, Any] = {
        "ranking_score": 0.9135986597,
        "iptm": 0.7,
        "ptm": 0.8,
        "actifptm": 0.75,
        "chain_pair_iptm": [[0.8, 0.6], [0.6, 0.3]],
        "chain_pair_pae_min": [[0.8, 4.0], [4.0, 0.8]],
    }
    (directory / "fixture_summary_confidences.json").write_text(json.dumps(summary))
    with (directory / "fixture_ranking_scores.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("seed", "sample", "ranking_score"))
        writer.writeheader()
        writer.writerow({"seed": 1, "sample": 0, "ranking_score": summary["ranking_score"]})
    return directory, names


def test_distribution_reports_counts_and_bins() -> None:
    result = _distribution([20.0, 50.0, 69.0, 70.0, 90.0])

    assert result["count"] == 5
    assert result["mean"] == pytest.approx(59.8)
    assert result["median"] == 69.0
    assert {name: item["count"] for name, item in result["bins"].items()} == {
        "lt50": 1, "50to70": 2, "ge70": 2
    }


def test_single_match_excludes_summary_confidence_file(tmp_path: Path) -> None:
    regular = tmp_path / "model_confidences.json"
    summary = tmp_path / "model_summary_confidences.json"
    regular.touch()
    summary.touch()

    assert _single_match(tmp_path, "_confidences.json") == regular


def test_poly_audit_reports_named_atoms_monads_pae_and_hashes(tmp_path: Path) -> None:
    model_dir, names = _write_model(tmp_path / "poly")

    result = audit_model(model_dir, "poly_xama_maxa")

    assert result["ligand_atom_confidence"]["count"] == len(names)
    assert result["reactive_atoms"] == {"C033": 83.84, "O034": 54.43, "N035": 60.16}
    assert result["monad_confidence"]["terminal_monads_1_and_8"]["mean"] == 20.0
    assert result["monad_confidence"]["internal_monads_2_through_7"]["mean"] > 60
    assert result["interchain_pae"]["ligand_to_protein"]["mean"] == 7.0
    assert result["interchain_pae"]["protein_to_ligand"]["mean"] == 4.0
    assert len(result["sha256"]["inputs"]) == 4
    assert "ranking_score" in result["summary_metrics"]


def test_pa6_audit_includes_measured_reactive_atom_names(tmp_path: Path) -> None:
    model_dir, _ = _write_model(tmp_path / "pa6", "pa6_l2")

    result = audit_model(model_dir, "pa6_l2")

    assert result["reactive_atom_names"] == ["C54", "O49", "N59"]
    assert set(result["reactive_atoms"]) == {"C54", "O49", "N59"}
    assert "monad_confidence" not in result


def test_audit_accepts_summary_ranking_score_rounded_to_two_decimals(
    tmp_path: Path,
) -> None:
    model_dir, _ = _write_model(tmp_path / "poly")
    summary_path = model_dir / "fixture_summary_confidences.json"
    summary = json.loads(summary_path.read_text())
    summary["ranking_score"] = 0.91
    summary_path.write_text(json.dumps(summary))

    result = audit_model(model_dir, "poly_xama_maxa")

    assert result["summary_metrics"]["ranking_score"] == pytest.approx(0.9135986597)
    assert result["summary_metrics"]["ranking_score_summary"] == 0.91
    assert "rounding to two decimals" in result["summary_metrics"]["ranking_score_precision_note"]


def test_audit_rejects_mismatched_csv_ranking_score(tmp_path: Path) -> None:
    model_dir, _ = _write_model(tmp_path / "poly")
    with (model_dir / "fixture_ranking_scores.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("seed", "sample", "ranking_score"))
        writer.writeheader()
        writer.writerow({"seed": 1, "sample": 0, "ranking_score": 0.1})

    with pytest.raises(ValueError, match="ranking scores disagree"):
        audit_model(model_dir, "poly_xama_maxa")


def test_audit_rejects_invalid_ligand_plddt(tmp_path: Path) -> None:
    model_dir, _ = _write_model(tmp_path / "poly")
    path = model_dir / "fixture_confidences.json"
    confidence = json.loads(path.read_text())
    confidence["atom_plddts"][1] = 120
    path.write_text(json.dumps(confidence))
    cif = model_dir / "fixture_model.cif"
    text = cif.read_text().replace("A PRO 45.00", "A PRO 45.00")
    cif.write_text(text.replace("I N001 20.00", "I N001 100.00"))

    with pytest.raises(ValueError, match="Invalid pLDDT"):
        audit_model(model_dir, "poly_xama_maxa")
