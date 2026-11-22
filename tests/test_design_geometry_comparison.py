"""Tests for design-stage PA6 geometry comparisons."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent import design_geometry_comparison as comparison
from scripts import screen_pa6_binding_orientation as screen

SEED_PDB = (
    screen.ROOT / "qmmm/pa6_static_es_seed/pa6_l4z_l2_nyl12_wt_model_23_canonical.pdb"
)
REFERENCE = screen.ROOT / "cif/n12c.cif"


def test_one_ligand_comparison_writes_four_site_metrics(tmp_path: Path) -> None:
    """Preserve real model23 geometry and select its closest catalytic site."""
    report = comparison.compare_models(
        [
            {
                "variant": "wt",
                "seed": 1,
                "stage": "raw",
                "path": str(SEED_PDB),
                "ligand_count": 1,
            }
        ],
        REFERENCE,
        tmp_path / "comparison",
    )
    pairs = list(csv.DictReader((tmp_path / "comparison/all_sites.csv").open()))
    occupied = list(csv.DictReader((tmp_path / "comparison/occupied_sites.csv").open()))
    summary = list(csv.DictReader((tmp_path / "comparison/model_summary.csv").open()))
    assert len(pairs) == 4
    assert len(occupied) == 1
    assert occupied[0]["site"] == "AE"
    assert float(occupied[0]["herschlag_phi"]) == pytest.approx(99.01, abs=0.02)
    assert occupied[0]["tier"] == "A"
    assert summary[0]["joint_pass_count"] == "1"
    assert report["model_count"] == 1
    assert "nested" in report["replicate_unit"]
    assert (
        json.loads((tmp_path / "comparison/report.json").read_text())["comparison_rows"]
        == 4
    )


def test_four_ligands_use_minimum_sum_one_to_one_assignment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep all 16 raw pairs while selecting the globally minimal assignment."""
    pdb = tmp_path / "synthetic.pdb"
    pdb.write_text("synthetic fixture marker\n")
    atoms = {}
    for index, chain in enumerate("IJKL"):
        for name in comparison.LIGAND_NAMES:
            atoms[(chain, "1", name)] = ("C", np.array([float(index), 0.0, 0.0]))
    monkeypatch.setattr(screen, "read_atoms", lambda _path: atoms)
    monkeypatch.setattr(screen, "layout_site", lambda _atoms, i: {"site_marker": i})

    def fake_measure(site: dict, ligand: dict, _reference: dict) -> dict:
        index = int(ligand["C54"][0])
        target = (index + 1) % 4
        return {
            "d1": float(abs(site["site_marker"] - target)),
            "tier": "A",
            "herschlag_phi": 90.0,
            "native_pass": True,
            "joint_pass": True,
        }

    monkeypatch.setattr(screen, "measure", fake_measure)
    report = comparison.compare_models(
        [
            {
                "variant": "variant_01",
                "seed": 2,
                "stage": "relaxed",
                "path": str(pdb),
                "ligand_count": 4,
            }
        ],
        REFERENCE,
        tmp_path / "four",
    )
    raw = list(csv.DictReader((tmp_path / "four/all_sites.csv").open()))
    occupied = list(csv.DictReader((tmp_path / "four/occupied_sites.csv").open()))
    assert len(raw) == 16
    assert len(occupied) == 4
    assert {row["site"] for row in occupied} == set(comparison.SITE_NAMES)
    assert all(row["is_assigned"] == "True" for row in occupied)
    assert report["occupied_site_count"] == 4
    summary = list(csv.DictReader((tmp_path / "four/model_summary.csv").open()))
    assignment = json.loads(summary[0]["assignment"])
    assert assignment == {"I:1": "BF", "J:1": "CG", "K:1": "DH", "L:1": "AE"}


def test_assignment_policy_handles_single_and_four_ligands() -> None:
    """Check closest-site assignment and the globally optimal permutation."""
    one = {"I:1": dict(zip(comparison.SITE_NAMES, (4.0, 1.0, 3.0, 2.0), strict=True))}
    assert comparison._assign_ligands(["I:1"], one) == {"I:1": "BF"}
    costs = {
        ligand: dict(zip(comparison.SITE_NAMES, row, strict=True))
        for ligand, row in zip(
            ("I:1", "J:1", "K:1", "L:1"),
            ((7, 0, 7, 7), (7, 7, 0, 7), (7, 7, 7, 0), (0, 7, 7, 7)),
            strict=True,
        )
    }
    assignment = comparison._assign_ligands(list(costs), costs)
    assert len(set(assignment.values())) == 4
    assert assignment == {"I:1": "BF", "J:1": "CG", "K:1": "DH", "L:1": "AE"}


def test_explicit_ligand_count_rejects_unexpected_inventory(tmp_path: Path) -> None:
    """Fail when structure ligand inventory disagrees with entry metadata."""
    with pytest.raises(ValueError, match="expected 4 PA6 ligands, found 1"):
        comparison.compare_models(
            [
                {
                    "variant": "wt",
                    "seed": 1,
                    "stage": "raw",
                    "path": str(SEED_PDB),
                    "ligand_count": 4,
                }
            ],
            REFERENCE,
            tmp_path / "bad",
        )
