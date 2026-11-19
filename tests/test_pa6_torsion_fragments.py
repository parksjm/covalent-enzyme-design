"""Test deterministic PA6 torsion benchmark fragment preparation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.pa6_torsion_fragments import FRAGMENTS, prepare_torsion_fragments


def test_prepare_fragments_writes_declared_chemistry(tmp_path: Path) -> None:
    output = tmp_path / "fragments"
    manifest = prepare_torsion_fragments(output)

    assert [record["label"] for record in manifest["fragments"]] == [
        "internal_neutral",
        "n_terminus_cation",
        "c_terminus_anion",
    ]
    assert [record["charge"] for record in manifest["fragments"]] == [0, 1, -1]
    assert json.loads((output / "manifest.json").read_text()) == manifest
    for spec, record in zip(FRAGMENTS, manifest["fragments"], strict=True):
        fragment = output / spec.label
        assert (fragment / "seed_mmff94s.sdf").is_file()
        assert (fragment / "seed_mmff94s.xyz").is_file()
        gaussian = (fragment / "opt_smd.gjf").read_text()
        assert "wB97XD/def2SVP" in gaussian
        assert "SCRF=(SMD,Solvent=Water)" in gaussian
        assert "CalcFC" not in gaussian
        assert f"\n{spec.charge} 1\n" in gaussian
        assert record["proper_torsions_gaussian_one_based"]


def test_fragment_preparation_refuses_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "fragments"
    output.mkdir()
    with pytest.raises(FileExistsError):
        prepare_torsion_fragments(output)
