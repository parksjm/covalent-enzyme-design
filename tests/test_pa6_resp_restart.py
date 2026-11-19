"""Test explicit, non-overwriting PA6 RESP geometry restart preparation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.pa6_resp_preparation import prepare_pa6_resp_conformers
from nyl12_covalent.pa6_resp_restart import GaussianSettings, prepare_pa6_resp_restart

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "provenance/xa_ma_capped_reference_20261022/pa6_l4_source.mol2"
MODEL = (
    ROOT
    / "models/nyl12_template_conditioned_es/02_amber_rescue_relaxation/model_7/LIG_0001.pdb"
)
TORSIONS = ROOT / "inputs/pa6_l4_resp_gaff2/frozen_backbone_torsions.json"


def _source_manifest(tmp_path: Path) -> Path:
    source_dir = tmp_path / "source"
    prepare_pa6_resp_conformers(
        SOURCE,
        {"template_model7": MODEL},
        TORSIONS,
        source_dir,
    )
    return source_dir / "manifest.json"


def _settings() -> GaussianSettings:
    return GaussianSettings(
        method="wB97XD",
        basis="def2SVP",
        solvent="SMD:Water",
        max_cycles=300,
        cores=8,
        memory_gb=16,
        calc_fc=False,
    )


def test_prepares_standalone_def2svp_no_calcfc_attempt(tmp_path: Path) -> None:
    output = tmp_path / "restart"
    manifest = prepare_pa6_resp_restart(_source_manifest(tmp_path), output, _settings())

    assert json.loads((output / "manifest.json").read_text()) == manifest
    assert manifest["gaussian"]["basis"] == "def2SVP"
    assert manifest["gaussian"]["calc_fc"] is False
    record = manifest["conformers"][0]
    assert (output / record["coordinates_xyz"]).is_file()
    gaussian = (output / record["gaussian_input"]).read_text()
    assert "wB97XD/def2SVP" in gaussian
    assert "CalcFC" not in gaussian
    assert "Opt=(Tight,ModRedundant,MaxCycles=300)" in gaussian
    assert "\n\nD 1 2 3 4 F" in gaussian
    assert gaussian.count(" F\n") == 25


def test_rejects_existing_output_without_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "restart"
    output.mkdir()
    with pytest.raises(FileExistsError, match="already exists"):
        prepare_pa6_resp_restart(_source_manifest(tmp_path), output, _settings())


@pytest.mark.parametrize(
    "settings",
    [
        GaussianSettings,
    ],
)
def test_settings_require_explicit_smd_syntax(settings: type[GaussianSettings]) -> None:
    with pytest.raises(ValueError, match="SMD"):
        settings("wB97XD", "def2SVP", "Water", 300, 8, 16, False)
