"""Test strict PA6 RESP conformer preparation without launching Gaussian."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.pa6_resp_preparation import prepare_pa6_resp_conformers
from nyl12_covalent.pa6_resp_restart import (
    GaussianSettings,
    prepare_pa6_resp_restart,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "provenance/xa_ma_capped_reference_20261022/pa6_l4_source.mol2"


def _write_lig_pdb(path: Path, omitted_name: str | None = None) -> None:
    substrate = load_pa6(SOURCE)
    conformer = substrate.mol.GetConformer()
    lines = []
    serial = 1
    for atom in substrate.mol.GetAtoms():
        name = atom.GetProp("name")
        if name == omitted_name:
            continue
        point = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"HETATM{serial:5d} {name:>4s} LIG A   1    "
            f"{point.x:8.3f}{point.y:8.3f}{point.z:8.3f}"
            f"  1.00  0.00          {atom.GetSymbol():>2s}"
        )
        serial += 1
    path.write_text("\n".join([*lines, "END", ""]))


def _write_torsions(path: Path) -> None:
    path.write_text(
        json.dumps(
            [
                {
                    "name": "n_terminal_backbone_01",
                    "atoms": ["N56", "C64", "C72", "C76"],
                },
                {
                    "name": "n_terminal_backbone_02",
                    "atoms": ["C64", "C72", "C76", "C68"],
                },
            ]
        )
    )


def test_prepares_canonical_conformer_and_hashed_manifest(tmp_path: Path) -> None:
    pdb = tmp_path / "bound.pdb"
    torsions = tmp_path / "torsions.json"
    output = tmp_path / "prepared"
    _write_lig_pdb(pdb)
    _write_torsions(torsions)

    manifest = prepare_pa6_resp_conformers(SOURCE, {"bound_01": pdb}, torsions, output)

    assert manifest["atom_count"] == 79
    assert manifest["atom_names_in_canonical_order"][:4] == [
        "N56",
        "C64",
        "C72",
        "C76",
    ]
    assert manifest["gaussian"] == {
        "route": manifest["gaussian"]["route"],
        "cores": 8,
        "memory_gb": 16,
        "charge": 0,
        "multiplicity": 1,
        "jobs_launched": False,
    }
    gaussian = (output / "bound_01/opt_smd.gjf").read_text()
    assert "%nprocshared=8" in gaussian
    assert "%mem=16GB" in gaussian
    assert "wB97XD/def2TZVP" in gaussian
    assert "SCRF=(SMD,Solvent=Water)" in gaussian
    assert "MaxCycles=300" in gaussian
    assert "\n0 1\n" in gaussian
    assert "D 1 2 3 4 F" in gaussian
    assert "Freq" not in gaussian
    assert len(manifest["conformers"][0]["gaussian_input_sha256"]) == 64
    assert json.loads((output / "manifest.json").read_text()) == manifest


def test_rejects_incomplete_lig_mapping_before_output(tmp_path: Path) -> None:
    pdb = tmp_path / "incomplete.pdb"
    torsions = tmp_path / "torsions.json"
    _write_lig_pdb(pdb, omitted_name="H46")
    _write_torsions(torsions)

    with pytest.raises(ValueError, match=r"missing=\['H46'\]"):
        prepare_pa6_resp_conformers(
            SOURCE, {"bound_01": pdb}, torsions, tmp_path / "prepared"
        )
    assert not (tmp_path / "prepared").exists()


def test_invalid_second_conformer_leaves_no_partial_output(tmp_path: Path) -> None:
    valid = tmp_path / "valid.pdb"
    invalid = tmp_path / "invalid.pdb"
    torsions = tmp_path / "torsions.json"
    output = tmp_path / "prepared"
    _write_lig_pdb(valid)
    _write_lig_pdb(invalid, omitted_name="H46")
    _write_torsions(torsions)

    with pytest.raises(ValueError, match=r"missing=\['H46'\]"):
        prepare_pa6_resp_conformers(
            SOURCE,
            {"conformer_01": valid, "conformer_02": invalid},
            torsions,
            output,
        )
    assert not output.exists()


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ([], "nonempty list"),
        (
            [{"name": "bad", "atoms": ["N56", "C64", "C72", "H6"]}],
            "only heavy atoms",
        ),
        (
            [{"name": "bad", "atoms": ["N56", "C64", "C72", "C68"]}],
            "not consecutively bonded",
        ),
    ],
)
def test_rejects_missing_or_invalid_torsions(
    tmp_path: Path, payload: list[dict[str, object]], message: str
) -> None:
    pdb = tmp_path / "bound.pdb"
    torsions = tmp_path / "torsions.json"
    _write_lig_pdb(pdb)
    torsions.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match=message):
        prepare_pa6_resp_conformers(
            SOURCE, {"bound_01": pdb}, torsions, tmp_path / "prepared"
        )


def test_rejects_existing_output_and_implicit_conformer_set(tmp_path: Path) -> None:
    output = tmp_path / "prepared"
    output.mkdir()
    with pytest.raises(FileExistsError):
        prepare_pa6_resp_conformers(SOURCE, {}, tmp_path / "torsions.json", output)

    output.rmdir()
    with pytest.raises(ValueError, match="explicit conformer"):
        prepare_pa6_resp_conformers(SOURCE, {}, tmp_path / "torsions.json", output)


def test_restart_inputs_preserve_order_constraints_and_explicit_provenance(
    tmp_path: Path,
) -> None:
    source_manifest = (
        ROOT / "parameters/pa6_l4_zwitterion_resp_gaff2/01_dft_geometry/manifest.json"
    )
    output = tmp_path / "def2svp_restart"
    settings = GaussianSettings("wB97XD", "def2SVP", "SMD:Water", 300, 8, 16, False)

    manifest = prepare_pa6_resp_restart(source_manifest, output, settings)

    assert manifest["atom_count"] == 79
    assert len(manifest["frozen_torsions"]) == 25
    assert manifest["gaussian"]["method"] == "wB97XD"
    assert manifest["gaussian"]["basis"] == "def2SVP"
    assert manifest["gaussian"]["solvent"] == "Water"
    assert manifest["gaussian"]["calc_fc"] is False
    assert manifest["gaussian"]["jobs_launched"] is False
    assert len(manifest["conformers"]) == 4
    gaussian = (output / "md_replica1_frame1161/opt_smd.gjf").read_text()
    assert "%nprocshared=8" in gaussian
    assert "%mem=16GB" in gaussian
    assert "wB97XD/def2SVP" in gaussian
    assert "SCRF=(SMD,Solvent=Water)" in gaussian
    assert "MaxCycles=300" in gaussian
    assert "CalcFC" not in gaussian
    assert "\n0 1\n" in gaussian
    assert gaussian.count(" F\n") == 25
    assert len(manifest["conformers"][0]["source_xyz_sha256"]) == 64

    with pytest.raises(FileExistsError):
        prepare_pa6_resp_restart(source_manifest, output, settings)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"solvent": "Water"}, "SMD:<Gaussian solvent>"),
        ({"max_cycles": 0}, "must all be positive"),
        ({"method": "wB97XD SCRF"}, "Invalid Gaussian method"),
    ],
)
def test_restart_settings_are_explicit_and_validated(
    kwargs: dict[str, object], message: str
) -> None:
    values: dict[str, object] = {
        "method": "wB97XD",
        "basis": "def2SVP",
        "solvent": "SMD:Water",
        "max_cycles": 300,
        "cores": 8,
        "memory_gb": 16,
        "calc_fc": False,
    }
    values.update(kwargs)
    with pytest.raises(ValueError, match=message):
        GaussianSettings(**values)  # type: ignore[arg-type]
