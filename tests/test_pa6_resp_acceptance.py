"""Test PA6 DFT acceptance and fixed-geometry ESP preparation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.pa6_resp_acceptance import prepare_pa6_hf_esp
from nyl12_covalent.pa6_resp_preparation import prepare_pa6_resp_conformers

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "provenance/xa_ma_capped_reference_20261022/pa6_l4_source.mol2"


def _prepare(tmp_path: Path) -> tuple[Path, np.ndarray, list[int]]:
    substrate = load_pa6(SOURCE)
    conformer = substrate.mol.GetConformer()
    lines = []
    for serial, atom in enumerate(substrate.mol.GetAtoms(), start=1):
        point = conformer.GetAtomPosition(atom.GetIdx())
        name = atom.GetProp("name")
        lines.append(
            f"HETATM{serial:5d} {name:>4s} LIG A   1    "
            f"{point.x:8.3f}{point.y:8.3f}{point.z:8.3f}"
            f"  1.00  0.00          {atom.GetSymbol():>2s}"
        )
    pdb = tmp_path / "source.pdb"
    pdb.write_text("\n".join([*lines, "END", ""]))
    torsions = tmp_path / "torsions.json"
    torsions.write_text(
        json.dumps([{"name": "backbone_01", "atoms": ["N56", "C64", "C72", "C76"]}])
    )
    prepared = tmp_path / "prepared"
    prepare_pa6_resp_conformers(SOURCE, {"conf_01": pdb}, torsions, prepared)
    xyz_lines = (prepared / "conf_01/coordinates.xyz").read_text().splitlines()[2:]
    periodic = Chem.GetPeriodicTable()
    elements = [periodic.GetAtomicNumber(line.split()[0]) for line in xyz_lines]
    xyz = np.asarray(
        [[float(value) for value in line.split()[1:]] for line in xyz_lines]
    )
    return prepared / "manifest.json", xyz, elements


def _gaussian_log(elements: list[int], xyz: np.ndarray) -> str:
    rows = [
        f"{index:6d}{number:11d}{0:12d}{x:16.6f}{y:12.6f}{z:12.6f}"
        for index, (number, (x, y, z)) in enumerate(
            zip(elements, xyz, strict=True), start=1
        )
    ]
    return "\n".join(
        [
            " #p wB97XD/def2TZVP SCRF=(SMD,Solvent=Water)",
            " Charge = 0 Multiplicity = 1",
            " Standard orientation:",
            " ---------------------------------------------------------------------",
            " Center Atomic Atomic Coordinates (Angstroms)",
            " Number Number Type X Y Z",
            " ---------------------------------------------------------------------",
            *rows,
            " ---------------------------------------------------------------------",
            " Optimization completed.",
            " Normal termination of Gaussian 16",
            "",
        ]
    )


def test_accepted_geometry_writes_fixed_hf_mk_input(tmp_path: Path) -> None:
    manifest, xyz, elements = _prepare(tmp_path)
    log = tmp_path / "gaussian.log"
    log.write_text(_gaussian_log(elements, xyz))
    acceptance = tmp_path / "scientific_acceptance.json"
    esp = tmp_path / "hf_mk.gjf"

    report = prepare_pa6_hf_esp(SOURCE, manifest, "conf_01", log, 2.0, acceptance, esp)

    assert report["accepted"] and report["esp_input_written"]
    text = esp.read_text()
    assert "HF/6-31G(d) Pop=MK IOp(6/33=2,6/42=6)" in text
    assert "Opt" not in text and "SCRF" not in text and "SMD" not in text
    assert "\n0 1\n" in text
    assert json.loads(acceptance.read_text()) == report


def test_acceptance_uses_explicit_method_basis_solvent_metadata(tmp_path: Path) -> None:
    manifest, xyz, elements = _prepare(tmp_path)
    payload = json.loads(manifest.read_text())
    payload["gaussian"].update(
        {
            "method": "wB97XD",
            "basis": "def2SVP",
            "solvation_model": "SMD",
            "solvent": "Water",
        }
    )
    manifest.write_text(json.dumps(payload))
    log = tmp_path / "gaussian.log"
    log.write_text(_gaussian_log(elements, xyz).replace("def2TZVP", "def2SVP"))
    report = prepare_pa6_hf_esp(
        SOURCE,
        manifest,
        "conf_01",
        log,
        2.0,
        tmp_path / "acceptance.json",
        tmp_path / "hf_mk.gjf",
    )
    assert report["accepted"]
    assert report["checks"]["method_basis_and_smd_solvent_match_manifest"]


def test_mismatched_and_missing_method_metadata_are_rejected(tmp_path: Path) -> None:
    manifest, xyz, elements = _prepare(tmp_path)
    payload = json.loads(manifest.read_text())
    payload["gaussian"].update(
        {
            "method": "wB97XD",
            "basis": "def2SVP",
            "solvation_model": "SMD",
            "solvent": "Water",
        }
    )
    manifest.write_text(json.dumps(payload))
    log = tmp_path / "gaussian.log"
    log.write_text(_gaussian_log(elements, xyz))
    with pytest.raises(RuntimeError, match="scientific acceptance"):
        prepare_pa6_hf_esp(
            SOURCE,
            manifest,
            "conf_01",
            log,
            2.0,
            tmp_path / "rejected.json",
            tmp_path / "unused.gjf",
        )
    assert not json.loads((tmp_path / "rejected.json").read_text())["accepted"]

    payload["gaussian"] = {}
    manifest.write_text(json.dumps(payload))
    with pytest.raises(RuntimeError, match="assessment failed"):
        prepare_pa6_hf_esp(
            SOURCE,
            manifest,
            "conf_01",
            log,
            2.0,
            tmp_path / "metadata_error.json",
            tmp_path / "unused2.gjf",
        )
    assert (
        "lacks explicit method/basis/solvent metadata"
        in json.loads((tmp_path / "metadata_error.json").read_text())["failure"]
    )


def test_terminal_proton_transfer_rejected_without_esp(tmp_path: Path) -> None:
    manifest, xyz, elements = _prepare(tmp_path)
    names = json.loads(manifest.read_text())["atom_names_in_canonical_order"]
    xyz[names.index("H1")] = xyz[names.index("O48")] + np.array([0.98, 0.0, 0.0])
    log = tmp_path / "gaussian.log"
    log.write_text(_gaussian_log(elements, xyz))
    acceptance = tmp_path / "scientific_acceptance.json"
    esp = tmp_path / "hf_mk.gjf"

    with pytest.raises(RuntimeError, match="scientific acceptance"):
        prepare_pa6_hf_esp(SOURCE, manifest, "conf_01", log, 2.0, acceptance, esp)

    report = json.loads(acceptance.read_text())
    assert not report["accepted"]
    assert not report["checks"]["terminal_proton_inventory_retained"]
    assert not esp.exists()


def test_malformed_log_writes_rejection_and_no_esp(tmp_path: Path) -> None:
    manifest, _, _ = _prepare(tmp_path)
    log = tmp_path / "gaussian.log"
    log.write_text("Error termination\n")
    acceptance = tmp_path / "scientific_acceptance.json"
    esp = tmp_path / "hf_mk.gjf"

    with pytest.raises(RuntimeError, match="assessment failed"):
        prepare_pa6_hf_esp(SOURCE, manifest, "conf_01", log, 2.0, acceptance, esp)
    assert not json.loads(acceptance.read_text())["accepted"]
    assert not esp.exists()


@pytest.mark.parametrize("tolerance", [0.0, 31.0, float("nan")])
def test_invalid_torsion_tolerance_fails_before_outputs(
    tmp_path: Path, tolerance: float
) -> None:
    with pytest.raises(ValueError, match="Torsion tolerance"):
        prepare_pa6_hf_esp(
            SOURCE,
            tmp_path / "manifest.json",
            "conf_01",
            tmp_path / "gaussian.log",
            tolerance,
            tmp_path / "scientific_acceptance.json",
            tmp_path / "hf_mk.gjf",
        )
