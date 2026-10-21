"""Exercise Gaussian audit provenance and portable report export paths."""

from __future__ import annotations

import hashlib
import json
import shlex
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from rdkit import Chem

from nyl12_covalent import gaussian_conformer_audit as audit
from nyl12_covalent.gaussian_conformer_geometry import sha256
from scripts import audit_nyl12_gaussian_conformers as cli


ATOMS = [
    {"index": 0, "name": "N", "residue": "THR", "element": "N"},
    {"index": 1, "name": "CA", "residue": "THR", "element": "C"},
    {"index": 2, "name": "C", "residue": "THR", "element": "C"},
    {"index": 3, "name": "N", "residue": "NME", "element": "N"},
]
XYZ = np.asarray([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.0, 1.1, 0.0], [2.8, 1.2, 1.0]])
FREQUENCIES = [float(index + 1) for index in range(57)]


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _write_molecule(path: Path) -> None:
    molecule = Chem.RWMol()
    for atomic_number in (7, 6, 6, 7):
        molecule.AddAtom(Chem.Atom(atomic_number))
    for first, second in ((0, 1), (1, 2), (2, 3)):
        molecule.AddBond(first, second, Chem.BondType.SINGLE)
    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    conformer = Chem.Conformer(result.GetNumAtoms())
    for index, point in enumerate(XYZ):
        conformer.SetAtomPosition(index, tuple(float(value) for value in point))
    result.AddConformer(conformer)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(result)


def _acceptance() -> dict[str, Any]:
    return {
        "accepted": True,
        "completed": True,
        "optimization_completed": True,
        "frequency_minimum_passed": True,
        "expected_frequency_count": 57,
        "frequencies_cm_inverse": FREQUENCIES,
        "electronic_state_passed": True,
        "esp_grid_present": True,
        "fitted_forcefield_ready": False,
        "final_coordinates_A": XYZ.tolist(),
        "geometry": {
            "accepted": True,
            "expected_stereocenters": {},
            "observed_stereocenters": {},
        },
    }


def _fixture_tree(root: Path) -> tuple[Path, Path, Path]:
    prepared = root / "prepared" / "conformer_01"
    calculation = root / "calculations" / "conformer_01"
    prepared.mkdir(parents=True)
    calculation.mkdir(parents=True)
    (prepared / "input.gjf").write_text("neutral test input\n", encoding="utf-8")
    _write_molecule(prepared / "molecule.sdf")
    _write_json(prepared / "atom_map.json", {"atoms": ATOMS})
    for name in ("input.gjf", "molecule.sdf", "atom_map.json"):
        (calculation / name).write_bytes((prepared / name).read_bytes())

    engine_root = root / "fake Gaussian"
    executable = engine_root / "g16" / "g16"
    profile = engine_root / "g16" / "bsd" / "g16.profile"
    espgen = root / "fake AmberTools" / "espgen"
    for path, contents in (
        (executable, "#!/bin/sh\nexit 0\n"),
        (profile, "# test profile\n"),
        (espgen, "#!/bin/sh\nexit 0\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents, encoding="utf-8")
    shell = calculation / "run_gaussian.sh"
    shell.write_text(
        "#!/usr/bin/env bash\n"
        f"source {shlex.quote(str(profile))}\n"
        f"{shlex.quote(str(executable))} < input.gjf > gaussian.log 2>&1\n",
        encoding="utf-8",
    )
    esp_command = calculation / "espgen_command.json"
    _write_json(esp_command, {"argv": [str(espgen), "-i", "gaussian.log"]})

    acceptance_source = Path(audit.__file__).resolve().parent / "gaussian_checks.py"
    source_dir = calculation / "source"
    source_dir.mkdir()
    (source_dir / acceptance_source.name).write_bytes(acceptance_source.read_bytes())
    input_hashes = {
        name: sha256(prepared / name)
        for name in ("input.gjf", "molecule.sdf", "atom_map.json")
    }
    identity = {
        "input_sha256": input_hashes,
        "source_sha256": {"gaussian_checks.py": sha256(acceptance_source)},
        "acceptance_source_sha256": sha256(acceptance_source),
        "engine_sha256": {
            "g16": sha256(executable),
            "profile": sha256(profile),
            "espgen": sha256(espgen),
        },
    }
    calculation_record = {
        "identity": identity,
        "source": str(prepared.resolve()),
        "source_sha256": identity["source_sha256"],
        "shell_sha256": sha256(shell),
    }
    _write_json(calculation / "calculation.json", calculation_record)
    _write_json(
        calculation / "process.json",
        {
            "exit_status": 0,
            "started_utc": "2026-10-20T00:00:00Z",
            "ended_utc": "2026-10-20T00:00:01Z",
        },
    )
    _write_json(calculation / "acceptance.json", _acceptance())
    (calculation / "gaussian.log").write_text(
        "SCF Done: E(RwB97XD) = -10.0 A.U.\n"
        "Optimization completed.\n"
        "Normal termination of Gaussian 16\n"
        "Frequencies -- 1.0\n"
        "Normal termination of Gaussian 16\n"
        "SCF Done: E(RHF) = -9.0 A.U.\n"
        "ESP Fit Center\nElectrostatic Properties\n"
        "Normal termination of Gaussian 16\n",
        encoding="utf-8",
    )
    (calculation / "esp.dat").write_text("4 1 0\n", encoding="utf-8")
    _write_json(calculation / "atom_map.json", {"atoms": ATOMS})
    return prepared, calculation, root


def _patch_science(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        audit,
        "assess_gaussian",
        lambda *_: {
            "accepted": True,
            "frequency_minimum_passed": True,
            "frequencies_cm_inverse": FREQUENCIES,
        },
    )
    monkeypatch.setattr(
        audit, "final_orientation", lambda *_: ([7, 6, 6, 7], XYZ.copy())
    )


def test_audit_conformer_checks_provenance_and_returns_energy_and_psi(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, calculation, _ = _fixture_tree(tmp_path)
    _patch_science(monkeypatch)

    result, xyz, elements = audit.audit_conformer("conformer_01", prepared, calculation)

    assert result["provenance"]["process_exit_status"] == 0
    assert result["optimized_wb97xd_energy_Eh"] == -10.0
    assert result["hf_esp_energy_Eh"] == -9.0
    assert np.array_equal(xyz, XYZ)
    assert elements == ["N", "C", "C", "N"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("input", "SHA-256 mismatch"),
        ("source_snapshot", "SHA-256 mismatch"),
        ("source_path", "source path mismatch"),
        ("incomplete_process", "failed or is incomplete"),
    ],
)
def test_provenance_rejects_changed_inputs_and_incomplete_records(
    tmp_path: Path, mutation: str, message: str
) -> None:
    prepared, calculation, _ = _fixture_tree(tmp_path)
    if mutation == "input":
        (prepared / "input.gjf").write_text("changed input\n", encoding="utf-8")
    elif mutation == "source_snapshot":
        (calculation / "source" / "gaussian_checks.py").write_text(
            "changed source\n", encoding="utf-8"
        )
    elif mutation == "source_path":
        record = json.loads((calculation / "calculation.json").read_text())
        record["source"] = str(tmp_path / "wrong-prepared-path")
        _write_json(calculation / "calculation.json", record)
    else:
        _write_json(calculation / "process.json", {"exit_status": 1})

    with pytest.raises(ValueError, match=message):
        audit.verify_provenance("conformer_01", prepared, calculation)


def test_provenance_rejects_missing_process_record(tmp_path: Path) -> None:
    prepared, calculation, _ = _fixture_tree(tmp_path)
    (calculation / "process.json").unlink()

    with pytest.raises(FileNotFoundError):
        audit.verify_provenance("conformer_01", prepared, calculation)


@pytest.mark.parametrize(
    "mutation", ["coordinates", "atom_map_order", "atom_map_element"]
)
def test_audit_conformer_rejects_final_coordinate_and_atom_mapping_mismatches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    prepared, calculation, _ = _fixture_tree(tmp_path)
    _patch_science(monkeypatch)
    if mutation == "coordinates":
        acceptance = json.loads((calculation / "acceptance.json").read_text())
        acceptance["final_coordinates_A"][0][0] += 0.01
        _write_json(calculation / "acceptance.json", acceptance)
        message = "coordinates disagree"
    else:
        atom_map = json.loads((calculation / "atom_map.json").read_text())
        if mutation == "atom_map_order":
            atom_map["atoms"][0]["index"] = 1
        else:
            atom_map["atoms"][0]["element"] = "O"
        for map_path in (prepared / "atom_map.json", calculation / "atom_map.json"):
            _write_json(map_path, atom_map)
        calculation_record = json.loads((calculation / "calculation.json").read_text())
        calculation_record["identity"]["input_sha256"]["atom_map.json"] = sha256(
            calculation / "atom_map.json"
        )
        _write_json(calculation / "calculation.json", calculation_record)
        message = "[Aa]tom map"

    with pytest.raises(ValueError, match=message):
        audit.audit_conformer("conformer_01", prepared, calculation)


def test_cli_writes_report_and_optimized_exports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calculation_root = tmp_path / "calculations"
    prepared_root = tmp_path / "prepared"
    rows = []
    for index, name in enumerate(cli.CONFORMERS):
        folder = prepared_root / name
        folder.mkdir(parents=True)
        _write_molecule(folder / "molecule.sdf")
        rows.append(
            {
                "name": name,
                "optimized_wb97xd_energy_Eh": -10.0 + 0.001 * index,
                "hf_esp_energy_Eh": -9.0,
                "optimized_coordinates_A": XYZ.tolist(),
                "atom_map": ATOMS,
            }
        )
    monkeypatch.setattr(
        cli,
        "audit_conformer",
        lambda name, *_: (
            next(row for row in rows if row["name"] == name),
            XYZ.copy(),
            ["N", "C", "C", "N"],
        ),
    )
    output = tmp_path / "report.json"
    structures = tmp_path / "optimized"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "audit_nyl12_gaussian_conformers.py",
            "--calculation-root",
            str(calculation_root),
            "--prepared-root",
            str(prepared_root),
            "--output",
            str(output),
            "--structures-dir",
            str(structures),
        ],
    )

    cli.main()

    report = json.loads(output.read_text())
    assert report["readiness"]["charge_model_or_fitted_force_field_ready"] is False
    assert report["conformers"][0]["relative_optimized_wb97xd_energy_Eh"] == 0.0
    assert report["conformers"][1][
        "relative_optimized_wb97xd_energy_kcal_mol"
    ] == pytest.approx(0.627509474)
    for name in cli.CONFORMERS:
        row = next(item for item in report["conformers"] if item["name"] == name)
        for kind in ("sdf", "xyz"):
            artifact = structures / row["optimized_structure_files"][kind]["path"]
            assert artifact.is_file()
            assert (
                hashlib.sha256(artifact.read_bytes()).hexdigest()
                == row["optimized_structure_files"][kind]["sha256"]
            )
