"""Test matched enzyme-path method benchmark preparation and parsing."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.enzyme_path_method_analysis import (
    HARTREE_TO_KCAL_MOL,
    analyze,
    correct_periodic_profile,
)
from nyl12_covalent.enzyme_path_method_benchmark import (
    ANGSTROM_TO_BOHR,
    WINDOWS,
    parse_gaussian_embedding,
    prepare_benchmark,
    write_xtb_embedding,
)


def _gaussian_input(path: Path, atom_count: int = 106) -> None:
    atoms = [f"H {index}.0 0.0 0.0" for index in range(atom_count)]
    path.write_text(
        "%chk=test.chk\n# force charge prop=(field,read)\n\nTitle\n\n-2 1\n"
        + "\n".join(atoms)
        + "\n\n1.0 2.0 3.0 -0.5\n-1.0 -2.0 -3.0 0.5\n\n"
        + "1.0 2.0 3.0\n\n"
    )


def test_gaussian_embedding_conversion_preserves_order_and_units(
    tmp_path: Path,
) -> None:
    source = tmp_path / "input.gjf"
    _gaussian_input(source)
    embedding = parse_gaussian_embedding(source)
    assert embedding.charge == -2
    assert embedding.multiplicity == 1
    assert len(embedding.atoms) == 106
    assert len(embedding.point_charges) == 2
    metadata = write_xtb_embedding(source, tmp_path)
    lines = (tmp_path / "point_charges.pc").read_text().splitlines()
    assert lines[0] == "2"
    fields = [float(value) for value in lines[1].split()]
    assert fields == pytest.approx(
        [
            -0.5,
            1.0 * ANGSTROM_TO_BOHR,
            2.0 * ANGSTROM_TO_BOHR,
            3.0 * ANGSTROM_TO_BOHR,
            999.0,
        ]
    )
    assert metadata["qm_link_atom_count"] == 106
    assert metadata["point_charge_count"] == 2
    assert metadata["point_charge_hardness"] == 999.0
    assert json.loads((tmp_path / "xtb_input_manifest.json").read_text()) == metadata


def test_gaussian_embedding_rejects_wrong_qm_atom_count(tmp_path: Path) -> None:
    source = tmp_path / "input.gjf"
    _gaussian_input(source, atom_count=105)
    with pytest.raises(ValueError, match="106"):
        parse_gaussian_embedding(source)


def test_gaussian_embedding_rejects_missing_charges(tmp_path: Path) -> None:
    source = tmp_path / "input.gjf"
    atoms = "\n".join(f"H {index}.0 0.0 0.0" for index in range(106))
    source.write_text(f"# route\n\nTitle\n\n-2 1\n{atoms}\n\n")
    with pytest.raises(ValueError, match="no external"):
        parse_gaussian_embedding(source)


def test_prepare_benchmark_copies_exact_five_windows(tmp_path: Path) -> None:
    coordinate_drive = tmp_path / "coordinate_drive"
    for index, window in enumerate(WINDOWS):
        folder = coordinate_drive / window
        folder.mkdir(parents=True)
        (folder / "optimized.ncrst").write_text(f"restart {index}\n")
        (folder / "acceptance.json").write_text(
            json.dumps(
                {
                    "native_converged": True,
                    "distances_A": {
                        "OG1_HG1_A": 1.0 + index,
                        "N_HG1_A": 2.0 - index / 10,
                        "proton_difference_A": -1.0 + index / 2,
                    },
                }
            )
        )
    topology = tmp_path / "system.parm7"
    topology.write_text("topology\n")
    output = tmp_path / "benchmark"
    prepare_benchmark(coordinate_drive, topology, output)
    manifest = json.loads((output / "benchmark_manifest.json").read_text())
    assert manifest["windows"] == list(WINDOWS)
    assert manifest["expected_link_atom_count"] == 8
    assert len(manifest["records"]) == 5
    assert (output / "windows" / "pt_19" / "input.ncrst").read_text() == "restart 4\n"
    dftb_input = (
        output / "windows" / "pt_01" / "dftb3_mm" / "minimize.mdin"
    ).read_text()
    assert "qm_ewald=0" in dftb_input
    assert "qm_pme=0" in dftb_input
    assert "qmcut=20.0" in dftb_input


def _completed_window(root: Path, window: str, index: int) -> None:
    folder = root / "windows" / window
    dftb = folder / "dftb3_mm"
    dft = folder / "wb97xd_mm"
    xtb = folder / "gfn2_xtb"
    for path in (dftb, dft, xtb):
        path.mkdir(parents=True)
    (dftb / "single_point.mdout").write_text(
        f"DFTBESCF= {-100.0 + index * 2}\n"
        f"EPtot = {-1000.0 + index * 3}\n"
        "Energy calculation finished\n"
    )
    (dft / "single_point.mdout").write_text(
        f"EXTERNESCF= {-200.0 + index}\n"
        f"EPtot = {-1100.0 + index * 2}\n"
        "Energy calculation finished\n"
    )
    xtb_energy = -50.0 + index * 1.5 / HARTREE_TO_KCAL_MOL
    (xtb / "gfn2.out").write_text(
        f"| TOTAL ENERGY {xtb_energy:.15f} Eh\nnormal termination of xtb\n"
    )


def test_analyze_writes_relative_profiles_and_metrics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    records = []
    for index, window in enumerate(WINDOWS):
        _completed_window(tmp_path, window, index)
        records.append(
            {
                "window": window,
                "acceptance": {
                    "distances_A": {
                        "OG1_HG1_A": 1.0 + index / 10,
                        "N_HG1_A": 2.0 - index / 10,
                        "proton_difference_A": -1.0 + index / 2,
                    }
                },
            }
        )
    (tmp_path / "benchmark_manifest.json").write_text(json.dumps({"records": records}))
    monkeypatch.setattr(
        "nyl12_covalent.enzyme_path_method_analysis._plot",
        lambda rows, output: output.write_bytes(b"PNG"),
    )
    output = tmp_path / "analysis"
    analyze(tmp_path, output)
    rows = json.loads((output / "matched_energies.json").read_text())
    metrics = json.loads((output / "method_errors.json").read_text())
    assert rows[-1]["relative_dftb3_qm_kcal_mol"] == pytest.approx(8.0)
    assert rows[-1]["relative_wb97xd_qm_kcal_mol"] == pytest.approx(4.0)
    assert rows[-1]["relative_gfn2_xtb_qm_kcal_mol"] == pytest.approx(6.0)
    assert metrics["dftb3"]["qm_endpoint_error_kcal_mol"] == pytest.approx(4.0)
    assert metrics["gfn2_xtb"]["qm_endpoint_error_kcal_mol"] == pytest.approx(2.0)
    assert (output / "embedded_method_comparison.png").stat().st_size > 0


def test_correct_periodic_profile_uses_double_difference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    matched = []
    coordinate_drive = tmp_path / "drive"
    for index, window in enumerate(WINDOWS):
        matched.append(
            {
                "window": window,
                "proton_difference_A": -1.0 + index / 2,
                "relative_dftb3_qm_kcal_mol": float(index * 3),
                "relative_wb97xd_qm_kcal_mol": float(index * 2),
            }
        )
        folder = coordinate_drive / window
        folder.mkdir(parents=True)
        (folder / "minimized.mdout").write_text(
            f"DFTBESCF= {-100.0 + index * 4}\n"
            f"EAMBER (non-restraint) = {-1000.0 + index * 5}\n"
        )
    matched_path = tmp_path / "matched.json"
    matched_path.write_text(json.dumps(matched))
    monkeypatch.setattr(
        "nyl12_covalent.enzyme_path_method_analysis._plot_correction",
        lambda rows, output: output.write_bytes(b"PNG"),
    )
    output = tmp_path / "corrected"
    correct_periodic_profile(matched_path, coordinate_drive, output)
    result = json.loads((output / "dft_corrected_periodic_profile.json").read_text())
    endpoint = result["rows"][-1]
    assert endpoint["finite_embedding_dftb3_minus_dft_kcal_mol"] == 4.0
    assert endpoint["relative_periodic_dftb3_hybrid_kcal_mol"] == 20.0
    assert endpoint["dft_corrected_hybrid_kcal_mol"] == 16.0
    assert result["raw_finite_cutoff_profile_accepted"] is False
    assert result["stationary_point_claim"] is False
