"""Tests for immutable native-Thr calibration and explicit data failures."""

from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.native_thr_calibration import (
    LIMITS,
    frozen_geometry,
    gate_failures,
    measure,
    read_sites,
    summarize_controls,
    write_table,
)

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "qmmm/pa6_flipped_model30/10_wt_protonation_md/analysis_after_md"
GEOMETRY = FROZEN / "source_snapshot/src/nyl12_covalent/catalytic_alignment_geometry.py"
CIF = FROZEN / "source_inputs/n12c.cif"
CONTROLS = (
    ROOT
    / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/analysis/occupied_sites.csv"
)


def test_gate_boundaries_and_each_failure() -> None:
    assert gate_failures(dict(LIMITS)) == []
    for key in LIMITS:
        row = dict(LIMITS)
        row[key] += 1e-8
        assert gate_failures(row) == [key]
    for key in ("psi_delta", "chi1_delta"):
        row = dict(LIMITS)
        row[key] *= -1
        assert gate_failures(row) == []


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1.0])
def test_invalid_displacements_fail(value: float) -> None:
    with pytest.raises(ValueError):
        gate_failures({**LIMITS, "og1_shift": value})


def test_cif_auth_label_mapping_and_rigid_invariance() -> None:
    pytest.importorskip(
        "Bio", reason="CIF calibration requires the analysis environment"
    )
    sites, metadata = read_sites(CIF)
    assert set(sites) == {"A", "F", "H", "V"}
    assert metadata["F"]["label_chain"] == "B"
    assert metadata["A"]["residue_map"]["233"] == "221"
    reference = sites["A"]
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    mobile = {
        key: (element, point @ rotation + 123)
        for key, (element, point) in reference.items()
    }
    result = measure(reference, mobile, frozen_geometry(GEOMETRY))
    assert result["native_pass"]
    for metric in (*LIMITS, "scaffold_rmsd"):
        assert result[metric] == pytest.approx(0, abs=1e-10)
    mobile[(233, "OG1")][1][0] += 2.0
    result = measure(reference, mobile, frozen_geometry(GEOMETRY))
    assert result["og1_shift"] == pytest.approx(2.0)
    assert "og1_shift" in result["failures"]


def test_all_control_flags_reproduced_and_csv_projection(tmp_path: Path) -> None:
    summaries = summarize_controls(CONTROLS)
    assert len(summaries) == 22
    assert all(row["sites"] == 20 for row in summaries)
    assert all(row["models"] == 5 for row in summaries)
    path = tmp_path / "controls.csv"
    write_table(path, summaries)
    assert len(path.read_text().splitlines()) == 23


@pytest.mark.parametrize("duplicate", [True, False])
def test_duplicate_or_disagreeing_control_rejected(
    tmp_path: Path, duplicate: bool
) -> None:
    import csv

    with CONTROLS.open() as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        row = next(reader)
    path = tmp_path / "controls.csv"
    if not duplicate:
        row["native_pass"] = "False" if row["native_pass"] == "True" else "True"
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([row, row] if duplicate else [row])
    with pytest.raises(ValueError):
        summarize_controls(path)


def test_empty_cif_fails(tmp_path: Path) -> None:
    pytest.importorskip(
        "Bio", reason="CIF calibration requires the analysis environment"
    )
    path = tmp_path / "empty.cif"
    path.write_text("data_empty\n")
    with pytest.raises(KeyError):
        read_sites(path)


def test_md_scalar_summary_rejects_empty_or_inconsistent(tmp_path: Path) -> None:
    from nyl12_covalent.native_thr_calibration import summarize_md

    with pytest.raises(ValueError, match="No MD"):
        summarize_md([])
    path = tmp_path / "md.csv"
    path.write_text(
        "case,native_og1_shift_A,native_nalpha_shift_A,psi_delta_deg,chi1_delta_deg,native_pass,scaffold_rmsd_A\na__replica_01,0,0,0,0,False,0.1\n"
    )
    with pytest.raises(ValueError, match="MD gate"):
        summarize_md([path])


@pytest.mark.parametrize(
    "script", ["calibrate_native_thr_gate", "summarize_native_thr_diagnostics"]
)
def test_cli_writes_scientific_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, script: str
) -> None:
    import json
    import runpy
    import sys

    if script == "calibrate_native_thr_gate":
        pytest.importorskip(
            "Bio", reason="CIF calibration requires the analysis environment"
        )
    output = tmp_path / "result"
    monkeypatch.setattr(
        sys, "argv", [script, "--root", str(ROOT), "--output", str(output)]
    )
    runpy.run_path(str(ROOT / "scripts" / f"{script}.py"), run_name="__main__")
    acceptance = json.loads((output / "scientific_acceptance.json").read_text())
    if script == "calibrate_native_thr_gate":
        assert acceptance["scientific_acceptance"]["crystal_native_passes"] == 14
        assert (output / "crystal_calibration.png").stat().st_size > 1000
    else:
        assert acceptance["existing_md_frames"] == 120000
        assert acceptance["trajectory_reprocessing"] is False
        assert (output / "known_active_raw_to_relaxed.png").stat().st_size > 1000


def test_missing_biopython_fails_explicitly(monkeypatch: pytest.MonkeyPatch) -> None:
    import builtins

    original_import = builtins.__import__

    def import_without_bio(name: str, *args: object, **kwargs: object) -> object:
        if name == "Bio.PDB.MMCIF2Dict":
            raise ModuleNotFoundError("No module named Bio")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_bio)
    with pytest.raises(
        RuntimeError, match="requires Biopython in the analysis environment"
    ):
        read_sites(CIF)
