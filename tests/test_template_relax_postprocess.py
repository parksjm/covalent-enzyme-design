"""Focused checks for matched Amber inspection orchestration."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = (
    ROOT
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/wt_protein_template_controls_20261127"
)
SCRIPT = BASE / "amber_relaxation_20261127/postprocess_template_relaxation.py"


@pytest.fixture
def module() -> ModuleType:
    """Load the retained artifact script without running its CLI."""
    spec = importlib.util.spec_from_file_location("template_postprocess", SCRIPT)
    assert spec is not None and spec.loader is not None
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_pipeline_real_control_and_reject_overwrite(
    module: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise real pairing, all48 geometry, diagnostics and portable publication."""
    source = (
        BASE
        / "calibration_20261127/template_canonical_roundtrip/af3_template_seed20261126/canonical_mature_local.pdb"
    )

    def export(topology: Path, restart: Path, conversion: Path, output: Path) -> dict:
        shutil.copy2(source, output)
        return {"test_export_only": True}

    monkeypatch.setattr(module, "export_canonical_heavy", export)
    args = argparse.Namespace(
        case="af3_template_seed20261126",
        topology=source,
        before_restart=source,
        after_restart=source,
        conversion=source,
        output_dir=tmp_path / "inspection",
    )
    module.process(args)
    root = args.output_dir
    rows = json.loads((root / "before_after_all48.json").read_text())
    assert len(rows) == 48
    assert all(row["delta_d1"] == 0 for row in rows)
    report = json.loads(
        (root / (args.case + "_after") / "geometry_any_amide.json").read_text()
    )
    assert len(report["assigned"]) == len(report["nearest_d1_per_site"]) == 4
    assert {row["site"] for row in report["assigned"]} == {"AE", "BF", "CG", "DH"}
    assert (root / "MANIFEST.sha256").is_file()
    assert "canonical_precursor.pdb" in (root / "view_controls.py").read_text()
    with pytest.raises(FileExistsError):
        module.process(args)


def test_reject_incomplete_or_changed_identity(module: ModuleType) -> None:
    """Fixed-identity deltas fail rather than pair rows by their changing order."""
    with pytest.raises(ValueError, match="48"):
        module.compare_rows([], [])


def test_derived_pair_and_biological_numbering(
    module: ModuleType, tmp_path: Path
) -> None:
    """A nonidentity physical chain assignment is explicit and preserves coordinates."""
    from nyl12_covalent.design_pa6_mapping import _pdb_atom_line

    source = tmp_path / "original.pdb"
    source.write_text(
        _pdb_atom_line("ATOM", 1, "CA", "THR", "F", 1, 1.0, 2.0, 3.0, "C") + "\nEND\n"
    )
    module.derived_views(source, tmp_path, {"A": "F", "B": "E", "C": "G", "D": "H"})
    local = (tmp_path / "canonical_mature_local.pdb").read_text().splitlines()[0]
    precursor = (tmp_path / "canonical_precursor.pdb").read_text().splitlines()[0]
    assert local[21] == precursor[21] == "E"
    assert int(local[22:26]) == 1
    assert int(precursor[22:26]) == 233
    assert local[30:54] == precursor[30:54]
