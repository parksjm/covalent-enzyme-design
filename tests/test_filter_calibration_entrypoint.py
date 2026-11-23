"""Checks for reproducible input preparation and complete cohort extraction."""

import csv
import json
from pathlib import Path
import sys

import pytest

from scripts import calibrate_nyl12_filters as cli


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write the explicit fixture schema."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def setup(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Create a complete two-variant, two-model frozen control panel."""
    sources = tmp_path / "sources"
    sources.mkdir()
    source_paths = [
        "scripts/calibrate_nyl12_filters.py",
        "src/nyl12_covalent/filter_calibration.py",
        "src/nyl12_covalent/filter_calibration_plots.py",
        "tests/test_filter_calibration.py",
        "tests/test_filter_calibration_plots.py",
        "tests/test_filter_calibration_entrypoint.py",
        "pyproject.toml",
        "uv.lock",
    ]
    for relative in source_paths:
        path = sources / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# retained fixture source\n")
    monkeypatch.setattr(cli, "ROOT", sources)
    models = []
    for variant in ["WT", "F55N"]:
        for index in range(2):
            models.append(
                {
                    "sequence_name": f"p_{variant.lower()}",
                    "model": f"p_{variant.lower()}_model_{index}",
                    "score_source": "amber_openmm",
                    "geometry_tier": "B",
                    "d1": 3.4,
                    "d2": 3.9,
                    "d3": 3.9,
                    "a1": 87,
                    "herschlag_phi": -115,
                    "interface_delta": -22.0,
                }
            )
    write_csv(tmp_path / "models.csv", models)
    write_csv(
        tmp_path / "activities.csv",
        [
            {
                "system": "PA6",
                "comparison": "L2 signal",
                "variant": "F55N",
                "relative_activity": 1.3,
            }
        ],
    )
    write_csv(
        tmp_path / "summary.csv",
        [
            {
                "substrate_register": "register",
                "variant": variant,
                "n_models": 2,
                "n_tier_ab": 2,
            }
            for variant in ["WT", "F55N"]
        ],
    )
    config = {
        "inputs": {
            "models": str(tmp_path / "models.csv"),
            "activity": str(tmp_path / "activities.csv"),
            "historical_summary": str(tmp_path / "summary.csv"),
        },
        "cohorts": [
            {
                "id": "pa6",
                "system": "PA6",
                "register": "register",
                "experimental_comparison": "L2 signal",
                "variants": ["WT", "F55N"],
                "expected_models": 2,
                "prefix": "p_",
                "score_source": "amber_openmm",
                "model_input": "models",
                "thresholds": {
                    "d1_max": 3.5,
                    "oxyanion_max": 4.0,
                    "a1_min": 85.0,
                    "a1_max": 125.0,
                    "herschlag_phi_min": 60.0,
                    "herschlag_phi_max": 125.0,
                },
            }
        ],
        "scope": "retrospective fixtures",
        "assay_caveats": ["no raw uncertainty"],
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    return config_path, tmp_path / "analysis"


def test_complete_freeze_and_analysis(setup: tuple[Path, Path]) -> None:
    """Run the real analysis and ensure provenance and sign survive export."""
    config, output = setup
    cli.prepare(config, output)
    cli.analyze(output / "prepared.json")
    report = json.loads((output / "report.json").read_text())
    assert (
        report["cohorts"][0]["per_variant"][1]["median_interface_delta_tier_ab"] == -22
    )
    assert (
        report["cohorts"][0]["per_variant"][1]["selected_endpoint_relative_activity"]
        == 1.3
    )
    assert (output / "plots/calibration.png").stat().st_size > 1000
    acceptance = json.loads((output / "scientific_acceptance.json").read_text())
    assert acceptance["calibration_for_available_derived_labels_complete"]
    assert not acceptance["activity_predictor_validated"]
    with pytest.raises(FileExistsError):
        cli.prepare(config, output)


def test_modified_input_stops_analysis(setup: tuple[Path, Path]) -> None:
    """Reject altered scientific inputs before producing a new report."""
    config, output = setup
    cli.prepare(config, output)
    (output / "inputs/models.csv").write_text("changed\n")
    with pytest.raises(ValueError, match="Frozen input/source changed"):
        cli.analyze(output / "prepared.json")
    assert not (output / "report.json").exists()


def test_historical_disagreement_stops(setup: tuple[Path, Path]) -> None:
    """Historical ensemble counts are independently reproduced, not trusted blindly."""
    config, output = setup
    raw_summary = config.parent / "summary.csv"
    raw_summary.write_text(raw_summary.read_text().replace(",2,2", ",2,1"))
    cli.prepare(config, output)
    with pytest.raises(ValueError, match="Historical count disagreement"):
        cli.analyze(output / "prepared.json")


def test_selection_and_invalid_source(setup: tuple[Path, Path]) -> None:
    """The fixed budget excludes later indices and refuses mixed preparation."""
    config, _ = setup
    cohort = json.loads(config.read_text())["cohorts"][0]
    path = config.parent / "models.csv"
    rows = cli.read_csv(path)
    rows.append({**rows[0], "model": "p_wt_model_99"})
    rows.append({**rows[0], "sequence_name": "different_register"})
    write_csv(path, rows)
    assert len(cli.load_models(path, cohort)) == 4
    rows[0]["score_source"] = "unrelaxed"
    write_csv(path, rows)
    with pytest.raises(ValueError, match="Mixed preparation"):
        cli.load_models(path, cohort)
    rows[0]["score_source"] = "amber_openmm"
    rows[0]["model"] = "p_wt_model_-1"
    write_csv(path, rows)
    with pytest.raises(ValueError, match="Negative model index"):
        cli.load_models(path, cohort)


def test_empty_tables_and_cli(
    setup: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject empty input/output tables and exercise the persistent entry point."""
    config, output = setup
    empty = config.parent / "empty.csv"
    empty.write_text("variant\n")
    with pytest.raises(ValueError, match="Empty source"):
        cli.read_csv(empty)
    with pytest.raises(ValueError, match="empty table"):
        cli.write_table(empty, [])
    monkeypatch.setattr(
        sys,
        "argv",
        ["calibrate", "prepare", "--config", str(config), "--output", str(output)],
    )
    cli.main()
    monkeypatch.setattr(
        sys,
        "argv",
        ["calibrate", "analyze", "--prepared", str(output / "prepared.json")],
    )
    cli.main()
