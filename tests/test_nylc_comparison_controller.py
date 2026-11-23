"""Check dated native comparison gates against malformed scientific artifacts."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[1] / (
    "analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/"
    "repaired_comparison_20261120/source_revision03/run_comparison.py"
)
spec = importlib.util.spec_from_file_location("nylc_dated_comparison", SOURCE)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def occupancy(tmp_path: Path, rows: list[tuple[str, str]]) -> Path:
    """Write a small native-shaped occupancy artifact."""
    path = tmp_path / "fort.38"
    path.write_text(
        "ph "
        + " ".join(map(str, range(15)))
        + "\n"
        + "".join(name + " " + " ".join([value] * 15) + "\n" for name, value in rows)
    )
    return path


def test_native_occupancy_normalization_at_output_precision(tmp_path: Path) -> None:
    """Accept rounding at native three-decimal output precision."""
    ids = [f"ASP0{i}A0040_00{i}" for i in range(3)]
    path = occupancy(tmp_path, [(name, "0.333") for name in ids])
    module.validate_occupancies(path, ids)


@pytest.mark.parametrize("value", ["nan", "inf", "-0.001", "1.001", "0.8"])
def test_reject_invalid_probabilities(tmp_path: Path, value: str) -> None:
    """Reject nonfinite, outside-unit and unnormalized occupancy values."""
    name = "ASP01A0040_001"
    with pytest.raises(ValueError):
        module.validate_occupancies(occupancy(tmp_path, [(name, value)]), [name])


def test_reject_missing_and_duplicated_conformers(tmp_path: Path) -> None:
    """Require exact native conformer inventory including multiplicity."""
    name = "ASP01A0040_001"
    with pytest.raises(ValueError, match="inventory"):
        module.validate_occupancies(
            occupancy(tmp_path, [(name, "1")]), [name, "ASP02A0040_002"]
        )
    with pytest.raises(ValueError, match="inventory"):
        module.validate_occupancies(
            occupancy(tmp_path, [(name, ".5"), (name, ".5")]), [name]
        )


def test_reject_wrong_grid(tmp_path: Path) -> None:
    """Reject a native run that changes the declared pH grid."""
    path = occupancy(tmp_path, [("ASP01A0040_001", "1")])
    path.write_text(path.read_text().replace("ph 0 1 2", "ph -1 1 2"))
    with pytest.raises(ValueError, match="grid"):
        module.validate_occupancies(path, ["ASP01A0040_001"])


def test_upstream_failure_stops_without_receipt(tmp_path: Path) -> None:
    """A nonzero upstream exit must never fall back to older structures."""
    controller = object.__new__(module.Controller)
    controller.base = tmp_path
    controller.profile = {"upstream_exit_key": "exit_status"}
    controller.journal = lambda *args: None
    status = tmp_path / "status.json"
    status.write_text(json.dumps({"exit_status": 1}))
    with pytest.raises(ValueError, match="Upstream failed"):
        controller.wait(status)


def test_reject_undeclared_legacy_adoption(tmp_path: Path) -> None:
    """Missing declared termini require the explicit reviewed legacy contract."""
    with pytest.raises(ValueError, match="Missing explicit"):
        module.validate_adoption(
            {"adoption_contract": "explicit_refinement_termini"},
            {"geometry": {"failures": []}},
            {},
            tmp_path,
        )


def test_reject_geometry_failures(tmp_path: Path) -> None:
    """A rejected refinement cannot proceed to MCCE atom naming."""
    with pytest.raises(ValueError, match="geometry failures"):
        module.validate_adoption(
            {}, {"geometry": {"failures": ["bad_peptide"]}}, {}, tmp_path
        )


def test_wait_missing_status_then_complete(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A not-yet-published verified status is a bounded waiting condition."""
    controller = object.__new__(module.Controller)
    controller.base = tmp_path
    controller.profile = {"upstream_exit_key": "exit_status"}
    controller.journal = lambda *args: None
    controller.status = lambda *args, **kwargs: None
    path = tmp_path / "status.json"
    monkeypatch.setattr(
        module.time,
        "sleep",
        lambda seconds: path.write_text(json.dumps({"exit_status": 0})),
    )
    assert controller.wait(path)["exit_status"] == 0


def test_explicit_custom_exit_contract(tmp_path: Path) -> None:
    """Historical custom exit_code receipts require an explicit profile field."""
    controller = object.__new__(module.Controller)
    controller.base = tmp_path
    controller.profile = {"upstream_exit_key": "exit_code"}
    controller.journal = lambda *args: None
    path = tmp_path / "status.json"
    path.write_text(json.dumps({"exit_code": 0}))
    assert controller.wait(path)["exit_code"] == 0
