"""Focused tests for MCCE titration-curve validation and charge grouping."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import plot_mcce_titration_audit as audit  # noqa: E402


def _native(
    tmp_path: Path, *, omit_row: str | None = None, bad_occupancy: str | None = None
) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    header = "iConf CONFORMER FL occ crg Em0 pKa0 ne nH\n"
    conformers = [
        ("NTR01E0001_001", 0.0, [0.8, 0.1]),
        ("NTR+1E0001_002", 1.0, [0.2, 0.9]),
        ("ASP01E0040_001", 0.0, [0.1, 0.1]),
        ("ASP-1E0040_002", -1.0, [0.9, 0.9]),
    ]
    # Keep all four target chains present; duplicate-site rows allow simple fixtures.
    for chain in "FGH":
        conformers.extend(
            [
                (f"NTR01{chain}0001_001", 0.0, [0.8, 0.1]),
                (f"NTR+1{chain}0001_002", 1.0, [0.2, 0.9]),
                (f"ASP01{chain}0040_001", 0.0, [0.1, 0.1]),
                (f"ASP-1{chain}0040_002", -1.0, [0.9, 0.9]),
            ]
        )
    # Curves require both target Asp positions, too.
    for chain in "EFGH":
        conformers.extend(
            [
                (f"ASP01{chain}0042_001", 0.0, [0.1, 0.1]),
                (f"ASP-1{chain}0042_002", -1.0, [0.9, 0.9]),
            ]
        )
    head = [header]
    fort = [" ph 0.0 14.0\n"]
    for index, (name, charge, values) in enumerate(conformers, 1):
        head.append(f"{index:05d} {name} f 0.00 {charge:.3f} 0 0 0 0\n")
        if name != omit_row:
            rendered = [str(v) for v in values]
            if name == bad_occupancy:
                rendered[0] = "oops"
            fort.append(f"{name} {' '.join(rendered)}\n")
    (tmp_path / "head3.lst").write_text("".join(head), encoding="utf-8")
    (tmp_path / "fort.38").write_text("".join(fort), encoding="utf-8")
    return tmp_path


def test_charge_grouping_distinguishes_base_and_acid(tmp_path: Path) -> None:
    ph, curves = audit.load_curves(_native(tmp_path))
    assert ph == [0.0, 14.0]
    assert curves["NTR1"]["chains"]["E"]["protonated_fraction"] == [0.2, 0.9]
    assert curves["ASP40"]["chains"]["E"]["protonated_fraction"] == [0.1, 0.1]
    assert curves["ASP42"]["chains"]["H"]["protonated_fraction"] == [0.1, 0.1]


def test_missing_conformer_row_fails(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="IDs differ"):
        audit.load_curves(_native(tmp_path, omit_row="NTR+1E0001_002"))


def test_nonnumeric_occupancy_fails(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="nonnumeric occupancy"):
        audit.load_curves(_native(tmp_path, bad_occupancy="NTR+1E0001_002"))


def test_occupancy_sum_outside_tolerance_fails(tmp_path: Path) -> None:
    folder = _native(tmp_path)
    fort = folder / "fort.38"
    text = fort.read_text(encoding="utf-8").replace(
        "NTR+1E0001_002 0.2 0.9", "NTR+1E0001_002 0.3 0.9"
    )
    fort.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="occupancy sum"):
        audit.load_curves(folder)


def test_nonincreasing_ph_grid_fails(tmp_path: Path) -> None:
    folder = _native(tmp_path)
    fort = folder / "fort.38"
    fort.write_text(
        fort.read_text(encoding="utf-8").replace("ph 0.0 14.0", "ph 2.0 2.0"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="strictly increasing"):
        audit.load_curves(folder)


def test_out_of_range_occupancy_fails(tmp_path: Path) -> None:
    folder = _native(tmp_path)
    fort = folder / "fort.38"
    text = fort.read_text(encoding="utf-8").replace(
        "NTR+1E0001_002 0.2 0.9", "NTR+1E0001_002 1.2 0.9"
    )
    fort.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=r"occupancies must be finite and in \[0, 1\]"):
        audit.load_curves(folder)


def test_main_writes_json_without_fitted_pka(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    native = _native(tmp_path / "native")
    output = tmp_path / "output"

    def skip_plot(ph: list[float], curves: dict[str, object], path: Path) -> None:
        del ph, curves, path

    monkeypatch.setattr(audit, "_plot", skip_plot)
    assert audit.main(["--native", str(native), "--output", str(output)]) == 0
    report = json.loads((output / "titration_curves.json").read_text(encoding="utf-8"))
    assert report["curves"]["ASP40"]["chains"]["E"]["protonated_fraction"] == [0.1, 0.1]
    assert "pKa" in report["interpretation"]
    assert "script" in report["source_sha256"]


def test_plot_uses_protonated_fractions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ph, curves = audit.load_curves(_native(tmp_path / "native"))
    matplotlib = ModuleType("matplotlib")
    matplotlib.use = Mock()
    pyplot = ModuleType("matplotlib.pyplot")
    figure = Mock()
    axes = [Mock(), Mock(), Mock()]
    pyplot.subplots = Mock(return_value=(figure, axes))
    pyplot.close = Mock()
    monkeypatch.setitem(sys.modules, "matplotlib", matplotlib)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", pyplot)

    output = tmp_path / "curves.png"
    audit._plot(ph, curves, output)

    assert axes[1].plot.call_args.args[1] == [0.1, 0.1]
    figure.savefig.assert_called_once_with(output, dpi=180)


def test_png_smoke_when_matplotlib_is_available(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    ph, curves = audit.load_curves(_native(tmp_path / "native"))
    output = tmp_path / "curves.png"
    audit._plot(ph, curves, output)
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
