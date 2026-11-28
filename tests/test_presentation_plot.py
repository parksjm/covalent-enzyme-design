"""Check source-bound presentation plotting and its stable entry point."""

from pathlib import Path
from typing import Any

import pytest

from PIL import Image

from scripts import plot_nyl12_presentation as plotting


def test_plot_retains_all_source_points(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "path.csv"
    source.write_text("image,relative_physical_kcal_mol\n0,0\n1,16.31\n2,-4.15\n")
    destination = tmp_path / "path.png"
    original = plotting.plt.subplots
    axes = []

    def capture(*args: Any, **kwargs: Any) -> tuple[Any, Any]:
        figure, axis = original(*args, **kwargs)
        axes.append(axis)
        return figure, axis

    monkeypatch.setattr(plotting.plt, "subplots", capture)
    plotting.plot_relay(source, destination)
    assert list(axes[0].lines[0].get_ydata()) == [0, 16.31, -4.15]
    assert axes[0].get_ylabel() == "ΔE (kcal/mol)"
    assert "<svg" in destination.with_suffix(".svg").read_text()
    assert destination.with_suffix(".pdf").read_bytes().startswith(b"%PDF")
    with Image.open(destination) as result:
        assert result.width >= 2800 and result.height >= 1000


def test_entry_point_uses_canonical_sources(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []
    monkeypatch.setattr(plotting, "plot_relay", lambda *paths: calls.append(paths))
    plotting.main()
    source, destination = calls[0]
    assert source.name == "connected_relay_pes.csv"
    assert destination.name == "relay_pes_large.png"


def test_occupancy_plot_uses_replica_means_and_errors(tmp_path: Path) -> None:
    source = tmp_path / "occupancy.csv"
    source.write_text(
        "state,observable,mean,sample_sd\n"
        "nh2_asp272_asp274,tier_a_or_b,0.277675,0.246419\n"
        "nh2_asp272_ash274,tier_a_or_b,0.0518,0.067988\n"
        "nh3_asp272_asp274,tier_a_or_b,0.2899,0.1299\n"
    )
    destination = tmp_path / "occupancy.png"
    plotting.plot_wt_occupancy(source, destination)
    assert destination.is_file()
    assert destination.with_suffix(".svg").is_file()
    assert destination.with_suffix(".pdf").read_bytes().startswith(b"%PDF")


def test_direct_neb_and_pmf_plots_retain_source_values(tmp_path: Path) -> None:
    neb_source = tmp_path / "neb.csv"
    neb_source.write_text(
        "image,relative_energy_kcal_mol\n1,0\n2,8.5\n3,35.81\n4,-2.1\n"
    )
    pmf_source = tmp_path / "pmf.csv"
    pmf_source.write_text(
        "q_acyl_A,F_relative_ES_bin_kcal_mol,bootstrap95_low,bootstrap95_high\n"
        "-1.0,0,-0.5,0.5\n-0.5,20,17,23\n0.0,47.05,42,51\n"
    )
    neb_destination = tmp_path / "neb.png"
    pmf_destination = tmp_path / "pmf.png"
    plotting.plot_direct_neb(neb_source, neb_destination)
    plotting.plot_direct_pmf(pmf_source, pmf_destination)
    assert neb_destination.is_file() and pmf_destination.is_file()
    assert neb_destination.with_suffix(".svg").is_file()
    assert pmf_destination.with_suffix(".pdf").read_bytes().startswith(b"%PDF")


def test_orientation_and_bonded_panel_crops_are_reproducible(tmp_path: Path) -> None:
    orientation = tmp_path / "orientation.png"
    Image.new("RGB", (400, 100), "white").save(orientation)
    orientation_output = tmp_path / "orientation_crop.png"
    plotting.crop_orientation_panel(orientation, orientation_output)
    with Image.open(orientation_output) as cropped:
        assert cropped.size == (352, 80)

    profiles = tmp_path / "profiles.png"
    Image.new("RGB", (400, 200), "white").save(profiles)
    profiles_output = tmp_path / "profiles_crop.png"
    plotting.crop_bonded_profiles(profiles, profiles_output)
    with Image.open(profiles_output) as cropped:
        assert cropped.width == 208
        assert cropped.height == 91
