"""Test strict parsing and comparison of Gaussian relaxed torsion scans."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.pa6_torsion_scan_analysis import (
    HARTREE_TO_KCAL_MOL,
    analyze_manifest,
)
from nyl12_covalent.pa6_gaussian_scan_parser import dihedral_degrees, parse_scan_log
from nyl12_covalent.pa6_torsion_scan_plot import unwrap_angles_degrees
from nyl12_covalent.pa6_torsion_scan_matching import match_profile_overlap


def _orientation(angle_degrees: float) -> str:
    angle = math.radians(angle_degrees)
    points = [
        (-1.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (math.cos(angle), 1.0, math.sin(angle)),
    ]
    rows = [
        f"{index:6d}{6:11d}{0:12d}{x:16.6f}{y:12.6f}{z:12.6f}"
        for index, (x, y, z) in enumerate(points, start=1)
    ]
    return "\n".join(
        [
            "Standard orientation:",
            "---------------------------------------------------------------------",
            "Center Atomic Atomic Coordinates (Angstroms)",
            "Number Number Type X Y Z",
            "---------------------------------------------------------------------",
            *rows,
            "---------------------------------------------------------------------",
        ]
    )


def _log(angles: list[float], energies: list[float]) -> str:
    blocks = []
    for angle, energy in zip(angles, energies, strict=True):
        blocks.append(
            "\n".join(
                [
                    f" SCF Done:  E(RwB97XD) = {energy:.10f} A.U. after 8 cycles",
                    _orientation(angle),
                    " Optimization completed.",
                ]
            )
        )
    return "\n".join([*blocks, "Normal termination of Gaussian 16", ""])


def test_dihedral_signed_angle_and_invalid_coordinates() -> None:
    coordinates = np.array(
        [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]
    )
    assert abs(abs(dihedral_degrees(coordinates)) - 90.0) < 1e-8
    with pytest.raises(ValueError, match="finite 4 by 3"):
        dihedral_degrees(np.zeros((3, 3)))


def test_unwrap_angles_across_signed_boundary() -> None:
    assert unwrap_angles_degrees([170.0, 179.0, -179.0, -170.0]) == pytest.approx(
        [170.0, 179.0, 181.0, 190.0]
    )
    with pytest.raises(ValueError, match="finite and in"):
        unwrap_angles_degrees([181.0])


def test_match_profile_overlap_handles_wrap_and_validates_arrays() -> None:
    matched = match_profile_overlap([(1, 179.0)], [(1, -179.5)], [1.25], [0.75], 2.0)
    assert matched == [
        {
            "forward_point": 1,
            "reverse_point": 1,
            "angle_degrees": 179.0,
            "hysteresis_kcal_mol": 0.5,
        }
    ]
    with pytest.raises(ValueError, match="matching lengths"):
        match_profile_overlap([(1, 0.0)], [], [], [], 1.0)


def test_parse_scan_log_extracts_each_completed_window(tmp_path: Path) -> None:
    log = tmp_path / "scan.log"
    log.write_text(_log([90.0, 120.0], [-100.0, -99.99]))

    points = parse_scan_log(log, [1, 2, 3, 4], 2)

    assert len(points) == 2
    assert points[0].energy_hartree == pytest.approx(-100.0)
    assert abs(points[0].angle_degrees) == pytest.approx(90.0)
    assert points[1].energy_hartree == pytest.approx(-99.99)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (
            _log([90.0], [-100.0]).replace("Normal termination of Gaussian 16", ""),
            "normal Gaussian termination",
        ),
        (
            _log([90.0], [-100.0]).replace("Optimization completed.", ""),
            "completed scan windows",
        ),
        (
            _log([90.0], [-100.0]).replace("SCF Done:", "SCF energy omitted:"),
            "No SCF energy",
        ),
    ],
)
def test_parse_scan_log_rejects_incomplete_windows(
    tmp_path: Path, text: str, message: str
) -> None:
    log = tmp_path / "incomplete.log"
    log.write_text(text)
    with pytest.raises(ValueError, match=message):
        parse_scan_log(log, [1, 2, 3, 4], 1)


def test_analyze_manifest_writes_profiles_csv_json_and_png(tmp_path: Path) -> None:
    scan_root = tmp_path / "scans"
    scan_root.mkdir()
    directions = {}
    for direction, angles, energies, increment in (
        ("forward", [90.0, 120.0], [-100.0, -99.99], 30.0),
        ("reverse", [90.0, 60.0], [-100.0, -99.98], -30.0),
    ):
        directory = scan_root / "fragment" / "torsion" / direction
        directory.mkdir(parents=True)
        (directory / "scan.log").write_text(_log(angles, energies))
        (directory / "scan.gjf").write_text("input placeholder")
        directions[direction] = {
            "input": str(directory / "scan.gjf"),
            "increment_degrees": increment,
        }
    manifest = {
        "scans": [
            {
                "fragment": "fragment",
                "label": "torsion",
                "coordinate_kind": "proper",
                "steps_per_direction": 1,
                "atom_indices_gaussian_one_based": [1, 2, 3, 4],
                "directions": directions,
            }
        ]
    }
    manifest_path = scan_root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / "analysis"

    report = analyze_manifest(manifest_path, output)

    analysis = report["scans"][0]
    assert len(report["source_manifest_sha256"]) == 64
    assert analysis["coordinate_kind"] == "proper"
    assert len(analysis["forward_log_sha256"]) == 64
    assert len(analysis["reverse_log_sha256"]) == 64
    assert analysis["matched_overlap_count"] == 1
    assert analysis[
        "matched_overlap_mean_abs_energy_difference_kcal_mol"
    ] == pytest.approx(0.0)
    assert analysis["unmatched_forward_count"] == 1
    assert len(analysis["raw_signed_profiles"]["forward"]) == 2
    assert (output / "scan_points.csv").is_file()
    assert (output / "analysis.json").is_file()
    assert (output / "fragment" / "torsion" / "energy_profile.svg").is_file()
    assert (output / "fragment" / "torsion" / "energy_profile.png").is_file()
    assert HARTREE_TO_KCAL_MOL > 627.0


def test_improper_scans_report_only_the_complementary_branch_overlap(
    tmp_path: Path,
) -> None:
    scan_root = tmp_path / "scans"
    directions = {}
    for direction, angles, energies, increment in (
        ("forward", [0.0, 5.0, 10.0], [-50.0, -49.99, -49.98], 5.0),
        ("reverse", [0.0, -5.0, -10.0], [-50.0, -49.99, -49.98], -5.0),
    ):
        directory = scan_root / "frag" / "improper" / direction
        directory.mkdir(parents=True)
        (directory / "scan.log").write_text(_log(angles, energies))
        directions[direction] = {
            "input": str(directory / "scan.gjf"),
            "increment_degrees": increment,
        }
    manifest_path = scan_root / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "scans": [
                    {
                        "fragment": "frag",
                        "label": "improper",
                        "coordinate_kind": "improper",
                        "steps_per_direction": 2,
                        "atom_indices_gaussian_one_based": [1, 2, 3, 4],
                        "directions": directions,
                    }
                ]
            }
        )
    )

    report = analyze_manifest(manifest_path, tmp_path / "analysis")

    result = report["scans"][0]
    assert result["coordinate_kind"] == "improper"
    assert result["matched_overlap_count"] == 1
    assert result["unmatched_forward_count"] == 2
    assert "complementary branches" in result["comparison_interpretation"]
    assert result[
        "matched_overlap_mean_abs_energy_difference_kcal_mol"
    ] == pytest.approx(0.0)
    assert result["raw_signed_profiles"]["reverse"][1][
        "angle_degrees"
    ] == pytest.approx(-175.0, abs=1e-4)


def test_analyze_manifest_requires_normal_gaussian_output(tmp_path: Path) -> None:
    scan_root = tmp_path / "scans"
    directory = scan_root / "fragment" / "torsion" / "forward"
    directory.mkdir(parents=True)
    (directory / "scan.log").write_text("partial output")
    reverse = scan_root / "fragment" / "torsion" / "reverse"
    reverse.mkdir(parents=True)
    (reverse / "scan.log").write_text("partial output")
    manifest = {
        "scans": [
            {
                "fragment": "fragment",
                "label": "torsion",
                "steps_per_direction": 1,
                "atom_indices_gaussian_one_based": [1, 2, 3, 4],
                "directions": {
                    "forward": {
                        "input": str(directory / "scan.gjf"),
                        "increment_degrees": 30,
                    },
                    "reverse": {
                        "input": str(reverse / "scan.gjf"),
                        "increment_degrees": -30,
                    },
                },
            }
        ]
    }
    path = scan_root / "manifest.json"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="normal Gaussian termination"):
        analyze_manifest(path, tmp_path / "analysis")
