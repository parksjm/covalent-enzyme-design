"""Tests for strict non-overwriting PA6 amide omega recovery preparation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from nyl12_covalent.amide_omega_recovery import prepare_amide_omega_recovery


def _source_tree(tmp_path: Path) -> tuple[Path, Path, Path]:
    scan_dir = tmp_path / "amide_omega"
    (scan_dir / "forward").mkdir(parents=True)
    manifest_path = tmp_path / "manifest.json"
    gjf = scan_dir / "forward" / "scan.gjf"
    gjf.write_text(
        "%nprocshared=4\n%mem=8GB\n%chk=seed.chk\n"
        "#p wB97XD/def2SVP SCRF=(SMD,Solvent=Water) "
        "Opt=(ModRedundant,Tight,MaxCycles=200) SCF=(Tight,XQC) NoSymm\n\n"
        "PA6 benchmark internal_neutral amide_omega forward\n\n"
        "0 1\nC 0.0 0.0 0.0\nC 1.0 0.0 0.0\nO 1.5 1.0 0.0\n"
        "N 2.0 1.0 0.0\nC 2.5 2.0 0.0\n\nD 1 2 4 5 S 24 15.000000\n"
    )
    lines = [
        "Standard orientation:",
        " ---------------------------------------------------------------------",
        " Center Atomic Atomic Coordinates (Angstroms)",
        " Number Number Type X Y Z",
        " ---------------------------------------------------------------------",
    ]
    for index, (atomic_number, point) in enumerate(
        [
            (6, (0.1, 0.0, 0.0)),
            (6, (1.1, 0.0, 0.0)),
            (8, (1.6, 1.0, 0.0)),
            (7, (2.1, 1.0, 0.0)),
            (6, (2.6, 2.0, 0.0)),
        ],
        start=1,
    ):
        lines.append(
            f"{index:6d}{atomic_number:11d}{0:12d}"
            f"{point[0]:16.6f}{point[1]:12.6f}{point[2]:12.6f}"
        )
    lines.extend(
        [
            " ---------------------------------------------------------------------",
            "Optimization completed.",
            "Standard orientation:",
            " ---------------------------------------------------------------------",
            " Center Atomic Atomic Coordinates (Angstroms)",
            " Number Number Type X Y Z",
            " ---------------------------------------------------------------------",
        ]
    )
    for index, (atomic_number, point) in enumerate(
        [
            (6, (9.1, 0.0, 0.0)),
            (6, (10.1, 0.0, 0.0)),
            (8, (10.6, 1.0, 0.0)),
            (7, (11.1, 1.0, 0.0)),
            (6, (11.6, 2.0, 0.0)),
        ],
        start=1,
    ):
        lines.append(
            f"{index:6d}{atomic_number:11d}{0:12d}"
            f"{point[0]:16.6f}{point[1]:12.6f}{point[2]:12.6f}"
        )
    lines.extend(
        [
            " ---------------------------------------------------------------------",
            "Optimization completed.",
            *["Optimization completed." for _ in range(8)],
            " Error termination request processed by link 9999.",
            " Error termination via Lnk1e in l9999.exe",
        ]
    )
    (scan_dir / "forward" / "scan.log").write_text("\n".join(lines) + "\n")
    manifest = {
        "schema_version": 1,
        "method": "wB97XD/def2SVP/SMD(water)",
        "scans": [
            {
                "fragment": "internal_neutral",
                "label": "amide_omega",
                "coordinate_kind": "proper",
                "atom_indices_zero_based": [0, 1, 3, 4],
                "atom_indices_gaussian_one_based": [1, 2, 4, 5],
                "steps_per_direction": 24,
                "directions": {
                    "forward": {
                        "input": str(gjf),
                        "sha256": hashlib.sha256(gjf.read_bytes()).hexdigest(),
                        "increment_degrees": 15.0,
                    },
                    "reverse": {"increment_degrees": -15.0},
                },
            }
        ],
    }
    manifest_path.write_text(json.dumps(manifest))
    return scan_dir, manifest_path, gjf


def test_recovery_schedules_and_uses_expected_geometries(tmp_path: Path) -> None:
    scan_dir, manifest_path, source_gjf = _source_tree(tmp_path)
    output = tmp_path / "recovery_passed"

    manifest = prepare_amide_omega_recovery(scan_dir, manifest_path, output)

    forward = (output / "forward" / "scan.gjf").read_text()
    reverse = (output / "reverse" / "scan.gjf").read_text()
    assert "Opt=(ModRedundant,Tight,CalcFC,MaxCycles=300)" in forward
    assert "SCF=(Tight,XQC) NoSymm" in forward
    assert "D 1 2 4 5 S 45 5.000000" in forward
    assert "D 1 2 4 5 S 72 -5.000000" in reverse
    assert "C   9.1000000000" in forward
    assert "C   0.0000000000" in reverse
    assert manifest["preparation_status"] == "passed"
    assert manifest["parent"]["optimization_completed_count"] == 10
    assert manifest["parent"]["error_termination_event_count"] == 1
    assert manifest["recovery_schedules"]["forward"]["expected_points"] == 46
    assert manifest["recovery_schedules"]["reverse"]["expected_points"] == 73
    assert manifest["recovery_schedules"]["total_expected_points"] == 119
    assert (
        manifest["source_hashes"]["forward_gjf_sha256"]
        == hashlib.sha256(source_gjf.read_bytes()).hexdigest()
    )


def test_existing_output_directory_is_never_overwritten(tmp_path: Path) -> None:
    scan_dir, manifest_path, _ = _source_tree(tmp_path)
    output = tmp_path / "recovery"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("preserve")

    with pytest.raises(FileExistsError):
        prepare_amide_omega_recovery(scan_dir, manifest_path, output)

    assert marker.read_text() == "preserve"


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("optimization_count", "exactly 10 completed optimizations"),
        ("second_error", "exactly one Gaussian error-termination event"),
        ("normal", "no normal terminations"),
    ],
)
def test_forward_log_must_pass_strict_termination_checks(
    tmp_path: Path, mutation: str, message: str
) -> None:
    scan_dir, manifest_path, _ = _source_tree(tmp_path)
    log_path = scan_dir / "forward" / "scan.log"
    log = log_path.read_text()
    if mutation == "optimization_count":
        log = log.replace("Optimization completed.", "", 1)
    elif mutation == "second_error":
        log += "\n Error termination via Lnk1e in second.exe\n"
    else:
        log += "\n Normal termination of Gaussian 16\n"
    log_path.write_text(log)

    with pytest.raises(RuntimeError, match=message):
        prepare_amide_omega_recovery(scan_dir, manifest_path, tmp_path / "recovery")


def test_parent_schedule_mismatch_fails_before_output_creation(tmp_path: Path) -> None:
    scan_dir, manifest_path, _ = _source_tree(tmp_path)
    manifest = json.loads(manifest_path.read_text())
    manifest["scans"][0]["steps_per_direction"] = 23
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / "recovery"

    with pytest.raises(ValueError, match="24 steps"):
        prepare_amide_omega_recovery(scan_dir, manifest_path, output)

    assert not output.exists()
