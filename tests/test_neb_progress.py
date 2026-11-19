"""Tests using the exact DL-FIND dlf_neb.F90 fixed-width report format."""

import csv
import hashlib
import json
from pathlib import Path

import pytest

from nyl12_covalent.neb_progress import (
    HARTREE_TO_KCAL_MOL,
    parse_last_complete_neb_report,
    write_neb_progress,
)


def native_report(energies: list[float], ci: tuple[int, float] | None = None) -> bytes:
    """Build the literal Fortran report format from dlf_neb.F90."""
    lines = [
        "NEB Report",
        "             Energy       F_tang    F_perp     Dist"
        "     Angle 1-3 Ang 1-2 Sum ",
    ]
    rows = [("Img ", i, energy) for i, energy in enumerate(energies, 1)]
    if ci is not None:
        rows.append(("Cimg", *ci))
    for label, image, energy in rows:
        lines.append(
            f"{label}{image:4d}{energy:15.7f}"
            f"{0.0:10.5f}{0.1:10.5f}{1.0:10.5f}"
            f"{0.0:8.2f}{0.0:8.2f}{0.0:8.2f}  frozen"
        )
    lines.append(f"Total path length {4.0:10.5f}")
    return ("\n".join(lines) + "\n").encode("ascii")


def test_last_complete_and_cycle_association() -> None:
    first = native_report([-10, -9, -9.5, -8, -11], (4, -7.5))
    payload = first + b"Testing convergence in cycle       7\nNEB Report\nImg "
    report = parse_last_complete_neb_report(payload, 5)
    assert report.cycle == 7
    assert report.report_index == 1
    assert report.report_end_byte == len(first)
    assert report.source_sha256 == hashlib.sha256(payload).hexdigest()
    assert report.provisional_minimum_images == (3,)
    assert report.images[1].relative_energy_kcal_mol == pytest.approx(
        HARTREE_TO_KCAL_MOL
    )
    assert report.climbing_image is not None
    assert report.climbing_image.climbing
    assert report.climbing_image.image == 4
    assert len(report.images) == 5


def test_second_report_does_not_reuse_previous_cycle() -> None:
    payload = native_report([-10, -9, -11]) + b"Testing convergence in cycle 1\n"
    payload += native_report([-10, -8, -11])
    report = parse_last_complete_neb_report(payload, 3)
    assert report.cycle is None
    assert report.report_index == 2
    assert report.climbing_image is None
    assert report.provisional_minimum_images == ()


@pytest.mark.parametrize(
    "payload,count",
    [
        (b"", 3),
        (native_report([-10, -9, -11])[:-1], 3),
        (native_report([-10, -9, -11]), 4),
        (native_report([-10, -9, -11]), 1),
        (native_report([-10, float("nan"), -11]), 3),
        (native_report([-10, -9, -11], (8, -8)), 3),
        (native_report([-10, -9, -11]).replace(b"Img    2", b"Img    1"), 3),
    ],
)
def test_invalid_reports_fail(payload: bytes, count: int) -> None:
    with pytest.raises(ValueError):
        parse_last_complete_neb_report(payload, count)


def test_plateaus_and_endpoints_are_not_minima() -> None:
    report = parse_last_complete_neb_report(native_report([-11, -9, -9, -10]), 4)
    assert report.provisional_minimum_images == ()


def test_exports_keep_ci_separate_and_record_source(tmp_path: Path) -> None:
    source = tmp_path / "neb.mdout"
    source.write_bytes(native_report([-10, -9, -11], (2, -8)))
    destination = tmp_path / "progress"
    report = write_neb_progress(source, destination, 3)
    metadata = json.loads((destination / "neb_progress.json").read_text())
    assert metadata["source_mtime_ns"] == source.stat().st_mtime_ns
    assert metadata["restart_files_associated"] is False
    assert metadata["source_sha256"] == report.source_sha256
    with (destination / "neb_progress.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    assert rows[-1]["climbing"] == "True"
    assert rows[-1]["provisional_minimum"] == "False"


def test_source_mutation_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "neb.mdout"
    source.write_bytes(native_report([-10, -9, -11]))
    original_read = Path.read_bytes

    def mutate_read(path: Path) -> bytes:
        payload = original_read(path)
        path.write_bytes(payload + b"new output\n")
        return payload

    monkeypatch.setattr(Path, "read_bytes", mutate_read)
    with pytest.raises(RuntimeError, match="changed"):
        write_neb_progress(source, tmp_path / "progress", 3)
