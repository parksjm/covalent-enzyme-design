"""Portable tests for the standalone Gaussian log-to-XYZ converter."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.gaussian_log_to_xyz import convert, parse_log, trajectory_text


def _orientation(kind: str, rows: tuple[str, ...], *, unit: str = "Angstroms") -> str:
    """Build one realistic Gaussian orientation table."""
    return "\n".join(
        [
            f" {kind} orientation:",
            " ---------------------------------------------------------------------",
            f" Center     Atomic      Atomic             Coordinates ({unit})",
            " Number     Number       Type             X           Y           Z",
            " ---------------------------------------------------------------------",
            *rows,
            " ---------------------------------------------------------------------",
        ]
    )


ROWS = (
    "      1          8           0        0.000000    0.000000    0.000000",
    "      2          1           0        0.900000    0.000000    0.000000",
)


def test_auto_prefers_standard_and_keeps_all_standard_snapshots() -> None:
    """Auto chooses the Standard stream and retains repeated snapshots."""
    text = (
        "NoSymm\n"
        + _orientation("Input", ROWS)
        + "\n"
        + _orientation("Standard", ROWS)
        + "\n"
        + _orientation("Standard", ROWS)
    )
    parsed = parse_log(text)
    assert len(parsed.frames) == 2
    assert all(frame.orientation == "Standard" for frame in parsed.frames)
    assert parsed.counts == {"Standard": 2, "Input": 1}
    xyz = trajectory_text(parsed, "model.log", "a" * 64)
    assert xyz.count('"orientation":"Standard"') == 2
    assert xyz.splitlines()[2].startswith("O ")
    assert xyz.splitlines()[3].startswith("H ")


def test_input_selection_and_no_symm_log() -> None:
    """Explicit Input selection works when Standard tables are absent."""
    parsed = parse_log("NoSymm\n" + _orientation("Input", ROWS), "input")
    assert len(parsed.frames) == 1
    assert parsed.frames[0].orientation == "Input"


def test_multilink_orientation_tables_are_not_deduplicated() -> None:
    """Identical tables from different Gaussian links remain separate frames."""
    parsed = parse_log(
        _orientation("Standard", ROWS)
        + "\n--Link1--\n"
        + _orientation("Standard", ROWS)
    )
    assert len(parsed.frames) == 2
    assert parsed.frames[0].coordinates == parsed.frames[1].coordinates


def test_truncated_trailing_orientation_requires_opt_in() -> None:
    """Only an EOF-truncated final table can be ignored with the option."""
    complete = _orientation("Standard", ROWS)
    truncated = " Standard orientation:\n -----\n Center Atomic Atomic Coordinates (Angstroms)\n"
    with pytest.raises(ValueError, match="Truncated Standard"):
        parse_log(complete + "\n" + truncated)
    parsed = parse_log(complete + "\n" + truncated, allow_incomplete=True)
    assert len(parsed.frames) == 1
    assert parsed.skipped_trailing == {"Standard": 1, "Input": 0}


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("      1          0           0        0.0 0.0 0.0", "atomic number 0"),
        ("      1        119           0        0.0 0.0 0.0", "atomic number 119"),
        ("      1          8           0        nan 0.0 0.0", "Nonfinite"),
    ],
)
def test_dummy_unsupported_and_nonfinite_atoms_fail(row: str, message: str) -> None:
    """Dummy/unsupported elements and nonfinite coordinates are rejected."""
    with pytest.raises(ValueError, match=message):
        parse_log(_orientation("Standard", (row,)))


def test_atom_center_order_and_element_order_must_match() -> None:
    """A later selected snapshot cannot silently change atom correspondence."""
    reordered = (
        "      2          8           0        0.000000    0.000000    0.000000",
        "      1          1           0        0.900000    0.000000    0.000000",
    )
    with pytest.raises(ValueError, match="Center numbers"):
        parse_log(_orientation("Standard", reordered))
    changed_element = (
        ROWS[0],
        ROWS[1].replace("          1           0", "          6           0"),
    )
    with pytest.raises(ValueError, match="Atom count/order changes"):
        parse_log(
            _orientation("Standard", ROWS)
            + "\n"
            + _orientation("Standard", changed_element)
        )


def test_units_header_must_explicitly_be_angstroms() -> None:
    """Coordinate values in another unit cannot be mislabeled as XYZ Angstroms."""
    with pytest.raises(ValueError, match="not in Angstroms"):
        parse_log(_orientation("Standard", ROWS, unit="Bohr"))


def test_convert_writes_trajectory_and_provenance_without_overwrite(
    tmp_path: Path,
) -> None:
    """The CLI conversion creates the XYZ plus input/source/output hashes."""
    log = tmp_path / "run.log"
    output = tmp_path / "run.xyz"
    log.write_text(
        _orientation("Standard", ROWS) + "\n" + _orientation("Standard", ROWS),
        encoding="utf-8",
    )
    report = convert(log, output, "auto", False)
    sidecar = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert report["selected_frame_count"] == 2
    assert sidecar["source_sha256"] == report["source_sha256"]
    assert sidecar["script_sha256"]
    assert len(output.read_text(encoding="utf-8").splitlines()) == 8
    with pytest.raises(FileExistsError):
        convert(log, output, "auto", False)


def test_sidecar_collision_is_detected_before_trajectory_write(tmp_path: Path) -> None:
    """An existing sidecar prevents either output from being created."""
    log = tmp_path / "run.log"
    output = tmp_path / "run.xyz"
    log.write_text(_orientation("Standard", ROWS), encoding="utf-8")
    output.with_suffix(".json").write_text("existing", encoding="utf-8")
    with pytest.raises(FileExistsError):
        convert(log, output, "auto", False)
    assert not output.exists()
