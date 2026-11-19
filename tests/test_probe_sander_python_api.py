"""Focused tests for Amber namelist parsing used by the Python-sander probe."""

from pathlib import Path

import pytest

from scripts.probe_sander_python_api import (
    _active_mask,
    _namelist,
    _validate_coordinates,
    _value,
)


def test_value_parses_fortran_scalars_and_strings() -> None:
    assert _value("1.0d-10,") == pytest.approx(1.0e-10)
    assert _value("-2,") == -2
    assert _value(".true.,") is True
    assert _value("'DFTB3',") == "DFTB3"


def test_namelist_retains_comma_delimited_mask(tmp_path: Path) -> None:
    mdin = tmp_path / "dimer.mdin"
    mdin.write_text("&qmmm\n qmmask='@1-3,5,7-9', qmcharge=-2, dftb_telec=300,\n/\n")

    parsed = _namelist(mdin, "qmmm")

    assert parsed == {"qmmask": "@1-3,5,7-9", "qmcharge": -2, "dftb_telec": 300}


def test_namelist_requires_requested_section(tmp_path: Path) -> None:
    mdin = tmp_path / "dimer.mdin"
    mdin.write_text("&cntrl\n cut=10.0,\n/\n")

    with pytest.raises(ValueError, match="Missing &qmmm"):
        _namelist(mdin, "qmmm")


def test_active_mask_survives_absolute_paths_before_it(tmp_path: Path) -> None:
    mdin = tmp_path / "dimer.mdin"
    mdin.write_text(
        "&dlfind\n refcrd='/root/path/reference.ncrst', active=':1-3,5,7-9',\n/\n"
    )

    assert _active_mask(mdin) == ":1-3,5,7-9"


def test_validate_coordinates_accepts_single_frame_batch() -> None:
    import numpy as np

    coordinates = np.zeros((1, 2, 3))

    assert _validate_coordinates(coordinates, 2).shape == (2, 3)
