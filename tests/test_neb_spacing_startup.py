"""First-frame prefix parsing avoids live-history tails and rejects wrong mapping."""

from pathlib import Path

import numpy as np
import pytest

from scripts.verify_neb_spacing_startup import PendingFrameError, read_first_frame


def frame() -> bytes:
    """Return the installed native fixed-width seven-decimal XYZ format."""
    return (
        "2\n\n"
        + f"{'N':2s}{1.0:15.7f}{2.0:15.7f}{3.0:15.7f}\n"
        + f"{'H':2s}{4.0:15.7f}{5.0:15.7f}{6.0:15.7f}\n"
    ).encode()


def test_reads_only_first_frame_before_incomplete_growing_tail(tmp_path: Path) -> None:
    path = tmp_path / "neb_0001.xyz"
    prefix = frame()
    path.write_bytes(prefix + b"2\n\nINCOMPLETE APPENDING TAIL")
    xyz, payload = read_first_frame(path, ["N", "H"])
    assert payload == prefix
    np.testing.assert_array_equal(xyz, [[1, 2, 3], [4, 5, 6]])
    path.write_bytes(prefix + frame())
    _, again = read_first_frame(path, ["N", "H"])
    assert again == prefix


@pytest.mark.parametrize("missing", [True, False])
def test_missing_or_incomplete_prefix_is_explicit_pending(
    tmp_path: Path, missing: bool
) -> None:
    path = tmp_path / "neb_0001.xyz"
    if not missing:
        path.write_bytes(frame()[:-1])
    with pytest.raises(PendingFrameError):
        read_first_frame(path, ["N", "H"])


@pytest.mark.parametrize("defect", ["count", "comment", "element", "columns", "nan"])
def test_bad_native_mapping_or_values_fail_explicitly(
    tmp_path: Path, defect: str
) -> None:
    path = tmp_path / "neb_0001.xyz"
    data = frame()
    if defect == "count":
        data = data.replace(b"2\n", b"3\n", 1)
    elif defect == "comment":
        data = data.replace(b"2\n\n", b"2\nnotblank\n")
    elif defect == "element":
        data = data.replace(b"N ", b"O ", 1)
    elif defect == "columns":
        data = data.replace(b"1.0000000", b"1.0", 1)
    elif defect == "nan":
        data = data.replace(b"      1.0000000", b"           nan", 1)
    path.write_bytes(data)
    with pytest.raises(ValueError):
        read_first_frame(path, ["N", "H"])
