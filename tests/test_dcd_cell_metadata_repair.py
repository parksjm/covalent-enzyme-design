"""Require a DCD cell-only derivative to preserve every noncell byte."""

import importlib.util
from pathlib import Path
import struct

import pytest

SOURCE = (
    Path(__file__).resolve().parents[1]
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/supplemental_movies_post5ps_prepare_20261127/repair_dcd_cell_metadata.py"
)
spec = importlib.util.spec_from_file_location("cell_only_repair", SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def rec(payload: bytes) -> bytes:
    return struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))


def read(handle):
    length = struct.unpack("<i", handle.read(4))[0]
    payload = handle.read(length)
    assert handle.read(4) == struct.pack("<i", length)
    return payload


def fixture(tmp_path: Path):
    controls = [0] * 20
    controls[0] = 2
    controls[10] = 1
    boxes = [
        [110.0, 111.0, 112.0, 109.0, 108.0, 107.0],
        [113.0, 114.0, 115.0, 106.0, 105.0, 104.0],
    ]
    data = (
        rec(b"CORD" + struct.pack("<20i", *controls))
        + rec(b"title")
        + rec(struct.pack("<i", 2))
    )
    for box in boxes:
        data += rec(struct.pack("<6d", *box)) + b"".join(
            rec(struct.pack("<2f", i, i + 1)) for i in range(3)
        )
    source = tmp_path / "original.dcd"
    source.write_bytes(data)
    return source, tmp_path / "derived.dcd", boxes


def test_cell_only_success(tmp_path: Path) -> None:
    source, dest, boxes = fixture(tmp_path)
    before = source.read_bytes()
    proof = m.repair_cell_metadata(source, dest, boxes, read)
    assert (
        source.read_bytes() == before
        and proof["header_title_count_coordinate_and_record_marker_bytes_unchanged"]
    )
    ranges = [(r["offset"], r["length"]) for r in proof["allowed_changed_byte_ranges"]]
    assert len(ranges) == 2
    for (offset, _), box in zip(ranges, boxes, strict=True):
        a, b, c, alpha, beta, gamma = box
        assert struct.unpack("<6d", dest.read_bytes()[offset : offset + 48]) == (
            a,
            gamma,
            b,
            beta,
            alpha,
            c,
        )
    with pytest.raises(FileExistsError):
        m.repair_cell_metadata(source, dest, boxes, read)


@pytest.mark.parametrize("kind", ["header", "coordinate"])
def test_reject_any_noncell_change(tmp_path: Path, kind: str) -> None:
    source, dest, boxes = fixture(tmp_path)
    proof = m.repair_cell_metadata(source, dest, boxes, read)
    ranges = [(r["offset"], r["length"]) for r in proof["allowed_changed_byte_ranges"]]
    data = bytearray(dest.read_bytes())
    data[8 if kind == "header" else -10] ^= 1
    dest.write_bytes(data)
    with pytest.raises(ValueError, match="Noncell"):
        m.verify_noncell_bytes(source, dest, ranges)


def test_unproven_source_box_rejected(tmp_path: Path) -> None:
    source, dest, boxes = fixture(tmp_path)
    boxes[0][0] += 1
    with pytest.raises(ValueError, match="not proven"):
        m.repair_cell_metadata(source, dest, boxes, read)


def test_no_original_overwrite(tmp_path: Path) -> None:
    source, _, boxes = fixture(tmp_path)
    with pytest.raises(ValueError, match="immutable"):
        m.repair_cell_metadata(source, source, boxes, read)


def test_wrong_frame_count(tmp_path: Path) -> None:
    source, dest, boxes = fixture(tmp_path)
    with pytest.raises(ValueError, match="layout"):
        m.repair_cell_metadata(source, dest, boxes[:1], read)


def test_extra_data_rejected(tmp_path: Path) -> None:
    source, dest, boxes = fixture(tmp_path)
    source.write_bytes(source.read_bytes() + b"extra")
    with pytest.raises(ValueError, match="extra"):
        m.repair_cell_metadata(source, dest, boxes, read)


def test_size_and_range_validation(tmp_path: Path) -> None:
    source, dest, boxes = fixture(tmp_path)
    m.repair_cell_metadata(source, dest, boxes, read)
    with pytest.raises(ValueError, match="range"):
        m.verify_noncell_bytes(source, dest, [(10, 3)])
    dest.write_bytes(dest.read_bytes()[:-1])
    with pytest.raises(ValueError, match="size"):
        m.verify_noncell_bytes(source, dest, [])
