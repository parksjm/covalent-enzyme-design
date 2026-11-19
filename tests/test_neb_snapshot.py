"""Tests for native save-loop band association and fail-closed snapshots."""

import json
from pathlib import Path

import pytest
from test_neb_progress import native_report

from nyl12_covalent import neb_snapshot


def setup_band(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Create native ASCII restart files and their closed unbuffered save log."""
    monkeypatch.setattr(neb_snapshot, "_require_safe_writer", lambda pid: "test writer")
    mdout = tmp_path / "neb.mdout"
    rows = []
    for image in range(1, 4):
        path = tmp_path / f"neb.{image:04d}.rst7"
        path.write_text(
            "test\n    1\n" + "".join(f"{float(image):12.7f}" for _ in range(3)) + "\n"
        )
        rows.append(f"Saving img {image:3d} to restart: {path}\n".encode())
    mdout.write_bytes(
        native_report([-10, -9, -11])
        + b"".join(rows)
        + b"Testing convergence in cycle 1\n"
    )
    return mdout, tmp_path / "neb.rst7"


def test_snapshot_copies_matching_closed_band(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mdout, base = setup_band(tmp_path, monkeypatch)
    destination = tmp_path / "snapshot"
    report = neb_snapshot.snapshot_complete_band(
        mdout, base, destination, 3, writer_pid=99
    )
    metadata = json.loads((destination / "band_snapshot.json").read_text())
    assert report.cycle == 1
    assert metadata["natom"] == 1
    assert len(metadata["restarts"]) == 3
    assert (destination / "neb.0002.rst7").read_bytes() == (
        tmp_path / "neb.0002.rst7"
    ).read_bytes()
    with pytest.raises(FileExistsError):
        neb_snapshot.snapshot_complete_band(mdout, base, destination, 3, writer_pid=99)


@pytest.mark.parametrize(
    "mutation", ["no_cycle", "missing_save", "new_report", "wrong_path", "truncated"]
)
def test_snapshot_refuses_unproven_band(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    mdout, base = setup_band(tmp_path, monkeypatch)
    payload = mdout.read_bytes()
    if mutation == "no_cycle":
        payload = payload[: payload.index(b"Testing convergence")]
    elif mutation == "missing_save":
        payload = b"\n".join(
            line for line in payload.split(b"\n") if b"Saving img   2" not in line
        )
    elif mutation == "new_report":
        payload += b"NEB Report\n"
    elif mutation == "wrong_path":
        payload = payload.replace(b"neb.0001.rst7", b"bad.0001.rst7")
    else:
        (tmp_path / "neb.0002.rst7").write_bytes(b"test\n1\n")
    mdout.write_bytes(payload)
    with pytest.raises(ValueError):
        neb_snapshot.snapshot_complete_band(
            mdout, base, tmp_path / "snapshot", 3, writer_pid=99
        )
    assert not (tmp_path / "snapshot").exists()


def test_concurrent_restart_mutation_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mdout, base = setup_band(tmp_path, monkeypatch)
    original = Path.read_bytes

    def mutate(path: Path) -> bytes:
        data = original(path)
        if path.name == "neb.0003.rst7":
            path.write_bytes(data + b"\n")
        return data

    monkeypatch.setattr(Path, "read_bytes", mutate)
    with pytest.raises(RuntimeError, match="changed"):
        neb_snapshot.snapshot_complete_band(
            mdout, base, tmp_path / "snapshot", 3, writer_pid=99
        )


def test_netcdf_restart_validation(tmp_path: Path) -> None:
    import numpy as np
    from parmed.amber import Rst7
    from nyl12_covalent.neb_snapshot import _restart_natom

    restart = Rst7(natom=2)
    restart.coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    path = tmp_path / "test.ncrst"
    restart.write(str(path), netcdf=True)
    assert _restart_natom(path.read_bytes()) == 2


def test_wrong_writer_rejected() -> None:
    import os
    from nyl12_covalent.neb_snapshot import _require_safe_writer
    with pytest.raises(ValueError, match='positive'):
        _require_safe_writer(0)
    with pytest.raises(ValueError, match='supported serial sander'):
        _require_safe_writer(os.getpid())
