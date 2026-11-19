"""Exact native XYZ mapping, anchor guards and coherent reconstruction tests."""

import hashlib
import json
from pathlib import Path

import numpy as np
import parmed as pmd
import pytest
from parmed.amber import Rst7
from test_neb_progress import native_report

from nyl12_covalent import neb_xyz_snapshot as subject


def fixture() -> tuple[pmd.Structure, Rst7, bytes]:
    """Make active indices with an index gap and a native mass-filter exclusion."""
    topology = pmd.Structure()
    for name, number, mass in [
        ("O", 8, 16),
        ("CA", 6, 12),
        ("H", 1, 0.05),
        ("CB", 6, 12),
    ]:
        topology.add_atom(
            pmd.Atom(name=name, atomic_number=number, mass=mass), "ALA", 1
        )
    topology.bonds.append(pmd.Bond(topology.atoms[1], topology.atoms[3]))
    coordinates = np.array(
        [[9.0, 9.0, 9.0], [0.0, 0.0, 0.0], [6.0, 6.0, 6.0], [1.53, 0.0, 0.0]]
    )
    topology.coordinates = coordinates.copy()
    reference = Rst7(natom=4)
    reference.coordinates = coordinates.copy()
    reference.box = [30, 30, 30, 90, 90, 90]
    lines = []
    for shift in [0.0, 0.1]:
        lines += ["           2", " "]
        for x in [shift, 1.53 + shift]:
            lines.append(f"{'C':2}{x:15.7f}{0.0:15.7f}{0.0:15.7f}")
    return topology, reference, ("\n".join(lines) + "\n").encode()


def test_exact_ascending_mapping_with_gap_mass_filter_and_frozen_exterior() -> None:
    topology, reference, payload = fixture()
    frames, active = subject.reconstruct_xyz_band(
        payload, topology, reference, "@2-4", 2, 2
    )
    assert active == [1, 3]
    assert len(frames) == 2
    assert frames[1][1, 0] == pytest.approx(0.1)
    assert frames[1][3, 0] == pytest.approx(1.63)
    np.testing.assert_array_equal(
        frames[1][[0, 2]], reference.coordinates.reshape(-1, 3)[[0, 2]]
    )


@pytest.mark.parametrize(
    "defect",
    [
        "missing_newline",
        "columns",
        "element",
        "count",
        "extra",
        "anchor",
        "nonfinite",
        "short",
    ],
)
def test_bad_native_frames_fail(defect: str) -> None:
    topology, reference, payload = fixture()
    if defect == "missing_newline":
        payload = payload[:-1]
    elif defect == "columns":
        payload = payload.replace(b"C ", b"C", 1)
    elif defect == "element":
        payload = payload.replace(b"C ", b"N ", 1)
    elif defect == "count":
        payload = payload.replace(b"           2", b"           3", 1)
    elif defect == "extra":
        payload += payload
    elif defect == "anchor":
        payload = payload.replace(b"      1.5300000", b"     10.5300000")
    elif defect == "nonfinite":
        payload = payload.replace(b"      1.5300000", b"           nan")
    else:
        payload = payload[:30] + b"\n"
    with pytest.raises(ValueError):
        subject.reconstruct_xyz_band(payload, topology, reference, "@2-4", 2, 2)


def test_mapping_and_anchor_prerequisites_fail() -> None:
    topology, reference, payload = fixture()
    with pytest.raises(ValueError, match="count mismatch"):
        subject.reconstruct_xyz_band(payload, topology, reference, "@2-4", 2, 3)
    topology.bonds.clear()
    with pytest.raises(ValueError, match="anchor"):
        subject.reconstruct_xyz_band(payload, topology, reference, "@2-4", 2, 2)


def setup_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, dict]:
    """Mock only topology loading and writer proof; retain real NetCDF I/O."""
    topology, reference, xyz = fixture()
    ref = tmp_path / "reference.ncrst"
    reference.write(str(ref), netcdf=True)
    top = tmp_path / "topology.parm7"
    top.write_bytes(b"original topology bytes")
    xyz_path = tmp_path / "nebpath.xyz"
    xyz_path.write_bytes(xyz)
    mdout = tmp_path / "minimized.mdout"
    output = native_report([-10.0, -9.0])
    for image in [1, 2]:
        output += (
            f"Saving img {image:3d} to restart: optimized.{image:04d}.ncrst\n".encode()
        )
    output += b"Testing convergence in cycle 10\n"
    mdout.write_bytes(output)
    monkeypatch.setattr(
        subject, "_require_safe_writer", lambda pid: "verified test writer"
    )
    monkeypatch.setattr(subject, "AmberParm", lambda *args, **kwargs: topology)
    return (
        mdout,
        xyz_path,
        {
            "reactant": str(ref),
            "topology": str(top),
            "active_mask": "@2-4",
            "active_atoms": 2,
        },
    )


def test_snapshot_reconstructs_real_netcdf_and_archives_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mdout, xyz, config = setup_snapshot(tmp_path, monkeypatch)
    destination = tmp_path / "snapshot"
    report = subject.snapshot_complete_band_from_xyz(
        mdout, xyz, destination, config, 2, writer_pid=99
    )
    assert report.cycle == 10
    metadata = json.loads((destination / "band_snapshot.json").read_text())
    assert metadata["active_indices_0based"] == [1, 3]
    assert "write_nc_restart" in metadata["native_restart_defect"]
    assert metadata["source_xyz_sha256"] == hashlib.sha256(xyz.read_bytes()).hexdigest()
    result = Rst7.open(str(destination / "optimized.0002.ncrst"))
    assert result.coordinates.reshape(-1, 3)[3, 0] == pytest.approx(1.63)
    np.testing.assert_array_equal(result.box, [30, 30, 30, 90, 90, 90])
    assert (destination / "source_nebpath.xyz").read_bytes() == xyz.read_bytes()


@pytest.mark.parametrize(
    "defect", ["no_cycle", "missing_save", "new_report", "mutation"]
)
def test_snapshot_rejects_unproven_or_changing_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, defect: str
) -> None:
    mdout, xyz, config = setup_snapshot(tmp_path, monkeypatch)
    output = mdout.read_bytes()
    if defect == "no_cycle":
        mdout.write_bytes(output[: output.index(b"Testing convergence")])
    elif defect == "missing_save":
        mdout.write_bytes(
            b"\n".join(
                line for line in output.split(b"\n") if b"Saving img   2" not in line
            )
        )
    elif defect == "new_report":
        mdout.write_bytes(output + b"NEB Report\n")
    else:
        original_read = Path.read_bytes
        mutated = False

        def mutate_during_capture(path: Path) -> bytes:
            nonlocal mutated
            result = original_read(path)
            if path == xyz and not mutated:
                mutated = True
                path.write_bytes(result + b"\n")
            return result

        monkeypatch.setattr(Path, "read_bytes", mutate_during_capture)
    with pytest.raises((ValueError, RuntimeError)):
        subject.snapshot_complete_band_from_xyz(
            mdout, xyz, tmp_path / "snapshot", config, 2, writer_pid=99
        )
    assert not (tmp_path / "snapshot").exists()


def test_native_writer_may_advance_during_offline_reconstruction(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mdout, xyz, config = setup_snapshot(tmp_path, monkeypatch)
    original = subject.reconstruct_xyz_band
    captured_output = mdout.read_bytes()
    captured_xyz = xyz.read_bytes()

    def advance_native(*args: object, **kwargs: object) -> object:
        mdout.write_bytes(captured_output + b"NEB Report\n")
        xyz.write_bytes(b"next in-progress band\n")
        return original(*args, **kwargs)

    monkeypatch.setattr(subject, "reconstruct_xyz_band", advance_native)
    destination = tmp_path / "snapshot"
    report = subject.snapshot_complete_band_from_xyz(
        mdout, xyz, destination, config, 2, writer_pid=99
    )
    assert report.cycle == 10
    assert (destination / "source_nebpath.xyz").read_bytes() == captured_xyz
    assert (
        destination / "source_output_snapshot.mdout"
    ).read_bytes() == captured_output


def test_fixed_inputs_must_remain_stable_during_offline_reconstruction(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mdout, xyz, config = setup_snapshot(tmp_path, monkeypatch)
    original = subject.reconstruct_xyz_band

    def mutate_fixed(*args: object, **kwargs: object) -> object:
        result = original(*args, **kwargs)
        Path(config["topology"]).write_bytes(b"changed topology")
        return result

    monkeypatch.setattr(subject, "reconstruct_xyz_band", mutate_fixed)
    with pytest.raises(RuntimeError, match="Fixed reference or topology"):
        subject.snapshot_complete_band_from_xyz(
            mdout, xyz, tmp_path / "snapshot", config, 2, writer_pid=99
        )
    assert not (tmp_path / "snapshot").exists()
