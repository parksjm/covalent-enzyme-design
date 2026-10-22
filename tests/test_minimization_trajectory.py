"""Validate actual Amber trajectory frames and compressed inspection export."""

import gzip
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import minimization_trajectory as trajectory
from nyl12_covalent.local_repair_run import write_stage_input


def test_future_minimization_inputs_record_trajectory(tmp_path: Path) -> None:
    mdin = write_stage_input(tmp_path, ":1", [[1, 3], [2, 3]], 100, 0.05, 5.0)
    assert "ntwx=1" in mdin.read_text()
    assert "ioutfm=0" in mdin.read_text()


def _write(path: Path, frames: list[np.ndarray]) -> None:
    with path.open("w") as stream:
        stream.write("actual test frames\n")
        for frame in frames:
            numbers = frame.ravel()
            for start in range(0, len(numbers), 10):
                stream.write(
                    "".join(f"{value:8.3f}" for value in numbers[start : start + 10])
                    + "\n"
                )


def test_fixed_width_negative_fields_and_repeated_frames(tmp_path: Path) -> None:
    path = tmp_path / "test.mdcrd"
    frame = np.array([[-123.456, -123.456, 0], [1, 2, 3], [4, 5, 6], [7, 8, 9]])
    _write(path, [frame, frame, frame + 1])
    frames = list(trajectory.read_frames(path, 4))
    assert len(frames) == 3
    assert frames[0] == pytest.approx(frame)
    assert np.array_equal(frames[0], frames[1])
    assert frames[2] == pytest.approx(frame + 1)


@pytest.mark.parametrize(
    "payload",
    [
        "",
        "title\n",
        "title\n 1.0\n",
        "title\n********\n",
        "title\n     nan\n",
        "title\n   1.000\n",
    ],
)
def test_rejects_bad_or_incomplete_frames(tmp_path: Path, payload: str) -> None:
    path = tmp_path / "bad.mdcrd"
    path.write_text(payload)
    with pytest.raises(ValueError):
        list(trajectory.read_frames(path, 1))


def test_atom_count_must_be_positive(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Positive"):
        list(trajectory.read_frames(tmp_path / "unused", 0))


def test_compressed_full_system_pdb_has_real_frames(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    top, coordinates, output = (
        tmp_path / "model.prmtop",
        tmp_path / "coords.mdcrd",
        tmp_path / "export",
    )
    top.write_text("topology fixture")
    frame = np.arange(12).reshape(4, 3).astype(float)
    _write(coordinates, [frame, frame + 0.1])
    atoms = [
        SimpleNamespace(
            idx=i,
            name=name,
            atomic_number=element,
            residue=SimpleNamespace(idx=residue, name=resname),
        )
        for i, (name, element, residue, resname) in enumerate(
            [
                ("N", 7, 0, "ALA"),
                ("OG1", 8, 232, "THR"),
                ("C54", 6, 328, "LIG"),
                ("O", 8, 329, "WAT"),
            ]
        )
    ]
    residues = [SimpleNamespace(name="ALA") for _ in range(330)]
    for index, name in ((232, "THR"), (328, "LIG"), (329, "WAT")):
        residues[index].name = name
    monkeypatch.setattr(
        trajectory,
        "AmberParm",
        lambda path: SimpleNamespace(atoms=atoms, residues=residues),
    )
    report = trajectory.export_trajectory(top, coordinates, output)
    assert report["frames"] == 2 and report["atoms_per_frame"] == 4
    assert report["full_system"]
    content = gzip.open(output / "optimization_trajectory.pdb.gz", "rt").read()
    assert content.count("MODEL ") == 2
    assert content.count("ENDMDL") == 2
    assert content.count("ATOM  ") == 8
    assert "THR E   1" in content
    with pytest.raises(FileExistsError):
        trajectory.export_trajectory(top, coordinates, output)
