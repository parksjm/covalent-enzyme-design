"""Permutation-invariant chemistry and full-precision periodic observables."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.io import netcdf_file

from nyl12_covalent import neb_umbrella_observables as obs

ROLES = {
    "carbon": 1,
    "og1": 2,
    "leaving_n": 3,
    "hg1": 4,
    "nalpha": 5,
    "nalpha_h1": 6,
    "nalpha_h2": 7,
}
XYZ = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [3, 0, 0],
        [1.987654, 0, 0],
        [0, 3, 0],
        [0, 4, 0],
        [0.8, 3.6, 0],
    ],
    dtype=float,
)
BOX = np.array([20, 20, 20, 109.471219, 109.471219, 109.471219])
OWNERS = {4: 2, 6: 5, 7: 5}


def measure(xyz: np.ndarray) -> dict:
    return obs.frame_observables(xyz, BOX, ROLES, [4, 6, 7], [2, 3, 5], OWNERS)


def test_proton_permutation_invariant_chemical_cv() -> None:
    first = measure(XYZ)
    swapped = XYZ.copy()
    swapped[[3, 5]] = swapped[[5, 3]]
    second = measure(swapped)
    assert first["q_h"] != pytest.approx(second["q_h"])
    for key in ["cn_progress", "cn_nalpha_hpool", "cn_nleave_hpool", "cn_og1_hpool"]:
        assert first[key] == pytest.approx(second[key], abs=1e-14)
    assert second["owner_changed"]


def test_truncated_octahedron_wrap_invariance_and_precision() -> None:
    cell = np.asarray(
        obs.box_lengths_and_angles_to_vectors(*BOX).value_in_unit(obs.unit.angstroms)
    )
    wrapped = XYZ.copy()
    wrapped[3] += cell[1]
    wrapped[2] -= cell[2]
    first, second = measure(XYZ), measure(wrapped)
    for key in ["q_h", "q_acyl", "cn_progress"]:
        assert first[key] == pytest.approx(second[key], abs=1e-13)
    assert abs(first["q_h"] - round(first["q_h"], 3)) > 1e-5


@pytest.mark.parametrize("bad", [np.full((7, 3), np.nan), np.ones((7, 2))])
def test_invalid_coordinates_rejected(bad: np.ndarray) -> None:
    with pytest.raises(ValueError):
        measure(bad)


def test_all_frames_retained_and_source_hashes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    top = tmp_path / "system.prmtop"
    top.write_text("mock topology")
    data = {
        "POINTERS": [7],
        "ATOMIC_NUMBER": [6, 8, 7, 1, 7, 1, 1],
        "BONDS_INC_HYDROGEN": [3, 9, 1, 12, 15, 1, 12, 18, 1],
        "BONDS_WITHOUT_HYDROGEN": [],
    }
    monkeypatch.setattr(obs, "AmberFormat", lambda _: SimpleNamespace(parm_data=data))
    path = tmp_path / "trajectory.nc"
    changed = XYZ.copy()
    changed[[3, 5]] = changed[[5, 3]]
    with netcdf_file(str(path), "w") as nc:
        for key, size in [("frame", 2), ("atom", 7), ("spatial", 3)]:
            nc.createDimension(key, size)
        nc.createVariable("coordinates", "d", ("frame", "atom", "spatial"))[:] = [
            XYZ,
            changed,
        ]
        nc.createVariable("cell_lengths", "d", ("frame", "spatial"))[:] = [BOX[:3]] * 2
        nc.createVariable("cell_angles", "d", ("frame", "spatial"))[:] = [BOX[3:]] * 2
        nc.createVariable("time", "d", ("frame",))[:] = [0.05, 0.1]
    report = obs.extract_observables(
        top, path, ROLES, list(range(1, 8)), tmp_path / "results"
    )
    assert report["frames"] == 2 and report["frames_censored"] == 0
    assert report["owner_changed_frames"] >= 1
    assert report["hydrogen_pool_1based"] == [4, 6, 7]
    assert len(report["csv_sha256"]) == 64
    assert (tmp_path / "results/observables.csv").read_text().count("\n") == 3
