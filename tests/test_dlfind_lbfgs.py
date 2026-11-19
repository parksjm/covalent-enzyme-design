"""Exercise the installed native optimizer and Python safety boundaries."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.dlfind_lbfgs import DLFindLBFGS

LIBRARY = Path("/data/jpc/code/ambertools26/venv/lib/libdlfind.so")
SOURCE = Path(
    "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/libdlfind/src/dlf_lbfgs.F90"
)
pytestmark = pytest.mark.skipif(
    not LIBRARY.is_file() or not SOURCE.is_file(),
    reason="Installed DL-FIND and ABI source required",
)


def optimizer(n: int = 3) -> DLFindLBFGS:
    """Construct the installed native implementation."""
    return DLFindLBFGS(n, LIBRARY, SOURCE)


def test_quadratic_converges_and_provenance() -> None:
    with optimizer() as opt:
        assert len(opt.provenance["source_sha256"]) == 64
        x = np.array([1.0, -2.0, 3.0])
        stiffness = np.array([1.0, 4.0, 10.0])
        for _ in range(100):
            g = stiffness * x
            if np.max(np.abs(g)) < 1e-7:
                break
            step = opt.step(x, g)
            x += step * min(1.0, 0.2 / np.max(np.abs(step)))
        assert np.max(np.abs(stiffness * x)) < 1e-7


@pytest.mark.parametrize(
    "x,g",
    [
        ([1.0, 2.0], [1.0, 2.0, 3.0]),
        ([1.0, 2.0, 3.0], [1.0, np.nan, 1.0]),
        ([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]),
        ([1.0, 2.0, 3.0], [1e308] * 3),
        ([[1.0, 2.0, 3.0]], [1.0, 2.0, 3.0]),
    ],
)
def test_invalid_vectors(x: list, g: list) -> None:
    with optimizer() as opt, pytest.raises(ValueError):
        opt.step(x, g)


def test_history_restart_and_exclusivity() -> None:
    with optimizer() as opt:
        x, g = np.ones(3), np.ones(3)
        first = opt.step(x, g)
        with pytest.raises(ValueError, match="curvature"):
            opt.step(x, g)
        with pytest.raises(ValueError, match="curvature"):
            opt.step(x + 1.0, g - 0.5)
        with pytest.raises(RuntimeError, match="global state"):
            optimizer()
        with ThreadPoolExecutor(1) as pool:
            with pytest.raises(RuntimeError, match="owning thread"):
                pool.submit(opt.step, x, g).result()
        opt.restart()
        np.testing.assert_array_equal(opt.step(x, g), first)
    opt.close()
    with pytest.raises(RuntimeError, match="closed"):
        opt.restart()
    with optimizer():
        pass


@pytest.mark.parametrize("n,memory", [(0, 5), (3, 0), (True, 5), (2**31, 5)])
def test_invalid_dimensions(n: int, memory: int) -> None:
    with pytest.raises(ValueError):
        DLFindLBFGS(n, LIBRARY, SOURCE, memory)


def test_missing_library(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        DLFindLBFGS(3, tmp_path / "absent.so", SOURCE)
    with optimizer():
        pass
