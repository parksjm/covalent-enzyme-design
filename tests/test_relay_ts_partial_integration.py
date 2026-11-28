"""Block-Hessian assembly and lifted direction checks against harmonic models."""

from concurrent.futures import Future
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import relay_ts_partial as module


class InlinePool:
    """Execute column tasks inline while preserving the future interface."""

    def __init__(self, max_workers: int) -> None:
        self.workers = max_workers

    def __enter__(self) -> "InlinePool":
        return self

    def __exit__(self, *args: object) -> None:
        return None

    def submit(self, function: object, argument: object) -> Future:
        future = Future()
        future.set_result(function(argument))
        return future


@pytest.mark.parametrize("negative", [True, False])
def test_complete_partial_hessian_workflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, negative: bool
) -> None:
    """Known harmonic curvature yields the correct partial spectrum and direction."""
    model = tmp_path / "model"
    model.mkdir()
    restart = model / "restart"
    topology = model / "topology"
    restart.write_text("fixture")
    topology.write_text("fixture")
    xyz = np.zeros((99, 3))
    active = np.array([True] * 98 + [False])
    diagonal = np.ones(294)
    diagonal[0] = -1 if negative else 1
    state = SimpleNamespace(
        xyz=xyz,
        active=active,
        indices=np.arange(98),
        inv_sqrt=np.ones(294),
        topology=topology,
        mdin_source="qmmask='@1-98'",
    )

    def gradient(coordinates: np.ndarray, evidence: Path) -> tuple[float, np.ndarray]:
        cartesian = coordinates[active].reshape(-1)
        value = -718 + float(np.dot(cartesian, diagonal * cartesian)) / 2
        return value, (cartesian * diagonal).reshape((-1, 3))

    state.gradient = gradient
    monkeypatch.setattr(module, "NativeState", lambda *args: state)
    monkeypatch.setattr(module, "AmberParm", lambda *args: None)
    monkeypatch.setattr(
        module, "AmberMask", lambda *args: SimpleNamespace(Selection=lambda: active)
    )
    monkeypatch.setattr(module, "ProcessPoolExecutor", InlinePool)

    def column(args: tuple[Path, Path, Path, int, float]) -> str:
        _, output, _, index, _ = args
        folder = output / "columns" / f"column_{index:05d}"
        folder.mkdir(parents=True)
        vector = np.zeros(294)
        vector[index] = diagonal[index]
        np.save(folder / "column.npy", vector)
        return str(folder / "column.npy")

    monkeypatch.setattr(module, "worker_column", column)
    if negative:
        result = module.partial_screen(model, restart, tmp_path / "output")
        assert result["partial_negative_count"] == 1
        assert result["partial_negative_direction_verified"]
        assert not result["full_active_index_verified"]
        np.testing.assert_allclose(
            result["directional_curvatures_Eh_per_A2"], [-1.0, -1.0]
        )
    else:
        with pytest.raises(RuntimeError, match="no negative mode"):
            module.partial_screen(model, restart, tmp_path / "output")


def test_invalid_resources_fail_before_mutation(tmp_path: Path) -> None:
    """Unauthorized worker limits fail without creating an output attempt."""
    with pytest.raises(ValueError):
        module.partial_screen(
            tmp_path, tmp_path / "restart", tmp_path / "output", workers=9
        )
    assert not (tmp_path / "output").exists()
