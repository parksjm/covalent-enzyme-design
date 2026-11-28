"""Known-matrix and strict-budget checks for bounded low-mode diagnostics."""

from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.relay_ts_davidson import davidson, orthogonalize


def test_block_solver_matches_known_dense_spectrum() -> None:
    """A dense rotated indefinite matrix yields independent low-mode residuals."""
    rng = np.random.default_rng(12)
    q, _ = np.linalg.qr(rng.normal(size=(20, 20)))
    expected = np.linspace(0.1, 2.0, 20)
    expected[0] = -0.5
    matrix = q @ np.diag(expected) @ q.T
    result = davidson(
        lambda v: matrix @ v,
        lambda v, e: v,
        rng.normal(size=(20, 4)),
        budget=20,
        tolerance=1e-10,
    )
    np.testing.assert_allclose(result["values"], expected[:3], atol=1e-12)
    assert max(result["residuals"]) < 1e-12
    assert result["products"] <= 20


def test_strict_budget_preserves_unconverged_residual() -> None:
    """Reaching an expensive-product budget records an unresolved calculation."""
    rng = np.random.default_rng(3)
    matrix = np.diag(np.linspace(-0.1, 5.0, 100))
    calls = []

    def operator(vector: np.ndarray) -> np.ndarray:
        calls.append(1)
        return matrix @ vector

    result = davidson(operator, lambda v, e: v, rng.normal(size=(100, 4)), budget=7)
    assert len(calls) == 7
    assert not result["residual_converged"]


def test_projection_removes_duplicate_and_existing_directions() -> None:
    """Twice-projected subspace additions preserve orthogonality."""
    basis = np.eye(4)[:, :1]
    result = orthogonalize(
        np.column_stack((basis[:, 0], np.eye(4)[:, 1], np.eye(4)[:, 1])), basis
    )
    assert result.shape == (4, 1)
    np.testing.assert_allclose(basis.T @ result, 0, atol=1e-15)
    assert orthogonalize(np.empty((4, 0)), basis).shape == (4, 0)


def test_invalid_block_fails() -> None:
    """An insufficient initial subspace fails explicitly."""
    with pytest.raises(ValueError):
        davidson(lambda v: v, lambda v, e: v, np.zeros((10, 4)))
    with pytest.raises(ValueError):
        davidson(lambda v: v, lambda v, e: v, np.ones(5))


def test_full_native_operator_orchestration_harmonic(
    tmp_path: "Path", monkeypatch: pytest.MonkeyPatch
) -> None:
    """Full operator products retain independent h/h2 mode probes within budget."""
    from pathlib import Path
    from types import SimpleNamespace
    from nyl12_covalent import relay_ts_davidson as module

    restart = tmp_path / "restart"
    restart.write_text("fixture")
    screen = tmp_path / "screen"
    screen.mkdir()
    np.save(screen / "partial_columns.npy", np.arange(3))
    np.save(screen / "hessian_raw.npy", np.diag([-0.5, 0.4, 0.5]))
    chemical = np.zeros(12)
    chemical[0] = 1.0
    np.save(screen / "lifted_cartesian_direction.npy", chemical)
    xyz = np.zeros((4, 3))
    diagonal = np.linspace(0.2, 1.4, 12)
    diagonal[0] = -0.5
    state = SimpleNamespace(
        xyz=xyz, active=np.ones(4, dtype=bool), inv_sqrt=np.ones(12)
    )

    def gradient(coordinates: np.ndarray, evidence: Path) -> tuple[float, np.ndarray]:
        vector = coordinates.reshape(-1)
        return float(vector @ (diagonal * vector)) / 2, (diagonal * vector).reshape(
            (-1, 3)
        )

    state.gradient = gradient
    monkeypatch.setattr(module, "NativeState", lambda *args: state)
    result = module.run_davidson(tmp_path, restart, screen, tmp_path / "out")
    np.testing.assert_allclose(
        result["lowest_Ritz_values"], sorted(diagonal)[:3], atol=1e-10
    )
    assert result["operator_products"] <= 80
    assert len(result["independent_mode_probes"]) == 3
    assert not result["full_active_index_verified"]
    native = module.NativeHessianProducts(state, tmp_path / "out")
    with pytest.raises(ValueError):
        native.product(np.ones(3))
    with pytest.raises(ValueError):
        native.product(np.zeros(12))
    native.count = 80
    with pytest.raises(ValueError):
        native.product(np.ones(12))
