"""Exercise explicit chemical modes and conditional both-sign downhill dispatch."""

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import run_relay_ts_followup as module
from scripts import run_relay_ts_partial as partial_cli


def state_fixture() -> SimpleNamespace:
    """Return a full-coordinate inventory with distinct reactive atom positions."""
    xyz = np.zeros((20000, 3))
    xyz[19915] = [2.0, 0.0, 0.0]
    xyz[13915] = [0.0, 0.0, 0.0]
    xyz[13916] = [0.5, 0.5, 0.0]
    xyz[13904] = [1.0, 1.0, 0.0]
    active = np.ones(20000, dtype=bool)
    return SimpleNamespace(
        xyz=xyz,
        active=active,
        indices=np.arange(20000),
        restart=None,
        mdin_source="fixture",
    )


@pytest.mark.parametrize("kind", ["attack", "proton_transfer"])
def test_chemical_seed_modes_are_unit_local_and_nonzero(kind: str) -> None:
    """Mode construction changes only declared reactive atoms, independent of MM drift."""
    state = state_fixture()
    direction = module.chemical_direction(state, kind).reshape((-1, 3))
    assert np.linalg.norm(direction) == pytest.approx(1.0)
    assert np.count_nonzero(np.linalg.norm(direction, axis=1)) <= 3
    state.active[:] = False
    with pytest.raises(ValueError):
        module.chemical_direction(state, kind)


@pytest.mark.parametrize("explicit_gate", [False, True])
@pytest.mark.parametrize(
    "kind,sign",
    [("attack", None), ("proton_transfer", None), ("downhill", -1), ("downhill", 1)],
)
def test_followup_dispatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    kind: str,
    sign: int | None,
    explicit_gate: bool,
) -> None:
    """Both uphill searches and conditional downhill signs retain original full arrays."""
    state = state_fixture()
    center = state.xyz.copy()
    direction = np.zeros(state.xyz.size)
    direction[0] = 1.0

    def gradient(
        xyz: np.ndarray, evidence: Path, **kwargs: object
    ) -> tuple[float, np.ndarray]:
        displacement = xyz - center
        return -718.0 - float(np.sum(displacement**2)) / 2, -displacement

    state.gradient = gradient
    monkeypatch.setattr(module, "NativeState", lambda *args: state)
    captured = {}

    def guarded(*args: object, **kwargs: object) -> dict[str, object]:
        captured["call"] = (args, kwargs)
        return {}

    monkeypatch.setattr(module, "run_guarded", guarded)

    def restart(path: Path, xyz: np.ndarray, original: object) -> None:
        path.write_text("fixture")
        assert xyz.shape == (20000, 3)

    monkeypatch.setattr(module, "write_restart", restart)
    monkeypatch.setattr(module, "_dimer_input", lambda *args: "maxstep=0.05")
    monkeypatch.setattr(
        module, "minimization_input", lambda *args, **kwargs: "maxstep=0.05"
    )
    screen = tmp_path / "screen"
    screen.mkdir()
    np.save(screen / "lifted_cartesian_direction.npy", direction)
    (screen / "scientific_acceptance.json").write_text(
        json.dumps(
            {
                "relative_half_step_difference": 0.0,
                "partial_negative_direction_verified": True,
            }
        )
    )
    args = [
        "runner",
        "--model",
        str(tmp_path),
        "--restart",
        str(tmp_path / "restart"),
        "--output",
        str(tmp_path / "out"),
        "--kind",
        kind,
    ]
    if sign is not None:
        args += ["--screen", str(screen), "--sign", str(sign)]
    if explicit_gate and kind == "downhill":
        gate = tmp_path / "chemical_gate.json"
        gate.write_text(
            json.dumps(
                {
                    "passed": True,
                    "restart_sha256": "fixture",
                    "mode_sha256": "fixture",
                    "source_precision_receipt_sha256": "fixture",
                }
            )
        )
        monkeypatch.setattr(module, "digest", lambda path: "fixture")
        args += ["--chemical-sign-gate", str(gate)]
    monkeypatch.setattr(sys, "argv", args)
    module.main()
    assert captured["call"]
    actual_args, _ = captured["call"]
    assert ("maxstep=0.02" if kind == "downhill" else "maxstep=0.01") in actual_args[3]


def test_partial_cli_dispatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The partial screen CLI passes absolute inputs and the explicit worker budget."""
    calls = []
    monkeypatch.setattr(partial_cli, "partial_screen", lambda *args: calls.append(args))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "runner",
            "--model",
            str(tmp_path),
            "--restart",
            str(tmp_path / "rst"),
            "--output",
            str(tmp_path / "out"),
            "--workers",
            "4",
        ],
    )
    partial_cli.main()
    assert calls[0][-1] == 4
