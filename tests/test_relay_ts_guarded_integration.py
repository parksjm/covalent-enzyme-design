"""Exercise real guard orchestration with a deterministic fake native process."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import relay_ts_guarded as module


@pytest.fixture
def native_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, str]:
    """Provide immutable input paths and a minimal two-atom native state."""
    model = tmp_path / "model"
    (model / "library_build").mkdir(parents=True)
    for p in (
        model / "topology",
        model / "sander",
        model / "start",
        model / "library_build/libdlfind.so",
    ):
        p.write_text("fixture")
    mdin = "&cntrl ntr=0,nmropt=0,cut=10.0 /\n&qmmm qmmask='@1-2',qm_theory='DFTB3',dftb_telec=0.0,qmcharge=-2,scfconv=1.d-12 /\n&dlfind active='@1-2',dimer=.true. /"
    state = SimpleNamespace(
        xyz=np.zeros((2, 3)),
        active=np.array([True, False]),
        mdin_source=mdin,
        topology=model / "topology",
        sander=model / "sander",
    )
    state.gradient = lambda xyz, evidence: (-718.0, np.zeros((2, 3)))
    monkeypatch.setattr(module, "NativeState", lambda *a: state)
    monkeypatch.setattr(module, "geometry", lambda xyz: {"fixture": True})
    monkeypatch.setattr(module, "hydrogen_audit", lambda *a: {"passed": True})
    monkeypatch.setattr(
        module.Rst7, "open", lambda *a: SimpleNamespace(coordinates=np.zeros((2, 3)))
    )
    return model, model / "start", mdin


def fake_native(
    monkeypatch: pytest.MonkeyPatch, text: str, *, active: bool = False, code: int = 0
) -> dict[str, bool]:
    """Install a fake exact-child process and expose termination evidence."""
    flags = {"terminated": False}

    class Child:
        returncode = code
        pid = 123

        def poll(self) -> int | None:
            return None if active else self.returncode

        def terminate(self) -> None:
            flags["terminated"] = True
            self.returncode = -15

        def wait(self, timeout: int = 0) -> int:
            return self.returncode

    def launch(command: list[str], **kwargs: object) -> Child:
        Path(kwargs["cwd"], "native.mdout").write_text(text)
        return Child()

    monkeypatch.setattr(module.subprocess, "Popen", launch)
    return flags


def test_native_success(
    native_fixture: tuple[Path, Path, str],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Independent force, electronic and exterior gates permit a force candidate."""
    model, start, mdin = native_fixture
    fake_native(monkeypatch, "Optimisation converged\nEtot = -450560.0")
    result = module.run_guarded(model, start, tmp_path / "run", mdin)
    assert result["accepted_force_candidate"]
    assert result["process_completed"]
    assert not result["full_active_index_verified"]


def test_active_invalid_trial_terminates_exact_child(
    native_fixture: tuple[Path, Path, str],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A running SCC-invalid child is terminated and its record preserved."""
    model, start, mdin = native_fixture
    flags = fake_native(monkeypatch, "QM calculation failed", active=True)
    with pytest.raises(RuntimeError, match="Guarded native stage failed"):
        module.run_guarded(model, start, tmp_path / "run", mdin)
    assert flags["terminated"]
    assert (tmp_path / "run/process.json").exists()


@pytest.mark.parametrize(
    "original,replacement",
    [("active='@1-2'", "active='@1'"), ("cut=10.0", "cut=8.0"), ("ntr=0", "ntr=1")],
)
def test_changed_scientific_controls_fail(
    native_fixture: tuple[Path, Path, str],
    tmp_path: Path,
    original: str,
    replacement: str,
) -> None:
    """Changed masks, Hamiltonian settings and positional restraints fail closed."""
    model, start, mdin = native_fixture
    with pytest.raises(ValueError):
        module.run_guarded(
            model, start, tmp_path / "run", mdin.replace(original, replacement)
        )


def test_invalid_budgets_fail(
    native_fixture: tuple[Path, Path, str], tmp_path: Path
) -> None:
    """Nonpositive divergence budgets cannot silently disable guards."""
    model, start, mdin = native_fixture
    with pytest.raises(ValueError):
        module.run_guarded(model, start, tmp_path / "run", mdin, hours=0)
