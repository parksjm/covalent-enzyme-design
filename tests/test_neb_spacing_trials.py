"""Check that numerical trial rendering cannot silently change the QM/MM setup."""

from pathlib import Path
import json

import pytest

from scripts import launch_neb_spacing_trials as launcher
from scripts.launch_neb_spacing_trials import render_trial_input


def source_input() -> str:
    """Return an input with distinguishable scientific and numerical controls."""
    return (
        "&cntrl ntr=0, nmropt=0, maxcyc=100000 /\n"
        "&qmmm qmcharge=-1, qmmask='@1-102', qm_theory='DFTB3' /\n"
        "&dlfind\n neb='FROZEN', neb_nimage=25, nebk=0.01,\n"
        " refcrd='/old/start', optalg='LBFGS', maxstep=0.05,\n"
        + "".join(f" coords2({i})='/old/{i}',\n" for i in range(1, 13))
        + " active=':1-500', tol=1.d-4, tole=1.d-5, maxcycle=16000,\n/\n"
    )


@pytest.mark.parametrize("images,optimizer", [(24, "LBFGS"), (48, "LBFGS"), (24, "CG")])
def test_trial_retains_hamiltonian_masks_and_endpoints(
    images: int, optimizer: str
) -> None:
    """Band density and optimizer changes preserve the scientific problem."""
    controls = [Path(f"/new/control_{i:02d}.ncrst") for i in range(13)]
    base = source_input()
    result = render_trial_input(base, controls, images, optimizer)
    assert result.partition("&dlfind")[0] == base.partition("&dlfind")[0]
    for unchanged in (
        "active=':1-500'",
        "neb='FROZEN'",
        "nebk=0.01",
        "tol=1.d-4",
        "tole=1.d-5",
        "maxcycle=16000",
    ):
        assert unchanged in result
    assert f"neb_nimage={images + 1}" in result
    assert f"optalg='{optimizer}'" in result
    assert "maxstep=0.02" in result
    assert f"refcrd='{controls[0]}'" in result
    assert f"coords2(12)='{controls[-1]}'" in result
    assert "/old/" not in result


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "restraint"])
def test_invalid_source_fails_before_launch(mutation: str) -> None:
    """Ambiguous controls or changed scientific restraints must be rejected."""
    base = source_input()
    if mutation == "missing":
        base = base.replace("coords2(12)", "absent(12)")
    elif mutation == "duplicate":
        base = base.replace("neb_nimage=25,", "neb_nimage=25, neb_nimage=49,")
    else:
        base = base.replace("ntr=0", "ntr=1")
    with pytest.raises(ValueError):
        render_trial_input(base, [Path(str(i)) for i in range(13)], 24, "LBFGS")


@pytest.mark.parametrize(
    "count,images,optimizer", [(12, 24, "LBFGS"), (13, 40, "LBFGS"), (13, 24, "FIRE")]
)
def test_unsupported_trial_is_explicit(count: int, images: int, optimizer: str) -> None:
    """Reject unsupported adapter inputs rather than substituting settings."""
    with pytest.raises(ValueError):
        render_trial_input(
            source_input(), [Path(str(i)) for i in range(count)], images, optimizer
        )


def trial_config(tmp_path: Path) -> Path:
    """Create a minimal prepared trial for launch-gate tests."""
    calculation = tmp_path / "calculation"
    calculation.mkdir()
    controls = [tmp_path / f"control_{i}.ncrst" for i in range(13)]
    (calculation / "minimize.mdin").write_text(
        render_trial_input(source_input(), controls, 24, "LBFGS")
    )
    (tmp_path / "preflight.json").write_text(
        json.dumps({"scientific_preparation_passed": True})
    )
    (tmp_path / "source_snapshot").mkdir()
    path = tmp_path / "config.json"
    path.write_text(
        json.dumps(
            {
                "calculation": str(calculation),
                "qm_charge": -1,
                "active_atoms": 4662,
                "regular_images": 24,
                "reactant": str(controls[0]),
                "images": [str(p) for p in controls],
                "sander": "/amber/bin/sander",
                "topology": "/source/topology.prmtop",
            }
        )
    )
    return path


def test_launch_keeps_explicit_single_thread_environment(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A valid trial uses the recorded runner and unbuffered single-thread jobs."""
    path = trial_config(tmp_path)
    captured = {}

    def fake_launch(*args: object, **kwargs: object) -> Path:
        captured.update({"args": args, **kwargs})
        return tmp_path / "job/run.sh"

    monkeypatch.setattr(launcher, "launch_job", fake_launch)
    assert launcher.launch_trial(path) == tmp_path / "job/run.sh"
    assert captured["env_overrides"]["OMP_NUM_THREADS"] == "1"
    assert captured["env_overrides"]["GFORTRAN_UNBUFFERED_ALL"] == "1"
    assert path in captured["input_paths"]
    assert tmp_path / "preflight.json" in captured["input_paths"]


@pytest.mark.parametrize(
    "field,value",
    [
        ("qm_charge", 0),
        ("active_atoms", 4000),
        ("regular_images", 40),
        ("regular_images", 48),
        ("reactant", "/other/start.ncrst"),
    ],
)
def test_bad_trial_configuration_cannot_launch(
    tmp_path: Path,
    field: str,
    value: object,
) -> None:
    """Wrong charge, movable region, seed or band counts fail before execution."""
    path = trial_config(tmp_path)
    data = json.loads(path.read_text())
    data[field] = value
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        launcher.launch_trial(path)


def test_failed_preparation_cannot_launch(tmp_path: Path) -> None:
    """A produced input is insufficient without explicit scientific acceptance."""
    path = trial_config(tmp_path)
    (tmp_path / "preflight.json").write_text("{}")
    with pytest.raises(ValueError, match="scientific gates"):
        launcher.launch_trial(path)
