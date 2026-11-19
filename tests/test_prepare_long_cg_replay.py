"""Budget-only CG replay input and preserved-input provenance contracts."""

import json
from pathlib import Path
import shlex

import pytest

from scripts.prepare_long_cg_replay import digest, prepare, replay_input


def pilot_input(controls: list[Path]) -> str:
    """Return explicit synthetic Amber input without running a chemistry engine."""
    rows = [
        "CG pilot original title",
        "&cntrl",
        " maxcyc=100000, ntr=0, nmropt=0,",
        "/",
        "&qmmm",
        " qmcharge=-1, qmmask='@1-102', dftb_telec=0.0,",
        "/",
        "&dlfind",
        " optalg='CG', maxstep=0.02, neb_nimage=25, neb='FROZEN',",
        " maxcycle=250, tol=1.d-4, tole=1.d-5, active=':1-99',",
        f" refcrd='{controls[0]}',",
    ]
    rows += [f" coords2({i})='{p}'," for i, p in enumerate(controls[1:], 1)]
    return "\n".join([*rows, "/", ""])


def fixture_pilot(root: Path) -> Path:
    """Create immutable-parent fixture inputs and fake binaries that are never run."""
    pilot = root / "pilot"
    (pilot / "calculation").mkdir(parents=True)
    controls = []
    for i in range(13):
        p = root / f"control_{i:02d}.ncrst"
        p.write_text(f"seed {i}")
        controls.append(p)
    case, topology, seed = (
        root / "case.json",
        root / "system.prmtop",
        root / "seed_manifest.json",
    )
    case.write_text("{}")
    topology.write_text("topology fixture")
    seed.write_text("{}")
    sander = root / "amber/bin/sander"
    sander.parent.mkdir(parents=True)
    sander.write_text("not executable fixture")
    library = root / "amber/lib/libdlfind.so"
    library.parent.mkdir()
    library.write_text("library fixture")
    cfg = {
        "regular_images": 24,
        "active_atoms": 4662,
        "qm_charge": -1,
        "optimizer": "CG",
        "qm_atom_ids": list(range(1, 103)),
        "images": list(map(str, controls)),
        "reactant": str(controls[0]),
        "source_seed_manifest": str(seed),
        "trial_source_sha256": {str(seed): digest(seed)},
        "sander": str(sander),
        "case_path": str(case),
        "topology": str(topology),
    }
    (pilot / "config.json").write_text(json.dumps(cfg))
    (pilot / "preflight.json").write_text('{"scientific_preparation_passed":true}')
    for name in ["startup_precision_review.json", "startup_verification.json"]:
        (pilot / name).write_text("{}")
    (pilot / "calculation/minimize.mdin").write_text(pilot_input(controls))
    return pilot


def test_only_two_budget_tokens_change() -> None:
    original = pilot_input([Path(f"/immutable/source/seed_{i:02d}") for i in range(13)])
    result = replay_input(original)
    assert (
        result.replace("maxcycle=32000", "maxcycle=250").replace(
            "maxcyc=300000", "maxcyc=100000"
        )
        == original
    )
    assert "optalg='CG', maxstep=0.02" in result
    assert "qmcharge=-1, qmmask='@1-102', dftb_telec=0.0" in result
    assert len(result.splitlines()) == len(original.splitlines())


@pytest.mark.parametrize(
    "old,new",
    [
        ("maxcycle=250", "maxcycle=251"),
        ("maxcycle=250", "maxcycle=250, maxcycle=250"),
        ("&dlfind", "&wrong"),
        ("coords2(12)", "coords2(13)"),
        ("neb_nimage=25", "neb_nimage=49"),
    ],
)
def test_invalid_budget_or_control_count_fails(old: str, new: str) -> None:
    original = pilot_input([Path(f"/seed_{i}") for i in range(13)])
    with pytest.raises(ValueError):
        replay_input(original.replace(old, new))


def test_prepare_snapshots_and_launcher_preserve_pilot(tmp_path: Path) -> None:
    pilot = fixture_pilot(tmp_path)
    destination = pilot / "long_run"
    original = (pilot / "calculation/minimize.mdin").read_bytes()
    config_path = prepare(pilot, destination)
    assert (pilot / "calculation/minimize.mdin").read_bytes() == original
    cfg = json.loads(config_path.read_text())
    preflight = json.loads((destination / "preflight.json").read_text())
    assert (
        cfg["native_energy_evaluation_budget"] == 224000
        and preflight["budget_multiplier"] == 128
    )
    assert preflight["budget_only_delta_verified"] is True
    assert preflight["force_converged"] is None
    assert cfg["images"] == json.loads((pilot / "config.json").read_text())["images"]
    assert (
        len(
            [
                p
                for p in (destination / "source_snapshot").iterdir()
                if p.name.endswith(".ncrst")
            ]
        )
        == 13
    )
    args = shlex.split((destination / "launch.sh").read_text().split("exec ", 1)[1])
    assert "nyl12_neb_cg_long_20261119" in args
    assert str(destination / "calculation") == args[args.index("--cwd") + 1]
    assert "-c" in args and args[args.index("-c") + 1] == cfg["images"][0]
    assert not (destination / "job").exists()
    with pytest.raises(FileExistsError):
        prepare(pilot, destination)


def test_wrong_config_image_count_fails_before_writes(tmp_path: Path) -> None:
    pilot = fixture_pilot(tmp_path)
    path = pilot / "config.json"
    cfg = json.loads(path.read_text())
    cfg["regular_images"] = 48
    path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="24images"):
        prepare(pilot, pilot / "long_run")
    assert not (pilot / "long_run").exists()
