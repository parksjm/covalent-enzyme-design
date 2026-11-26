"""Validate in-place continuation preparation and preserved historical identities."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.neb_umbrella_campaign import digest
from nyl12_covalent.relay_umbrella_extension import (
    patch_guard_source,
    patch_runtime_source,
    prepare_extension,
)

ROOT = Path(__file__).resolve().parents[1]


def write_json(path: Path, value: dict) -> None:
    """Write fixture input with predictable serialization."""
    path.write_text(json.dumps(value))


def fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    """Build a small 57-window six-block preparation without any native engine."""
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    original = campaign / "original"
    source = original / "source_snapshot"
    (source / "src/nyl12_covalent").mkdir(parents=True)
    (source / "scripts").mkdir()
    runtime = source / "src/nyl12_covalent/relay_umbrella_runtime.py"
    runtime.write_text(
        (ROOT / "src/nyl12_covalent/relay_umbrella_runtime.py").read_text()
    )
    script = source / "scripts/run_relay_umbrella.py"
    script.write_text("# retainedqueue\n")
    repair = tmp_path / "repair"
    (repair / "src/nyl12_covalent").mkdir(parents=True)
    (repair / "src/nyl12_covalent/relay_umbrella_diagnostics.py").write_text(
        "# validatedreader\n"
    )
    guard = tmp_path / "guard.py"
    guard.write_text("# noexecutionduringpreparation\n")
    windows, frozen = [], {}
    for index in range(57):
        label = f"window_{index}"
        stages = []
        for number in range(1, 9):
            name = f"{number:02d}_production_5ps"
            path = campaign / label / name
            path.mkdir(parents=True)
            mdin = path / "run.mdin"
            mdin.write_text(
                f"&cntrl irest=1,ntx=5,dt=0.0005,nstlim=10000,ig={1000 + index * 20 + number},/\n"
            )
            (path / "umbrella.RST").write_text("unchangedrestraint\n")
            frozen[str(mdin.relative_to(campaign))] = digest(mdin)
            stages.append(
                dict(
                    label=name,
                    input=str(mdin.relative_to(campaign)),
                    kind="production" if number >= 3 else "equil",
                    seed=1000 + index * 20 + number,
                )
            )
        windows.append(dict(label=label, stages=stages))
    prepared = campaign / "prepared.json"
    write_json(
        prepared,
        dict(
            planned_production_blocks_per_window=6,
            windows=windows,
            frozen_files_sha256=frozen,
        ),
    )
    plan = original / "plan.json"
    write_json(
        plan,
        dict(
            production_blocks=6,
            workers=57,
            initial_workers=28,
            prepared=str(prepared),
            prepared_sha256=digest(prepared),
            source_sha256={str(runtime): digest(runtime), str(script): digest(script)},
            runtime_sha256={},
        ),
    )
    plan.with_suffix(".sha256").write_text(digest(plan))
    return plan, campaign / "continuation_60ps", guard, repair


def test_patches_are_strict_and_guard_loads_plan_before_preparation() -> None:
    runtime = patch_runtime_source(
        (ROOT / "src/nyl12_covalent/relay_umbrella_runtime.py").read_text()
    )
    assert 'Path(plan["stage_plan_paths"][spec["label"]])' in runtime
    assert 'limit <= plan["production_blocks"] <= 12' in runtime
    assert 'base = Path(plan["campaign"])' in runtime
    compile(runtime, "runtime", "exec")
    guardian = patch_guard_source(
        (ROOT / "scripts/guard_relay_umbrella.py").read_text()
    )
    assert guardian.index("plan = json.loads(args.plan.read_text())") < guardian.index(
        'cfg = json.loads(Path(plan["prepared"])'
    )
    compile(guardian, "guardian", "exec")
    with pytest.raises(ValueError, match="anchor"):
        patch_runtime_source("changedunsupportedsource")
    with pytest.raises(ValueError, match="anchor"):
        patch_guard_source("changedunsupportedsource")


def test_prepares_exact_continuation_with_only_new_seed_in_inputs(
    tmp_path: Path,
) -> None:
    args = fixture(tmp_path)
    parent, output, guard, repair = args
    parent_hash = digest(parent)
    original_cfg = Path(json.loads(parent.read_text())["prepared"])
    cfg_hash = digest(original_cfg)
    result = prepare_extension(*args)
    assert result["windows"] == 57 and result["new_stage_paths"] == 684
    assert result["native_launched"] is False
    assert digest(parent) == parent_hash and digest(original_cfg) == cfg_hash
    new = json.loads((output / "prepared.json").read_text())
    plan = json.loads((output / "plan.json").read_text())
    assert plan["initial_workers"] == 57 and plan["production_blocks"] == 12
    assert plan["stage_plan_paths"]["08_production_5ps"] == str(parent)
    assert plan["stage_plan_paths"]["09_production_5ps"] == str(output / "plan.json")
    seeds = []
    for row in new["windows"]:
        for stage in row["stages"][8:]:
            mdin = original_cfg.parent / stage["input"]
            assert f"ig={stage['seed']}" in mdin.read_text()
            assert stage["starting_restart"].endswith(
                f"{int(stage['label'][:2]) - 1:02d}_production_5ps/restart.ncrst"
            )
            seeds.append(stage["seed"])
    assert len(seeds) == len(set(seeds)) == 342
    with pytest.raises(FileExistsError):
        prepare_extension(*args)


def test_changed_original_fails_before_new_stage_creation(tmp_path: Path) -> None:
    parent, output, guard, repair = fixture(tmp_path)
    campaign = Path(json.loads(parent.read_text())["prepared"]).parent
    (campaign / "window_0/08_production_5ps/run.mdin").write_text("changed")
    with pytest.raises(ValueError, match="frozen"):
        prepare_extension(parent, output, guard, repair)
    assert not (campaign / "window_0/09_production_5ps").exists()


def test_existing_unverified_source_clone_is_rejected(tmp_path: Path) -> None:
    parent, output, guard, repair = fixture(tmp_path)
    (output / "source_snapshot/src/nyl12_covalent").mkdir(parents=True)
    (output / "source_snapshot/scripts").mkdir()
    (
        output / "source_snapshot/src/nyl12_covalent/relay_umbrella_runtime.py"
    ).write_text("changed")
    with pytest.raises(ValueError, match="clone mismatch"):
        prepare_extension(parent, output, guard, repair)
