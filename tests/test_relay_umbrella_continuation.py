"""Verify source-preserving 60-to100ps continuation preparation."""

from pathlib import Path
import json
import re

import pytest

from nyl12_covalent.relay_umbrella_continuation import (
    digest,
    patch_block_limit,
    prepare_continuation,
)


def write(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


@pytest.fixture
def campaign(tmp_path: Path) -> tuple[Path, Path]:
    base = tmp_path / "campaign"
    old = base / "continuation60"
    source = old / "source_snapshot/src/nyl12_covalent"
    source.mkdir(parents=True)
    runtime = source / "relay_umbrella_runtime.py"
    runtime.write_text(
        'if not 1 <= limit <= plan["production_blocks"] <= 12:\n    raise ValueError()\n'
    )
    engine = base / "sander"
    engine.write_text("validated host binary")
    rows = []
    frozen = {}
    for index in range(57):
        folder = base / f"w{index:02d}"
        stages = []
        for number in range(3, 15):
            label = f"{number:02d}_production_5ps"
            stage = folder / label
            stage.mkdir(parents=True)
            mdin = stage / "run.mdin"
            mdin.write_text(
                f'&cntrl irest=1,ntx=5,nstlim=10000,ig={index * 20 + number}, /\n&wt type="END" /\nDISANG=umbrella.RST\n'
            )
            bias = stage / "umbrella.RST"
            bias.write_text("same bias25.0")
            for p in [mdin, bias]:
                frozen[str(p.relative_to(base))] = digest(p)
            stages.append(
                {
                    "label": label,
                    "input": str(mdin.relative_to(base)),
                    "kind": "production",
                    "seed": index * 20 + number,
                    "steps": 10000,
                }
            )
        checkpoint = folder / "14_production_5ps/restart.ncrst"
        checkpoint.write_text("accepted coordinates velocities box time")
        write(
            folder / "segment_receipt_limit_12.json",
            {
                "status": "accepted",
                "production_blocks": 12,
                "checkpoint_sha256": digest(checkpoint),
            },
        )
        write(
            folder / "14_production_5ps/acceptance.json",
            {"accepted": True, "checkpoint_sha256": digest(checkpoint)},
        )
        write(folder / "14_production_5ps/process.json", {"exit_status": 0})
        rows.append({"label": folder.name, "stages": stages})
    prepared = old / "prepared.json"
    write(
        prepared,
        {
            "windows": rows,
            "planned_production_blocks_per_window": 12,
            "frozen_files_sha256": frozen,
            "system": {"unchanged": "chemistry"},
        },
    )
    plan = old / "plan.json"
    write(
        plan,
        {
            "workers": 57,
            "production_blocks": 12,
            "campaign": str(base),
            "prepared": str(prepared),
            "prepared_sha256": digest(prepared),
            "sander": str(engine),
            "runtime_sha256": {str(engine): digest(engine)},
            "source_sha256": {str(runtime): digest(runtime)},
            "stage_plan_paths": {
                f"{i:02d}_production_5ps": str(plan) for i in range(3, 15)
            },
            "environment": {},
        },
    )
    plan.with_suffix(".sha256").write_text(digest(plan) + "\n")
    return plan, base / "continuation100"


def test_full57_prefix_and_only_ig_changes(campaign: tuple[Path, Path]) -> None:
    old, out = campaign
    before = old.read_bytes()
    oldplan = json.loads(before)
    original = json.loads(Path(oldplan["prepared"]).read_text())
    result = prepare_continuation(old, out)
    plan = json.loads((out / "plan.json").read_text())
    new = json.loads((out / "prepared.json").read_text())
    assert (
        old.read_bytes() == before
        and result["windows"] == 57
        and result["workers"] == 57
    )
    assert plan["stage_plan_paths"]["14_production_5ps"] == str(old)
    assert plan["stage_plan_paths"]["15_production_5ps"] == str(out / "plan.json")
    seeds = []
    for first, last in zip(original["windows"], new["windows"], strict=True):
        assert first["stages"] == last["stages"][:12]
        base = Path(oldplan["campaign"])
        template = (base / first["stages"][-1]["input"]).read_text()
        for stage in last["stages"][12:]:
            assert re.sub(
                r"ig=\d+", "ig=SEED", (base / stage["input"]).read_text()
            ) == re.sub(r"ig=\d+", "ig=SEED", template)
            assert (base / stage["restraint"]).read_text() == "same bias25.0"
            seeds.append(stage["seed"])
    assert len(seeds) == len(set(seeds)) == 456
    assert new["system"] == original["system"]
    assert plan["environment"] == {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "CUDA_VISIBLE_DEVICES": "",
    }
    with pytest.raises(ValueError):
        prepare_continuation(old, out)


@pytest.mark.parametrize(
    "what",
    ["checkpoint", "frozen", "source", "runtime", "future", "acceptance", "parenthash"],
)
def test_invalid_prefix_never_prepares(campaign: tuple[Path, Path], what: str) -> None:
    old, out = campaign
    plan = json.loads(old.read_text())
    base = Path(plan["campaign"])
    if what == "checkpoint":
        (base / "w00/14_production_5ps/restart.ncrst").write_text("changed")
    elif what == "frozen":
        (base / "w00/14_production_5ps/run.mdin").write_text("changed")
    elif what == "source":
        Path(next(iter(plan["source_sha256"]))).write_text("changed")
    elif what == "runtime":
        Path(plan["sander"]).write_text("changed")
    elif what == "future":
        (base / "w00/15_production_5ps").mkdir()
    elif what == "acceptance":
        write(base / "w00/14_production_5ps/process.json", {"exit_status": 1})
    elif what == "parenthash":
        old.write_text("{}")
    with pytest.raises((ValueError, FileExistsError)):
        prepare_continuation(old, out)
    assert not out.exists()


@pytest.mark.parametrize(
    "text,target",
    [
        ("no anchor", 20),
        ('if not 1 <= limit <= plan["production_blocks"] <= 12:' * 2, 20),
        ("anything", 12),
    ],
)
def test_source_patch_strict(text: str, target: int) -> None:
    with pytest.raises(ValueError):
        patch_block_limit(text, target)


def test_budget_refuses_extra_workers(campaign: tuple[Path, Path]) -> None:
    with pytest.raises(ValueError):
        prepare_continuation(*campaign, workers=58)
