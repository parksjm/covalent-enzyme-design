"""Synthetic process checks for the frozen outer comparison controller."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import run_original_matched_comparison as module


def setup_config(tmp_path: Path) -> dict:
    """Create synthetic immutable sources and completed search statuses."""
    source = tmp_path / "source"
    source.write_text("frozen")
    prediction = tmp_path / "plan.json"
    prediction.write_text(
        json.dumps(
            {
                "seeds": list(module.SEEDS),
                "candidates": [
                    {"candidate": "WT", "gpu_uuid": "gpu0"},
                    {"candidate": "relay_prot101_06", "gpu_uuid": "gpu0"},
                    {"candidate": "relay_prot101_03", "gpu_uuid": "gpu1"},
                ],
            }
        )
    )
    amber = tmp_path / "amber.json"
    amber.write_text(
        json.dumps(
            {
                "input_sha256": {str(source): module.digest(source)},
                "relaxation_config": {"gpu": "gpu0"},
            }
        )
    )
    for name in ("relay_prot101_03", "relay_prot101_06"):
        folder = tmp_path / "jobs" / f"msa_{name}"
        folder.mkdir(parents=True)
        (folder / "status.json").write_text('{"exit_status":0}')
    return {
        "repository_root": str(tmp_path),
        "campaign_root": str(tmp_path),
        "python": "/fake/python",
        "source_snapshot": str(tmp_path / "snapshot"),
        "prediction_plan": str(prediction),
        "amber_template": str(amber),
        "source_sha256": {str(source): module.digest(source)},
    }


@pytest.mark.parametrize("failure", [None, "native", "amber", "missing"])
def test_full_controller(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    """All native calls precede Amber; independent failures retain six-case census."""
    config = setup_config(tmp_path)
    monkeypatch.setattr(module, "journal", lambda *args: None)
    calls = []

    def fake_invoke(config: dict, argv: list[str], name: str) -> int:
        calls.append(name)
        if name.startswith("prepare_"):
            candidate = argv[argv.index("--candidate") + 1]
            for seed in module.SEEDS:
                work = tmp_path / "prediction/cases" / module.case_name(candidate, seed)
                work.mkdir(parents=True)
                (work / "case.json").write_text("{}")
        if (
            name.startswith("predict_")
            and failure == "native"
            and "WT_seed_20261126" in name
        ):
            return 1
        if name.startswith("amber_"):
            assert sum(x.startswith("predict_") for x in calls) == 6
            case = argv[-1]
            if failure == "missing" and case == "WT_seed_20261126":
                return 1
            out = tmp_path / "relaxation"
            out.mkdir(exist_ok=True)
            failed = failure == "amber" and case == "WT_seed_20261126"
            (out / f"failure_census_{case}.json").write_text(
                json.dumps(
                    {
                        "cases": {
                            case: {
                                "status": "failed" if failed else "completed",
                                "positive_ab_before": 2,
                                "positive_ab_after": 4,
                            }
                        }
                    }
                )
            )
            return int(failed)
        return 0

    monkeypatch.setattr(module, "invoke", fake_invoke)
    report = module.run(config)
    assert len(report["cases"]) == 6
    assert report["failure_count"] == int(failure is not None)
    assert (tmp_path / "controller/comparison.csv").is_file()
    for worker in range(2):
        p = json.loads(
            (tmp_path / f"controller/worker{worker}_prediction_plan.json").read_text()
        )
        assert {x["gpu_uuid"] for x in p["candidates"]} == {f"gpu{worker}"}


def test_search_wait_and_fail(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Pending statuses wait once; failed searches stop immediately."""
    config = setup_config(tmp_path)
    path = tmp_path / "jobs/msa_relay_prot101_03/status.json"
    path.write_text('{"exit_status":null}')
    monkeypatch.setattr(
        module.time, "sleep", lambda seconds: path.write_text('{"exit_status":0}')
    )
    module.wait_searches(config)
    path.write_text('{"exit_status":7}')
    with pytest.raises(RuntimeError, match="Search failed"):
        module.wait_searches(config)


def test_pins_and_immutable_outputs(tmp_path: Path) -> None:
    """Invalid pins and overwritten evidence fail explicitly."""
    with pytest.raises(ValueError):
        module.verify_pins({})
    path = tmp_path / "file"
    path.write_text("test")
    with pytest.raises(ValueError):
        module.verify_pins({str(path): "wrong"})
    with pytest.raises(FileExistsError):
        module.write_json(path, {})


def test_invoke_and_journal(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A stage captures command/status and makes required journal subprocesses."""
    config = setup_config(tmp_path)
    (tmp_path / "controller").mkdir()
    seen = []

    def process(argv: list[str], **kwargs: object) -> SimpleNamespace:
        seen.append(argv)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(module.subprocess, "run", process)
    assert module.invoke(config, ["/fake/command"], "stage") == 0
    assert len(seen) == 3
    assert (
        json.loads((tmp_path / "controller/stage.json").read_text())["exit_status"] == 0
    )


def test_invalid_absolute_config(tmp_path: Path) -> None:
    """Relative roots are never accepted."""
    config = setup_config(tmp_path)
    config["repository_root"] = "relative"
    with pytest.raises(ValueError, match="absolute"):
        module.run(config)
