"""Test the run-specific host launcher without reserving real GPUs or tmux."""

import importlib.util
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from nyl12_covalent import jobs, prepared_md_queue_plan

SOURCE = Path(__file__).resolve().parents[1] / (
    "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/launch_host.py"
)


def load_launcher(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> object:
    spec = importlib.util.spec_from_file_location("design_md_host", SOURCE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "__file__", str(tmp_path / "launch_host.py"))
    engine = tmp_path / "pmemd.cuda_SPFP"
    engine.write_bytes(b"validated CUDA build")
    monkeypatch.setattr(
        module,
        "PROFILES",
        {
            "violet": {
                "engine": str(engine),
                "sha256": hashlib.sha256(engine.read_bytes()).hexdigest(),
                "config_prefix": "",
                "queue_suffix": "",
                "receipt": "",
            },
        },
    )
    (tmp_path / "configs").mkdir()
    for label in ("WT", "broader_mask07", "focused_mask03", "V111C"):
        (tmp_path / "configs" / f"{label}.json").write_text(
            json.dumps({"cuda_engine": str(engine)})
        )
    monkeypatch.setattr(module.platform, "node", lambda: "violet.example")
    monkeypatch.setattr(
        module.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout="GPU0/GPU1")
    )
    return module


@pytest.mark.parametrize("prepare", [False, True])
def test_explicit_host_and_two_serial_queues(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    prepare: bool,
) -> None:
    module = load_launcher(tmp_path, monkeypatch)
    plans, launches = [], []
    monkeypatch.setattr(
        prepared_md_queue_plan,
        "runtime",
        lambda *a: SimpleNamespace(plan=lambda *a: plans.append(a)),
    )
    monkeypatch.setattr(jobs, "launch_job", lambda *a, **k: launches.append((a, k)))
    argv = ["launch_host", "--host", "violet"]
    if prepare:
        argv.append("--prepare-plans")
    monkeypatch.setattr(module.sys, "argv", argv)
    module.main()
    assert len(plans) == (4 if prepare else 0)
    assert len(launches) == 2
    assert {row[1]["env_overrides"]["CUDA_VISIBLE_DEVICES"] for row in launches} == {
        "0",
        "1",
    }
    assert (tmp_path / "queues/violet_environment.json").is_file()
    with pytest.raises(FileExistsError):
        module.main()


def test_wrong_execution_host_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_launcher(tmp_path, monkeypatch)
    monkeypatch.setattr(module.sys, "argv", ["launch", "--host", "azalea"])
    with pytest.raises(ValueError, match="execution workstation"):
        module.main()


@pytest.mark.parametrize("wrong_hash", [False, True])
def test_wrong_host_binary_rejected_before_planning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    wrong_hash: bool,
) -> None:
    module = load_launcher(tmp_path, monkeypatch)
    monkeypatch.setattr(module.sys, "argv", ["launch", "--host", "violet"])
    if wrong_hash:
        (tmp_path / "pmemd.cuda_SPFP").write_bytes(b"another build")
    else:
        other = tmp_path / "generic_other_cuda"
        other.write_bytes(b"wrong architecture")
        (tmp_path / "configs/WT.json").write_text(
            json.dumps({"cuda_engine": str(other)})
        )
    with pytest.raises(ValueError, match="hash differs|wrong workstation"):
        module.main()
