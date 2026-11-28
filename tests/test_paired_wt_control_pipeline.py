"""Check paired WT staging, resource barriers and stock command preservation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import pytest

CAMPAIGN = (
    Path(__file__).resolve().parents[1]
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/wt_paired_precursor_controls_20261127"
)


def load(name: str) -> ModuleType:
    sys.path.insert(0, str(CAMPAIGN))
    spec = importlib.util.spec_from_file_location(
        "paired_test_" + name, CAMPAIGN / (name + ".py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def queue_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[ModuleType, dict, list, list]:
    m = load("queue_worker")
    monkeypatch.setattr(m, "BASE", tmp_path)
    monkeypatch.setattr(sys, "argv", ["queue", "--gpu", "GPU2"])
    record = tmp_path / "output/boltz_results_WT/processed/records/WT.json"
    write(record, {})
    native = tmp_path / "native.py"
    native.write_text("unchanged native feature source")
    input_path = tmp_path / "WT.yaml"
    input_path.write_text("stock WT input")
    paths = []
    for i in range(3):
        job = tmp_path / f"old_job_{i}"
        write(job / "status.json", {"exit_status": None})
        paths.append(str(job))
    csvs = {}
    for part in ("alpha", "beta"):
        path = tmp_path / "msas" / f"{part}_paired.csv"
        path.parent.mkdir(exist_ok=True)
        path.write_text("sequence,key\nA,-1\nC,1\n")
        csvs[path.name] = sha(path)
    receipt = tmp_path / "msas/precursor_pairing_receipt.json"
    write(receipt, {"output_sha256": csvs})
    checkpoint = tmp_path / "checkpoint"
    checkpoint.write_text("validated checkpoint placeholder")
    case = {
        "case": "WT",
        "gpu_uuid": "GPU2",
        "wait_for_existing_jobs": paths,
        "remote_input": str(input_path),
        "remote_output": str(tmp_path / "output"),
        "remote_job": str(tmp_path / "new_job"),
        "cpu_affinity": "16-23",
        "seed": 20261126,
    }
    write(
        tmp_path / "campaign.json",
        {"cases": [case, {"case": "other", "gpu_uuid": "GPU3"}]},
    )
    write(
        tmp_path / "cpu_acceptance.json",
        {
            "status": "accepted",
            "stock_native_entry_only": True,
            "paired_source_receipt_sha256": sha(receipt),
            "native_featurizer_source_path": str(native),
            "native_featurizer_source_sha256": sha(native),
            "models": [
                {
                    "case": "WT",
                    "input_sha256": sha(input_path),
                    "record_sha256": sha(record),
                }
            ],
        },
    )
    write(
        tmp_path / "baseline_boltz_status.json",
        {
            "command": [
                "taskset",
                "-c",
                "16-23",
                "/home/jpc/code/miniforge3/envs/boltz-env/bin/boltz",
                "predict",
                "old.yaml",
                "--checkpoint",
                str(checkpoint),
                "--seed",
                "1",
                "--out_dir",
                "old",
            ],
            "env_overrides": {"OMP_NUM_THREADS": "2", "CUDA_VISIBLE_DEVICES": "old"},
        },
    )
    real_sha = m.sha
    monkeypatch.setattr(
        m,
        "sha",
        lambda p: (
            "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"
            if p == checkpoint
            else real_sha(p)
        ),
    )
    monkeypatch.setattr(
        m.shutil, "disk_usage", lambda p: SimpleNamespace(free=11 * 1024**3)
    )
    sleeps = []

    def advance(seconds: int) -> None:
        sleeps.append(seconds)
        for path in paths:
            write(Path(path) / "status.json", {"exit_status": 0})

    monkeypatch.setattr(m.time, "sleep", advance)
    monkeypatch.setattr(m, "require_idle", lambda gpu: None)
    commands = []

    def run(argv: list[str], **kwargs: object) -> SimpleNamespace:
        commands.append((argv, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(m.subprocess, "run", run)
    return m, case, sleeps, commands


def test_queue_waits_existing_three_and_keeps_stock_entry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    m, case, sleeps, commands = queue_fixture(tmp_path, monkeypatch)
    m.main()
    assert sleeps == [60]
    assert len(commands) == 1
    argv, options = commands[0]
    assert argv[3].endswith("/bin/boltz")
    assert "isolated_boltz_entry.py" not in " ".join(argv)
    assert options["env"]["CUDA_VISIBLE_DEVICES"] == "GPU2"
    status = json.loads((Path(case["remote_job"]) / "status.json").read_text())
    assert status["exit_status"] == 0
    assert "scientific_acceptance" in status


def test_queue_stops_on_existing_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    m, case, _, commands = queue_fixture(tmp_path, monkeypatch)
    write(Path(case["wait_for_existing_jobs"][0]) / "status.json", {"exit_status": 1})
    with pytest.raises(RuntimeError, match="Existing repair queue failed"):
        m.main()
    assert not commands


def test_queue_detects_csv_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    m, _, _, commands = queue_fixture(tmp_path, monkeypatch)
    (tmp_path / "msas/alpha_paired.csv").write_text("changed")
    with pytest.raises(ValueError, match="paired CSV changed"):
        m.main()
    assert not commands


def test_queue_records_native_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    m, case, _, _ = queue_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(
        m.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=7)
    )
    with pytest.raises(RuntimeError, match="exit7"):
        m.main()
    assert (
        json.loads((Path(case["remote_job"]) / "status.json").read_text())[
            "exit_status"
        ]
        == 7
    )


def test_idle_guard_waits_without_signaling(monkeypatch: pytest.MonkeyPatch) -> None:
    m = load("queue_worker")
    results = iter(["GPU2,42\n", "GPU3,55\n"])
    monkeypatch.setattr(
        m.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout=next(results))
    )
    waits = []
    monkeypatch.setattr(m.time, "sleep", waits.append)
    m.require_idle("GPU2")
    assert waits == [60]


def test_orchestrator_full_qualified_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    m = load("orchestrate_controls")
    root = tmp_path / "repo"
    base = tmp_path / "campaign"
    base.mkdir()
    package = root / "src/nyl12_covalent"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("# frozen source")
    source = base / "source_acceptance.json"
    write(source, {"status": "accepted", "precursor_a3m": "valid.a3m"})
    monkeypatch.setattr(m, "BASE", base)
    monkeypatch.setattr(m, "ROOT", root)
    monkeypatch.setattr(m, "SSH", ["ssh", "host"])
    monkeypatch.setattr(
        sys, "argv", ["orchestrator", "--source-acceptance", str(source)]
    )
    monkeypatch.setattr(m, "journal", lambda *args: None)
    commands = []

    def run(argv: list[str], **kwargs: object) -> SimpleNamespace:
        commands.append(argv)
        if "--precursor-a3m" in argv:
            write(
                base / "campaign.json",
                {"cases": [{"gpu_uuid": "GPU2"}, {"gpu_uuid": "GPU3"}]},
            )
            write(
                base / "cpu_acceptance.json", {"status": "accepted", "models": [{}, {}]}
            )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(m.subprocess, "run", run)
    m.main()
    state = json.loads((base / "pipeline_state.json").read_text())
    assert state["stock_native_entry_only"] and state["target_phi_sign"] == 1
    assert state["status"].endswith("postprocessing_active")
    assert len([a for a in commands if a[:2] == ["tmux", "new-session"]]) == 2
    assert (
        len(
            [
                a
                for a in commands
                if "nyl12_paired_WT_gpu2" in a or "nyl12_paired_WT_gpu3" in a
            ]
        )
        == 2
    )
    assert (base / "source_hashes.json").is_file()


def test_stage_preserves_wt_graph_and_uses_one_precursor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import shutil

    m = load("stage_campaign")
    pilot = tmp_path / "pilot"
    base = pilot / "wt_paired_precursor_controls_20261127"
    base.mkdir(parents=True)
    reference = CAMPAIGN.parent
    for relative in (
        "boltz_masked_loop95_designs_20261127/campaign.json",
        "boltz_masked_loop95_designs_20261127/template_alignment_repair01/campaign.json",
        "boltz_masked_loop95_designs_20261127/baseline_boltz_status.json",
        "boltz_masked_loop95_designs_20261127/target_pose_convention.json",
        "prediction_preparation/ligand_graph_manifest.json",
    ):
        target = pilot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(reference / relative, target)
    wt = json.loads(
        (reference / "boltz_masked_loop95_designs_20261127/campaign.json").read_text()
    )["cases"][0]
    a, b = wt["sequences"]["alpha"], wt["sequences"]["beta"]
    source = tmp_path / "single_precursor.a3m"
    source.write_text(
        f"#328\t1\n>query\n{a + b}\n>one_actual_homolog\nC{a[1:]}A{b[1:]}\n"
    )
    provenance = tmp_path / "source_acceptance.json"
    write(provenance, {"status": "accepted", "precursor_a3m_sha256": sha(source)})
    monkeypatch.setattr(m, "BASE", base)
    monkeypatch.setattr(m, "PILOT", pilot)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "stage",
            "--precursor-a3m",
            str(source),
            "--source-provenance",
            str(provenance),
        ],
    )
    m.main()
    plan = json.loads((base / "campaign.json").read_text())
    assert plan["prediction_budget"] == 2
    assert "mask" not in plan and "negative_gate" not in plan
    assert plan["target_phi_sign"] == 1
    assert plan["paired_precursor_msa"]["paired_homolog_row_count"] == 1
    for case in plan["cases"]:
        data = json.loads(Path(case["input"]).read_text())
        assert set(data) == {"version", "sequences"}
        assert len(case["wait_for_existing_jobs"]) == 3
        assert data["sequences"][0]["protein"]["sequence"] == a
        assert data["sequences"][1]["protein"]["sequence"] == b
        assert data["sequences"][2]["ligand"]["id"] == list("IJKL")
    with pytest.raises(FileExistsError):
        m.main()


def test_paired_publication_exposes_only_verified_subset(tmp_path: Path) -> None:
    m = load("finalize_publication")
    output = tmp_path / "inspection_export_20261127"
    output.mkdir()
    query = output / "query.yaml"
    query.write_text("verified two WT inputs")
    manifest = output / "MANIFEST.sha256"
    manifest.write_text(f"{sha(query)}  query.yaml\n")
    write(
        tmp_path / "publication_acceptance.json",
        {
            "status": "READY_TO_PULL",
            "matched_cases": 2,
            "target_phi_sign": 1,
            "manifest_sha256": sha(manifest),
        },
    )
    (tmp_path / "pull_after_acceptance.block").write_text(
        "# paired WT unique block\necho valid\n"
    )
    pull = tmp_path / "pull.sh"
    pull.write_text("#!/bin/bash\necho existing\n")
    m.verify_and_wire(tmp_path, pull)
    assert "# paired WT unique block" in pull.read_text()
    assert (tmp_path / "filtered_pull_acceptance.json").is_file()
    with pytest.raises(FileExistsError):
        m.verify_and_wire(tmp_path, pull)
    write(tmp_path / "publication_acceptance.json", {"status": "not accepted"})
    with pytest.raises(ValueError):
        m.verify_and_wire(tmp_path, pull)


def test_paired_publication_source_binding(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    m = load("finalize_publication")
    source = tmp_path / "publish_campaign.py"
    source.write_text("pass")
    write(
        tmp_path / "followthrough_source_binding.json",
        {"publisher_sha256": sha(source), "qualification_status": "accepted"},
    )
    monkeypatch.setattr(m, "BASE", tmp_path)
    monkeypatch.setattr(m, "ROOT", tmp_path)
    calls = []
    monkeypatch.setattr(m.subprocess, "run", lambda *a, **k: calls.append(a))
    monkeypatch.setattr(m, "verify_and_wire", lambda *a: calls.append(a))
    m.main()
    assert len(calls) == 2
    source.write_text("changed")
    with pytest.raises(ValueError):
        m.main()


def test_publication_waits_for_explicit_schema_qualification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    m = load("finalize_publication")
    source = tmp_path / "publish_campaign.py"
    source.write_text("pass")
    binding = tmp_path / "followthrough_source_binding.json"
    write(
        binding,
        {
            "publisher_sha256": sha(source),
            "qualification_status": "pending_publication_schema_v2",
        },
    )
    monkeypatch.setattr(m, "BASE", tmp_path)
    monkeypatch.setattr(m, "ROOT", tmp_path)
    calls = []
    waits = []

    def qualify(seconds: int) -> None:
        assert not calls
        waits.append(seconds)
        write(
            binding,
            {"publisher_sha256": sha(source), "qualification_status": "accepted"},
        )

    monkeypatch.setattr(m.time, "sleep", qualify)
    monkeypatch.setattr(m.subprocess, "run", lambda *a, **k: calls.append(a))
    monkeypatch.setattr(m, "verify_and_wire", lambda *a: calls.append(a))
    m.main()
    assert waits == [15] and len(calls) == 2
    write(binding, {"publisher_sha256": sha(source), "qualification_status": "invalid"})
    with pytest.raises(ValueError, match="Qualified publication source required"):
        m.main()
