"""Portable student MD input/protocol/provenance and no-fallback contracts."""

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import parmed
import pytest

TOOLS = Path(__file__).resolve().parents[1] / "student_handoff/nyl12_ash274_md/tools"
sys.path.insert(0, str(TOOLS))
import prepared_md_support as support  # noqa:E402
import run_prepared_md as runner  # noqa:E402

ROOT = Path(__file__).resolve().parents[1]


def mapped_topology() -> tuple[SimpleNamespace, list[dict]]:
    """Build four explicit protected-copy records without inferring residue offsets."""
    residues = []
    rows = []
    atoms = []
    for site in support.SITES:
        for sequence in (113, 156, 186, 233, 272, 274):
            resname = support.PROTECTED.get(
                sequence, "ASH" if sequence == 274 else "ASP"
            )
            names = ["CA"] + (
                ["H1", "H2"] if sequence == 233 else ["HD2"] if resname == "ASH" else []
            )
            members = [SimpleNamespace(name=n) for n in names]
            residues.append(SimpleNamespace(name=resname, atoms=members))
            atoms.extend(members)
            rows.append(
                {
                    "site": site,
                    "chain": site[0] if sequence < 233 else site[1],
                    "sequence_residue": sequence,
                    "structure_residue": sequence if sequence < 233 else sequence - 232,
                    "amber_residue": len(residues),
                    "resname": resname,
                }
            )
    return SimpleNamespace(
        residues=residues,
        atoms=atoms,
        coordinates=np.zeros((len(atoms), 3)),
        box=np.array([80.0, 80.0, 80.0, 90.0, 90.0, 90.0]),
    ), rows


def dataset(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, SimpleNamespace]:
    """Pin tiny fixture artifacts while mocking only the topology file parser."""
    top, rows = mapped_topology()
    monkeypatch.setattr(parmed, "load_file", lambda *a, **k: top)
    files = {}
    for key in ("topology", "restart", "preparation_acceptance"):
        p = tmp_path / f"{key}.fixture"
        p.write_text(key)
        files[key] = {"path": p.name, "sha256": support.digest(p)}
    receipt = tmp_path / "receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "process_complete": True,
                "restart_readable": True,
                "coordinates_sha256": files["restart"]["sha256"],
                "topology_sha256": files["topology"]["sha256"],
            }
        )
    )
    files["minimization_acceptance"] = {
        "path": receipt.name,
        "sha256": support.digest(receipt),
    }
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "datasets": {
                    "wt_asp272_ash274": {
                        **files,
                        "microstate": "ASP272/ASH274",
                        "restart_stage": "unrestrained_minimized",
                        "protected_residues": rows,
                    }
                }
            }
        )
    )
    home = tmp_path / "amber"
    (home / "lib").mkdir(parents=True)
    (home / "bin").mkdir()
    engine = home / "bin/pmemd.cuda"
    engine.write_text("fake CUDA binary, never executed")
    engine.chmod(0o750)
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "dataset_manifest": manifest.name,
                "dataset": "wt_asp272_ash274",
                "amber_home": str(home),
                "cuda_engine": str(engine),
                "output": "runs",
                "session_prefix": "student_md",
                "replicas": [
                    {"name": "replica_01", "seed": 100, "gpu": "0"},
                    {"name": "replica_02", "seed": 101, "gpu": "1"},
                ],
            }
        )
    )
    return manifest, config, top


def test_artifact_and_independent_seeds(tmp_path: Path) -> None:
    p = tmp_path / "file"
    p.write_text("original")
    record = {"path": "file", "sha256": support.digest(p)}
    assert support.artifact(record, tmp_path) == p
    p.write_text("changed")
    with pytest.raises(ValueError, match="SHA-256"):
        support.artifact(record, tmp_path)
    assert (
        len(
            support.replicas(
                {
                    "replicas": [
                        {"name": "r1", "seed": 10, "gpu": "0"},
                        {"name": "r2", "seed": 11, "gpu": "1"},
                    ]
                }
            )
        )
        == 2
    )
    with pytest.raises(ValueError, match="overlap"):
        support.replicas(
            {
                "replicas": [
                    {"name": "r1", "seed": 10, "gpu": "0"},
                    {"name": "r2", "seed": 10010, "gpu": "1"},
                ]
            }
        )


@pytest.mark.parametrize(
    "rows",
    [
        [],
        [{"name": "bad/name", "seed": 1, "gpu": "0"}],
        [{"name": "r", "seed": 0, "gpu": "0"}],
        [{"name": "r", "seed": 2147483647, "gpu": "0"}],
        [{"name": "r", "seed": 1, "gpu": "cuda"}],
        [{"name": "r", "seed": 1, "gpu": "0"}, {"name": "r", "seed": 2, "gpu": "0"}],
    ],
)
def test_invalid_replica_declarations(rows: list) -> None:
    with pytest.raises(ValueError):
        support.replicas({"replicas": rows})


def test_four_copy_mapping_and_protected_mutation_rejection() -> None:
    top, rows = mapped_topology()
    support.validate_protected(top, rows)
    with pytest.raises(ValueError, match="24mapped"):
        support.validate_protected(top, rows[:-1])
    rows[0]["chain"] = "B"
    with pytest.raises(ValueError, match="chain"):
        support.validate_protected(top, rows)
    rows[0]["chain"] = "A"
    top.residues[0].name = "ALA"
    with pytest.raises(ValueError, match="identity"):
        support.validate_protected(top, rows)
    top.residues[0].name = "PHE"
    top.residues[3].atoms.append(SimpleNamespace(name="H3"))
    with pytest.raises(ValueError, match="neutral"):
        support.validate_protected(top, rows)


def test_receipt_matches_exact_coordinates_and_microstate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, _, _ = dataset(tmp_path, monkeypatch)
    data, top, rst, hashes = support.load_dataset(manifest, "wt_asp272_ash274")
    assert len(data["protected_residues"]) == 24 and hashes[str(rst)] == support.digest(
        rst
    )
    doc = json.loads(manifest.read_text())
    entry = doc["datasets"]["wt_asp272_ash274"]
    receipt = tmp_path / entry["minimization_acceptance"]["path"]
    record = json.loads(receipt.read_text())
    record["coordinates_sha256"] = "wrong"
    receipt.write_text(json.dumps(record))
    entry["minimization_acceptance"]["sha256"] = support.digest(receipt)
    manifest.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="these exact"):
        support.load_dataset(manifest, "wt_asp272_ash274")


@pytest.mark.parametrize(
    "change",
    [
        {"restart_stage": "solvated"},
        {"microstate": "silent_other_state"},
        {"microstate": "ASP272/ASP274"},
    ],
)
def test_no_silent_start_or_microstate_substitution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change: dict
) -> None:
    manifest, _, _ = dataset(tmp_path, monkeypatch)
    doc = json.loads(manifest.read_text())
    doc["datasets"]["wt_asp272_ash274"].update(change)
    manifest.write_text(json.dumps(doc))
    with pytest.raises(ValueError):
        support.load_dataset(manifest, "wt_asp272_ash274")


def test_explicit_gpu_gate_has_no_cpu_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, _ = dataset(tmp_path, monkeypatch)
    cfg = json.loads(config.read_text())
    assert support.engines(cfg, tmp_path)[1].name == "pmemd.cuda"
    cpu = tmp_path / "sander"
    cpu.write_text("cpu")
    cpu.chmod(0o750)
    cfg["cuda_engine"] = str(cpu)
    with pytest.raises(ValueError, match="fallback"):
        support.engines(cfg, tmp_path)
    monkeypatch.setattr(
        support.subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(stdout="0, GPU model, GPU-testuuid\n"),
    )
    assert support.require_gpu("0")["gpu_index"] == "0"
    with pytest.raises(ValueError, match="unavailable"):
        support.require_gpu("1")


def test_plan_uses_exact_original_10ns_plus_90ns_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, _ = dataset(tmp_path, monkeypatch)
    out = runner.plan(config, ROOT)
    from nyl12_covalent.pilot_md_protocol import mdin_files
    from scripts.continue_model30_md import continuation_input

    expected = mdin_files(100)
    folder = out / "replica_01/planned_inputs"
    for name, text in expected.items():
        assert (folder / f"{name}.in").read_text() == text
    assert (folder / "06_production_90ns.in").read_text() == continuation_input(
        expected["05_production_10ns"], 50100
    )
    assert " ntr=0," in (folder / "06_production_90ns.in").read_text()
    assert not (out / "replica_01/job").exists()
    cfg = json.loads(config.read_text())
    engine = Path(cfg["cuda_engine"])
    engine.write_text("mutated")
    with pytest.raises(ValueError, match="changed"):
        runner.frozen(config, ROOT)


def test_launch_uses_existing_tmux_logger_and_explicit_devices(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, _ = dataset(tmp_path, monkeypatch)
    out = runner.plan(config, ROOT)
    import nyl12_covalent.jobs as jobs

    calls = []
    monkeypatch.setattr(jobs, "launch_job", lambda *a, **k: calls.append((a, k)))
    monkeypatch.setattr(runner, "require_gpu", lambda gpu: {"gpu_index": gpu})
    runner.launch(config, ROOT)
    assert (
        len(calls) == 2 and calls[0][1]["env_overrides"]["CUDA_VISIBLE_DEVICES"] == "0"
    )
    assert calls[1][1]["env_overrides"]["CUDA_VISIBLE_DEVICES"] == "1"
    assert calls[0][0][1] == out / "replica_01/job"
    assert "execute" in calls[0][0][3]


def test_execute_reuses_existing_stage_runners_and_stops_on_gpu_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, _ = dataset(tmp_path, monkeypatch)
    out = runner.plan(config, ROOT)
    monkeypatch.setattr(runner, "require_gpu", lambda gpu: {"gpu_index": gpu})
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "1")
    with pytest.raises(ValueError, match="differs"):
        runner.execute(config, ROOT, "replica_01")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    calls = []
    monkeypatch.setattr(
        runner.subprocess, "run", lambda argv, **kwargs: calls.append(argv)
    )
    runner.execute(config, ROOT, "replica_01")
    assert calls[0][1] == str(ROOT / "scripts/run_model30_md_replica.py")
    assert calls[1][1] == str(ROOT / "scripts/continue_model30_md.py")
    assert (out / "replica_01/100ns_process_complete.json").exists()


def test_explicit_engine_environment_is_pinned_to_plan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, _ = dataset(tmp_path, monkeypatch)
    cfg = json.loads(config.read_text())
    home = cfg["amber_home"]
    engine = cfg["cuda_engine"]
    cfg.update(amber_home="$STUDENT_TEST_AMBERHOME", cuda_engine="$STUDENT_TEST_CUDA")
    config.write_text(json.dumps(cfg))
    monkeypatch.delenv("STUDENT_TEST_CUDA", raising=False)
    monkeypatch.setenv("STUDENT_TEST_AMBERHOME", home)
    with pytest.raises(ValueError, match="Unresolved"):
        runner.plan(config, ROOT)
    monkeypatch.setenv("STUDENT_TEST_CUDA", engine)
    runner.plan(config, ROOT)
    other = Path(engine).parent / "pmemd.cuda_other"
    other.write_text("other GPU binary")
    other.chmod(0o750)
    monkeypatch.setenv("STUDENT_TEST_CUDA", str(other))
    with pytest.raises(ValueError, match="environment differs"):
        runner.frozen(config, ROOT)


def test_all_gpus_are_checked_before_any_replica_launch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, _ = dataset(tmp_path, monkeypatch)
    runner.plan(config, ROOT)
    import nyl12_covalent.jobs as jobs

    launched = []
    monkeypatch.setattr(jobs, "launch_job", lambda *a, **k: launched.append(a))

    def check(gpu: str) -> dict:
        if gpu == "1":
            raise ValueError("GPU1 unavailable")
        return {"gpu_index": gpu}

    monkeypatch.setattr(runner, "require_gpu", check)
    with pytest.raises(ValueError, match="GPU1"):
        runner.launch(config, ROOT)
    assert not launched
