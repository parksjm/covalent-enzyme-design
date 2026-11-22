"""Verify the explicit-runtime adapter to the established Boltz worker."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import queue
import sys
from unittest.mock import patch

import pytest

from scripts.run_design_boltz_predictions import (
    Job,
    canonical_dir,
    command,
    merge_output,
    parse_args,
    read_jobs,
    run_job,
)


def _args(tmp_path: Path) -> argparse.Namespace:
    return argparse.Namespace(
        output_root=tmp_path / "out",
        source_dir=tmp_path / "source",
        cache=tmp_path / "cache",
        boltz_python=Path("/existing/env/bin/python"),
    )


def test_native_runtime_and_seeded_command(tmp_path: Path) -> None:
    """Keep the existing science settings while using the exact supplied runtime."""
    job = Job("pilot", "variant_01", 3, 2026112204, tmp_path / "x.yaml")
    args = _args(tmp_path)
    argv = command(args, job)
    assert argv[:3] == [
        str(args.boltz_python),
        str(args.source_dir / "run_boltz_cli.py"),
        "predict",
    ]
    for flag, value in [
        ("--seed", "2026112204"),
        ("--sampling_steps", "200"),
        ("--recycling_steps", "3"),
        ("--diffusion_samples", "1"),
    ]:
        assert argv[argv.index(flag) + 1] == value
    assert "pixi" not in argv and "--use_potentials" in argv
    assert "--affinity_mw_correction" not in argv


def test_job_manifest_rejects_duplicate_or_missing_inputs(tmp_path: Path) -> None:
    """A repeated sequence/seed index or absent input must stop dispatch."""
    path = tmp_path / "jobs.csv"
    template = tmp_path / "WT.yaml"
    template.write_text("{}")
    header = "campaign,sequence_name,model_index,seed,yaml_path\n"
    row = f"pilot,WT,0,1,{template}\n"
    path.write_text(header + row)
    assert len(read_jobs(path)) == 1
    path.write_text(header + row + row)
    with pytest.raises(ValueError, match="Duplicate"):
        read_jobs(path)
    template.unlink()
    with pytest.raises(FileNotFoundError):
        read_jobs(path)


def test_dispatch_merges_exact_outputs_and_returns_gpu(tmp_path: Path) -> None:
    """Successful native output is exported once with its confidence metadata."""
    args = _args(tmp_path)
    job = Job("pilot", "variant_01", 2, 3, tmp_path / "x.yaml")
    gpus: queue.Queue[str] = queue.Queue()
    gpus.put("1")

    def fake_run(argv: list[str], **kwargs: object) -> None:
        assert kwargs["env"]["CUDA_VISIBLE_DEVICES"] == "1"
        out = (
            Path(argv[argv.index("--out_dir") + 1])
            / "boltz_results_variant_01/predictions/variant_01"
        )
        out.mkdir(parents=True)
        (out / "variant_01_model_0.pdb").write_text("EXACT_NATIVE_BYTES\n")
        (out / "confidence_variant_01_model_0.json").write_text(
            json.dumps({"confidence_score": 0.8})
        )

    with patch(
        "scripts.run_design_boltz_predictions.subprocess.run", side_effect=fake_run
    ) as native:
        run_job(args, job, gpus)
        run_job(args, job, gpus)
    assert native.call_count == 1
    assert gpus.get_nowait() == "1"
    assert (
        canonical_dir(args.output_root, job) / "variant_01_model_2.pdb"
    ).read_text() == "EXACT_NATIVE_BYTES\n"
    # A missing raw artifact must never be marked merged.
    missing = Job("pilot", "WT", 0, 4, tmp_path / "WT.yaml")
    with pytest.raises(FileNotFoundError):
        merge_output(args.output_root, missing)


def test_native_failure_propagates_and_releases_gpu(tmp_path: Path) -> None:
    """Dispatch must not swallow a failed scientific process or lose the GPU slot."""
    args = _args(tmp_path)
    job = Job("pilot", "WT", 0, 4, tmp_path / "WT.yaml")
    gpus: queue.Queue[str] = queue.Queue()
    gpus.put("0")
    with (
        patch(
            "scripts.run_design_boltz_predictions.subprocess.run",
            side_effect=RuntimeError("native failed"),
        ),
        pytest.raises(RuntimeError, match="native failed"),
    ):
        run_job(args, job, gpus)
    assert gpus.get_nowait() == "0"


def _write_one_job(tmp_path: Path) -> tuple[Path, Path]:
    yaml_path = tmp_path / "WT.yaml"
    yaml_path.write_text("{}")
    jobs_path = tmp_path / "jobs.csv"
    jobs_path.write_text(
        "campaign,sequence_name,model_index,seed,yaml_path\n"
        f"pilot,WT,0,2026112201,{yaml_path}\n"
    )
    return jobs_path, yaml_path


def _cli_args(tmp_path: Path, jobs_path: Path, gpu_devices: str = "0") -> list[str]:
    return [
        "run_design_boltz_predictions.py",
        "--jobs",
        str(jobs_path),
        "--output-root",
        str(tmp_path / "out"),
        "--source-dir",
        str(tmp_path / "source"),
        "--cache",
        str(tmp_path / "cache"),
        "--boltz-python",
        "/existing/env/bin/python",
        "--gpu-devices",
        gpu_devices,
    ]


def test_parse_args_accepts_explicit_runtime_and_devices(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    jobs_path, _ = _write_one_job(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_args(tmp_path, jobs_path, "0, 1"))
    args = parse_args()
    assert args.boltz_python == Path("/existing/env/bin/python")
    assert args.gpu_devices == "0, 1"
    assert args.jobs == jobs_path


def test_main_writes_completion_marker_from_canonical_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    jobs_path, _ = _write_one_job(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_args(tmp_path, jobs_path))

    def fake_run_job(
        args: argparse.Namespace, job: Job, _gpu_queue: queue.Queue[str]
    ) -> None:
        destination = canonical_dir(args.output_root, job)
        destination.mkdir(parents=True)
        (destination / f"{job.sequence_name}_model_{job.model_index}.pdb").write_text(
            "PDB\n"
        )
        (
            destination / f"confidence_{job.sequence_name}_model_{job.model_index}.json"
        ).write_text("{}\n")

    with patch(
        "scripts.run_design_boltz_predictions.run_job", side_effect=fake_run_job
    ):
        from scripts.run_design_boltz_predictions import main

        main()
    marker = json.loads((tmp_path / "out/WORKER_COMPLETE.json").read_text())
    assert marker == {"jobs": 1, "complete": 1}


def test_main_rejects_empty_gpu_device_list(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    jobs_path, _ = _write_one_job(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_args(tmp_path, jobs_path, " , "))
    from scripts.run_design_boltz_predictions import main

    with pytest.raises(ValueError, match="At least one GPU device"):
        main()
