"""Tests for the preflight-gated frozen LigandMPNN pilot executor."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from scripts import run_nyl12_ligandmpnn_pilot as runner


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def setup_campaign(tmp_path: Path) -> tuple[Path, Path, Path, dict[str, object]]:
    engine = tmp_path / "LigandMPNN"
    engine.mkdir()
    (engine / "run.py").write_text("fixture engine\n")
    input_pdb = tmp_path / "input.pdb"
    input_pdb.write_text("fixture structure\n")
    source_file = tmp_path / "source.dat"
    source_file.write_text("source fixture\n")
    source_manifest = tmp_path / "input_sources.json"
    source_manifest.write_text(
        json.dumps({"files": {"source.dat": {"sha256": _hash(source_file)}}})
    )
    jobs = []
    for temperature, replicate, seed in (
        (0.1, "r1", 20261122),
        (0.1, "r2", 20261123),
        (0.2, "r1", 20261222),
        (0.2, "r2", 20261223),
    ):
        job_id = f"model_23_t{temperature}_{replicate}"
        output_folder = tmp_path / "ligandmpnn_outputs" / job_id
        jobs.append(
            {
                "job_id": job_id,
                "argv": [
                    "/usr/bin/python",
                    str(engine / "run.py"),
                    "--pdb_path",
                    str(input_pdb),
                    "--out_folder",
                    str(output_folder),
                    "--temperature",
                    str(temperature),
                    "--seed",
                    str(seed),
                    "--batch_size",
                    "1",
                    "--number_of_batches",
                    "6",
                ],
            }
        )
    manifest_data: dict[str, object] = {
        "expected_sequences_before_deduplication": 24,
        "design_mask": {"alpha_positions": [50], "beta_positions": [39]},
        "config": {
            "models": [23],
            "protected_alpha_positions": [113, 186, *range(222, 233)],
            "protected_beta_positions": [1, 2, 3, 40, 42],
        },
        "jobs": jobs,
        "backbones": [
            {"design_pdb": str(input_pdb), "design_pdb_sha256": _hash(input_pdb)}
        ],
    }
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(manifest_data))
    receipt = tmp_path / "preflight.json"
    receipt_data = {
        "accepted": True,
        "parent_models": [23],
        "sample_budget": 24,
        "manifest_sha256": _hash(manifest),
        "input_source_manifest": str(source_manifest),
        "input_source_manifest_sha256": _hash(source_manifest),
        "execution_sha256": {
            str(engine / "run.py"): _hash(engine / "run.py"),
            str(source_file): _hash(source_file),
            str(Path(runner.__file__).resolve()): _hash(
                Path(runner.__file__).resolve()
            ),
            str(
                Path(runner.__file__).resolve().parents[1]
                / "src/nyl12_covalent/ligandmpnn_pilot_audit.py"
            ): _hash(
                Path(runner.__file__).resolve().parents[1]
                / "src/nyl12_covalent/ligandmpnn_pilot_audit.py"
            ),
        },
        "environment": runner.REQUIRED_ENVIRONMENT,
        "host": "violet",
    }
    receipt.write_text(json.dumps(receipt_data))
    return manifest, receipt, tmp_path / "results", manifest_data


def _samples_report() -> dict[str, object]:
    samples = []
    for index in range(24):
        samples.append(
            {
                "sample_id": f"job{index // 6}:{index % 6 + 1}",
                "parent_model": f"model_{40 + index // 12}",
                "job_id": f"job{index // 6}",
                "temperature": 0.1 if index // 6 % 2 == 0 else 0.2,
                "seed": 4000 + index // 6,
                "header": f"id={index % 6 + 1}",
                "alpha_sequence": "A" * 232,
                "beta_sequence": "T" * 96,
                "mutations": [],
                "mutation_count": 0,
            }
        )
    return {
        "manifest": "fixture",
        "samples": samples,
        "duplicate_groups": [],
        "summary": {
            "job_count": 4,
            "generated_sample_count": 24,
            "native_records_excluded": 4,
            "unique_sequence_count": 1,
            "duplicate_sample_count": 23,
            "mutation_count_min": 0,
            "mutation_count_max": 0,
        },
    }


@pytest.mark.parametrize(
    "failure", ["unaccepted", "manifest_hash", "execution_hash", "source_hash"]
)
def test_preflight_failures_block_all_launches(
    tmp_path: Path,
    setup_campaign: tuple[Path, Path, Path, dict[str, object]],
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
) -> None:
    manifest, receipt, output, _ = setup_campaign
    data = json.loads(receipt.read_text())
    if failure == "unaccepted":
        data["accepted"] = False
    elif failure == "manifest_hash":
        data["manifest_sha256"] = "0" * 64
    elif failure == "execution_hash":
        engine_path = Path(next(iter(data["execution_sha256"])))
        data["execution_sha256"][str(engine_path)] = "0" * 64
    else:
        (tmp_path / "source.dat").write_text("mutated\n")
    receipt.write_text(json.dumps(data))
    called = False

    def forbidden_run(
        *args: object, **kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        nonlocal called
        called = True
        raise AssertionError("subprocess must not run before preflight passes")

    monkeypatch.setattr(runner.subprocess, "run", forbidden_run)
    with pytest.raises(ValueError):
        runner.run_campaign(manifest, receipt, output)
    assert called is False
    assert not output.exists()


def test_native_failure_is_logged_statused_and_reraised(
    setup_campaign: tuple[Path, Path, Path, dict[str, object]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, receipt, output, _ = setup_campaign

    def failing_run(
        argv: list[str], **kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        stdout = kwargs["stdout"]
        stdout.write("native failure evidence\n")
        stdout.flush()
        raise subprocess.CalledProcessError(7, argv)

    monkeypatch.setattr(runner.subprocess, "run", failing_run)
    with pytest.raises(subprocess.CalledProcessError) as error:
        runner.run_campaign(manifest, receipt, output)
    assert error.value.returncode == 7
    assert (
        output / "jobs/model_23_t0.1_r1.log"
    ).read_text() == "native failure evidence\n"
    assert json.loads((output / "status.json").read_text())["status"] == "failed"
    assert (
        json.loads((output / "jobs/model_23_t0.1_r1.status.json").read_text())[
            "returncode"
        ]
        == 7
    )


def test_successful_campaign_uses_frozen_argv_and_writes_exports(
    setup_campaign: tuple[Path, Path, Path, dict[str, object]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, receipt, output, manifest_data = setup_campaign
    calls: list[tuple[list[str], dict[str, object]]] = []

    def successful_run(
        argv: list[str], **kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        calls.append((argv, kwargs))
        output_path = Path(argv[argv.index("--out_folder") + 1])
        output_path.mkdir(parents=True)
        kwargs["stdout"].write(f"ran {argv[1]}\n")
        return subprocess.CompletedProcess(argv, 0)

    audited_jobs: list[str] = []
    monkeypatch.setattr(runner.subprocess, "run", successful_run)
    monkeypatch.setattr(
        runner,
        "audit_job",
        lambda _manifest, job: audited_jobs.append(job["job_id"]) or [{}] * 6,
    )
    report = _samples_report()
    monkeypatch.setattr(runner, "audit_campaign", lambda _path: report)

    result = runner.run_campaign(manifest, receipt, output)
    assert result == report
    assert [call[0] for call in calls] == [job["argv"] for job in manifest_data["jobs"]]
    assert audited_jobs == [job["job_id"] for job in manifest_data["jobs"]]
    assert all(call[1]["timeout"] == 1200 for call in calls)
    assert all(call[1]["env"]["CUDA_VISIBLE_DEVICES"] == "" for call in calls)
    assert len(list(csv.DictReader((output / "samples.csv").open()))) == 24
    assert len((output / "all_samples.fasta").read_text().splitlines()) == 48
    assert len((output / "unique_sequences.fasta").read_text().splitlines()) == 2
    assert json.loads((output / "status.json").read_text())["status"] == "completed"
    acceptance = json.loads((output / "scientific_acceptance.json").read_text())
    assert acceptance["copy_ties_validated"] is True
    assert acceptance["catalytic_activity_validated"] is False
