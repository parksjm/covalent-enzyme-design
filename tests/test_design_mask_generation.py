"""Manifest and execution checks for the mask comparison runner."""

from __future__ import annotations

import hashlib
import csv
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

import pytest

from src.nyl12_covalent import design_mask_generation as runner


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest(tmp_path: Path) -> tuple[Path, Path]:
    engine = tmp_path / "engine" / "run.py"
    checkpoint = tmp_path / "weights.pt"
    engine.parent.mkdir()
    engine.write_text("engine")
    checkpoint.write_bytes(b"checkpoint")
    jobs: list[dict[str, Any]] = []
    sources = {str(engine): _hash(engine), str(checkpoint): _hash(checkpoint)}
    for parent in range(2):
        pdb = tmp_path / f"parent_{parent}.pdb"
        pdb.write_text(f"pdb {parent}")
        sources[str(pdb)] = _hash(pdb)
        for temperature_index, temperature in enumerate((0.1, 0.2)):
            for arm in ("wide", "focused"):
                job_id = f"parent{parent}_t{temperature_index}_{arm}"
                output = tmp_path / "job_outputs" / job_id
                argv = [
                    sys.executable,
                    str(engine),
                    "--pdb_path",
                    str(pdb),
                    "--checkpoint_ligand_mpnn",
                    str(checkpoint),
                    "--out_folder",
                    str(output),
                    "--temperature",
                    str(temperature),
                    "--seed",
                    str(parent * 2 + temperature_index + 1),
                    "--batch_size",
                    "1",
                    "--number_of_batches",
                    "6",
                ]
                jobs.append(
                    {
                        "job_id": job_id,
                        "parent": str(parent),
                        "arm": arm,
                        "argv": argv,
                        "audit_manifest": {
                            "design_mask": {
                                "alpha_positions": [1],
                                "beta_positions": [1],
                            },
                            "config": {
                                "protected_alpha_positions": [],
                                "protected_beta_positions": [],
                            },
                        },
                    }
                )
    manifest = {
        "schema_version": 1,
        "expected_jobs": 8,
        "expected_samples": 48,
        "environment": runner.REQUIRED_ENVIRONMENT,
        "source_sha256": sources,
        "jobs": jobs,
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    return path, tmp_path / "comparison"


def _native_sequences(_: Path) -> dict[str, str]:
    return {**dict.fromkeys("ABCD", "A" * 232), **dict.fromkeys("EFGH", "A" * 96)}


def _samples(_: dict[str, Any], job: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "sample_id": f"{job['job_id']}:{index}",
            "job_id": job["job_id"],
            "parent_model": str(job["parent"]),
            "temperature": float(job["argv"][job["argv"].index("--temperature") + 1]),
            "seed": int(job["argv"][job["argv"].index("--seed") + 1]),
            "header": f"id={index}, seed=1, T=0.1",
            "alpha_sequence": "A" * 232,
            "beta_sequence": "A" * 96,
            "mutation_count": 0,
            "mutations": [],
        }
        for index in range(1, 7)
    ]


def test_run_comparison_exports_all_jobs_and_cross_arm_lineage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, output = _manifest(tmp_path)
    monkeypatch.setattr(runner, "_pdb_sequences", _native_sequences)
    monkeypatch.setattr(runner, "audit_job", _samples)
    monkeypatch.setattr(
        runner.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0, "ok", ""),
    )

    report = runner.run_comparison(manifest, output)

    assert report["summary"] == {
        "job_count": 8,
        "raw_sequence_count": 48,
        "unique_sequence_count": 1,
        "by_arm": {"wide": 24, "focused": 24},
    }
    assert len(report["duplicate_groups"][0]["sample_ids"]) == 48
    assert (output / "samples.csv").is_file()
    with (output / "samples.csv").open() as stream:
        exported = list(csv.DictReader(stream))
    assert len(exported) == 48
    assert exported[0]["job_id"] == report["samples"][0]["job_id"]
    assert exported[0]["header"] == report["samples"][0]["header"]
    fasta = (output / "all_samples.fasta").read_text().splitlines()
    assert len(fasta) == 96
    assert fasta[1].count(":") == 7
    assert len(list((output / "jobs").glob("*/process.json"))) == 8


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda value: value.update(expected_samples=47), "expected_samples"),
        (lambda value: value["jobs"].pop(), "exactly 8 jobs"),
        (
            lambda value: value["jobs"][1].update(job_id=value["jobs"][0]["job_id"]),
            "Duplicate job_id",
        ),
        (
            lambda value: value["jobs"][0]["argv"].__setitem__(
                value["jobs"][0]["argv"].index("--seed") + 1, "0"
            ),
            "seed must be positive",
        ),
    ],
)
def test_manifest_contract_fails_before_subprocess(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    change: Any,
    message: str,
) -> None:
    manifest, output = _manifest(tmp_path)
    data = json.loads(manifest.read_text())
    change(data)
    manifest.write_text(json.dumps(data))
    monkeypatch.setattr(runner, "_pdb_sequences", _native_sequences)
    monkeypatch.setattr(
        runner.subprocess, "run", lambda *args, **kwargs: pytest.fail("ran")
    )
    with pytest.raises(ValueError, match=message):
        runner.run_comparison(manifest, output)
    assert not output.exists()


def test_source_hash_and_existing_output_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, output = _manifest(tmp_path)
    checkpoint = tmp_path / "weights.pt"
    checkpoint.write_bytes(b"changed")
    monkeypatch.setattr(
        runner.subprocess, "run", lambda *args, **kwargs: pytest.fail("ran")
    )
    with pytest.raises(ValueError, match="Source SHA256 mismatch"):
        runner.run_comparison(manifest, output)
    output.mkdir()
    with pytest.raises(FileExistsError):
        runner.run_comparison(manifest, output)


def test_failed_job_is_recorded_and_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, output = _manifest(tmp_path)
    monkeypatch.setattr(runner, "_pdb_sequences", _native_sequences)
    monkeypatch.setattr(
        runner.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 7, "out", "bad"),
    )
    with pytest.raises(subprocess.CalledProcessError):
        runner.run_comparison(manifest, output)
    processes = list((output / "jobs").glob("*/process.json"))
    assert processes
    assert json.loads(processes[0].read_text())["status"] == "failed"
    assert (processes[0].parent / "stderr.log").read_text() == "bad"
