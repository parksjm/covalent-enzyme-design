"""Matched exact-precursor searches and stock Boltz prediction orchestration."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

from nyl12_covalent.precursor_paired_msa import (
    read_precursor_a3m,
    split_precursor_records,
    write_paired_precursor_msas,
)

ROOT = Path(
    os.environ.get("NYL12_PROJECT_ROOT", Path(__file__).resolve().parents[2])
).resolve()
SEARCH_FLAGS = [
    "--gpu",
    "1",
    "--prefilter-mode",
    "1",
    "--threads",
    "2",
    "--use-env",
    "1",
    "--use-templates",
    "0",
]
DATABASE_NAMES = ("uniref30_2302_db", "colabfold_envdb_202108_db")


def sha(path: Path) -> str:
    """Return the SHA256 of one explicit input."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Publish finite JSON atomically within the destination directory."""
    temporary = path.with_suffix(".pending")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def journal(phase: str, action: str) -> None:
    """Append an action record before dependent work or after its result."""
    if not (ROOT / "scripts/record_agent_action.py").is_file():
        raise ValueError("NYL12_PROJECT_ROOT must identify the active project journal")
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "matched_original_boltz",
            "--phase",
            phase,
            "--permission-mode",
            "escalated",
            "--action",
            action,
            "--reason",
            "PI-authorized original-backbone matched precursor controls",
            "--rollback",
            "Retain immutable inputs and receipts; new attempt for failures",
        ],
        check=True,
    )


def candidate(plan: dict[str, Any], name: str) -> dict[str, Any]:
    """Require one exact mature-chain query and the approved positive-phi target."""
    if plan["target_phi_sign"] != 1:
        raise ValueError("Explicit approved target_phi_sign=+1 required")
    matches = [c for c in plan["candidates"] if c["candidate"] == name]
    if (
        len(matches) != 1
        or not name
        or any(
            c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-"
            for c in name
        )
    ):
        raise ValueError("Unique safe candidate identifier required")
    item = matches[0]
    a, b = item["sequences"]["alpha"], item["sequences"]["beta"]
    if (
        len(a) != 232
        or len(b) != 96
        or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in a + b)
    ):
        raise ValueError("Exact canonical 232+96 mature sequence required")
    return item


def require_idle(gpu: str) -> None:
    """Fail on an occupied assigned GPU without touching other processes."""
    result = subprocess.run(
        [
            "nvidia-smi",
            "--query-compute-apps=gpu_uuid,pid",
            "--format=csv,noheader,nounits",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if any(line.split(",")[0].strip() == gpu for line in result.stdout.splitlines()):
        raise RuntimeError(f"Assigned GPU occupied: {gpu}")


def search(plan: dict[str, Any], output: Path, name: str) -> Path:
    """Run one exact query once, preserving native process and paired acceptance."""
    item = candidate(plan, name)
    query = item["sequences"]["alpha"] + item["sequences"]["beta"]
    settings = plan["search"]
    database = Path(settings["database"])
    executable = Path(settings["executable"])
    identity = {
        "query_sha256": hashlib.sha256(query.encode()).hexdigest(),
        "command_flags": SEARCH_FLAGS,
        "executable": str(executable),
        "executable_sha256": sha(executable),
        "database_files": {
            n: {
                "path": str(database / n),
                "bytes": (database / n).stat().st_size,
                "mtime_ns": (database / n).stat().st_mtime_ns,
            }
            for n in DATABASE_NAMES
        },
    }
    work = output / "searches" / name
    acceptance = work / "source_acceptance.json"
    if acceptance.exists():
        record = json.loads(acceptance.read_text())
        if record["identity"] != identity or record["status"] != "accepted":
            raise ValueError("Exact query/database/settings cache identity changed")
        source = Path(record["precursor_a3m"])
        if sha(source) != record["precursor_a3m_sha256"]:
            raise ValueError("Accepted MSA source changed")
        return source
    journal("intent", f"Run exact precursor search {name} on {item['gpu_uuid']}")
    work.mkdir(parents=True, exist_ok=False)
    fasta = work / "query.fasta"
    fasta.write_text(f">{name}\n{query}\n")
    require_idle(item["gpu_uuid"])
    argv = [
        str(executable),
        *SEARCH_FLAGS,
        str(fasta),
        str(database),
        str(work / "local_attempt01"),
    ]
    overrides = {
        "PATH": str(executable.parent) + ":" + os.environ["PATH"],
        "CUDA_VISIBLE_DEVICES": item["gpu_uuid"],
        "OMP_NUM_THREADS": "2",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "2",
    }
    process = {
        "status": "running",
        "command": argv,
        "environment_overrides": overrides,
        "identity": identity,
        "source_sha256": sha(Path(__file__)),
        "query_fasta_sha256": sha(fasta),
    }
    write_json(work / "process.json", process)
    with (work / "search.log").open("w") as handle:
        result = subprocess.run(
            argv,
            env={**os.environ, **overrides},
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    process.update(
        status="process_completed" if result.returncode == 0 else "failed",
        exit_status=result.returncode,
    )
    write_json(work / "process.json", process)
    journal("outcome", f"Exact precursor search {name} exited {result.returncode}")
    if result.returncode:
        raise RuntimeError(f"Search failed for {name}; retain attempt evidence")
    source = work / "local_attempt01" / f"{name}.a3m"
    paired = split_precursor_records(
        read_precursor_a3m(source), alpha_query=query[:232], beta_query=query[232:]
    )
    write_json(
        acceptance,
        {
            "status": "accepted",
            "identity": identity,
            "precursor_a3m": str(source.resolve()),
            "precursor_a3m_sha256": sha(source),
            "paired_homolog_rows_after_native_dedup": paired[
                "paired_homolog_row_count"
            ],
            "source_row_count": paired["source_row_count"],
            "cleavage_after_precursor_residue1": 232,
        },
    )
    return source


def prepare(plan: dict[str, Any], output: Path, name: str) -> list[Path]:
    """Stage both seeds from a validated exact-query precursor cache."""
    item = candidate(plan, name)
    if "precursor_a3m" in item:
        source = Path(item["precursor_a3m"])
        accepted = json.loads(Path(item["source_acceptance"]).read_text())
        if (
            accepted["status"] != "accepted"
            or sha(source) != accepted["precursor_a3m_sha256"]
        ):
            raise ValueError("Reused MSA source requires accepted unchanged provenance")
    else:
        acceptance = json.loads(
            (output / "searches" / name / "source_acceptance.json").read_text()
        )
        source = Path(acceptance["precursor_a3m"])
        if (
            acceptance["status"] != "accepted"
            or sha(source) != acceptance["precursor_a3m_sha256"]
        ):
            raise ValueError("Accepted unchanged exact-query search required")
    graph_path = Path(plan["ligand_graph_manifest"])
    graph = json.loads(graph_path.read_text())
    if (
        graph["heavy_atoms"] != 33
        or graph["formal_charge"] != 0
        or graph["covalent_ES_bonds"]
    ):
        raise ValueError("Whole zwitterionic noncovalent PA6L4 graph required")
    Path(plan["approved_reference_manifest"]).read_bytes()
    msas = output / "msas" / name
    write_paired_precursor_msas(
        source,
        msas,
        alpha_query=item["sequences"]["alpha"],
        beta_query=item["sequences"]["beta"],
    )
    if plan["seeds"] != [20261126, 20261127]:
        raise ValueError("Fixed matched seed pair required")
    cases = []
    for seed in plan["seeds"]:
        case_id = f"{name}_seed_{seed}"
        work = output / "cases" / case_id
        work.mkdir(parents=True, exist_ok=False)
        path = work / f"{case_id}.yaml"
        document = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": list("ABCD"),
                        "sequence": item["sequences"]["alpha"],
                        "msa": str(msas / "alpha_paired.csv"),
                    }
                },
                {
                    "protein": {
                        "id": list("EFGH"),
                        "sequence": item["sequences"]["beta"],
                        "msa": str(msas / "beta_paired.csv"),
                    }
                },
                {"ligand": {"id": list("IJKL"), "smiles": graph["smiles"]}},
            ],
        }
        write_json(path, document)
        write_json(
            work / "case.json",
            {
                **item,
                "case": case_id,
                "seed": seed,
                "input": str(path),
                "remote_input": str(path),
                "remote_output": str(work / "runs"),
                "remote_job": str(work / "job"),
                "msas": str(msas),
                "source_sha256": {
                    str(source): sha(source),
                    str(graph_path): sha(graph_path),
                    plan["approved_reference_manifest"]: sha(
                        Path(plan["approved_reference_manifest"])
                    ),
                },
            },
        )
        cases.append(work)
    return cases
