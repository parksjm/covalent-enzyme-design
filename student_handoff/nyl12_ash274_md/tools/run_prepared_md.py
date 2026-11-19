#!/usr/bin/env python3
"""Plan/launch the unchanged10+90ns protocol from a validated minimized dataset."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from prepared_md_support import (  # noqa:E402
    digest,
    engines,
    load_dataset,
    replicas,
    require_gpu,
)

LOGGER = logging.getLogger(__name__)


def context(repo: Path) -> None:
    """Require the supplied repository and expose its unchanged protocol/launcher."""
    for name in (
        "scripts/run_model30_md_replica.py",
        "scripts/continue_model30_md.py",
        "src/nyl12_covalent/pilot_md_protocol.py",
        "src/nyl12_covalent/jobs.py",
    ):
        if not (repo / name).is_file():
            raise FileNotFoundError(repo / name)
    sys.path.insert(0, str(repo))
    sys.path.insert(0, str(repo / "src"))


def write_json(path: Path, value: Any) -> None:
    """Write a new immutable plan/process record, refusing an existing artifact."""
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def plan(config_path: Path, repo: Path) -> Path:
    """Freeze a portable100ns plan and exact original protocol inputs without MD.

    Args:
        config_path: Explicit dataset, engine, output and independent replicas.
        repo: Repository checkout providing the authoritative protocol modules.

    Returns:
        New output directory with plan, environment and stage-input snapshots.
    """
    context(repo)
    from nyl12_covalent.pilot_md_protocol import mdin_files, protocol_manifest
    from scripts.continue_model30_md import continuation_input

    config = json.loads(config_path.read_text())
    rows = replicas(config)
    manifest = (config_path.parent / config["dataset_manifest"]).resolve(strict=True)
    dataset, topology, restart, hashes = load_dataset(manifest, config["dataset"])
    home, engine = engines(config, config_path.parent)
    output = (config_path.parent / config["output"]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    sources = [
        Path(__file__).resolve(),
        HERE / "prepared_md_support.py",
        config_path,
        repo / "scripts/run_model30_md_replica.py",
        repo / "scripts/continue_model30_md.py",
        repo / "src/nyl12_covalent/pilot_md_protocol.py",
        repo / "src/nyl12_covalent/jobs.py",
        engine,
    ]
    hashes.update({str(p): digest(p) for p in sources})
    snapshot = output / "source"
    snapshot.mkdir()
    for index, source in enumerate(sources):
        if source != engine:
            shutil.copy2(source, snapshot / f"{index:02d}_{source.name}")
    for row in rows:
        folder = output / row["name"]
        folder.mkdir()
        inputs = folder / "planned_inputs"
        inputs.mkdir()
        text = mdin_files(row["seed"])
        text["06_production_90ns"] = continuation_input(
            text["05_production_10ns"], row["seed"] + 50000
        )
        for name, value in text.items():
            (inputs / f"{name}.in").write_text(value)
        record = protocol_manifest(row["seed"])
        record.update(
            production_total_ns=100,
            continuation_seed=row["seed"] + 50000,
            supplied_minimized_restart=str(restart),
            minimization="Accepted upstream unrestrained MMminimized coordinates; not rerun",
            preparation_restraints="Heating5,density2,weak0.5kcal/mol/A²; equilibration/productionntr0",
            planned_input_sha256={str(p): digest(p) for p in inputs.iterdir()},
        )
        write_json(folder / "protocol_plan.json", record)
    write_json(
        output / "environment.json",
        {
            "python": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "numpy": importlib.metadata.version("numpy"),
            "parmed": importlib.metadata.version("ParmEd"),
            "AMBERHOME": str(home),
            "cuda_engine": str(engine),
            "GPU_check": "Deferred until explicit launch/execute; plan launches nothing",
        },
    )
    write_json(
        output / "plan.json",
        {
            "config": config,
            "repo": str(repo),
            "dataset": dataset,
            "topology": str(topology),
            "restart": str(restart),
            "source_sha256": hashes,
            "amber_home": str(home),
            "cuda_engine": str(engine),
            "production_ns": 100,
            "production_restraints": "ntr=0 nmropt=0",
            "acceptance_scope": "Process/temperature gates are not sampling convergence or catalytic acceptance",
        },
    )
    LOGGER.info("Prepared %s; no MD launched", output)
    return output


def frozen(
    config_path: Path, repo: Path
) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    """Require a current immutable plan before any CUDA check or launch."""
    context(repo)
    config = json.loads(config_path.read_text())
    output = (config_path.parent / config["output"]).resolve(strict=True)
    record = json.loads((output / "plan.json").read_text())
    if record["config"] != config or record["repo"] != str(repo):
        raise ValueError(
            "Config/repository differs from frozen plan; prepare a new output"
        )
    home, engine = engines(config, config_path.parent)
    if str(home) != record["amber_home"] or str(engine) != record["cuda_engine"]:
        raise ValueError("Resolved Amber/CUDA environment differs from frozen plan")
    for path, expected in record["source_sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Frozen source/input changed: {path}")
    return config, output, record


def execute(config_path: Path, repo: Path, name: str) -> None:
    """Execute only the named prepared replica using existing10ns/90ns entrypoints."""
    config, output, record = frozen(config_path, repo)
    row = next(r for r in config["replicas"] if r["name"] == name)
    gpu = require_gpu(str(row["gpu"]))
    if os.environ.get("CUDA_VISIBLE_DEVICES") != str(row["gpu"]):
        raise ValueError("CUDA_VISIBLE_DEVICES differs from planned explicit GPU")
    folder = output / name
    write_json(folder / "gpu_check.json", gpu)
    engine = engines(config, config_path.parent)[1]
    command = [
        sys.executable,
        str(repo / "scripts/run_model30_md_replica.py"),
        "--topology",
        record["topology"],
        "--coordinates",
        record["restart"],
        "--engine",
        str(engine),
        "--output",
        str(folder / "dynamics"),
        "--seed",
        str(row["seed"]),
        "--label",
        f"{config['dataset']}_{name}",
    ]
    continuation = [
        sys.executable,
        str(repo / "scripts/continue_model30_md.py"),
        "--replica-dir",
        str(folder / "dynamics"),
        "--source-job-dir",
        str(folder / "job"),
        "--engine",
        str(engine),
        "--seed",
        str(row["seed"] + 50000),
    ]
    for label, argv in (("pilot_10ns", command), ("continuation_90ns", continuation)):
        write_json(
            folder / f"{label}.command.json",
            {"command": argv, "started_utc": datetime.now(timezone.utc).isoformat()},
        )
        subprocess.run(argv, check=True, cwd=repo)
    write_json(
        folder / "100ns_process_complete.json",
        {"completed": True, "sampling_convergence": "Not assessed"},
    )


def launch(config_path: Path, repo: Path) -> None:
    """Launch explicit independent GPU replicas through the unchanged tmux job logger."""
    config, output, record = frozen(config_path, repo)
    from nyl12_covalent.jobs import launch_job

    home, _ = engines(config, config_path.parent)
    for row in config["replicas"]:
        require_gpu(str(row["gpu"]))
    for row in config["replicas"]:
        launch_job(
            f"{config['session_prefix']}_{row['name']}",
            output / row["name"] / "job",
            repo,
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "execute",
                "--repo",
                str(repo),
                "--config",
                str(config_path),
                "--replica",
                row["name"],
            ],
            env_overrides={
                "AMBERHOME": str(home),
                "LD_LIBRARY_PATH": str(home / "lib"),
                "CUDA_VISIBLE_DEVICES": str(row["gpu"]),
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
            },
            input_paths=[Path(p) for p in record["source_sha256"]],
        )


def main() -> None:
    """Require an explicit plan/launch/execute action; default never launches MD."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "launch", "execute"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--replica")
    args = parser.parse_args()
    repo = args.repo.resolve(strict=True)
    config = args.config.resolve(strict=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.action == "plan":
        plan(config, repo)
    elif args.action == "launch":
        launch(config, repo)
    else:
        if not args.replica:
            parser.error("execute requires --replica")
        execute(config, repo, args.replica)


if __name__ == "__main__":
    main()
