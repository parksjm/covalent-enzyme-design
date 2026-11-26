"""Prepare an immutable, bounded 30-to-60 ps direct-US continuation."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from . import neb_umbrella_campaign as campaign
from . import neb_umbrella_segments as segments
from . import direct_us_continuation_runtime as runtime_identity

TARGET_PS = 60
ORIGINAL_PS = 30
BLOCK_PS = 5
MAX_WORKERS = 40
FIRST_BLOCK = 7
LAST_BLOCK = 12
PLAN_SCHEMA = 1
FIXED_BURNIN = Path("analysis/fixed_burnin_20261126")


def _sha256(path: Path) -> str:
    """Hash a file without loading large native artifacts into memory."""
    return campaign.digest(path)


def _write_new_json(path: Path, value: dict[str, Any]) -> None:
    """Create one JSON file exclusively; never replace prior evidence."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, allow_nan=False) + "\n"
    with path.open("x", encoding="utf-8") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())


def replace_unique_seed(text: str, seed: int) -> str:
    """Replace exactly one Amber `ig` seed while preserving every other input byte."""
    if type(seed) is not int or not 1 <= seed <= 2_147_483_647:
        raise ValueError("Amber random seed must be a positive 32-bit integer")
    updated, count = re.subn(r"(?im)(\big\s*=\s*)\d+", rf"\g<1>{seed}", text)
    if count != 1:
        raise ValueError(f"Expected exactly one explicit ig seed; found {count}")
    if not re.search(r"(?i)\birest\s*=\s*1\b", updated) or not re.search(
        r"(?i)\bntx\s*=\s*5\b", updated
    ):
        raise ValueError("Continuation must retain restart coordinates and velocities")
    return updated


def _source_paths(prepared: Path) -> list[Path]:
    """Return the source files whose identities govern continuation behavior."""
    root = campaign.ROOT
    paths = [
        *campaign.source_paths(),
        root / "src/nyl12_covalent/neb_umbrella_segments.py",
        root / "src/nyl12_covalent/direct_us_coordinates.py",
        root / "src/nyl12_covalent/direct_us_continuation.py",
        root / "src/nyl12_covalent/direct_us_continuation_worker.py",
        root / "src/nyl12_covalent/direct_us_continuation_scheduler.py",
        root / "src/nyl12_covalent/direct_us_continuation_runtime.py",
        root / "src/nyl12_covalent/direct_us_continuation_io.py",
        root / "scripts/continue_direct_us.py",
        prepared,
    ]
    return list(dict.fromkeys(path.resolve(strict=True) for path in paths))


def _validate_last_stage(
    cfg: dict[str, Any],
    row: dict[str, Any],
    topology: Any,
    label: str,
) -> dict[str, Any]:
    """Validate final accepted 30 ps checkpoint against the immutable original seed."""
    folder = Path(cfg["output"]) / label
    last_label = "03_production_block_06_5ps"
    stage = folder / last_label
    process_path = stage / "process.json"
    acceptance_path = stage / "acceptance.json"
    restart = stage / "restart.ncrst"
    if not all(path.is_file() for path in (process_path, acceptance_path, restart)):
        raise FileNotFoundError(
            f"Final 30 ps checkpoint evidence is incomplete: {label}"
        )
    status_path = folder / "window_status.json"
    status = json.loads(status_path.read_text())
    if (
        status.get("status") != "completed"
        or status.get("production_blocks") != 6
        or status.get("production_ps") != ORIGINAL_PS
        or last_label not in status.get("completed", [])
    ):
        raise ValueError(f"Window is not accepted at exactly 30 production ps: {label}")
    process = json.loads(process_path.read_text())
    acceptance = json.loads(acceptance_path.read_text())
    restart_hash = _sha256(restart)
    if process.get("exit_status") != 0:
        raise ValueError(f"Final native process did not exit successfully: {label}")
    if acceptance.get("checkpoint_sha256") != restart_hash:
        raise ValueError(f"Final accepted checkpoint digest mismatch: {label}")
    if (
        acceptance.get("completed_steps") != 10_000
        or acceptance.get("finite_energy_temperature") is not True
        or acceptance.get("SCC_error_detected") is not False
        or not isinstance(acceptance.get("final_temperature_K"), (int, float))
        or not math.isfinite(acceptance["final_temperature_K"])
    ):
        raise ValueError(f"Final production acceptance is invalid: {label}")

    seed = Rst7.open(str(row["seed"]))
    baseline = np.asarray(seed.coordinates).reshape(-1, 3).copy()
    if baseline.shape != (len(topology.atoms), 3) or not np.isfinite(baseline).all():
        raise ValueError(f"Original immutable seed coordinates are invalid: {label}")
    active = np.asarray(list(AmberMask(topology, cfg["active_mask"]).Selected()))
    if len(active) != 4165 or len(np.unique(active)) != 4165:
        raise ValueError("Frozen-boundary active atom count changed from 4165")
    movable = np.zeros(len(baseline), dtype=bool)
    movable[active] = True
    checkpoint, xyz = segments._checkpoint(
        restart, baseline, movable, np.asarray(seed.box)
    )
    return {
        "stage": str(stage.resolve()),
        "process_sha256": _sha256(process_path),
        "acceptance_sha256": _sha256(acceptance_path),
        "restart": str(restart.resolve()),
        "restart_sha256": restart_hash,
        "checkpoint": checkpoint,
        "atom_count": int(len(xyz)),
    }


def prepare_continuation(
    prepared: Path, control: Path, workers: int = MAX_WORKERS
) -> Path:
    """Pin a 60 ps continuation and materialize blocks 07–12 without launching MD.

    Args:
        prepared: Immutable 63-window campaign ``prepared.json``.
        control: New control directory, normally ``continuation_60ps_20261126``.
        workers: Maximum number of dynamic stage workers, from one through 40.

    Returns:
        Path to the exclusively created continuation plan.

    Raises:
        FileExistsError: If plan, lock, stage, or continuation input already exists.
        ValueError: If source identities or accepted checkpoints fail validation.
    """
    prepared = Path(prepared).resolve(strict=True)
    control = Path(control).resolve()
    if type(workers) is not int or not 1 <= workers <= MAX_WORKERS:
        raise ValueError(f"workers must be an integer from 1 through {MAX_WORKERS}")
    cfg = campaign.frozen(prepared)
    if len(cfg["windows"]) != 63:
        raise ValueError(
            "Direct-US continuation requires the prepared 63-window campaign"
        )
    if cfg["temperature_K"] != 300 or cfg["temperature_gate_K"] != [250, 350]:
        raise ValueError(
            "Prepared temperature or gate differs from the reviewed campaign"
        )
    if (
        cfg["bias_convention"]
        != "Amberrk:centralE=k*delta²;nativeouterlineartails;nohalf"
    ):
        raise ValueError("Prepared umbrella bias convention changed")
    if control.exists():
        raise FileExistsError(
            f"Continuation control directory already exists: {control}"
        )

    output = Path(cfg["output"]).resolve(strict=True)
    if control.parent != output:
        raise ValueError(
            "Continuation control directory must be directly under campaign root"
        )
    augmented_root = output / FIXED_BURNIN / "augmented"
    if not augmented_root.is_dir():
        raise FileNotFoundError(
            f"Existing fixed-burnin augmentation root is missing: {augmented_root}"
        )

    topology = parmed.load_file(cfg["topology"])
    identities: dict[str, dict[str, Any]] = {}
    stage_records: list[dict[str, Any]] = []
    used_seeds = {
        int(stage["seed"])
        for source_window in cfg["windows"]
        for stage in source_window["stages"]
    }
    for window_index, row in enumerate(cfg["windows"]):
        label = row["label"]
        identities[label] = _validate_last_stage(cfg, row, topology, label)
        source_stage = output / label / "03_production_block_06_5ps"
        source_input = source_stage / "run.mdin"
        source_restraint = source_stage / "umbrella.RST"
        source_text = source_input.read_text()
        source_seeds = re.findall(r"(?im)\big\s*=\s*(\d+)", source_text)
        original_block = next(
            spec
            for spec in row["stages"]
            if spec["label"] == "03_production_block_06_5ps"
        )
        if len(source_seeds) != 1 or int(source_seeds[0]) != int(
            original_block["seed"]
        ):
            raise ValueError(
                f"Final source input seed differs from prepared block06: {label}"
            )
        for block in range(FIRST_BLOCK, LAST_BLOCK + 1):
            seed = int(cfg["seed_base"] + window_index * 100 + block + 2)
            if seed > 2_147_483_647 or seed in used_seeds:
                raise ValueError(
                    "Continuation random seeds are duplicate or out of range"
                )
            used_seeds.add(seed)
            stage_label = f"03_production_block_{block:02d}_5ps"
            stage = output / label / stage_label
            if stage.exists():
                raise FileExistsError(f"Continuation stage already exists: {stage}")
            if not source_restraint.is_file():
                raise FileNotFoundError(source_restraint)
            mdin = replace_unique_seed(source_text, seed)
            stage_records.append(
                {
                    "window": label,
                    "block": block,
                    "label": stage_label,
                    "seed": seed,
                    "stage": str(stage.resolve()),
                    "source_mdin_sha256": _sha256(source_input),
                    "source_restraint_sha256": _sha256(source_restraint),
                    "mdin_sha256": hashlib.sha256(mdin.encode()).hexdigest(),
                    "restraint_sha256": _sha256(source_restraint),
                }
            )

    sources = {str(path): _sha256(path) for path in _source_paths(prepared)}
    external_runtime = runtime_identity.capture_external_runtime(cfg)
    plan = {
        "schema_version": PLAN_SCHEMA,
        "campaign_root": str(output),
        "control_dir": str(control),
        "prepared": str(prepared),
        "prepared_sha256": _sha256(prepared),
        "original_production_ps_per_window": ORIGINAL_PS,
        "target_production_ps_per_window": TARGET_PS,
        "continuation_production_blocks": LAST_BLOCK - FIRST_BLOCK + 1,
        "block_ps": BLOCK_PS,
        "window_count": len(cfg["windows"]),
        "workers": workers,
        "max_workers": MAX_WORKERS,
        "temperature_K": cfg["temperature_K"],
        "temperature_gate_K": cfg["temperature_gate_K"],
        "bias_convention": cfg["bias_convention"],
        "execution_sha256": cfg["execution_sha256"],
        "external_runtime_sha256": external_runtime,
        "engine_binary_sha256": external_runtime[str(Path(cfg["sander"]).resolve())],
        "topology_sha256": external_runtime[str(Path(cfg["topology"]).resolve())],
        "qm_charge": cfg["qm_charge"],
        "qm_atoms": cfg["qm_atoms"],
        "active_mask": cfg["active_mask"],
        "scc_step_verification": "not_claimed; native failure warnings are checked, full-step SCC output is absent from source input",
        "prepared_input_identity": identities,
        "source_sha256": sources,
        "windows": [
            {
                "label": row["label"],
                "stages": [
                    {
                        "label": f"03_production_block_{block:02d}_5ps",
                        "kind": "production",
                        "block": block,
                    }
                    for block in range(FIRST_BLOCK, LAST_BLOCK + 1)
                ],
            }
            for row in cfg["windows"]
        ],
        "stages": stage_records,
        "augmentation_root": str(augmented_root.resolve()),
        "sampling_convergence": "Not established by completion",
    }

    control.mkdir(parents=True, exist_ok=False)
    for record in stage_records:
        stage = Path(record["stage"])
        stage.mkdir(parents=False, exist_ok=False)
        source_stage = output / record["window"] / "03_production_block_06_5ps"
        mdin = replace_unique_seed(
            (source_stage / "run.mdin").read_text(), record["seed"]
        )
        with (stage / "run.mdin").open("x", encoding="utf-8") as stream:
            stream.write(mdin)
        with (stage / "umbrella.RST").open("xb") as stream:
            stream.write((source_stage / "umbrella.RST").read_bytes())
        if (
            _sha256(stage / "run.mdin") != record["mdin_sha256"]
            or _sha256(stage / "umbrella.RST") != record["restraint_sha256"]
        ):
            raise ValueError(f"Materialized continuation input hash mismatch: {stage}")
    plan_path = control / "plan.json"
    _write_new_json(plan_path, plan)
    with (control / "plan.sha256").open("x", encoding="utf-8") as stream:
        stream.write(_sha256(plan_path) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    return plan_path


def verify_plan(plan_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Verify plan, runtime sources, continuation inputs, and optionally all frozen inputs."""
    plan_path = Path(plan_path).resolve(strict=True)
    expected = (plan_path.parent / "plan.sha256").read_text().strip()
    if _sha256(plan_path) != expected:
        raise ValueError("Continuation plan digest mismatch")
    plan = json.loads(plan_path.read_text())
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError("Unsupported continuation plan schema")
    prepared = Path(plan["prepared"])
    if _sha256(prepared) != plan["prepared_sha256"]:
        raise ValueError("Original prepared campaign changed")
    cfg = campaign.frozen(prepared)
    for name, digest in plan["source_sha256"].items():
        if _sha256(Path(name)) != digest:
            raise ValueError(f"Continuation source changed: {name}")
    for name, digest in cfg["execution_sha256"].items():
        if _sha256(Path(name)) != digest:
            raise ValueError(f"Prepared execution binary changed: {name}")
    runtime_identity.verify_external_runtime(plan["external_runtime_sha256"])
    if len(plan["stages"]) != 63 * 6:
        raise ValueError("Continuation stage manifest must contain 378 blocks")
    for record in plan["stages"]:
        stage = Path(record["stage"])
        if _sha256(stage / "run.mdin") != record["mdin_sha256"]:
            raise ValueError(f"Continuation run.mdin changed: {stage}")
        if _sha256(stage / "umbrella.RST") != record["restraint_sha256"]:
            raise ValueError(f"Continuation umbrella restraint changed: {stage}")
    return plan, cfg
