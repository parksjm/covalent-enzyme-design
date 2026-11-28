"""Prepare immutable, bounded 300 K relay replicas without launching jobs."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import re
import shutil
from typing import Any

from .neb_umbrella_campaign import digest, write_json
from .relay_umbrella_cv import RelayBiasWindow
from .relay_umbrella_prepare import native_stage_input

STAGES = (
    ("01_equilibration_1ps", "equil", 2000),
    ("02_equilibration_4ps", "equil", 8000),
    ("03_production_5ps", "production", 10000),
    ("04_production_5ps", "production", 10000),
)


def target_stage_input(
    source: str,
    kind: str,
    seed: int,
    active_mask: str,
    fresh: bool = False,
    steps: int = 10000,
) -> str:
    """Render bounded native NVT input, optionally drawing fresh 300 K velocities.

    Args:
        source: Audited source MDIN with its original QM and active masks.
        kind: Native equil or production stage.
        seed: Explicit positive Langevin random seed.
        active_mask: Unchanged source active mask.
        fresh: Initialize velocities at 300 K for the first equilibration only.
        steps: Explicit positive integration-step count at 0.5 fs per step.

    Returns:
        Native input preserving the source Hamiltonian and frozen boundary.
    """
    if kind not in {"equil", "production"} or (fresh and kind != "equil"):
        raise ValueError("Fresh velocities are allowed only for equilibration")
    text = native_stage_input(
        source, kind, seed, Path("umbrella.RST"), active_mask, steps=steps
    )
    if fresh:
        marker = "irest=1, ntx=5"
        if text.count(marker) != 1:
            raise ValueError("Expected exactly one native restart control")
        text = text.replace(marker, "irest=0, ntx=1, tempi=300.0", 1)
    return text


def verify_hashes(files: dict[str, str]) -> None:
    """Fail on a missing, relative, or changed pinned file."""
    if not files:
        raise ValueError("Explicit file hash manifest is required")
    for name, expected in files.items():
        path = Path(name)
        if not path.is_absolute() or digest(path) != expected:
            raise ValueError(f"Pinned file changed or path is relative: {name}")


def _absolute(value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise ValueError(f"Absolute path required: {value}")
    return path.resolve(strict=True)


def prepare_targets(plan_path: Path) -> Path:
    """Freeze one new campaign child from an explicit reviewed JSON plan.

    Args:
        plan_path: JSON containing parent_prepared, output, topology,
            source_mdin, sander, environment, runtime_sha256, source_manifest,
            and targets. Each target declares label, source_restart and its
            SHA256, feature_names, centers_A, k=25 and initial_ig. Optional
            source_lineage and validation slots are retained verbatim.

    Returns:
        Absolute child prepared.json. Existing destinations are never reused.
    """
    plan_path = plan_path.resolve(strict=True)
    plan = json.loads(plan_path.read_text())
    parent_path = _absolute(plan["parent_prepared"])
    parent = json.loads(parent_path.read_text())
    if digest(parent_path) != parent_path.with_suffix(".sha256").read_text().strip():
        raise ValueError("Parent prepared hash mismatch")
    topology = _absolute(plan["topology"])
    source_mdin = _absolute(plan["source_mdin"])
    sander = _absolute(plan["sander"])
    verify_hashes(plan["runtime_sha256"])
    verify_hashes(plan["source_manifest"])
    if str(sander) not in plan["runtime_sha256"]:
        raise ValueError("Runtime manifest must pin the engine")
    for name, expected in parent["frozen_files_sha256"].items():
        if digest(parent_path.parent / name) != expected:
            raise ValueError(f"Parent scientific input changed: {name}")
    if digest(topology) != digest(parent_path.parent / parent["topology"]):
        raise ValueError("Topology differs from accepted parent topology")
    output = Path(plan["output"])
    if (
        not output.is_absolute()
        or output.name != "additional_sampling_20261126"
        or output.parent.resolve() != parent_path.parent
    ):
        raise ValueError("Output must be the existing campaign's named child")
    source = source_mdin.read_text()
    active = parent["system"]["active_mask"]
    rows: list[dict[str, Any]] = []
    labels: set[str] = set()
    seeds: set[int] = set()
    for target in plan["targets"]:
        label = target["label"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", label) or label in labels:
            raise ValueError("Target labels must be unique safe directory names")
        labels.add(label)
        seed = target["initial_ig"]
        if type(seed) is not int or not 0 < seed <= 2147483647 - len(STAGES) + 1:
            raise ValueError("Seed and all stage offsets must be positive int32")
        stage_seeds = set(range(seed, seed + len(STAGES)))
        if seeds & stage_seeds:
            raise ValueError("Stage seeds must be unique within and across targets")
        seeds.update(stage_seeds)
        if target["k"] != 25:
            raise ValueError("Target force constant must be explicitly 25")
        bias = RelayBiasWindow(
            tuple(target["feature_names"]), tuple(target["centers_A"]), (25.0, 25.0)
        )
        restart = _absolute(target["source_restart"])
        if digest(restart) != target["source_restart_sha256"]:
            raise ValueError("Source restart hash mismatch")
        stages = []
        for index, (name, kind, steps) in enumerate(STAGES):
            stage_seed = seed + index
            text = target_stage_input(
                source, kind, stage_seed, active, fresh=index == 0, steps=steps
            )
            stages.append(
                dict(
                    label=name,
                    kind=kind,
                    steps=steps,
                    pmf_included=kind == "production",
                    sampling_role="production"
                    if kind == "production"
                    else "equilibration_excluded_from_pmf",
                    seed=stage_seed,
                    input=f"{label}/{name}/run.mdin",
                    text=text,
                )
            )
        rows.append(
            dict(
                label=label,
                seed=f"{label}/seed.ncrst",
                bias=asdict(bias),
                biased_cv_names=list(bias.feature_names),
                stages=stages,
                source_restart=str(restart),
                source_restart_sha256=digest(restart),
                source_lineage=target.get("source_lineage"),
                validation=target.get("validation", {}),
            )
        )
    if not rows:
        raise ValueError("At least one target is required")
    output.mkdir(exist_ok=False)
    frozen: dict[str, str] = {}
    for row in rows:
        folder = output / row["label"]
        folder.mkdir()
        seed_path = output / row["seed"]
        shutil.copyfile(row["source_restart"], seed_path)
        frozen[row["seed"]] = digest(seed_path)
        for spec in row["stages"]:
            mdin = output / spec["input"]
            mdin.parent.mkdir()
            mdin.write_text(spec.pop("text"))
            rst = mdin.with_name("umbrella.RST")
            rst.write_text(RelayBiasWindow(**row["bias"]).restraint_text())
            for path in (mdin, rst):
                frozen[str(path.relative_to(output))] = digest(path)
    # Preserve diagnostic schema and the original physical Hamiltonian identity.
    cfg = dict(parent)
    cfg.update(
        windows=rows,
        topology=str(topology),
        temperature_K=300.0,
        parent_prepared=str(parent_path),
        targets_plan=plan,
        target_plan_sha256=digest(plan_path),
        validation=plan.get("validation", {}),
        frozen_files_sha256=frozen,
        external_files_sha256={
            str(parent_path): digest(parent_path),
            str(topology): digest(topology),
            str(source_mdin): digest(source_mdin),
            **{
                str(parent_path.parent / name): value
                for name, value in parent["frozen_files_sha256"].items()
            },
            **plan["runtime_sha256"],
            **plan["source_manifest"],
        },
        temperature_gate_K=plan.get("temperature_gate_K", [150.0, 450.0]),
        production_horizon_ps=10.0,
        automatic_primary_analysis=False,
    )
    analysis = parent_path.parent / "analysis_request.json"
    shutil.copyfile(analysis, output / analysis.name)
    cfg["frozen_files_sha256"][analysis.name] = digest(output / analysis.name)
    prepared = output / "prepared.json"
    write_json(prepared, cfg)
    prepared.with_suffix(".sha256").write_text(digest(prepared) + "\n")
    return prepared


def verify_target(prepared: Path, expected_sha256: str) -> dict[str, Any]:
    """Verify the caller-pinned preparation and all frozen execution inputs."""
    if digest(prepared) != expected_sha256:
        raise ValueError("Target prepared hash mismatch")
    if prepared.with_suffix(".sha256").read_text().strip() != expected_sha256:
        raise ValueError("Target sidecar mismatch")
    cfg = json.loads(prepared.read_text())
    verify_hashes(cfg["external_files_sha256"])
    verify_hashes(
        {
            str(prepared.parent / name): value
            for name, value in cfg["frozen_files_sha256"].items()
        }
    )
    return cfg
