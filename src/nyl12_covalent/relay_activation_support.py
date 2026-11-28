"""Isolated preparation for the approved33-state physical activation support grid."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict
import json
from pathlib import Path
import re
import shutil
from typing import Any

from .neb_umbrella_campaign import digest, write_json
from .relay_targeted import _absolute, target_stage_input, verify_hashes
from .relay_umbrella_cv import RelayBiasWindow

STAGES = (
    ("01_equilibration_1ps", "equil", 2000),
    ("02_equilibration_4ps", "equil", 8000),
    *((f"{i:02d}_production_5ps", "production", 10000) for i in range(3, 9)),
)


def validate_support_grid(targets: list[dict[str, Any]]) -> None:
    """Require exactly the approved33physicaltwo-CVcenters and k25."""
    expected = {
        (x, round(-1 + 0.2 * i, 8)) for x in (-1.25, -1.5, -1.75) for i in range(11)
    }
    actual = []
    for target in targets:
        if target["feature_names"] != ["q_acyl", "q_act"] or target["k"] != 25:
            raise ValueError("Supportrequiresq_acyl/q_actandk25")
        actual.append(tuple(round(float(x), 8) for x in target["centers_A"]))
    if len(actual) != 33 or set(actual) != expected:
        raise ValueError("Supportmustcontainexactlytheapproved33centers")


def prepare_support_targets(plan_path: Path) -> Path:
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
    return prepare_reviewed_targets(
        plan_path, "activation_support_20261127", validate_support_grid
    )


def prepare_reviewed_targets(
    plan_path: Path, output_name: str, validator: Callable[[list[dict[str, Any]]], None]
) -> Path:
    """Prepare a named child after its explicit scientific plan validator passes."""
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
        or output.name != output_name
        or output.parent.resolve() != parent_path.parent
    ):
        raise ValueError("Output must be the existing campaign's named child")
    validator(plan["targets"])
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
        raise ValueError("All33supporttargetsarerequired")
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
        production_horizon_ps=30.0,
        automatic_primary_analysis=False,
    )
    analysis = parent_path.parent / "analysis_request.json"
    shutil.copyfile(analysis, output / analysis.name)
    cfg["frozen_files_sha256"][analysis.name] = digest(output / analysis.name)
    prepared = output / "prepared.json"
    write_json(prepared, cfg)
    prepared.with_suffix(".sha256").write_text(digest(prepared) + "\n")
    return prepared
