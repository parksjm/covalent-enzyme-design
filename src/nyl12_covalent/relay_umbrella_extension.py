"""Prepare immutable restart continuation inputs for the relay umbrella campaign."""

from __future__ import annotations

import json
from pathlib import Path
import re
import shutil
from typing import Any

from .neb_umbrella_campaign import digest


def _write_json(path: Path, value: dict[str, Any]) -> None:
    """Write one new JSON artifact without replacing an existing file."""
    with path.open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def patch_runtime_source(source: str) -> str:
    """Return a strictly patched copy of the pinned runner for 12 blocks.

    The active 6-block source snapshot is never changed. Historical stage input
    identities are checked against their original plan; extension stages use
    the new extension plan.
    """
    replacements = (
        (
            "base = prepared.parent",
            'base = Path(plan["campaign"]).resolve(strict=True)',
        ),
        (
            "if not 1 <= limit <= 6:",
            'if not 1 <= limit <= plan["production_blocks"] <= 12:',
        ),
        (
            "            input_plan = plan_path",
            '            input_plan = Path(plan["stage_plan_paths"][spec["label"]])',
        ),
        (
            'status="completed" if count == 6 else "paused_at_block_boundary",',
            'status="completed" if count == plan["production_blocks"] else "paused_at_block_boundary",',
        ),
    )
    for old, new in replacements:
        if source.count(old) != 1:
            raise ValueError(f"Pinned runner patch anchor must occur once: {old}")
        source = source.replace(old, new, 1)
    return source


def patch_guard_source(source: str) -> str:
    """Return a guardian copy that scans the plan-selected prepared config."""
    old = (
        'cfg = json.loads((args.campaign / "prepared.json").read_text())\n'
        "    plan = json.loads(args.plan.read_text())"
    )
    new = (
        "plan = json.loads(args.plan.read_text())\n"
        '    cfg = json.loads(Path(plan["prepared"]).read_text())'
    )
    if source.count(old) != 1:
        raise ValueError("Pinned guardian prepared-config patch anchor must occur once")
    return source.replace(old, new, 1)


def prepare_extension(
    original_plan_path: Path,
    output_directory: Path,
    extension_guard_source: Path,
    derived_repair_snapshot: Path,
) -> dict[str, Any]:
    """Create the isolated 30-to-60 ps plan, stage inputs and pinned source copy.

    No native process is launched. Existing campaign plans, accepted stages,
    trajectories, restarts, receipts and workers are read-only inputs.

    Args:
        original_plan_path: Validated six-block execution plan for the campaign.
        output_directory: New ``continuation_60ps_20261126`` control directory.
        extension_guard_source: Versioned guarded-launch script to pin into the
            copied runtime snapshot.
        derived_repair_snapshot: Previously validated float32 trajectory-time
            reader source, with no changes to native QM/MM calculations.

    Returns:
        A compact provenance summary for the prepared extension artifacts.
    """
    original_plan_path = original_plan_path.resolve(strict=True)
    output_directory = output_directory.resolve()
    plan_hash_path = original_plan_path.with_suffix(".sha256")
    original_plan = json.loads(original_plan_path.read_text())
    if digest(original_plan_path) != plan_hash_path.read_text().strip():
        raise ValueError("Original execution plan hash mismatch")
    if (
        original_plan.get("production_blocks") != 6
        or original_plan.get("workers") != 57
    ):
        raise ValueError("Expected the active 57-window, six-block plan")
    campaign = Path(original_plan["prepared"]).resolve(strict=True).parent
    prepared_path = campaign / "prepared.json"
    prepared = json.loads(prepared_path.read_text())
    if digest(prepared_path) != original_plan["prepared_sha256"]:
        raise ValueError("Original prepared input changed")
    for name, expected in prepared["frozen_files_sha256"].items():
        if digest(campaign / name) != expected:
            raise ValueError(f"Original frozen input changed: {name}")
    if prepared.get("planned_production_blocks_per_window") != 6:
        raise ValueError("Original prepared campaign is not the six-block baseline")
    if len(prepared.get("windows", [])) != 57:
        raise ValueError("Expected 57 frozen windows")
    if (output_directory / "plan.json").exists():
        raise FileExistsError(f"Extension plan already exists: {output_directory}")
    if (output_directory / "prepared.json").exists():
        raise FileExistsError("Prepared extension config already exists")
    for row in prepared["windows"]:
        for stage in row["stages"]:
            if stage["kind"] == "production":
                stage_dir = campaign / row["label"] / stage["label"]
                if not (stage_dir / "run.mdin").is_file():
                    raise FileNotFoundError(
                        f"Missing frozen production input: {stage_dir}"
                    )
        for block in range(9, 15):
            if (campaign / row["label"] / f"{block:02d}_production_5ps").exists():
                raise FileExistsError(
                    f"Extension stage already exists: {row['label']} block {block}"
                )

    source_paths = list(original_plan["source_sha256"])
    source_root = next(
        (
            Path(name.split("/source_snapshot/")[0]) / "source_snapshot"
            for name in source_paths
            if "/source_snapshot/" in name
        ),
        None,
    )
    if source_root is None or not source_root.is_dir():
        raise FileNotFoundError("Pinned live source snapshot unavailable for safe copy")
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    output_directory.mkdir(exist_ok=True)
    extension_source = output_directory / "source_snapshot"
    if extension_source.exists():
        for name, expected in original_plan["source_sha256"].items():
            clone = extension_source / Path(name).relative_to(source_root)
            if digest(clone) != expected:
                raise ValueError(f"Unmodified pinned source clone mismatch: {clone}")
    else:
        shutil.copytree(source_root, extension_source)
    runtime_copy = extension_source / "src/nyl12_covalent/relay_umbrella_runtime.py"
    guard_copy = extension_source / "scripts/guard_relay_umbrella.py"
    runtime_copy.write_text(patch_runtime_source(runtime_copy.read_text()))
    # Use the already tested derived-only guardian and float32 time-reader repair.
    project = Path(__file__).resolve().parents[2]
    guard_copy.write_text(
        patch_guard_source((project / "scripts/guard_relay_umbrella.py").read_text())
    )
    for name in ("relay_umbrella_diagnostics.py",):
        candidate = derived_repair_snapshot / "src/nyl12_covalent" / name
        if not candidate.is_file():
            raise FileNotFoundError(
                f"Validated derived-reader repair missing: {candidate}"
            )
        shutil.copy2(candidate, extension_source / "src/nyl12_covalent" / name)
    shutil.copy2(
        extension_guard_source,
        extension_source / "scripts/guard_relay_umbrella_extension.py",
    )

    extension_prepared = json.loads(json.dumps(prepared))
    extension_prepared["planned_production_blocks_per_window"] = 12
    extension_prepared["planned_production_ps_per_window"] = 60
    extension_prepared["extension_parent_prepared_sha256"] = digest(prepared_path)
    extension_prepared["extension_parent_plan_sha256"] = digest(original_plan_path)
    frozen = extension_prepared["frozen_files_sha256"]
    frozen["prepared.json"] = digest(prepared_path)
    extended_stage_paths: list[str] = []
    first_seed = 1 + max(
        stage["seed"] for row in prepared["windows"] for stage in row["stages"]
    )
    for window_index, row in enumerate(extension_prepared["windows"]):
        previous = next(
            stage for stage in row["stages"] if stage["label"] == "08_production_5ps"
        )
        template_text = (campaign / previous["input"]).read_text()
        restraint_text = (
            campaign / row["label"] / previous["label"] / "umbrella.RST"
        ).read_text()
        for block in range(9, 15):
            stage_label = f"{block:02d}_production_5ps"
            stage_directory = campaign / row["label"] / stage_label
            stage_directory.mkdir()
            restraint = stage_directory / "umbrella.RST"
            restraint.write_text(restraint_text)
            seed = first_seed + 6 * window_index + (block - 9)
            if seed >= 2**31:
                raise ValueError("Thermostat seed exceeds signed32bit range")
            mdin = stage_directory / "run.mdin"
            text, changes = re.subn(r"\big\s*=\s*\d+", f"ig={seed}", template_text)
            if changes != 1:
                raise ValueError("Expected one explicit production thermostat seed")
            mdin.write_text(text)
            frozen[str(restraint.relative_to(campaign))] = digest(restraint)
            frozen[str(mdin.relative_to(campaign))] = digest(mdin)
            extended_stage_paths.extend((str(restraint), str(mdin)))
            row["stages"].append(
                {
                    "input": str(mdin.relative_to(campaign)),
                    "kind": "production",
                    "label": stage_label,
                    "output_restart": str(
                        (Path(row["label"]) / stage_label / "restart.ncrst")
                    ),
                    "restraint": str(restraint.relative_to(campaign)),
                    "seed": seed,
                    "starting_restart": str(
                        Path(row["label"]) / previous["label"] / "restart.ncrst"
                    ),
                    "steps": 10_000,
                }
            )
            previous = row["stages"][-1]

    prepared_extension_path = output_directory / "prepared.json"
    _write_json(prepared_extension_path, extension_prepared)
    extension_plan = dict(original_plan)
    extension_plan.update(
        prepared=str(prepared_extension_path),
        prepared_sha256=digest(prepared_extension_path),
        production_blocks=12,
        initial_workers=57,
        campaign=str(campaign),
        session_prefix="relay_us_cont60_20261126",
        original_plan=str(original_plan_path),
        original_plan_sha256=digest(original_plan_path),
        stage_plan_paths={
            stage["label"]: str(original_plan_path)
            for row in prepared["windows"]
            for stage in row["stages"]
        },
        extension_gate_requires_root_review=True,
        extension_from_ps=30,
        extension_to_ps=60,
        extension_new_blocks=list(range(7, 13)),
    )
    extension_plan_path = output_directory / "plan.json"
    extension_plan["stage_plan_paths"].update(
        {
            f"{index:02d}_production_5ps": str(extension_plan_path)
            for index in range(9, 15)
        }
    )
    source_prefix = str(source_root.resolve())
    new_source_prefix = str(extension_source.resolve())
    new_sources = {}
    for name in original_plan["source_sha256"]:
        if not name.startswith(source_prefix + "/"):
            raise ValueError(f"Pinned source path is outside source snapshot: {name}")
        new_path = Path(name.replace(source_prefix, new_source_prefix, 1))
        new_sources[str(new_path)] = digest(new_path)
    new_sources[str(extension_source / "scripts/guard_relay_umbrella_extension.py")] = (
        digest(extension_source / "scripts/guard_relay_umbrella_extension.py")
    )
    for path in extension_source.rglob("*.py"):
        new_sources[str(path)] = digest(path)
    extension_plan["source_sha256"] = new_sources
    extension_plan["runtime_sha256"][str(original_plan_path)] = digest(
        original_plan_path
    )
    _write_json(extension_plan_path, extension_plan)
    (output_directory / "plan.sha256").write_text(digest(extension_plan_path) + "\n")
    summary = {
        "status": "prepared_waiting_for_root_review",
        "campaign": str(campaign),
        "original_plan": str(original_plan_path),
        "original_plan_sha256": digest(original_plan_path),
        "prepared_extension": str(prepared_extension_path),
        "prepared_extension_sha256": digest(prepared_extension_path),
        "extension_plan": str(extension_plan_path),
        "extension_plan_sha256": digest(extension_plan_path),
        "source_snapshot": str(extension_source),
        "source_snapshot_file_count": len(new_sources),
        "windows": 57,
        "new_stage_count": 6,
        "new_stage_paths": len(extended_stage_paths),
        "new_stage_labels": [f"{block:02d}_production_5ps" for block in range(9, 15)],
        "native_launched": False,
        "baseline_blocks_unchanged": 6,
    }
    _write_json(output_directory / "preparation_provenance.json", summary)
    return summary
