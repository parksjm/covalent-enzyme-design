"""Prepare a source-pinned continuation of the accepted 57-window 60-ps prefix."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Any


def digest(path: Path) -> str:
    """Return the exact file SHA256 identity."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def patch_block_limit(source: str, target_blocks: int) -> str:
    """Change only the copied accepted-prefix runner's maximum block count.

    Args:
        source: Validated 60-ps copied runtime source.
        target_blocks: Authorized total production blocks, exactly20.

    Returns:
        Source differing only in the maximum block-count guard.
    """
    old = 'if not 1 <= limit <= plan["production_blocks"] <= 12:'
    if target_blocks != 20 or source.count(old) != 1:
        raise ValueError("Expected exactly one validated12-block guard and target20")
    return source.replace(
        old, 'if not 1 <= limit <= plan["production_blocks"] <= 20:', 1
    )


def prepare_continuation(
    original_plan_path: Path,
    output_directory: Path,
    *,
    target_blocks: int = 20,
    workers: int = 57,
) -> dict[str, Any]:
    """Append eight unchanged-protocol blocks using the existing native queue.

    Args:
        original_plan_path: Accepted host-specific60-ps execution plan.
        output_directory: New control directory inside the same native campaign.
        target_blocks: Authorized total20 production blocks (100ps).
        workers: Explicit one-thread workers, at most57 under PI authorization.

    Returns:
        Immutable preparation/source/checkpoint receipt; launches no processes.
    """
    original_plan_path = original_plan_path.resolve(strict=True)
    old = json.loads(original_plan_path.read_text())
    if (
        digest(original_plan_path)
        != original_plan_path.with_suffix(".sha256").read_text().strip()
    ):
        raise ValueError("Parent execution plan changed")
    if old["production_blocks"] != 12 or target_blocks != 20 or not 1 <= workers <= 57:
        raise ValueError("Require accepted60ps to authorized100ps, 1..57 workers")
    base = Path(old["campaign"]).resolve(strict=True)
    output_directory = output_directory.resolve()
    if output_directory.parent != base or output_directory.exists():
        raise ValueError("Require a new continuation directory in the same campaign")
    prepared = Path(old["prepared"])
    cfg = json.loads(prepared.read_text())
    if digest(prepared) != old["prepared_sha256"] or len(cfg["windows"]) != 57:
        raise ValueError("Parent prepared identity/window inventory changed")
    for relative, value in cfg["frozen_files_sha256"].items():
        if digest(base / relative) != value:
            raise ValueError("Frozen scientific input changed: " + relative)
    for name, value in old["source_sha256"].items():
        if digest(Path(name)) != value:
            raise ValueError("Parent source changed: " + name)
    if digest(Path(old["sander"])) != old["runtime_sha256"][old["sander"]]:
        raise ValueError("Host-specific validated Sander changed")
    if cfg.get("planned_production_blocks_per_window") != 12:
        raise ValueError("Parent is not the complete12-block preparation")
    for row in cfg["windows"]:
        labels = [s["label"] for s in row["stages"] if s["kind"] == "production"]
        if labels != [f"{i:02d}_production_5ps" for i in range(3, 15)]:
            raise ValueError(
                "Parent production stages are not the complete accepted prefix"
            )
    checkpoints = {}
    for row in cfg["windows"]:
        folder = base / row["label"]
        receipt = json.loads((folder / "segment_receipt_limit_12.json").read_text())
        last = folder / "14_production_5ps"
        accept = json.loads((last / "acceptance.json").read_text())
        process = json.loads((last / "process.json").read_text())
        checkpoint = last / "restart.ncrst"
        if (
            receipt["status"] != "accepted"
            or receipt["production_blocks"] != 12
            or process["exit_status"] != 0
            or accept.get("accepted") is not True
        ):
            raise ValueError("Window lacks accepted60ps prefix: " + row["label"])
        if (
            digest(checkpoint) != receipt["checkpoint_sha256"]
            or digest(checkpoint) != accept["checkpoint_sha256"]
        ):
            raise ValueError("Accepted endpoint checkpoint changed")
        checkpoints[row["label"]] = {
            "path": str(checkpoint),
            "sha256": digest(checkpoint),
        }
        text = (last / "run.mdin").read_text()
        for key, value in [("irest", "1"), ("ntx", "5"), ("nstlim", "10000")]:
            if not re.search(r"\b" + key + r"\s*=\s*" + value + r"\b", text, re.I):
                raise ValueError("Expected unchanged restart production controls")
        if len(re.findall(r"\big\s*=\s*\d+", text)) != 1:
            raise ValueError("Expected one explicit production thermostat seed")
        for number in range(15, 23):
            if (folder / f"{number:02d}_production_5ps").exists():
                raise FileExistsError("Future stage already exists")
    old_source = original_plan_path.parent / "source_snapshot"
    if not old_source.is_dir():
        raise FileNotFoundError(old_source)
    output_directory.mkdir()
    source = output_directory / "source_snapshot"
    shutil.copytree(old_source, source)
    runtime = source / "src/nyl12_covalent/relay_umbrella_runtime.py"
    runtime.write_text(patch_block_limit(runtime.read_text(), target_blocks))
    newcfg = copy.deepcopy(cfg)
    newcfg.update(
        planned_production_blocks_per_window=20,
        planned_production_ps_per_window=100,
        extension_parent_prepared_sha256=digest(prepared),
        extension_parent_plan_sha256=digest(original_plan_path),
    )
    frozen = newcfg["frozen_files_sha256"]
    frozen[str(prepared.relative_to(base))] = digest(prepared)
    first_seed = 1 + max(s["seed"] for row in cfg["windows"] for s in row["stages"])
    for index, row in enumerate(newcfg["windows"]):
        previous = row["stages"][-1]
        if previous["label"] != "14_production_5ps":
            raise ValueError("Parent does not end at60ps")
        text = (base / previous["input"]).read_text()
        bias = (base / row["label"] / previous["label"] / "umbrella.RST").read_bytes()
        for offset, number in enumerate(range(15, 23)):
            label = f"{number:02d}_production_5ps"
            directory = base / row["label"] / label
            directory.mkdir()
            seed = first_seed + index * 8 + offset
            if seed >= 2**31:
                raise ValueError("Thermostat seed exceeds signed32bit range")
            mdin = directory / "run.mdin"
            mdin.write_text(re.sub(r"\big\s*=\s*\d+", f"ig={seed}", text))
            restraint = directory / "umbrella.RST"
            restraint.write_bytes(bias)
            for p in [mdin, restraint]:
                frozen[str(p.relative_to(base))] = digest(p)
            stage = dict(previous)
            stage.update(
                label=label,
                input=str(mdin.relative_to(base)),
                output_restart=str((directory / "restart.ncrst").relative_to(base)),
                restraint=str(restraint.relative_to(base)),
                seed=seed,
                starting_restart=str(
                    Path(row["label"]) / previous["label"] / "restart.ncrst"
                ),
                steps=10000,
            )
            row["stages"].append(stage)
            previous = stage
    newprepared = output_directory / "prepared.json"
    newprepared.write_text(json.dumps(newcfg, indent=2) + "\n")
    plan = copy.deepcopy(old)
    plan.update(
        prepared=str(newprepared),
        prepared_sha256=digest(newprepared),
        production_blocks=20,
        workers=workers,
        initial_workers=workers,
        parent_worker_count=old["workers"],
        session_prefix="relay_us_cont100_20261127",
        original_plan=str(original_plan_path),
        original_plan_sha256=digest(original_plan_path),
        extension_from_ps=60,
        extension_to_ps=100,
        extension_new_blocks=list(range(13, 21)),
        root_review_authorized=True,
    )
    planpath = output_directory / "plan.json"
    plan["stage_plan_paths"].update(
        {f"{i:02d}_production_5ps": str(planpath) for i in range(15, 23)}
    )
    plan["runtime_sha256"][str(original_plan_path)] = digest(original_plan_path)
    plan["source_sha256"] = {str(p): digest(p) for p in source.rglob("*.py")}
    for key in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        plan["environment"][key] = "1"
    plan["environment"]["CUDA_VISIBLE_DEVICES"] = ""
    planpath.write_text(json.dumps(plan, indent=2) + "\n")
    planpath.with_suffix(".sha256").write_text(digest(planpath) + "\n")
    result = {
        "status": "prepared_native_not_launched",
        "parent_plan": str(original_plan_path),
        "parent_plan_sha256": digest(original_plan_path),
        "plan": str(planpath),
        "plan_sha256": digest(planpath),
        "windows": 57,
        "workers": workers,
        "from_ps": 60,
        "to_ps": 100,
        "new_production_labels": [f"{i:02d}_production_5ps" for i in range(15, 23)],
        "accepted60ps_checkpoints": checkpoints,
        "scientific_changes": "Onlynew5psblocks and disjointig seeds underexistingblockprotocol; sameirest1ntx5, velocities/time/box retained; noextraequilibration orHamiltonian/SCF/bias changes",
        "runtime_copy_change": "Only12to20maximumblockguard; existingaccepted-prefixlogicunchanged",
    }
    (output_directory / "preparation_provenance.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result
