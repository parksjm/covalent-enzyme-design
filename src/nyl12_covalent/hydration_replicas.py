"""Preparation and gates for independent short hydration branches."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any

from .hydration_continuation_101ps import validate_accepted_segment
from .deacylation_hydration import build_input
from .hydration_continuation_sources_101ps import replace_once

ROOT = Path(__file__).resolve().parents[2]
BASE_RELATIVE = Path(
    "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
)
SPECS = (
    {
        "name": "replica_02_from_26ps_extend04",
        "checkpoint": 4,
        "first": 5,
        "last": 9,
        "seeds": (2026112721, 2026112722, 2026112723, 2026112724, 2026112725),
    },
    {
        "name": "replica_03_from_41ps_extend07",
        "checkpoint": 7,
        "first": 8,
        "last": 12,
        "seeds": (2026112731, 2026112732, 2026112733, 2026112734, 2026112735),
    },
)


def digest(path: Path) -> str:
    """Hash one input or receipt without modifying it."""
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, value: dict[str, Any], exclusive: bool = False) -> None:
    """Write finite JSON, optionally refusing any pre-existing artifact."""
    mode = "x" if exclusive else "w"
    with path.open(mode) as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def verify_hashes(hashes: dict[str, str]) -> None:
    """Fail closed when a pinned branch or historical input changes."""
    if not hashes:
        raise ValueError("Empty replica source/input manifest")
    for filename, expected in hashes.items():
        if digest(Path(filename)) != expected:
            raise ValueError(f"Replica input identity changed: {filename}")


def _copy_snapshot(source: Path, destination: Path) -> None:
    """Copy one immutable source snapshot, excluding Python caches."""
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("__pycache__"))


def _patch_native(
    branch: Path, spec: dict[str, Any], source_base: Path
) -> dict[str, str]:
    """Copy the accepted native core and specialize its five segment constants."""
    source = source_base / "continuation_101ps_source"
    destination = branch / "native_source"
    _copy_snapshot(source, destination)
    extension = destination / "nyl12_covalent/deacylation_hydration_extension.py"
    text = extension.read_text()
    old_names = (
        'SEGMENT_NAMES = tuple(f"extend_{index:02d}" for index in range(10, 20))'
    )
    names = f'SEGMENT_NAMES = tuple(f"extend_{{index:02d}}" for index in range({spec["first"]}, {spec["last"] + 1}))'
    text = replace_once(text, old_names, names, f"{spec['name']} segment range")
    old_seeds = "SEEDS = (2026112610, 2026112611, 2026112612, 2026112613, 2026112614, 2026112615, 2026112616, 2026112617, 2026112618, 2026112619)"
    seed_source = "SEEDS = (" + ", ".join(map(str, spec["seeds"])) + ")"
    text = replace_once(text, old_seeds, seed_source, f"{spec['name']} Langevin seeds")
    text = text.replace("continuation_101ps_config.json", "runtime_config.json")
    text = text.replace("continuation_101ps_source", "native_source")
    extension.write_text(text)
    # The frozen runner and gate modules also resolve their config/source paths.
    for module in sorted((destination / "nyl12_covalent").glob("*.py")):
        if module == extension:
            continue
        module_text = module.read_text()
        module_text = module_text.replace(
            "continuation_101ps_config.json", "runtime_config.json"
        )
        module_text = module_text.replace("continuation_101ps_source", "native_source")
        module.write_text(module_text)
    return {
        str(path.relative_to(destination)): digest(path)
        for path in sorted(destination.rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts
    }


def _patch_postprocess(
    branch: Path, spec: dict[str, Any], source_base: Path
) -> dict[str, str]:
    """Copy reviewed analysis/VMD code and cap its receipt scan at this fork end."""
    destination = branch / "postprocess_source"
    _copy_snapshot(source_base / "continuation_101ps_postprocess_source", destination)
    adapter = destination / "src/nyl12_covalent/hydration_receipt_adapter.py"
    text = replace_once(
        adapter.read_text(),
        "range(1, 20)",
        f"range(1, {spec['last'] + 1})",
        "fork receipt cap",
    )
    text = replace_once(
        text,
        "extend_01 through extend_19",
        f"extend_01 through extend_{spec['last']:02d}",
        "fork receipt range",
    )
    adapter.write_text(text)
    retained = branch / "native_source/postprocess"
    for folder in (destination / "scripts", destination / "src/nyl12_covalent"):
        for source in folder.glob("*.py"):
            shutil.copy2(source, retained / source.name)
    return {
        str(path.resolve()): digest(path)
        for folder in (destination / "scripts", destination / "src/nyl12_covalent")
        for path in sorted(folder.glob("*.py"))
    }


def _link_shared_history(branch: Path, original: Path, checkpoint: int) -> None:
    """Make workflow-owned read-only aliases to original accepted history."""
    for name in (
        "config.json",
        "hydration_endpoint_status_20261124.json",
        "extension_20ps_config.json",
        "followon_runtime_config.json",
        "inputs",
        "heat",
        "equil",
    ):
        (branch / name).symlink_to(
            (original / name).resolve(), target_is_directory=(original / name).is_dir()
        )
    for index in range(1, checkpoint + 1):
        name = f"extend_{index:02d}"
        (branch / name).symlink_to(
            (original / name).resolve(), target_is_directory=True
        )


def validate_checkpoint(original: Path, index: int) -> Path:
    """Verify every accepted block through a fork checkpoint and return its restart."""
    if index not in {4, 7}:
        raise ValueError(f"Unsupported accepted fork checkpoint: {index}")
    for prior in range(1, index + 1):
        validate_accepted_segment(original, prior)
    block = original / f"extend_{index:02d}"
    status = json.loads((block / "status.json").read_text())
    restart = block / "restart.ncrst"
    if status.get("restart_sha256") != digest(restart):
        raise ValueError(
            f"Fork restart does not match accepted receipt: extend_{index:02d}"
        )
    return restart


def prepare_branch(original: Path, spec: dict[str, Any], python: str) -> dict[str, Any]:
    """Prepare one isolated fork with independent outputs and shared history links."""
    original = original.resolve(strict=True)
    checkpoint = spec["checkpoint"]
    fork_restart = validate_checkpoint(original, checkpoint)
    branch = original / "replicas" / spec["name"]
    branch.mkdir(parents=True, exist_ok=False)
    _link_shared_history(branch, original, checkpoint)
    native_hashes = _patch_native(branch, spec, original)
    post_hashes = _patch_postprocess(branch, spec, original)
    native_hashes = {
        str(path.relative_to(branch / "native_source")): digest(path)
        for path in sorted((branch / "native_source").rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts
    }
    checkpoint_dir = branch / f"extend_{checkpoint:02d}"
    names = [f"extend_{i:02d}" for i in range(spec["first"], spec["last"] + 1)]
    seeds = dict(zip(names, spec["seeds"], strict=True))
    template = json.loads((original / "continuation_101ps_template.json").read_text())
    template.update(
        scope=f"25 ps independent hydration fork from accepted extend_{checkpoint:02d}",
        endpoint_restart_path=str((checkpoint_dir / "restart.ncrst").resolve()),
        endpoint_restart_sha256=digest(fork_restart),
        segment_names=names,
        seeds=seeds,
        source_snapshot_sha256=native_hashes,
        postprocess_source_sha256=post_hashes,
        postprocess_argv=[
            python,
            str(branch / "postprocess_source/scripts/refresh_deacylation_hydration.py"),
            "--base",
            str(branch),
            "--project-root",
            str(branch / "postprocess_source"),
        ],
        postprocess_environment={
            "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        },
        segment_input_sha256={},
    )
    source = (branch / "inputs/source.mdin").read_text()
    for offset, name in enumerate(names):
        block = branch / name
        block.mkdir(exist_ok=False)
        inputs = block / "inputs"
        inputs.mkdir()
        (inputs / "source.mdin").symlink_to((branch / "inputs/source.mdin").resolve())
        (inputs / "system.prmtop").symlink_to(
            (branch / "inputs/system.prmtop").resolve()
        )
        mdin = build_input(source, template["active_mask"], "equil", seeds[name])
        (inputs / "segment.mdin").write_text(
            mdin.replace(
                "Acyl enzyme hydration only: equil",
                f"Acyl enzyme replica {spec['name']}: {name}",
                1,
            )
        )
        predecessor_index = checkpoint + offset
        if offset == 0:
            (inputs / "restart.ncrst").symlink_to(
                (branch / f"extend_{checkpoint:02d}/restart.ncrst").resolve()
            )
        write_json(
            block / "status.json",
            {
                "status": "prepared",
                "segment": spec["first"] + offset,
                "duration_ps": 5.0,
                "predecessor": f"extend_{predecessor_index:02d}",
            },
            exclusive=True,
        )
        template["segment_input_sha256"][name] = {
            path.name: digest(path) for path in sorted(inputs.iterdir())
        }
    runtime = branch / "runtime_config.json"
    write_json(runtime, template, exclusive=True)
    pinned_paths = [
        runtime,
        original / "config.json",
        original / "hydration_endpoint_status_20261124.json",
        original / "continuation_101ps_template.json",
        original / "inputs/source.mdin",
        original / "inputs/system.prmtop",
        fork_restart,
    ]
    for index in range(1, checkpoint + 1):
        block = original / f"extend_{index:02d}"
        pinned_paths.extend(
            block / part
            for part in (
                "status.json",
                "process.json",
                "scientific_acceptance.json",
                "restart.ncrst",
                "inputs/restart.ncrst",
                "inputs/source.mdin",
                "inputs/system.prmtop",
                "inputs/segment.mdin",
            )
        )
    pinned_paths.extend(
        path for path in (branch / "native_source").rglob("*") if path.is_file()
    )
    pinned_paths.extend(Path(path) for path in post_hashes)
    pinned_paths.extend(
        path for name in names for path in (branch / name / "inputs").iterdir()
    )
    pinned = {str(path.resolve()): digest(path) for path in pinned_paths}
    plan = {
        "status": "PREPARED_ONLY_NOT_LAUNCHED",
        "requires_root_review": True,
        "branch": spec["name"],
        "base": str(branch),
        "original_base": str(original),
        "root": str(ROOT),
        "checkpoint_segment": f"extend_{checkpoint:02d}",
        "checkpoint_restart": str(fork_restart),
        "checkpoint_restart_sha256": digest(fork_restart),
        "shared_prefix_segments": [f"extend_{i:02d}" for i in range(1, checkpoint + 1)],
        "shared_prefix_is_not_independent_sampling": True,
        "fork_time_ps": 6 + checkpoint * 5,
        "new_duration_ps": 25,
        "target_horizon_ps": 6 + spec["last"] * 5,
        "segment_names": names,
        "seeds": seeds,
        "cpu_threads": 1,
        "native_source_sha256": native_hashes,
        "postprocess_source_sha256": post_hashes,
        "pinned_sha256": pinned,
        "candidate_gate": {
            "O_C_A_max": 3.5,
            "attack_angle_deg": [95, 125],
            "stop_on_candidate": True,
        },
        "statistics_scope": "Use only branch-local extend blocks for independent sampling statistics; include the shared prefix only for chronological geometry screening.",
        "water_handling": "All waters remain MM. A candidate is only a position observation; chemistry requires separate QM promotion and reoptimization.",
    }
    write_json(branch / "replica_plan.json", plan, exclusive=True)
    plan_hash = digest(branch / "replica_plan.json")
    (branch / "replica_plan.sha256").write_text(plan_hash + "\n")
    return {
        "branch": str(branch),
        "plan": str(branch / "replica_plan.json"),
        "plan_sha256": plan_hash,
        "checkpoint": checkpoint,
        "shared_prefix_ps": plan["fork_time_ps"],
        "segments": names,
    }


def candidate_check(branch: Path, native_source: Path, count: int, python: str) -> bool:
    """Run the frozen candidate gate in a fresh process using branch code."""
    command = [
        python,
        "-c",
        "import sys; from pathlib import Path; "
        f"sys.path.insert(0, {str(native_source)!r}); "
        "from nyl12_covalent.hydration_followon_gates import report_candidates; "
        f"print('candidate' if report_candidates(Path({str(branch)!r}), {count}) else 'clear')",
    ]
    env = dict(
        os.environ,
        PYTHONPATH=str(native_source),
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    result = subprocess.run(
        command, cwd=ROOT, env=env, capture_output=True, text=True, check=True
    )
    return result.stdout.strip().splitlines()[-1] == "candidate"
