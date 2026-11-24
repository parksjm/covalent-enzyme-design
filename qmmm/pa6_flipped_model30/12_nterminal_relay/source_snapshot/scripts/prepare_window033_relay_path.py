"""Freeze the two relay branches without executing new quantum calculations."""

from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import importlib.util
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.relay_minimum import verify_frozen  # noqa: E402
from nyl12_covalent.stationary_native import NativeState, digest, write_json  # noqa: E402

LOGGER = logging.getLogger(__name__)


def freeze_sources(destination: Path) -> dict[str, str]:
    """Copy the static local import closure of both persistent entrypoints.

    Args:
        destination: New private source snapshot root.

    Returns:
        Absolute snapshot paths with their content SHA-256 values.
    """
    pending = [
        ROOT / "scripts/watch_window033_relay.py",
        ROOT / "scripts/run_window033_relay_path.py",
        Path(__file__),
        ROOT / "scripts/record_agent_action.py",
    ]
    seen: set[Path] = set()
    while pending:
        source = pending.pop()
        if source in seen:
            continue
        if not source.is_file():
            raise FileNotFoundError(source)
        seen.add(source)
        relative = source.relative_to(ROOT)
        package = ".".join(relative.with_suffix("").parts[:-1])
        if package.startswith("src."):
            package = package[4:]
        tree = ast.parse(source.read_text())
        modules = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                module = "." * node.level + (node.module or "")
                modules.append(
                    importlib.util.resolve_name(module, package)
                    if node.level
                    else module
                )
            elif isinstance(node, ast.Import):
                modules.extend(item.name for item in node.names)
        for module in modules:
            if module.startswith("nyl12_covalent."):
                candidate = ROOT / "src" / (module.replace(".", "/") + ".py")
            elif module.startswith("scripts."):
                candidate = ROOT / (module.replace(".", "/") + ".py")
            else:
                continue
            if candidate.is_file():
                pending.append(candidate)
            else:
                raise FileNotFoundError(f"Unresolved local source import: {module}")
    seen.add(ROOT / "src/nyl12_covalent/__init__.py")
    result = {}
    for source in sorted(seen):
        target = destination / source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        result[str(target)] = digest(target)
    return result


def prepare(seed: Path, output: Path) -> dict[str, Any]:
    """Freeze exact source state, display mapping and guarded relay execution.

    Args:
        seed: Existing force-accepted window033 release.
        output: Fresh logically numbered relay campaign.

    Returns:
        Written configuration; quantum execution remains forbidden until gating.
    """
    seed = seed.resolve(strict=True)
    verification = seed / "minimum_verification"
    frozen = verification / "01_frozen"
    verify_frozen(frozen)
    release = json.loads((frozen / "release_acceptance.json").read_text())
    if release["accepted"] is not True or release["force_converged"] is not True:
        raise ValueError("Original window033 release is not force-accepted")
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    model = output / "model"
    shutil.copytree(frozen, model)
    state = NativeState(model)
    viewer = output / "viewer_template"
    viewer.mkdir()
    for name in ("stripped.parm7", "view_relay.vmd", "manifest.json"):
        shutil.copy2(seed / "vmd_progress_verified" / name, viewer / name)
    viewer_record = json.loads((viewer / "manifest.json").read_text())
    if viewer_record["topology_sha256"] != digest(state.topology):
        raise ValueError("Viewer uses a different source topology")
    sources = freeze_sources(output / "source_snapshot")
    for directory in (model, viewer):
        sources.update({str(p): digest(p) for p in directory.rglob("*") if p.is_file()})
    source_campaign = seed.parent
    parameters = json.loads((source_campaign / "prepared.json").read_text())[
        "parameter_sha256"
    ]
    amber = Path("/data/jpc/code/ambertools26/venv")
    for name, expected in parameters.items():
        path = amber / "dat/slko/3ob-3-1" / name
        if digest(path) != expected:
            raise ValueError(f"Original 3ob parameter changed: {path}")
        sources[str(path)] = expected
    for binary in (state.sander, amber / "bin/cpptraj"):
        sources[str(binary)] = digest(binary)
    config = {
        "verification_root": str(verification),
        "model": "model",
        "viewer_template": "viewer_template",
        "cpptraj": str(amber / "bin/cpptraj"),
        "source_topology_sha256": digest(state.topology),
        "reference_energy_hartree": release["final_energy_hartree"],
        "source_sha256": sources,
        "proton": 13917,
        "branches": {
            "reverse_activation": "01_activation_reverse",
            "forward_acylation": "02_acylation_forward",
        },
        "path_final_A": 4.2,
        "path_maxincrement_A": 0.15,
        "path_maxcycles": 4000,
        "scope": "Two native single-CPU branches only after qualified full-shell minimum acceptance",
        "reaction_bias": "One asynchronous composite per branch; no positional or proton-identity restraints",
        "fixed_exterior": "Identical window033 coordinates in both branches and endpoints",
        "ts_search": "Strict interior physical-PES maximum; native dimer up to5000cycles, then independent force check; index verification separate",
    }
    write_json(output / "config.json", config)
    write_json(
        output / "preparation.json",
        {
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
            "source_seed": str(seed),
            "source_frozen_manifest_sha256": digest(frozen / "manifest.json"),
            "config_sha256": digest(output / "config.json"),
            "repository_revision": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
            "native_calculations_launched": False,
            "minimum_classification": "Pending actual full-Hessian scientific/process gate",
        },
    )
    LOGGER.info("Prepared %s; no quantum calculation launched", output)
    return config


def main() -> None:
    """Prepare explicit original seed into a fresh campaign without dispatch."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    prepare(args.seed, args.output)


if __name__ == "__main__":
    main()
