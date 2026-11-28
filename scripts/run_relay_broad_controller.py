#!/usr/bin/env python3
"""Persist a finite, fail-closed RFD3-to-broad-pilot continuation in tmux."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys
import time
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.relay_broad_runtime import (  # noqa: E402
    digest,
    journal,
    native_run,
    verify_parent,
    write_new,
)

LOGGER = logging.getLogger(__name__)


def read_json(path: Path) -> dict[str, Any]:
    """Read one required JSON object."""
    return json.loads(path.read_text())


def wait_for_gate(
    gate: Path, terminal: Path, deadline: float, sleep_seconds: float = 60
) -> dict[str, Any]:
    """Wait for acceptance, stopping on final rejection or finite timeout."""
    while time.monotonic() < deadline:
        value = read_json(gate)
        if value["accepted"] is True:
            if len(value["parents"]) != 2:
                raise ValueError(
                    "Initial qualification must accept exactly two parents"
                )
            for parent in value["parents"]:
                verify_parent(parent)
            return value
        if terminal.exists():
            end = read_json(terminal)
            if end["terminal"] is True and end["accepted"] is False:
                raise RuntimeError(
                    "Final qualification rejected; no continuation authorized"
                )
        time.sleep(sleep_seconds)
    raise TimeoutError("Finite qualification waiting budget exhausted")


def check_hashes(hashes: dict[str, str]) -> None:
    """Reject changed frozen sources or immutable inputs."""
    for path, expected in hashes.items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Frozen controller source/input changed: {path}")


def execute(config: dict[str, Any]) -> None:
    """Run expansion, bounded generation/crossscores and exactly eight predictions."""
    output = Path(config["output"])
    check_hashes(config["source_sha256"])
    journal(ROOT, "intent", "Wait for qualified broad-pilot gate", config["gate"])
    deadline = time.monotonic() + config["max_wait_seconds"]
    accepted = wait_for_gate(Path(config["gate"]), Path(config["terminal"]), deadline)
    check_hashes(config["source_sha256"])
    release_code = native_run(
        ROOT,
        config["release_only_argv"],
        output / "verify_expansion_release",
        config["rfd3_environment"],
    )
    if release_code:
        raise RuntimeError(
            "Independent expansion execution release verification failed"
        )
    release_path = Path(config["expansion_release"])
    while not release_path.exists() and time.monotonic() < deadline:
        time.sleep(60)
    release = read_json(release_path)
    if (
        release["execution_release_checked"] is not True
        or release["qualification_gate_verified"] is not True
        or release["output_scientific_acceptance"] is not False
    ):
        raise ValueError(
            "Expansion execution release rejected or misstates unrun outputs"
        )
    manifest_path = Path(release["manifest"])
    if digest(manifest_path) != release["manifest_sha256"]:
        raise ValueError("Expansion release manifest changed")
    expansion_manifest = read_json(manifest_path)
    # A newly frozen release is permitted only for the same six bounded targets.
    if len(expansion_manifest["jobs"]) != 6:
        raise ValueError("Expansion must contain exactly six remaining jobs")
    runner = Path(expansion_manifest["runner"])
    if digest(runner) != release["runner_sha256"]:
        raise ValueError("Expansion runner differs from reviewed release")
    write_new(output / "qualification_gate_acceptance.json", accepted)
    write_new(output / "expansion_release_frozen.json", release)
    journal(
        ROOT,
        "outcome",
        "Wait for qualified broad-pilot gate",
        "Two native parents accepted",
    )
    argv = [
        config["rfd3_python"],
        str(runner),
        str(manifest_path),
        "--phase",
        "expansion",
    ]
    code = native_run(ROOT, argv, output / "expansion", config["rfd3_environment"])
    expanded_path = Path(config["expanded_parents"])
    if not expanded_path.exists():
        raise RuntimeError(
            f"Expansion ended without terminal parent manifest: exit={code}"
        )
    expanded = read_json(expanded_path)
    if expanded["status"] not in {"accepted", "scientific_failed"}:
        raise RuntimeError(f"Expansion operational failure: {expanded['status']}")
    if code and expanded["status"] != "scientific_failed":
        raise RuntimeError(f"Expansion process failed: exit={code}")
    for parent in expanded["parents"]:
        verify_parent(parent)
    check_hashes(config["source_sha256"])
    native = read_json(Path(config["native_manifest"]))
    request = read_json(Path(config["request"]))
    broad_manifest = {
        "request": request,
        "native": native,
        "accepted_manifests": [accepted, expanded],
        "local_pocket_mask": config["local_pocket_mask"],
        "source_sha256": {
            **config["source_sha256"],
            config["gate"]: digest(Path(config["gate"])),
            str(expanded_path): digest(expanded_path),
            str(manifest_path): digest(manifest_path),
            str(runner): digest(runner),
        },
    }
    broad_path = output / "broad_execution_manifest.json"
    write_new(broad_path, broad_manifest)
    generation = output / "generation_run01"
    code = native_run(
        ROOT,
        [native["python"], config["broad_runner"], str(broad_path), str(generation)],
        output / "broad_pilot",
        native["environment"],
    )
    if (
        code
        or read_json(generation / "scientific_acceptance.json")["accepted"] is not True
    ):
        raise RuntimeError(f"Broad native pilot did not pass: exit={code}")
    check_hashes(config["source_sha256"])
    prediction_argv = [
        part.replace("{design_json}", str(generation / "preselected_candidate.json"))
        for part in config["prediction_argv"]
    ]
    code = native_run(
        ROOT, prediction_argv, output / "predictions", config["prediction_environment"]
    )
    if code:
        raise RuntimeError(f"Prediction comparison failed: exit={code}")
    write_new(
        output / "continuation_scientific_acceptance.json",
        {
            "accepted": True,
            "bounded_generation_crossscores_completed": True,
            "prediction_process_completed": True,
            "geometry_acceptance": "separate_receipt",
            "activity_validated": False,
        },
    )


def main() -> None:
    """Run once with immutable configuration and preserve a terminal receipt."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configuration", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    config = read_json(args.configuration)
    try:
        execute(config)
    except Exception as error:
        journal(ROOT, "outcome", "Broad continuation terminated", repr(error))
        write_new(
            Path(config["output"]) / "continuation_terminal.json",
            {"terminal": True, "accepted": False, "error": repr(error)},
        )
        raise
    write_new(
        Path(config["output"]) / "continuation_terminal.json",
        {"terminal": True, "accepted": True},
    )


if __name__ == "__main__":
    main()
