"""Explicit, archived continuation of an existing direct-acyl diagnostic."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import shutil
from typing import Any


def prepare_resume(
    args: argparse.Namespace,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Validate a retained checkpoint, then archive the rejected suffix.

    Args:
        args: Runner arguments, including an explicit checkpoint digest.

    Returns:
        The unchanged remaining schedule and status with the accepted prefix.

    Raises:
        ValueError: If checkpoint, status or schedule contradict the continuation.
        FileExistsError: If the rejected-attempt archive already exists.
    """
    output = args.output
    inputs_path, status_path = output / "inputs.json", output / "status.json"
    inputs = json.loads(inputs_path.read_text())
    previous = json.loads(status_path.read_text())
    rows = inputs["schedule"]
    labels = [row["label"] for row in rows]
    index = labels.index(args.resume_after) + 1
    drive = output / "01_coupled_coordinate_drive"
    checkpoint = drive / args.resume_after / "optimized.ncrst"
    if args.start.resolve() != checkpoint.resolve():
        raise ValueError("Resume start must be the explicit accepted-window checkpoint")
    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    if not args.checkpoint_sha256 or digest != args.checkpoint_sha256:
        raise ValueError("Resume checkpoint SHA256 mismatch or missing digest")
    if previous["state"] != "stopped":
        raise ValueError("Resume requires a stopped scientific status")
    if previous["completed_windows"][:index] != labels[:index]:
        raise ValueError("Existing completed windows do not match the schedule prefix")
    for label in labels[:index]:
        acceptance = json.loads((drive / label / "acceptance.json").read_text())
        if not acceptance["native_converged"] or not acceptance["chemistry_retained"]:
            raise ValueError(f"Rejected checkpoint prefix window: {label}")
    if not rows[index:]:
        raise ValueError("No remaining scheduled windows")
    if args.final_composite != rows[-1]["composite_target_A"]:
        raise ValueError("Resumption must preserve the original final composite target")
    archive = output.parent / "clutter" / "initial_window21_proton_detour"
    if archive.exists():
        raise FileExistsError(archive)
    if not (output.parent / "source_snapshot").is_dir():
        raise FileNotFoundError("Missing immutable original source snapshot")
    suffix = [drive / label for label in labels[index:] if (drive / label).exists()]
    # All validation precedes mutations. No files from accepted windows are changed.
    archive.mkdir(parents=True)
    shutil.copy2(status_path, archive / "status.json")
    shutil.copy2(inputs_path, archive / "inputs.json")
    shutil.copytree(output.parent / "source_snapshot", archive / "source_snapshot")
    for folder in suffix:
        shutil.move(str(folder), str(archive / folder.name))
    now = datetime.now(timezone.utc).isoformat()
    manifest = {
        "continued_utc": now,
        "checkpoint": str(checkpoint.resolve()),
        "checkpoint_sha256": digest,
        "resume_after": args.resume_after,
        "original_schedule_retained": True,
        "substrate_proton_guard": args.substrate_proton_guard,
        "guard": {
            "atoms_one_based": [19929, 13917],
            "r1_r2_r3_r4_A": [0.8, 0.95, 1.25, 1.5],
            "rk2_rk3_kcal_mol_A2": [50, 50],
        }
        if args.substrate_proton_guard
        else None,
        "bias_protocol": "Piecewise changed from resumed window onward"
        if args.substrate_proton_guard
        else "Original composite bias unchanged",
        "rejected_attempt": str(archive.resolve()),
        "rejected_suffix_excluded_from_continuous_path": True,
        "release": "All reaction biases removed together at acyl-like candidate",
        "schedule": rows[index:],
    }
    (output / "continuation.json").write_text(json.dumps(manifest, indent=2) + "\n")
    status = {
        "state": "running",
        "completed_windows": labels[:index],
        "started_utc": now,
        "continuation": manifest,
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n")
    return rows[index:], status
