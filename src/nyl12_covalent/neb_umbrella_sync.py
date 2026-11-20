"""Validate and publish compact remote production results with acceptance last."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import os
from pathlib import Path
import re
import shutil

import numpy as np

from .neb_umbrella_dispatch import digest, read, validate_plan, write

CRITICAL = (
    "samples_and_cross_bias.npz",
    "analysis_request.json",
    "observables/observables.json",
    "observables/observables.csv",
    "process.json",
    "acceptance.json",
)


def context(path: Path, remote: str) -> tuple[dict, dict, Path, dict[str, list[str]]]:
    """Verify immutable plan/config and expose only declared remote production stages."""
    if digest(path) != (path.parent / "plan.sha256").read_text().strip():
        raise ValueError("Dispatch plan changed")
    plan = read(path)
    validate_plan(plan)
    if remote != plan["remote_host"] or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9.-]*", remote
    ):
        raise ValueError("Remote must be the explicit remote owner in the plan")
    prepared = Path(plan["hosts"][plan["local_host"]]["prepared"])
    if digest(prepared) != plan["hosts"][plan["local_host"]]["sha256"]:
        raise ValueError("Master prepared identity changed")
    if Path(plan["hosts"][remote]["prepared"]).parent != prepared.parent:
        raise ValueError(
            "Sync requires identical canonical campaign paths on both hosts"
        )
    cfg = read(prepared)
    owned = {task["label"] for tasks in plan["assignments"][remote] for task in tasks}
    blocks = {
        row["label"]: [s["label"] for s in row["stages"] if s["kind"] == "production"]
        for row in cfg["windows"]
        if row["label"] in owned
    }
    if set(blocks) != owned or any(not re.fullmatch(r"window_\d{3}", w) for w in owned):
        raise ValueError("Remote window ownership differs from master")
    if any(
        not re.fullmatch(r"03_production_block_\d{2}_5ps", s)
        for ss in blocks.values()
        for s in ss
    ):
        raise ValueError("Unexpected production stage identity")
    return plan, cfg, prepared, blocks


def verify_hashes(folder: Path, hashes: dict[str, str]) -> None:
    """Reject missing, linked, changed, absolute or traversing compact payload paths."""
    for name, expected in hashes.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise ValueError("Unsafe transferred relative path")
        path = folder / relative
        if (
            any(
                p.is_symlink() for p in [path, *path.parents[: len(relative.parts) - 1]]
            )
            or not path.is_file()
            or digest(path) != expected
        ):
            raise ValueError(f"Transferred payload hash mismatch: {path}")


def validate_block(folder: Path, hashes: dict[str, str], cfg: dict) -> None:
    """Validate all frames, CSV identity, finite NPZ arrays and complete native process."""
    if set(hashes) != set(CRITICAL):
        raise ValueError("Completed block must contain exactly the critical payloads")
    verify_hashes(folder, hashes)
    receipt, request = (
        read(folder / "observables/observables.json"),
        read(folder / "analysis_request.json"),
    )
    accepted, process = read(folder / "acceptance.json"), read(folder / "process.json")
    if (
        accepted["finite_energy_temperature"] is not True
        or accepted["completed_steps"] <= 0
        or process["exit_status"] != 0
        or process["ended_utc"] is None
    ):
        raise ValueError("Production process/acceptance is incomplete or failed")
    csv_path = folder / "observables/observables.csv"
    if digest(csv_path) != receipt["csv_sha256"] or receipt["frames_censored"] != 0:
        raise ValueError("Observable CSV identity or censoring contract changed")
    with csv_path.open() as stream:
        rows = list(csv.DictReader(stream))
    samples = np.array(
        [
            [float(row[k]) for k in ("q_acyl", "q_h", "cn_progress", "cn_nalpha_hpool")]
            for row in rows
        ]
    )
    if (
        not rows
        or len(rows) != receipt["frames"]
        or request["samples"] != len(rows)
        or not np.isfinite(samples).all()
    ):
        raise ValueError("Invalid observable frame count or nonfinite CVs")
    shape = (len(cfg["windows"]), len(rows))
    if request["bias_matrix_shape"] != list(shape) or request[
        "all_window_bias_definitions"
    ] != [w["bias"] for w in cfg["windows"]]:
        raise ValueError(
            "Remote bias definitions or matrix dimensions differ from master"
        )
    with np.load(folder / "samples_and_cross_bias.npz", allow_pickle=False) as archive:
        cv = archive["cv_values"]
        if cv.shape != (len(rows), 2) or not np.allclose(
            cv, samples[:, :2], rtol=0, atol=1e-8
        ):
            raise ValueError("NPZ/observable frame identity mismatch")
        for name in ("bias_kcal_mol", "u_kn_bias_only"):
            values = archive[name]
            if values.shape != shape or not np.isfinite(values).all():
                raise ValueError("Invalid/nonfinite cross-bias matrix")


def identity(plan_path: Path, remote: str, source: str, hashes: dict) -> dict:
    """Describe immutable remote identity separately from observation time."""
    plan = read(plan_path)
    return {
        "source_host": remote,
        "source_path": source,
        "plan_sha256": digest(plan_path),
        "remote_prepared_sha256": plan["hosts"][remote]["sha256"],
        "file_sha256": hashes,
    }


def published(
    folder: Path, hashes: dict, expected: dict, gate: str = "acceptance.json"
) -> bool:
    """Require identical previously accepted payloads and their transfer provenance."""
    if not (folder / gate).exists():
        return False
    verify_hashes(folder, hashes)
    receipt = read(folder / "transfer_receipt.json")
    if any(receipt.get(key) != value for key, value in expected.items()):
        raise ValueError("Accepted remote transfer identity changed")
    return True


def atomic_copy(source: Path, destination: Path) -> None:
    """Publish one verified compact file through an adjacent atomic replacement."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    pending = destination.with_name(f".{destination.name}.sync_{os.getpid()}")
    shutil.copyfile(source, pending)
    pending.replace(destination)


def publish_block(
    incoming: Path, target: Path, hashes: dict, cfg: dict, expected: dict
) -> bool:
    """Commit acceptance only after every validated payload and provenance are present."""
    validate_block(incoming, hashes, cfg)
    if published(target, hashes, expected):
        return False
    for name in CRITICAL:
        existing = target / name
        if existing.exists() and digest(existing) != hashes[name]:
            raise ValueError(f"Existing partial block payload changed: {existing}")
    for name in CRITICAL[:-1]:
        atomic_copy(incoming / name, target / name)
    write(
        target / "transfer_receipt.json",
        {**expected, "transferred_utc": datetime.now(timezone.utc).isoformat()},
    )
    atomic_copy(incoming / "acceptance.json", target / "acceptance.json")
    return True


def manifest_hashes(folder: Path) -> dict[str, str]:
    """Parse a complete immutable endpoint-view bundle manifest."""
    raw = (folder / "MANIFEST.sha256").read_text()
    if not raw or not raw.endswith("\n"):
        raise ValueError("Incomplete view bundle manifest")
    hashes = {}
    for line in raw.splitlines():
        expected, name = line.split("  ", 1)
        if (
            not re.fullmatch(r"[0-9a-f]{64}", expected)
            or Path(name).name != name
            or name in hashes
            or name in (".", "..", "MANIFEST.sha256")
        ):
            raise ValueError("Invalid view bundle manifest")
        hashes[name] = expected
    if not {"manifest.json", "view.pdb", "view_umbrella.vmd"} <= set(hashes):
        raise ValueError("Incomplete view bundle payload listing")
    return hashes


def publish_view(incoming: Path, target: Path, hashes: dict, expected: dict) -> bool:
    """Publish a whole validated immutable compact view directory atomically."""
    verify_hashes(incoming, hashes)
    entries = manifest_hashes(incoming)
    if set(hashes) != {*entries, "MANIFEST.sha256"}:
        raise ValueError("View inventory/manifest payload sets differ")
    verify_hashes(incoming, entries)
    if published(target, hashes, expected, "MANIFEST.sha256"):
        return False
    if target.exists():
        raise FileExistsError("Refuse to replace an incomplete existing view directory")
    target.parent.mkdir(parents=True, exist_ok=True)
    pending = target.parent / f".{target.name}.sync_{os.getpid()}"
    pending.mkdir(exist_ok=False)
    for name in hashes:
        shutil.copyfile(incoming / name, pending / name)
    write(
        pending / "transfer_receipt.json",
        {**expected, "transferred_utc": datetime.now(timezone.utc).isoformat()},
    )
    pending.rename(target)
    return True


def publish_latest(generation: Path) -> None:
    """Point latest at a fully published compact view using an atomic symlink swap."""
    if not (generation / "MANIFEST.sha256").is_file():
        raise ValueError("Cannot point latest at an incomplete view")
    pointer = generation.parent / "latest"
    if pointer.exists() and not pointer.is_symlink():
        raise ValueError("Refuse to replace a non-symlink latest artifact")
    pending = generation.parent / f".latest_sync_{os.getpid()}"
    pending.symlink_to(generation.name)
    pending.replace(pointer)
