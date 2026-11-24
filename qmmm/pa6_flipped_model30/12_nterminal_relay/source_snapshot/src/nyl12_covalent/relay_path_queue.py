"""Fail-closed gate for a relay pathway conditioned on full-shell inertia."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from .relay_minimum import verify_frozen


def digest(path: Path) -> str:
    """Return the SHA-256 digest of an exact input or receipt.

    Args:
        path: Existing file to hash.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_minimum(root: Path) -> dict[str, Any] | None:
    """Require consistent full-shell scientific and native-process acceptance.

    Args:
        root: Existing window033 minimum-verification directory.

    Returns:
        Hash-bound receipt on acceptance, or None while work is pending.

    Raises:
        ValueError: If a finished or scientific result fails the explicit gate.
        FileNotFoundError: If required accepted-state evidence is missing.
    """
    frozen = root / "01_frozen"
    verify_frozen(frozen)
    frozen_manifest = json.loads((frozen / "manifest.json").read_text())
    required = {
        "inputs/optimized.ncrst",
        "inputs/system.prmtop",
        "inputs/dimer.mdin",
        "library_build/libdlfind.so",
    }
    if not required <= frozen_manifest.keys():
        raise ValueError("Frozen manifest omits required native inputs")
    stage = root / "06_full_hessian"
    process_path = stage / "job/status.json"
    summary_path = root / "verification_summary.json"
    process = json.loads(process_path.read_text()) if process_path.is_file() else None
    if (
        process
        and process.get("ended_utc") is not None
        and process.get("exit_status") != 0
    ):
        raise ValueError("Full-Hessian native job failed")
    if not summary_path.is_file():
        if process and process.get("ended_utc") is not None:
            raise ValueError("Full-Hessian job ended without a scientific summary")
        return None
    summary = json.loads(summary_path.read_text())
    status = summary["status"]
    if status in {"incomplete", "not_verified", "inconclusive_preflight"}:
        raise ValueError(f"Minimum verification failed: {status}")
    if status == "full_hessian_running":
        if process and process.get("ended_utc") is not None:
            raise ValueError(
                "Full-Hessian job ended without a terminal scientific result"
            )
        return None
    if (
        status != "numerical_local_minimum"
        or summary.get("numerically_verified") is not True
    ):
        raise ValueError(f"Unexpected or rejected verification status: {status}")
    inertia = summary["literal_inertia"]
    expected = {
        "passed": True,
        "negative": 0,
        "unresolved": 0,
        "expected_negative": 0,
        "dimension": 12495,
    }
    if any(
        inertia.get(k) != v or type(inertia.get(k)) is not type(v)
        for k, v in expected.items()
    ):
        raise ValueError("Full-Hessian inertia fails minimum gate")
    if not process or process.get("ended_utc") is None:
        return None
    acceptance_path = stage / "scientific_acceptance.json"
    if json.loads(acceptance_path.read_text()) != inertia:
        raise ValueError("Summary inertia differs from native scientific acceptance")
    specification_path = stage / "specification.json"
    spec = json.loads(specification_path.read_text())
    restart = frozen / "inputs/optimized.ncrst"
    topology = frozen / "inputs/system.prmtop"
    if (
        spec.get("dimension") != 12495
        or spec.get("restart_sha256") != digest(restart)
        or spec.get("library_sha256") != digest(frozen / "library_build/libdlfind.so")
    ):
        raise ValueError(
            "Hessian specification differs from frozen restart or dimension"
        )
    return {
        "summary_sha256": digest(summary_path),
        "hessian_acceptance_sha256": digest(acceptance_path),
        "specification_sha256": digest(specification_path),
        "native_job_status_sha256": digest(process_path),
        "frozen_manifest_sha256": digest(frozen / "manifest.json"),
        "restart_sha256": digest(restart),
        "topology_sha256": digest(topology),
        "raw_hessian_sha256": inertia["raw_hessian_sha256"],
        "symmetric_hessian_sha256": inertia["symmetric_hessian_sha256"],
        "literal_inertia": inertia,
    }
