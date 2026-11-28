"""Pinned campaign source identities and stable JSON snapshot utilities."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


BRANCHES = {
    "01_activation_reverse": "reverse_activation",
    "02_acylation_forward": "forward_acylation",
}


ROLES = {
    "carbon": 19916,
    "leaving_n": 19929,
    "og1": 13916,
    "hg1": 13917,
    "nalpha": 13905,
    "nalpha_h1": 13906,
    "nalpha_h2": 13907,
}


TRACKED_H = {13906, 13907, 13917, 19937}


DISCONTINUITY_A = 0.35


class SnapshotChangedError(RuntimeError):
    """Raised when a live campaign record changes during an audit pass."""


def digest(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    """Read an authoritative JSON record and require an object root."""
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _read_snapshot(path: Path) -> tuple[dict[str, Any], str]:
    """Read one JSON document and return its in-memory snapshot digest."""
    payload = path.read_bytes()
    snapshot_hash = hashlib.sha256(payload).hexdigest()
    try:
        value = json.loads(payload)
    except json.JSONDecodeError as exc:
        if digest(path) != snapshot_hash:
            raise SnapshotChangedError(f"Live JSON changed while read: {path}") from exc
        raise
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value, snapshot_hash


def _verify_source(campaign: Path, config: dict[str, Any]) -> dict[str, Any]:
    model = campaign / config["model"]
    manifest_path = model / "manifest.json"
    manifest = read_json(manifest_path)
    preparation_path = campaign / "preparation.json"
    preparation = read_json(preparation_path)
    if preparation.get("source_frozen_manifest_sha256") != digest(manifest_path):
        raise ValueError("Preparation manifest does not match frozen model manifest")
    config_hash = digest(campaign / "config.json")
    optimized_hash = manifest["inputs/optimized.ncrst"]
    checked: list[dict[str, str]] = []
    checked.append({"path": str(manifest_path), "sha256": digest(manifest_path)})
    checked.append({"path": str(preparation_path), "sha256": digest(preparation_path)})
    required_inputs = {
        "inputs/system.prmtop",
        "inputs/dimer.mdin",
        "inputs/optimized.ncrst",
        "inputs/scientific_acceptance.json",
        "source_identity.json",
    }
    if not required_inputs.issubset(manifest):
        raise ValueError("Frozen model manifest is missing a required source identity")
    for relative, expected in manifest.items():
        path = model / relative
        observed = digest(path)
        if observed != expected:
            raise ValueError(f"Frozen model artifact changed: {relative}")
        checked.append({"path": str(path), "sha256": observed})
    topology_hash = manifest["inputs/system.prmtop"]
    source_identity_path = model / "source_identity.json"
    source_identity = read_json(source_identity_path)
    source_identity_hash = digest(source_identity_path)
    if source_identity_hash != manifest["source_identity.json"]:
        raise ValueError("Source identity record differs from frozen model manifest")
    checked.append(
        {
            "path": str(source_identity_path),
            "sha256": source_identity_hash,
            "coordinate_sha256_semantics": "pre-optimization trajectory seed; distinct from optimized restart",
            "coordinate_sha256": source_identity["coordinate_sha256"],
        }
    )
    for directory in BRANCHES:
        plan_path = campaign / directory / "exploratory_plan_20261124.json"
        plan = read_json(plan_path)
        if plan.get("config_sha256") != config_hash:
            raise ValueError(f"Branch plan config hash mismatch: {plan_path}")
        if plan.get("source_restart_sha256") != optimized_hash:
            raise ValueError(
                f"Branch plan optimized restart hash mismatch: {plan_path}"
            )
        if plan.get("topology_sha256") != topology_hash:
            raise ValueError(f"Branch plan topology hash mismatch: {plan_path}")
        authorization_path = campaign / directory / plan["authorization_file"]
        if digest(authorization_path) != plan["authorization_sha256"]:
            raise ValueError(
                f"Branch authorization hash mismatch: {authorization_path}"
            )
        checked.append({"path": str(plan_path), "sha256": digest(plan_path)})
        checked.append(
            {"path": str(authorization_path), "sha256": digest(authorization_path)}
        )
    snapshot_hashes = config["source_sha256"]
    for path_text, expected in snapshot_hashes.items():
        path = Path(path_text)
        if not (
            "source_snapshot" in path.parts
            or path.suffix == ".skf"
            or path.name == "sander"
        ):
            continue
        observed = digest(path)
        if observed != expected:
            raise ValueError(f"Frozen source snapshot changed: {path}")
        checked.append({"path": str(path), "sha256": observed})
    return {
        "manifest_path": str(manifest_path),
        "checked": checked,
        "base_hamiltonian_input_sha256": manifest["inputs/dimer.mdin"],
        "coordinate_sha256_semantics": "source_identity coordinate hash describes the pre-optimization trajectory seed; optimized restart identity is in model/manifest.json and branch plans",
    }
