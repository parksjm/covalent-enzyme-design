"""Cheap, bounded checks for regenerated portable viewers with pinned chemistry."""

from pathlib import Path
import hashlib
import json
from typing import Any

from .vmd_movie_policy import POLICY_VERSION
from .vmd_style_manifests import digest, verify_checksums
from .vmd_style_refresh import BEGIN, END, AmberParm


def inspect_viewer(row: dict[str, Any], state: dict[str, Any]) -> str:
    """Inspect one logical movie, allowing refresh only after a stable poll.

    Args:
        row: Explicit refresh mapping with logical viewer/topology and pinned hash.
        state: Mutable per-movie observed/accepted signature state.

    Returns:
        ``unchanged``, ``unstable``, ``styled``, ``waiting_checksums`` or
        ``refresh``. Cheap metadata checks avoid loading unchanged topologies.

    Raises:
        ValueError: If a stable topology no longer matches its validated hash.
    """
    viewer, topology = Path(row["viewer"]), Path(row["topology"])
    if not viewer.is_file() or not topology.is_file():
        return "waiting_files"
    signature = (
        str(viewer.resolve()),
        viewer.stat().st_size,
        viewer.stat().st_mtime_ns,
        str(topology.resolve()),
        topology.stat().st_size,
        topology.stat().st_mtime_ns,
    )
    if state.get("accepted") == signature:
        return "unchanged"
    if state.get("observed") != signature:
        state["observed"] = signature
        return "unstable"
    actual = digest(topology)
    if actual != row["topology_sha256"]:
        identity = topology_identity(AmberParm(str(topology)))
        if identity != row["topology_identity_sha256"]:
            raise ValueError(
                "HOLD: logical movie topology identity changed; physical mapping must be revalidated"
            )
        state["equivalent_topology"] = dict(
            pinned_raw_sha256=row["topology_sha256"],
            actual_raw_sha256=actual,
            identity_sha256=identity,
            reason="Raw export bytes changed; exact atom identities, charges, types and native bond graph retained",
        )
    text = viewer.read_text()
    if (
        BEGIN in text
        and END in text
        and f"# Display policy version: {POLICY_VERSION}" in text
    ):
        state["accepted"] = signature
        return "styled"
    if verify_checksums(viewer.parent):
        return "waiting_checksums"
    return "refresh"


def topology_identity(topology: Any) -> str:
    """Hash ordered physical atom identities, charges, types and native bonds.

    Args:
        topology: Parsed physical ParmEd topology.

    Returns:
        Exact identity digest insensitive to exported header dates.
    """
    atoms = [
        (
            a.name,
            a.atomic_number,
            a.residue.idx,
            a.residue.name,
            a.residue.chain,
            a.type,
            a.charge,
            a.mass,
        )
        for a in topology.atoms
    ]
    bonds = sorted(tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in topology.bonds)
    return hashlib.sha256(
        json.dumps(dict(atoms=atoms, bonds=bonds), separators=(",", ":")).encode()
    ).hexdigest()
