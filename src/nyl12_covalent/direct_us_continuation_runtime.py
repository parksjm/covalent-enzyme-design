"""Capture and verify DFTB3 parameter, engine, and topology identities."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import neb_umbrella_campaign as campaign

ELEMENTS = ("C", "H", "N", "O")
RELAY_PREPARED = Path(
    "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/"
    "full_path_300K/prepared.json"
)


def capture_external_runtime(cfg: dict[str, Any]) -> dict[str, str]:
    """Pin the current 3ob-3-1 SKF set against the relay campaign reference.

    Args:
        cfg: Frozen direct-US campaign configuration.

    Returns:
        SHA-256 map for all 16 SKF files, Sander executable, and topology.

    Raises:
        FileNotFoundError: If any runtime prerequisite is missing.
        ValueError: If any current SKF identity differs from the relay pin.
    """
    relay = campaign.ROOT / RELAY_PREPARED
    relay_cfg = json.loads(relay.read_text())
    expected = relay_cfg.get("external_runtime_sha256")
    if not isinstance(expected, dict):
        raise ValueError("Relay prepared campaign lacks external_runtime_sha256")
    parameter_dir = Path(cfg["sander"]).parent.parent / "dat/slko/3ob-3-1"
    identities = {}
    for first in ELEMENTS:
        for second in ELEMENTS:
            path = parameter_dir / f"{first}-{second}.skf"
            name = str(path.resolve(strict=True))
            digest = campaign.digest(Path(name))
            if expected.get(name) != digest:
                raise ValueError(
                    f"Current 3ob-3-1 parameter differs from relay pin: {name}"
                )
            identities[name] = digest
    for path in (Path(cfg["sander"]), Path(cfg["topology"])):
        resolved = path.resolve(strict=True)
        identities[str(resolved)] = campaign.digest(resolved)
    return identities


def verify_external_runtime(identities: dict[str, str]) -> None:
    """Fail if any pinned SKF, engine, or topology file changed."""
    for name, expected in identities.items():
        path = Path(name)
        if campaign.digest(path) != expected:
            raise ValueError(f"Pinned external runtime file changed: {path}")
