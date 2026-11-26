"""Verify explicitly reviewed native heating-checkpoint adoption receipts."""

from __future__ import annotations
import json
from pathlib import Path
from typing import Any
from .neb_umbrella_campaign import digest


def heat_adoption(
    plan: dict[str, Any], label: str, stage: Path
) -> tuple[Path, Path] | None:
    """Return an accepted heat receipt and original plan only when pinned exactly.

    Args:
        plan: New immutable execution plan, optionally with explicit adoptions.
        label: Existing scientific window identity.
        stage: Original heat-stage directory containing preserved native outputs.

    Returns:
        Separate reviewed acceptance path and original input-plan path, or None
        when this window has no explicit adoption.
    """
    if "adopted_heat" not in plan or label not in plan["adopted_heat"]:
        return None
    entry = plan["adopted_heat"][label]
    path = Path(entry["receipt"])
    if digest(path) != entry["sha256"]:
        raise ValueError("Reviewed heat adoption receipt changed")
    receipt = json.loads(path.read_text())
    if receipt["window"] != label or Path(receipt["stage"]) != stage:
        raise ValueError("Reviewed heat adoption window/stage mismatch")
    if receipt["temperature_acceptance_role"] != "completed_nonstationary_heating_ramp":
        raise ValueError("Heat adoption must not assert thermal equilibrium")
    for filename, expected in receipt["native_source_sha256"].items():
        if digest(Path(filename)) != expected:
            raise ValueError(f"Reviewed heat native source changed: {filename}")
    original = Path(receipt["original_input_plan"])
    if digest(original) != receipt["original_input_plan_sha256"]:
        raise ValueError("Original heat input plan changed")
    return path, original
