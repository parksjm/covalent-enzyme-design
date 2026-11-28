"""Frozen70ps main100 plus six independent reactant-boundary histories."""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from collect import collect
from nyl12_covalent.relay_pooled_sources import load_cohort, checkpoint_receipts
from nyl12_covalent.relay_umbrella_analysis_io import digest

WORK = Path(__file__).resolve().parent.parent


def registration() -> dict:
    """Load explicit immutable baseline and boundary identities."""
    return json.loads((WORK / "registration.json").read_text())


def ready(horizon: int) -> bool:
    """Require every consecutive03 onward production receipt for all six starts."""
    config = registration()
    if horizon not in config["minimum_new_production_checkpoints_ps"]:
        raise ValueError("Unauthorized boundary horizon")
    prepared = Path(config["boundary_prepared"])
    if digest(prepared) != config["boundary_prepared_sha256"]:
        raise ValueError("Boundary manifest identity changed")
    windows = json.loads(prepared.read_text())["windows"]
    if len(windows) != 6 or len({w["label"] for w in windows}) != 6:
        raise ValueError("Require six distinct histories")
    return all(
        checkpoint_receipts(prepared.parent, w["label"], horizon) for w in windows
    )


def combine(old: dict, new: dict, horizon: int) -> dict:
    """Append six histories without merging duplicate bias centers or pruning data."""
    if len(old["states"]) != 100 or len(new["states"]) != 6:
        raise ValueError("Require exact100+6 history inventory")
    if new["counts"].tolist() != [(horizon - 5) * 20] * 6:
        raise ValueError("Expected once-only5ps production discard")
    states = old["states"] + new["states"]
    proofs = old["proofs"] + new["proofs"]
    csv = [
        next(v for k, v in p["source_sha256"].items() if k.endswith(".csv"))
        for p in proofs
    ]
    if len(set(csv)) != len(csv) or len({s["state_id"] for s in states}) != 106:
        raise ValueError("Duplicate physical block or history")
    biases = {json.dumps(s["bias"], sort_keys=True) for s in states}
    if len(biases) != 95:
        raise ValueError("Expected95 distinct physical bias potentials")
    return dict(
        old,
        q=np.concatenate([old["q"], new["q"]]),
        counts=np.concatenate([old["counts"], new["counts"]]),
        states=states,
        proofs=proofs,
        bootstrap_eligible=False,
        frozen_baseline_main_horizon_ps=70,
        new_horizon_ps=horizon,
    )


def collect_boundary(horizon: int, output: Path) -> dict:
    """Reuse unchanged native-parity loader and frozen100-history collector."""
    config = registration()
    selected = json.loads(Path(config["baseline_inventory"]).read_text())
    old = collect(selected, output)
    new = load_cohort(
        Path(config["boundary_prepared"]),
        horizon,
        old["config"],
        "reactant_boundary6",
        6,
    )
    if new is None:
        raise ValueError("Complete accepted boundary prefix disappeared")
    return combine(old, new, horizon)
