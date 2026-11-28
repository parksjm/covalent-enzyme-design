"""Collect native-accepted lineage prefixes through the existing parity loader."""

from pathlib import Path
import json
import numpy as np
from nyl12_covalent.relay_pooled_sources import BASE, FIXED, NEW, OLD, load_cohort
from nyl12_covalent.relay_final_sources import MIXING
from nyl12_covalent.relay_umbrella_analysis_io import digest

EXTENDED = BASE / "analysis/main100_continuation_20261127"

SUPPLEMENTS = (
    (OLD, "old_targets4"),
    (NEW, "activation_grid33"),
    (MIXING, "geranium_mixing6"),
)


def accepted_horizon(root: Path, label: str, maximum: int = 30) -> int:
    """Return only consecutive complete, explicitly accepted production blocks."""
    duration = 0
    for number in range(3, maximum // 5 + 3):
        receipt = root / label / f"{number:02d}_production_5ps/analysis_acceptance.json"
        if not receipt.exists():
            break
        data = json.loads(receipt.read_text())
        if data["accepted"] is not True or data["window"] != label:
            raise ValueError("Rejected or wrong-lineage native receipt")
        duration += 5
    return duration


def inventory() -> dict:
    """Select latest published main horizon and every eligible supplement separately."""
    reports = [(h, FIXED / f"fixed_{h}ps_analysis.json") for h in (50, 55, 60)]
    reports += [(h, EXTENDED / f"evidence/fixed_{h}ps_analysis.json") for h in range(65, 101, 5)]
    main, report = max((h, p) for h, p in reports if p.exists())
    lineages = []
    for prepared, group in SUPPLEMENTS:
        data = json.loads(prepared.read_text())
        for window in data["windows"]:
            label = window["label"]
            horizon = accepted_horizon(prepared.parent, label)
            proofs = [
                prepared.parent
                / label
                / f"{i:02d}_production_5ps/analysis_acceptance.json"
                for i in range(3, horizon // 5 + 3)
            ]
            lineages.append(
                dict(
                    group=group,
                    prepared=str(prepared),
                    window=label,
                    horizon_ps=horizon,
                    eligible=horizon >= 10,
                    receipt_sha256={str(p): digest(p) for p in proofs},
                )
            )
    return dict(
        main_horizon_ps=main,
        primary_report=str(report),
        primary_sha256=digest(report),
        supplemental_lineages=lineages,
    )


def collect(selected: dict, work: Path) -> dict:
    """Use unchanged per-lineage native parity checks and one-time 5 ps discard."""
    main = selected["main_horizon_ps"]
    reference = json.loads(
        (
            (FIXED / f"automatic_progression_20261126/checkpoint_{main}ps/analysis_request.json"
             if main <= 60 else EXTENDED / f"checkpoint_{main}ps/analysis_request.json")
        ).read_text()
    )
    groups = [load_cohort(BASE / "prepared.json", main, reference, "main57", 57)]
    for lineage in selected["supplemental_lineages"]:
        if not lineage["eligible"]:
            continue
        original = Path(lineage["prepared"])
        prepared = json.loads(original.read_text())
        window = next(x for x in prepared["windows"] if x["label"] == lineage["window"])
        adapter = work / "source_adapters" / lineage["group"] / lineage["window"]
        adapter.mkdir(parents=True)
        (adapter / window["label"]).symlink_to(
            original.parent / window["label"], target_is_directory=True
        )
        prepared["windows"] = [window]
        prepared["derived_adapter_source"] = dict(
            path=str(original), sha256=digest(original), native_inputs_changed=False
        )
        manifest = adapter / "prepared.json"
        manifest.write_text(json.dumps(prepared, indent=2) + "\n")
        groups.append(
            load_cohort(manifest, lineage["horizon_ps"], reference, lineage["group"], 1)
        )
    if any(x is None for x in groups):
        raise ValueError("Accepted native source disappeared")
    states = [s for group in groups for s in group["states"]]
    proofs = [p for group in groups for p in group["proofs"]]
    csv = [
        next(v for k, v in p["source_sha256"].items() if k.endswith(".csv"))
        for p in proofs
    ]
    if len(set(csv)) != len(csv) or len({s["state_id"] for s in states}) != len(states):
        raise ValueError("Duplicated native block or independent state")
    return dict(
        q=np.concatenate([g["q"] for g in groups]),
        counts=np.concatenate([g["counts"] for g in groups]),
        states=states,
        proofs=proofs,
        config=reference,
        bootstrap_eligible=False,
        inventory=selected,
    )
