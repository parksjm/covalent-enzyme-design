"""Actual most-extensive100-lineage source gate; no frozen45ps substitution."""

from __future__ import annotations

import json

import numpy as np

from .relay_pooled_sources import (
    BASE,
    FIXED,
    NEW,
    OLD,
    load_cohort,
)

MIXING = (
    BASE
    / "analysis/final_candidate_20261127/geranium_derived/mixing_checks_20261127/prepared.json"
)
COHORTS = (
    (BASE / "prepared.json", 60, "main60", 57),
    (NEW, 30, "grid30", 33),
    (OLD, 30, "old_targets30", 4),
    (MIXING, 30, "mixing_checks30", 6),
)


def readiness() -> dict:
    """Inspect exact acceptance prefixes cheaply before expensive source validation.

    Returns:
        Every cohort's actual minimum/maximum complete5psproduction prefix.
    """
    cohorts = []
    for path, horizon, group, count in COHORTS:
        prepared = json.loads(path.read_text())
        if len(prepared["windows"]) != count:
            raise ValueError("Finalcohortinventorychanged:" + group)
        durations = []
        for window in prepared["windows"]:
            duration = 0
            for block in range(3, horizon // 5 + 3):
                receipt = (
                    path.parent
                    / window["label"]
                    / f"{block:02d}_production_5ps/analysis_acceptance.json"
                )
                if not receipt.is_file():
                    break
                proof = json.loads(receipt.read_text())
                if proof["accepted"] is not True or proof["window"] != window["label"]:
                    raise ValueError("Nonacceptedorwronglineagereceipt")
                duration += 5
            durations.append(duration)
        cohorts.append(
            dict(
                group=group,
                required_production_ps=horizon,
                lineages=count,
                minimum_accepted_ps=min(durations),
                maximum_accepted_ps=max(durations),
                complete=all(x == horizon for x in durations),
            )
        )
    return dict(
        ready=all(c["complete"] for c in cohorts),
        cohorts=cohorts,
        required_lineages=100,
        expected_postdiscard_frames=84200,
    )


def final_data() -> dict | None:
    """Load all actual60psmain and30pssupplementary trajectories after individual5psdiscards.

    Returns:
       100separateMBARstates, allphysicalfeatures/proofs, orNone until complete.
    """
    if not readiness()["ready"]:
        return None
    reference = json.loads(
        (
            FIXED
            / "automatic_progression_20261126/checkpoint_45ps/analysis_request.json"
        ).read_text()
    )
    groups = []
    for path, horizon, group, count in COHORTS:
        cohort = load_cohort(path, horizon, reference, group, count)
        if cohort is None:
            raise ValueError("Acceptedfinalprefixdisappeared")
        groups.append(cohort)
    states = [s for g in groups for s in g["states"]]
    counts = np.concatenate([g["counts"] for g in groups])
    if len(states) != 100 or counts.tolist() != [1100] * 57 + [500] * 43:
        raise ValueError("Finalanalysisrequiresactualmain60plusall43supplement30")
    if len({s["state_id"] for s in states}) != 100:
        raise ValueError("Independentlineagesmerged")
    proofs = [p for g in groups for p in g["proofs"]]
    csv_hashes = [
        next(v for k, v in p["source_sha256"].items() if k.endswith(".csv"))
        for p in proofs
    ]
    if len(set(csv_hashes)) != len(csv_hashes):
        raise ValueError("Duplicatedphysicalsourceblock")
    return dict(
        q=np.concatenate([g["q"] for g in groups]),
        counts=counts,
        states=states,
        proofs=proofs,
        samples=[s for g in groups for s in g["samples"]],
        config=reference,
        main_actual_horizon_ps=60,
        supplemental_horizon_ps=30,
        bootstrap_eligible=True,
        prepared_sources=[
            dict(path=g["prepared_path"], sha256=g["prepared_sha256"]) for g in groups
        ],
        comparison_scope="Mostextensive100lineageFINALCANDIDATE;notautomaticallyconvergedorprimary",
    )
