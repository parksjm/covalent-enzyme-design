"""Deterministic selection of diverse, audited mask candidates."""

from __future__ import annotations
from collections import defaultdict
import math
from typing import Any

AA = set("ACDEFGHIKLMNPQRSTVWY")


def _sequence_key(sample: dict[str, Any]) -> tuple[str, str]:
    return sample["alpha_sequence"], sample["beta_sequence"]


def _distance(left: tuple[str, str], right: tuple[str, str]) -> int:
    return sum(
        a != b
        for x, y in zip(left, right, strict=True)
        for a, b in zip(x, y, strict=True)
    )


def select_candidates(samples: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Select up to four diverse samples per arm with complete lineage.

    The first parent/temperature stratum chooses its lexicographically smallest
    sample ID. Later strata greedily maximize minimum sequence Hamming distance
    from earlier selections within that arm. Identical sequences selected in
    both arms become one record with both arm lineages.

    Args:
        samples: Audited records containing sample, arm, parent, temperature,
            alpha_sequence, and beta_sequence fields.

    Returns:
        Sequence-deduplicated candidates with source arms and selection lineage.

    Raises:
        ValueError: A sample field is missing/invalid or the two-by-two-by-two
            arm/parent/temperature design is incomplete.
    """
    required = {
        "sample_id",
        "arm",
        "parent",
        "temperature",
        "alpha_sequence",
        "beta_sequence",
    }
    if not samples:
        raise ValueError("At least one audited sample is required")
    ids: set[str] = set()
    arms, parents, temperatures = set(), set(), set()
    lengths: tuple[int, int] | None = None
    strata: dict[tuple[str, str, float], list[dict[str, Any]]] = defaultdict(list)
    for sample in samples:
        missing = required - sample.keys()
        if missing:
            raise ValueError(f"Sample missing fields: {sorted(missing)}")
        sample_id, arm, parent = (
            sample[name] for name in ("sample_id", "arm", "parent")
        )
        alpha, beta = sample["alpha_sequence"], sample["beta_sequence"]
        if not all(
            isinstance(value, str) and value
            for value in (sample_id, arm, str(parent), alpha, beta)
        ):
            raise ValueError(
                "Sample identifiers, parent/arm and sequences must be nonempty"
            )
        if sample_id in ids:
            raise ValueError(f"Duplicate sample_id: {sample_id}")
        ids.add(sample_id)
        if not set(alpha) <= AA or not set(beta) <= AA:
            raise ValueError(
                f"Sample {sample_id} contains an invalid amino-acid symbol"
            )
        if lengths is None:
            lengths = (len(alpha), len(beta))
        if lengths != (len(alpha), len(beta)):
            raise ValueError(f"Sample {sample_id} has inconsistent chain lengths")
        try:
            temperature = float(sample["temperature"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Sample {sample_id} has invalid temperature") from exc
        if not math.isfinite(temperature) or temperature <= 0:
            raise ValueError(f"Sample {sample_id} has invalid temperature")
        arms.add(arm)
        parents.add(str(parent))
        temperatures.add(temperature)
        strata[(arm, str(parent), temperature)].append(sample)
    if len(arms) != 2 or len(parents) != 2 or len(temperatures) != 2:
        raise ValueError(
            "Samples must cover exactly two arms, parents, and temperatures"
        )
    expected = {
        (arm, parent, temp)
        for arm in arms
        for parent in parents
        for temp in temperatures
    }
    if set(strata) != expected:
        raise ValueError("Samples do not cover every arm/parent/temperature stratum")

    selected: dict[tuple[str, str], dict[str, Any]] = {}
    for arm in sorted(arms):
        chosen: list[tuple[str, str]] = []
        for _, parent, temperature in sorted(
            (key for key in strata if key[0] == arm), key=lambda x: (x[1], x[2])
        ):
            choices = sorted(
                strata[(arm, parent, temperature)], key=lambda row: row["sample_id"]
            )
            unique = {_sequence_key(row): row for row in reversed(choices)}
            available = [(key, row) for key, row in unique.items() if key not in chosen]
            if not available:
                continue
            if not chosen:
                key, row = min(available, key=lambda item: item[1]["sample_id"])
            else:
                key, row = min(
                    available,
                    key=lambda item: (
                        -min(_distance(item[0], previous) for previous in chosen),
                        item[1]["sample_id"],
                    ),
                )
            chosen.append(key)
            lineage = {
                "sample_id": row["sample_id"],
                "parent": str(row["parent"]),
                "temperature": row["temperature"],
                "arm": arm,
            }
            if key not in selected:
                selected[key] = {
                    "alpha_sequence": key[0],
                    "beta_sequence": key[1],
                    "arms": [],
                    "selection_lineage": [],
                }
            selected[key]["arms"].append(arm)
            selected[key]["selection_lineage"].append(lineage)
    return sorted(
        selected.values(), key=lambda row: (row["alpha_sequence"], row["beta_sequence"])
    )
