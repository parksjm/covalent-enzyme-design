"""Preserve raw sampling lineage while separating replacement pilot identities."""

from __future__ import annotations

from collections.abc import Callable
import re
from typing import Any


def namespace_pool(
    samples: list[dict[str, Any]],
    *,
    pool: Callable[[list[dict[str, Any]]], list[dict[str, Any]]],
    namespace: str,
) -> list[dict[str, Any]]:
    """Prefix raw identifiers before deduplication and rename pooled candidates.

    Args:
        samples: Audited native samples, modified to retain source and new IDs.
        pool: Existing whole-sequence deduplicator, unchanged computational logic.
        namespace: Explicit alphanumeric replacement campaign identifier.

    Returns:
        Deduplicated candidates with coherent sample IDs and complete lineage.

    Raises:
        ValueError: Namespace, raw identities, or pooled lineage is inconsistent.
    """
    if re.fullmatch(r"[a-z][a-z0-9]*", namespace) is None:
        raise ValueError("Require a lowercase alphanumeric namespace")
    identifiers = [sample["sample_id"] for sample in samples]
    if len(set(identifiers)) != len(identifiers) or any(
        "source_sample_id" in sample for sample in samples
    ):
        raise ValueError("Raw sample identities must be unique and unnamespaced")
    for sample in samples:
        sample["source_sample_id"] = sample["sample_id"]
        sample["sample_id"] = f"{namespace}_{sample['sample_id']}"
    candidates = pool(samples)
    observed = []
    for index, candidate in enumerate(candidates, 1):
        candidate["source_candidate_id"] = candidate["candidate_id"]
        candidate["candidate_id"] = f"relay_{namespace}_{index:02d}"
        if candidate["sample_id"] not in candidate["lineage"]:
            raise ValueError("Candidate representative is absent from its lineage")
        observed.extend(candidate["lineage"])
    if sorted(observed) != sorted(sample["sample_id"] for sample in samples):
        raise ValueError("Pooling changed or lost raw sampling lineage")
    return candidates
