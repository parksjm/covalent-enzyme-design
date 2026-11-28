"""Validate the sampling horizon represented by an umbrella movie export."""

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from typing import Any


def validated_sampling_scope(
    campaign: Path, prepared: dict[str, Any], manifest: dict[str, Any]
) -> dict[str, Any]:
    """Derive movie coverage from its actual accepted production blocks.

    Args:
        campaign: Root containing the prepared windows and accepted receipts.
        prepared: Authoritative campaign preparation manifest.
        manifest: Export manifest recording each sampled frame and block.

    Returns:
        Equal per-window horizon, frame count, block count and receipt hashes.

    Raises:
        ValueError: If coverage is inconsistent, nonconsecutive or unaccepted.
    """
    labels = [row["label"] for row in prepared["windows"]]
    if manifest["windows"] != labels:
        raise ValueError("Movie window order differs from the prepared campaign")
    records = manifest["frame_plan"]
    counts = Counter(record["window"] for record in records)
    if set(counts) != set(labels) or len(set(counts.values())) != 1:
        raise ValueError("Movie must cover every window with equal frame counts")
    if len(records) != manifest["frames"]:
        raise ValueError("Movie frame plan length does not match its frame count")
    blocks: dict[str, set[str]] = defaultdict(set)
    for record in records:
        blocks[record["window"]].add(record["block"])
    horizons, block_counts, hashes = [], [], {}
    for row in prepared["windows"]:
        label = row["label"]
        stages = [s["label"] for s in row["stages"] if s["kind"] == "production"]
        selected = blocks[label]
        if selected != set(stages[: len(selected)]):
            raise ValueError(f"Movie must use a consecutive production prefix: {label}")
        duration = 0.0
        for block in sorted(selected):
            path = campaign / label / block / "analysis_acceptance.json"
            payload = path.read_bytes()
            receipt = json.loads(payload)
            if (
                receipt["accepted"] is not True
                or receipt["window"] != label
                or receipt["physical_hamiltonian_sha256"]
                != prepared["physical_hamiltonian_sha256"]
                or receipt["sampled_duration_ps"] <= 0
            ):
                raise ValueError(f"Invalid production receipt for movie: {path}")
            duration += float(receipt["sampled_duration_ps"])
            hashes[str(path.resolve())] = hashlib.sha256(payload).hexdigest()
        horizons.append(duration)
        block_counts.append(len(selected))
    if len(set(horizons)) != 1 or len(set(block_counts)) != 1:
        raise ValueError("Movie windows have different sampled horizons")
    return {
        "sampled_horizon_per_window_ps": horizons[0],
        "sampled_frames_per_window": next(iter(counts.values())),
        "production_blocks_per_window": block_counts[0],
        "acceptance_receipt_sha256": hashes,
    }
