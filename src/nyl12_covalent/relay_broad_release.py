"""Explicit four-template allocation and fail-closed broad-pilot release guards."""

from __future__ import annotations

from typing import Any

from .relay_broad_pilot import symmetry_manifest

FOUR_TEMPLATE_PROTOCOL = "two_original_two_remodeled_six_each_20261127"


def four_template_schedule(
    original: list[str], remodeled: list[str]
) -> list[dict[str, Any]]:
    """Allocate exactly six samples to each of two original and two remodeled TSs.

    Args:
        original: Two distinct original TS source IDs.
        remodeled: Two independently accepted remodeled output IDs.

    Returns:
        Four explicit allocations totaling 24 raw samples.

    Raises:
        ValueError: Parent counts or global identities differ from the protocol.
    """
    if len(original) != 2 or len(remodeled) != 2 or len(set(original + remodeled)) != 4:
        raise ValueError(
            "Four-template protocol requires two distinct parents per kind"
        )
    return [
        {"parent_id": parent, "parent_kind": kind, "raw_samples": 6}
        for names, kind in ((original, "original_TS"), (remodeled, "remodeled_TS"))
        for parent in names
    ]


def validate_generation_request(
    request: dict[str, Any], protection: dict[str, Any]
) -> None:
    """Require the reviewed 229/99 mask, exact ties and fixed 24-sample settings.

    Args:
        request: Concrete four-template generation request.
        protection: Existing canonical sequence-protection census.

    Raises:
        ValueError: Budget, protection mask, ties or sampling method differs.
    """
    generation = request["generation"]
    mask = generation["mask"]
    if mask != {
        "alpha": protection["mutable"]["alpha"],
        "beta": protection["mutable"]["beta"],
    }:
        raise ValueError("Broad mutable mask differs from canonical protection receipt")
    if (len(mask["alpha"]), len(mask["beta"]), protection["counts"]["fixed_total"]) != (
        172,
        57,
        99,
    ):
        raise ValueError("Expected exactly 229 mutable and 99 protected positions")
    expected = symmetry_manifest(mask)
    if any(generation[key] != value for key, value in expected.items()):
        raise ValueError("Four-copy symmetry or redesigned positions differ")
    if (
        generation["temperature"],
        generation["native_identity_bias"],
        generation["max_raw_samples"],
    ) != (0.1, 0, 24):
        raise ValueError("Sampling settings differ from authorized method")


def validate_original_lineage(parents: list[dict[str, Any]]) -> None:
    """Require independently accepted remodeled outputs from both original TSs.

    Args:
        parents: Exact-hash accepted parent records, checked separately on disk.

    Raises:
        ValueError: Either original is absent, duplicated or uses a rejected arm.
    """
    if len(parents) != 2 or {p["parent_source_id"] for p in parents} != {
        "medoid01",
        "medoid03",
    }:
        raise ValueError("Both distinct original TS sources must qualify")
    if any(p["noise_A"] != 5 or p["sequence_scope"] != "native_WT" for p in parents):
        raise ValueError("This release permits only accepted native-WT 5A candidates")


def validate_idle_gpu(gpu_rows: str, compute_rows: str) -> None:
    """Require idle Violet GPU0 immediately before bounded generation.

    Args:
        gpu_rows: nvidia-smi GPU0 memory.used and utilization.gpu CSV output.
        compute_rows: nvidia-smi GPU0 compute PID CSV output.

    Raises:
        RuntimeError: GPU0 has a compute process, material memory use or activity.
    """
    memory, utilization = [int(value.strip()) for value in gpu_rows.strip().split(",")]
    if compute_rows.strip() or memory > 256 or utilization != 0:
        raise RuntimeError(
            "Violet GPU0 is occupied; generation remains held without retry"
        )
