"""Build and independently audit selected cohort records."""

from __future__ import annotations
from collections import defaultdict
from typing import Any
from nyl12_covalent.design_boltz_inputs import _mutations
from nyl12_covalent.design_mask_validation_inputs import _candidate_mutations


def _audit_samples_and_build_records(
    samples: list[dict[str, Any]],
    selected: list[dict[str, Any]],
    masks: dict[str, dict[str, set[int]]],
    protected: dict[str, set[int]],
    native_alpha: str,
    native_beta: str,
    source_strata: set[tuple[str, str, float]],
) -> tuple[list[dict[str, Any]], dict[str, set[int]]]:
    """Audit all source samples and build WT, controls, and selected records."""
    sample_counts: dict[tuple[str, str, float], int] = defaultdict(int)
    for sample in samples:
        key = (sample["arm"], str(sample["parent"]), float(sample["temperature"]))
        if key not in source_strata:
            raise ValueError(
                f"Sample {sample['sample_id']} is absent from the frozen generation manifest"
            )
        sample_counts[key] += 1
        expected_mutations = _mutations(
            sample["alpha_sequence"],
            sample["beta_sequence"],
            native_alpha,
            native_beta,
            masks[sample["arm"]]["alpha"],
            masks[sample["arm"]]["beta"],
            protected["alpha"],
            protected["beta"],
        )
        if expected_mutations != sample.get("mutations") or len(
            expected_mutations
        ) != sample.get("mutation_count"):
            raise ValueError(
                f"Sample {sample['sample_id']} failed an independent sequence audit"
            )
    if set(sample_counts) != source_strata or any(
        count != 6 for count in sample_counts.values()
    ):
        raise ValueError(
            "Audited report must retain six samples for every frozen stratum"
        )

    records = [
        {
            "id": "WT",
            "alpha_sequence": native_alpha,
            "beta_sequence": native_beta,
            "arms": [],
            "selection_lineage": [],
            "cohort_tag": "baseline",
        }
    ]
    for control, pos, target in (
        ("control_K59D", 59, "D"),
        ("control_V111C", 111, "C"),
    ):
        if native_alpha[pos - 1] != ("K" if pos == 59 else "V"):
            raise ValueError(
                f"Native alpha sequence does not support {control} at position {pos}"
            )
        records.append(
            {
                "id": control,
                "alpha_sequence": native_alpha[: pos - 1] + target + native_alpha[pos:],
                "beta_sequence": native_beta,
                "arms": [],
                "selection_lineage": [],
                "cohort_tag": "matched_positive_control",
            }
        )
    sample_by_id = {sample["sample_id"]: sample for sample in samples}
    union = {
        family: set().union(*(mask[family] for mask in masks.values()))
        for family in ("alpha", "beta")
    }
    for index, item in enumerate(selected, 1):
        alpha, beta = item["alpha_sequence"], item["beta_sequence"]
        mutations = _candidate_mutations(alpha, beta, native_alpha, native_beta)
        for lineage in item["selection_lineage"]:
            sample = sample_by_id[lineage["sample_id"]]
            if (sample["alpha_sequence"], sample["beta_sequence"]) != (
                alpha,
                beta,
            ) or sample["arm"] != lineage["arm"]:
                raise ValueError("Selected lineage does not match its source sample")
            for mutation in mutations:
                if (
                    mutation["physical_position"]
                    not in masks[lineage["arm"]][mutation["family"]]
                ):
                    raise ValueError(
                        f"Mutation falls outside source arm mask: {mutation}"
                    )
        checked = _mutations(
            alpha,
            beta,
            native_alpha,
            native_beta,
            union["alpha"],
            union["beta"],
            protected["alpha"],
            protected["beta"],
        )
        if checked != mutations:
            raise ValueError("Candidate mutation audit disagrees with source sequence")
        records.append(
            {
                "id": f"mask_{index:02d}",
                "alpha_sequence": alpha,
                "beta_sequence": beta,
                "arms": item["arms"],
                "selection_lineage": item["selection_lineage"],
                "mutations": mutations,
                "mutation_count": len(mutations),
                "cohort_tag": "masked_design",
            }
        )
    for row in records[:3]:
        row.setdefault(
            "mutations",
            _candidate_mutations(
                row["alpha_sequence"], row["beta_sequence"], native_alpha, native_beta
            ),
        )
        row.setdefault("mutation_count", len(row["mutations"]))
    return records, union
