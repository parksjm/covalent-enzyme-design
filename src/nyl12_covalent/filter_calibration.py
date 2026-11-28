"""Retrospectively compare geometry filters with measured activity labels."""

from __future__ import annotations

import logging
import math
from statistics import median
from typing import Any

LOGGER = logging.getLogger(__name__)
GATES = ("d1", "d2", "d3", "a1", "abs_phi")


def _number(value: Any, field: str, *, nonnegative: bool = False) -> float:
    """Return a finite number, rejecting invalid scientific inputs."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be numeric")
    result = float(value)
    if not math.isfinite(result) or (nonnegative and result < 0):
        raise ValueError(f"{field} must be finite and nonnegative when required")
    return result


def _validate_config(config: dict[str, Any]) -> tuple[list[str], int, dict[str, float]]:
    """Validate the cohort definition and inclusive gate thresholds."""
    required = (
        "id",
        "system",
        "register",
        "experimental_comparison",
        "variants",
        "expected_models",
        "thresholds",
    )
    missing = [key for key in required if key not in config]
    if missing:
        raise ValueError(f"Missing config fields: {missing}")
    variants = config["variants"]
    if (
        not isinstance(variants, list)
        or not variants
        or len(set(variants)) != len(variants)
        or "WT" not in variants
    ):
        raise ValueError("variants must be a nonempty unique list")
    n_models = config["expected_models"]
    if isinstance(n_models, bool) or not isinstance(n_models, int) or n_models < 1:
        raise ValueError("expected_models must be a positive integer")
    thresholds = config["thresholds"]
    keys = (
        "d1_max",
        "oxyanion_max",
        "a1_min",
        "a1_max",
        "herschlag_phi_min",
        "herschlag_phi_max",
    )
    if not isinstance(thresholds, dict) or any(k not in thresholds for k in keys):
        raise ValueError("thresholds lack required gate values")
    values = {key: _number(thresholds[key], f"thresholds.{key}") for key in keys}
    if (
        values["a1_min"] > values["a1_max"]
        or values["herschlag_phi_min"] > values["herschlag_phi_max"]
    ):
        raise ValueError("threshold bounds are reversed")
    if (
        values["d1_max"] <= 0
        or values["oxyanion_max"] <= 0
        or not 0 <= values["a1_min"] <= values["a1_max"] <= 180
        or not 0 <= values["herschlag_phi_min"] <= values["herschlag_phi_max"] <= 180
    ):
        raise ValueError(
            "distance maxima must be positive and angles must be within 0..180"
        )
    return variants, n_models, values


def _pearson(rows: list[dict[str, Any]]) -> dict[str, float | int] | None:
    """Compute descriptive Pearson correlations across non-WT variants."""
    pairs = [
        (row["tier_ab_fraction"], row["selected_endpoint_relative_activity"])
        for row in rows
        if row["variant"] != "WT"
        and row["selected_endpoint_relative_activity"] is not None
    ]
    if len(pairs) < 3:
        return None
    xs, ys = zip(*pairs)
    xbar, ybar = sum(xs) / len(xs), sum(ys) / len(ys)
    xss = sum((x - xbar) ** 2 for x in xs)
    yss = sum((y - ybar) ** 2 for y in ys)
    if xss == 0 or yss == 0:
        return None
    correlation = sum((x - xbar) * (y - ybar) for x, y in pairs) / math.sqrt(xss * yss)
    return {"n_variants": len(pairs), "pearson_r": correlation}


def analyze_cohort(
    records: list[dict[str, Any]],
    activities: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    """Analyze a fixed structural cohort against explicitly selected activities.

    WT endpoint activity is normalized to 1.0 by convention. It is not treated
    as a raw replicate. Other endpoints retain their measured fold values;
    positive reports in any same-system comparison are descriptive evidence
    only and are never pooled into a total-activity estimate.

    Args:
        records: Normalized per-variant model rows with measured geometry.
        activities: Original activity rows, including system and comparison.
        config: Cohort identity, variants, model count, and gate thresholds.

    Returns:
        A JSON-compatible cohort summary with per-variant and per-model rows.

    Raises:
        ValueError: If cohort rows, activity labels, or thresholds are invalid.
    """
    variants, expected_models, t = _validate_config(config)
    system, comparison = config["system"], config["experimental_comparison"]
    indexed: dict[str, dict[int, dict[str, Any]]] = {
        variant: {} for variant in variants
    }
    per_model: list[dict[str, Any]] = []
    mismatches: list[dict[str, Any]] = []
    for source in records:
        variant, model = source.get("variant"), source.get("model_index")
        if variant not in indexed:
            raise ValueError(f"Unexpected variant: {variant}")
        if (
            isinstance(model, bool)
            or not isinstance(model, int)
            or not 0 <= model < expected_models
        ):
            raise ValueError(f"Invalid model_index for {variant}: {model}")
        if model in indexed[variant]:
            raise ValueError(f"Duplicate model {model} for {variant}")
        tier = source.get("geometry_tier")
        if tier not in {"A", "B", "C", "reject"}:
            raise ValueError(f"Invalid geometry_tier for {variant}/{model}: {tier}")
        vals = {
            key: _number(
                source.get(key),
                f"{variant}/{model}.{key}",
                nonnegative=key in ("d1", "d2", "d3"),
            )
            for key in ("d1", "d2", "d3", "a1", "herschlag_phi", "interface_delta")
        }
        gates = {
            "d1": vals["d1"] <= t["d1_max"],
            "d2": vals["d2"] <= t["oxyanion_max"],
            "d3": vals["d3"] <= t["oxyanion_max"],
            "a1": t["a1_min"] <= vals["a1"] <= t["a1_max"],
            "abs_phi": t["herschlag_phi_min"]
            <= abs(vals["herschlag_phi"])
            <= t["herschlag_phi_max"],
        }
        gate_pass = all(gates.values())
        stored_retained = tier in {"A", "B"}
        mismatch = gate_pass != stored_retained
        row = {
            "variant": variant,
            "model_index": model,
            "geometry_tier": tier,
            "gate_flags": gates,
            "gate_pass": gate_pass,
            "tier_mismatch": mismatch,
            "d1": vals["d1"],
            "d2": vals["d2"],
            "d3": vals["d3"],
            "a1": vals["a1"],
            "herschlag_phi": vals["herschlag_phi"],
            "interface_delta": vals["interface_delta"],
        }
        indexed[variant][model] = row
        per_model.append(row)
        if mismatch:
            mismatches.append(
                {
                    "variant": variant,
                    "model_index": model,
                    "stored_tier": tier,
                    "gate_pass": gate_pass,
                }
            )
    for variant in variants:
        actual = set(indexed[variant])
        expected = set(range(expected_models))
        if actual != expected:
            raise ValueError(
                f"Incomplete models for {variant}: missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
            )

    selected: dict[str, float] = {}
    any_positive = {variant: False for variant in variants}
    for activity in activities:
        if activity.get("system") != system:
            continue
        variant = activity.get("variant")
        if variant not in any_positive:
            continue
        value = _number(
            activity.get("relative_activity"),
            f"activity {variant}/{activity.get('comparison')}",
            nonnegative=True,
        )
        if value > 0:
            any_positive[variant] = True
        if activity.get("comparison") == comparison:
            if variant in selected:
                raise ValueError(f"Duplicate selected activity mapping for {variant}")
            selected[variant] = value

    per_variant: list[dict[str, Any]] = []
    for variant in variants:
        models = [indexed[variant][i] for i in range(expected_models)]
        retained = [row for row in models if row["geometry_tier"] in {"A", "B"}]
        selected_value = 1.0 if variant == "WT" else selected.get(variant)
        label = (
            "wt_normalized_reference"
            if variant == "WT"
            else "missing"
            if selected_value is None
            else "zero_reported_not_confirmed_inactive"
            if selected_value == 0
            else "reported_positive_selected_endpoint"
        )
        gate_fail_counts = {
            gate: sum(not row["gate_flags"][gate] for row in models) for gate in GATES
        }
        gate_only = {
            gate: sum(
                not row["gate_flags"][gate]
                and all(row["gate_flags"][other] for other in GATES if other != gate)
                for row in models
            )
            for gate in GATES
        }
        phis = [row["herschlag_phi"] for row in models]
        per_variant.append(
            {
                "variant": variant,
                "n_models": expected_models,
                "n_tier_a": sum(row["geometry_tier"] == "A" for row in models),
                "n_tier_ab": len(retained),
                "tier_ab_fraction": len(retained) / expected_models,
                "any_model_retained": bool(retained),
                "gate_fail_counts": gate_fail_counts,
                "gate_only_failure_counts": gate_only,
                "n_phi_positive": sum(value > 0 for value in phis),
                "n_phi_negative": sum(value < 0 for value in phis),
                "n_phi_zero": sum(value == 0 for value in phis),
                "selected_endpoint_relative_activity": selected_value,
                "endpoint_label": label,
                "any_product_reported_positive": any_positive[variant],
                "median_interface_delta_tier_ab": median(
                    [row["interface_delta"] for row in retained]
                )
                if retained
                else None,
            }
        )
    non_wt = [row for row in per_variant if row["variant"] != "WT"]
    known_positive_retained = sum(
        row["any_model_retained"] and row["any_product_reported_positive"]
        for row in non_wt
    )
    known_positive_count = sum(row["any_product_reported_positive"] for row in non_wt)
    summary = {
        "n_variants": len(variants),
        "known_positive_variant_count": known_positive_count,
        "known_positive_variants_with_retained_model": known_positive_retained,
        "known_positive_variant_count_basis": "variants, excluding WT; positivity in any same-system product comparison",
        "zero_activity_caveat": "A reported zero is not evidence of confirmed inactivity.",
        "predictor_claim": "No validated predictor claim; retrospective descriptive calibration only.",
        "missing_endpoint_labels": [
            row["variant"] for row in per_variant if row["endpoint_label"] == "missing"
        ],
        "descriptive_activity_pearson": _pearson(per_variant),
    }
    LOGGER.info(
        "Analyzed %d variants and %d models; %d stored-tier/gate mismatches",
        len(variants),
        len(per_model),
        len(mismatches),
    )
    return {
        "id": config["id"],
        "system": system,
        "register": config["register"],
        "experimental_comparison": comparison,
        "per_variant": per_variant,
        "per_model": per_model,
        "summary": summary,
        "tier_ab_mismatches": mismatches,
    }
