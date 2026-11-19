"""Mapped physical gates for a diagnostic composite-coordinate path."""

from __future__ import annotations


def chemistry_retained(values: dict[str, float]) -> bool:
    """Check identities while allowing shared HG1 and asynchronous progress."""
    return (
        all(0.85 <= values[f"N_{h}_A"] <= 1.25 for h in ("H1", "H2"))
        and all(values[f"OG1_{h}_A"] >= 1.30 for h in ("H1", "H2"))
        and 0.75 <= min(values["N_HG1_A"], values["OG1_HG1_A"]) <= 1.70
        and 1.30 <= values["OG1_carbonyl_C_A"] <= 4.00
        and 1.10 <= values["carbonyl_C_O_A"] <= 1.55
        and 1.20 <= values["carbonyl_C_N_A"] <= 1.80
    )


def candidate_geometry(values: dict[str, float]) -> bool:
    """Identify a mapped TI1-like candidate for an unbiased release test."""
    return (
        chemistry_retained(values)
        and values["OG1_carbonyl_C_A"] <= 1.80
        and values["carbonyl_C_O_A"] >= 1.28
        and 0.85 <= values["N_HG1_A"] <= 1.25
        and values["OG1_HG1_A"] >= 1.30
    )
