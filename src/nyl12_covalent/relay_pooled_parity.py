"""Native per-frame parity contract reused from relay_umbrella_analysis.analyze."""

from __future__ import annotations

import numpy as np

from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies


def validate_parity(block: dict, bias: RelayBiasWindow) -> dict:
    """Validate independent Amber CV and bias-energy parity for every eligible frame.

    Args:
        block: Hash-validated accepted block returned by existing load_block.
        bias: This independent trajectory's declared native umbrella bias.

    Returns:
        Exact frame coverage and maximum errors/bounds; no thinning or fit.
    """
    values = np.asarray(block["values"], float)
    native = np.asarray(block["native_bias_kcal_mol"], float)
    features = dict(zip(("q_acyl", "q_act", "q_del"), values[:, :3].T, strict=True))
    energy = cross_bias_energies(features, [bias])[0]
    mask = np.asarray(block["native_bias_evaluated"], bool)
    if mask.shape != energy.shape or not mask[:-1].all() or not mask.any():
        raise ValueError("Nativeparitycanbeunevaluatedonlyforterminalframe")
    eps = float(block["native_cv_quantization_A"]) + float(
        block["coordinate_serialization_tolerance_A"]
    )
    if not np.isfinite(eps) or eps < 0:
        raise ValueError("InvalidnativeCVbound")
    allowed = np.zeros_like(energy)
    for name, center, force in zip(
        bias.feature_names, bias.centers_A, bias.force_constants, strict=True
    ):
        column = features[name]
        allowed += (
            2
            * force
            * np.minimum(np.abs(column - center) + eps, bias.outer_extent_A)
            * eps
        )
    rounding = float(block["native_energy_quantization_kcal_mol"])
    if not np.isfinite(rounding) or rounding < 0:
        raise ValueError("Invalidnativeenergybound")
    allowed += rounding
    numerical = 32 * np.finfo(float).eps * np.maximum(1, np.abs(energy))
    if (
        native.shape != energy.shape
        or not np.isfinite(native[mask]).all()
        or np.any(abs(native[mask] - energy[mask]) > allowed[mask] + numerical[mask])
    ):
        raise ValueError("NativeAmber/computedper-framebiasparityfailed")
    actual = np.column_stack([features[name] for name in bias.feature_names])
    emitted = np.asarray(block["native_cv_values_A"], float)
    if (
        emitted.shape != actual.shape
        or not np.isfinite(emitted[mask]).all()
        or np.any(abs(emitted[mask] - actual[mask]) > eps + 1e-12)
    ):
        raise ValueError("IndependentnativeCVparityfailed")
    return dict(
        evaluated_frames=int(mask.sum()),
        unevaluated_terminal_frames=int((~mask).sum()),
        native_bias_source=block["native_bias_source"],
        maximum_energy_error_kcal_mol=float(np.max(abs(native[mask] - energy[mask]))),
        maximum_energy_bound_kcal_mol=float(np.max(allowed[mask])),
        CV_bound_A=eps,
    )
