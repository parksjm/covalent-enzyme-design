"""Cumulative production-block selection and common-reference PMF overlays."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

PRODUCTION_STAGE = re.compile(r"^(?P<index>\d+)_production_5ps$")


def select_horizon(
    receipts: list[Path], windows: list[dict[str, Any]], horizon_ps: float
) -> tuple[list[Path], dict[str, float]] | None:
    """Select complete accepted receipts ending at one shared duration.

    Args:
        receipts: Sorted immutable accepted production-block receipts.
        windows: Ordered analysis manifest windows.
        horizon_ps: Requested cumulative duration per window.

    Returns:
        Selected receipts and per-window durations, or ``None`` until all
        windows have data through the exact common horizon.
    """
    labels = {window["label"] for window in windows}
    selected: list[Path] = []
    durations = dict.fromkeys(labels, 0.0)
    block_indices: dict[str, list[int]] = {label: [] for label in labels}
    for path in receipts:
        proof = json.loads(path.read_text())
        label = proof["window"]
        if label not in durations:
            raise ValueError(f"Unexpected accepted receipt window: {label}")
        duration = float(proof["sampled_duration_ps"])
        if proof.get("accepted") is not True or duration <= 0:
            raise ValueError(f"Invalid accepted production receipt: {path}")
        stage = next(
            (
                PRODUCTION_STAGE.fullmatch(part)
                for part in path.parts
                if PRODUCTION_STAGE.fullmatch(part)
            ),
            None,
        )
        if stage is None:
            raise ValueError(
                f"Receipt is not from a named 5 ps production block: {path}"
            )
        index = int(stage.group("index"))
        if index in block_indices[label]:
            raise ValueError(f"Duplicate production block index for {label}: {index}")
        block_indices[label].append(index)
        if durations[label] + duration <= horizon_ps + 1e-8:
            durations[label] += duration
            selected.append(path)
    for label, indices in block_indices.items():
        if indices and indices != list(range(3, 3 + len(indices))):
            raise ValueError(f"Noncontiguous production blocks for {label}: {indices}")
    if any(abs(duration - horizon_ps) > 1e-6 for duration in durations.values()):
        return None
    return selected, durations


def plot_progression(
    reports: list[tuple[float, dict[str, Any]]], output: Path
) -> dict[str, Any]:
    """Plot only successful connected PMFs on their shared declared grid.

    Args:
        reports: Duration/report pairs for successful connected analyses.
        output: New figure path; existing outputs are not overwritten.
    """
    if output.exists():
        raise FileExistsError(output)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    reference: dict[str, Any] | None = None
    plotted: list[tuple[Any, float, dict[str, Any]]] = []
    for duration, report in reports:
        pmf = report.get("pmf")
        if report.get("status") != "preliminary_connected_not_converged" or not pmf:
            continue
        identity = {
            key: report.get(key)
            for key in (
                "physical_hamiltonian_sha256",
                "temperature_K",
                "path_3cv_sha256",
                "projection",
                "es_search_bounds",
                "acyl_search_bounds",
            )
        }
        identity["edges"] = pmf.get("edges")
        if reference is None:
            reference = identity
        elif reference != identity:
            raise ValueError(
                "Progression PMFs do not share physics, projection, and reference bounds"
            )
        values = [np.nan if value is None else value for value in pmf["pmf_kcal_mol"]]
        (line,) = axis.plot(pmf["centers"], values, ".-", label=f"{duration:g} ps")
        plotted.append((line, duration, report))
    if reference is None:
        axis.text(0.5, 0.5, "No connected PMF checkpoints yet", ha="center")
    else:
        common_reference_s = plotted[0][2]["pmf"]["reference_s"]
        included: list[float] = []
        omitted: list[float] = []
        offsets: list[tuple[float, float]] = []
        for line, duration, report in plotted:
            pmf = report.get("pmf")
            if not pmf or report.get("status") != "preliminary_connected_not_converged":
                continue
            index = pmf["centers"].index(common_reference_s)
            reference_value = pmf["pmf_kcal_mol"][index]
            if reference_value is None:
                line.set_label(f"{duration:g} ps (no support at common reference)")
                line.set_visible(False)
                omitted.append(duration)
            else:
                offsets.append((duration, float(reference_value)))
                line.set_ydata(
                    [
                        None if value is None else value - reference_value
                        for value in pmf["pmf_kcal_mol"]
                    ]
                )
                included.append(duration)
        axis.set(
            xlabel="s (raw normalized path arclength)",
            ylabel="F(s), shifted to common E-S bin (kcal/mol)",
        )
        axis.legend(title="Cumulative production")
    if reference is None:
        axis.set_title(
            "Preliminary MBAR progression (not converged; no uncertainty shown)"
        )
    else:
        axis.set_title(
            "Preliminary MBAR progression (not converged; "
            f"common E-S reference s={common_reference_s:.3f}; no uncertainty shown)"
        )
    figure.savefig(output, dpi=160)
    plt.close(figure)
    return {
        "common_reference_s": common_reference_s if reference is not None else None,
        "included_horizons_ps": included if reference is not None else [],
        "omitted_horizons_ps": omitted if reference is not None else [],
        "reference_offsets_kcal_mol": (offsets if reference is not None else []),
        "uncertainty": "not calculated",
    }
