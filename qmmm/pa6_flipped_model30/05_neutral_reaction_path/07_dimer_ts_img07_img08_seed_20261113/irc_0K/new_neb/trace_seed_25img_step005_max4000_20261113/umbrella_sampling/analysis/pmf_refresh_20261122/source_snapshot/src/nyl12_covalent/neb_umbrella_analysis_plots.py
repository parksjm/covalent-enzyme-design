"""Matplotlib scientific plots of immutable umbrella-analysis JSON artifacts."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def export(output: Path, result: dict[str, Any]) -> None:
    """Write PNG/SVG coverage, overlap and interim FES using standard Matplotlib."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def save(fig: Any, name: str) -> None:
        fig.suptitle(
            "INTERIM umbrella diagnostics; no barrier or convergence acceptance",
            fontsize=10,
        )
        fig.tight_layout()
        for extension in ("png", "svg"):
            fig.savefig(output / f"{name}.{extension}", dpi=160)
        plt.close(fig)

    if "coverage_points" not in result:
        fig, axis = plt.subplots(figsize=(7, 4))
        axis.axis("off")
        axis.text(0.5, 0.5, "Waiting for completed production blocks", ha="center")
        save(fig, "waiting_for_production")
        return
    points = np.asarray(result["coverage_points"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for axis, y, label in [
        (axes[0], 1, "qH (Å)"),
        (axes[1], 2, "Proton PI (coordination difference)"),
    ]:
        axis.scatter(points[:, 0], points[:, y], s=2, alpha=0.2, rasterized=True)
        axis.set(xlabel="qAcyl (Å)", ylabel=label, title="All valid production frames")
    save(fig, "sample_coverage")
    for key in ("histogram_overlap", "mbar_overlap"):
        if key not in result:
            continue
        fig, axis = plt.subplots(figsize=(7, 6))
        image = axis.imshow(
            np.asarray(result[key]), origin="lower", vmin=0, vmax=1, cmap="viridis"
        )
        axis.set(
            xlabel="Window index", ylabel="Window index", title=key.replace("_", " ")
        )
        fig.colorbar(image, ax=axis, label="Overlap")
        save(fig, key)
    for key, ylabel in [
        ("primary_fes", "Proton PI (coordination difference)"),
        ("diagnostic_fes", "qH (Å)"),
    ]:
        if key not in result:
            continue
        record = result[key]
        values = np.array(
            [
                [np.nan if v is None else v for v in row]
                for row in record["free_energy_kcal_mol"]
            ]
        )
        fig, axis = plt.subplots(figsize=(7, 6))
        axis.set_facecolor("#eeeeee")
        image = axis.pcolormesh(
            record["x_edges"],
            record["y_edges"],
            np.ma.masked_invalid(values.T),
            shading="flat",
            cmap="viridis",
        )
        axis.set(
            xlabel="qAcyl (Å)",
            ylabel=ylabel,
            title="Interim unbiased histogram FES; unsampled bins gray",
        )
        fig.colorbar(image, ax=axis, label="Relative free energy (kcal/mol)")
        save(fig, key)
