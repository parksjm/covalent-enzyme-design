"""Plot snapshotted presentation evidence without recomputing the PMF."""

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    """Render a readable PMF panel with the qualified original-lineage band."""
    base = Path(__file__).resolve().parent
    pooled = json.loads((base / "evidence/relay_point.json").read_text())
    original = json.loads((base / "evidence/relay_original57.json").read_text())
    x = np.array(original["primary"]["centers"])
    y = np.array(pooled["projections"]["hard"]["pmf"], dtype=float)
    bootstrap = original["bootstrap"]["5ps"]
    low = np.array(bootstrap["lower_95_kcal_mol"], dtype=float)
    high = np.array(bootstrap["upper_95_kcal_mol"], dtype=float)
    plt.rcParams.update(
        {
            "font.size": 22,
            "axes.labelsize": 26,
            "xtick.labelsize": 22,
            "ytick.labelsize": 22,
        }
    )
    fig, ax = plt.subplots(figsize=(12.8, 5.1), layout="constrained")
    ax.fill_between(
        x,
        low,
        high,
        color="#9babb7",
        alpha=0.45,
        label="Original 57: pointwise 95% interval",
    )
    ax.plot(x, y, color="#008c95", lw=3, label="Pooled 61: point estimate")
    ax.set(
        xlabel="Path arclength, s",
        ylabel="PMF (kcal mol⁻¹)",
        xlim=(-0.08, 1.1),
        ylim=(-5, 31),
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper left", fontsize=17)
    ax.text(
        0.98,
        0.92,
        "300 K; sampling not converged",
        transform=ax.transAxes,
        ha="right",
        fontsize=21,
    )
    fig.savefig(base / "assets/relay_pmf_readable.png", dpi=200)
    fig.savefig(base / "assets/relay_pmf_readable.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
