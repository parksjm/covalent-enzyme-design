"""Plot archived NEB profiles and peak evolution without changing run data."""

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    """Render early/recent band profiles and sampled maximum energies."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plots", type=Path)
    args = parser.parse_args()
    profiles = {}
    sources = []
    for path in sorted((args.plots / "reports").glob("*/neb_progress.csv")):
        raw = path.read_bytes()
        rows = list(csv.DictReader(raw.decode().splitlines()))
        regular = [r for r in rows if r["climbing"].lower() == "false"]
        if not regular:
            raise ValueError(f"No regular band images: {path}")
        cycle = int(regular[0]["cycle"])
        x = [int(r["image"]) for r in regular]
        y = [float(r["relative_energy_kcal_mol"]) for r in regular]
        if x != list(range(1, len(x) + 1)) or abs(y[0]) > 1e-6:
            raise ValueError(f"Invalid image ordering or energy reference: {path}")
        profiles[cycle] = (x, y)
        sources.append({"path": str(path), "sha256": hashlib.sha256(raw).hexdigest()})
    cycles = sorted(profiles)
    if not cycles:
        raise ValueError("No archived profiles")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    groups = [[c for c in cycles if c <= 43], cycles[-6:]]
    for ax, group, title in zip(axes[:2], groups, ["Early relaxation", "Latest sampled profiles"], strict=True):
        for color, cycle in zip(plt.cm.viridis([i / max(1, len(group) - 1) for i in range(len(group))]), group, strict=True):
            ax.plot(*profiles[cycle], marker="o", markersize=3, color=color, label=f"Iteration {cycle}")
        ax.set(title=title, xlabel="NEB image", ylabel="Relative potential energy (kcal/mol)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)
    peaks = [max(profiles[c][1]) for c in cycles]
    axes[2].plot(cycles, peaks, "o-", color="#a33236", markersize=4)
    axes[2].set(title="Band maximum during optimization", xlabel="NEB iteration", ylabel="Maximum relative energy (kcal/mol)")
    axes[2].grid(alpha=0.2)
    axes[2].annotate(f"{peaks[-1]:.1f}", (cycles[-1], peaks[-1]), xytext=(-30, 12), textcoords="offset points")
    fig.suptitle("DFTB3/MM acylenzyme NEB — evolution of the potential-energy profile\nEach profile referenced to its own E–S endpoint; unconverged band maxima are not activation barriers", fontsize=12)
    for suffix in ("png", "svg"):
        fig.savefig(args.plots / f"neb_evolution.{suffix}", dpi=170)
    metadata = {"sources": sources, "cycles": cycles, "peak_kcal_mol": peaks, "reference": "same-band image 1", "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.plots / "neb_evolution.json").write_text(json.dumps(metadata, indent=2) + "\n")
    logging.info("Plotted %d profiles; latest cycle %d, maximum %.3f kcal/mol", len(cycles), cycles[-1], peaks[-1])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
