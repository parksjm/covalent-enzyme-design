"""Render the retained5ps MBAR overlap matrix without rerunning analysis."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

root = Path(__file__).resolve().parent
data = json.loads((root / "mbar_overlap_05ps_source.json").read_text())
matrix = np.asarray(data["overlap"], dtype=float)
plt.rcParams.update({"font.size": 15, "axes.titlesize": 19, "axes.labelsize": 17})
figure, axis = plt.subplots(figsize=(10, 8.5), constrained_layout=True, facecolor="white")
masked = np.ma.masked_equal(matrix, 0.0)
color = plt.get_cmap("viridis").copy()
color.set_bad("white")
artist = axis.imshow(masked, origin="upper", cmap=color, norm=LogNorm(vmin=1e-6, vmax=1, clip=True), interpolation="nearest", extent=(0.5, 57.5, 57.5, 0.5))
axis.axvline(31.5, color="white", linewidth=1.1, linestyle="--")
axis.axhline(31.5, color="white", linewidth=1.1, linestyle="--")
axis.set(xticks=[1, 10, 20, 31, 40, 50, 57], yticks=[1, 10, 20, 31, 40, 50, 57], xlabel="Umbrella state j", ylabel="Umbrella state i", title="MBAR overlap at 5 ps per window\n57 states · 300 K · preliminary")
bar = figure.colorbar(artist, ax=axis, shrink=0.86, pad=0.035, ticks=[1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
bar.set_label(r"Overlap $O_{ij}$ (logarithmic color scale)")
figure.supxlabel("States 1–31: OG1 → N-terminal N; 32–57: N-terminal N → leaving-group N", fontsize=12)
figure.savefig(root / "mbar_overlap_05ps.png", dpi=200, metadata={"Source_SHA256": data["source_sha256"], "Scope": "5ps/window; connected graph does not establish convergence"})
plt.close(figure)
