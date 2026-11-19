"""Plot saved directional derivatives with distinct energy definitions."""

from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    """Write a static two-CV comparison without mixing gradient definitions."""
    here = Path(__file__).resolve().parent
    data = json.loads((here / "results.json").read_text())
    rows = data["projections"]
    states = ["01_es", "02_activated", "03_attack"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for col, (cv, unit) in enumerate(
        (("q_occ", "kcal/mol per q"), ("attack_hold_q_occ", "kcal/mol/Å"))
    ):
        for row, methods in enumerate(
            (("DFT_fixed_PC_QM", "GFN2_fixed_PC_QM"), ("DFTB3_native_TOTAL",))
        ):
            ax = axes[row, col]
            for method in methods:
                y = [
                    next(
                        r["derivative"]
                        for r in rows
                        if r["state"] == s
                        and r["method"] == method
                        and r["direction"] == cv
                    )
                    for s in states
                ]
                ax.plot(range(3), y, "o-", label=method)
            ax.axhline(0, color="gray", linewidth=0.6)
            ax.set_title(cv)
            ax.set_ylabel(unit)
            ax.set_xticks(range(3), ["E–S", "activated", "attack"])
            ax.legend(fontsize=8)
    fig.suptitle(
        "Local energy slopes; positive = uphill in increasing CV\nTop: fixed-PC QM/link; bottom: native total QM/MM (separate definition)"
    )
    fig.tight_layout()
    fig.savefig(here / "directional_derivatives.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
