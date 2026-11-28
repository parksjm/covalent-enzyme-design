"""Analyze the complete matched WT cohort after all twelve MD runs finish."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.matched_wt_occupancy import (  # noqa: E402
    summarize_states,
    validate_config,
)
from nyl12_covalent.matched_wt_occupancy_io import analyze_case, digest  # noqa: E402
from nyl12_covalent.modular_es_snapshot_selection import TIERS  # noqa: E402

LOGGER = logging.getLogger(__name__)


def plot_occupancy(states: dict, output: Path) -> None:
    """Plot equally weighted replica means with replica SD bars.

    Args:
        states: Per-state independent-replica summaries.
        output: New PNG output path.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(9, 4.5))
    keys = ("tier_a", "tier_b_exclusive", "tier_a_or_b", "tier_a_or_b_and_native")
    labels = ("Tier A", "Tier B exclusive", "Tier A or B", "A or B + native Thr")
    x = np.arange(len(states))
    for ordinal, (key, label) in enumerate(zip(keys, labels, strict=True)):
        values = [s["observables"][key] for s in states.values()]
        ax.bar(
            x + (ordinal - 1.5) * 0.19,
            [v["mean"] for v in values],
            0.19,
            yerr=[v["sample_sd"] for v in values],
            capsize=2,
            label=label,
        )
    ax.set_xticks(x, [name.replace("_", "\n", 1) for name in states])
    ax.set_ylabel("Fraction of all saved production frames")
    ax.set_ylim(0, 1.05)
    ax.set_title("Matched WT: four independent 100 ns replicas per state")
    ax.legend(fontsize=8)
    fig.text(
        0.5,
        0.01,
        "Bars: replica mean; error bars: sample SD (n=4), not frame confidence intervals.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    """Validate the cohort and write a new immutable occupancy result directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    config = json.loads(args.config.read_text())
    validate_config(config)
    config_hash = digest(args.config)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    replicas = []
    for case in config["cases"]:
        result = analyze_case(case, config, args.output)
        replicas.append(result)
        (args.output / f"{case['name']}_summary.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        LOGGER.info("Completed %s: %s", case["name"], result["fractions"])
    sources = [
        Path(__file__),
        *[
            Path(__file__).resolve().parents[1] / "src" / "nyl12_covalent" / name
            for name in (
                "matched_wt_occupancy.py",
                "matched_wt_occupancy_io.py",
                "leaving_group_relays.py",
                "leaving_group_relay_mapping.py",
                "nalpha_water_relay.py",
                "catalytic_alignment.py",
                "catalytic_alignment_geometry.py",
                "active_site_geometry.py",
                "modular_es_snapshot_selection.py",
            )
        ],
    ]
    states = summarize_states(replicas)
    summary = {
        "config_sha256": config_hash,
        "source_sha256": {str(p.resolve()): digest(p) for p in sources},
        "tier_thresholds": {key: asdict(v) for key, v in TIERS.items()},
        "native_gates": config["native_gates"],
        "replicas": replicas,
        "states": states,
        "interpretation": "Fixed-protonation geometric occupancies, not kinetics or protonation probabilities. n=4 independent replicas; sample SD/range describe between-replica spread. Correlated frames are not independent observations. Ten-ns blocks describe stability, not sampling convergence. Tier B is exclusive; A or B is inclusive. Phi signs are separate and do not affect abs-phi tier gates.",
    }
    if digest(args.config) != config_hash:
        raise ValueError("Configuration changed during analysis")
    plot_occupancy(states, args.output / "occupancy.png")
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
