#!/usr/bin/env python3
"""Compare complete native NEB trial reports; never launch or signal engines."""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path
import shutil
import sys
import time
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa:E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_spacing_monitor import (  # noqa:E402
    atomic_json,
    digest,
    update_trial,
    write_csv,
)
from scripts.record_agent_action import append_action  # noqa:E402

LOGGER = logging.getLogger(__name__)
DEFAULT = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/04_neb/spacing_tests"
)


def source_snapshot(output: Path) -> None:
    """Retain immutable numbered source milestones when this monitor changes."""
    paths = [
        Path(__file__).resolve(),
        ROOT / "src/nyl12_covalent/neb_spacing_monitor.py",
        ROOT / "src/nyl12_covalent/neb_progress.py",
    ]
    hashes = {str(p): digest(p) for p in paths}
    folder = output / "source_milestones"
    folder.mkdir(exist_ok=True)
    milestones = sorted(folder.glob("milestone_*"))
    if (
        milestones
        and json.loads((milestones[-1] / "manifest.json").read_text())["source_sha256"]
        == hashes
    ):
        return
    destination = folder / f"milestone_{len(milestones) + 1:04d}"
    destination.mkdir()
    for path in paths:
        shutil.copy2(path, destination / path.name)
    atomic_json(
        destination / "manifest.json", {"source_sha256": hashes, "command": sys.argv}
    )


def plot_comparison(output: Path, profiles: dict[str, dict[str, Any]]) -> None:
    """Plot current physical profiles, maxima and Max/RMS target ratios by cycle."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 11))
    histories = []
    for index, (name, profile) in enumerate(profiles.items()):
        color = f"C{index}"
        images = profile["images"]
        n = len(images)
        axes[0].plot(
            [i / (n - 1) for i in range(n)],
            [i["relative_energy_kcal_mol"] for i in images],
            marker=".",
            color=color,
            label=f"{name} ({n}; cycle {profile['cycle']})",
        )
        path = output / name / "history.jsonl"
        history = [json.loads(line) for line in path.read_text().splitlines()]
        histories.extend({"trial": name, **r} for r in history)
        cycles = [r["cycle"] for r in history]
        axes[1].plot(
            cycles,
            [r["physical_max_kcal_mol"] for r in history],
            color=color,
            label=name,
        )
        axes[2].plot(
            cycles,
            [r["max_gradient_ratio"] for r in history],
            color=color,
            label=f"{name}: Max",
        )
        axes[2].plot(
            cycles,
            [r["rms_gradient_ratio"] for r in history],
            color=color,
            ls="--",
            label=f"{name}: RMS",
        )
    axes[0].set(
        xlabel="Regular image fraction (E–S → acyl)", ylabel="Physical ΔE (kcal/mol)"
    )
    axes[1].set(xlabel="Native convergence cycle", ylabel="Profile maximum (kcal/mol)")
    axes[2].set(
        xlabel="Native convergence cycle",
        ylabel="Gradient / native target",
        yscale="log",
    )
    axes[2].axhline(1, color="0.5", lw=0.8)
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    fig.suptitle(
        "NEB spacing trials — unconverged path candidates; no barrier acceptance"
    )
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        pending = output / f"comparison.pending.{suffix}"
        fig.savefig(pending, dpi=160)
        pending.replace(output / f"comparison.{suffix}")
    plt.close(fig)
    write_csv(output / "comparison_history.csv", histories)


def refresh(root: Path, output: Path) -> dict[str, Any]:
    """Discover trials explicitly and update only this isolated monitor's outputs.

    Args:
        root: Spacing-test directory containing immediate trial/config.json files.
        output: Dedicated plots/history/source milestone destination.

    Returns:
        Per-trial complete-report/waiting states and checkpoint comparison.
    """
    output.mkdir(parents=True, exist_ok=True)
    source_snapshot(output)
    configs = sorted(root.glob("*/config.json"))
    statuses, profiles = [], {}
    for config in configs:
        name = config.parent.name
        statuses.append(update_trial(config, output / name))
        path = output / name / "current_profile.json"
        if path.is_file():
            profiles[name] = json.loads(path.read_text())
    comparison = {
        "status": "waiting_for_configs"
        if not configs
        else "monitoring_complete_reports",
        "trials": statuses,
        "checkpoints": [50, 100, 250],
        "interpretation": "Native physical energies, NEB springs excluded; no validated barriers or automatic acceptance",
        "coordinate_policy": "Chemical gaps only from independently coherent same-cycle exported bands",
        "vmd_policy": "Existing watch_relay_acyl_neb.py supports24regular images only;48requires separate coherent export",
    }
    comparison["profile_identity"] = {
        name: p["native_block_sha256"] for name, p in profiles.items()
    }
    prior = output / "comparison.json"
    changed = (
        not prior.exists()
        or json.loads(prior.read_text()).get("profile_identity")
        != comparison["profile_identity"]
    )
    atomic_json(prior, comparison)
    if profiles and changed:
        plot_comparison(output, profiles)
    LOGGER.info("trials=%d profiles=%d", len(configs), len(profiles))
    return comparison


def main(argv: list[str] | None = None) -> None:
    """Run a bounded read-only monitor; root launches scientific jobs separately."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=int, choices=(30, 60), default=30)
    parser.add_argument("--hours", type=float, default=12)
    parser.add_argument(
        "--journal", type=Path, default=ROOT / "output_agent_actions/20261021.jsonl"
    )
    args = parser.parse_args(argv)
    if not math.isfinite(args.hours) or args.hours <= 0:
        parser.error("--hours must be finite positive")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    root = args.root.resolve()
    output = args.output.resolve() if args.output else root / "plots"
    deadline = time.monotonic() + args.hours * 3600
    while True:
        append_action(
            args.journal,
            "neb_spacing_monitor",
            "intent",
            f"Refresh isolated NEB spacing comparison {root}",
            "tool",
            "Authorized complete-report analysis only",
            "Generated monitor outputs can be regenerated; native files untouched",
        )
        try:
            status = refresh(root, output)
        except Exception as error:
            append_action(
                args.journal,
                "neb_spacing_monitor",
                "outcome",
                f"Comparison failed: {type(error).__name__}: {error}",
                "tool",
                "Explicit failure; no method substitution",
                "Native scientific jobs untouched",
            )
            raise
        append_action(
            args.journal,
            "neb_spacing_monitor",
            "outcome",
            f"Refreshed {len(status['trials'])} trial states",
            "tool",
            "Native output consumed incrementally; no engine calls",
            "Generated monitor outputs only",
        )
        if not args.watch:
            return
        if time.monotonic() >= deadline:
            LOGGER.info("Bounded monitor wall time reached; scientific jobs continue")
            return
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
