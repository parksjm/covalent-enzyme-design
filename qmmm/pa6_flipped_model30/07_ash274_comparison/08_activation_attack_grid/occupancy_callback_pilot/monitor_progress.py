"""Plot only accepted pilot iterations and export stable accepted snapshots."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from export_support import (
    HERE,
    ROOT,
    csv_table,
    digest,
    export_movie,
    snapshot_exporter,
    snapshot_frames,
)

LOGGER = logging.getLogger(__name__)
JOURNAL = ROOT / "output_agent_actions/20261021.jsonl"


def journal(phase: str, action: str) -> None:
    """Record exporter intent/outcome; stop before writes if journaling fails."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(JOURNAL),
            "--agent",
            "occupancy_pilot_monitor",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "tool",
            "--reason",
            "Accepted-only pilot inspection export",
            "--rollback",
            "Generated progress plots/movie can be regenerated; retain source snapshots",
        ],
        check=True,
    )


def accepted_rows(pilot: Path) -> list[dict[str, Any]]:
    """Read complete authoritative lines, excluding an in-progress final line."""
    rows = []
    for stage in ("point_A", "point_B", "point_C", "point_D"):
        path = pilot / stage / "accepted_iterations.jsonl"
        if not path.exists():
            continue
        for line in path.read_text().splitlines(keepends=True):
            if not line.endswith("\n"):
                continue
            record = json.loads(line)
            rows.append({"stage": stage, **record})
    return rows


def plot_progress(output: Path, rows: list[dict[str, Any]]) -> None:
    """Plot physical/total/bias energies, both occupancies and attack distance."""
    csv_table(output / "accepted_iterations.csv", rows)
    x = range(len(rows))
    baseline = rows[0]["physical_energy_kcal_mol"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    for key, label in [
        ("physical_energy_kcal_mol", "Physical ΔE"),
        ("total_energy_kcal_mol", "Total ΔE"),
    ]:
        axes[0].plot(x, [r[key] - baseline for r in rows], label=label)
    axes[0].plot(x, [r["bias_energy_kcal_mol"] for r in rows], label="Bias E", ls=":")
    axes[0].set_ylabel("Energy (kcal/mol)")
    for key, label in [
        ("q", "q = CN_N − CN_OG1"),
        ("CN_N", "CN_N"),
        ("CN_OG1", "CN_OG1"),
    ]:
        axes[1].plot(x, [r[key] for r in rows], label=label)
    axes[1].set_ylabel("Smooth coordination")
    axes[2].plot(x, [r["r_A"] for r in rows], label="OG1–substrate C")
    axes[2].plot(
        x, [r["r_A"] - r["r_target_lag_A"] for r in rows], label="Target", ls=":"
    )
    axes[2].set(
        ylabel="Distance (Å)",
        xlabel="Recorded accepted states (initial/final included)",
    )
    for index in range(1, len(rows)):
        if rows[index]["stage"] != rows[index - 1]["stage"]:
            for axis in axes:
                axis.axvline(index, c="0.7", lw=0.7)
    for axis in axes:
        axis.legend(loc="best")
        axis.grid(alpha=0.2)
    fig.suptitle("Occupancy-bias diagnostic minimization; not dynamics or PMF")
    fig.tight_layout()
    fig.savefig(output / "accepted_progress.png", dpi=170)
    fig.savefig(output / "accepted_progress.svg")
    plt.close(fig)


def refresh(pilot: Path, output: Path, previous: dict[str, tuple[int, int]]) -> bool:
    """Refresh plots and snapshots stable over two checkpoints; return terminal state."""
    rows = accepted_rows(pilot)
    output.mkdir(parents=True, exist_ok=True)
    snapshot_exporter(output)
    stability = output / "snapshot_stability.json"
    if not previous and stability.is_file():
        previous.update(
            {
                key: tuple(value)
                for key, value in json.loads(stability.read_text()).items()
            }
        )
    if rows:
        plot_progress(output, rows)
        frames = snapshot_frames(pilot, rows)
        stable, observed = [], {}
        terminal = (pilot / "completion.json").is_file() or any(
            (pilot / f"point_{stage}" / "process_failure.json").is_file()
            for stage in "ABCD"
        )
        for path, row in frames:
            stat = path.stat()
            stamp = stat.st_size, stat.st_mtime_ns
            observed[str(path)] = stamp
            if terminal or previous.get(str(path)) == stamp:
                stable.append((path, row))
        previous.clear()
        previous.update(observed)
        stability.write_text(json.dumps(observed, indent=2) + "\n")
        export_movie(output, stable)
    terminal = (pilot / "completion.json").is_file() or any(
        (pilot / f"point_{stage}" / "process_failure.json").is_file()
        for stage in "ABCD"
    )
    manifest = {
        "accepted_records": len(rows),
        "terminal": terminal,
        "energy_reference": "First recorded physical energy; bias plotted absolute",
        "source_jsonl_sha256": {
            str(p): digest(p)
            for p in sorted(pilot.glob("point_*/accepted_iterations.jsonl"))
        },
        "movie_frames": "Stable accepted initial/every25/final only; no trial evaluations",
    }
    (output / "progress.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("accepted_records=%d terminal=%s", len(rows), terminal)
    return terminal


def main() -> None:
    """Refresh once or monitor at bounded checkpoints; never launch a scientific job."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, default=HERE / "pilot")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=float, default=120.0)
    parser.add_argument("--hours", type=float, default=48.0)
    args = parser.parse_args()
    if args.interval <= 0 or args.hours <= 0:
        parser.error("interval and hours must be positive")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    pilot = args.pilot.resolve()
    output = args.output.resolve() if args.output else pilot / "progress"
    deadline = time.monotonic() + args.hours * 3600
    previous = {}
    while True:
        journal(
            "intent", f"Refresh accepted pilot plots/movie from {pilot} into {output}"
        )
        try:
            terminal = refresh(pilot, output, previous)
        except Exception as error:
            journal("outcome", f"Pilot export failed: {type(error).__name__}: {error}")
            raise
        journal("outcome", f"Refreshed accepted pilot plots/movie; terminal={terminal}")
        if terminal or not args.watch:
            return
        if time.monotonic() >= deadline:
            raise TimeoutError("Pilot export monitor reached bounded deadline")
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
