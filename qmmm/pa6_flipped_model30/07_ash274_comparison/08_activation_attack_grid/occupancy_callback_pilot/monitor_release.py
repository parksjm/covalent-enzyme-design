"""Monitor accepted unbiased D-release records and export stable snapshots."""

from __future__ import annotations

import argparse
import json
import logging
import shutil
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

TITLE = "Unbiased point D release; both biases 0"
LOGGER = logging.getLogger(__name__)


def journal(phase: str, action: str) -> None:
    """Record monitor actions before and after writes; fail on logging errors."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "occupancy_D_release_monitor",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "tool",
            "--reason",
            TITLE,
            "--rollback",
            "Regenerate progress outputs; retain exporter sources and inputs",
        ],
        check=True,
    )


def accepted_rows(release: Path) -> list[dict[str, Any]]:
    """Read complete accepted ledger lines and reject any nonzero bias."""
    path = release / "accepted_iterations.jsonl"
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().splitlines(keepends=True):
        if not line.endswith("\n"):
            continue
        row = json.loads(line)
        for key in ["bias_energy_kcal_mol", "bias_q_kcal_mol", "bias_r_kcal_mol"]:
            if row[key] != 0.0:
                raise ValueError(f"Nonzero release bias: {key}={row[key]}")
        if row["physical_energy_kcal_mol"] != row["total_energy_kcal_mol"]:
            raise ValueError("Release total energy differs from physical energy")
        rows.append({"stage": "unbiased_release", **row})
    return rows


def preserve_sources(output: Path) -> None:
    """Retain the established exporter sources plus this release monitor."""
    snapshot_exporter(output)
    folder = output / "exporter_source"
    identity = {str(Path(__file__).resolve()): digest(Path(__file__))}
    manifest = folder / "release_monitor_manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text())["source_sha256"] != identity:
            raise ValueError("Release monitor changed; use a new progress directory")
        return
    shutil.copy2(Path(__file__), folder / Path(__file__).name)
    manifest.write_text(
        json.dumps(
            {"source_sha256": identity, "command": sys.argv, "state": TITLE}, indent=2
        )
        + "\n"
    )


def plot_progress(output: Path, rows: list[dict[str, Any]]) -> None:
    """Plot physical energy relative to the first release state and geometry."""
    csv_table(output / "accepted_iterations.csv", rows)
    baseline = rows[0]["physical_energy_kcal_mol"]
    x = [row["accepted_iteration"] for row in rows]
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    axes[0].plot(
        x,
        [r["physical_energy_kcal_mol"] - baseline for r in rows],
        label="Physical ΔE from release initial",
    )
    axes[0].plot(
        x, [r["bias_energy_kcal_mol"] for r in rows], label="Bias E = 0", ls=":"
    )
    axes[0].set_ylabel("Energy (kcal/mol)")
    for key in ["q", "CN_N", "CN_OG1"]:
        axes[1].plot(x, [r[key] for r in rows], label=key)
    axes[1].set_ylabel("Smooth coordination")
    axes[2].plot(x, [r["r_A"] for r in rows], label="OG1–substrate C")
    axes[2].set(ylabel="Distance (Å)", xlabel="Accepted release iteration")
    for axis in axes:
        axis.legend(loc="best")
        axis.grid(alpha=0.2)
    fig.suptitle(TITLE + "\nDiagnostic minimization; not dynamics or PMF")
    fig.tight_layout()
    for suffix in ["png", "svg"]:
        fig.savefig(output / f"accepted_progress.{suffix}", dpi=170)
    plt.close(fig)


def label_bundle(output: Path, baseline: float) -> None:
    """Set release-specific VMD state and regenerate affected checksums."""
    bundle = output / "vmd"
    if not bundle.exists():
        return
    path = bundle / "view_path.vmd"
    text = path.read_text().replace(
        "Occupancy pilot accepted snapshots; not dynamics", TITLE
    )
    path.write_text(text)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.update(
        state=TITLE,
        energy_reference="First recorded physical release energy",
        energy_reference_kcal_mol=baseline,
        both_biases_zero=True,
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    checksum = bundle / "MANIFEST.sha256"
    names = [line.split("  ", 1)[1] for line in checksum.read_text().splitlines()]
    checksum.write_text("".join(f"{digest(bundle / name)}  {name}\n" for name in names))


def refresh(release: Path, output: Path, previous: dict[str, tuple[int, int]]) -> bool:
    """Refresh accepted records and snapshots stable across checkpoints."""
    rows = accepted_rows(release)
    terminal = any(
        (release / name).is_file()
        for name in ["release_assessment.json", "process_failure.json"]
    )
    output.mkdir(parents=True, exist_ok=True)
    preserve_sources(output)
    stability = output / "snapshot_stability.json"
    if not previous and stability.exists():
        previous.update(
            {
                key: tuple(value)
                for key, value in json.loads(stability.read_text()).items()
            }
        )
    stable, observed = [], {}
    if rows:
        plot_progress(output, rows)
        for path, row in snapshot_frames(release.parent, rows):
            stat = path.stat()
            stamp = stat.st_size, stat.st_mtime_ns
            observed[str(path)] = stamp
            if terminal or previous.get(str(path)) == stamp:
                stable.append((path, row))
        previous.clear()
        previous.update(observed)
        stability.write_text(json.dumps(observed, indent=2) + "\n")
        export_movie(output, stable)
        label_bundle(output, rows[0]["physical_energy_kcal_mol"])
    manifest = {
        "state": TITLE,
        "accepted_records": len(rows),
        "terminal": terminal,
        "both_biases_zero": True,
        "energy_reference": "First recorded physical release energy",
        "energy_reference_kcal_mol": rows[0]["physical_energy_kcal_mol"]
        if rows
        else None,
        "source_jsonl_sha256": digest(release / "accepted_iterations.jsonl")
        if rows
        else None,
        "movie_frames": "Stable accepted initial/every25/final; no trial evaluations",
        "stable_snapshots": len(stable),
    }
    (output / "progress.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info(
        "accepted_records=%d stable_snapshots=%d terminal=%s",
        len(rows),
        len(stable),
        terminal,
    )
    return terminal


def main() -> None:
    """Refresh once or watch every 120 seconds for at most 48 hours."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--release", type=Path, default=HERE / "pilot/point_D/unbiased_release"
    )
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
    release = args.release.resolve()
    output = args.output.resolve() if args.output else release / "progress"
    deadline = time.monotonic() + args.hours * 3600
    previous: dict[str, tuple[int, int]] = {}
    while True:
        journal("intent", f"Refresh {TITLE}: {release} into {output}")
        try:
            terminal = refresh(release, output, previous)
        except Exception as error:
            journal(
                "outcome", f"Release export failed: {type(error).__name__}: {error}"
            )
            raise
        journal("outcome", f"Refreshed release progress; terminal={terminal}")
        if terminal or not args.watch:
            return
        if time.monotonic() >= deadline:
            raise TimeoutError("Release monitor reached bounded deadline")
        time.sleep(min(args.interval, max(0, deadline - time.monotonic())))


if __name__ == "__main__":
    main()
