#!/usr/bin/env python3
"""Analyze completed production blocks once or watch; never launch scientific jobs."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
import math
import subprocess
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_umbrella_analysis import analyze, collect, digest  # noqa:E402


def refresh(
    prepared: Path, output: Path, bins: int = 20, plot_python: Path | None = None
) -> dict:
    """Publish a new analysis generation only when completed production inputs change."""
    cfg, windows = collect(prepared)
    source_paths = [
        Path(__file__),
        *[
            ROOT / "src/nyl12_covalent" / name
            for name in (
                "neb_umbrella_analysis.py",
                "neb_umbrella_analysis_plots.py",
                "neb_umbrella_bias.py",
            )
        ],
    ]
    source_paths.append(ROOT / "scripts/plot_neb_umbrella_analysis.py")
    source_identity = {str(p): digest(p) for p in source_paths}
    identity = hashlib.sha256(
        json.dumps(
            {
                "window_sources": {w["label"]: w["source_sha256"] for w in windows},
                "prepared_sha256": digest(prepared),
                "bins": bins,
                "analysis_sha256": source_identity,
                "plot_python": str(plot_python) if plot_python else None,
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()
    prior = output / "latest.json"
    if prior.exists() and json.loads(prior.read_text())["input_identity"] == identity:
        return json.loads(prior.read_text())
    result = analyze(cfg, windows, bins)
    result.update(
        input_identity=identity,
        as_of_utc=datetime.now(timezone.utc).isoformat(),
        prepared_sha256=digest(prepared),
    )
    output.mkdir(parents=True, exist_ok=True)
    generation = output / f"update_{len(list(output.glob('update_*'))) + 1:05d}"
    generation.mkdir()
    result["analysis_source_sha256"] = source_identity
    result["generation"] = str(generation)
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    (generation / "analysis.json").write_text(text)
    if plot_python is not None:
        subprocess.run(
            [
                str(plot_python),
                str(ROOT / "scripts/plot_neb_umbrella_analysis.py"),
                "--analysis",
                str(generation / "analysis.json"),
            ],
            check=True,
        )
    pending = output / "latest.pending.json"
    pending.write_text(text)
    pending.replace(prior)
    link = output / "latest"
    pending_link = output / "latest.pending"
    pending_link.symlink_to(generation.name)
    pending_link.replace(link)
    logging.info(
        "Analyzed %d/%dwindows; %s", len(windows), len(cfg["windows"]), generation
    )
    return result


def main() -> None:
    """Provide a stable post-block one-shot hook or bounded read-only watcher."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prepared", required=True, type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--bins", type=int, default=20)
    p.add_argument("--watch", action="store_true")
    p.add_argument(
        "--plot-python",
        type=Path,
        help="Explicit existing Matplotlib Python; omit for JSON only",
    )
    p.add_argument("--interval", type=int, choices=(30, 60), default=60)
    p.add_argument("--hours", type=float, default=48)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO)
    if not math.isfinite(args.hours) or args.hours <= 0:
        p.error("--hours must be positive")
    deadline = time.monotonic() + 3600 * args.hours
    while True:
        refresh(
            args.prepared.resolve(),
            args.output or args.prepared.parent / "analysis",
            args.bins,
            args.plot_python,
        )
        if not args.watch or time.monotonic() >= deadline:
            return
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
