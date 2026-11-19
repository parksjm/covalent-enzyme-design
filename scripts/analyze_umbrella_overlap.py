#!/usr/bin/env python3
"""Generate umbrella-window overlap data and plots for one campaign stage."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.umbrella_overlap import (  # noqa: E402
    campaign_samples,
    overlap_summary,
    write_overlap_data,
)

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Parse a campaign, calculate overlap, and render its diagnostic plot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--stage", default="03_production_2ps")
    parser.add_argument("--discard-fraction", type=float, default=0.0)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--plot-python", type=Path, required=True)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    windows = campaign_samples(
        args.campaign,
        stage=args.stage,
        discard_fraction=args.discard_fraction,
        allow_incomplete=args.allow_incomplete,
    )
    summaries, edges = overlap_summary(windows)
    args.output.mkdir(parents=True, exist_ok=True)
    summary_path, _ = write_overlap_data(summaries, edges, args.output)
    samples_path = args.output / "overlap_samples.json"
    samples_path.write_text(
        json.dumps(
            [
                {
                    **{key: value for key, value in row.items() if key != "samples_A"},
                    "samples_A": row["samples_A"].tolist(),
                }
                for row in windows
            ],
            indent=2,
        )
        + "\n"
    )
    plot = args.output / "window_overlap.png"
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    subprocess.run(
        [
            str(args.plot_python),
            str(Path(__file__).with_name("plot_umbrella_overlap.py")),
            str(samples_path),
            str(summary_path),
            str(plot),
            "--title",
            args.title,
        ],
        env=environment,
        check=True,
    )
    files = sorted(args.output.glob("overlap_*")) + [plot]
    manifest = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}"
        for path in files
    ]
    (args.output / "MANIFEST.sha256").write_text("\n".join(manifest) + "\n")
    LOGGER.info("Analyzed %d windows in %s", len(windows), args.output)


if __name__ == "__main__":
    main()
