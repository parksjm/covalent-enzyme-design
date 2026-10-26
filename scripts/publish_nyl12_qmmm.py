"""Maintain the intuitive QM/MM inspection layout alongside a running path."""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.qmmm_inspection_package import build_qmmm_package  # noqa: E402


def main() -> None:
    """Publish once or follow finalized path inspection updates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    config = json.loads(args.config.read_text())
    run = Path(config["run"])
    package_stamp = None
    while True:
        generation = (run / "live_inspection/latest").resolve(strict=True)
        data = json.loads((generation / "reaction_distances.json").read_text())
        output = build_qmmm_package(args.config.resolve())
        stamp = (output / "package_sources.json").stat().st_mtime_ns
        if stamp != package_stamp:
            logging.info(
                "QM/MM inspection updated: %s (%d frames)", output, len(data["rows"])
            )
            package_stamp = stamp
        state = json.loads((run / "status.json").read_text())
        # A controller updates terminal status before its final inspection export.
        if not args.watch or (
            state["state"] in {"completed", "stopped"} and data["run_state"] == state
        ):
            break
        time.sleep(15)


if __name__ == "__main__":
    main()
