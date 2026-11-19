#!/usr/bin/env python3
"""Prepare and run one modular-PA6 E-S MM minimization when a GPU is idle."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.es_classical_relaxation import prepare_relaxation  # noqa: E402
from scripts.run_when_gpu_idle import wait_and_run  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Prepare immutable inputs, wait non-preemptively, and run pmemd.cuda."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--engine", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case-label", required=True)
    parser.add_argument("--max-cycles", type=int, default=10_000)
    parser.add_argument("--gpu", required=True, type=int)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--nvidia-smi", type=Path, default=Path("/usr/bin/nvidia-smi"))
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    metadata = prepare_relaxation(
        args.topology.resolve(),
        args.coordinates.resolve(),
        args.engine.resolve(),
        args.output.resolve(),
        args.case_label,
        args.max_cycles,
    )
    status = wait_and_run(
        args.nvidia_smi.resolve(),
        args.gpu,
        args.poll_seconds,
        metadata["command"],
    )
    if status:
        raise RuntimeError(f"MM minimization exited with status {status}")
    LOGGER.info("MM minimization completed: %s", args.output)


if __name__ == "__main__":
    main()
