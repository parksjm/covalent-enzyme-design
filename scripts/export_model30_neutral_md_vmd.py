#!/usr/bin/env python3
"""Export stripped 100 ns VMD movies for model30 neutral MD replicas."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.md_vmd_export import export_replica  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Parse paths and export every declared neutral replica."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--md-root", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--stride", type=int, default=25)
    parser.add_argument(
        "--replica", action="append", default=["neutral_replica_01", "neutral_replica_02"]
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    seen: set[str] = set()
    for name in args.replica:
        if name in seen:
            continue
        seen.add(name)
        replica = args.md_root / name
        output = replica / "vmd_stripped_100ns"
        manifest = export_replica(replica, args.cpptraj, output, args.stride)
        LOGGER.info("Verified %s with %d frames", name, manifest["frames"])


if __name__ == "__main__":
    main()
