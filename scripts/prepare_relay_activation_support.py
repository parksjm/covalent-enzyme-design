"""Prepare the reviewed33-stateactivation-supportchild; never launch MD."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.relay_activation_support import prepare_support_targets  # noqa: E402


def main() -> None:
    """Read an immutable reviewed plan and write new native inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    args = parser.parse_args()
    prepare_support_targets(args.plan)


if __name__ == "__main__":
    main()
