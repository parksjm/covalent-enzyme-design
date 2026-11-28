"""Execute exactly one reviewed33-state support window through30productionps."""

from __future__ import annotations

import argparse
from pathlib import Path

from nyl12_covalent.relay_activation_support_runtime import run_support_target


def main() -> None:
    """Run one hash-pinned native support task; caller owns launch authorization."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    run_support_target(args.prepared, args.label, args.sha256)


if __name__ == "__main__":
    main()
