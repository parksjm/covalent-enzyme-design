"""Run exactly one approved boundary history after separate release authorization."""

from __future__ import annotations

import argparse
from pathlib import Path

from nyl12_covalent.relay_activation_support_runtime import run_reviewed_target


def main() -> None:
    """Invoke existing hash-verifying runtime with the exact six-history count."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    run_reviewed_target(args.prepared, args.label, args.sha256, 6)


if __name__ == "__main__":
    main()
