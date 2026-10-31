#!/usr/bin/env python3
"""Prepare validated modular PA6 L4 Nyl12 MD systems and replicas."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.modular_pa6_md_setup import (  # noqa: E402
    SystemSpec,
    prepare_campaign,
)

LOGGER = logging.getLogger("prepare_nyl12_modular_pa6_md")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tleap", required=True, type=Path)
    parser.add_argument("--pmemd-cuda", required=True, type=Path)
    parser.add_argument("--replicas", type=int, default=3)
    parser.add_argument(
        "--system",
        action="append",
        required=True,
        metavar="LABEL:PDB:MANIFEST",
        help="Accepted converted PDB and conversion manifest",
    )
    return parser


def _system_spec(value: str) -> SystemSpec:
    label, separator, remainder = value.partition(":")
    pdb_text, separator_two, manifest_text = remainder.partition(":")
    if (
        not label
        or not separator
        or not pdb_text
        or not separator_two
        or not manifest_text
    ):
        raise ValueError(f"Invalid --system value: {value!r}")
    return SystemSpec(label, Path(pdb_text), Path(manifest_text))


def main(argv: list[str] | None = None) -> int:
    """Run the modular PA6 campaign builder.

    Args:
        argv: Optional command arguments.

    Returns:
        Process exit status.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    manifest = prepare_campaign(
        [_system_spec(value) for value in args.system],
        args.release_dir,
        args.output_dir,
        args.tleap,
        args.pmemd_cuda,
        replicas=args.replicas,
    )
    LOGGER.info(
        "prepared %d systems and %d replicas in %s",
        len(manifest["systems"]),
        len(manifest["replica_directories"]),
        args.output_dir,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
