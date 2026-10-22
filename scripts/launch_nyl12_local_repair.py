"""Snapshot a fresh local TI1 repair and launch it persistently in tmux."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.jobs import launch_job  # noqa: E402
from nyl12_covalent.ti1_attempt import _snapshot  # noqa: E402


def main() -> None:
    """Retain actual inputs/source and launch a new preparation or engine job."""
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("topology", "coordinates", "cif", "amber-bin", "output"):
        parser.add_argument(f"--{argument}", type=Path, required=True)
    parser.add_argument("--session", required=True)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    root = args.output.absolute()
    root.mkdir(parents=True, exist_ok=False)
    sources = _snapshot(root)
    (root / "source_hashes.json").write_text(json.dumps(sources, indent=2) + "\n")
    inputs = {}
    for name, path in (
        ("source.prmtop", args.topology),
        ("source.rst7", args.coordinates),
        ("af3.cif", args.cif),
    ):
        destination = root / name
        shutil.copy2(path, destination)
        inputs[name] = destination
    command = [
        sys.executable,
        str(root / "source/scripts/run_nyl12_local_repair.py"),
        "--topology",
        str(inputs["source.prmtop"]),
        "--coordinates",
        str(inputs["source.rst7"]),
        "--cif",
        str(inputs["af3.cif"]),
        "--amber-bin",
        str(args.amber_bin.absolute()),
        "--output",
        str(root / "calculation"),
    ]
    if args.prepare_only:
        command.append("--prepare-only")
    launch_job(
        args.session,
        root / "job",
        root,
        command,
        env_overrides={"OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "1"},
        input_paths=[*inputs.values(), root / "source_hashes.json"],
    )
    logging.info("Launched %s; retained artifacts at %s", args.session, root)


if __name__ == "__main__":
    main()
