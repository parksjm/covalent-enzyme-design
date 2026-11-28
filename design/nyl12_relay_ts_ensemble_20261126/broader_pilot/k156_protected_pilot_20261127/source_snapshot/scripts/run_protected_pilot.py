#!/usr/bin/env python3
"""Run the root-released protected pilot through the existing broad runner."""

from __future__ import annotations

import argparse
from functools import partial
import importlib.util
import json
import logging
from pathlib import Path
import sys

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE / "src"))

from nyl12_covalent.relay_pilot_namespace import namespace_pool  # noqa: E402


def main() -> None:
    """Namespace raw lineage before pooling and execute the unchanged broad run."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    spec = importlib.util.spec_from_file_location(
        "protected_broad_runner", SOURCE / "scripts/run_relay_broad_pilot.py"
    )
    if spec is None or spec.loader is None:
        raise ImportError("Frozen broad runner cannot be loaded")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.pool_samples = partial(
        namespace_pool, pool=runner.pool_samples, namespace="prot101"
    )
    runner.run(manifest, args.output)


if __name__ == "__main__":
    main()
