"""Measure a stable unrestrained restart without calling it an optimized endpoint."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import tempfile

from parmed.amber import AmberParm

from run_optimization import geometry


def main() -> None:
    """Save a timestamped geometry observation with the checkpoint identity."""
    root = Path(__file__).resolve().parents[2]
    run = (
        root
        / "output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/phase1_unrestrained"
    )
    metadata = json.loads((run / "preparation.json").read_text())
    source = run / "optimized.ncrst"
    before = source.stat()
    data = source.read_bytes()
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError("Restart changed while reading; inspect a later snapshot")
    with tempfile.TemporaryDirectory(prefix="nyl12-checkpoint-") as scratch:
        snapshot = Path(scratch) / "checkpoint.ncrst"
        snapshot.write_bytes(data)
        parm = AmberParm(metadata["topology"], xyz=str(snapshot))
        report = {
            "observed_utc": datetime.now(timezone.utc).isoformat(),
            "source": str(source),
            "source_sha256": hashlib.sha256(data).hexdigest(),
            "source_mtime_ns": after.st_mtime_ns,
            "scope": "Latest restart checkpoint, not necessarily final or converged",
            "distances_A": geometry(parm),
        }
    output = run / "checkpoint_geometry.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    logging.info("Checkpoint geometry: %s", report["distances_A"])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
