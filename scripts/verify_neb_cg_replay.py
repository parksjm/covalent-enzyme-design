"""Read-only printed-profile comparison of an identical-start higherbudget CG replay."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from nyl12_covalent.neb_replay import compare_replay, parse_replay_profiles  # noqa:E402

LOGGER = logging.getLogger(__name__)


def save(path: Path, record: dict, *, exclusive: bool = False) -> None:
    """Write immutable checkpoints or atomically replace the latest monitor report."""
    target = path if exclusive else path.with_suffix(".pending")
    with target.open("x" if exclusive else "w") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")
    if not exclusive:
        target.replace(path)


def watch_replay(
    pilot_path: Path,
    replay_path: Path,
    output: Path,
    *,
    watch: bool = False,
    interval: float = 60.0,
) -> dict:
    """Pin immutable pilot bytes and compare incremental complete replay profiles.

    Args:
        pilot_path: Finished pilot minimized.mdout, never modified.
        replay_path: Higherbudget replay minimized.mdout, read-only.
        output: New evidence directory; may contain an outer launcher job folder.
        watch: Continue until a divergence or first completed post-pilot cycle.
        interval: Explicit positive monitor interval in seconds.

    Returns:
        Last comparison with independent replay/force/scientific states.

    Raises:
        ValueError: If pilot, input precision or monitoring interval is invalid.
        FileExistsError: If immutable pilot provenance already exists.
    """
    if interval <= 0:
        raise ValueError("Monitor interval must be positive")
    pilot_path, replay_path, output = (
        p.resolve() for p in (pilot_path, replay_path, output)
    )
    output.mkdir(parents=True, exist_ok=True)
    pilot_bytes = pilot_path.read_bytes()
    pilot_hash = hashlib.sha256(pilot_bytes).hexdigest()
    pilot = parse_replay_profiles(pilot_bytes)
    if not pilot:
        raise ValueError("Pilot contains no complete gradient/profile cycles")
    with (output / "pilot_reference.mdout").open("xb") as stream:
        stream.write(pilot_bytes)
    sources = {str(pilot_path): pilot_hash}
    recipes = output / "source_snapshot"
    recipes.mkdir()
    for path in (
        Path(__file__),
        ROOT / "src/nyl12_covalent/neb_replay.py",
        ROOT / "src/nyl12_covalent/neb_progress.py",
    ):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        copied = recipes / path.name
        shutil.copy2(path, copied)
        if hashlib.sha256(copied.read_bytes()).hexdigest() != digest:
            raise ValueError("Replay monitor snapshot mismatch")
        sources[str(path)] = digest
    save(
        output / "provenance.json",
        dict(
            pilot=str(pilot_path),
            pilot_sha256=pilot_hash,
            pilot_last_cycle=max(pilot),
            source_sha256=sources,
            comparison_policy="Same cycle;24absoluteenergies at1.1e-7Eh; gradient half-lastdigit sum; no scientific acceptance",
            prefix_read_limit_bytes=len(pilot_bytes) + 1048576,
        ),
        exclusive=True,
    )
    checkpoints = set()
    while True:
        if hashlib.sha256(pilot_path.read_bytes()).hexdigest() != pilot_hash:
            raise ValueError("Immutable pilot changed")
        if replay_path.exists():
            with replay_path.open("rb") as stream:
                payload = stream.read(len(pilot_bytes) + 1048576)
        else:
            payload = b""
        parsed = parse_replay_profiles(payload)
        result = compare_replay(pilot, parsed)
        result.update(
            observed_utc=datetime.now(timezone.utc).isoformat(),
            pilot_sha256=pilot_hash,
            replay_source=str(replay_path),
            replay_prefix_sha256=hashlib.sha256(payload).hexdigest(),
            replay_prefix_bytes=len(payload),
            new_energy_or_force_evaluations=False,
        )
        save(output / "last.json", result)
        with (output / "checks.jsonl").open("a") as stream:
            stream.write(json.dumps(result, allow_nan=False) + "\n")
        for cycle in (1, 50, max(pilot)):
            if cycle in parsed and cycle not in checkpoints:
                folder = output / f"checkpoint_cycle_{cycle:03d}"
                folder.mkdir()
                end = parsed[cycle].end_byte
                (folder / "replay_prefix.mdout").write_bytes(payload[:end])
                row = next((r for r in result["cycles"] if r["cycle"] == cycle), None)
                save(
                    folder / "comparison.json",
                    dict(
                        cycle=cycle,
                        comparison=row,
                        pilot_prefix_sha256=hashlib.sha256(
                            pilot_bytes[: pilot[cycle].end_byte]
                        ).hexdigest(),
                        replay_prefix_sha256=hashlib.sha256(payload[:end]).hexdigest(),
                        status_at_capture=result["status"],
                    ),
                    exclusive=True,
                )
                checkpoints.add(cycle)
        LOGGER.info(
            "Replay%s compared=%d latest=%d energymaxdiff=%.3gEh",
            result["status"],
            result["compared_cycles"],
            result["replay_last_complete_cycle"],
            result["max_energy_difference_Eh"],
        )
        if not watch or result["failed_cycles"] or result["replay_passed_pilot"]:
            return result
        time.sleep(interval)


def main() -> None:
    """Start explicit finishedpilot versus replay comparison; never launch engines."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot", type=Path)
    parser.add_argument("replay", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=float, default=60.0)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    result = watch_replay(
        args.pilot, args.replay, args.output, watch=args.watch, interval=args.interval
    )
    if result["failed_cycles"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
