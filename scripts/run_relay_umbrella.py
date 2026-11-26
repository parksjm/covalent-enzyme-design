#!/usr/bin/env python3
"""Launch bounded relay workers with a global first-production-block gate."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import json
import logging
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.jobs import launch_job  # noqa: E402
from nyl12_covalent.neb_umbrella_campaign import write_json  # noqa: E402
from nyl12_covalent.relay_umbrella_runtime import run_segment, verify_plan  # noqa: E402


def event(folder: Path, action: str, **details: object) -> None:
    """Append durable scheduler intent/outcome records."""
    with (folder / "events.jsonl").open("a") as stream:
        stream.write(
            json.dumps(
                dict(
                    utc=datetime.now(timezone.utc).isoformat(), action=action, **details
                )
            )
            + "\n"
        )
        stream.flush()
        os.fsync(stream.fileno())


def claim_task(path: Path, index: int) -> tuple[str, int] | str | None:
    """Claim one least-sampled window block under an exclusive host queue lock."""
    plan, cfg, base = verify_plan(path)
    folder = path.parent
    claims = folder / "claims"
    claims.mkdir(exist_ok=True)
    with (folder / "queue.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        counts = {}
        busy = set()
        for row in cfg["windows"]:
            label = row["label"]
            status = base / label / "window_status.json"
            if status.exists() and json.loads(status.read_text())["status"] == "failed":
                raise RuntimeError(f"Campaign peer failed: {label}")
            count = 0
            for limit in range(1, plan["production_blocks"] + 1):
                receipt = base / label / f"segment_receipt_limit_{limit:02d}.json"
                if receipt.exists():
                    data = json.loads(receipt.read_text())
                    if (
                        data["production_blocks"] != limit
                        or data["status"] != "accepted"
                        or limit != count + 1
                    ):
                        raise ValueError("Invalid non-prefix block receipt")
                    count = limit
            counts[label] = count
        for claim in claims.glob("*.json"):
            data = json.loads(claim.read_text())
            if data["status"] == "failed":
                raise RuntimeError(f"Failed queue claim requires review: {claim}")
            if data["status"] == "running":
                try:
                    os.kill(data["pid"], 0)
                except ProcessLookupError as error:
                    raise RuntimeError(
                        f"Unresolved stale queue claim: {claim}"
                    ) from error
                busy.add(data["window"])
        minimum = min(counts.values())
        if minimum == plan["production_blocks"]:
            return None
        if minimum >= 1:
            for row in cfg["windows"]:
                first = next(
                    spec for spec in row["stages"] if spec["kind"] == "production"
                )
                outcome = base / Path(first["input"]).parent / "analysis_outcome.json"
                if (
                    not outcome.is_file()
                    or json.loads(outcome.read_text())["status"] != "completed"
                ):
                    event(folder, "extensions_held_for_analysis", window=row["label"])
                    return "analysis_hold"
        candidates = [
            label
            for label, count in counts.items()
            if count == minimum and label not in busy
        ]
        if not candidates:
            return "wait"
        label = candidates[0]
        limit = minimum + 1
        claim = claims / f"{label}_block_{limit:02d}.json"
        with claim.open("x") as stream:
            json.dump(
                dict(
                    status="running",
                    window=label,
                    limit=limit,
                    worker=index,
                    pid=os.getpid(),
                    started_utc=datetime.now(timezone.utc).isoformat(),
                ),
                stream,
            )
            stream.flush()
            os.fsync(stream.fileno())
        event(folder, "block_claimed", window=label, limit=limit, worker=index)
        return label, limit


def finish_claim(path: Path, label: str, limit: int, error: str | None) -> None:
    """Finalize one owned claim durably without changing its scientific receipt."""
    folder = path.parent
    with (folder / "queue.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        claim = folder / "claims" / f"{label}_block_{limit:02d}.json"
        data = json.loads(claim.read_text())
        if data["pid"] != os.getpid() or data["status"] != "running":
            raise ValueError("Claim ownership changed")
        data.update(
            status="failed" if error else "completed",
            error=error,
            ended_utc=datetime.now(timezone.utc).isoformat(),
        )
        write_json(claim, data)
        event(folder, "block_claim_finalized", window=label, limit=limit, error=error)


def worker(path: Path, index: int) -> None:
    """Work a bounded shared queue, prioritizing minimum sampled block count."""
    plan, _, _ = verify_plan(path)
    folder = path.parent / "workers" / f"worker_{index:02d}"
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "claim").open("x") as stream:
        stream.write(str(os.getpid()))
    deadline = time.monotonic() + plan["coverage_gate_timeout_hours"] * 3600
    try:
        while True:
            task = claim_task(path, index)
            if task is None:
                break
            if task == "analysis_hold":
                write_json(
                    folder / "status.json",
                    dict(status="held_for_analysis", sampling_converged=False),
                )
                return
            if task == "wait":
                if time.monotonic() > deadline:
                    raise TimeoutError("Queue integration deadline reached")
                time.sleep(30)
                continue
            label, limit = task
            event(folder, "segment_intent", window=label, limit=limit)
            try:
                run_segment(path, label, limit)
            except Exception as error:
                finish_claim(path, label, limit, f"{type(error).__name__}: {error}")
                raise
            finish_claim(path, label, limit, None)
            event(folder, "segment_accepted", window=label, limit=limit)
        write_json(
            folder / "status.json", dict(status="completed", sampling_converged=False)
        )
    except Exception as error:
        event(folder, "worker_failed", error=f"{type(error).__name__}: {error}")
        write_json(folder / "status.json", dict(status="failed", error=str(error)))
        raise


def main() -> None:
    """Launch an explicitly reviewed worker group or run one worker."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("launch", "worker", "expand"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--index", type=int)
    parser.add_argument("--first", type=int)
    parser.add_argument("--count", type=int)
    args = parser.parse_args()
    plan, cfg, base = verify_plan(args.plan)
    if args.action == "worker":
        if args.index is None or not 0 <= args.index < plan["workers"]:
            parser.error("worker requires a valid index")
        worker(args.plan, args.index)
    else:
        if plan["root_review_authorized"] is not True or not 1 <= plan["workers"] <= 57:
            raise ValueError("Explicit root review and 1..57 workers required")
        folder = args.plan.parent
        if args.action == "expand":
            if (
                args.first is None
                or args.count is None
                or args.count < 1
                or args.first < 0
                or args.first + args.count > plan["workers"]
            ):
                parser.error("expand requires an explicit in-budget first/count range")
            indices = range(args.first, args.first + args.count)
            claim_name = f"expansion_{args.first:02d}_{args.count:02d}_claim"
        else:
            initial = (
                plan["initial_workers"]
                if "initial_workers" in plan
                else plan["workers"]
            )
            if not 1 <= initial <= plan["workers"]:
                raise ValueError("Initial worker count exceeds reviewed budget")
            indices = range(initial)
            claim_name = "launch_claim"
        with (folder / claim_name).open("x") as stream:
            stream.write(str(os.getpid()))
        for index in indices:
            event(folder, "worker_launch_intent", index=index)
            launch_job(
                f"{plan['session_prefix']}_{index:02d}",
                folder / "workers" / f"worker_{index:02d}" / "job",
                base,
                [
                    plan["python"],
                    str(Path(__file__).resolve()),
                    "worker",
                    "--plan",
                    str(args.plan),
                    "--index",
                    str(index),
                ],
                env_overrides=plan["environment"],
                input_paths=[args.plan],
            )
            event(folder, "worker_launch_outcome", index=index)
        write_json(
            folder / "status.json",
            dict(
                status="workers_launched",
                workers=plan["workers"],
                window_count=len(cfg["windows"]),
                sampling_converged=False,
            ),
        )


if __name__ == "__main__":
    main()
