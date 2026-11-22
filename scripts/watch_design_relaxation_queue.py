#!/usr/bin/env python3
"""Consume explicit, hashed model requests on one assigned Amber GPU."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.design_model_relaxation import relax_model  # noqa: E402


def _publish_status(path: Path, text: str) -> None:
    """Publish a terminal/readiness receipt without exposing a truncated JSON file."""
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(text)
    os.replace(temp, path)


def main() -> None:
    """Wait for disjoint requests, run minimization once and retain terminal state."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--gpu", choices=["0", "1"], required=True)
    args = parser.parse_args()
    base = args.campaign.resolve()
    plan = json.loads((base / "prepared/candidate_manifest.json").read_text())
    config = json.loads((base / "configuration/relaxation.json").read_text())
    config["gpu"] = args.gpu
    jobs = [j for i, j in enumerate(plan["jobs"]) if i % 2 == int(args.gpu)]
    deadline = time.monotonic() + 43200
    while jobs and time.monotonic() < deadline:
        pending = []
        for job in jobs:
            case = (
                base / "models" / job["variant"] / f"seed_{job['model_index'] + 1:02d}"
            )
            request = case / "request.json"
            if not request.exists():
                pending.append(job)
                continue
            payload = json.loads(request.read_text())
            source = Path(payload["source"])
            if hashlib.sha256(source.read_bytes()).hexdigest() != payload["sha256"]:
                raise ValueError(f"Normalized source changed: {source}")
            status = case / "relaxation_worker_status.json"
            if status.exists():
                old = json.loads(status.read_text())
                if old["status"] != "completed":
                    raise RuntimeError(f"Explicit restart review required: {status}")
                continue
            _publish_status(
                status,
                json.dumps(
                    {
                        "status": "running",
                        "gpu": args.gpu,
                        "source_sha256": payload["sha256"],
                    },
                    indent=2,
                )
                + "\n",
            )
            try:
                result = relax_model(source, case / "relaxation", config)
            except Exception as exc:
                _publish_status(
                    status,
                    json.dumps(
                        {
                            "status": "failed",
                            "gpu": args.gpu,
                            "error": f"{type(exc).__name__}: {exc}",
                        },
                        indent=2,
                    )
                    + "\n",
                )
                raise
            _publish_status(
                status,
                json.dumps(
                    {"status": "completed", "gpu": args.gpu, "result": result}, indent=2
                )
                + "\n",
            )
            logging.info(
                "Completed %s seed %s on GPU%s",
                job["variant"],
                job["model_index"] + 1,
                args.gpu,
            )
        jobs = pending
        if jobs:
            time.sleep(30)
    if jobs:
        raise TimeoutError(f"{len(jobs)} requests did not complete within12h")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
