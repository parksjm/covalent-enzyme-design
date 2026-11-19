"""Wait for grid and reference gates, then prepare or run bounded releases."""

from __future__ import annotations

import argparse
import logging
import os

from release_execution import (
    SOURCE_MDIN,
    TOPOLOGY,
    journal,
    prepare_releases,
    run_release,
    write_status,
)
from release_queue_support import (
    GRID,
    PROTOCOL,
    QUEUE_GATE,
    STATUS_PATH,
    digest,
    load_json,
    select_candidates,
    wait_until_ready,
)

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Wait at most 48 hours and optionally run gated releases in root tmux."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--wait", action="store_true", required=True)
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run sequential releases in root-launched tmux",
    )
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.run and not os.environ.get("TMUX"):
        raise RuntimeError(
            "--run is allowed only inside the root-launched tmux session"
        )
    if STATUS_PATH.exists():
        raise FileExistsError(f"Queue controller has already run: {STATUS_PATH}")
    protocol = load_json(PROTOCOL)
    expected_mdin_hash = protocol["source_sha256"].get(str(SOURCE_MDIN))
    expected_topology_hash = protocol["source_sha256"].get(str(TOPOLOGY))
    if (
        digest(SOURCE_MDIN) != expected_mdin_hash
        or digest(TOPOLOGY) != expected_topology_hash
    ):
        raise ValueError(
            "Source unrestrained mdin or topology differs from the frozen grid protocol"
        )
    mdin = SOURCE_MDIN.read_text().replace(" ", "")
    if "nmropt=0" not in mdin or "ntr=0" not in mdin or "qmcharge=-1" not in mdin:
        raise ValueError(
            "Unbiased release source does not retain required QM/MM controls"
        )
    write_status(
        "waiting", "waiting up to 48 h for terminal grid and reference comparison"
    )
    ready, reason, evidence = wait_until_ready()
    if not ready:
        write_status("pending", reason, evidence=evidence)
        journal("outcome", "Wait for grid/reference release gate", f"pending: {reason}")
        LOGGER.warning("Candidate releases remain pending: %s", reason)
        return
    candidates, decisions = select_candidates()
    journal(
        "intent",
        "Prepare up to two unique accepted candidate release inputs after root-reviewed gate",
    )
    release_dirs = prepare_releases(candidates, decisions, evidence["gate"])
    journal(
        "outcome",
        "Prepare up to two unique accepted candidate release inputs after root-reviewed gate",
        f"prepared {len(release_dirs)} candidate(s); classes={decisions['classes']}",
    )
    if not args.run:
        LOGGER.info(
            "Prepared %d release candidate(s); no engine launched", len(release_dirs)
        )
        return
    write_status(
        "running",
        "root-launched sequential unbiased releases",
        selected=decisions,
        reference_gate_sha256=digest(QUEUE_GATE),
        releases=[str(p.relative_to(GRID)) for p in release_dirs],
    )
    outcomes = []
    for folder, candidate in zip(release_dirs, candidates, strict=True):
        outcomes.append(run_release(folder, candidate))
    write_status(
        "completed",
        "all prepared releases reached terminal classification",
        selected=decisions,
        outcomes=outcomes,
        releases=[str(p.relative_to(GRID)) for p in release_dirs],
    )


if __name__ == "__main__":
    main()
