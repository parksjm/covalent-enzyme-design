"""Bounded automatic65–100ps fixed-protocol PMF and movie checkpoint updates."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/data/jpc/covalent-enzyme-design")


def journal(phase: str, action: str) -> None:
    """Record derived-only operations and fail if the append-only journal fails."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_auto_checkpoints",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "PIstandingautomatic5psPMFandmovieupdates;4CPU;no35refit;nativejobsunchanged",
            "--rollback",
            "Stoponlyderivedmonitor;preservesourcesresultsandarchivedmovies",
        ],
        check=True,
    )


def run(base: Path, output: Path, source: Path, hours: float) -> None:
    """Publish exact common future horizons once, then terminate after100ps.

    Args:
        base: Existing native physical campaign, read-only to this monitor.
        output: New derived monitor and source-bound checkpoint directory.
        source: Pinned analysis/movie/discovery sources.
        hours: Positive finite monitoring time bound.
    """
    from checkpoint_sources import bind_checkpoint
    from movie_update import sha
    from watch_relay_fixed_pmf import publish

    manifest = json.loads((output / "source_manifest.json").read_text())
    for name, expected in manifest.items():
        if sha(output / name) != expected:
            raise ValueError("Pinnedautomaticmonitorsourcechanged:" + name)
    baseline = base / "analysis/fixed_burnin_20261126"
    evidence = output / "evidence"
    evidence.mkdir(exist_ok=True)
    (output / "plots").mkdir(exist_ok=True)
    initial = json.loads((baseline / "fixed_35ps_publication.json").read_text())
    if initial["horizon_ps"] != 35:
        raise ValueError("Required35pspublicationmissing")
    state = dict(
        status="watching",
        horizons=list(range(65, 101, 5)),
        workers=4,
        bootstrap_replicates_per_width=200,
        CI="unchangedpointwiseconditional95%",
        start_after_completed35ps=True,
        native_jobs_changed=False,
        completed=[],
    )
    deadline = time.monotonic() + hours * 3600
    for horizon in state["horizons"]:
        while True:
            if time.monotonic() > deadline:
                raise TimeoutError("Boundedautomaticcheckpointmonitorreacheddeadline")
            state.update(waiting_for_common_horizon_ps=horizon)
            pending = output / "monitor_state.pending"
            pending.write_text(json.dumps(state, indent=2) + "\n")
            pending.replace(output / "monitor_state.json")
            if (base / ".compact_sync_in_progress").exists():
                time.sleep(60)
                continue
            checkpoint = bind_checkpoint(base, output, horizon)
            if checkpoint is not None:
                break
            time.sleep(60)
        journal(
            "intent", f"Publishfixed{horizon}ps200bootstrap4workersanalysis_only"
        )
        publish(checkpoint, evidence, output / "plots", source, ROOT, workers=4)
        # Main60ps movie remains immutable; this sidecar publishes analysis only.
        state["completed"].append(horizon)
        journal(
            "outcome",
            f"Publishedfixed{horizon}psPMF95percentconditionalCIandvalidatedmovie",
        )
    state.update(status="complete", waiting_for_common_horizon_ps=None)
    (output / "monitor_state.json").write_text(json.dumps(state, indent=2) + "\n")


def main() -> None:
    """Use pinned discovery/analysis/movie code and a finite72hourmonitor."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--hours", type=float, required=True)
    args = parser.parse_args()
    if not 0 < args.hours <= 72:
        raise ValueError("Explicitpositivefiniteboundatmost72hoursrequired")
    sys.path.insert(0, str(args.source / "scripts"))
    sys.path.insert(0, str(args.source / "monitor"))
    journal("intent", "Startpinnedautomatic65through100psderivedcheckpointmonitor")
    run(args.base, args.output, args.source, args.hours)
    journal("outcome", "Automaticderivedcheckpointmonitorcompletedthrough100ps")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
