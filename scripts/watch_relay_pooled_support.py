"""Bounded supplemental94state diagnostic monitor; primary aliases are never written."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from nyl12_covalent.relay_pooled_analysis import analyze
from nyl12_covalent.relay_pooled_sources import BASE, NEW, OLD, pooled_data
from nyl12_covalent.relay_umbrella_analysis_io import digest

ROOT = Path("/data/jpc/covalent-enzyme-design")


def journal(phase: str, action: str) -> None:
    """Record exact derived checkpoint intent/outcome; logging failure stops analysis."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_pooled_watcher",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "PIauthorizedsupplemental94statepoint/conditionalCIchecks;4CPU;leadreviewbeforeprimary",
            "--rollback",
            "Retainpooledresultsandfrozenevidence;stoponlyderivedmonitorifneeded",
        ],
        check=True,
    )


def watch(output: Path, hours: float, poll: float) -> None:
    """Wait for exact common10/15/20/25/30ps supplemental checkpoints.

    Args:
        output: Dedicated nonprimary diagnostic namespace.
        hours: Explicit positive runtime bound, atmost72hours.
        poll: Atleast60seconds between readiness checks.
    """
    if not 0 < hours <= 72 or poll < 60:
        raise ValueError("Requirebounded72hoursand>=60secondpoll")
    source_paths = [NEW, OLD, OLD.parent / "extension_20261127/manifest.json"]
    identity = {str(p): digest(p) for p in source_paths}
    status = dict(
        status="watching",
        horizons=[10, 15, 20, 25, 30],
        completed=[],
        independent_states=94,
        main_horizon_ps=45,
        discard_ps_once_per_lineage=5,
        workers=4,
        replicates_per_block_width=200,
        primary_replacement=False,
        source_sha256=identity,
    )
    deadline = time.monotonic() + hours * 3600
    for horizon in status["horizons"]:
        checkpoint = output / f"pooled94_supplement{horizon:02d}ps"
        if (checkpoint / "acceptance.json").exists():
            prior = json.loads((checkpoint / "acceptance.json").read_text())
            if prior["states"] != 94 or prior["horizon_ps"] != horizon:
                raise ValueError("Existingpooledcheckpointidentitymismatch")
            status["completed"].append(horizon)
            continue
        while True:
            if time.monotonic() > deadline:
                raise TimeoutError("Pooledmonitorboundeddeadline")
            for p in source_paths:
                if digest(p) != identity[str(p)]:
                    raise ValueError("Pinnedcohortmanifestchanged")
            status["waiting_for_common_supplemental_ps"] = horizon
            (output / "monitor_state.json").write_text(
                json.dumps(status, indent=2) + "\n"
            )
            if not (BASE / ".compact_sync_in_progress").exists():
                data = pooled_data(output / "baseline_cache", horizon, True)
                if data is not None:
                    break
            time.sleep(poll)
        journal(
            "intent",
            f"Analyze94independentstatesmain45supplement{horizon}ps;5psdiscardonce;4CPU;10pointonly15pluspairedCI",
        )
        analyze(data, checkpoint, workers=4)
        status["completed"].append(horizon)
        journal("outcome", f"Completednonprimarypooled94supplement{horizon}ps")
    status.update(status="complete", waiting_for_common_supplemental_ps=None)
    (output / "monitor_state.json").write_text(json.dumps(status, indent=2) + "\n")


def main() -> None:
    """Run only the user-authorized isolated diagnostic monitor."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hours", type=float, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--poll-seconds", type=float, default=60)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    pinned = json.loads(args.source_manifest.read_text())
    for relative, expected in pinned["files"].items():
        if digest(args.source_manifest.parent / relative) != expected:
            raise ValueError("Pinnedanalysiscodechanged:" + relative)
    for filename, expected in pinned["external_inputs"].items():
        if digest(Path(filename)) != expected:
            raise ValueError("Pinnedanalysisinputchanged:" + filename)
    journal("intent", "Startbounded94statepooleddiagnosticmonitor")
    watch(args.output, args.hours, args.poll_seconds)
    journal("outcome", "Pooleddiagnosticmonitorcompletedthrough30ps")


if __name__ == "__main__":
    main()
