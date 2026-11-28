"""Bounded isolated106-history MBAR checkpoints against frozen70ps main100."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
from boundary_collect import WORK, registration, ready, collect_boundary
from diagnostics import summarize
from nyl12_covalent.relay_pooled_estimator import point
from nyl12_covalent.relay_pooled_analysis import write
from nyl12_covalent.relay_umbrella_analysis_io import digest

ROOT = Path("/data/jpc/covalent-enzyme-design")


def journal(phase: str, action: str) -> None:
    """Fail closed on unavailable append-only action logging."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "boundary_mbar_queue",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Approved isolated106-history fixed-protocol MBAR checkpoint",
            "--rollback",
            "Preserve all attempts; stop only this derived monitor",
        ],
        check=True,
    )


def verify() -> None:
    """Validate all pinned sources, prepared manifests and frozen baseline files."""
    for path, expected in json.loads(
        (WORK / "source_manifest.json").read_text()
    ).items():
        if digest(Path(path)) != expected:
            raise ValueError("Pinned input changed: " + path)


def publish(horizon: int) -> None:
    """Preserve attempted fit and expose no current symlink before acceptance."""
    folder = WORK / f"checkpoint_boundary{horizon:02d}ps_frozen_main70ps"
    folder.mkdir()
    data = collect_boundary(horizon, folder)
    report, _, _, lw = point(data)
    write(folder / "point.json", report)
    write(
        folder / "sources.json",
        {k: v for k, v in data.items() if k not in ("q", "counts")},
    )
    np.savez_compressed(
        folder / "weights.npz",
        q=data["q"],
        counts=data["counts"],
        logweights=lw,
        f_k=report["f_k"],
    )
    diagnostics = summarize(data, lw)
    write(folder / "diagnostics.json", diagnostics)
    baseline = Path(registration()["baseline"])
    prior = json.loads((baseline / "point.json").read_text())
    old = np.array(prior["projections"]["hard"]["pmf"], float)
    new = np.array(report["projections"]["hard"]["pmf"], float)
    write(
        folder / "comparison.json",
        dict(
            frozen_baseline=str(baseline),
            frozen_main_horizon_ps=70,
            old_histories=100,
            new_histories=106,
            unique_bias_potentials=95,
            new_production_horizon_ps=horizon,
            PMF_difference_kcal_mol=[
                float(b - a) if np.isfinite(a + b) else None
                for a, b in zip(old, new, strict=True)
            ],
            baseline_support=prior["projections"]["hard"]["ESS"],
            checkpoint_support=report["projections"]["hard"]["ESS"],
            baseline_physical_diagnostics=str(WORK / "baseline_diagnostics.json"),
            CI_status="Point only; uncertainty unresolved",
            retained_5ps_chunks_per_new_history=(horizon - 5) // 5,
            pooled_CI_qualified=False,
            mixing_diagnostic_only=True,
            sampling_converged=False,
            primary_replacement=False,
        ),
    )
    write(
        folder / "acceptance.json",
        dict(
            process_complete=True,
            native_parity_validated=True,
            direct_CV_identities_verified=3,
            exact_weight_order=True,
            connected=report["connected"],
            accepted=bool(report["connected"]),
            scientific_convergence_accepted=False,
            primary_replacement=False,
        ),
    )
    write(
        folder / "payload_sha256.json",
        {p.name: digest(p) for p in folder.iterdir() if p.is_file()},
    )


def main() -> None:
    """Wait for common10ps then15ps and5ps checkpoints without changing active jobs."""
    config = registration()
    deadline = time.monotonic() + config["timeout_hours"] * 3600
    completed = []
    for horizon in config["minimum_new_production_checkpoints_ps"]:
        while True:
            journal(
                "intent",
                f"Verify frozen sources and inspect six accepted prefixes for{horizon}ps",
            )
            verify()
            if ready(horizon):
                journal("outcome", f"All six accepted{horizon}ps prefixes available")
                break
            write(
                WORK / "monitor_state.json",
                dict(
                    status="waiting_for_accepted_boundary_stages",
                    waiting_for_common_new_production_ps=horizon,
                    frozen_baseline_main_ps=70,
                    completed=completed,
                    poll_seconds=config["poll_seconds"],
                    native_changes=False,
                ),
            )
            journal(
                "outcome",
                f"Waiting normally for six common{horizon}ps production histories",
            )
            if time.monotonic() > deadline:
                raise TimeoutError("Authorized72h derived monitor bound reached")
            time.sleep(config["poll_seconds"])
        journal(
            "intent",
            f"Fit isolated106-history global MBAR at new{horizon}ps; frozen70ps main100",
        )
        publish(horizon)
        completed.append(horizon)
        journal(
            "outcome",
            f"Retained isolated{horizon}ps checkpoint; no primary/current symlink changes",
        )
    write(
        WORK / "monitor_state.json",
        dict(status="complete", completed=completed, primary_replacement=False),
    )


if __name__ == "__main__":
    main()
