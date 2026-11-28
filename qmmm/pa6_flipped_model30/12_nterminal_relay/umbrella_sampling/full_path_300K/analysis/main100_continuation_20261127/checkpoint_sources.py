"""Source-bind complete manifest-window checkpoints without preliminary fits."""

from __future__ import annotations

import json
from pathlib import Path

from nyl12_covalent.relay_mbar_progression import select_horizon
from nyl12_covalent.relay_umbrella_analysis_io import digest, load_block
from watch_relay_umbrella_analysis import discover


def bind_checkpoint(base: Path, output: Path, horizon: int) -> Path | None:
    """Bind an exact complete continuous-production horizon to accepted sources.

    Args:
        base: Existing native campaign root; only read.
        output: Derived monitor directory, outside legacy checkpoint discovery.
        horizon: One authorized future horizon,65 through100 in5ps steps.

    Returns:
        Immutable checkpoint directory, orNone when common coverage is incomplete.
    """
    if horizon not in range(65, 101, 5):
        raise ValueError("Only authorized65through100ps checkpoints accepted")
    previous = base / "analysis/fixed_burnin_20261126/checkpoint_35ps_refresh_20261126"
    config = json.loads((previous / "analysis_request.json").read_text())
    receipts, _ = discover(config, base, "analysis_acceptance.json")
    selected = select_horizon(receipts, config["windows"], horizon)
    if selected is None:
        return None
    paths, durations = selected
    config["accepted_block_receipts"] = [str(p.resolve()) for p in paths]
    blocks = [
        load_block(
            p,
            temperature=config["temperature_K"],
            hamiltonian=config["physical_hamiltonian_sha256"],
            diagnostic_columns=config["diagnostic_columns"],
        )
        for p in paths
    ]
    checkpoint = output / f"checkpoint_{horizon:02d}ps"
    request = checkpoint / "analysis_request.json"
    if checkpoint.exists():
        if json.loads(request.read_text()) != config:
            raise ValueError("Existingcheckpointsourceconfigurationchanged")
        parent = json.loads((checkpoint / "analysis.json").read_text())
        if parent["request_sha256"] != digest(request):
            raise ValueError("Existingcheckpointrequesthashchanged")
        return checkpoint
    checkpoint.mkdir()
    request.write_text(json.dumps(config, indent=2) + "\n")
    prior = json.loads((previous / "analysis.json").read_text())
    parent = dict(
        request_sha256=digest(request),
        accepted_sources=[
            dict(receipt_sha256=b["receipt_sha256"], source_sha256=b["source_sha256"])
            for b in blocks
        ],
        path_3cv_sha256=prior["path_3cv_sha256"],
        checkpoint=dict(
            prior["checkpoint"],
            actual_horizon_ps=horizon,
            scope="Exactmanifest57windowcontinuousproduction; sourcebindingonly",
        ),
        window_durations_ps=durations,
        source_sha256={
            str(p): digest(p)
            for p in [
                Path(__file__),
                previous / "analysis.json",
                previous / "analysis_request.json",
            ]
        },
        preliminary_MBAR_fit=False,
    )
    (checkpoint / "analysis.json").write_text(json.dumps(parent, indent=2) + "\n")
    return checkpoint
