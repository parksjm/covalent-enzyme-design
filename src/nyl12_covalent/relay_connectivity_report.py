"""Human-readable output and bounded watcher for relay connectivity audits."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import time
from typing import Any

from .relay_connectivity_audit import (
    BRANCHES,
    SnapshotChangedError,
    audit_campaign,
    digest,
    read_json,
)

LOGGER = logging.getLogger(__name__)


def render_markdown(report: dict[str, Any]) -> str:
    """Render a compact human-readable audit from authoritative JSON."""
    lines = [
        "# N-terminal relay connectivity audit",
        "",
        f"Status: **{report['audit_status']}**. Connectivity: **not established**.",
        "",
        "| Branch | Status | Accepted path frames | Last geometry | C–OG1 / C–Nleave / C–O (Å) | Physical energy (kcal/mol) | Bias (kcal/mol) |",
        "|---|---|---:|---|---:|---:|---:|",
    ]
    for branch, data in report["branches"].items():
        last = data["accepted_path_rows"][-1] if data["accepted_path_rows"] else None
        if last:
            distances = last["distances_A"]
            lines.append(
                f"| {branch} | {data['status']} | {len(data['accepted_path_rows'])} | "
                f"{last['geometry_label']} | {distances['C_OG']:.3f} / "
                f"{distances['C_NL']:.3f} / {distances['C_O']:.3f} | "
                f"{last['physical_kcal_mol']:.4f} | {last['bias_kcal_mol']:.4f} |"
            )
        else:
            lines.append(f"| {branch} | {data['status']} | 0 | none | — | — | — |")
    lines.extend(
        [
            "",
            "Tracked H ownership uses original atom IDs; all QM N/O-bound hydrogens use the existing PBC-aware pool observable.",
            "Distance discontinuities over 0.35 Å are listed in the JSON. Torsions are unavailable in these records.",
            "",
            "Shared physical topology/exterior/Hamiltonian hashes are checked. Branch-specific coordinate biases are reported separately from physical energies.",
            "",
            report["required_next_gate"],
            "",
        ]
    )
    return "\n".join(lines)


def refresh_watch(
    campaign: Path,
    output: Path,
    interval_s: float = 30.0,
    max_seconds: float = 43200.0,
    once: bool = False,
) -> None:
    """Refresh only audit outputs when progress or status receipts change.

    Args:
        campaign: Existing campaign directory.
        output: Dedicated connectivity_audit output directory.
        interval_s: Positive polling interval in seconds.
        max_seconds: Positive maximum watch duration.
        once: Run one snapshot and return.
    """
    if interval_s <= 0 or max_seconds <= 0:
        raise ValueError("Watch interval and maximum duration must be positive")
    output.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + max_seconds
    prior_signature = None
    while True:
        watched = [
            campaign / branch / name
            for branch in BRANCHES
            for name in ("progress.json", "branch_status.json")
        ]
        signature = tuple((str(path), digest(path)) for path in watched)
        if signature != prior_signature:
            try:
                report = audit_campaign(campaign)
            except SnapshotChangedError:
                LOGGER.warning("Campaign changed during snapshot; deferring refresh")
                if once:
                    return
                time.sleep(interval_s)
                continue
            json_pending = output / "connectivity_audit.json.pending"
            markdown_pending = output / "connectivity_audit.md.pending"
            json_pending.write_text(
                json.dumps(report, indent=2, allow_nan=False) + "\n"
            )
            markdown_pending.write_text(render_markdown(report))
            json_pending.replace(output / "connectivity_audit.json")
            markdown_pending.replace(output / "connectivity_audit.md")
            prior_signature = signature
            LOGGER.info("Refreshed connectivity audit in %s", output)
        if once or time.monotonic() >= deadline:
            return
        statuses = [
            read_json(campaign / branch / "branch_status.json").get("status")
            for branch in BRANCHES
        ]
        if all(status not in {"running", "queued", "pending"} for status in statuses):
            if signature == prior_signature:
                return
        time.sleep(min(interval_s, max(0.0, deadline - time.monotonic())))
