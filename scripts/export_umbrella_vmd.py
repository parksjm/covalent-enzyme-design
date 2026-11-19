#!/usr/bin/env python3
"""Export sequential, downsampled umbrella-window frames for VMD inspection."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.path_inspection import dcd_counts  # noqa: E402
from nyl12_covalent.umbrella_vmd import (  # noqa: E402
    campaign_frame_plan,
    cpptraj_trajectory_input,
)
from scripts.export_direct_pt_progress import (  # noqa: E402
    STRIP_MASK,
    run_cpptraj,
    sha256,
    vmd_script,
    write_checksum_manifest,
)

LOGGER = logging.getLogger(__name__)


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def export(
    campaign: Path,
    topology: Path,
    cpptraj: Path,
    output: Path,
    maximum_frames_per_window: int,
    allow_incomplete: bool = False,
) -> dict[str, Any]:
    """Create and verify one immutable umbrella visualization bundle.

    Args:
        campaign: Umbrella campaign directory.
        topology: Full Amber topology matching every trajectory.
        cpptraj: AmberTools cpptraj executable.
        output: New output directory.
        maximum_frames_per_window: Independent downsampling ceiling per window.
        allow_incomplete: Skip windows without completed production blocks.

    Returns:
        Bundle manifest.
    """
    for path in (campaign / "inputs.json", topology, cpptraj):
        if not path.is_file():
            raise FileNotFoundError(path)
    output.mkdir(parents=True, exist_ok=False)
    plan = campaign_frame_plan(campaign, maximum_frames_per_window, allow_incomplete)
    trajectory = output / "umbrella_windows.dcd"
    trajectory_input = output / "trajectory.cpptraj"
    trajectory_input.write_text(
        cpptraj_trajectory_input(topology, plan, trajectory, STRIP_MASK)
    )
    topology_input = output / "topology.cpptraj"
    topology_input.write_text(
        f"parm {topology.resolve()}\nparmstrip {STRIP_MASK}\n"
        f"parmwrite out {(output / 'stripped.parm7').resolve()}\nquit\n"
    )
    run_cpptraj(
        cpptraj,
        trajectory_input,
        output,
        output / "trajectory_export.log",
    )
    run_cpptraj(
        cpptraj,
        topology_input,
        output,
        output / "topology_export.log",
    )
    frames, atoms = dcd_counts(trajectory)
    if frames != len(plan):
        raise ValueError("Exported DCD frame count differs from frame plan")
    viewer = output / "view_umbrella.vmd"
    viewer.write_text(
        vmd_script(
            frames,
            trajectory_name=trajectory.name,
            description="sequential umbrella-sampling",
        )
    )
    serializable_plan = [
        {**row, "source": str(row["source"].resolve())} for row in plan
    ]
    _write_json(output / "frame_manifest.json", serializable_plan)
    sources = sorted({row["source"].resolve() for row in plan})
    window_counts: dict[str, int] = {}
    for row in plan:
        window_counts[row["window"]] = window_counts.get(row["window"], 0) + 1
    manifest = {
        "campaign": str(campaign.resolve()),
        "campaign_inputs_sha256": sha256(campaign / "inputs.json"),
        "topology": str(topology.resolve()),
        "topology_sha256": sha256(topology),
        "maximum_frames_per_window": maximum_frames_per_window,
        "allow_incomplete": allow_incomplete,
        "windows": window_counts,
        "frames": frames,
        "stripped_atoms": atoms,
        "strip_mask": STRIP_MASK,
        "source_trajectories": [
            {"path": str(path), "sha256": sha256(path)} for path in sources
        ],
    }
    _write_json(output / "manifest.json", manifest)
    write_checksum_manifest(output)
    return manifest


def main() -> None:
    """Parse paths and export an umbrella-window VMD bundle."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--maximum-frames-per-window", type=int, default=100)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    manifest = export(
        args.campaign,
        args.topology,
        args.cpptraj,
        args.output,
        args.maximum_frames_per_window,
        args.allow_incomplete,
    )
    LOGGER.info(
        "Published %d frames from %d windows in %s",
        manifest["frames"],
        len(manifest["windows"]),
        args.output,
    )


if __name__ == "__main__":
    main()
