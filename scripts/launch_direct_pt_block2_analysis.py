"""Reuse recorded analysis commands for the complete second production block."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.jobs import launch_job  # noqa: E402


def main() -> None:
    """Launch MBAR and chemistry audits in separate provenance-recorded sessions."""
    campaign = (
        ROOT
        / "qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/07_direct_pt_umbrella_30ps_violet"
    )
    for window in range(20):
        process = (
            campaign / f"window_{window:02d}/02_production_block_02_5ps/process.json"
        )
        record = json.loads(process.read_text())
        if record["exit_status"] != 0:
            raise ValueError(f"Block 2 incomplete: {process}")
    previous = json.loads(
        (campaign / "analysis/mbar_7ps_complete/job/command.json").read_text()
    )
    command = previous["command"]
    index = command.index("--output")
    command[index:index] = ["--source", str(campaign), "02_production_block_02_5ps"]
    command[command.index("--output") + 1] = str(
        campaign / "analysis/mbar_12ps_complete/result"
    )
    command[command.index("--title") + 1] = (
        "Direct HG1 transfer: interim 12 ps/window PMF"
    )
    launch_job(
        "nyl12_direct_pt_mbar_12ps_20261110",
        campaign / "analysis/mbar_12ps_complete/job",
        ROOT,
        command,
        env_overrides={"OPENBLAS_NUM_THREADS": "1"},
        input_paths=[Path(p) for p in previous["input_sha256"]],
    )
    previous = json.loads(
        (
            campaign / "analysis/proton_state_block01_5ps/job_retry/command.json"
        ).read_text()
    )
    command = previous["command"]
    command[command.index("--stage") + 1] = "02_production_block_02_5ps/trajectory.nc"
    command[command.index("--output") + 1] = str(
        campaign / "analysis/proton_state_block02_5ps/result"
    )
    launch_job(
        "nyl12_direct_pt_audit_block02_20261110",
        campaign / "analysis/proton_state_block02_5ps/job",
        Path(previous["cwd"]),
        command,
        env_overrides=previous["env_overrides"],
        input_paths=[Path(p) for p in previous["input_sha256"]],
    )


if __name__ == "__main__":
    main()
