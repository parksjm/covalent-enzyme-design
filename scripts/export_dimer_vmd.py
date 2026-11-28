"""Export saved dimer search states as an immutable stripped VMD bundle."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from export_direct_pt_progress import (  # noqa: E402
    QM_INDICES,
    STRIP_MASK,
    run_cpptraj,
    sha256,
    vmd_script,
    write_checksum_manifest,
)
from nyl12_covalent.path_inspection import dcd_counts  # noqa: E402
from scipy.io import netcdf_file  # noqa: E402


from nyl12_covalent.vmd_style_refresh import append_file_policy  # noqa: E402


def export(base: Path) -> Path:
    """Snapshot the native trajectory and export its saved search states.

    Args:
        base: Prepared native dimer campaign directory.

    Returns:
        New immutable visualization bundle directory.
    """
    base = base.resolve(strict=True)
    source = base / "calculation/optimization.nc"
    output = (
        base
        / "vmd_progress"
        / datetime.now(timezone.utc).strftime("update_%Y%m%dT%H%M%S")
    )
    output.mkdir(parents=True, exist_ok=False)
    snapshot = output / "source_snapshot"
    snapshot.mkdir()
    shutil.copy2(source, snapshot / source.name)
    shutil.copy2(base / "calculation/dimer.mdout", snapshot / "dimer.mdout")
    shutil.copy2(Path(__file__), snapshot / Path(__file__).name)
    with netcdf_file(snapshot / source.name, "r", mmap=True) as nc:
        shape = nc.variables["coordinates"].shape
        saved_frames, original_atoms = int(shape[0]), int(shape[1])
    if saved_frames < 1 or original_atoms != 104922:
        raise ValueError(f"Unexpected native trajectory shape: {shape}")
    topology = base / "inputs/system.prmtop"
    initial = base / "inputs/ts_guess.ncrst"
    cpptraj = Path("/data/jpc/code/ambertools26/venv/bin/cpptraj")
    trajectory_input = output / "trajectory.cpptraj"
    trajectory_input.write_text(
        f"parm {topology}\ntrajin {initial}\n"
        f"trajin {snapshot / source.name} 1 {saved_frames}\n"
        f"strip {STRIP_MASK}\ntrajout {output / 'ts_optimization.dcd'} dcd\nrun\nquit\n"
    )
    topology_input = output / "topology.cpptraj"
    topology_input.write_text(
        f"parm {topology}\nparmstrip {STRIP_MASK}\n"
        f"parmwrite out {output / 'stripped.parm7'}\nquit\n"
    )
    for script in (trajectory_input, topology_input):
        run_cpptraj(cpptraj, script, output, output / f"{script.stem}_export.log")
    frames, atoms = dcd_counts(output / "ts_optimization.dcd")
    if frames != saved_frames + 1 or atoms != 19967:
        raise ValueError(f"Unexpected stripped inventory: {frames}, {atoms}")
    (output / "view_ts_optimization.vmd").write_text(
        vmd_script(
            frames,
            trajectory_name="ts_optimization.dcd",
            description="saved dimer search states",
        )
    )
    append_file_policy(
        output / "view_ts_optimization.vmd",
        output / "stripped.parm7",
        [i - 1 for i in QM_INDICES],
    )
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "source": str(source),
                "snapshot_sha256": sha256(snapshot / source.name),
                "topology_sha256": sha256(topology),
                "initial_sha256": sha256(initial),
                "frames": frames,
                "stripped_atoms": atoms,
                "strip_mask": STRIP_MASK,
                "frame_zero": "Exact dimer starting guess",
                "remaining_frames": "Native saved optimizer states in file order; may include dimer-displaced states; not MD time or a converged reaction pathway",
                "native_output_snapshot_sha256": sha256(snapshot / "dimer.mdout"),
                "stationary_point_validation": "Not established by this visualization",
            },
            indent=2,
        )
        + "\n"
    )
    (output / "README.md").write_text(
        "Load view_ts_optimization.vmd with VMD. Frame0 is the exact starting "
        "guess; later frames are saved native dimer search states, not MD or a "
        "validated reaction path. Protein/substrate retained; solvent/ions stripped.\n"
    )
    write_checksum_manifest(output)
    latest = output.parent / "latest"
    pending = output.parent / ".latest.pending"
    pending.symlink_to(output.name)
    pending.replace(latest)
    logging.info("Exported %d frames, %d atoms to %s", frames, atoms, output)
    return output


def main() -> None:
    """Parse the campaign path and create a visualization checkpoint."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    export(args.base)


if __name__ == "__main__":
    main()
