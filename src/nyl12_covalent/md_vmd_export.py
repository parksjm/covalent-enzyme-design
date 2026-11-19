"""Create compact, provenance-rich VMD bundles from Amber MD trajectories."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import subprocess
from typing import Any

from nyl12_covalent.path_inspection import dcd_counts
from nyl12_covalent.umbrella_vmd import netcdf_frame_count

LOGGER = logging.getLogger(__name__)
STRIP_MASK = ":WAT,Na+,Cl-,NA,CL"
SUBSTRATE_MASK = ":1313-1316"


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest for an artifact.

    Args:
        path: File to hash.

    Returns:
        Lowercase hexadecimal digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cpptraj_trajectory_input(
    topology: Path,
    trajectories: list[Path],
    frame_counts: list[int],
    stride: int,
    output: Path,
) -> str:
    """Build an ordered, imaged, stripped cpptraj export input.

    Args:
        topology: Amber topology shared by every trajectory.
        trajectories: Production segments in chronological order.
        frame_counts: Available frames corresponding to ``trajectories``.
        stride: Positive source-frame stride.
        output: Destination DCD path.

    Returns:
        Complete cpptraj input text.
    """
    if not trajectories or len(trajectories) != len(frame_counts):
        raise ValueError("Trajectory and frame-count inventories must match")
    if stride < 1 or any(count < 1 for count in frame_counts):
        raise ValueError("Frame counts and stride must be positive")
    lines = [f"parm {topology.resolve()}"]
    lines.extend(
        f"trajin {path.resolve()} 1 {count} {stride}"
        for path, count in zip(trajectories, frame_counts, strict=True)
    )
    final_frame = frame_counts[-1]
    if (final_frame - 1) % stride != 0:
        lines.append(
            f"trajin {trajectories[-1].resolve()} {final_frame} {final_frame}"
        )
    lines.extend(
        [
            f"autoimage anchor {SUBSTRATE_MASK}",
            f"strip {STRIP_MASK}",
            f"trajout {output.resolve()} dcd",
            "run",
            "quit",
        ]
    )
    return "\n".join(lines) + "\n"


def cpptraj_topology_input(topology: Path, output: Path) -> str:
    """Build cpptraj input for the topology matching the stripped DCD."""
    return (
        f"parm {topology.resolve()}\n"
        f"parmstrip {STRIP_MASK}\n"
        f"parmwrite out {output.resolve()}\nquit\n"
    )


def vmd_script(frame_count: int) -> str:
    """Return a self-contained VMD loader for a neutral-replica movie."""
    if frame_count < 1:
        raise ValueError("frame_count must be positive")
    return f"""set here [file dirname [file normalize [info script]]]
mol new [file join $here stripped.parm7] type parm7 waitfor all
mol addfile [file join $here production_100ns_stripped.dcd] type dcd waitfor all
set molid [molinfo top]
if {{[molinfo $molid get numframes] != {frame_count}}} {{
    error "Expected {frame_count} frames"
}}
mol delrep 0 $molid

# Retain the full tetramer cartoon but hide it by default.
mol representation NewCartoon 0.3 10.0 4.1 0
mol color ColorID 7
mol selection "protein"
mol material Opaque
mol addrep $molid
mol showrep $molid 0 0

# PA6 L4 uses normal element colors, with only its carbons overlaid in orange.
mol representation Licorice 0.18 12.0 12.0
mol color Element
mol selection "resid 1313 to 1316"
mol material Opaque
mol addrep $molid
mol representation Licorice 0.22 12.0 12.0
mol color ColorID 3
mol selection "resid 1313 to 1316 and element C"
mol material Opaque
mol addrep $molid

# Display the changing local protein environment; update it every frame.
mol representation Licorice 0.14 10.0 10.0
mol color Element
mol selection "protein and same residue as within 5.0 of resid 1313 to 1316"
mol material Opaque
mol addrep $molid
mol selupdate 3 $molid on

display projection Orthographic
color Display Background black
animate goto 0
puts "Loaded neutral-replica 100 ns movie: {frame_count} frames"
puts "Representation 0 is the full protein cartoon and is hidden by default."
"""


def _run_cpptraj(executable: Path, input_path: Path, log_path: Path) -> None:
    """Run cpptraj with one retained input and combined log."""
    with log_path.open("w") as log:
        subprocess.run(
            [str(executable), "-i", str(input_path)],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def _write_checksum_manifest(output: Path, names: list[str]) -> None:
    """Write checksums for the declared portable bundle members."""
    (output / "MANIFEST.sha256").write_text(
        "".join(f"{sha256(output / name)}  {name}\n" for name in names)
    )


def export_replica(
    replica: Path,
    cpptraj: Path,
    output: Path,
    stride: int = 25,
) -> dict[str, Any]:
    """Export and validate one neutral-replica 100 ns VMD bundle.

    Args:
        replica: Directory containing the matching Amber topology and segments.
        cpptraj: AmberTools cpptraj executable.
        output: New bundle directory; it must not already exist.
        stride: Source-frame stride. Production frames are separated by 10 ps.

    Returns:
        Serializable export manifest.
    """
    topology = replica / "system.prmtop"
    trajectories = [
        replica / "05_production_10ns.nc",
        replica / "06_production_90ns.nc",
    ]
    for path in [topology, cpptraj, *trajectories]:
        if not path.is_file():
            raise FileNotFoundError(path)
    output.mkdir(parents=True, exist_ok=False)
    frame_counts = [netcdf_frame_count(path) for path in trajectories]
    if frame_counts != [1000, 9000]:
        raise ValueError(f"Expected 1000 + 9000 production frames: {frame_counts}")
    trajectory = output / "production_100ns_stripped.dcd"
    trajectory_input = output / "trajectory.cpptraj"
    topology_input = output / "topology.cpptraj"
    trajectory_input.write_text(
        cpptraj_trajectory_input(
            topology, trajectories, frame_counts, stride, trajectory
        )
    )
    topology_input.write_text(
        cpptraj_topology_input(topology, output / "stripped.parm7")
    )
    _run_cpptraj(cpptraj, trajectory_input, output / "trajectory_export.log")
    _run_cpptraj(cpptraj, topology_input, output / "topology_export.log")
    frames, atoms = dcd_counts(trajectory)
    expected_frames = sum((count - 1) // stride + 1 for count in frame_counts)
    if (frame_counts[-1] - 1) % stride != 0:
        expected_frames += 1
    if frames != expected_frames:
        raise ValueError(f"Expected {expected_frames} frames, found {frames}")
    viewer = output / "view_production.vmd"
    viewer.write_text(vmd_script(frames))
    sources = [
        {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in trajectories
    ]
    manifest: dict[str, Any] = {
        "replica": replica.name,
        "source_topology": str(topology.resolve()),
        "source_topology_sha256": sha256(topology),
        "source_trajectories": sources,
        "source_frame_counts": frame_counts,
        "source_frame_interval_ps": 10.0,
        "stride": stride,
        "display_frame_interval_ps": 10.0 * stride,
        "exact_final_frame_appended": (frame_counts[-1] - 1) % stride != 0,
        "frames": frames,
        "stripped_atoms": atoms,
        "strip_mask": STRIP_MASK,
        "imaging": f"autoimage anchor {SUBSTRATE_MASK}",
        "scientific_scope": "visualization only; downsampled coordinates are not an analysis trajectory",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    (output / "README.md").write_text(
        "# Neutral Nyl12–PA6 production trajectory\n\n"
        "Run `vmd -e view_production.vmd` from this directory. The 100 ns "
        "production trajectory is downsampled to approximately 0.25 ns spacing "
        "for visual inspection. Protein and PA6 L4 are retained; bulk water and "
        "ions are removed. The full cartoon exists but is hidden by default. "
        "Use the original NetCDF files for quantitative analysis.\n"
    )
    members = [
        "README.md",
        "manifest.json",
        "production_100ns_stripped.dcd",
        "stripped.parm7",
        "topology.cpptraj",
        "topology_export.log",
        "trajectory.cpptraj",
        "trajectory_export.log",
        "view_production.vmd",
    ]
    _write_checksum_manifest(output, members)
    LOGGER.info(
        "Exported %s: %d frames, %d atoms, %.1f MiB",
        replica.name,
        frames,
        atoms,
        trajectory.stat().st_size / 1024**2,
    )
    return manifest
