#!/usr/bin/env python3
"""Export finalized direct-proton-transfer windows for immediate VMD inspection."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.path_inspection import dcd_counts  # noqa: E402

LOGGER = logging.getLogger(__name__)
STRIP_MASK = ":WAT,Na+,Cl-,NA,CL"
QM_INDICES = [
    *range(1621, 1656),
    *range(2727, 2735),
    *range(13905, 13932),
    *range(14555, 14561),
    *range(14574, 14580),
    *range(19914, 19918),
    *range(19925, 19932),
    *range(19937, 19942),
]
PHE113_N = 1637
ASN186_ND2 = 2731
THR233_N = 13904
THR233_OG1 = 13915
THR233_HG1 = 13916
PA6_CARBONYL_C = 19915
PA6_CARBONYL_O = 19916
PA6_SCISSILE_N = 19928


def sha256(path: Path) -> str:
    """Return a file's SHA-256 identity."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_checksum_manifest(directory: Path) -> None:
    """Write a portable SHA-256 manifest for every current bundle file."""
    files = sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.name != "MANIFEST.sha256"
    )
    (directory / "MANIFEST.sha256").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in files)
    )


def completed_restarts(run: Path, coordinate: str = "direct_pt") -> list[Path]:
    """Return converged reaction-path endpoints in recorded order.

    Args:
        run: Coordinate-driving run directory.
        coordinate: Direct proton transfer or asynchronous E--S-to-TI1.

    Returns:
        Ordered Amber restart paths for finalized windows.
    """
    layouts = {
        "direct_pt": ("01_coordinate_drive", "02_unrestrained_endpoint"),
        "asynchronous": (
            "01_coupled_coordinate_drive",
            "02_unrestrained_ti1_endpoint",
        ),
        "attack_only": ("windows", "02_unrestrained_endpoint"),
    }
    if coordinate not in layouts:
        raise ValueError(f"Unsupported reaction-path coordinate: {coordinate}")
    drive_directory, release_directory = layouts[coordinate]
    status = json.loads((run / "status.json").read_text())
    paths = [
        run / drive_directory / label / "optimized.ncrst"
        for label in status["completed_windows"]
    ]
    if coordinate == "attack_only":
        start = Path(status["starting_restart"])
        paths.insert(0, start)
    release = run / release_directory / "optimized.ncrst" if release_directory else None
    extension_status_path = run / "extension_status.json"
    if coordinate == "asynchronous" and extension_status_path.is_file():
        extension_status = json.loads(extension_status_path.read_text())
        paths.extend(
            run / "03_coordinate_drive_extension_to_ti1" / label / "optimized.ncrst"
            for label in extension_status["completed_windows"]
        )
        extended_release = (
            run / "04_unrestrained_extended_ti1_endpoint" / "optimized.ncrst"
        )
        if extension_status["state"] == "completed" and extended_release.is_file():
            paths.append(extended_release)
    elif release is not None and status["state"] == "completed" and release.is_file():
        paths.append(release)
    if not paths or not all(path.is_file() for path in paths):
        raise ValueError("No complete, internally consistent endpoint set")
    return paths


def vmd_script(
    frame_count: int,
    trajectory_name: str = "reaction_path.dcd",
    description: str = "finalized pathway",
) -> str:
    """Return the VMD loader and reactive-site display configuration."""
    if frame_count < 1:
        raise ValueError("frame_count must be positive")
    if not trajectory_name.endswith(".dcd") or any(
        character in trajectory_name for character in "/\\"
    ):
        raise ValueError("trajectory_name must be a local DCD filename")
    if not description or any(character in description for character in '"\n'):
        raise ValueError("description must be nonempty single-line text")
    qm = " ".join(str(index - 1) for index in QM_INDICES)
    return f"""set here [file dirname [file normalize [info script]]]
mol new [file join $here stripped.parm7] type parm7 waitfor all
mol addfile [file join $here {trajectory_name}] type dcd waitfor all
set molid [molinfo top]
mol delrep 0 $molid

# Remove static topology bonds whose connectivity changes along the path.
proc remove_reactive_bond {{molid first second}} {{
    foreach {{a b}} [list $first $second $second $first] {{
        set selected [atomselect $molid "index $a"]
        set neighbors [lindex [$selected getbonds] 0]
        set position [lsearch -exact $neighbors $b]
        if {{$position >= 0}} {{
            $selected setbonds [list [lreplace $neighbors $position $position]]
        }}
        $selected delete
    }}
}}
remove_reactive_bond $molid {THR233_OG1} {THR233_HG1}
remove_reactive_bond $molid {THR233_N} {THR233_HG1}
remove_reactive_bond $molid {PA6_SCISSILE_N} {THR233_HG1}
remove_reactive_bond $molid {THR233_OG1} {PA6_CARBONYL_C}
remove_reactive_bond $molid {PA6_CARBONYL_C} {PA6_SCISSILE_N}

# Full protein cartoon is retained but hidden by default.
mol representation NewCartoon 0.3 10.0 4.1 0
mol color ColorID 7
mol selection "protein"
mol material Opaque
mol addrep $molid
mol showrep $molid 0 0

# QM subsystem.
mol representation Licorice 0.2 12.0 12.0
mol color Name
mol selection "index {qm}"
mol material Opaque
mol addrep $molid

# Non-QM PA6 atoms use element colors.
mol representation Licorice 0.18 12.0 12.0
mol color Element
mol selection "resid 1313 to 1316 and not index {qm}"
mol material Opaque
mol addrep $molid

# Only non-QM substrate carbons receive an orange overlay.
mol representation Licorice 0.22 12.0 12.0
mol color ColorID 3
mol selection "resid 1313 to 1316 and element C and not index {qm}"
mol material Opaque
mol addrep $molid

# Nearby non-QM pocket atoms.
mol representation Lines 1.0
mol color Name
mol selection "protein and within 5.0 of resid 1313 to 1316 and not index {qm}"
mol material Opaque
mol addrep $molid

# Separate dynamic representations prevent one permissive cutoff from drawing
# a broken proton bond. Proton-transfer bonds use 1.40 A, allowing HG1 to
# float between heavy atoms; the forming heavy-atom bond uses 1.90 A.
foreach {{first second cutoff}} {{
    {THR233_OG1} {THR233_HG1} 1.40
    {THR233_N} {THR233_HG1} 1.40
    {PA6_SCISSILE_N} {THR233_HG1} 1.40
    {THR233_OG1} {PA6_CARBONYL_C} 1.90
    {PA6_CARBONYL_C} {PA6_SCISSILE_N} 1.90
}} {{
    mol representation DynamicBonds $cutoff 0.12 12.0
    mol color Name
    mol selection "index $first $second"
    mol material Opaque
    mol addrep $molid
}}

# Transferring HG1 as a small sphere.
mol representation VDW 0.28 18.0
mol color ColorID 4
mol selection "index {THR233_HG1}"
mol material Opaque
mol addrep $molid

label add Bonds $molid/{THR233_OG1} $molid/{THR233_HG1}
label add Bonds $molid/{THR233_N} $molid/{THR233_HG1}
label add Bonds $molid/{THR233_HG1} $molid/{PA6_SCISSILE_N}
label add Bonds $molid/{THR233_OG1} $molid/{PA6_CARBONYL_C}
label add Bonds $molid/{PA6_CARBONYL_O} $molid/{PHE113_N}
label add Bonds $molid/{PA6_CARBONYL_O} $molid/{ASN186_ND2}
label add Angles $molid/{THR233_OG1} $molid/{PA6_CARBONYL_C} $molid/{PA6_CARBONYL_O}
label textsize 1.0
animate goto 0
display resetview
puts "Loaded {frame_count} {description} frames."
puts "Labels: OG1-HG1, Thr233(N)-HG1, HG1-substrate(N), OG1-C, O-Phe113(N), O-Asn186(ND2), and OG1-C-O Burgi-Dunitz angle."
"""


def run_cpptraj(cpptraj: Path, input_path: Path, cwd: Path, log: Path) -> None:
    """Run cpptraj with an explicit Amber environment and retained log."""
    environment = os.environ.copy()
    amberhome = cpptraj.parent.parent
    environment.update(AMBERHOME=str(amberhome), OMP_NUM_THREADS="1")
    library = str(amberhome / "lib")
    environment["LD_LIBRARY_PATH"] = (
        f"{library}:{environment['LD_LIBRARY_PATH']}"
        if environment.get("LD_LIBRARY_PATH")
        else library
    )
    with log.open("w") as handle:
        subprocess.run(
            [str(cpptraj.resolve()), "-i", str(input_path.resolve())],
            cwd=cwd,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def export(
    run: Path,
    topology: Path,
    cpptraj: Path,
    output: Path,
    coordinate: str = "direct_pt",
) -> Path:
    """Publish one immutable progress bundle and update the ``latest`` link."""
    restarts = completed_restarts(run, coordinate)
    output.mkdir(parents=True, exist_ok=True)
    generation = output / f"update_{len(restarts):04d}"
    if generation.exists() and coordinate == "asynchronous":
        generation = output / f"update_{len(restarts):04d}_extended_path"
    if generation.exists():
        raise FileExistsError(generation)
    generation.mkdir()
    trajectory_input = generation / "trajectory.cpptraj"
    trajectory_input.write_text(
        "\n".join(
            [f"parm {topology.resolve()}"]
            + [f"trajin {path.resolve()}" for path in restarts]
            + [
                f"strip {STRIP_MASK}",
                f"trajout {(generation / 'reaction_path.dcd').resolve()} dcd",
                "run",
                "quit",
            ]
        )
        + "\n"
    )
    topology_input = generation / "topology.cpptraj"
    topology_input.write_text(
        f"parm {topology.resolve()}\nparmstrip {STRIP_MASK}\n"
        f"parmwrite out {(generation / 'stripped.parm7').resolve()}\nquit\n"
    )
    run_cpptraj(
        cpptraj, trajectory_input, generation, generation / "trajectory_export.log"
    )
    run_cpptraj(cpptraj, topology_input, generation, generation / "topology_export.log")
    frames, atoms = dcd_counts(generation / "reaction_path.dcd")
    if frames != len(restarts):
        raise ValueError("Exported DCD frame count differs from finalized windows")
    (generation / "reaction_path.vmd").write_text(vmd_script(frames))
    manifest: dict[str, Any] = {
        "frames": frames,
        "stripped_atoms": atoms,
        "source_run": str(run.resolve()),
        "coordinate": coordinate,
        "source_status_sha256": sha256(run / "status.json"),
        "source_restarts": [
            {"path": str(path.resolve()), "sha256": sha256(path)} for path in restarts
        ],
        "strip_mask": STRIP_MASK,
        "files": {},
    }
    for path in sorted(generation.iterdir()):
        if path.is_file():
            manifest["files"][path.name] = sha256(path)
    (generation / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    write_checksum_manifest(generation)
    latest = output / "latest"
    pending = output / ".latest.pending"
    if pending.exists() or pending.is_symlink():
        pending.unlink()
    pending.symlink_to(generation.name)
    pending.replace(latest)
    return generation


def main() -> None:
    """Parse paths and export the finalized progress frames."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--coordinate",
        choices=("direct_pt", "asynchronous", "attack_only"),
        default="direct_pt",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    generation = export(
        args.run, args.topology, args.cpptraj, args.output, args.coordinate
    )
    LOGGER.info("Published %s", generation)


if __name__ == "__main__":
    main()
