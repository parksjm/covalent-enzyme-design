#!/usr/bin/env python3
"""Export a compact, mapped VMD trace for both DFTB3/MM downhill endpoints."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import tempfile
import sys
from pathlib import Path
from typing import Any

import numpy as np
import parmed as pmd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.vmd_style_refresh import append_file_policy  # noqa: E402
from scripts.export_direct_pt_progress import QM_INDICES  # noqa: E402

CAMPAIGN = Path(
    "qmmm/pa6_flipped_model30/05_neutral_reaction_path/"
    "07_dimer_ts_img07_img08_seed_20261113"
)
DISPLAY_STRIP_MASK = "!:1-1316,8602,22254,23286"
DISPLAY_RESIDUES = set(range(1, 1317)) | {8602, 22254, 23286}
SAMPLE_FRAMES_PER_BRANCH = 49
ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl"}
REACTIVE_PAIRS = {
    "Thr233_OG1_substrate_C": (13916, 19916),
    "Thr233_OG1_HG1": (13916, 13917),
    "HG1_substrate_leaving_N": (13917, 19929),
    "substrate_C_leaving_N": (19916, 19929),
}


def digest(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_residue_mask(mask: str) -> set[int]:
    """Convert a one-based comma/range residue mask to zero-based indices."""
    residues: set[int] = set()
    for token in mask.lstrip(":").split(","):
        bounds = [int(value) for value in token.split("-")]
        first = bounds[0]
        last = bounds[-1]
        residues.update(range(first - 1, last))
    return residues


def read_selected_xyz(
    path: Path,
    selected: set[int],
    expected_symbols: list[str],
) -> tuple[dict[int, np.ndarray], int]:
    """Read selected frames from an XYZ trajectory and validate atom counts."""
    atom_count = len(expected_symbols)
    frames: dict[int, np.ndarray] = {}
    frame_index = 0
    with path.open() as handle:
        while line := handle.readline():
            try:
                frame_atoms = int(line.strip())
            except ValueError as exc:
                raise ValueError(f"Invalid XYZ atom count in {path}") from exc
            if frame_atoms != atom_count:
                raise ValueError(f"Unexpected atom count {frame_atoms} in {path}")
            if not handle.readline():
                raise ValueError(f"Missing XYZ comment line in {path}")
            coordinates = np.empty((atom_count, 3), dtype=np.float64)
            for atom_index in range(atom_count):
                fields = handle.readline().split()
                if len(fields) != 4:
                    raise ValueError(f"Incomplete XYZ frame in {path}")
                if fields[0] != expected_symbols[atom_index]:
                    raise ValueError(
                        f"Active XYZ atom order mismatch at {atom_index} in {path}"
                    )
                if frame_index in selected:
                    coordinates[atom_index] = [float(value) for value in fields[1:]]
            if frame_index in selected:
                frames[frame_index] = coordinates
            frame_index += 1
    if set(frames) != selected:
        raise ValueError(f"Selected XYZ frames missing from {path}")
    return frames, frame_index


def dcd_counts(path: Path) -> tuple[int, int]:
    """Read frame and atom counts from a little-endian coordinate DCD."""
    with path.open("rb") as handle:
        if handle.read(4) != struct.pack("<i", 84):
            raise ValueError("Unsupported DCD record format")
        header = handle.read(84)
        if len(header) != 84 or header[:4] != b"CORD":
            raise ValueError("Invalid DCD header")
        frames = struct.unpack("<i", header[4:8])[0]
        if handle.read(4) != struct.pack("<i", 84):
            raise ValueError("Invalid DCD header trailer")
        title_size = struct.unpack("<i", handle.read(4))[0]
        handle.seek(title_size + 4, os.SEEK_CUR)
        if handle.read(4) != struct.pack("<i", 4):
            raise ValueError("Invalid DCD atom count record")
        atoms = struct.unpack("<i", handle.read(4))[0]
    if frames < 1 or atoms < 1:
        raise ValueError("Empty DCD")
    return frames, atoms


def write_mdcrd(path: Path, frames: list[np.ndarray]) -> None:
    """Write full-system coordinates in Amber ASCII trajectory format."""
    with path.open("w") as handle:
        handle.write(
            "Endpoint minimization display coordinates; XYZ frames in Angstrom\n"
        )
        for coordinates in frames:
            values = coordinates.reshape(-1)
            for offset in range(0, len(values), 10):
                handle.write(
                    "".join(f"{value:8.3f}" for value in values[offset : offset + 10])
                )
                handle.write("\n")


def vmd_script() -> str:
    """Return the viewer script with active-site representations and labels."""
    qm_atoms = [
        1621,
        1655,
        2727,
        2734,
        13905,
        13931,
        14555,
        14560,
        14574,
        14579,
        19914,
        19917,
        19925,
        19931,
        19937,
        19941,
    ]
    qm_indices: list[int] = []
    for first, last in zip(qm_atoms[::2], qm_atoms[1::2]):
        qm_indices.extend(range(first - 1, last))
    qm_selection = " ".join(str(index) for index in qm_indices)
    return f"""set here [file dirname [file normalize [info script]]]
mol new [file join $here stripped.parm7] type parm7 waitfor all
mol addfile [file join $here trajectory.dcd] type dcd waitfor all
set molid [molinfo top]
mol delrep 0 $molid

proc remove_reactive_bond {{molid first second}} {{
    foreach {{a b}} [list $first $second $second $first] {{
        set selected [atomselect $molid "index $a"]
        set neighbors [lindex [$selected getbonds] 0]
        set position [lsearch -exact $neighbors $b]
        if {{$position >= 0}} {{$selected setbonds [list [lreplace $neighbors $position $position]]}}
        $selected delete
    }}
}}
foreach {{first second}} {{13915 13916 13904 13916 13916 19928 13915 19915 19915 19928}} {{
    remove_reactive_bond $molid $first $second
}}

mol representation NewCartoon 0.3 10.0 4.1 0
mol color ColorID 7
mol selection "protein"
mol material Opaque
mol addrep $molid
mol showrep $molid 0 0

mol representation Licorice 0.2 12.0 12.0
mol color Name
mol selection "index {qm_selection}"
mol material Opaque
mol addrep $molid

mol representation Licorice 0.18 12.0 12.0
mol color Element
mol selection "resid 1313 to 1316 and not index {qm_selection}"
mol material Opaque
mol addrep $molid

mol representation Licorice 0.20 12.0 12.0
mol color Name
mol selection "resname WAT"
mol material Opaque
mol addrep $molid

mol representation Lines 1.0
mol color Name
mol selection "protein and within 5.0 of resid 1313 to 1316 and not index {qm_selection}"
mol material Opaque
mol addrep $molid

foreach {{first second cutoff}} {{
    13915 13916 1.40
    13904 13916 1.40
    13916 19928 1.40
    13915 19915 1.90
    19915 19928 1.90
}} {{
    mol representation DynamicBonds $cutoff 0.12 12.0
    mol color Name
    mol selection "index $first $second"
    mol material Opaque
    mol addrep $molid
}}

label add Bonds $molid/13915 $molid/19915
label add Bonds $molid/13915 $molid/13916
label add Bonds $molid/13916 $molid/19928
label add Bonds $molid/19915 $molid/19928
label textsize 1.0
animate goto 0
display resetview
puts "Endpoint minimization trace: E-S frames 0-48; explicit TS frame 49; acylenzyme frames 50-98."
puts "This is not an IRC or MEP. Both branches are unrestrained 0 K DFTB3/MM endpoint minimizations."
"""


def main() -> None:
    """Export the immutable, downsampled DFTB3/MM minimization movie."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument(
        "--cpptraj",
        type=Path,
        default=Path("/data/jpc/code/ambertools26/venv/bin/cpptraj"),
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    campaign = repo / CAMPAIGN
    irc = campaign / "irc_0K"
    endpoint = irc / "endpoint_minimizations"
    ts_restart = campaign / "continuation_cycle200_20261113/calculation/optimized.ncrst"
    full_topology = campaign / "inputs/system.prmtop"
    source_topology = campaign / "vmd_progress/latest/stripped.parm7"
    analysis_path = irc / "plots/endpoint_minimization_analysis.json"
    distance_path = irc / "inputs/distance_pairs.json"
    preparation_path = campaign / "preparation.json"
    output = (
        args.output_dir.resolve()
        if args.output_dir
        else irc / "vmd_progress/update_20261113_endpoint_minimizations"
    )
    latest = output.parent / "latest"
    if output.exists() or latest.exists():
        raise FileExistsError(f"Refusing to overwrite {output} or {latest}")
    required = [
        ts_restart,
        full_topology,
        source_topology,
        analysis_path,
        distance_path,
        preparation_path,
    ]
    required.extend(
        endpoint / mode / "calculation/path.xyz" for mode in ("mode_minus", "mode_plus")
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    analysis = json.loads(analysis_path.read_text())
    preparation = json.loads(preparation_path.read_text())
    active_mask = preparation["active_mask"]
    for mode in ("mode_minus", "mode_plus"):
        mdin = (endpoint / mode / "calculation/dlfind.mdin").read_text()
        match = re.search(r"active='([^']+)'", mdin)
        if not match or match.group(1) != active_mask:
            raise ValueError(f"Active mask mismatch with preparation.json: {mode}")
    if (
        not analysis["all_endpoints_complete"]
        or not analysis["final_frames_match_minimized_restarts"]
    ):
        raise ValueError(
            "Endpoint analysis does not certify both complete minimizations"
        )
    for mode in ("mode_minus", "mode_plus"):
        state = json.loads((endpoint / mode / "job_run/status.json").read_text())
        if state["exit_status"] != 0:
            raise ValueError(f"Endpoint minimization failed: {mode}")

    full_parm = pmd.load_file(str(full_topology))
    display_parm = pmd.load_file(str(source_topology))
    ts = pmd.load_file(str(ts_restart))
    ts_coordinates = np.asarray(ts.coordinates, dtype=np.float64).reshape((-1, 3))
    if len(full_parm.atoms) != len(ts_coordinates):
        raise ValueError("TS restart atom count differs from full topology")

    retained_original = [
        atom.idx for atom in full_parm.atoms if atom.residue.idx + 1 in DISPLAY_RESIDUES
    ]
    if len(retained_original) != len(display_parm.atoms):
        raise ValueError("Retained atom map differs from the stripped VMD topology")
    for original_idx, stripped_atom in zip(retained_original, display_parm.atoms):
        original_atom = full_parm.atoms[original_idx]
        if (
            original_atom.name != stripped_atom.name
            or original_atom.residue.name != stripped_atom.residue.name
        ):
            raise ValueError(
                f"Stripped topology atom identity mismatch at {stripped_atom.idx}"
            )
    active_residue_indices = parse_residue_mask(active_mask)
    active_original = [
        atom.idx
        for atom in full_parm.atoms
        if atom.residue.idx in active_residue_indices
    ]
    if len(active_original) != 4165:
        raise ValueError(f"Expected 4165 active atoms, found {len(active_original)}")
    active_symbols = [
        ELEMENTS.get(full_parm.atoms[index].atomic_number) for index in active_original
    ]
    if any(symbol is None for symbol in active_symbols):
        raise ValueError("Unsupported active-atom element in topology")

    branch_frames: dict[str, dict[int, np.ndarray]] = {}
    source_counts: dict[str, int] = {}
    selected_by_mode: dict[str, list[int]] = {}
    for mode in ("mode_minus", "mode_plus"):
        path = endpoint / mode / "calculation/path.xyz"
        frame_count = int(analysis["branches"][mode]["frame_count"])
        if frame_count < SAMPLE_FRAMES_PER_BRANCH:
            raise ValueError(f"Not enough frames in {path}")
        if mode == "mode_minus":
            selected = np.rint(
                np.linspace(frame_count - 1, 0, SAMPLE_FRAMES_PER_BRANCH)
            ).astype(int)
        else:
            selected = np.rint(
                np.linspace(0, frame_count - 1, SAMPLE_FRAMES_PER_BRANCH)
            ).astype(int)
        frames, parsed_count = read_selected_xyz(
            path, set(map(int, selected)), active_symbols
        )
        if parsed_count != frame_count:
            raise ValueError(f"XYZ frame count changed while reading {path}")
        branch_frames[mode] = frames
        source_counts[mode] = frame_count
        selected_by_mode[mode] = list(map(int, selected))

    records: list[tuple[str, str, int | None, np.ndarray]] = []
    for source_frame in selected_by_mode["mode_minus"]:
        records.append(
            (
                "E-S",
                "mode_minus",
                source_frame,
                branch_frames["mode_minus"][source_frame],
            )
        )
    records.append(
        ("TS", "optimized_ts_restart", None, ts_coordinates[active_original])
    )
    for source_frame in selected_by_mode["mode_plus"]:
        records.append(
            (
                "acylenzyme",
                "mode_plus",
                source_frame,
                branch_frames["mode_plus"][source_frame],
            )
        )
    if len(records) != 99:
        raise ValueError(f"Expected 99 display records, found {len(records)}")

    dcd_frames: list[np.ndarray] = []
    distances: list[dict[str, float]] = []
    topology_pairs = {
        label: tuple(index - 1 for index in pair)
        for label, pair in REACTIVE_PAIRS.items()
    }
    for _, _, _, active_coordinates in records:
        coordinates = ts_coordinates.copy()
        coordinates[active_original] = active_coordinates
        dcd_frames.append(coordinates)
        distances.append(
            {
                label: float(np.linalg.norm(coordinates[first] - coordinates[second]))
                for label, (first, second) in topology_pairs.items()
            }
        )

    es_distances = distances[0]
    acyl_distances = distances[-1]
    if not (
        es_distances["Thr233_OG1_substrate_C"] > 2.5
        and es_distances["Thr233_OG1_HG1"] < 1.2
        and es_distances["substrate_C_leaving_N"] < 1.7
    ):
        raise ValueError(f"Mode-minus endpoint fails E-S geometry gate: {es_distances}")
    if not (
        acyl_distances["Thr233_OG1_substrate_C"] < 1.6
        and acyl_distances["HG1_substrate_leaving_N"] < 1.3
        and acyl_distances["substrate_C_leaving_N"] > 2.0
    ):
        raise ValueError(
            f"Mode-plus endpoint fails acylenzyme geometry gate: {acyl_distances}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    shutil.copy2(source_topology, output / "stripped.parm7")
    (output / "view_downhill_minimizations.vmd").write_text(vmd_script())
    append_file_policy(
        output / "view_downhill_minimizations.vmd",
        output / "stripped.parm7",
        [i - 1 for i in QM_INDICES],
    )
    (output / "README.md").write_text(
        "# PA6 DFTB3/MM endpoint minimization trace\n\n"
        "Load with `source /absolute/path/view_downhill_minimizations.vmd`. "
        "The 99-frame movie shows the stripped 19,976-atom display system, including "
        "the complete enzyme, substrate, and three local waters. Each frame is first "
        "reconstructed in the 104,922-atom system from all 4,165 active coordinates; "
        "the existing display strip is then applied.\n\n"
        "Frames 0–48 traverse the reversed mode-minus endpoint minimization from the "
        "E-S minimum toward the TS. Frame 49 is the optimized TS restart. Frames "
        "50–98 traverse the forward mode-plus minimization to the acylenzyme minimum. "
        "Reactive bonds are dynamic and key distances are labeled.\n\n"
        "This is an endpoint minimization trace, not an IRC or MEP and not a dynamical trajectory. "
        "The endpoint minimizations are unrestrained 0 K DFTB3/MM DL-FIND L-BFGS runs. "
        "See `frame_manifest.json` and `../../plots/endpoint_minimization_analysis.json` "
        "for source mapping and acceptance evidence.\n"
    )

    with tempfile.TemporaryDirectory(prefix="downhill-vmd-") as temp_dir:
        temp = Path(temp_dir)
        mdcrd = temp / "combined.mdcrd"
        control = temp / "export.cpptraj"
        log_path = output / "trajectory_export.log"
        write_mdcrd(mdcrd, dcd_frames)
        control.write_text(
            f"parm {full_topology}\n"
            f"trajin {mdcrd}\n"
            f"strip {DISPLAY_STRIP_MASK}\n"
            f"trajout {output / 'trajectory.dcd'} dcd\n"
            "run\nquit\n"
        )
        result = subprocess.run(
            [str(args.cpptraj), "-i", str(control)],
            env={
                **os.environ,
                "AMBERHOME": str(args.cpptraj.resolve().parent.parent),
                "LD_LIBRARY_PATH": str(args.cpptraj.resolve().parent.parent / "lib"),
                "OMP_NUM_THREADS": "1",
            },
            capture_output=True,
            text=True,
            check=False,
        )
        log_path.write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"cpptraj export failed; see {log_path}")
        frame_count, atom_count = dcd_counts(output / "trajectory.dcd")
        if (frame_count, atom_count) != (99, len(display_parm.atoms)):
            raise ValueError(f"DCD count mismatch: {(frame_count, atom_count)}")

    with (output / "reactive_distances.csv").open("w", newline="") as handle:
        columns = [
            "frame",
            "segment",
            "source_mode",
            "source_frame_1based",
            *REACTIVE_PAIRS,
        ]
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for frame_index, ((segment, mode, source_frame, _), row) in enumerate(
            zip(records, distances)
        ):
            writer.writerow(
                {
                    "frame": frame_index,
                    "segment": segment,
                    "source_mode": mode,
                    "source_frame_1based": ""
                    if source_frame is None
                    else source_frame + 1,
                    **row,
                }
            )

    source_paths = [
        ts_restart,
        full_topology,
        source_topology,
        analysis_path,
        distance_path,
        preparation_path,
    ]
    source_paths.extend(
        endpoint / mode / "calculation/path.xyz" for mode in ("mode_minus", "mode_plus")
    )
    manifest: dict[str, Any] = {
        "description": "Downsampled endpoint minimization trace; not an IRC, MEP, or dynamics trajectory.",
        "frame_count": 99,
        "sequence": {"E-S": [0, 48], "optimized_TS": 49, "acylenzyme": [50, 98]},
        "display_atom_count": len(display_parm.atoms),
        "original_system_atom_count": len(full_parm.atoms),
        "active_atom_count": len(active_original),
        "active_mask_residues_1based": active_mask,
        "active_atom_order": "Topology order within residues selected by the recorded DL-FIND active mask.",
        "display_strip_mask": DISPLAY_STRIP_MASK,
        "retained_local_water_residues_1based": [8602, 22254, 23286],
        "full_system_reconstruction": "Insert all 4,165 active XYZ coordinates into the 104,922-atom TS restart before applying the existing display strip.",
        "source_xyz_frame_counts": source_counts,
        "selected_source_frames_1based": {
            mode: [frame + 1 for frame in indices]
            for mode, indices in selected_by_mode.items()
        },
        "source_sha256": {str(path): digest(path) for path in source_paths},
        "endpoint_geometry_angstrom": {
            "E-S": es_distances,
            "acylenzyme": acyl_distances,
        },
        "export": {
            "cpptraj": str(args.cpptraj),
            "topology": "stripped.parm7",
            "trajectory": "trajectory.dcd",
        },
    }
    (output / "frame_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.copy2(Path(__file__).resolve(), output / Path(__file__).name)
    checksum_files = sorted(
        path
        for path in output.iterdir()
        if path.is_file() and path.name != "MANIFEST.sha256"
    )
    (output / "MANIFEST.sha256").write_text(
        "".join(f"{digest(path)}  {path.name}\n" for path in checksum_files)
    )
    os.symlink(output.name, latest)
    print(f"Wrote {output}: {frame_count} frames x {atom_count} atoms")
    print(f"E-S endpoint distances: {es_distances}")
    print(f"Acylenzyme endpoint distances: {acyl_distances}")


if __name__ == "__main__":
    main()
