"""Build a stripped VMD package for one Amber optimizer trajectory."""

from __future__ import annotations

import gzip
import json
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm, NetCDFTraj

from .vmd_export_io import select_display_atoms, sha, validate_dcd

LOGGER = logging.getLogger(__name__)
REACTIVE_LABELS = (
    ("OG1_HG1", 13915, "OG1", 13916, "HG1"),
    ("OG1_C54", 13915, "OG1", 19902, "C54"),
    ("C54_O49", 19902, "C54", 19903, "O49"),
    ("C54_N59", 19902, "C54", 19904, "N59"),
    ("Asp272_OD1_HG1", 14558, "OD1", 13916, "HG1"),
    ("Thr_N_HG1", 13904, "N", 13916, "HG1"),
    ("O49_Phe113_N", 19903, "O49", 1637, "N"),
    ("O49_Asn186_ND2", 19903, "O49", 2731, "ND2"),
)
DYNAMIC_BONDS = (
    (13915, 13916, 1.25),
    (14558, 13916, 1.25),
    (13904, 13916, 1.25),
    (13915, 19902, 1.9),
)
BURGI_DUNITZ_ANGLE = (13915, "OG1", 19902, "C54", 19903, "O49")


def parse_atom_mask(mask: str, atom_count: int) -> list[int]:
    """Parse an explicit Amber numeric atom mask into zero-based indices."""
    if not mask.startswith("@"):
        raise ValueError("QM mask must use explicit numeric atom indices")
    indices: list[int] = []
    for token in mask[1:].split(","):
        bounds = token.split("-")
        if len(bounds) not in {1, 2} or not all(part.isdigit() for part in bounds):
            raise ValueError(f"Invalid numeric atom-mask token: {token}")
        start = int(bounds[0])
        stop = int(bounds[-1])
        if start < 1 or stop < start or stop > atom_count:
            raise ValueError(f"Atom-mask token is outside topology: {token}")
        indices.extend(range(start - 1, stop))
    if not indices or len(indices) != len(set(indices)):
        raise ValueError("QM atom mask must select unique atoms")
    return indices


def residue_mask(topology: AmberParm, atom_indices: np.ndarray) -> str:
    """Convert selected atom indices to a compact 1-based residue mask."""
    residues = sorted(
        {topology.atoms[int(index)].residue.idx + 1 for index in atom_indices}
    )
    if not residues:
        raise ValueError("Display atom selection is empty")
    ranges: list[str] = []
    start = previous = residues[0]
    for value in residues[1:] + [None]:
        if value is not None and value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        if value is not None:
            start = previous = value
    return ":" + ",".join(ranges)


def read_frames(
    trajectory_path: Path,
    atom_count: int,
    fixed_box: np.ndarray,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Load finite optimizer frames and require the declared fixed box."""
    trajectory = NetCDFTraj.open_old(str(trajectory_path))
    try:
        coordinates = np.asarray(trajectory.coordinates)
        boxes = np.asarray(trajectory.box)
    finally:
        trajectory.close()
    if coordinates.ndim != 3 or coordinates.shape[1:] != (atom_count, 3):
        raise ValueError(f"Unexpected trajectory coordinate shape: {coordinates.shape}")
    if boxes.shape != (len(coordinates), 6) or len(coordinates) == 0:
        raise ValueError("Trajectory must contain frames with six box values")
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for xyz, raw_box in zip(coordinates, boxes, strict=True):
        if not np.isfinite(xyz).all() or not np.isfinite(raw_box).all():
            raise ValueError("Trajectory contains nonfinite coordinates or box")
        if np.all(raw_box == 0):
            box = fixed_box.copy()
        elif np.allclose(raw_box, fixed_box, atol=1e-8, rtol=0):
            box = raw_box.copy()
        else:
            raise ValueError("Trajectory box conflicts with the declared fixed box")
        states.append((xyz, box))
    return states


def _viewer_text(
    atom_count: int,
    qm_indices: list[int],
    remap: dict[int, int],
    topology: AmberParm,
) -> str:
    """Create VMD loader with explicit QM licorice and current reactive labels."""
    labels = REACTIVE_LABELS
    for _, a, name_a, b, name_b in labels:
        for original, expected in ((a, name_a), (b, name_b)):
            if original not in remap or topology.atoms[original].name != expected:
                raise ValueError(
                    f"Reactive label identity mismatch at {original}/{expected}"
                )
    mapped = [remap[index] for index in qm_indices]
    lines = [
        "# Stage2 PA6 coordinate-drive optimizer; optimizer states are not MD time.",
        "set ::stage2_dir [file dirname [file normalize [info script]]];",
        "set _topology [file join $::stage2_dir stripped.parm7]",
        'if {![file isfile $_topology]} {exec gzip -dc "$_topology.gz" > $_topology}',
        "set ::stage2_mol [mol new $_topology type parm7 waitfor all]",
        "mol addfile [file join $::stage2_dir optimization.dcd] type dcd waitfor all molid $::stage2_mol",
        f'if {{[molinfo $::stage2_mol get numatoms] != {atom_count}}} {{error "Unexpected stripped atom count"}}',
        f"set ::stage2_qm {{{' '.join(map(str, mapped))}}}",
        "mol delrep 0 $::stage2_mol",
        'set _qm_selection "index $::stage2_qm"',
        'set _water_selection "same residue as (water and within 5.0 of ($_qm_selection))"',
        "foreach {selection style} [list protein NewCartoon {resname LIG} Licorice $_qm_selection Licorice $_water_selection Licorice] {",
        "    mol representation $style; mol color Element; mol selection $selection; mol addrep $::stage2_mol",
        "}",
        'mol representation Licorice 0.31 12.0 12.0; mol color ColorID 3; mol selection "resname LIG and element C and not ($_qm_selection)"; mol addrep $::stage2_mol',
        "mol showrep $::stage2_mol 0 off",
        "proc ::stage2_remove_bond {molid a b} {",
        "    foreach {first second} [list $a $b $b $a] {",
        '        set selected [atomselect $molid "index $first"]',
        "        set neighbors [lindex [$selected getbonds] 0]",
        "        set position [lsearch -exact $neighbors $second]",
        "        if {$position >= 0} {$selected setbonds [list [lreplace $neighbors $position $position]]}",
        "        $selected delete",
        "    }",
        "}",
    ]
    # The transferred proton is always included in QM and has an explicit small-sphere highlight.
    proton_index = remap[13916]
    lines[-1] = (
        "mol representation VDW 0.32 16; mol color ColorID 7; "
        f'mol selection "index {proton_index}"; mol addrep $::stage2_mol'
    )
    for a, b, cutoff in DYNAMIC_BONDS:
        lines.extend(
            (
                f"::stage2_remove_bond $::stage2_mol {remap[a]} {remap[b]}",
                f"mol representation DynamicBonds {cutoff} 0.12 12",
                "mol color ColorID 1",
                f'mol selection "index {remap[a]} {remap[b]}"',
                "mol addrep $::stage2_mol",
            )
        )
    for title, a, name_a, b, name_b in labels:
        mapped_a, mapped_b = remap[a], remap[b]
        lines.extend(
            [
                f'set _a [atomselect $::stage2_mol "index {mapped_a} and name {name_a}"]',
                f'set _b [atomselect $::stage2_mol "index {mapped_b} and name {name_b}"]',
                f'if {{[$_a num] != 1 || [$_b num] != 1}} {{error "Reactive label identity mismatch: {title}"}}',
                "$_a delete; $_b delete",
                "set _label [llength [label list Bonds]]",
                f"label add Bonds $::stage2_mol/{mapped_a} $::stage2_mol/{mapped_b}",
                "label show Bonds $_label",
            ]
        )
    a, name_a, b, name_b, c, name_c = BURGI_DUNITZ_ANGLE
    for original, expected in ((a, name_a), (b, name_b), (c, name_c)):
        if original not in remap or topology.atoms[original].name != expected:
            raise ValueError(
                f"Burgi-Dunitz atom identity mismatch at {original}/{expected}"
            )
    angle_id = "[llength [label list Angles]]"
    lines.extend(
        (
            f"set _angle {angle_id}",
            f"label add Angles $::stage2_mol/{remap[a]} $::stage2_mol/{remap[b]} $::stage2_mol/{remap[c]}",
            "label show Angles $_angle",
        )
    )
    lines.append(
        'puts "Loaded [molinfo $::stage2_mol get numframes] stage2 optimizer frames; reactive bond labels update with frame."'
    )
    return "\n".join(lines) + "\n"


def export_stage2(
    run_dir: Path,
    system_dir: Path,
    inputs_path: Path,
    output: Path,
    cpptraj: Path,
) -> dict[str, Any]:
    """Export and verify one stage2 trajectory without touching source files.

    Args:
        run_dir: Stage directory containing ``minimize.mdin`` and ``optimization.nc``.
        system_dir: Directory containing the source ``md_solvated.prmtop``.
        inputs_path: Original run inputs JSON containing the fixed periodic box.
        output: New destination directory for VMD files.
        cpptraj: AmberTools cpptraj executable.

    Returns:
        Export manifest data with frame count, selected atoms, waters and hashes.
    """
    if output.exists():
        raise FileExistsError(output)
    topology_path = system_dir / "md_solvated.prmtop"
    trajectory_path = run_dir / "optimization.nc"
    mdin_path = run_dir / "minimize.mdin"
    for path in (topology_path, trajectory_path, mdin_path, inputs_path, cpptraj):
        if not path.is_file():
            raise FileNotFoundError(path)
    topology = AmberParm(str(topology_path))
    mdin = mdin_path.read_text()
    match = re.search(r"qmmask\s*=\s*'([^']+)'", mdin)
    if match is None:
        raise ValueError("Could not find explicit qmmask in stage input")
    qm_indices = parse_atom_mask(match.group(1), len(topology.atoms))
    fixed_box = np.asarray(json.loads(inputs_path.read_text())["box"], dtype=float)
    if fixed_box.shape != (6,) or not np.isfinite(fixed_box).all():
        raise ValueError("inputs.json must contain six finite fixed-box values")
    states = read_frames(trajectory_path, len(topology.atoms), fixed_box)
    kept, water_ids = select_display_atoms(topology, states, qm_indices)
    keep_mask = residue_mask(topology, kept)
    remap = {int(old): new for new, old in enumerate(kept)}
    mapped_qm = [remap[index] for index in qm_indices]
    sources = (topology_path, trajectory_path, mdin_path, inputs_path)
    source_hashes = {str(path): sha(path) for path in sources}

    output.mkdir(parents=True)
    parm_script = output / "stripped_topology.cpptraj"
    parm_script.write_text(
        f"parm {topology_path}\nparmstrip !{keep_mask}\nparmwrite out {output / 'stripped.parm7'}\n"
    )
    trajectory_script = output / "export.cpptraj"
    box_names = ("x", "y", "z", "alpha", "beta", "gamma")
    box_line = "box " + " ".join(
        f"{key} {value:.12g}" for key, value in zip(box_names, fixed_box, strict=True)
    )
    trajectory_script.write_text(
        "\n".join(
            (
                f"parm {topology_path}",
                f"trajin {trajectory_path}",
                f"strip !{keep_mask}",
                box_line,
                f"trajout {output / 'optimization.dcd'} dcd namdcell",
                "run",
            )
        )
        + "\n"
    )
    environment = os.environ.copy()
    environment["LD_LIBRARY_PATH"] = str(cpptraj.parent.parent / "lib")
    with (output / "cpptraj.log").open("w") as log:
        subprocess.run(
            [str(cpptraj), "-i", str(parm_script)],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
        )
        subprocess.run(
            [str(cpptraj), "-i", str(trajectory_script)],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
        )
    validation = validate_dcd(
        output / "optimization.dcd", [(xyz[kept], box) for xyz, box in states]
    )
    with (output / "stripped.parm7").open("rb") as incoming:
        with (output / "stripped.parm7.gz").open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                shutil.copyfileobj(incoming, zipped)
    (output / "qm_indices.tcl").write_text(
        "set ::stage2_qm_indices {" + " ".join(map(str, mapped_qm)) + "}\n"
    )
    (output / "view_stage2.vmd").write_text(
        _viewer_text(len(kept), qm_indices, remap, topology)
    )
    manifest: dict[str, Any] = {
        "frame_count": len(states),
        "source_topology_atoms": len(topology.atoms),
        "stripped_atom_count": len(kept),
        "qm_atom_count": len(qm_indices),
        "qm_indices_source_0based": qm_indices,
        "qm_indices_stripped_0based": mapped_qm,
        "retained_water_residues_1based": water_ids,
        "water_selection": "Complete WAT residues with oxygen within 4 A of any QM heavy atom in any source frame; Euclidean coordinates, no reimaging",
        "optimizer_states_are_not_md_times": True,
        "source_sha256": source_hashes,
        "validation": validation,
    }
    (output / "frame_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if any(sha(path) != source_hashes[str(path)] for path in sources):
        raise ValueError("Source input changed during export")
    (output / "README.md").write_text(
        "# Failed stage-II optimizer trajectory\n\n"
        "In VMD, run `source /absolute/path/view_stage2.vmd`. Frames are saved "
        "optimizer states, not equally spaced MD time. The protein cartoon is loaded "
        "but hidden; QM atoms use licorice, non-QM ligand carbons are orange, and "
        "reactive distances and the Burgi-Dunitz angle are shown.\n"
    )
    payload = (
        "README.md",
        "frame_manifest.json",
        "optimization.dcd",
        "qm_indices.tcl",
        "stripped.parm7.gz",
        "view_stage2.vmd",
    )
    (output / "MANIFEST.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in payload)
    )
    LOGGER.info("Exported %d stage2 optimizer frames to %s", len(states), output)
    return manifest
