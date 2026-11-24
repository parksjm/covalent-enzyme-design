"""Export mapped, stripped native relay optimizer snapshots for VMD."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from nyl12_covalent.relay_release import ROLES  # noqa: E402
from nyl12_covalent.stationary_native import digest, write_json  # noqa: E402

STRIP_MASK = ":WAT,Na+,Cl-,NA,CL"


def select_optimizer_frames(count: int, limit: int = 47) -> list[int]:
    """Select unique saved optimizer indices including both endpoints.

    Args:
        count: Positive number of available saved frames.
        limit: Positive maximum selection size, at least two for multiple frames.

    Returns:
        Ordered original frame indices including zero and count minus one.

    Raises:
        ValueError: Counts are invalid or the limit cannot retain both endpoints.
    """
    if type(count) is not int or type(limit) is not int or count < 1 or limit < 1:
        raise ValueError("Frame count and limit must be positive integers")
    if count > 1 and limit < 2:
        raise ValueError("At least two slots are required to retain both endpoints")
    return np.linspace(0, count - 1, num=min(count, limit), dtype=int).tolist()


def export(seed: Path, output: Path) -> None:
    """Snapshot available saved iterations and verify their mapped DCD export.

    Args:
        seed: Prepared/running release seed directory.
        output: Fresh immutable destination for this visualization checkpoint.
    """
    from netCDF4 import Dataset
    from parmed.amber import AmberParm, Rst7
    from parmed.amber.mask import AmberMask
    from nyl12_covalent.vmd_export_io import validate_dcd

    output.mkdir(parents=True, exist_ok=False)
    source = seed / "optimization.nc"
    snapshot = output / "source_snapshot"
    snapshot.mkdir()
    original_source = source
    with (
        Dataset(original_source) as incoming,
        Dataset(snapshot / source.name, "w", format="NETCDF3_64BIT_OFFSET") as outgoing,
    ):
        original_count = len(incoming.dimensions["frame"])
        if original_count < 1:
            raise ValueError("Native trajectory has no saved optimizer frame")
        original_selected = select_optimizer_frames(original_count)
        outgoing.setncatts(incoming.__dict__)
        for key, dimension in incoming.dimensions.items():
            outgoing.createDimension(
                key, None if dimension.isunlimited() else len(dimension)
            )
        for key, variable in incoming.variables.items():
            target = outgoing.createVariable(key, variable.dtype, variable.dimensions)
            target.setncatts(variable.__dict__)
            if "frame" in variable.dimensions:
                if variable.dimensions[0] != "frame":
                    raise ValueError("Unsupported frame dimension order")
                for new, old in enumerate(original_selected):
                    target[new] = variable[old]
            else:
                target[:] = variable[:]
    shutil.copy2(Path(__file__), snapshot / Path(__file__).name)
    source = snapshot / source.name
    topology = seed / "inputs/system.prmtop"
    initial = seed / "inputs/optimized.ncrst"
    parm = AmberParm(str(topology))
    kept = np.flatnonzero(
        ~np.asarray(AmberMask(parm, STRIP_MASK).Selection(), dtype=bool)
    )
    remap = {int(original): new for new, original in enumerate(kept)}
    restart = Rst7.open(str(initial))
    xyz = np.asarray(restart.coordinates).reshape((-1, 3))
    states = [(xyz[kept], np.asarray(restart.box))]
    final_restart = seed / "optimized.ncrst"
    final_state = None
    restart_comparison = None
    if final_restart.exists():
        final = Rst7.open(str(final_restart))
        final_xyz = np.asarray(final.coordinates).reshape((-1, 3))
        with Dataset(original_source) as native:
            last_xyz = np.asarray(native["coordinates"][-1])
        restart_comparison = {
            "maximum_coordinate_difference_A": float(
                np.max(np.abs(last_xyz - final_xyz))
            ),
            "all_atom_float32_equal": bool(
                np.array_equal(
                    last_xyz.astype(np.float32), final_xyz.astype(np.float32)
                )
            ),
            "restart_sha256": digest(final_restart),
        }
        if not restart_comparison["all_atom_float32_equal"]:
            final_state = (final_xyz[kept], np.asarray(final.box))
            shutil.copy2(final_restart, snapshot / "final_optimized.ncrst")
    with Dataset(source) as nc:
        count = len(nc.dimensions["frame"])
        if count < 1:
            raise ValueError("Native trajectory has no saved optimizer frame")
        selected = select_optimizer_frames(count, 46 if final_state is not None else 47)
        for i in selected:
            states.append(
                (
                    np.asarray(nc["coordinates"][i, kept, :]),
                    np.r_[nc["cell_lengths"][i], nc["cell_angles"][i]],
                )
            )
    optimizer_indices = [original_selected[index] for index in selected]
    frame_sources = [{"kind": "initial_seed", "index_0based": None}]
    frame_sources += [
        {"kind": "optimizer_nc", "index_0based": index} for index in optimizer_indices
    ]
    trajectories = "".join(
        f"trajin {source} {index + 1} {index + 1}\n" for index in selected
    )
    if final_state is not None:
        states.append(final_state)
        frame_sources.append({"kind": "optimized_restart", "index_0based": None})
        trajectories += f"trajin {snapshot / 'final_optimized.ncrst'}\n"
    script = output / "export.cpptraj"
    script.write_text(
        f"parm {topology}\ntrajin {initial}\n"
        f"{trajectories}strip {STRIP_MASK}\n"
        f"trajout {output / 'optimization.dcd'} dcd namdcell\nrun\n"
        "quit\n"
    )
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = "/data/jpc/code/ambertools26/venv/lib"
    with (output / "cpptraj.log").open("w") as log:
        subprocess.run(
            ["/data/jpc/code/ambertools26/venv/bin/cpptraj", "-i", str(script)],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=True,
        )
    topology_script = output / "topology.cpptraj"
    topology_script.write_text(
        f"parm {topology}\nparmstrip {STRIP_MASK}\n"
        f"parmwrite out {output / 'stripped.parm7'}\nquit\n"
    )
    with (output / "topology_cpptraj.log").open("w") as log:
        subprocess.run(
            [
                "/data/jpc/code/ambertools26/venv/bin/cpptraj",
                "-i",
                str(topology_script),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=True,
        )
    validation = validate_dcd(output / "optimization.dcd", states)
    mapped = {key: remap[value - 1] for key, value in ROLES.items()}
    qm_match = re.search(
        r"qmmask\s*=\s*['\"]([^'\"]+)['\"]",
        (seed / "inputs/dimer.mdin").read_text(),
        re.IGNORECASE,
    )
    if qm_match is None:
        raise ValueError("Native input has no explicit qmmask")
    qm_original = np.flatnonzero(
        np.asarray(AmberMask(parm, qm_match.group(1)).Selection(), dtype=bool)
    ).tolist()
    if len(qm_original) != 98 or not set(qm_original) <= remap.keys():
        raise ValueError("Expected all 98 real QM atoms in the stripped viewer")
    qm_set = set(qm_original)
    substrate_original = [
        atom.idx
        for atom in parm.atoms
        if atom.residue.name in {"P6N", "P6I", "P6C"} and atom.idx not in qm_set
    ]
    if not substrate_original or not set(substrate_original) <= remap.keys():
        raise ValueError("Expected retained non-QM modular substrate atoms")
    qm_selection = " ".join(str(remap[index]) for index in qm_original)
    qm_carbons = " ".join(
        str(remap[index])
        for index in qm_original
        if parm.atoms[index].atomic_number == 6
    )
    substrate_selection = " ".join(str(remap[index]) for index in substrate_original)
    substrate_carbons = " ".join(
        str(remap[index])
        for index in substrate_original
        if parm.atoms[index].atomic_number == 6
    )
    lines = [
        "set here [file dirname [file normalize [info script]]]",
        "mol new [file join $here stripped.parm7] type parm7 waitfor all",
        "mol addfile [file join $here optimization.dcd] type dcd waitfor all",
        "animate goto 0",
        "mol delrep 0 top",
        "mol representation NewCartoon; mol selection protein; mol color Chain; mol addrep top",
        "mol showrep top 0 off",
        "color Element H white; color Element N blue; color Element O red",
        f"mol representation Licorice; mol selection {{index {qm_selection}}}; mol color Element; mol addrep top",
        f"mol representation Licorice; mol selection {{index {qm_carbons}}}; mol color ColorID 10; mol addrep top",
        f"mol representation Licorice; mol selection {{index {substrate_selection}}}; mol color Element; mol addrep top",
        f"mol representation Licorice; mol selection {{index {substrate_carbons}}}; mol color ColorID 3; mol addrep top",
        f"mol representation VDW 0.28 12; mol selection {{index {mapped['H1']} {mapped['H2']} {mapped['HG']} {mapped['HN']}}}; mol color Element; mol addrep top",
    ]
    lines += [
        "proc relay_remove_bond {a b} {",
        "  foreach {first second} [list $a $b $b $a] {",
        '    set selected [atomselect top "index $first"]',
        "    set neighbors [lindex [$selected getbonds] 0]",
        "    set position [lsearch -exact $neighbors $second]",
        "    if {$position >= 0} {$selected setbonds [list [lreplace $neighbors $position $position]]}",
        "    $selected delete",
        "  }",
        "}",
    ]
    for a, b in [("HG", "OG"), ("H1", "N"), ("H2", "N"), ("HN", "NL"), ("C", "NL")]:
        lines.append(f"relay_remove_bond {mapped[a]} {mapped[b]}")
    alternatives = [
        (h, a, 1.2) for h in ("H1", "H2", "HG", "HN") for a in ("N", "OG", "NL", "O")
    ]
    alternatives += [("C", "OG", 1.85), ("C", "NL", 1.85)]
    for a, b, cutoff in alternatives:
        lines.append(
            f"mol representation DynamicBonds {cutoff} 0.12 12; mol selection {{index {mapped[a]} {mapped[b]}}}; mol color Element; mol addrep top"
        )
    for a, b in [
        ("C", "OG"),
        ("C", "NL"),
        ("C", "O"),
        ("HG", "N"),
        ("H1", "N"),
        ("H2", "N"),
        ("HG", "OG"),
        ("HG", "NL"),
    ]:
        lines.append(
            f"label add Bonds [molinfo top]/{mapped[a]} [molinfo top]/{mapped[b]}"
        )
    for residue, name in ((113, "N"), (186, "ND2")):
        donors = [a.idx for a in parm.residues[residue - 1].atoms if a.name == name]
        if len(donors) != 1 or donors[0] not in remap:
            raise ValueError("Oxyanion-hole donor atom mapping unavailable")
        lines.append(
            f"label add Bonds [molinfo top]/{mapped['O']} [molinfo top]/{remap[donors[0]]}"
        )
    lines += [
        f"label add Angles [molinfo top]/{mapped['OG']} [molinfo top]/{mapped['C']} [molinfo top]/{mapped['O']}",
        "label show Bonds all; label show Angles all",
        f"set visible [atomselect top {{index {qm_selection} {substrate_selection}}} frame 0]",
        "if {[$visible num] == 0} {$visible delete; error {No visible active-site atoms}}",
        "if {[catch {measure minmax $visible} bounds]} {$visible delete; error $bounds}",
        "$visible delete",
        "set center [vecscale 0.5 [vecadd [lindex $bounds 0] [lindex $bounds 1]]]",
        "set span [vecsub [lindex $bounds 1] [lindex $bounds 0]]",
        "set extent [tcl::mathfunc::max {*}$span]",
        "if {$extent <= 0} {error {Visible atom extent must be positive}}",
        "set identity {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}",
        "set centering $identity; set scaling $identity",
        "for {set i 0} {$i < 3} {incr i} {lset centering $i 3 [expr {-[lindex $center $i]}]; lset scaling $i $i [expr {2.6/$extent}]}",
        "molinfo top set center_matrix [list $centering]",
        "molinfo top set rotate_matrix [list $identity]",
        "molinfo top set global_matrix [list $identity]",
        "molinfo top set scale_matrix [list $scaling]; display update",
    ]
    lines.append(
        'puts "Frame0: exact sampled seed. Later frames: optimizer iterations, not MD time or reaction path."'
    )
    (output / "view_relay.vmd").write_text("\n".join(lines) + "\n")
    write_json(
        output / "manifest.json",
        {
            "seed": str(seed),
            "frame_count": len(states),
            "source_optimizer_frames_0based": optimizer_indices,
            "snapshot_optimizer_frames_0based": original_selected,
            "frame_sources": frame_sources,
            "final_restart_comparison": restart_comparison,
            "original_source_frame_count_at_snapshot": original_count,
            "stripped_to_original_atom_indices_0based": kept.tolist(),
            "role_indices_stripped_0based": mapped,
            "qm_atom_indices_original_0based": qm_original,
            "qm_atom_indices_stripped_0based": [remap[index] for index in qm_original],
            "mm_substrate_atom_indices_stripped_0based": [
                remap[index] for index in substrate_original
            ],
            "validation": validation,
            "topology_sha256": digest(topology),
            "source_snapshot_sha256": digest(source),
            "seed_coordinate_sha256": digest(initial),
            "strip_mask": STRIP_MASK,
            "interpretation": "Optimizer states; not MD time, a timed reaction path, or evidence of convergence",
        },
    )
    (output / "README.md").write_text(
        "Source view_relay.vmd in VMD. Frame0 is the exact trajectory seed; later frames are saved native optimizer iterations. Numerical/scientific acceptance is separate.\n"
    )
    (output / "MANIFEST.sha256").write_text(
        "".join(
            f"{digest(p)}  {p.name}\n" for p in sorted(output.iterdir()) if p.is_file()
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    export(args.seed.resolve(), args.output.resolve())
