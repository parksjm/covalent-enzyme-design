"""Build self-contained VMD Tcl for the deacylation hydration movie."""

from __future__ import annotations

from .vmd_movie_policy import render_movie_policy


def viewer_text(
    frame_count: int,
    atom_count: int,
    qm_indices: list[int],
    water_atom_indices: list[int],
    proton_index: int,
    atom_pairs: dict[str, tuple[int, int]],
    angle_indices: tuple[int, int, int],
    center_indices: list[int],
    *,
    substrate_indices: list[int],
    cc_bonds: list[tuple[int, int]],
    chain_groups: dict[str, list[int]],
) -> str:
    """Return a portable viewer script with mapped reactive bonds and labels."""
    qm = " ".join(str(value) for value in qm_indices)
    waters = " ".join(str(value) for value in water_atom_indices)
    pair_commands = "\n".join(
        f"removebond $top {first} {second}\npairbond $top {first} {second} {cutoff}"
        for (first, second), cutoff in (
            (atom_pairs["OG1_C199"], 1.85),
            (atom_pairs["C199_Nleave"], 1.85),
            (atom_pairs["OG1_HG1"], 1.4),
            (atom_pairs["Nalpha_HG1"], 1.4),
            (atom_pairs["Nleave_HG1"], 1.4),
        )
    )
    labels = "\n".join(
        f"label add Bonds $top/{a} $top/{b}"
        for a, b in (
            atom_pairs["OG1_C199"],
            atom_pairs["C199_Nleave"],
            atom_pairs["O199_Phe113N"],
            atom_pairs["O199_Asn186ND2"],
            atom_pairs["candidate_O_C199"],
            atom_pairs["candidate_O_Nalpha"],
        )
    )
    angle = " ".join(f"$top/{value}" for value in angle_indices)
    center = " ".join(str(value) for value in center_indices)
    policy = render_movie_policy(
        qm_indices,
        substrate_indices,
        cc_bonds,
        molecule="$top",
        chain_groups=chain_groups,
    )
    return f"""# Portable stripped deacylation hydration movie; run: vmd -e view_hydration.vmd
set here [file dirname [file normalize [info script]]]
set top [mol new [file join $here stripped.parm7] type parm7 waitfor all]
mol addfile [file join $here hydration.dcd] type dcd waitfor all molid $top
if {{[molinfo $top get numframes] != {frame_count}}} {{error "Unexpected frame count"}}
if {{[molinfo $top get numatoms] != {atom_count}}} {{error "Unexpected atom count"}}
mol delrep 0 $top
color Element H white
color Element N blue
color Element O red
color change rgb 3 1.0 0.5 0.0
mol representation NewCartoon
mol color Structure
mol selection protein
mol addrep $top
mol showrep $top 0 off
mol representation Licorice 0.18 12 12
mol color Element
mol selection "resname P6N P6I P6C"
mol addrep $top
mol representation Licorice 0.28 12 12
mol color ColorID 3
mol selection "resname P6N P6I P6C and element C"
mol addrep $top
mol representation Licorice 0.22 12 12
mol color Element
mol selection "index {qm}"
mol addrep $top
mol representation Licorice 0.18 12 12
mol color Element
mol selection "index {waters}"
mol addrep $top
mol representation VDW 0.28 12
mol color Element
mol selection "index {proton_index}"
mol addrep $top
# Remove labile topology edges in the viewer only; the native Amber topology is unchanged.
proc removebond {{molid a b}} {{
    foreach {{first second}} [list $a $b $b $a] {{
        set selected [atomselect $molid "index $first"]
        set neighbors [lindex [$selected getbonds] 0]
        set position [lsearch -exact $neighbors $second]
        if {{$position >= 0}} {{$selected setbonds [list [lreplace $neighbors $position $position]]}}
        $selected delete
    }}
}}
# Pair-specific dynamic representations redraw reactive bonds at the stated cutoffs.
proc pairbond {{mol a b cutoff}} {{
    mol representation DynamicBonds $cutoff 0.1 12
    mol color Element
    mol selection "index $a $b"
    mol addrep $mol
    set rep [expr {{[molinfo $mol get numreps] - 1}}]
    mol selupdate $rep $mol on
}}
{pair_commands}
{policy}
mol rename $top "PA6 deacylation hydration reconnaissance"
{labels}
label add Angles {angle}
display depthcue off
display resetview
set focus [atomselect $top "index {center}"]
molinfo $top set center_matrix [trans origin [measure center $focus]]
$focus delete
scale by 2.5
animate goto 0
puts "Loaded {frame_count} chronological frames; initial endpoint, then heating and equilibration."
"""
