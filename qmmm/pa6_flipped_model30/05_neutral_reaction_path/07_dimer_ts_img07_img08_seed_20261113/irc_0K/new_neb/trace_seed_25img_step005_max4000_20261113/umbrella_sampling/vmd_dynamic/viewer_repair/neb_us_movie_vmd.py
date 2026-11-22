"""Portable VMD rendering for concatenated NEB umbrella trajectories."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from nyl12_covalent.relay_progress_vmd import substrate_indices

_REQUIRED_ROLES = {
    "carbon", "og1", "leaving_n", "hg1", "thr_n", "oxygen", "phe_n", "asn_nd2"
}


def _numbers(indices: Iterable[int]) -> str:
    """Return sorted, space-separated VMD atom indices."""
    return " ".join(str(index) for index in sorted(set(indices)))


def _validate_mapping(
    topology: Any,
    original_atom_ids_1based: list[int],
    qm_atom_ids_1based: list[int],
    roles_1based: dict[str, int],
    frames: int,
) -> dict[int, int]:
    """Validate source identities and build original-to-stripped index map."""
    atom_count = len(topology.atoms)
    kept = list(original_atom_ids_1based)
    if not kept or kept != sorted(set(kept)):
        raise ValueError("Kept original atom IDs must be nonempty, unique, and sorted")
    if kept[0] < 1 or kept[-1] > atom_count:
        raise ValueError("Kept original atom ID is outside the topology")
    if frames < 1:
        raise ValueError("Expected frame count must be positive")
    qm = set(qm_atom_ids_1based)
    if not qm or not qm <= set(kept):
        raise ValueError("QM atom IDs must be nonempty and included in kept atoms")
    missing = sorted(_REQUIRED_ROLES - roles_1based.keys())
    if missing:
        raise ValueError(f"Missing required role mappings: {', '.join(missing)}")
    for role, atom_id in roles_1based.items():
        if not isinstance(atom_id, int) or atom_id < 1 or atom_id > atom_count:
            raise ValueError(f"Role {role!r} has an invalid original atom ID")
        if atom_id not in kept:
            raise ValueError(f"Role {role!r} atom must be included in kept atoms")
    return {original_id: index for index, original_id in enumerate(kept)}


def render_movie_vmd(
    topology: Any,
    original_atom_ids_1based: list[int],
    qm_atom_ids_1based: list[int],
    roles_1based: dict[str, int],
    frames: int,
) -> str:
    """Render a portable, mapped VMD script for a NEB-US trajectory.

    The input trajectory is expected to contain concatenated umbrella windows;
    it is a biased-window movie, not a continuous reaction trajectory.

    Args:
        topology: Original full ParmEd topology with atom identities and bonds.
        original_atom_ids_1based: Sorted original IDs retained in stripped files.
        qm_atom_ids_1based: Original one-based IDs in the QM region.
        roles_1based: Explicit role-to-original-one-based-atom mapping.
        frames: Expected number of frames in the local path.dcd file.

    Returns:
        Tcl script that loads local path.parm7 and path.dcd files.

    Raises:
        ValueError: If atom mappings or required role identities are invalid.
    """
    mapping = _validate_mapping(
        topology, original_atom_ids_1based, qm_atom_ids_1based, roles_1based, frames
    )
    substrate = substrate_indices(
        topology, roles_1based["carbon"] - 1, roles_1based["leaving_n"] - 1
    )
    qm_original = set(qm_atom_ids_1based)
    qm_local = {mapping[atom_id] for atom_id in qm_original}
    substrate_local = {mapping[index + 1] for index in substrate}
    nonqm_substrate_carbons = {
        mapping[index + 1]
        for index in substrate
        if index + 1 not in qm_original
        and topology.atoms[index].atomic_number == 6
    }
    reactive = {
        frozenset((roles_1based["carbon"], roles_1based["og1"])),
        frozenset((roles_1based["carbon"], roles_1based["leaving_n"])),
    }
    removed_bonds: set[tuple[int, int]] = set()
    for atom_id in sorted(qm_original):
        atom = topology.atoms[atom_id - 1]
        if atom.atomic_number == 1:
            for partner in atom.bond_partners:
                partner_id = partner.idx + 1
                if partner_id in qm_original and partner.atomic_number in {7, 8}:
                    removed_bonds.add(tuple(sorted((mapping[atom_id], mapping[partner_id]))))
        for partner in atom.bond_partners:
            partner_id = partner.idx + 1
            if frozenset((atom_id, partner_id)) in reactive:
                removed_bonds.add(tuple(sorted((mapping[atom_id], mapping[partner_id]))))

    atom = {name: mapping[value] for name, value in roles_1based.items()}
    qm_n_o = {
        mapping[atom_id]
        for atom_id in qm_original
        if topology.atoms[atom_id - 1].atomic_number in {7, 8}
    }
    qm_h = [
        mapping[atom_id]
        for atom_id in sorted(qm_original)
        if topology.atoms[atom_id - 1].atomic_number == 1
    ]
    visible = qm_local | substrate_local
    lines = [
        "set title {NEB-US movie}",
        "foreach filename {path.parm7 path.dcd} {",
        ' if {![file readable [file join $here $filename]]} {error "Missing readable movie file: [file join $here $filename]"}',
        "}",
        "proc ::neb_us_movie::remove {molid} {",
        " foreach kind {Atoms Bonds Angles Dihedrals} {",
        "  set labels [label list $kind]",
        "  for {set i [expr {[llength $labels]-1}]} {$i >= 0} {incr i -1} {",
        "   set owned 0",
        "   foreach spec [lindex $labels $i] {if {[llength $spec] == 2 && [lindex $spec 0] eq $molid} {set owned 1}}",
        "   if {$owned} {label delete $kind $i}",
        "  }",
        " }",
        " mol delete $molid",
        "}",
        "set previous $::neb_us_movie::molecule",
        "if {$previous ne {} && [lsearch -exact [molinfo list] $previous] >= 0} {::neb_us_movie::remove $previous}",
        "set molid [mol new [file join $here path.parm7] type parm7 waitfor all]",
        "mol rename $molid $title",
        "mol addfile [file join $here path.dcd] type dcd waitfor all molid $molid",
        'if {[molinfo $molid get numframes] != EXPECTED_FRAMES} {error {Unexpected trajectory frame count}}'.replace("EXPECTED_FRAMES", str(frames)),
        f'if {{[molinfo $molid get numatoms] != {len(mapping)}}} {{error {{Unexpected topology atom count}}}}',
        "molinfo $molid set frame 0",
        'if {[molinfo $molid get frame] != 0} {error {Could not activate trajectory frame zero}}',
        "set ::neb_us_movie::molecule $molid",
        "mol delrep 0 $molid",
        "color Element H white; color Element N blue; color Element O red; label textsize 0.6",
        "color change rgb 3 1.0 0.5 0.0",
        "mol representation NewCartoon; mol color Structure; mol selection protein; mol addrep $molid; mol showrep $molid 0 off",
        f"mol representation Licorice 0.18 12 12; mol color Element; mol selection {{index {_numbers(visible)}}}; mol addrep $molid",
        f"mol representation Licorice 0.25 12 12; mol color ColorID 3; mol selection {{index {_numbers(nonqm_substrate_carbons)}}}; mol addrep $molid",
        "set sel [atomselect $molid all]; set bonds [$sel getbonds]",
        "foreach pair {" + " ".join(f"{{{a} {b}}}" for a, b in sorted(removed_bonds)) + "} {",
        " set a [lindex $pair 0]; set b [lindex $pair 1]",
        " set row [lindex $bonds $a]; set pos [lsearch -exact $row $b]; if {$pos >= 0} {set row [lreplace $row $pos $pos]}; lset bonds $a $row",
        " set row [lindex $bonds $b]; set pos [lsearch -exact $row $a]; if {$pos >= 0} {set row [lreplace $row $pos $pos]}; lset bonds $b $row",
        "}",
        "$sel setbonds $bonds; $sel delete",
        f"foreach pair {{{{{atom['og1']} {atom['carbon']}}} {{{atom['carbon']} {atom['leaving_n']}}}}} {{mol representation DynamicBonds 1.85 0.1 12; mol color Element; mol selection [format {{index %s %s}} [lindex $pair 0] [lindex $pair 1]]; mol addrep $molid}}",
    ]
    for hydrogen in qm_h:
        lines.append(
            f"mol representation DynamicBonds 1.4 0.1 12; mol color Element; mol selection {{index {hydrogen} {_numbers(qm_n_o)}}}; mol addrep $molid"
        )
    lines.extend(
        [
            f"mol representation VDW 0.28 12; mol color Element; mol selection {{index {atom['hg1']}}}; mol addrep $molid",
            "proc ::neb_us_movie::add_label {kind ids offset} {",
            " set molid $::neb_us_movie::molecule; set specs {}",
            " foreach id $ids {lappend specs $molid/$id}",
            " ::label add $kind {*}$specs",
            " set index [expr {[llength [::label list $kind]]-1}]",
            " ::label show $kind $index; ::label textoffset $kind $index $offset",
            "}",
        ]
    )
    pairs = [
        ("og1", "hg1"), ("hg1", "leaving_n"), ("hg1", "thr_n"),
        ("og1", "carbon"), ("carbon", "leaving_n"),
        ("oxygen", "phe_n"), ("oxygen", "asn_nd2"),
    ]
    offsets = [
        (0.04, 0.08), (0.04, -0.08), (-0.13, 0.07), (0.08, -0.05),
        (-0.12, -0.09), (0.06, 0.13), (-0.11, 0.14),
    ]
    for (first, second), (dx, dy) in zip(pairs, offsets, strict=True):
        lines.append(
            f"::neb_us_movie::add_label Bonds [list {atom[first]} {atom[second]}] {{{dx} {dy}}}"
        )
    lines.append(
        f"::neb_us_movie::add_label Angles [list {atom['og1']} {atom['carbon']} {atom['oxygen']}] {{0.04 -0.14}}"
    )
    lines.extend(
        [
            f"set visible [atomselect $molid {{index {_numbers(visible)}}} frame 0]",
            'if {[$visible num] == 0} {$visible delete; error {No visible active-site atoms}}',
            'if {[catch {measure minmax $visible} bounds]} {$visible delete; error "Cannot measure movie frame zero: $bounds"}',
            "$visible delete",
            "set center [vecscale 0.5 [vecadd [lindex $bounds 0] [lindex $bounds 1]]]",
            "set span [vecsub [lindex $bounds 1] [lindex $bounds 0]]; set extent [tcl::mathfunc::max {*}$span]",
            "if {$extent <= 0} {error {Visible atom extent must be positive}}",
            "set identity {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}; set centering $identity; set scaling $identity",
            "for {set i 0} {$i < 3} {incr i} {lset centering $i 3 [expr {-[lindex $center $i]}]; lset scaling $i $i [expr {2.6/$extent}]}",
            "molinfo $molid set center_matrix [list $centering]; molinfo $molid set rotate_matrix [list $identity]; molinfo $molid set global_matrix [list $identity]; molinfo $molid set scale_matrix [list $scaling]",
            "mol top $molid; mol on $molid; display update",
            "return $molid",
        ]
    )
    header = [
        "# Concatenated biased windows; not a continuous reaction trajectory.",
        "namespace eval ::neb_us_movie {variable molecule; if {![info exists molecule]} {set molecule {}}}",
        "proc ::neb_us_movie::load {here} {",
    ]
    return "\n".join(
        header + lines + ["}", "::neb_us_movie::load [file dirname [file normalize [info script]]]"]
    ) + "\n"
