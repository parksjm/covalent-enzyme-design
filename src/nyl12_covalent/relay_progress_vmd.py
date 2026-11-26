"""Portable, molecule-scoped VMD views of accepted relay snapshots."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

from .vmd_chain_groups import display_chain_groups
from .vmd_movie_policy import mapped_cc_bonds, render_movie_policy


VMD_API_REFERENCES = (
    "https://www.ks.uiuc.edu/Research/vmd/doxygen/cmd__label_8C-source.html",
    "https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node140.html",
    "https://www.ks.uiuc.edu/Research/vmd/doxygen/DrawMolItem_8C-source.html",
)


def substrate_indices(topology: Any, carbon: int, leaving_n: int) -> set[int]:
    """Find the mapped substrate component using zero-based bonded atom indices.

    Args:
        topology: ParmEd topology with atom bond partners and residue names.
        carbon: Mapped scissile carbon index.
        leaving_n: Mapped leaving nitrogen index.

    Returns:
        Connected P6N/P6I/P6C atoms, excluding protein and other polymer molecules.

    Raises:
        ValueError: If mapped reactive atoms are not in one polymer component.
    """
    allowed = {a.idx for a in topology.atoms if a.residue.name in {"P6N", "P6I", "P6C"}}
    if carbon not in allowed or leaving_n not in allowed:
        raise ValueError("Mapped scissile atoms must belong to P6N/P6I/P6C substrate")
    found, pending = set(), [carbon]
    while pending:
        current = pending.pop()
        if current in found:
            continue
        found.add(current)
        pending.extend(
            a.idx
            for a in topology.atoms[current].bond_partners
            if a.idx in allowed and a.idx not in found
        )
    if leaving_n not in found:
        raise ValueError(
            "Mapped carbon and leaving N are in different substrate components"
        )
    return found


# Helpers never change visibility or labels belonging to another molecule.
_SCOPED_TCL = r"""
namespace eval ::relay_progress {variable molecules; if {![info exists molecules]} {set molecules {}}}
proc ::relay_progress::remove {molid} {
    foreach kind {Atoms Bonds Angles Dihedrals} {
        set labels [label list $kind]
        for {set i [expr {[llength $labels]-1}]} {$i >= 0} {incr i -1} {
            set owned 0
            foreach spec [lindex $labels $i] {
                if {[llength $spec] == 2 && [lindex $spec 0] eq $molid} {set owned 1}
            }
            if {$owned} {label delete $kind $i}
        }
    }
    mol delete $molid
}
proc ::relay_progress::add_label {kind atoms {offset {0.03 0.03}}} {
    set index 0
    foreach item [label list $kind] {
        if {[lrange $item 0 [expr {[llength $atoms]-1}]] eq $atoms} {
            label show $kind $index
            label textoffset $kind $index $offset
            return $index
        }
        incr index
    }
    set specs {}
    foreach atom $atoms {lappend specs "[lindex $atom 0]/[lindex $atom 1]"}
    label add $kind {*}$specs
    label show $kind $index
    label textoffset $kind $index $offset
    return $index
}
"""


def vmd_script(
    manifest: dict[str, Any],
    case: dict[str, Any],
    topology: Any,
    *,
    distance_pairs: Sequence[tuple[str, str]] | None = None,
) -> str:
    """Render mapped substrate/QM atoms, route-specific labels and a local camera.

    Args:
        manifest: Bundle's explicit stripped-to-original zero-based index map.
        case: Original one-based QM atom list and named reactive atom mapping.
        topology: Original ParmEd topology with residue identities/connectivity.
        distance_pairs: Explicit route atom-map keys; None uses six acyl-route labels.

    Returns:
        Portable Tcl; repeated sourcing replaces only this campaign's molecule.
    """
    mapping = {
        old + 1: new
        for new, old in enumerate(manifest["stripped_to_original_atom_indices_0based"])
    }
    original_qm = set(map(int, case["qm_atom_ids"]))
    substrate = substrate_indices(
        topology,
        int(case["atom_map"]["carbon"]) - 1,
        int(case["atom_map"]["leaving_n"]) - 1,
    )
    qm = [mapping[i] for i in sorted(original_qm)]
    hydrogens = [
        mapping[i]
        for i in sorted(original_qm)
        if topology.atoms[i - 1].atomic_number == 1
    ]
    on_indices = [
        mapping[i]
        for i in sorted(original_qm)
        if topology.atoms[i - 1].atomic_number in (7, 8)
    ]
    nonqm_c = {
        mapping[i + 1]
        for i in substrate
        if topology.atoms[i].atomic_number == 6 and i + 1 not in original_qm
    }
    visible = sorted(set(qm) | {mapping[i + 1] for i in substrate})
    element_atoms = visible
    atom = {key: mapping[int(value)] for key, value in case["atom_map"].items()}

    def numbers(ids: Iterable[int]) -> str:
        return " ".join(map(str, sorted(ids)))

    lines = [
        "# Accepted diagnostic states; not dynamics or a validated relay.",
        "set here [file dirname [file normalize [info script]]]",
        _SCOPED_TCL,
        "set owner [file dirname $here]",
        "if {[regexp {^update_[0-9]+$} [file tail $owner]]} {set owner [file dirname $owner]}",
        'set title "Relay campaign: $owner"',
        "set previous {}",
        "if {[dict exists $::relay_progress::molecules $owner]} {",
        " set registered [dict get $::relay_progress::molecules $owner]",
        " if {[lsearch -exact [molinfo list] $registered] >= 0 && [lindex [molinfo $registered get name] 0] eq $title} {set previous $registered}",
        "}",
        "if {$previous eq {}} {",
        " set expected_top [file normalize [file join $here path.parm7]]",
        " set expected_traj [file normalize [file join $here path.dcd]]",
        " set matches {}",
        " foreach candidate [molinfo list] {",
        "  set source_files {}",
        "  foreach group [lindex [molinfo $candidate get filename] 0] {",
        "   if {[file extension $group] in {.parm7 .dcd}} {set group [list $group]}",
        "   foreach filename $group {lappend source_files [file normalize $filename]}",
        "  }",
        "  if {[lsearch -exact $source_files $expected_top] >= 0 && [lsearch -exact $source_files $expected_traj] >= 0} {lappend matches $candidate}",
        " }",
        " if {[llength $matches] > 1} {error {Multiple exact topology+trajectory matches; remove duplicate current-bundle molecules explicitly before sourcing this script}}",
        " if {[llength $matches] == 1} {set previous [lindex $matches 0]}",
        "}",
        "if {$previous ne {}} {::relay_progress::remove $previous}",
        "set molid [mol new [file join $here path.parm7] type parm7 waitfor all]",
        "mol rename $molid $title",
        "dict set ::relay_progress::molecules $owner $molid",
        "mol addfile [file join $here path.dcd] type dcd waitfor all molid $molid",
        "mol delrep 0 $molid",
        f"set proton_indices {{{numbers(hydrogens)}}}",
        f"set hetero_indices {{{numbers(on_indices)}}}",
        f"set water_hydrogens {{{atom['water_h1']} {atom['water_h2']}}}",
        f"set attack_o {atom['thr_og']}; set scissile_c {atom['carbon']}; set scissile_n {atom['leaving_n']}",
        "set sel [atomselect $molid all]",
        "set bonds [$sel getbonds]",
        "for {set i 0} {$i < [llength $bonds]} {incr i} {",
        " set keep {}",
        " foreach j [lindex $bonds $i] {",
        "  set labileH [expr {([lsearch -exact $proton_indices $i] >= 0 && [lsearch -exact $hetero_indices $j] >= 0) || ([lsearch -exact $proton_indices $j] >= 0 && [lsearch -exact $hetero_indices $i] >= 0)}]",
        "  set reactive [expr {($i == $scissile_c && ($j == $scissile_n || $j == $attack_o)) || ($j == $scissile_c && ($i == $scissile_n || $i == $attack_o))}]",
        "  set waterHH [expr {[lsearch -exact $water_hydrogens $i] >= 0 && [lsearch -exact $water_hydrogens $j] >= 0}]",
        "  if {!$labileH && !$reactive && !$waterHH} {lappend keep $j}",
        " }",
        " lset bonds $i $keep",
        "}",
        "$sel setbonds $bonds",
        "$sel delete",
        "color Element H white",
        "color Element N blue",
        "color Element O red",
        "color change rgb 3 1.0 0.5 0.0",
        "mol representation NewCartoon; mol color Structure; mol selection protein; mol addrep $molid",
        "mol showrep $molid 0 off",
    ]
    selections = [("Licorice 0.18 12 12", "Element", numbers(element_atoms))]
    if nonqm_c:
        selections.append(("Licorice 0.25 12 12", "ColorID 3", numbers(nonqm_c)))
    selections.extend(
        ("DynamicBonds 1.85 0.1 12", "Element", f"{atom[a]} {atom[b]}")
        for a, b in (("thr_og", "carbon"), ("carbon", "leaving_n"))
    )
    selections.extend(
        ("DynamicBonds 1.4 0.1 12", "Element", f"{h} {numbers(on_indices)}")
        for h in hydrogens
    )
    for representation, color, selection in selections:
        lines.extend(
            [
                f"mol representation {representation}",
                f"mol color {color}",
                f'mol selection "index {selection}"',
                "mol addrep $molid",
            ]
        )
    defaults = [
        ("thr_og", "carbon"),
        ("carbon", "leaving_n"),
        ("thr_og", "water_h1"),
        ("leaving_n", "water_h1"),
        ("oxygen", "phe_n"),
        ("oxygen", "asn_nd2"),
    ]
    if distance_pairs is not None:
        defaults = list(distance_pairs)
    if not defaults or any(
        a not in atom or b not in atom or a == b for a, b in defaults
    ):
        raise ValueError("Distance labels require nonempty distinct mapped atom pairs")
    offsets = [
        (0.04, 0.08),
        (0.04, -0.08),
        (-0.13, 0.07),
        (0.08, -0.05),
        (-0.12, -0.09),
        (0.06, 0.13),
        (-0.11, 0.14),
    ]
    for index, (first, second) in enumerate(defaults):
        dx, dy = offsets[index % len(offsets)]
        lines.append(
            f"::relay_progress::add_label Bonds [list [list $molid {atom[first]}] [list $molid {atom[second]}]] {{{dx} {dy}}}"
        )
    lines.append(
        f"::relay_progress::add_label Angles [list [list $molid {atom['thr_og']}] [list $molid {atom['carbon']}] [list $molid {atom['oxygen']}]] {{0.04 -0.14}}"
    )
    network = [
        ("carbon", "oxygen"),
        ("thr_og", "thr_n"),
        ("thr_n", "water_o"),
        ("water_o", "leaving_n"),
        ("thr_og", "proton"),
        ("thr_n", "proton"),
        ("water_o", "water_h1"),
        ("water_o", "water_h2"),
        ("thr_n", "water_h1"),
    ]
    lines += [
        "# Optional: ::relay_progress::show_network $molid",
        "proc ::relay_progress::show_network {molid} {",
    ]
    for first, second in network:
        lines.append(
            f" ::relay_progress::add_label Bonds [list [list $molid {atom[first]}] [list $molid {atom[second]}]]"
        )
    lines += [
        "}",
        f"set visible [atomselect $molid {{index {numbers(visible)}}}]",
        "set bounds [measure minmax $visible]",
        "$visible delete",
        "set center [vecscale 0.5 [vecadd [lindex $bounds 0] [lindex $bounds 1]]]",
        "set span [vecsub [lindex $bounds 1] [lindex $bounds 0]]",
        "set extent [tcl::mathfunc::max {*}$span]",
        "if {$extent <= 0} {error {Visible atom extent must be positive}}",
        "set identity {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}",
        "set centering $identity",
        "for {set i 0} {$i < 3} {incr i} {lset centering $i 3 [expr {-[lindex $center $i]}]}",
        "set scaling $identity",
        "for {set i 0} {$i < 3} {incr i} {lset scaling $i $i [expr {2.6/$extent}]}",
        "molinfo $molid set center_matrix [list $centering]",
        "molinfo $molid set rotate_matrix [list $identity]",
        "molinfo $molid set global_matrix [list $identity]",
        "molinfo $molid set scale_matrix [list $scaling]",
    ]
    lines.append(
        render_movie_policy(
            qm,
            [mapping[index + 1] for index in substrate],
            mapped_cc_bonds(topology, {old - 1: new for old, new in mapping.items()}),
            chain_groups=display_chain_groups(topology, mapping),
        )
    )
    return "\n".join(lines) + "\n"
