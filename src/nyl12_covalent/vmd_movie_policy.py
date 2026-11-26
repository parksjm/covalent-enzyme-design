"""Portable, display-only VMD styling for physical reaction-movie atoms."""

from collections.abc import Iterable, Mapping
from numbers import Integral
import re
from typing import Any


POLICY_VERSION = "20261126-typeprefix3"


_POLICY = r"""
# Display policy version: 20261126-typeprefix3
namespace eval ::nyl12_movie_policy { variable original_types }
proc ::nyl12_movie_policy::palette {} {
    # VMD's 33 named solid colors include tan; brown is not a valid VMD color.
    # VMD Type coloring uses the first nonnumeric character, not the full type.
    set palette {C cyan Q orange H white N blue O red S yellow P tan F green I purple Z silver}
    set available [colorinfo colors]
    foreach {kind shade} $palette {
        if {$shade ni $available} {error "Unsupported VMD color '$shade' for $kind"}
    }
    return $palette
}
proc ::nyl12_movie_policy::restore_types {m} {
    variable original_types
    if {![info exists original_types($m)]} {error "No original types for molecule $m"}
    set all [atomselect $m all]
    $all set type $original_types($m)
    $all delete
    unset original_types($m)
}
proc ::nyl12_movie_policy::apply {m qm substrate cc_pairs chain_groups} {
    variable original_types
    set palette [::nyl12_movie_policy::palette]
    set all [atomselect $m all]
    set indices [$all get index]
    set elements [$all get element]
    if {![info exists original_types($m)]} {set original_types($m) [$all get type]}
    set visible [lsort -integer -unique [concat $qm $substrate]]
    foreach i $visible {
        if {[lsearch -exact $indices $i] < 0} {error "Invalid physical atom index $i"}
    }
    array set element {}
    array set sub {}
    foreach i $substrate {set sub($i) 1}
    set buckets {C C SUB_C Q H H N N O O S S P P F F CL F BR P I I NA N MG F CA Z ZN Z FE Q}
    set types {}
    foreach i $indices e $elements {
        set element($i) [string toupper $e]
        set kind $element($i)
        if {$element($i) eq "C" && [info exists sub($i)]} {set kind SUB_C}
        if {![dict exists $buckets $kind]} {error "Unsupported display element '$e' at atom $i"}
        lappend types "[dict get $buckets $kind]_NYL12_MOVIE_$kind"
    }
    # Type is display-only; original Amber types remain available to restore_types.
    $all set type $types
    foreach {kind shade} $palette {
        if {[lsearch -exact [colorinfo category Type] $kind] < 0} {
            color add item Type $kind $shade
        }
        color Type $kind $shade
    }
    set anchors [atomselect $m "protein and name CA CB and same residue as (index [join $qm { }])"]
    set visible [lsort -integer -unique [concat $visible [$anchors get index]]]
    # Only named CA/CB pairs within the same protein residue supplement supplied edges.
    array set ca {}
    array set cb {}
    foreach i [$anchors get index] r [$anchors get residue] n [$anchors get name] {
        if {$n eq "CA"} {set ca($r) $i}
        if {$n eq "CB"} {set cb($r) $i}
    }
    foreach r [array names ca] {
        if {[info exists cb($r)]} {lappend cc_pairs [list $ca($r) $cb($r)]}
    }
    $anchors delete
    array set bonds {}
    foreach i $indices neighbors [$all getbonds] {set bonds($i) $neighbors}
    foreach pair $cc_pairs {
        lassign $pair a b
        if {![info exists element($a)] || ![info exists element($b)]} {
            error "C-C pair references a nonphysical atom: $pair"
        }
        if {$element($a) ne "C" || $element($b) ne "C" || $a == $b} {
            error "Expected a physical C-C pair: $pair"
        }
        if {$a ni $visible || $b ni $visible} {continue}
        foreach i [list $a $b] j [list $b $a] {
            if {$j ni $bonds($i)} {lappend bonds($i) $j}
        }
    }
    set repaired {}
    foreach i $indices {lappend repaired $bonds($i)}
    $all setbonds $repaired
    $all delete
    # Preserve original water/pocket selections, including their frame updates.
    set static_selections [list "index [join $visible { }]"]
    # Remove split-color static reps; a single selection colors both bond halves.
    for {set r [expr {[molinfo $m get numreps] - 1}]} {$r >= 0} {incr r -1} {
        set style [lindex [molinfo $m get [list "rep $r"]] 0]
        set name [lindex $style 0]
        if {$name in {Licorice NewCartoon}} {
            if {$name eq "Licorice"} {
                lappend static_selections [lindex [molinfo $m get [list "selection $r"]] 0]
            }
            mol delrep $r $m
        } elseif {$name eq "DynamicBonds"} {
            mol modstyle $r $m DynamicBonds [lindex $style 1] 0.2 [lindex $style 3]
            mol modcolor $r $m Type
        } elseif {$name eq "VDW"} {
            mol modstyle $r $m VDW 0.2 [lindex $style 2]
            mol modcolor $r $m Type
        }
    }
    mol representation Licorice 0.2 12 12
    mol selection "([join $static_selections {) or (}])"
    mol color Type
    mol addrep $m
    mol selupdate [expr {[molinfo $m get numreps] - 1}] $m on
    foreach chain {A B C D E F G H} shade {7 10 27 4 7 10 27 4} {
        set cartoon_selection "protein and chain $chain"
        if {[dict exists $chain_groups $chain]} {
            set group [dict get $chain_groups $chain]
            foreach i $group {
                if {$i ni $indices} {error "Invalid cartoon physical atom index $i"}
            }
            set cartoon_selection "protein and index [join $group { }]"
        }
        set sel [atomselect $m $cartoon_selection]
        set count [$sel num]
        $sel delete
        if {$count == 0} {continue}
        mol representation NewCartoon 0.3 10 4.1 0
        mol selection $cartoon_selection
        mol color ColorID $shade
        mol addrep $m
        mol showrep $m [expr {[molinfo $m get numreps] - 1}] off
    }
    display depthcue off
    set ::nyl12_movie_policy::completed($m) 20261126-typeprefix3
    puts "NYL12 DISPLAY OK version=20261126-typeprefix3 molecule=$m qm=[llength $qm] substrate=[llength $substrate] visible=[llength $visible] radius=0.2"
}
"""


def render_movie_policy(
    qm_indices: Iterable[int],
    substrate_indices: Iterable[int],
    cc_bonds: Iterable[tuple[int, int]] = (),
    *,
    molecule: str = "$molid",
    chain_groups: Mapping[str, Iterable[int]] | None = None,
) -> str:
    """Render self-contained Tcl applying the common movie display policy.

    Args:
        qm_indices: Zero-based physical QM atom indices in the displayed topology.
        substrate_indices: All substrate atom indices in that topology.
        cc_bonds: Known original topology C-C edges; never inferred by distance.
        molecule: Integer molecule ID or simple Tcl variable reference.
        chain_groups: Validated A-H physical atom groups for missing native chain
            IDs. Only display selections change; native chain IDs are preserved.

    Returns:
        Inline Tcl; no installed helper file is needed by an exported movie.

    Raises:
        ValueError: If indices, bonds, or the molecule expression are invalid.

    Notes:
        Atom types change only inside VMD. Original types are retained in
        ``::nyl12_movie_policy::original_types`` and can be restored with
        ``::nyl12_movie_policy::restore_types $molid``. No files are written.
    """

    def checked(values: Iterable[int]) -> list[int]:
        result = list(values)
        if any(
            isinstance(i, bool) or not isinstance(i, Integral) or i < 0 for i in result
        ):
            raise ValueError("Expected nonnegative integer physical atom indices")
        return sorted({int(i) for i in result})

    qm = checked(qm_indices)
    substrate = checked(substrate_indices)
    if not qm or not substrate:
        raise ValueError("QM and substrate selections must both be nonempty")
    if molecule != "[molinfo top]" and not re.fullmatch(
        r"(?:[0-9]+|\$(?:::)?[A-Za-z_]\w*(?:::\w+)*)", molecule
    ):
        raise ValueError("Molecule must be an integer ID or simple Tcl variable")
    pairs = set()
    for pair in cc_bonds:
        values = list(pair)
        if len(values) != 2 or len(checked(values)) != 2:
            raise ValueError("Each C-C bond needs two distinct physical indices")
        pairs.add(tuple(checked(values)))
    qm_tcl = " ".join(map(str, qm))
    substrate_tcl = " ".join(map(str, substrate))
    pairs_tcl = " ".join("{" + f"{a} {b}" + "}" for a, b in sorted(pairs))
    groups = []
    for chain, indices in (chain_groups or {}).items():
        if chain not in tuple("ABCDEFGH"):
            raise ValueError("Cartoon groups must use known chain labels A-H")
        group = checked(indices)
        if not group:
            raise ValueError("Cartoon groups must be nonempty")
        groups.append(chain + " {" + " ".join(map(str, group)) + "}")
    groups_tcl = " ".join(groups)
    return (
        _POLICY
        + "\nif {[catch {\n"
        + f"    ::nyl12_movie_policy::apply {molecule} "
        + f"{{{qm_tcl}}} {{{substrate_tcl}}} {{{pairs_tcl}}} {{{groups_tcl}}}\n"
        + "} ::nyl12_movie_policy::message ::nyl12_movie_policy::options]} {\n"
        + '    puts stderr "NYL12 DISPLAY FAILED version=20261126-typeprefix3: $::nyl12_movie_policy::message"\n'
        + "    return -options $::nyl12_movie_policy::options $::nyl12_movie_policy::message\n"
        + "}\n"
    )


def render_palette_preflight() -> str:
    """Render the palette definition and validation without applying styling.

    Returns:
        Tcl that defines and validates the palette2 colors before scene edits.
    """
    marker = "proc ::nyl12_movie_policy::restore_types"
    definitions, separator, _ = _POLICY.partition(marker)
    if not separator:
        raise ValueError("Palette policy definition boundary is missing")
    return definitions + "\n::nyl12_movie_policy::palette\n"


def mapped_cc_bonds(topology: Any, mapping: Mapping[int, int]) -> list[tuple[int, int]]:
    """Map retained original carbon-carbon topology bonds into movie indices.

    Args:
        topology: ParmEd-compatible structure exposing atoms and bond_partners.
        mapping: Original zero-based index to displayed zero-based index.

    Returns:
        Unique sorted physical C-C bonds; omitted atoms never create edges.

    Raises:
        ValueError: If the mapping contains invalid or duplicate indices.
    """
    keys, values = list(mapping), list(mapping.values())
    if any(
        isinstance(i, bool) or not isinstance(i, Integral) or i < 0
        for i in keys + values
    ):
        raise ValueError("Mapping requires nonnegative integer physical indices")
    if len(set(values)) != len(values):
        raise ValueError("Mapping must be one-to-one")
    mapping = {int(key): int(value) for key, value in mapping.items()}
    pairs = set()
    for original, displayed in mapping.items():
        atom = topology.atoms[original]
        if atom.atomic_number != 6:
            continue
        for partner in atom.bond_partners:
            if partner.atomic_number == 6 and partner.idx in mapping:
                pairs.add(tuple(sorted((displayed, mapping[partner.idx]))))
    return sorted(pairs)
