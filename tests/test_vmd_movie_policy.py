"""Execute the portable policy in Tcl with a deterministic VMD command fixture."""

from types import SimpleNamespace
import tkinter

import pytest

from nyl12_covalent.vmd_movie_policy import mapped_cc_bonds, render_movie_policy


@pytest.fixture
def vmd() -> tkinter.Tcl:
    """Provide physical atoms, asymmetric boundary edges and scoped VMD commands."""
    tcl = tkinter.Tcl()
    tcl.eval(r"""
        set elements {C C H O C C H N C C}
        set names {C1 C2 H1 O1 CB CA HA N CA CB}
        set residues {0 0 0 0 1 1 1 1 2 2}
        set types {CT CT HC O CT CT H N CT CT}
        set bonds {{1 2} {3} {0} {1} {} {} {5} {5} {} {}}
        set reps {{Licorice 0.12 12 12} {Licorice 0.3 12 12} {DynamicBonds 1.7 0.12 16} {NewCartoon 0.3 10 4.1 0}}
        set selections {{index 0 1 2 3} {water and within 4 of index 0} {index 1 3} {protein}}
        set colors {Element Element Element Chain}
        set hidden {}
        # Canonical solid names from VMD User Guide, colorinfo colors (33).
        set solid_colors {blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3}
        set commands {}
        set nsel 0
        proc atomselect {m expression} {
            if {$m != 9} {error "wrong molecule"}
            set selected {}
            if {$expression eq "all"} {
                set selected {0 1 2 3 4 5 6 7 8 9}
            } elseif {[string match "protein and name CA CB*" $expression]} {
                set selected {4 5}
            } elseif {[regexp {^protein and chain ([A-H])$} $expression -> chain]} {
                if {$chain eq "A"} {set selected {4 5 6 7}}
                if {$chain eq "E"} {set selected {8 9}}
            } elseif {[regexp {^protein and index (.*)$} $expression -> indices]} {
                set selected $indices
            } else {error "unexpected selection: $expression"}
            set command sel[incr ::nsel]
            interp alias {} $command {} selected $selected
            return $command
        }
        proc selected {indices operation args} {
            if {$operation eq "delete"} {return}
            if {$operation eq "num"} {return [llength $indices]}
            if {$operation eq "setbonds"} {set ::bonds [lindex $args 0]; return}
            if {$operation eq "getbonds"} {return $::bonds}
            if {$operation eq "set"} {
                set ::types [lindex $args 1]
                # Real VMD registers/colors Type by its leading character.
                foreach kind $::types {
                    set prefix [string toupper [string index [string trimleft $kind 0123456789] 0]]
                    set ::registered($prefix) 1
                    if {![info exists ::palette($prefix)]} {set ::palette($prefix) blue}
                }
                return
            }
            set field [lindex $args 0]
            if {$field eq "index"} {return $indices}
            set values {}
            foreach i $indices {
                switch $field {
                    element {lappend values [lindex $::elements $i]}
                    name {lappend values [lindex $::names $i]}
                    residue {lappend values [lindex $::residues $i]}
                    type {lappend values [lindex $::types $i]}
                    default {error "unsupported property $field"}
                }
            }
            return $values
        }
        proc color {args} {
            if {[lindex $args end] ni $::solid_colors} {
                error "Cannot add color item: invalid color name '[lindex $args end]'"
            }
            if {[lrange $args 0 2] eq {add item Type}} {
                if {[info exists ::registered([lindex $args 3])]} {
                    error "Duplicate Type color registration"
                }
                set ::registered([lindex $args 3]) 1
                return
            }
            lassign $args category type shade
            if {$category ne "Type" || ![info exists ::registered($type)]} {
                error "Unregistered color category item: $category $type"
            }
            set ::palette($type) $shade
        }
        proc colorinfo {category {name ""}} {
            if {$category eq "colors" && $name eq ""} {return $::solid_colors}
            if {$category ne "category" || $name ne "Type"} {
                error "Unsupported colorinfo query"
            }
            return [array names ::registered]
        }
        proc display {args} {set ::display $args}
        proc molinfo {m args} {
            if {$m eq "top"} {return 9}
            if {$m != 9} {error "wrong molecule"}
            set field [lindex $args 1]
            if {$field eq "numreps"} {return [llength $::reps]}
            lassign [lindex $field 0] property r
            if {$property eq "rep"} {return [list [lindex $::reps $r]]}
            if {$property eq "selection"} {return [list [lindex $::selections $r]]}
            error "unexpected molinfo $field"
        }
        proc mol {operation args} {
            lappend ::commands [linsert $args 0 $operation]
            switch $operation {
                representation {set ::style $args}
                selection {set ::selection [lindex $args 0]}
                color {set ::color $args}
                addrep {
                    if {[lindex $args 0] != 9} {error "wrong molecule"}
                    lappend ::reps $::style
                    lappend ::selections $::selection
                    lappend ::colors $::color
                }
                delrep {
                    lassign $args r m
                    if {$m != 9} {error "wrong molecule"}
                    set ::reps [lreplace $::reps $r $r]
                    set ::selections [lreplace $::selections $r $r]
                    set ::colors [lreplace $::colors $r $r]
                }
                modstyle {
                    lassign $args r m
                    set ::reps [lreplace $::reps $r $r [lrange $args 2 end]]
                }
                modcolor {
                    lassign $args r m
                    set ::colors [lreplace $::colors $r $r [lrange $args 2 end]]
                }
                showrep {lappend ::hidden $args}
                selupdate {set ::selupdate $args}
                default {error "unexpected mol operation $operation"}
            }
        }
    """)
    return tcl


def test_executed_tcl_policy(vmd: tkinter.Tcl) -> None:
    original_types = vmd.eval("set types")
    script = render_movie_policy([0, 1, 4, 7], [0, 1, 2, 3], [(0, 1)])
    vmd.eval("set molid 9")
    vmd.eval(script)
    assert vmd.eval("lindex $bonds 1") == "3 0"
    assert vmd.eval("lindex $bonds 4") == "5"
    assert vmd.eval("lindex $bonds 5") == "4"
    assert vmd.eval("lindex $bonds 8") == ""  # Different protein residue unchanged.
    assert vmd.eval("set palette(Q)") == "orange"
    assert vmd.eval("set palette(C)") == "cyan"
    assert vmd.eval("set palette(H)") == "white"
    assert vmd.eval("set palette(N)") == "blue"
    assert vmd.eval("set palette(O)") == "red"
    assert vmd.eval("lrange $types 0 4") == (
        "Q_NYL12_MOVIE_SUB_C Q_NYL12_MOVIE_SUB_C H_NYL12_MOVIE_H O_NYL12_MOVIE_O C_NYL12_MOVIE_C"
    )
    assert vmd.eval("lindex $reps 0") == "DynamicBonds 1.7 0.2 16"
    assert vmd.eval("lindex $colors 0") == "Type"
    assert vmd.eval("lindex $reps 1") == "Licorice 0.2 12 12"
    assert vmd.eval("lindex $colors 1") == "Type"
    selection = vmd.eval("lindex $selections 1")
    assert "index 0 1 2 3 4 5 7" in selection
    assert "water and within 4 of index 0" in selection
    assert vmd.eval("set display") == "depthcue off"
    assert vmd.eval("set hidden") == "{9 2 off} {9 3 off}"
    assert vmd.eval("lrange $colors 2 end") == "{ColorID 7} {ColorID 7}"
    vmd.eval(script)  # Reapplying cannot duplicate restored edges or replace originals.
    assert vmd.eval("lindex $bonds 1") == "3 0"
    assert vmd.eval("lindex $bonds 4") == "5"
    vmd.eval("::nyl12_movie_policy::restore_types 9")
    assert vmd.eval("set types") == original_types
    with pytest.raises(tkinter.TclError, match="No original types"):
        vmd.eval("::nyl12_movie_policy::restore_types 9")


def test_real_vmd_type_prefix_colors_do_not_collapse_to_blue(vmd: tkinter.Tcl) -> None:
    """VMD ignores the suffix of a Type name when looking up its color."""
    vmd.eval("set palette(N) blue")
    legacy = "NYL12_MOVIE_C NYL12_MOVIE_H NYL12_MOVIE_N NYL12_MOVIE_O"
    vmd.eval(f"set legacy {{{legacy}}}")
    assert (
        vmd.eval("lmap kind $legacy {set palette([string index $kind 0])}")
        == "blue blue blue blue"
    )
    vmd.eval("set molid 9")
    vmd.eval(render_movie_policy([0, 1, 4, 7], [0, 1, 2, 3], [(0, 1)]))
    assert vmd.eval("lmap kind $types {set palette([string index $kind 0])}") == (
        "orange orange white red cyan cyan white blue cyan cyan"
    )


def test_unknown_elements_fail_before_display_types_change(vmd: tkinter.Tcl) -> None:
    original = vmd.eval("set types")
    vmd.eval("lset elements 0 X")
    vmd.eval("set molid 9")
    with pytest.raises(tkinter.TclError, match="Unsupported display element"):
        vmd.eval(render_movie_policy([0, 1, 4, 7], [0, 1, 2, 3], [(0, 1)]))
    assert vmd.eval("set types") == original


def test_indexed_chain_groups(vmd: tkinter.Tcl) -> None:
    vmd.eval(
        render_movie_policy(
            [0, 4],
            [0, 1, 2, 3],
            molecule="9",
            chain_groups={"A": [4, 5, 6, 7], "E": [8, 9]},
        )
    )
    assert vmd.eval("lrange $selections 2 end") == (
        "{protein and index 4 5 6 7} {protein and index 8 9}"
    )
    with pytest.raises(tkinter.TclError, match="Invalid cartoon"):
        vmd.eval(render_movie_policy([0], [1], molecule="9", chain_groups={"A": [99]}))


@pytest.mark.parametrize(
    "qm,sub,pairs", [([99], [0], []), ([0], [0], [(0, 99)]), ([0], [0], [(0, 2)])]
)
def test_tcl_rejects_nonphysical_or_noncarbon_pairs(
    vmd: tkinter.Tcl, qm: list[int], sub: list[int], pairs: list[tuple[int, int]]
) -> None:
    with pytest.raises(tkinter.TclError):
        vmd.eval(render_movie_policy(qm, sub, pairs, molecule="[molinfo top]"))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"qm_indices": []},
        {"substrate_indices": []},
        {"qm_indices": [-1]},
        {"qm_indices": [True]},
        {"molecule": "$m; exit"},
        {"cc_bonds": [(0, 0)]},
        {"cc_bonds": [(0,)]},
        {"cc_bonds": [(0, -1)]},
        {"chain_groups": {"Z": [0]}},
        {"chain_groups": {"A": []}},
    ],
)
def test_python_validation(kwargs: dict) -> None:
    values = {"qm_indices": [0], "substrate_indices": [1]} | kwargs
    with pytest.raises(ValueError):
        render_movie_policy(**values)


def test_all_eight_chain_counterparts(vmd: tkinter.Tcl) -> None:
    groups = {chain: [i] for i, chain in enumerate("ABCDEFGH")}
    vmd.eval(render_movie_policy([0], [1], molecule="9", chain_groups=groups))
    assert vmd.eval("llength $hidden") == "8"
    assert vmd.eval("lrange $colors 2 end") == (
        "{ColorID 7} {ColorID 10} {ColorID 27} {ColorID 4} "
        "{ColorID 7} {ColorID 10} {ColorID 27} {ColorID 4}"
    )


def test_mapped_topology_edges() -> None:
    atoms = [
        SimpleNamespace(idx=i, atomic_number=z, bond_partners=[])
        for i, z in enumerate([6, 6, 1, 6])
    ]
    atoms[0].bond_partners = [atoms[1], atoms[2], atoms[3]]
    atoms[1].bond_partners = [atoms[0]]
    topology = SimpleNamespace(atoms=atoms)
    assert mapped_cc_bonds(topology, {0: 3, 1: 2, 2: 1}) == [(2, 3)]
    for mapping in ({0: 1, 1: 1}, {-1: 0}, {0: -1}, {True: 0}):
        with pytest.raises(ValueError):
            mapped_cc_bonds(topology, mapping)


def test_existing_proton_spheres_use_requested_scale_and_valid_vdw_arity(
    vmd: tkinter.Tcl,
) -> None:
    vmd.eval(
        "lappend reps {VDW 0.28 12}; lappend selections {index 2}; lappend colors Element"
    )
    vmd.eval(render_movie_policy([0, 1], [2, 3], [(0, 1)], molecule="9"))
    assert vmd.eval("lindex $reps 1") == "VDW 0.2 12"


def test_palette_uses_real_vmd_names_and_rejects_brown(vmd: tkinter.Tcl) -> None:
    """Reproduce the actual GUI error that a registration-only mock missed."""
    with pytest.raises(tkinter.TclError, match="invalid color name 'brown'"):
        vmd.eval("color add item Type NYL12_MOVIE_BR brown")
    script = render_movie_policy([0, 1], [2, 3], [(0, 1)], molecule="9")
    vmd.eval(script)
    assert vmd.eval("set palette(P)") == "tan"
    assert vmd.eval("set ::nyl12_movie_policy::completed(9)") == "20261126-typeprefix3"


def test_palette_failure_is_explicit_and_precedes_mutation(vmd: tkinter.Tcl) -> None:
    original = {key: vmd.eval(f"set {key}") for key in ("types", "bonds", "reps")}
    vmd.eval("set solid_colors [lsearch -all -inline -not -exact $solid_colors tan]")
    with pytest.raises(tkinter.TclError, match="Unsupported VMD color 'tan'"):
        vmd.eval(render_movie_policy([0], [1], molecule="9"))
    assert {key: vmd.eval(f"set {key}") for key in original} == original
    assert not int(vmd.eval("info exists ::nyl12_movie_policy::completed(9)"))
