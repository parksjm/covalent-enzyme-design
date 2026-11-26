"""Topology/rendered-coordinate and stateful headless Tcl regressions."""

from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.relay_progress_vmd import substrate_indices, vmd_script
from nyl12_covalent import relay_progress_vmd


def scene() -> tuple[dict, dict, SimpleNamespace, np.ndarray]:
    """Create a reordered strip map, remote PRO and unrelated polymer component."""
    names = [
        "PRO",
        "PRO",
        "THR",
        "THR",
        "P6N",
        "P6N",
        "P6I",
        "P6I",
        "P6C",
        "P6C",
        "WAT",
        "WAT",
        "WAT",
        "PHE",
        "ASN",
        "THR",
        "P6N",
    ]
    elements = [6, 1, 8, 7, 6, 1, 6, 8, 7, 1, 8, 1, 1, 7, 7, 1, 6]
    atoms = [
        SimpleNamespace(
            idx=i,
            atomic_number=z,
            residue=SimpleNamespace(name=name, chain="A", idx=i),
            bond_partners=[],
        )
        for i, (name, z) in enumerate(zip(names, elements, strict=True))
    ]
    for first, second in [
        (0, 1),
        (4, 5),
        (4, 6),
        (6, 7),
        (6, 8),
        (8, 9),
        (2, 6),
        (10, 11),
        (10, 12),
        (11, 12),
    ]:
        atoms[first].bond_partners.append(atoms[second])
        atoms[second].bond_partners.append(atoms[first])
    mapping = [15, 1, 10, 7, 8, 6, 0, 13, 3, 2, 11, 12, 5, 4, 14, 9, 16]
    case = {
        "qm_atom_ids": [3, 4, 7, 8, 9, 11, 12, 13, 14, 15, 16],
        "atom_map": {
            "thr_og": 3,
            "thr_n": 4,
            "carbon": 7,
            "oxygen": 8,
            "leaving_n": 9,
            "water_o": 11,
            "water_h1": 12,
            "water_h2": 13,
            "proton": 16,
            "phe_n": 14,
            "asn_nd2": 15,
        },
    }
    xyz = np.arange(51, dtype=float).reshape(17, 3) / 10
    xyz[[0, 1, 16]] = 10000.0
    return (
        {"stripped_to_original_atom_indices_0based": mapping},
        case,
        SimpleNamespace(atoms=atoms),
        xyz,
    )


def test_substrate_component_and_invalid_mapping() -> None:
    _, _, topology, _ = scene()
    assert substrate_indices(topology, 6, 8) == {4, 5, 6, 7, 8, 9}
    with pytest.raises(ValueError, match="P6N/P6I/P6C"):
        substrate_indices(topology, 0, 8)
    with pytest.raises(ValueError, match="different substrate components"):
        substrate_indices(topology, 16, 8)


@pytest.mark.parametrize("adopt_old", [False, True])
def test_headless_view_scoping_labels_and_coordinates(
    tmp_path: Path, adopt_old: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    tclsh = shutil.which("tclsh")
    if tclsh is None:
        pytest.skip("Headless Tcl interpreter unavailable")
    # The shared Tcl policy has its own stateful tests; this harness isolates
    # ownership, labels and camera commands from the policy's richer VMD API.
    monkeypatch.setattr(
        relay_progress_vmd, "render_movie_policy", lambda *args, **kwargs: ""
    )
    manifest, case, topology, xyz = scene()
    text = vmd_script(manifest, case, topology)
    assert "residue 1311" not in text and "label show Bonds all" not in text
    assert text.count("mol representation DynamicBonds 1.4 0.1 12") == 3
    assert "DynamicBonds 1.85 0.1 12" in text
    folder = tmp_path / "bundle with spaces/progress/update_000001/vmd"
    folder.mkdir(parents=True)
    script = folder / "view_path.vmd"
    script.write_text(text)
    keep = manifest["stripped_to_original_atom_indices_0based"]
    reverse = {old: new for new, old in enumerate(keep)}
    visible = sorted(
        reverse[i] for i in (set(np.array(case["qm_atom_ids"]) - 1) | set(range(4, 10)))
    )
    coordinates = " ".join("{" + " ".join(map(str, xyz[i])) + "}" for i in keep)
    bonds = " ".join(
        "{"
        + " ".join(str(reverse[a.idx]) for a in topology.atoms[i].bond_partners)
        + "}"
        for i in keep
    )
    center = (
        xyz[[keep[i] for i in visible]].min(axis=0)
        + xyz[[keep[i] for i in visible]].max(axis=0)
    ) / 2
    harness = tmp_path / "stub.tcl"
    harness.write_text(
        _STUB
        + f"""\nset xyz {{{coordinates}}}
set initial_bonds {{{bonds}}}
dict set files 77 [list [list /foreign/path.parm7] [list /foreign/path.dcd]]
if {{{int(adopt_old)}}} {{
 dict set names 88 legacy_current_bundle
 dict set files 88 [list [list [file join [file dirname [lindex $argv 0]] path.parm7]] [list [file join [file dirname [lindex $argv 0]] path.dcd]]]
 dict lappend labels Bonds {{{{88 0}} {{88 1}} 100 show}}
}}
source [lindex $argv 0]
set first $molid
source [lindex $argv 0]
check {{[dict size $names] == 2}} {{Repeated source duplicated campaign molecule}}
check {{![dict exists $names 88]}} {{Unregistered exact-bundle molecule was not adopted}}
check {{![dict exists $names $first]}} {{Old owned molecule was not removed}}
check {{[llength [dict get $labels Bonds]] == 7}} {{Expected six local distances plus untouched foreign label}}
check {{[lindex [dict get $labels Bonds] 0] eq {{{{77 0}} {{77 1}} 100 hide}}}} {{Foreign label changed}}
check {{[llength [dict get $labels Angles]] == 1}} {{BD angle missing/duplicated}}
check {{[dict get $transforms 77] eq {{user_transform}}}} {{Foreign view changed}}
foreach i $last_visible {{
 foreach j [lindex $current_bonds $i] {{
  if {{[lsearch -exact $last_visible $j] >= 0}} {{
   set covered 0
   foreach selection [dict get $licorice_reps $molid] {{
    if {{[lsearch -exact $selection $i] >= 0 && [lsearch -exact $selection $j] >= 0}} {{set covered 1}}
   }}
   check {{$covered}} {{Retained bond endpoints split between Licorice selections}}
  }}
 }}
}}
check {{$last_visible eq {{{" ".join(map(str, visible))}}}}} {{Visible camera selected PRO or unrelated polymer}}
set center {{}}
foreach row [lrange [dict get $transforms $molid center_matrix] 0 2] {{lappend center [expr {{-[lindex $row 3]}}]}}
foreach actual $center expected {{{" ".join(map(str, center))}}} {{
 check {{abs($actual-$expected) < 1.e-10}} {{View centered on hidden/remote coordinates}}
}}
set activeH {reverse[11]}
check {{[lindex [lindex [dict get $labels Bonds] 3] 1] eq [list $molid $activeH]}} {{Wrong transferring H in OG1 label}}
check {{[lsearch -exact [lindex $current_bonds {reverse[6]}] {reverse[8]}] < 0}} {{Static scissile bond retained}}
check {{[lsearch -exact [lindex $current_bonds {reverse[2]}] {reverse[6]}] < 0}} {{Static attack bond duplicated}}
check {{[lsearch -exact [lindex $current_bonds {reverse[10]}] {reverse[11]}] < 0}} {{Static labile OH retained}}
check {{[lsearch -exact [lindex $current_bonds {reverse[11]}] {reverse[12]}] < 0}} {{Water HH constraint rendered as bond}}
check {{![dict exists $offsets Bonds 0]}} {{Foreign label offset changed}}
check {{[lsearch -exact [lindex $current_bonds {reverse[4]}] {reverse[5]}] >= 0}} {{CH bond lost}}
::relay_progress::show_network $molid
set count [llength [dict get $labels Bonds]]
::relay_progress::show_network $molid
check {{[llength [dict get $labels Bonds]] == $count}} {{Optional labels duplicated}}
puts {{Headless scoped VMD regression passed}}
"""
    )
    result = subprocess.run(
        [tclsh, str(harness), str(script)], check=True, capture_output=True, text=True
    )
    assert "regression passed" in result.stdout


def test_headless_ambiguous_adoption_fails_without_deleting(tmp_path: Path) -> None:
    tclsh = shutil.which("tclsh")
    if tclsh is None:
        pytest.skip("Headless Tcl interpreter unavailable")
    manifest, case, topology, _ = scene()
    folder = tmp_path / "progress/update_000001/vmd"
    folder.mkdir(parents=True)
    script = folder / "view_path.vmd"
    script.write_text(vmd_script(manifest, case, topology))
    harness = tmp_path / "ambiguous.tcl"
    harness.write_text(
        _STUB
        + r"""
foreach i {88 89} {
 dict set names $i duplicate_current
 dict set files $i [list [list [file join [file dirname [lindex $argv 0]] path.parm7]] [list [file join [file dirname [lindex $argv 0]] path.dcd]]]
}
set before $names
set status [catch {source [lindex $argv 0]} message]
check {$status == 1} {Ambiguous adoption did not fail}
check {[string match {*Multiple exact topology+trajectory*} $message]} {Unclear ambiguity error}
check {$names eq $before} {Molecules changed before ambiguous-adoption failure}
puts {Ambiguous adoption rejected safely}
"""
    )
    result = subprocess.run(
        [tclsh, str(harness), str(script)], check=True, capture_output=True, text=True
    )
    assert "rejected safely" in result.stdout


_STUB = r"""
set names [dict create 77 other_molecule]
set labels [dict create Atoms {} Bonds {{{77 0} {77 1} 100 hide}} Angles {} Dihedrals {}]
set transforms [dict create 77 user_transform]
set files [dict create 77 {{unrelated.parm7}}]
set offsets {}
set licorice_reps {}
set current_rep {}
set current_selection {}
set counter 0
proc check {condition message} {if {![uplevel 1 [list expr $condition]]} {error $message}}
proc mol {command args} {
 global counter names files licorice_reps current_rep current_selection
 switch $command {
  new {incr counter; dict set names $counter new; dict set files $counter [list [list [lindex $args 0]]]; return $counter}
  addfile {dict lappend files [lindex $args end] [list [lindex $args 0]]}
  representation {set current_rep [lindex $args 0]}
  selection {set current_selection [lindex $args 0]}
  addrep {
   if {$current_rep eq "Licorice"} {dict lappend licorice_reps [lindex $args 0] [lrange $current_selection 1 end]}
  }
  rename {dict set names [lindex $args 0] [lindex $args 1]}
  delete {dict unset names [lindex $args 0]; dict unset files [lindex $args 0]}
 }
}
proc molinfo {molid args} {
 global names transforms files
 if {$molid eq "list"} {return [dict keys $names]}
 if {[lindex $args 0] eq "get"} {
  if {[lindex $args 1] eq "filename"} {return [list [dict get $files $molid]]}
  return [list [dict get $names $molid]]
 }
 set attrs [lindex $args 1]; set values [lindex $args 2]
 check {[llength $attrs] == [llength $values]} {molinfo attribute/value count mismatch}
 foreach attr $attrs value $values {
  check {[llength $value] == 4} {molinfo matrix must have four rows}
  foreach row $value {check {[llength $row] == 4} {molinfo matrix row must have four values}}
  dict set transforms $molid $attr $value
 }
}
proc color {args} {}
proc label {command kind args} {
 global labels offsets
 set items [dict get $labels $kind]
 switch $command {
  list {return $items}
  add {
   set pairs {}
   foreach spec $args {
    if {![regexp {^([0-9]+)/([0-9]+)$} $spec -> molid atomid]} {error {label add requires molID/atomID}}
    lappend pairs [list $molid $atomid]
   }
   lappend items [concat $pairs {1.0 show}]; dict set labels $kind $items; return {}
  }
  delete {set index [lindex $args 0]; dict set labels $kind [lreplace $items $index $index]}
  show {
   set index [lindex $args 0]
   check {[string is integer -strict $index] && $index >= 0 && $index < [llength $items]} {label show requires a valid integer index}
  }
  textoffset {
   set index [lindex $args 0]; set offset [lindex $args 1]
   check {[string is integer -strict $index] && $index >= 0 && $index < [llength $items]} {textoffset requires valid label index}
   check {[llength $offset] == 2} {textoffset requires two values}
   dict set offsets $kind $index $offset
  }
 }
}
proc atomselect {molid selection} {
 global selection_text
 set selection_text $selection
 return atom_selection
}
proc atom_selection {command args} {
 global initial_bonds current_bonds
 switch $command {
  getbonds {return $initial_bonds}
  setbonds {set current_bonds [lindex $args 0]}
 }
}
proc measure {command selection} {
 global selection_text xyz last_visible
 set last_visible [lrange $selection_text 1 end]
 set lower {Inf Inf Inf}; set upper {-Inf -Inf -Inf}
 foreach i $last_visible {
  set point [lindex $xyz $i]
  set lo {}; set hi {}
  foreach x $point a $lower b $upper {lappend lo [expr {min($a,$x)}]; lappend hi [expr {max($b,$x)}]}
  set lower $lo; set upper $hi
 }
 return [list $lower $upper]
}
proc vecadd {a b} {set out {}; foreach x $a y $b {lappend out [expr {$x+$y}]}; return $out}
proc vecsub {a b} {set out {}; foreach x $a y $b {lappend out [expr {$x-$y}]}; return $out}
proc vecscale {scale a} {set out {}; foreach x $a {lappend out [expr {$scale*$x}]}; return $out}
proc trans {command args} {return [list $command {*}$args]}
"""


def test_policy_receives_mapped_qm_substrate_and_original_cc_edges() -> None:
    manifest, case, topology, _ = scene()
    text = vmd_script(manifest, case, topology)
    assert (
        "::nyl12_movie_policy::apply $molid {0 2 3 4 5 7 8 9 10 11 14} {3 4 5 12 13 15} {{5 13}}"
        in text
    )
