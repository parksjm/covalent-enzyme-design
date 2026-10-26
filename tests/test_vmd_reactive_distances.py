"""Validate reactive label selection and VMD frame-stage Tcl behavior."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from nyl12_covalent.vmd_reactive_distances import append_reactive_labels

ROOT = Path(__file__).resolve().parents[1]


def _config() -> dict:
    return {
        "measurement_atoms": {
            "O": {"index": 1, "name": "OG1"},
            "H": {"index": 2, "name": "HG1"},
            "D": {"index": 3, "name": "OD1"},
            "C": {"index": 4, "name": "C54"},
            "N": {"index": 5, "name": "N"},
        },
        "stage_measurements": {
            "1": [["O", "H"], ["H", "D"], ["O", "C"]],
            "2": [["D", "H"], ["H", "N"], ["O", "C"]],
        },
    }


def test_measurement_loader_maps_frames_and_checks_atoms(tmp_path: Path) -> None:
    viewer = tmp_path / "view.vmd"
    viewer.write_text("# existing viewer\n")
    append_reactive_labels(viewer, [{"stage": 1}, {"stage": 2}], _config())
    text = viewer.read_text()
    assert text.startswith("# existing viewer")
    assert "{1 2}" in text
    assert "{1 OG1 2 HG1}" in text
    assert "{3 OD1 2 HG1}" in text
    assert "{2 HG1 5 N}" in text


def test_optional_dynamic_pairs_and_persistent_labels_are_exported(
    tmp_path: Path,
) -> None:
    config = _config()
    config["measurement_atoms"].update(
        {
            "O49": {"index": 10, "name": "O49"},
            "PheN": {"index": 11, "name": "N"},
            "AsnND2": {"index": 12, "name": "ND2"},
        }
    )
    config["stage_dynamic_bonds"] = {
        "1": [["O", "H"], ["O", "C"]],
        "2": [["D", "H"], ["H", "N"]],
    }
    config["dynamic_bond_cutoffs"] = {
        "O-H": 1.25,
        "O-C": 1.9,
        "D-H": 1.25,
        "H-N": 1.25,
    }
    config["persistent_measurements"] = {
        "O49_to_Phe113_N": ["O49", "PheN"],
        "O49_to_Asn186_ND2": ["O49", "AsnND2"],
    }
    config["proton_atom"] = "H"
    viewer = tmp_path / "view.vmd"
    viewer.write_text("# existing viewer\n")
    append_reactive_labels(viewer, [{"stage": 1}, {"stage": 2}], config)
    text = viewer.read_text()
    assert "{1 OG1 2 HG1 1.25}" in text
    assert "{1 OG1 4 C54 1.9}" in text
    assert "{O49_to_Phe113_N 10 O49 11 N}" in text
    assert "{O49_to_Asn186_ND2 10 O49 12 ND2}" in text
    assert text.rstrip().endswith("2 HG1")
    tcl = (ROOT / "scripts/nyl12_reactive_distances.tcl").read_text()
    assert "mol representation DynamicBonds $cutoff 0.12 12" in tcl
    assert "$first getbonds" in tcl and "$first setbonds" in tcl


@pytest.mark.parametrize("bad", ["stage", "empty", "same", "index", "name"])
def test_invalid_labels_fail_before_viewer_changes(tmp_path: Path, bad: str) -> None:
    config = _config()
    if bad == "stage":
        config["stage_measurements"].pop("1")
    elif bad == "empty":
        config["stage_measurements"]["1"] = []
    elif bad == "same":
        config["stage_measurements"]["1"] = [["O", "O"]]
    elif bad == "index":
        config["measurement_atoms"]["O"]["index"] = -1
    else:
        config["measurement_atoms"]["O"]["name"] = "OG1;exit"
    viewer = tmp_path / "view.vmd"
    viewer.write_text("unchanged")
    with pytest.raises(ValueError):
        append_reactive_labels(viewer, [{"stage": 1}], config)
    assert viewer.read_text() == "unchanged"


def test_native_labels_switch_stage_without_deleting_other_labels(
    tmp_path: Path,
) -> None:
    # A Tcl execution test mocks only VMD APIs; native label values are VMD's job.
    mock = """
set labels {{unrelated}}
set shown [dict create 0 1]
set bad_identity 0
proc label {operation category args} {
    global labels shown
    switch $operation {
        list {return $labels}
        add {lappend labels $args}
        hide {dict set shown [lindex $args 0] 0}
        show {dict set shown [lindex $args 0] 1}
        default {error "Unexpected label operation: $operation"}
    }
}
proc selection {command} {
    global bad_identity
    if {$command eq "num"} {return [expr {!$bad_identity}]}
    if {$command ne "delete"} {error "Unexpected selection operation"}
}
proc atomselect {molid selection_text} {return selection}
proc molinfo {molid operation property} {
    if {$property eq "numframes"} {return 3}
    return $::vmd_frame($molid)
}
set ::vmd_frame(7) 0
set ::nyl12_mol 7
set ::nyl12_dir .
"""
    loader = tmp_path / "loader.tcl"
    loader.write_text(mock)
    append_reactive_labels(
        loader, [{"stage": 1}, {"stage": 1}, {"stage": 2}], _config()
    )
    # Replace only the generated relative source command for this test fixture.
    text = loader.read_text().replace(
        "[file join $::nyl12_dir nyl12_reactive_distances.tcl]",
        "{" + str(ROOT / "scripts/nyl12_reactive_distances.tcl") + "}",
    )
    text += """
if {[llength $labels] != 5} {error "Duplicate reverse-order bond label"}
if {$shown ne {0 1 1 1 2 1 3 1 4 0}} {error "Wrong stage1 labels: $shown"}
set ::vmd_frame(7) 1
if {[llength $labels] != 5} {error "Same-stage animation added labels"}
set ::vmd_frame(7) 2
if {$shown ne {0 1 1 0 2 1 3 1 4 1}} {error "Wrong stage2 labels: $shown"}
set ::vmd_frame(7) 0
if {$shown ne {0 1 1 1 2 1 3 1 4 0}} {error "Backward-stage switching failed"}
if {![catch {::nyl12_reactive::install 7 {1 1 2} {}}]} {error "Duplicate install accepted"}
if {![catch {::nyl12_reactive::install 8 {1} {}}]} {error "Frame mismatch accepted"}
set ::bad_identity 1
if {![catch {::nyl12_reactive::install 8 {1 1 1} {1 {{1 OG1 2 HG1}}}}]} {
    error "Atom identity mismatch accepted"
}
set ::vmd_frame(7) 9
"""
    # Tcl top-level errors must propagate via an explicit wrapper exit status.
    loader.write_text(
        "if {[catch {\n" + text + "\n} error]} {puts stderr $error; exit 1}\n"
    )
    # Out-of-range frame writes intentionally exercise the callback error.
    result = subprocess.run(["tclsh", str(loader)], capture_output=True, text=True)
    assert result.returncode == 1
    assert "No reaction stage for frame 9" in result.stderr
