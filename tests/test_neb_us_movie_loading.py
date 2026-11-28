"""Tclsh regression tests for the generated NEB-US VMD loader."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from test_neb_us_movie_vmd import _render


_VMD_STUB = r"""
namespace eval ::mock {variable frame -1; variable frames 1890; variable labels {}; variable fail_dcd 0; variable bounds_seen 0; variable display_updates 0; variable hidden 0; variable deleted 0; variable shown 0}
proc mol {command args} {
    if {$command eq "new"} {return 7}
    if {$command eq "off"} {set ::mock::hidden 1; return}
    if {$command eq "delete"} {incr ::mock::deleted; return}
    if {$command eq "on"} {incr ::mock::shown; return}
    if {$command eq "addfile"} {
        if {$::mock::fail_dcd && [string match *path.dcd [lindex $args 0]]} {error "DCD load failed"}
        return
    }
    if {$command eq "representation" && [lindex $args 0] eq "VDW" && [llength $args] != 3} {error "invalid VDW arity"}
    return
}
proc molinfo {molid {command ""} args} {
    if {$molid eq "list"} {return 7}
    if {$command eq "get"} {
        if {[lindex $args 0] eq "numframes"} {return $::mock::frames}
        if {[lindex $args 0] eq "frame"} {return $::mock::frame}
        return 10
    }
    if {$command eq "set" && [lindex $args 0] eq "frame"} {set ::mock::frame [lindex $args 1]}
    return
}
set ::solid_colors {blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3}
array set ::registered {}
proc colorinfo {kind {name ""}} {
    if {$kind eq "colors"} {return $::solid_colors}
    if {$kind eq "category" && $name eq "Type"} {return [array names ::registered]}
    error "unsupported colorinfo query"
}
proc color {args} {
    if {[lrange $args 0 2] eq {add item Type}} {set ::registered([lindex $args 3]) 1}
    return
}
proc ::label {command args} {
    if {$command eq "add" && $::mock::frame != 0} {error "label requested without active coordinates"}
    if {$command eq "add"} {lappend ::mock::labels $args; return}
    if {$command eq "list"} {return $::mock::labels}
    return
}
proc atomselect {molid selection args} {
    if {$selection eq "all"} {set indices {0 1 2 3 4 5 6 7 8 9}}
    if {[string match "protein and name CA CB*" $selection]} {set indices {}}
    if {[regexp {index ([0-9 ]+)} $selection -> chosen]} {set indices $chosen}
    if {[regexp {^protein and chain ([A-H])$} $selection]} {set indices {0 1 2 3 4 5 6 7 8 9}}
    set token ::mock::sel[incr ::mock::next_sel]
    interp alias {} $token {} ::mock::selected $indices
    return $token
}
proc ::mock::selected {indices command args} {
    if {$command eq "get"} {
        set field [lindex $args 0]
        if {$field eq "index"} {return $indices}
        if {$field eq "element"} {return [lrepeat [llength $indices] C]}
        if {$field eq "name"} {return [lrepeat [llength $indices] CA]}
        if {$field eq "residue"} {return [lrepeat [llength $indices] 0]}
        if {$field eq "type"} {return [lrepeat [llength $indices] CT]}
        if {$field eq "bonds"} {return [lrepeat 10 {}]}
    }
    if {$command eq "num"} {return [llength $indices]}
    if {$command eq "setbonds" || $command eq "set" || $command eq "delete"} {return}
    if {$command eq "getbonds"} {return [lrepeat 10 {}]}
    return
}
namespace eval ::mock {variable next_sel 0}
proc measure {command selection} {
    if {$command ne "minmax"} {error "unexpected measure command"}
    if {$::mock::frame != 0} {error "bounds requested at frame $::mock::frame"}
    set ::mock::bounds_seen 1
    return {{0 0 0} {4 5 6}}
}
proc vecadd {first second} {lmap a $first b $second {expr {$a + $b}}}
proc vecsub {first second} {lmap a $first b $second {expr {$a - $b}}}
proc vecscale {scale vector} {lmap value $vector {expr {$scale * $value}}}
proc display {args} {incr ::mock::display_updates}
"""


def _run_tcl(
    script: str, *, fail_dcd: bool = False, frames: int = 1890
) -> subprocess.CompletedProcess[str]:
    """Run a generated VMD script against a small Tcl VMD stub."""
    tclsh = shutil.which("tclsh")
    if tclsh is None:
        pytest.skip("tclsh is unavailable")
    with tempfile.TemporaryDirectory() as directory:
        workdir = Path(directory)
        (workdir / "path.parm7").touch()
        (workdir / "path.dcd").touch()
        driver = (
            _VMD_STUB
            + f"set ::mock::fail_dcd {int(fail_dcd)}\n"
            + f"set ::mock::frames {frames}\n"
            + "set code [catch {\n"
            + script
            + "\n} message options]\n"
            + "puts [list CODE $code MESSAGE $message FRAME $::mock::frame "
            + "LABELS [llength $::mock::labels] BOUNDS $::mock::bounds_seen "
            + "DISPLAY $::mock::display_updates HIDDEN $::mock::hidden "
            + "DELETED $::mock::deleted SHOWN $::mock::shown]\n"
        )
        driver_path = workdir / "driver.tcl"
        driver_path.write_text(driver)
        return subprocess.run(
            [tclsh, str(driver_path)], text=True, capture_output=True, check=False
        )


def test_generated_loader_initializes_frame_before_labels_and_camera() -> None:
    """The loader uses frame zero and valid VDW arguments before measuring."""
    result = _run_tcl(_render())
    assert result.returncode == 0, result.stderr
    assert "CODE 0" in result.stdout
    assert "FRAME 0" in result.stdout
    assert "BOUNDS 1" in result.stdout
    assert "DISPLAY 2" in result.stdout
    assert "invalid VDW arity" not in result.stdout
    assert "VDW 0.28 12 12" not in _render()


def test_missing_dcd_fails_before_labels_or_camera() -> None:
    """A DCD load error is reported before labels and camera setup proceed."""
    result = _run_tcl(_render(), fail_dcd=True)
    assert result.returncode == 0, result.stderr
    assert "CODE 1" in result.stdout
    assert "DCD load failed" in result.stdout
    assert "LABELS 0" in result.stdout
    assert "BOUNDS 0" in result.stdout
    assert "DISPLAY 0" in result.stdout
    assert "HIDDEN 1 DELETED 1 SHOWN 0" in result.stdout


def test_stale_download_cannot_leave_default_whole_protein_scene() -> None:
    """A size-skipped old DCD fails visibly and removes its incomplete molecule."""
    result = _run_tcl(_render(), frames=945)
    assert result.returncode == 0, result.stderr
    assert "CODE 1" in result.stdout
    assert "expected 1890, found 945" in result.stdout
    assert "Rerun the updated pull script" in result.stdout
    assert "NYL12 MOVIE LOAD FAILED" in result.stderr
    assert "HIDDEN 1 DELETED 1 SHOWN 0" in result.stdout
    assert "LABELS 0" in result.stdout


def test_loader_rejects_a_nonactivated_frame_before_labels_or_camera() -> None:
    """A regression that removes frame initialization fails at its frame guard."""
    script = _render().replace("molinfo $molid set frame 0\n", "")
    result = _run_tcl(script)
    assert result.returncode == 0, result.stderr
    assert "CODE 1" in result.stdout
    assert "Could not activate trajectory frame zero" in result.stdout
    assert "LABELS 0" in result.stdout
    assert "BOUNDS 0" in result.stdout
