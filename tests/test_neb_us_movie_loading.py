"""Tclsh regression tests for the generated NEB-US VMD loader."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from test_neb_us_movie_vmd import _render


_VMD_STUB = r"""
namespace eval ::mock {variable frame -1; variable frames 1890; variable labels {}; variable fail_dcd 0; variable bounds_seen 0; variable display_updates 0}
proc mol {command args} {
    if {$command eq "new"} {return 7}
    if {$command eq "addfile"} {
        if {$::mock::fail_dcd && [string match *path.dcd [lindex $args 0]]} {error "DCD load failed"}
        return
    }
    if {$command eq "representation" && [lindex $args 0] eq "VDW" && [llength $args] != 3} {error "invalid VDW arity"}
    return
}
proc molinfo {molid command args} {
    if {$molid eq "list"} {return 7}
    if {$command eq "get"} {
        if {[lindex $args 0] eq "numframes"} {return $::mock::frames}
        if {[lindex $args 0] eq "frame"} {return $::mock::frame}
        return 10
    }
    if {$command eq "set" && [lindex $args 0] eq "frame"} {set ::mock::frame [lindex $args 1]}
    return
}
proc color {args} {return}
proc ::label {command args} {
    if {$command eq "add" && $::mock::frame != 0} {error "label requested without active coordinates"}
    if {$command eq "add"} {lappend ::mock::labels $args; return}
    if {$command eq "list"} {return $::mock::labels}
    return
}
proc atomselect {molid selection args} {
    set token ::mock::sel[incr ::mock::next_sel]
    proc $token {command args} {
        if {$command eq "getbonds"} {return [lrepeat 10 {}]}
        if {$command eq "num"} {return 1}
        return
    }
    return $token
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
proc display {command} {incr ::mock::display_updates}
"""


def _run_tcl(
    script: str, *, fail_dcd: bool = False
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
            + "set code [catch {\n"
            + script
            + "\n} message options]\n"
            + "puts [list CODE $code MESSAGE $message FRAME $::mock::frame "
            + "LABELS [llength $::mock::labels] BOUNDS $::mock::bounds_seen "
            + "DISPLAY $::mock::display_updates]\n"
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
    assert "DISPLAY 1" in result.stdout
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


def test_loader_rejects_a_nonactivated_frame_before_labels_or_camera() -> None:
    """A regression that removes frame initialization fails at its frame guard."""
    script = _render().replace("molinfo $molid set frame 0\n", "")
    result = _run_tcl(script)
    assert result.returncode == 0, result.stderr
    assert "CODE 1" in result.stdout
    assert "Could not activate trajectory frame zero" in result.stdout
    assert "LABELS 0" in result.stdout
    assert "BOUNDS 0" in result.stdout
