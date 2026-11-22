"""Tests for mapped NEB umbrella movie VMD scripts."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any

import pytest

from nyl12_covalent.neb_us_movie_vmd import render_movie_vmd


@dataclass
class _Atom:
    idx: int
    atomic_number: int
    residue: Any
    bond_partners: list[Any] = field(default_factory=list)


def _topology() -> Any:
    """Build a small topology with an omitted atom before kept atoms."""
    atoms: list[_Atom] = []
    for index in range(11):
        residue_name = "P6N" if 1 <= index <= 6 or index == 10 else "ALA"
        atoms.append(_Atom(index, 6, SimpleNamespace(name=residue_name)))
    atoms[0].atomic_number = 6  # Omitted unrelated atom shifts every local index.
    atoms[1].atomic_number = 6  # Scissile carbon, original ID 2.
    atoms[2].atomic_number = 8  # Thr OG1, ID 3.
    atoms[3].atomic_number = 7  # Leaving N, ID 4.
    atoms[4].atomic_number = 8  # Carbonyl O, ID 5.
    atoms[5].atomic_number = 1  # QM H bonded to QM carbonyl O, ID 6.
    atoms[6].atomic_number = 1  # HG1, ID 7.
    atoms[7].atomic_number = 7  # Thr backbone N, ID 8.
    atoms[8].atomic_number = 7  # Phe backbone N, ID 9.
    atoms[9].atomic_number = 7  # Asn ND2, ID 10.
    for first, second in [(1, 2), (1, 3), (1, 4), (4, 5), (2, 6), (1, 10)]:
        atoms[first].bond_partners.append(atoms[second])
        atoms[second].bond_partners.append(atoms[first])
    return SimpleNamespace(atoms=atoms)


def _render() -> str:
    return render_movie_vmd(
        _topology(),
        list(range(2, 12)),
        [2, 3, 4, 5, 6, 7],
        {
            "carbon": 2,
            "og1": 3,
            "leaving_n": 4,
            "oxygen": 5,
            "hg1": 7,
            "thr_n": 8,
            "phe_n": 9,
            "asn_nd2": 10,
        },
        1890,
    )


def test_render_maps_original_ids_and_builds_reactive_and_qm_bonds() -> None:
    script = _render()
    assert "mol selection {index 0 1 2 3 4 5 9}" in script
    assert "Unexpected trajectory frame count" in script
    assert "get numframes] != 1890" in script
    assert "foreach pair {{1 0} {0 2}}" in script
    assert "{3 4}" in script  # Static QM H–O bond is removed.
    assert "DynamicBonds 1.85" in script
    assert "DynamicBonds 1.4" in script
    assert "VDW 0.28" in script


def test_render_sets_requested_colors_labels_angle_and_local_camera() -> None:
    script = _render()
    assert "color Element H white" in script
    assert "color Element N blue" in script
    assert "color Element O red" in script
    assert "color change rgb 3 1.0 0.5 0.0" in script
    assert "ColorID 3" in script
    assert script.count("::neb_us_movie::add_label Bonds") == 7
    assert "::neb_us_movie::add_label Angles [list 1 0 3]" in script
    assert "Licorice 0.25 12 12; mol color ColorID 3; mol selection {index 9}" in script
    assert "set visible [atomselect $molid {index 0 1 2 3 4 5 9}]" in script
    assert "center_matrix" in script and "scale_matrix" in script
    assert "Concatenated biased windows" in script


def test_script_cleans_only_labels_owned_by_its_molecule() -> None:
    script = _render()
    assert "if {$owned} {label delete $kind $i}" in script
    assert "foreach spec [lindex $labels $i]" in script
    assert "label delete $kind" in script
    assert "label delete $kind all" not in script
    assert "mol delete $molid" in script


def test_label_procedure_calls_global_vmd_label_command() -> None:
    tclsh = shutil.which("tclsh")
    if tclsh is None:
        pytest.skip("tclsh is unavailable")
    script = _render()
    start = script.index("proc ::neb_us_movie::add_label")
    end = script.index("\n}", start) + 2
    stub = """namespace eval ::neb_us_movie {variable molecule 23}
set ::calls {}
proc ::label {command args} {
    lappend ::calls [linsert $args 0 $command]
    if {$command eq "list"} {return {existing}}
}
"""
    run = subprocess.run(
        [tclsh],
        input=stub + script[start:end] + "\n::neb_us_movie::add_label Bonds {1 2} {0.04 0.08}\nputs $::calls\n",
        text=True,
        capture_output=True,
        check=True,
    )
    assert "add Bonds 23/1 23/2" in run.stdout
    assert "show Bonds 0" in run.stdout
    assert "textoffset Bonds 0" in run.stdout


@pytest.mark.parametrize(
    ("kept", "qm", "roles", "frames", "message"),
    [
        ([2, 3], [2], {"carbon": 2}, 4, "Missing required role"),
        ([3, 2], [2], {}, 4, "sorted"),
        ([2, 3], [2, 9], {}, 4, "QM atom IDs"),
        ([2, 3], [2], {"carbon": 99}, 4, "Missing required role"),
        ([2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 3], {
            "carbon": 2, "og1": 3, "leaving_n": 4, "hg1": 7,
            "thr_n": 8, "oxygen": 5, "phe_n": 9, "asn_nd2": 10,
        }, 0, "positive"),
    ],
)
def test_invalid_inputs_fail_explicitly(
    kept: list[int],
    qm: list[int],
    roles: dict[str, int],
    frames: int,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        render_movie_vmd(_topology(), kept, qm, roles, frames)


def test_role_must_be_retained_and_in_topology_range() -> None:
    roles = {
        "carbon": 2, "og1": 3, "leaving_n": 4, "hg1": 7,
        "thr_n": 8, "oxygen": 5, "phe_n": 9, "asn_nd2": 10,
    }
    with pytest.raises(ValueError, match="included in kept"):
        render_movie_vmd(_topology(), [2, 3, 4, 5, 6, 7, 8, 9], [2], roles, 1)
    roles["asn_nd2"] = 12
    with pytest.raises(ValueError, match="invalid original atom ID"):
        render_movie_vmd(_topology(), list(range(2, 12)), [2], roles, 1)
