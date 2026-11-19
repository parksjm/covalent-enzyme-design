"""Route-specific legacy07 viewer regression against strict VMD Tcl contracts."""

import importlib.util
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from test_relay_progress_vmd import _STUB, scene

ROOT = Path(__file__).resolve().parents[1]
SUPPORT = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/water_path_vmd_support.py"
)


def context() -> tuple:
    """Create a QM102/54H topology retaining a separate PRO and polymer component."""
    manifest, case, topology, _ = scene()
    original_count = len(topology.atoms)
    for i in range(91):
        topology.atoms.append(
            SimpleNamespace(
                idx=original_count + i,
                atomic_number=1 if i < 51 else (8 if i == 90 else 6),
                residue=SimpleNamespace(name="THR"),
                bond_partners=[],
            )
        )
    keep = manifest["stripped_to_original_atom_indices_0based"]
    keep.extend(range(original_count, len(topology.atoms)))
    case["qm_atom_ids"].extend(range(original_count + 1, len(topology.atoms) + 1))
    case["atom_map"]["asp_od1"] = len(topology.atoms)
    mapping = {old + 1: new for new, old in enumerate(keep)}
    keys = {
        "og1": "thr_og",
        "thr_n": "thr_n",
        "hg1": "proton",
        "water_o": "water_o",
        "water_h1": "water_h1",
        "water_h2": "water_h2",
        "asp_od1": "asp_od1",
        "carbonyl_c": "carbon",
        "carbonyl_o": "oxygen",
        "leaving_n": "leaving_n",
        "phe113_n": "phe_n",
        "asn_nd2": "asn_nd2",
    }
    atoms = {key: mapping[case["atom_map"][value]] for key, value in keys.items()}
    qm = [mapping[i] for i in case["qm_atom_ids"]]
    manifest.update(frames=5, qm_indices_stripped_0based=qm)
    spec = importlib.util.spec_from_file_location("legacy07_test_support", SUPPORT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, manifest, case, topology, qm, atoms


def test_legacy_adapter_seven_route_labels_and_strict_tcl(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    support, manifest, case, topology, qm, atoms = context()
    text = support.vmd_script(5, qm, atoms, topology=topology, manifest=manifest)
    assert "mol representation DynamicBonds 1.4" in text
    assert text.count("mol representation DynamicBonds 1.4") == 54
    assert "ColorID 7" not in text and "ColorID 4" not in text
    assert "residue 1312" not in text and "CPK" not in text
    named = {
        (929, "OG1"): "og1",
        (929, "N"): "thr_n",
        (929, "HG1"): "hg1",
        (1576, "O"): "water_o",
        (1576, "H1"): "water_h1",
        (1576, "H2"): "water_h2",
        (968, "OD1"): "asp_od1",
        (1314, "C"): "carbonyl_c",
        (1314, "O"): "carbonyl_o",
        (1315, "N"): "leaving_n",
        (113, "N"): "phe113_n",
        (186, "ND2"): "asn_nd2",
    }
    keep = manifest["stripped_to_original_atom_indices_0based"]
    monkeypatch.setattr(
        support,
        "atom_index",
        lambda top, resid, name: keep[atoms[named[(resid, name)]]],
    )
    assert support.render_manifest(topology, manifest) == text
    tclsh = shutil.which("tclsh")
    if tclsh is None:
        pytest.skip("Headless Tcl unavailable")
    script = tmp_path / "legacy07/vmd_progress/view_path.vmd"
    script.parent.mkdir(parents=True)
    script.write_text(text)
    xyz = " ".join(f"{{{i / 100} 0 0}}" for i in range(len(keep)))
    harness = tmp_path / "legacy07_stub.tcl"
    harness.write_text(
        _STUB
        + f"""
set xyz {{{xyz}}}
set initial_bonds [lrepeat {len(keep)} {{}}]
source [lindex $argv 0]
check {{[llength [dict get $labels Bonds]] == 8}} {{Expected seven07labels plus foreign label}}
check {{[llength [dict get $labels Angles]] == 1}} {{Expected BDangle}}
check {{[lindex [lindex [dict get $labels Bonds] 4] 0] eq [list $molid {atoms["water_h2"]}]}} {{Wrong07transferring hydrogen}}
check {{[lindex [lindex [dict get $labels Bonds] 4] 1] eq [list $molid {atoms["asp_od1"]}]}} {{Wrong07acceptor}}
check {{[dict size [dict get $offsets Bonds]] == 7}} {{Missing per-label offsets}}
check {{![dict exists $offsets Bonds 0]}} {{Foreign label altered}}
set status [catch {{label add Bonds {{77 0}} {{77 1}}}} error]
check {{$status == 1}} {{Stub incorrectly accepts pair lists forlabeladd}}
set status [catch {{molinfo $molid set center_matrix {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}}} error]
check {{$status == 1}} {{Stub incorrectly accepts barematrix}}
puts {{Legacy07 strict API and route regression passed}}
"""
    )
    result = subprocess.run(
        [tclsh, str(harness), str(script)], check=True, capture_output=True, text=True
    )
    assert "regression passed" in result.stdout


def test_legacy_adapter_missing_context_or_inventory_fails() -> None:
    support, manifest, _, topology, qm, atoms = context()
    with pytest.raises(TypeError):
        support.vmd_script(5, qm, atoms)
    with pytest.raises(ValueError, match="accepted frame"):
        support.vmd_script(0, qm, atoms, topology=topology, manifest=manifest)
    with pytest.raises(ValueError, match="QM102"):
        support.vmd_script(5, qm[:-1], atoms, topology=topology, manifest=manifest)
