"""Check input-mask and residue-mask handling for the generic stage exporter."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.stage2_vmd_export import (  # noqa: E402
    BURGI_DUNITZ_ANGLE,
    DYNAMIC_BONDS,
    REACTIVE_LABELS,
    _viewer_text,
    export_stage2,
    parse_atom_mask,
    read_frames,
    residue_mask,
)
from nyl12_covalent import stage2_vmd_export as export_module  # noqa: E402


@pytest.mark.parametrize(
    ("mask", "expected"),
    [
        ("@1,3-5,8", [0, 2, 3, 4, 7]),
        ("@2-2", [1]),
    ],
)
def test_parse_atom_mask_converts_amber_one_based_indices(
    mask: str, expected: list[int]
) -> None:
    assert parse_atom_mask(mask, 10) == expected


@pytest.mark.parametrize("mask", ["1-2", "@0", "@2-1", "@1-11", "@1,1", "@x"])
def test_parse_atom_mask_rejects_ambiguous_or_invalid_selection(mask: str) -> None:
    with pytest.raises(ValueError):
        parse_atom_mask(mask, 10)


@dataclass
class _Residue:
    idx: int
    name: str = "LIG"
    chain: str = "A"


@dataclass
class _Atom:
    residue: _Residue
    name: str = "X"
    idx: int = 0
    atomic_number: int = 6
    bond_partners: list[Any] = field(default_factory=list)


class _Topology:
    def __init__(self, residues: list[int]) -> None:
        self.atoms = [_Atom(_Residue(index)) for index in residues]
        self.atoms[0].residue.name = "ALA"
        for index, atom in enumerate(self.atoms):
            atom.idx = index


def test_residue_mask_compacts_selected_atom_residues() -> None:
    topology = _Topology([0, 0, 1, 2, 4, 4, 6])
    assert residue_mask(topology, np.array([0, 1, 3, 4, 6])) == ":1,3,5,7"


def test_residue_mask_rejects_empty_display_selection() -> None:
    with pytest.raises(ValueError, match="empty"):
        residue_mask(_Topology([0]), np.array([], dtype=int))


def test_viewer_uses_verified_zero_based_reactive_atoms_and_displays() -> None:
    """Keep the trajectory viewer aligned with established catalytic displays."""
    identities = {
        index: name
        for row in REACTIVE_LABELS
        for index, name in ((row[1], row[2]), (row[3], row[4]))
    }
    for index, name in zip(
        BURGI_DUNITZ_ANGLE[::2], BURGI_DUNITZ_ANGLE[1::2], strict=True
    ):
        identities[index] = name
    topology = _Topology([0] * (max(identities) + 1))
    for index, name in identities.items():
        topology.atoms[index].name = name
    remap = {index: offset for offset, index in enumerate(sorted(identities))}
    topology.atoms[1637].residue.name = "PHE"
    text = _viewer_text(len(remap), [13916], remap, topology)
    assert "mol showrep $::stage2_mol 0 off" in text
    assert "resname LIG and element C" in text
    assert "DynamicBonds" in text
    assert "label show Bonds" in text
    assert "label show Angles" in text
    assert "VDW 0.32" in text
    assert len(DYNAMIC_BONDS) == 4
    assert {row[2] for row in DYNAMIC_BONDS} == {1.25, 1.9}
    assert "::nyl12_movie_policy::apply $::stage2_mol" in text
    assert "    }\n}\nmol representation VDW" in text


def test_read_frames_accepts_zero_boxes_as_fixed_box(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Replace omitted NetCDF boxes with the declared fixed periodic cell."""

    class _Trajectory:
        coordinates = np.zeros((2, 4, 3))
        box = np.zeros((2, 6))

        def close(self) -> None:
            return None

    monkeypatch.setattr(export_module.NetCDFTraj, "open_old", lambda _: _Trajectory())
    fixed_box = np.array([10.0, 11.0, 12.0, 90.0, 90.0, 90.0])
    states = read_frames(Path("trajectory.nc"), 4, fixed_box)
    assert len(states) == 2
    assert all(np.array_equal(box, fixed_box) for _, box in states)


def test_export_stage2_writes_validated_transfer_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise export orchestration without invoking production cpptraj."""
    run = tmp_path / "stage2_00"
    system = tmp_path / "system"
    output = tmp_path / "vmd"
    run.mkdir()
    system.mkdir()
    topology_path = system / "md_solvated.prmtop"
    topology_path.write_bytes(b"topology")
    (run / "optimization.nc").write_bytes(b"trajectory")
    (run / "minimize.mdin").write_text(
        "&qmmm qmmask='@13905-13917,14559,19903-19905', /\n"
    )
    inputs = tmp_path / "inputs.json"
    fixed_box = [10.0, 11.0, 12.0, 90.0, 90.0, 90.0]
    inputs.write_text(json.dumps({"box": fixed_box}))
    cpptraj = tmp_path / "cpptraj"
    cpptraj.write_text("executable fixture")

    atom_count = 20000
    topology = _Topology([index // 10 for index in range(atom_count)])
    identities = {
        index: name
        for row in REACTIVE_LABELS
        for index, name in ((row[1], row[2]), (row[3], row[4]))
    }
    for index, name in identities.items():
        topology.atoms[index].name = name
    monkeypatch.setattr(export_module, "AmberParm", lambda _: topology)
    states = [
        (np.zeros((atom_count, 3)), np.asarray(fixed_box)),
        (np.ones((atom_count, 3)), np.asarray(fixed_box)),
    ]
    monkeypatch.setattr(export_module, "read_frames", lambda *args: states)
    kept = np.arange(atom_count)
    monkeypatch.setattr(
        export_module,
        "select_display_atoms",
        lambda *args: (kept, [101, 202]),
    )

    def fake_cpptraj(command: list[str], **kwargs: Any) -> None:
        assert kwargs["check"] is True
        input_path = Path(command[-1])
        if input_path.name == "stripped_topology.cpptraj":
            (output / "stripped.parm7").write_bytes(b"stripped topology")
        else:
            (output / "optimization.dcd").write_bytes(b"validated dcd fixture")

    monkeypatch.setattr(export_module.subprocess, "run", fake_cpptraj)
    monkeypatch.setattr(
        export_module,
        "validate_dcd",
        lambda *args: {"frames": 2, "atoms": atom_count},
    )
    manifest = export_stage2(run, system, inputs, output, cpptraj)
    assert manifest["frame_count"] == 2
    assert manifest["retained_water_residues_1based"] == [101, 202]
    assert (output / "MANIFEST.sha256").is_file()
    assert "DynamicBonds" in (output / "view_stage2.vmd").read_text()
    assert "optimizer states" in (output / "README.md").read_text()
