"""Tests for modular PA6 L4 campaign preparation."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import nyl12_covalent.modular_pa6_md_setup as md_setup
from nyl12_covalent.modular_pa6_md_protocol import (
    cpptraj_analysis,
    mdin_files,
    replica_runner,
)
from nyl12_covalent.modular_pa6_md_setup import (
    SystemSpec,
    _ion_counts,
    _require_clean_leap,
    _write_leap_input,
    _write_leap_pdb,
    _write_queues,
    prepare_campaign,
)
from nyl12_covalent.modular_pa6_topology import (
    PA6_SEQUENCE,
    sha256,
    validate_release,
    validate_topology,
)


def test_protocol_preserves_historical_lengths_and_modular_masks() -> None:
    files = mdin_files(810_117)
    assert len(files) == 7
    assert "nstlim=50000000" in files["07_prod_100ns.in"]
    assert "nstlim=50000" in files["03_heat_nvt.in"]
    assert "restraintmask='!:WAT,Na+,Cl- & !@H='" in files["03_heat_nvt.in"]
    analysis = cpptraj_analysis()
    assert ":1313-1316" in analysis
    assert ":1314@C :1314@O" in analysis
    assert ":1314@C5 :1315@N :1314@C :1314@O" in analysis


def test_protocol_rejects_bad_seed_and_relative_binary() -> None:
    with pytest.raises(ValueError, match="positive"):
        mdin_files(0)
    with pytest.raises(ValueError, match="absolute"):
        replica_runner(Path("pmemd.cuda"))
    runner = replica_runner(Path("/opt/amber/bin/pmemd.cuda"))
    assert runner.count("if [[ ! -s") == 7
    assert "07_prod_100ns.nc" in runner


def test_leap_inputs_select_mbondi3_release_and_explicit_ions(tmp_path: Path) -> None:
    count = tmp_path / "count.in"
    final = tmp_path / "final.in"
    _write_leap_input(count, None)
    _write_leap_input(final, (42, 44))
    assert "source leaprc.pa6_shared_p6i" in count.read_text()
    assert "source leaprc.gaff\n" not in count.read_text()
    assert "addIonsRand complex Na+ 42 Cl- 44 6.0" in final.read_text()


def test_pdb_preprocessing_restores_boundaries_and_n_terminal_name(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.pdb"
    source.write_text(
        "ATOM      1  H   ALA A   1       0.000   0.000   0.000  1.00  0.00           H  \n"
        "ATOM      2  OXT PRO A1312       1.000   0.000   0.000  1.00  0.00           O  \n"
        "HETATM    3  OXT P6C A1316       2.000   0.000   0.000  1.00  0.00           O  \n"
        "HETATM    4  H52 P6C A1316       3.000   0.000   0.000  1.00  0.00           H  \n"
        "END\n"
    )
    destination = tmp_path / "leap.pdb"
    _write_leap_pdb(source, destination)
    text = destination.read_text()
    assert " H1 " in text
    assert text.count("TER\n") == 2
    assert text.index("H52 P6C") < text.rindex("TER\n")
    assert text.endswith("END\n")


def test_ion_counts_and_leap_error_gate() -> None:
    sodium, chloride, pairs = _ion_counts(-2.0, 1_000_000.0)
    assert pairs == 90
    assert (sodium, chloride) == (92, 90)
    assert _ion_counts(-2.009, 1_000_000.0) == (92, 90, 90)
    with pytest.raises(ValueError, match="not integral"):
        _ion_counts(0.2, 1_000_000.0)
    _require_clean_leap("Errors = 0")
    with pytest.raises(ValueError, match="LEaP"):
        _require_clean_leap("Errors = 2")


def test_queue_assignment_is_balanced_and_portable(tmp_path: Path) -> None:
    replicas = [f"systems/wt/replica_{index:02d}" for index in range(1, 13)]
    _write_queues(tmp_path, replicas)
    assignments = json.loads((tmp_path / "host_assignment.json").read_text())
    assert set(assignments) == {
        "violet_gpu0",
        "violet_gpu1",
        "azalea_gpu0",
        "azalea_gpu1",
    }
    assert {len(entries) for entries in assignments.values()} == {3}
    assert all(
        not entry.startswith("/")
        for entries in assignments.values()
        for entry in entries
    )


def test_prepare_campaign_writes_manifest_and_replicas(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    release = tmp_path / "release"
    release.mkdir()
    (release / "manifest.json").write_text("{}\n")
    tleap = tmp_path / "tleap"
    tleap.write_text("binary\n")
    pdb = tmp_path / "wt.pdb"
    pdb.write_text("END\n")
    conversion = tmp_path / "conversion.json"
    conversion.write_text("{}\n")
    monkeypatch.setattr(
        md_setup,
        "validate_release",
        lambda _path: {"release_id": "test_release"},
    )
    monkeypatch.setattr(
        md_setup,
        "_prepare_system",
        lambda spec, _release, _output, _tleap: {"label": spec.label},
    )
    output = tmp_path / "campaign"
    manifest = prepare_campaign(
        [SystemSpec("wt", pdb, conversion)],
        release,
        output,
        tleap,
        Path("/opt/amber/bin/pmemd.cuda"),
        replicas=2,
    )
    assert manifest["release_id"] == "test_release"
    assert manifest["replica_directories"] == [
        "systems/wt/replica_01",
        "systems/wt/replica_02",
    ]
    assert (output / "systems/wt/replica_01/07_prod_100ns.in").is_file()
    saved = json.loads((output / "manifest.json").read_text())
    assert saved["replicas_per_system"] == 2
    with pytest.raises(FileExistsError):
        prepare_campaign(
            [SystemSpec("wt", pdb, conversion)],
            release,
            output,
            tleap,
            Path("/opt/amber/bin/pmemd.cuda"),
        )


def test_release_hash_validation(tmp_path: Path) -> None:
    artifact = tmp_path / "P6N.lib"
    artifact.write_text("immutable\n")
    manifest = {
        "md_ready": True,
        "chemistry": {"sequence": list(PA6_SEQUENCE)},
        "artifact_sha256": {artifact.name: sha256(artifact)},
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    assert validate_release(tmp_path) == manifest
    artifact.write_text("changed\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        validate_release(tmp_path)


def test_topology_validation_accepts_exact_release(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    residues = [SimpleNamespace(name="ALA", atoms=[]) for _ in range(1312)]
    pa6_residues = []
    for template in PA6_SEQUENCE:
        atoms = []
        atom_lines = [
            ("N", "n3", -0.25),
            ("C", "c", 0.25),
        ]
        (tmp_path / f"{template}.mol2").write_text(
            "@<TRIPOS>ATOM\n"
            + "\n".join(
                f"{index} {name} 0 0 0 {atom_type} 1 {template} {charge}"
                for index, (name, atom_type, charge) in enumerate(atom_lines, start=1)
            )
            + "\n@<TRIPOS>BOND\n"
        )
        residue = SimpleNamespace(name=template, atoms=atoms)
        for name, atom_type, charge in atom_lines:
            atoms.append(
                SimpleNamespace(
                    name=name, type=atom_type, charge=charge, residue=residue
                )
            )
        pa6_residues.append(residue)
    residues.extend(pa6_residues)
    bonds = [
        SimpleNamespace(
            atom1=pa6_residues[index].atoms[1], atom2=pa6_residues[index + 1].atoms[0]
        )
        for index in range(3)
    ]
    topology = SimpleNamespace(
        residues=residues,
        atoms=[atom for residue in residues for atom in residue.atoms],
        bonds=bonds,
        parm_data={"RADIUS_SET": ["modified Bondi radii (mbondi3)"]},
    )
    topology_path = tmp_path / "test.prmtop"
    topology_path.write_text("topology\n")
    monkeypatch.setattr(
        "nyl12_covalent.modular_pa6_topology.parmed.load_file",
        lambda _path: topology,
    )
    report = validate_topology(topology_path, tmp_path, solvated=True)
    assert report["total_charge"] == pytest.approx(0.0)
    assert report["radius_set"].endswith("(mbondi3)")
