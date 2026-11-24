"""Tests for derivation of portable design MD inputs."""

import json
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import design_md_inputs as adapter


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    source = tmp_path / "source"
    for directory in ("system", "minimization", "inputs/parameters"):
        (source / directory).mkdir(parents=True)
    release = source / "inputs/parameters"
    (release / "manifest.json").write_text('{"release_id":"release"}')
    for name in (
        "system/complex.prmtop",
        "minimization/min.rst7",
        "minimization/endpoint_complex.prmtop",
        "minimization/endpoint_complex.rst7",
    ):
        (source / name).write_text(name)
    mdin = source / "minimization/min.in"
    mdin.write_text("ntr=0, nmropt=0")
    (source / "minimization/min.out").write_text("native output")
    top = source / "system/complex.prmtop"
    rst = source / "minimization/min.rst7"
    argv = ["engine", "-p", str(top), "-r", str(rst), "-i", str(mdin)]
    command = {
        "event": "planned",
        "argv": argv,
        "input_sha256": {str(p): adapter._digest(p) for p in (top, mdin)},
    }
    (source / "command_log.jsonl").write_text(json.dumps(command) + "\n")
    process = {"status": "completed", "events": [{"argv": argv, "returncode": 0}]}
    (source / "process_acceptance.json").write_text(json.dumps(process))
    acceptance = {
        "status": "minimization_complete",
        "force_converged": True,
        "process_exit_code": 0,
        "parameter_release_id": "release",
        "final_rms_gradient": 0.04,
        "method": {
            "protein_force_field": "ff19SB",
            "water_model": "TIP3P",
            "ligand_force_field": "released modular PA6 RESP/GAFF2",
            "protonation": "Neutral beta Thr; ASP272/ASP274 anionic at all sites.",
            "restraints": "none; ntr=0, nmropt=0, ntc=1, ntf=1",
            "drms_target_kcal_mol_angstrom": 0.05,
        },
        "heavy_export": {
            "topology_sha256": adapter._digest(
                source / "minimization/endpoint_complex.prmtop"
            ),
            "restart_sha256": adapter._digest(
                source / "minimization/endpoint_complex.rst7"
            ),
        },
    }
    (source / "acceptance.json").write_text(json.dumps(acceptance))
    residues = [SimpleNamespace(name="ALA") for _ in range(1328)]
    for i in range(4):
        for index, name in (
            (113 + 232 * i, "PHE"),
            (156 + 232 * i, "LYS"),
            (186 + 232 * i, "ASN"),
            (929 + 96 * i, "THR"),
            (968 + 96 * i, "ASP"),
            (970 + 96 * i, "ASP"),
        ):
            residues[index - 1].name = name
    atoms = [SimpleNamespace(name="CA", residue=residues[0], charge=0.0, idx=0)]
    structure = SimpleNamespace(
        residues=residues,
        atoms=atoms,
        coordinates=np.array([[1.0, 2.0, 3.0]]),
        box=np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0]),
    )
    monkeypatch.setattr(adapter.parmed, "load_file", lambda *args, **kwargs: structure)
    monkeypatch.setattr(
        adapter,
        "validate_four_ligand_topology",
        lambda *args, **kwargs: {"release_id": "release"},
    )
    monkeypatch.setattr(adapter, "_final_results", lambda text: (-1.0, 0.04))
    return source, release, tmp_path / "output"


def test_derived_dataset_preserves_receipts(fixture):
    source, release, output = fixture
    result = adapter.prepare_design_dataset(source, output, "design", release)
    entry = result["datasets"]["design"]
    assert len(entry["protected_residues"]) == 24
    assert entry["terminal_protonation"] == "NH2"
    assert (output / "design/evidence/acceptance.json").read_bytes() == (
        source / "acceptance.json"
    ).read_bytes()
    receipt = json.loads((output / "design/minimization_acceptance.json").read_text())
    assert receipt["derived"] and receipt["coordinates_sha256"] == adapter._digest(
        output / "design/start.rst7"
    )
    assert json.loads((output / "datasets.json").read_text()) == result


@pytest.mark.parametrize(
    "mutation", ["hash", "force", "microstate", "process", "endpoint"]
)
def test_rejects_invalid_native_evidence(fixture, mutation):
    source, release, output = fixture
    acceptance_path = source / "acceptance.json"
    acceptance = json.loads(acceptance_path.read_text())
    if mutation == "hash":
        (source / "system/complex.prmtop").write_text("changed")
    elif mutation == "endpoint":
        (source / "minimization/endpoint_complex.rst7").write_text("changed")
    elif mutation == "process":
        (source / "process_acceptance.json").write_text(
            '{"status":"completed","events":[{"returncode":1}]}'
        )
    else:
        if mutation == "force":
            acceptance["force_converged"] = False
        else:
            acceptance["method"]["protonation"] = "ASH274"
        acceptance_path.write_text(json.dumps(acceptance))
    with pytest.raises(ValueError):
        adapter.prepare_design_dataset(source, output, "design", release)
    assert not output.exists()
