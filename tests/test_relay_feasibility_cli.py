"""CLI integration checks for the read-only relay feasibility report."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

import scripts.audit_nyl12_relay_feasibility as cli
from test_nyl12_relay_feasibility import _topology


def test_main_loads_arguments_writes_exclusively_and_hashes_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The CLI loads its paths, creates JSON once, and hashes source and inputs."""
    topology, coordinates = _topology()
    monkeypatch.setattr(
        cli,
        "AmberParm",
        lambda topology_path, xyz: SimpleNamespace(
            atoms=topology.atoms, residues=topology.residues, coordinates=coordinates
        ),
    )
    topology_path = tmp_path / "model.prmtop"
    coordinate_path = tmp_path / "model.rst7"
    output_path = tmp_path / "relay_feasibility.json"
    topology_path.write_bytes(b"synthetic topology input")
    coordinate_path.write_bytes(b"synthetic coordinate input")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "audit_nyl12_relay_feasibility.py",
            "--topology",
            str(topology_path),
            "--coordinates",
            str(coordinate_path),
            "--output",
            str(output_path),
        ],
    )

    cli.main()

    report = json.loads(output_path.read_text(encoding="utf-8"))
    provenance = report["provenance"]
    expected_paths = (
        topology_path,
        coordinate_path,
        Path(cli.__file__),
        cli.ROOT / "src/nyl12_covalent/relay_water.py",
        cli.ROOT / "src/nyl12_covalent/relay_feasibility.py",
    )
    for path in expected_paths:
        resolved = path.resolve()
        assert (
            provenance[str(resolved)]
            == hashlib.sha256(resolved.read_bytes()).hexdigest()
        )
    assert report["inventory"]["finite_coordinates"] is True
    assert set(report["scenarios"]) == {
        "original_all_atoms",
        "optimistic_asn_sidechain_and_beta_n_h_removed",
    }

    original = output_path.read_bytes()
    with pytest.raises(FileExistsError):
        cli.main()
    assert output_path.read_bytes() == original
