"""Focused checks for the read-only relay connectivity audit."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent import relay_connectivity_audit as audit
from nyl12_covalent import relay_connectivity_report as reporting
from nyl12_covalent import relay_connectivity_records as records
from nyl12_covalent.relay_connectivity_records import _endpoint_chemistry_passed
from nyl12_covalent.neb_umbrella_observables import frame_observables


def test_frame_observables_reports_all_h_owners_and_ambiguity() -> None:
    """Keep the pool-level owner alternatives and close competing contacts."""
    xyz = np.zeros((20, 3), dtype=float)
    xyz[1] = [1, 1, 1]
    xyz[2] = [2, 1, 1]
    xyz[3] = [1, 2, 1]
    xyz[4] = [1.6, 2, 1]
    xyz[5] = [5, 5, 5]
    xyz[6] = [5, 5, 6]
    xyz[7] = [5, 5, 4]
    roles = {
        "carbon": 2,
        "leaving_n": 3,
        "og1": 4,
        "hg1": 5,
        "nalpha": 6,
        "nalpha_h1": 7,
        "nalpha_h2": 8,
    }
    result = frame_observables(
        xyz,
        np.array([20, 20, 20, 90, 90, 90], dtype=float),
        roles,
        [5, 7, 8],
        [3, 4, 6],
        {5: 4, 7: 6, 8: 6},
    )
    assert result["hydrogen_owners"][5]["owner"] == 4
    assert result["hydrogen_owner_distances_A"][5][3] < 1.3
    assert result["ambiguous_hydrogens_within_1_30_A"] == [5]


def _row_geometry() -> dict:
    distances = {
        key: 1.0
        for key in (
            "C_OG",
            "C_NL",
            "C_O",
            "H1_N",
            "H1_OG",
            "H1_NL",
            "H2_N",
            "H2_OG",
            "H2_NL",
            "HG_N",
            "HG_OG",
            "HG_NL",
            "HN_N",
            "HN_OG",
            "HN_NL",
        )
    }
    return {
        "distances_A": distances,
        "geometry_label": "unassigned",
        "ownership": {
            name: [owner]
            for name, owner in (("H1", "N"), ("H2", "N"), ("HG", "N"), ("HN", "NL"))
        },
    }


def _accepted_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, list[dict]]:
    campaign = tmp_path / "campaign"
    model_inputs = campaign / "model/inputs"
    model_inputs.mkdir(parents=True)
    source = model_inputs / "optimized.ncrst"
    source.write_bytes(b"seed restart")
    topology = model_inputs / "system.prmtop"
    topology.write_bytes(b"topology")
    branch = campaign / "01_activation_reverse"
    folder = branch / "01_coordinate_drive/window_001"
    folder.mkdir(parents=True)
    restart = folder / "optimized.ncrst"
    restart.write_bytes(b"accepted restart")
    geometry = _row_geometry()
    process_path = folder / "process.json"
    process_path.write_text(
        json.dumps(
            {
                "exit_status": 0,
                "input_sha256": {
                    str(source.resolve()): audit.digest(source),
                    str(topology.resolve()): audit.digest(topology),
                },
            }
        )
    )
    acceptance_path = folder / "acceptance.json"
    acceptance_path.write_text(
        json.dumps(
            {
                "restart": str(restart.resolve()),
                "accepted": True,
                "native_converged": True,
                "electronic_history": {"passed": True},
                "nonrelay_hydrogen_audit": {"passed": True},
                "frozen_max_abs_change_A": 0.0,
                "geometry": geometry,
                "total_kcal_mol": 10.0,
                "bias_kcal_mol": 1.0,
                "physical_kcal_mol": 9.0,
            }
        )
    )
    seed_row = {
        "label": "window033_start",
        "restart": str(source.resolve()),
        "geometry": geometry,
        "total_kcal_mol": 8.0,
        "bias_kcal_mol": 0.0,
        "physical_kcal_mol": 8.0,
    }
    row = {
        "label": "window_001",
        "accepted": True,
        "restart": str(restart.resolve()),
        "geometry": geometry,
    }
    monkeypatch.setattr(
        records,
        "_frame_observables",
        lambda *args: {"hydrogen_pool_1based": [13906]},
    )
    return campaign, branch, [seed_row, row]


def test_accepted_rows_checks_predecessor_and_energy_decomposition(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Only a receipt-backed restart with an intact lineage enters the report."""
    campaign, branch, rows = _accepted_fixture(tmp_path, monkeypatch)
    topology_hash = audit.digest(campaign / "model/inputs/system.prmtop")
    accepted, receipts = audit._accepted_rows(
        campaign, branch, rows, topology_hash, [1], "reverse_activation"
    )
    assert [row["label"] for row in accepted] == ["window033_start", "window_001"]
    assert accepted[-1]["predecessor_restart_sha256"] == audit.digest(
        campaign / "model/inputs/optimized.ncrst"
    )
    assert accepted[-1]["physical_kcal_mol"] == 9.0
    assert len(receipts) == 1


def test_accepted_rows_rejects_wrong_predecessor_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A syntactically accepted row cannot break restart ancestry silently."""
    campaign, branch, rows = _accepted_fixture(tmp_path, monkeypatch)
    process_path = branch / "01_coordinate_drive/window_001/process.json"
    process = json.loads(process_path.read_text())
    process["input_sha256"][
        str((campaign / "model/inputs/optimized.ncrst").resolve())
    ] = "wrong"
    process_path.write_text(json.dumps(process))
    with pytest.raises(ValueError, match="Predecessor restart hash mismatch"):
        audit._accepted_rows(
            campaign,
            branch,
            rows,
            audit.digest(campaign / "model/inputs/system.prmtop"),
            [1],
            "reverse_activation",
        )


def test_accepted_unbiased_endpoint_requires_force_and_chemistry_gate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The endpoint receipt combines independent force and relay-geometry gates."""
    campaign, branch, rows = _accepted_fixture(tmp_path, monkeypatch)
    window_restart = Path(rows[-1]["restart"])
    topology = campaign / "model/inputs/system.prmtop"
    endpoint = branch / "02_unrestrained_endpoint"
    endpoint.mkdir()
    restart = endpoint / "optimized.ncrst"
    restart.write_bytes(b"unbiased endpoint restart")
    geometry = _row_geometry()
    geometry["distances_A"].update(C_OG=2.7, C_NL=1.45, C_O=1.25)
    geometry["ownership"] = {
        "H1": ["N"],
        "H2": ["N"],
        "HG": ["OG"],
        "HN": ["NL"],
    }
    process_path = endpoint / "process.json"
    process_path.write_text(
        json.dumps(
            {
                "exit_status": 0,
                "input_sha256": {
                    str(window_restart): audit.digest(window_restart),
                    str(topology.resolve()): audit.digest(topology),
                },
            }
        )
    )
    acceptance_path = endpoint / "acceptance.json"
    acceptance_path.write_text(
        json.dumps(
            {
                "restart": str(restart.resolve()),
                "accepted": True,
                "native_converged": True,
                "electronic_history": {"passed": True},
                "nonrelay_hydrogen_audit": {"passed": True},
                "frozen_max_abs_change_A": 0.0,
                "geometry": geometry,
                "total_kcal_mol": 2.0,
                "bias_kcal_mol": 0.0,
                "physical_kcal_mol": 2.0,
                "independent_energy_hartree": -700.0,
                "max_gradient_hartree_per_bohr": 8e-6,
            }
        )
    )
    rows.append(
        {
            "label": "unrestrained_endpoint",
            "accepted": True,
            "restart": str(restart.resolve()),
            "geometry": geometry,
        }
    )
    accepted, _ = records._accepted_rows(
        campaign,
        branch,
        rows,
        audit.digest(topology),
        [1],
        "reverse_activation",
    )
    endpoint_row = accepted[-1]
    assert endpoint_row["unbiased_endpoint_candidate_passed"] is True
    assert endpoint_row["maximum_gradient_hartree_per_bohr"] == 8e-6
    assert endpoint_row["independent_energy_hartree"] == -700.0


def test_read_snapshot_returns_hash_of_in_memory_document(tmp_path: Path) -> None:
    """Snapshot identity matches the bytes used to parse the JSON object."""
    path = tmp_path / "record.json"
    path.write_text('{"status":"running"}')
    value, fingerprint = audit._read_snapshot(path)
    assert value == {"status": "running"}
    assert fingerprint == audit.digest(path)


def test_refresh_watch_writes_only_audit_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One-shot refresh creates the two intended review files."""
    campaign = tmp_path / "campaign"
    for branch in audit.BRANCHES:
        folder = campaign / branch
        folder.mkdir(parents=True)
        (folder / "progress.json").write_text("{}")
        (folder / "branch_status.json").write_text('{"status":"running"}')
    report_value = {
        "audit_status": "provisional_incomplete",
        "connectivity_claim": "not established",
        "required_next_gate": "accepted unbiased endpoint",
        "branches": {},
    }
    monkeypatch.setattr(reporting, "audit_campaign", lambda _: report_value)
    output = tmp_path / "connectivity_audit"
    reporting.refresh_watch(campaign, output, once=True)
    assert sorted(path.name for path in output.iterdir()) == [
        "connectivity_audit.json",
        "connectivity_audit.md",
    ]


def test_refresh_watch_defers_changed_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A live-write race does not publish a partial or stale audit."""
    campaign = tmp_path / "campaign"
    for branch in audit.BRANCHES:
        folder = campaign / branch
        folder.mkdir(parents=True)
        (folder / "progress.json").write_text("{}")
        (folder / "branch_status.json").write_text('{"status":"running"}')

    def changed(_: Path) -> dict:
        raise audit.SnapshotChangedError("test race")

    monkeypatch.setattr(reporting, "audit_campaign", changed)
    output = tmp_path / "connectivity_audit"
    reporting.refresh_watch(campaign, output, once=True)
    assert list(output.iterdir()) == []


def test_endpoint_candidate_requires_unbiased_branch_geometry() -> None:
    """Endpoint acceptance includes bond changes and relay H ownership."""
    geometry = _row_geometry()
    geometry["distances_A"].update(C_OG=2.7, C_NL=1.45, C_O=1.25)
    geometry["ownership"] = {"H1": ["N"], "H2": ["N"], "HG": ["OG"], "HN": ["NL"]}
    assert _endpoint_chemistry_passed(geometry, "reverse_activation", 0.0)
    assert not _endpoint_chemistry_passed(geometry, "reverse_activation", 0.1)
