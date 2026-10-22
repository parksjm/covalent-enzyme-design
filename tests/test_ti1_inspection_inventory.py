from pathlib import Path

import pytest

from scripts import prepare_nyl12_inspection_download as package


def test_ti1_reference_files_returns_only_explicit_inventory(tmp_path: Path) -> None:
    audit = tmp_path / "audit"
    files = package.ti1_reference_files(audit)

    assert (
        files["qm/ti1_bonded_reference_progress/report.json"] == audit / "report.json"
    )
    assert files["qm/ti1_bonded_reference_progress/README.md"] == audit / "README.md"
    case_names = (
        "attack_R_cap_phi_minus120",
        "attack_R_cap_phi_minus60",
        "attack_S_cap_phi_minus120",
        "attack_S_cap_phi_minus60",
    )
    artifact_names = (
        "trajectory.xyz",
        "trajectory.json",
        "latest_geometry.xyz",
        "latest_geometry.sdf",
        "atom_map.json",
        "geometry_checks.json",
        "trajectory_metrics.csv",
        "view_progress.pml",
        "input.gjf",
    )
    expected = {
        f"qm/ti1_bonded_reference_progress/{case}/{artifact}"
        for case in case_names
        for artifact in artifact_names
    }
    expected.update(
        {
            "qm/ti1_bonded_reference_progress/report.json",
            "qm/ti1_bonded_reference_progress/README.md",
        }
    )
    assert set(files) == expected
    assert not any(path.endswith((".log", ".chk")) for path in files)


def test_selected_files_adds_optional_ti1_inventory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(package, "MODEL_ALIASES", {})
    protein_source = tmp_path / "protein_source"
    protein_source.mkdir()
    audit = tmp_path / "audit"

    files = package.selected_files(
        tmp_path / "repo",
        tmp_path / "audit_source",
        tmp_path / "structures",
        tmp_path / "trajectories",
        protein_source=protein_source,
        ti1_reference_audit=audit,
    )

    assert (
        files["qm/ti1_bonded_reference_progress/report.json"] == audit / "report.json"
    )
    assert files[
        "qm/ti1_bonded_reference_progress/attack_R_cap_phi_minus120/input.gjf"
    ] == (audit / "attack_R_cap_phi_minus120" / "input.gjf")


def test_build_package_rejects_missing_requested_ti1_artifact(
    tmp_path: Path,
) -> None:
    audit = tmp_path / "audit"
    output = tmp_path / "package"

    with pytest.raises(FileNotFoundError):
        package.build_package(package.ti1_reference_files(audit), output)

    assert not output.exists()
