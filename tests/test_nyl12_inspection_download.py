"""Check curated downloads, byte integrity, and laptop command safety."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import prepare_nyl12_inspection_download as export

ROOT = Path(__file__).resolve().parents[1]
PULL = ROOT / "scripts/pull_nyl12_inspection.sh"


def test_local_repair_diagnostic_inventory(tmp_path: Path) -> None:
    protein = tmp_path / "protein_source"
    protein.mkdir()
    for name in ("pa6_l4", "maxa_maxa"):
        (protein / f"{name}_ti1_starting.pdb").write_text("fixture")
    calculation = tmp_path / "calculation"
    files = export.selected_files(
        tmp_path,
        tmp_path / "audit",
        tmp_path / "structures",
        tmp_path / "trajectories",
        protein_source=protein,
        local_repair=calculation,
    )
    prefix = "diagnostics/local_repair/pa6_l4_alpha_beta/"
    assert (
        files[prefix + "pa6_l4_dftb3_diagnostic.pdb"]
        == calculation / "01_hydrogens_water/optimized.pdb"
    )
    assert (
        files[prefix + "dftb3_acceptance.json"]
        == calculation / "01_hydrogens_water/acceptance.json"
    )
    assert (
        files[prefix + "view_local_repair.pml"]
        == tmp_path / "scripts/nyl12_local_repair_view.pml"
    )
    assert (
        files[prefix + "optimized.pdb"]
        == calculation / "01_hydrogens_water/optimized.pdb"
    )
    assert sum(name.startswith(prefix) for name in files) == 9
    files = export.selected_files(
        tmp_path,
        tmp_path / "audit",
        tmp_path / "structures",
        tmp_path / "trajectories",
        protein_source=protein,
        local_repair=calculation,
        optimization_trajectory=tmp_path / "movie",
    )
    assert (
        files[prefix + "optimization_trajectory.pdb.gz"]
        == tmp_path / "movie/optimization_trajectory.pdb.gz"
    )
    assert (
        files[prefix + "view_optimization.pml"]
        == tmp_path / "scripts/nyl12_optimization_view.pml"
    )
    with pytest.raises(ValueError, match="corresponding"):
        export.selected_files(
            tmp_path,
            tmp_path / "audit",
            tmp_path / "structures",
            tmp_path / "trajectories",
            protein_source=protein,
            optimization_trajectory=tmp_path / "movie",
        )


def test_package_integrity_and_collision(tmp_path: Path) -> None:
    """Copied bytes and both manifests agree, with no silent overwrites."""
    source = tmp_path / "original.xyz"
    source.write_text("1\nexample\nH 0 0 0\n")
    output = tmp_path / "download"
    report = export.build_package({"qm/final.xyz": source}, output)
    assert report["payload_bytes"] == source.stat().st_size
    assert (output / "qm/final.xyz").read_bytes() == source.read_bytes()
    assert json.loads((output / "provenance.json").read_text()) == report
    for line in (output / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    with pytest.raises(FileExistsError):
        export.build_package({"qm/final.xyz": source}, output)


def test_package_preserves_script_executable_mode(tmp_path: Path) -> None:
    """A copied pull script stays runnable when rsync preserves permissions."""
    source = tmp_path / "pull_nyl12_inspection.sh"
    source.write_text("#!/bin/sh\nexit 0\n")
    source.chmod(0o755)
    output = tmp_path / "download"
    export.build_package({source.name: source}, output)
    assert os.access(output / source.name, os.X_OK)


@pytest.mark.parametrize(
    "name", ["../escape", "/escape", "bad\nname", "MANIFEST.sha256", "provenance.json"]
)
def test_unsafe_paths(tmp_path: Path, name: str) -> None:
    source = tmp_path / "source"
    source.write_text("data")
    with pytest.raises(ValueError):
        export.build_package({name: source}, tmp_path / "download")
    assert not (tmp_path / "download").exists()


def test_limits_and_missing_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    output = tmp_path / "download"
    with pytest.raises(FileNotFoundError):
        export.build_package({"a": source}, output)
    source.write_text("12345")
    monkeypatch.setattr(export, "MAX_FILE_BYTES", 4)
    with pytest.raises(ValueError, match="small regular"):
        export.build_package({"a": source}, output)
    monkeypatch.setattr(export, "MAX_FILE_BYTES", 10)
    monkeypatch.setattr(export, "MAX_TOTAL_BYTES", 4)
    with pytest.raises(ValueError, match="25 MiB"):
        export.build_package({"a": source}, output)
    assert not output.exists()


def test_inventory_and_cli(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Only small inspection types enter the explicit inventory."""
    old = tmp_path / "output_nyl12_resume_20261020/nyl12_pymol_starting_models"
    old.mkdir(parents=True)
    for name in (
        "pa6_l2_ti1_starting.pdb",
        "poly_xama_maxa_ti1_starting.pdb",
        "view.pml",
        "scratch.chk",
        "raw.log",
    ):
        (old / name).write_text("fixture")
    for folder in ("docs", "manuscript"):
        (tmp_path / folder).mkdir()
        (tmp_path / folder / "note.md").write_text("note")
    inventory = export.selected_files(
        tmp_path, tmp_path / "audit", tmp_path / "optimized", tmp_path / "trajectories"
    )
    assert "protein/pa6_l4_ti1_starting.pdb" in inventory
    assert "protein/maxa_maxa_ti1_starting.pdb" in inventory
    assert (
        inventory["protein/pa6_l4_ti1_starting.pdb"].name == "pa6_l2_ti1_starting.pdb"
    )
    assert "protein/scratch.chk" not in inventory
    assert "protein/raw.log" not in inventory
    assert inventory["protein/view_models.pml"] == export.VIEWER_SOURCE
    assert (
        inventory["qm/neutral_thr_nme/conformer_01/optimized.sdf"]
        == tmp_path / "optimized/conformer_01.sdf"
    )
    source = tmp_path / "source"
    source.write_text("source")
    monkeypatch.setattr(export, "selected_files", lambda *args: {"small.txt": source})
    output = tmp_path / "export"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "export",
            "--repo",
            str(tmp_path),
            "--qm-audit",
            str(source),
            "--qm-structures",
            str(tmp_path),
            "--qm-trajectories",
            str(tmp_path),
            "--output",
            str(output),
        ],
    )
    export.main()
    assert (output / "small.txt").read_text() == "source"


def test_alternate_protein_source_and_canonical_viewer(tmp_path: Path) -> None:
    """An explicit source directory is used and incomplete sources fail."""
    protein_source = tmp_path / "repaired_models"
    protein_source.mkdir()
    for name in (
        "pa6_l4_ti1_starting.pdb",
        "maxa_maxa_ti1_starting.pdb",
    ):
        (protein_source / name).write_text(name)
    inventory = export.selected_files(
        tmp_path,
        tmp_path / "audit",
        tmp_path / "optimized",
        tmp_path / "trajectories",
        protein_source,
    )
    assert inventory["protein/pa6_l4_ti1_starting.pdb"] == (
        protein_source / "pa6_l4_ti1_starting.pdb"
    )
    assert inventory["protein/maxa_maxa_ti1_starting.pdb"] == (
        protein_source / "maxa_maxa_ti1_starting.pdb"
    )
    viewer = inventory["protein/view_models.pml"].read_text()
    assert "load pa6_l4_ti1_starting.pdb, pa6_l4" in viewer
    assert "load maxa_maxa_ti1_starting.pdb, maxa_maxa" in viewer
    assert "set valence, 0" in viewer
    assert "show sticks, resi 113+186 and chain A" in viewer
    assert "resi 113 +186" not in viewer
    empty_source = tmp_path / "empty"
    empty_source.mkdir()
    with pytest.raises(FileNotFoundError, match="No source model"):
        export.selected_files(
            tmp_path,
            tmp_path / "audit",
            tmp_path / "optimized",
            tmp_path / "trajectories",
            empty_source,
        )


def _selected_files_inputs(repo: Path) -> tuple[Path, Path, Path, Path, Path]:
    """Build the two required protein-source files and inventory arguments."""
    proteins = repo / "proteins"
    proteins.mkdir(parents=True)
    (proteins / "pa6_l4_ti1_starting.pdb").write_text("pa6", encoding="utf-8")
    (proteins / "maxa_maxa_ti1_starting.pdb").write_text("maxa", encoding="utf-8")
    return (
        repo,
        repo / "audit.json",
        repo / "optimized",
        repo / "trajectories",
        proteins,
    )


def test_selected_files_adds_exact_af3_water_artifacts_per_attempt(
    tmp_path: Path,
) -> None:
    """Each explicit AF3 attempt contributes only its report, model, and viewer."""
    args = _selected_files_inputs(tmp_path / "repo")
    result = tmp_path / "attempt_01" / "result"
    result.mkdir(parents=True)
    names = (
        "report.json",
        "pa6_l4_alpha_beta_water_test.cif",
        "view_af3_water_test.pml",
    )

    inventory = export.selected_files(
        *args[:4], protein_source=args[4], af3_audits=[result]
    )

    expected = {
        f"diagnostics/af3_water/attempt_01/{name}": result / name for name in names
    }
    assert {name: inventory[name] for name in expected} == expected
    assert not any(
        name.startswith("diagnostics/af3_water/attempt_01/")
        for name in inventory
        if name not in expected
    )


def test_selected_files_default_inventory_has_no_af3_water_additions(
    tmp_path: Path,
) -> None:
    """Omitting optional AF3 audits preserves the previous inventory exactly."""
    args = _selected_files_inputs(tmp_path / "repo")

    default = export.selected_files(*args[:4], protein_source=args[4])
    explicit_empty = export.selected_files(
        *args[:4], protein_source=args[4], af3_audits=[]
    )

    assert default == explicit_empty
    assert not any(name.startswith("diagnostics/af3_water/") for name in default)


def test_selected_files_rejects_duplicate_af3_attempt_labels(
    tmp_path: Path,
) -> None:
    """Different directories with the same attempt label cannot collide."""
    args = _selected_files_inputs(tmp_path / "repo")
    first = tmp_path / "first" / "attempt_01" / "result"
    second = tmp_path / "second" / "attempt_01" / "result"
    first.mkdir(parents=True)
    second.mkdir(parents=True)

    with pytest.raises(ValueError, match="Duplicate AF3 audit label: attempt_01"):
        export.selected_files(
            *args[:4], protein_source=args[4], af3_audits=[first, second]
        )


def test_alias_provenance_keeps_substrate_and_product_distinct(
    tmp_path: Path,
) -> None:
    """Export aliases retain legacy IDs and measured product identities."""
    source = tmp_path / "source.txt"
    source.write_text("source")
    report = export.build_package({"source.txt": source}, tmp_path / "package")
    aliases = report["model_aliases"]
    assert aliases["pa6_l4"]["source_ids"] == ["pa6_l2", "pa6_l4"]
    assert aliases["pa6_l4"]["product"] == "L2"
    assert aliases["maxa_maxa"]["source_ids"] == [
        "poly_xama_maxa",
        "maxa_maxa",
    ]
    assert aliases["maxa_maxa"]["product"] == "MAXA"


def test_rsync_dry_run_and_arguments(tmp_path: Path) -> None:
    """A laptop destination with spaces stays one argument; no deletion flag."""
    binary = tmp_path / "bin"
    binary.mkdir()
    captured = tmp_path / "args"
    for name, body in {
        "rsync": '#!/bin/sh\nprintf "%s\\n" "$@" > "$CAPTURE_ARGS"\n',
        "ssh": "#!/bin/sh\nexit 0\n",
        "shasum": "#!/bin/sh\nexit 0\n",
    }.items():
        path = binary / name
        path.write_text(body)
        path.chmod(0o755)
    environment = {
        **os.environ,
        "PATH": f"{binary}:{os.environ['PATH']}",
        "CAPTURE_ARGS": str(captured),
    }
    destination = tmp_path / "laptop files"
    for dry_run in (True, False):
        argv = [
            "bash",
            str(PULL),
            *(["--dry-run"] if dry_run else []),
            "jpc@violet",
            str(destination),
        ]
        subprocess.run(argv, env=environment, check=True, capture_output=True)
        arguments = captured.read_text().splitlines()
        assert ("--dry-run" in arguments) == dry_run
        assert "--delete" not in arguments
        assert arguments[-1] == f"{destination}/"
        assert arguments[-2].endswith("output_laptop_inspection_20261021/")
        assert destination.exists() is not dry_run


def test_rsync_executable_uses_current_user_and_violet_by_default(
    tmp_path: Path,
) -> None:
    """Direct invocation needs no args and retains the checked default path."""
    binary = tmp_path / "bin"
    binary.mkdir()
    captured = tmp_path / "args"
    for name, body in {
        "rsync": '#!/bin/sh\nprintf "%s\\n" "$@" > "$CAPTURE_ARGS"\n',
        "ssh": "#!/bin/sh\nexit 0\n",
        "shasum": "#!/bin/sh\nexit 0\n",
    }.items():
        path = binary / name
        path.write_text(body)
        path.chmod(0o755)
    environment = {
        **os.environ,
        "PATH": f"{binary}:{os.environ['PATH']}",
        "CAPTURE_ARGS": str(captured),
    }
    current_user = subprocess.run(
        ["id", "-un"], check=True, text=True, capture_output=True
    ).stdout.strip()
    subprocess.run(
        [str(PULL)],
        cwd=tmp_path,
        env=environment,
        check=True,
        capture_output=True,
    )
    arguments = captured.read_text().splitlines()
    assert arguments[-2] == (
        f"{current_user}@violet:/data/jpc/covalent-enzyme-design/"
        "output_laptop_inspection_20261021/"
    )
    assert arguments[-1] == "./nyl12-inspection-20261021/"
    assert "--dry-run" not in arguments
    assert "-rtpvz" in arguments
    assert (tmp_path / "nyl12-inspection-20261021").is_dir()


@pytest.mark.parametrize("args", [["-evil"], ["host;touch foo"], ["goodhost", "-bad"]])
def test_rsync_rejects_invalid_arguments(args: list[str]) -> None:
    result = subprocess.run(["bash", str(PULL), *args], capture_output=True)
    assert result.returncode == 2
