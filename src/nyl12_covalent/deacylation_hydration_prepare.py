"""Freeze inputs and source for the bounded hydration pilot."""

from __future__ import annotations
import json
from pathlib import Path
import shutil
from .deacylation_hydration import build_input, digest, write_json


def prepare(campaign: Path, audit: Path, repo: Path, seed: int) -> None:
    """Freeze audited input files, implementation closure and stage inputs.

    Args:
        campaign: New prepared-campaign directory.
        audit: Immutable endpoint-water audit JSON.
        repo: Repository containing source and journal script.
        seed: Explicit Langevin random seed.
    """
    record = json.loads(audit.read_text())
    files = record["files"]
    acceptance = record["independent_acceptance"]
    if (
        acceptance["force_converged"] is not True
        or acceptance["electronic_history_passed"] is not True
        or not 0 <= acceptance["force_max_hartree_per_bohr"] <= 1e-5
    ):
        raise ValueError("Source endpoint lacks independent force/SCC acceptance")
    model = record["native_model"]
    expected = {
        "atom_count": 104922,
        "qm_atoms": 98,
        "qmcharge": -2,
        "link_atoms": 8,
        "active_atoms": 4165,
        "active_shell_A": 12,
    }
    if any(model[key] != value for key, value in expected.items()):
        raise ValueError("Source model identity differs from audited pilot")
    for item in files.values():
        if digest(Path(item["path"])) != item["sha256"]:
            raise ValueError("Audit source hash mismatch")
    active = record["native_model"]["active_mask"]
    source = Path(files["endpoint_native_mdin"]["path"]).read_text()
    if f"qmmask='{model['qmmask']}'" not in source:
        raise ValueError("Native input differs from audited QM mask")
    texts = {
        stage: build_input(source, active, stage, seed) for stage in ("heat", "equil")
    }
    campaign.mkdir(parents=True, exist_ok=False)
    inputs = campaign / "inputs"
    inputs.mkdir()
    copied = {}
    for key in ("endpoint_restart", "source_topology", "endpoint_native_mdin"):
        dst = (
            inputs
            / {
                "endpoint_restart": "endpoint.ncrst",
                "source_topology": "system.prmtop",
                "endpoint_native_mdin": "source.mdin",
            }[key]
        )
        shutil.copy2(files[key]["path"], dst)
        copied[key] = {"path": str(dst.resolve()), "sha256": digest(dst)}
    for stage, text in texts.items():
        (inputs / f"{stage}.mdin").write_text(text)
    snapshot = campaign / "source_snapshot"
    package = snapshot / "nyl12_covalent"
    package.mkdir(parents=True)
    for name in (
        "__init__.py",
        "deacylation_hydration.py",
        "deacylation_hydration_prepare.py",
        "neb_umbrella_inputs.py",
        "dftb_scc_validation.py",
    ):
        shutil.copy2(repo / "src/nyl12_covalent" / name, package / name)
    shutil.copy2(
        repo / "scripts/run_deacylation_hydration.py",
        snapshot / "run_deacylation_hydration.py",
    )
    snapshot_hashes = {
        str(p.relative_to(snapshot)): digest(p) for p in snapshot.rglob("*.py")
    }
    input_hashes = {str(p.relative_to(campaign)): digest(p) for p in inputs.iterdir()}
    config = {
        "audit_path": str(audit.resolve()),
        "audit_sha256": digest(audit),
        "sources": copied,
        "input_sha256": input_hashes,
        "source_snapshot_sha256": snapshot_hashes,
        "sander": files["sander"],
        "active_mask": active,
        "qmmask": model["qmmask"],
        "seed": seed,
        "method": "original 98-QM DFTB3/MM, charge -2; no water promotion; restrained only by frozen original exterior",
        "scope": "1 ps heating and 5 ps NVT hydration reconnaissance; not reaction or converged sampling",
        "journal_repo": str(repo.resolve()),
        "dftb_parameters_sha256": {
            str(path): digest(path)
            for a in "CHNO"
            for b in "CHNO"
            for path in [
                Path("/data/jpc/code/ambertools26/venv/dat/slko/3ob-3-1")
                / f"{a}-{b}.skf"
            ]
        },
    }
    write_json(campaign / "config.json", config)
