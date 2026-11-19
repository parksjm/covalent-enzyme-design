"""Prepare matched numerical trials from explicitly accepted shared seed manifests."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = next(
    p for p in Path(__file__).resolve().parents if (p / "src/nyl12_covalent").is_dir()
)
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts.launch_neb_spacing_trials import render_trial_input  # noqa: E402

LOGGER = logging.getLogger(__name__)
TRIALS = {
    "01_uniform_24_lbfgs": ("uniform", 24, "LBFGS"),
    "02_uniform_48_lbfgs": ("uniform", 48, "LBFGS"),
    "03_adaptive_24_lbfgs": ("adaptive", 24, "LBFGS"),
    "04_uniform_24_cg": ("uniform", 24, "CG"),
}


def prepare(name: str) -> Path:
    """Publish one accepted trial configuration without launching its engine.

    Args:
        name: Explicit comparison variant in the approved trial table.

    Returns:
        Prepared configuration path; existing destinations are an error.
    """
    home = Path(__file__).resolve().parent
    seed, count, optimizer = TRIALS[name]
    manifest_name = "manifest_24_48.json" if count == 48 else "manifest.json"
    manifest_path = home / "inputs" / seed / manifest_name
    manifest = json.loads(manifest_path.read_text())
    prediction = manifest["predictions"][str(count)]
    review = None
    if prediction["preflight_passed"] is not True:
        if count != 48:
            raise ValueError(f"Native starting-image gates failed for {name}")
        review = json.loads(
            (home / "inputs/uniform/shared_proton_review.json").read_text()
        )
        if (
            review["approved_for_seed_trial"] is not True
            or review["manifest_sha256"]
            != hashlib.sha256(manifest_path.read_bytes()).hexdigest()
            or not prediction["endpoint_within_1e_minus6_A"]
            or prediction["dropped_controls_zero_based"]
        ):
            raise ValueError("Shared-proton review does not approve this exact seed")
        for row in prediction["audit"]["images"]:
            if not row["no_compressed_H_contacts"] or row["spectator_parent_changes"]:
                raise ValueError("Overlap or spectator proton change is not approved")
            if (
                not row["strict_H_gate_passed"]
                and not row["shared_transfer20318_only_requires_root_decision"]
            ):
                raise ValueError("A different starting-image failure is not approved")
    if (
        manifest["endpoint_source_exact"] is not True
        or manifest["fixed_exterior_exact"] is not True
    ):
        raise ValueError("Seed endpoints/exterior have changed")
    controls = [Path(p) for p in manifest["images"]]
    if len(controls) != 13:
        raise ValueError("Expected thirteen seed coordinate controls")
    for path in controls:
        if (
            hashlib.sha256(path.read_bytes()).hexdigest()
            != manifest["image_sha256"][str(path)]
        ):
            raise ValueError(f"Seed input changed: {path}")
    original = home.parent / "config.json"
    config = json.loads(original.read_text())
    base_path = home.parent / "calculation/minimize.mdin"
    text = render_trial_input(base_path.read_text(), controls, count, optimizer)
    cycles = 250 if optimizer == "CG" else 16000
    if optimizer == "CG":
        if text.count("maxcycle=16000") != 1:
            raise ValueError("Expected the declared original iteration budget")
        text = text.replace("maxcycle=16000", "maxcycle=250")
    destination = home / name
    destination.mkdir()
    calculation = destination / "calculation"
    calculation.mkdir()
    (calculation / "minimize.mdin").write_text(text)
    source_folder = destination / "source_snapshot"
    source_folder.mkdir()
    sources = [
        Path(__file__),
        ROOT / "scripts/launch_neb_spacing_trials.py",
        ROOT / "scripts/prepare_neb_spacing_tests.py",
        ROOT / "src/nyl12_covalent/neb_spacing.py",
        ROOT / "scripts/run_nyl12_job.py",
        ROOT / "src/nyl12_covalent/jobs.py",
        ROOT / "scripts/watch_neb_spacing_trial.py",
        ROOT / "src/nyl12_covalent/relay_neb_progress.py",
        ROOT / "src/nyl12_covalent/neb_xyz_snapshot.py",
        ROOT / "src/nyl12_covalent/neb_progress.py",
        original,
        manifest_path,
        base_path,
    ]
    if review is not None:
        sources.append(home / "inputs/uniform/shared_proton_review.json")
    hashes = {}
    for i, path in enumerate(sources):
        copied = source_folder / f"{i:02d}_{path.name}"
        shutil.copy2(path, copied)
        hashes[str(path)] = hashlib.sha256(copied.read_bytes()).hexdigest()
    preflight = {
        "scientific_preparation_passed": True,
        "source_cycle": manifest["source_cycle"],
        "seed_manifest": str(manifest_path),
        "native_prediction": prediction,
        "shared_seed": seed,
        "regular_images": count,
        "optimizer": optimizer,
        "maxstep_bohr": 0.02,
        "qm_atoms": 102,
        "qm_charge": -1,
        "active_atoms": 4662,
        "qm_links": 8,
        "qm_real_H": 54,
        "electronic_temperature_K": 0,
        "unchanged": [
            "Hamiltonian",
            "proton_inventory",
            "12A_movable_region",
            "frozen_exterior",
            "force_tolerance",
            "fixed_endpoints",
        ],
        "source_sha256": hashes,
        "strict_geometry_preflight_passed": prediction["preflight_passed"],
        "explicit_shared_proton_seed_review": review,
        "maxcycle": cycles,
        "native_energy_evaluation_budget": 7 * cycles,
        "actual_native_startup_verification_required": True,
        "force_converged": False,
        "barrier_verified": False,
    }
    (destination / "preflight.json").write_text(json.dumps(preflight, indent=2) + "\n")
    config.update(
        parent_neb_config=str(original),
        trial_name=name,
        shared_seed=seed,
        source_seed_manifest=str(manifest_path),
        source_cycle=manifest["source_cycle"],
        reactant=str(controls[0]),
        start_restart=str(controls[0]),
        regular_images=count,
        expected_regular_images=count,
        native_neb_nimage=count + 1,
        images=[str(p) for p in controls],
        calculation=str(calculation),
        base_input=str(calculation / "minimize.mdin"),
        trial_source_sha256=hashes,
        maxstep_bohr=0.02,
        optimizer=optimizer,
    )
    launch = destination / "launch.sh"
    launch.write_text(
        "#!/bin/bash\nset -euo pipefail\ncd /data/jpc/covalent-enzyme-design\n"
        "exec /data/jpc/code/ambertools26/venv/bin/python "
        f"scripts/launch_neb_spacing_trials.py '{destination / 'config.json'}'\n"
    )
    launch.chmod(0o755)
    pending = destination / "config.pending"
    pending.write_text(json.dumps(config, indent=2) + "\n")
    pending.rename(destination / "config.json")
    LOGGER.info("Prepared %s (%d images, %s)", name, count, optimizer)
    return destination / "config.json"


def main() -> None:
    """Prepare only explicitly selected variants, failing on invalid prerequisites."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trials", choices=TRIALS, nargs="+")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    for name in args.trials:
        prepare(name)


if __name__ == "__main__":
    main()
