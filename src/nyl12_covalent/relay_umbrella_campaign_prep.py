"""Freeze a small, preparation-only 300 K N-terminal relay umbrella pilot."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Any

from nyl12_covalent.relay_umbrella_artifacts import (
    copy_checked as _copy_checked,
    digest,
    validated_parameters,
    verify_prepared as _verify_prepared,
)
from nyl12_covalent.relay_umbrella_path import (
    cumulative_arclength_A,
    equal_arclength_reference,
    full_path_windows,
    maximum_center_steps,
)
from nyl12_covalent.relay_umbrella_prepare import (
    ATOM_IDS,
    MOBILE_MASK,
    QM_MASK,
    PILOT_SEED_BASE,
    STAGES,
    pilot_stage_input,
    pilot_window_records,
    native_stage_input,
    write_json,
)

ROOT = Path(__file__).resolve().parents[2]
PARAMETER_RELEASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/umbrella_sampling/nterminal_relay_releases_20261123"
)
AMBERHOME = Path("/data/jpc/code/ambertools26/venv")
SANDER = AMBERHOME / "bin/sander"
SK_PARAMETER_DIR = AMBERHOME / "dat/slko/3ob-3-1"
EXPECTED_SANDER_SHA256 = "66323629c6d540a8038f8f8fe7557acf391f6797d88898844da708989de6a97e"
SOURCE_SNAPSHOT_FILES = (
    "scripts/prepare_relay_umbrella_pilot.py",
    "src/nyl12_covalent/__init__.py",
    "src/nyl12_covalent/neb_umbrella_inputs.py",
    "src/nyl12_covalent/restrained_qmmm_md.py",
    "src/nyl12_covalent/relay_umbrella_cv.py",
    "src/nyl12_covalent/relay_umbrella_prepare.py",
    "src/nyl12_covalent/relay_umbrella_campaign_prep.py",
    "src/nyl12_covalent/relay_umbrella_path.py",
    "src/nyl12_covalent/relay_umbrella_artifacts.py",
    "src/nyl12_covalent/relay_umbrella_diagnostics.py",
    "src/nyl12_covalent/neb_umbrella_observables.py",
    "src/nyl12_covalent/dftb_scc_validation.py",
    "src/nyl12_covalent/relay_umbrella_runtime.py",
)
EVIDENCE_SOURCES = (
    "config.json",
    "model/manifest.json",
    "model/release_acceptance.json",
    "model/source_identity.json",
    "01_activation_reverse/exploratory_plan_20261124.json",
    "02_acylation_forward/exploratory_plan_20261124.json",
)


def _validated_parameters() -> dict[str, str]:
    """Validate installed/released parameters using this module's configured paths."""
    from nyl12_covalent.relay_umbrella_artifacts import PARAMETER_RELEASE as release_root
    from nyl12_covalent.relay_umbrella_artifacts import SK_PARAMETER_DIR as parameter_root

    release = json.loads((PARAMETER_RELEASE / "prepared.json").read_text())
    expected = release.get("parameter_sha256", {})
    if set(expected) != {f"{first}-{second}.skf" for first in "HCNO" for second in "HCNO"}:
        raise ValueError("Validated release must declare all sixteen 3ob SKF hashes")
    if PARAMETER_RELEASE == release_root and SK_PARAMETER_DIR == parameter_root:
        return validated_parameters()
    observed = {}
    for name, value in expected.items():
        if digest(PARAMETER_RELEASE / "parameters_3ob" / name) != value:
            raise ValueError(f"Released parameter hash mismatch: {name}")
        if digest(SK_PARAMETER_DIR / name) != value:
            raise ValueError(f"Installed parameter hash mismatch: {name}")
        observed[name] = value
    return observed


def verify_prepared(directory: Path) -> dict[str, Any]:
    """Verify a preparation bundle's immutable local inputs and source snapshots."""
    return _verify_prepared(directory)


def _validate_model(base: Path) -> tuple[str, str]:
    """Validate exact retained QM/MM inputs and return their source strings."""
    reverse = json.loads(
        (base / "01_activation_reverse/exploratory_plan_20261124.json").read_text()
    )
    forward = json.loads(
        (base / "02_acylation_forward/exploratory_plan_20261124.json").read_text()
    )
    fixed = reverse.get("fixed_scientific_inputs", {})
    if fixed != forward.get("fixed_scientific_inputs"):
        raise ValueError("Activation/acylation branches do not share identical fixed inputs")
    if (
        fixed.get("full_atoms") != 104922
        or fixed.get("active_atoms") != 4165
        or fixed.get("qm_atoms") != 98
        or fixed.get("qm_charge") != -2
        or fixed.get("qm_method") != "DFTB3, 3ob-3-1 parameters; 0 K electronic temperature"
        or fixed.get("active_mask") != MOBILE_MASK.removeprefix(":")
        or fixed.get("qm_mask") != QM_MASK
        or fixed.get("protons_tracked") != [13906, 13907, 13917, 19937]
    ):
        raise ValueError("Retained scientific input inventory differs from reviewed model")
    source = Path(reverse["source_restart"]).resolve()
    topology = Path(reverse["topology"]).resolve()
    if source != (base / "model/inputs/optimized.ncrst").resolve():
        raise ValueError("Relay preparation points to an unexpected source restart")
    if topology != (base / "model/inputs/system.prmtop").resolve():
        raise ValueError("Relay preparation points to an unexpected topology")
    if digest(source) != reverse["source_restart_sha256"]:
        raise ValueError("Relay source restart differs from approved branch input")
    if digest(topology) != reverse["topology_sha256"]:
        raise ValueError("Relay topology differs from approved branch input")
    identity = json.loads((base / "model/source_identity.json").read_text())
    if (
        identity.get("full_atom_count") != 104922
        or identity.get("active_atom_count") != 4165
        or identity.get("qm_atom_count") != 98
        or identity.get("qm_mm_link_count") != 8
        or identity.get("qm_charge") != -2
        or identity.get("topology_sha256") != digest(topology)
    ):
        raise ValueError("Model source identity differs from accepted branch model")
    release = json.loads((base / "model/release_acceptance.json").read_text())
    if release.get("accepted") is not True or release.get("force_converged") is not True:
        raise ValueError("Source model release acceptance is missing")
    return str(source), str(topology)


def prepare_pilot(base: Path, output: Path, *, full_path: bool = False) -> Path:
    """Render 18 bounded-pilot window inputs with immutable local provenance.

    Args:
        base: Existing 12_nterminal_relay logical campaign directory.
        output: New output directory under its ``umbrella_sampling`` child.

    Returns:
        Path to the frozen preparation manifest. No native program is launched.
    """
    base = base.resolve(strict=True)
    output = output.resolve()
    umbrella_root = (base / "umbrella_sampling").resolve()
    if not output.is_relative_to(umbrella_root):
        raise ValueError("Pilot output must remain under the relay umbrella directory")
    if output.exists():
        raise FileExistsError(f"Prepared campaign already exists: {output}")
    source_restart, topology = _validate_model(base)
    parameters = _validated_parameters()
    old_mdin = (base / "model/inputs/dimer.mdin").read_text()
    windows = full_path_windows(base) if full_path else pilot_window_records(base)
    if full_path and len(windows) not in {55, 57}:
        raise ValueError(f"Unexpected densified path state count: {len(windows)}")
    if not full_path and len(windows) != 18:
        raise ValueError("Expected two reviewed 3x3 local overlap charts")
    stages = (
        [("01_heating_1ps", "heat", 2_000), ("02_equilibration_5ps", "equil", 10_000)]
        + [(f"{index:02d}_production_5ps", "production", 10_000) for index in range(3, 9)]
        if full_path
        else list(STAGES)
    )

    output.mkdir(parents=True)
    input_dir = output / "inputs"
    input_dir.mkdir()
    copied_hashes: dict[str, str] = {}
    external_hashes: dict[str, str] = {}
    copied_hashes["inputs/system.prmtop"] = _copy_checked(
        Path(topology), input_dir / "system.prmtop"
    )
    copied_hashes["inputs/source_dimer.mdin"] = _copy_checked(
        base / "model/inputs/dimer.mdin", input_dir / "source_dimer.mdin"
    )
    copied_hashes["inputs/source_seed.ncrst"] = _copy_checked(
        Path(source_restart), input_dir / "source_seed.ncrst"
    )
    copied_hashes["inputs/scientific_acceptance.json"] = _copy_checked(
        base / "model/inputs/scientific_acceptance.json",
        input_dir / "scientific_acceptance.json",
    )
    parameter_output = input_dir / "parameters_3ob"
    for name, expected in sorted(parameters.items()):
        external = SK_PARAMETER_DIR / name
        external_hashes[str(external)] = digest(external)
        copied = parameter_output / name
        observed = _copy_checked(external, copied)
        if observed != expected:
            raise ValueError(f"Copied parameter hash mismatch: {name}")
        copied_hashes[str(copied.relative_to(output))] = observed

    seed_files: dict[str, str] = {}
    seed_evidence: dict[str, str] = {}
    rows: list[dict[str, Any]] = []
    for index, window in enumerate(windows):
        seed_key = window.seed_id
        seed_relative = Path("inputs/seeds") / f"{seed_key}.ncrst"
        evidence_relative = Path("inputs/seed_evidence") / f"{seed_key}.json"
        if seed_key not in seed_files:
            seed_hash = _copy_checked(window.seed_source, output / seed_relative)
            if seed_hash != window.seed_sha256:
                raise ValueError(f"Accepted seed changed during preparation: {seed_key}")
            evidence_hash = _copy_checked(window.seed_evidence, output / evidence_relative)
            seed_files[seed_key] = str(seed_relative)
            seed_evidence[seed_key] = str(evidence_relative)
            copied_hashes[str(seed_relative)] = seed_hash
            copied_hashes[str(evidence_relative)] = evidence_hash

        stage_records = []
        prior_restart = seed_relative
        for stage_index, (stage_label, stage_kind, step_count) in enumerate(stages):
            seed = PILOT_SEED_BASE + index * 10 + stage_index
            stage_dir = output / window.label / stage_label
            stage_dir.mkdir(parents=True)
            restraint = window.bias.restraint_text()
            (stage_dir / "umbrella.RST").write_text(restraint)
            builder = native_stage_input if full_path else pilot_stage_input
            mdin = builder(
                old_mdin, stage_kind, seed, Path("umbrella.RST"), MOBILE_MASK,
                **({"steps": step_count} if full_path else {}),
            )
            (stage_dir / "run.mdin").write_text(mdin)
            copied_hashes[str((stage_dir / "umbrella.RST").relative_to(output))] = digest(
                stage_dir / "umbrella.RST"
            )
            copied_hashes[str((stage_dir / "run.mdin").relative_to(output))] = digest(
                stage_dir / "run.mdin"
            )
            stage_records.append(
                {
                    "label": stage_label,
                    "kind": stage_kind,
                    "steps": step_count,
                    "seed": seed,
                    "input": str((stage_dir / "run.mdin").relative_to(output)),
                    "restraint": str((stage_dir / "umbrella.RST").relative_to(output)),
                    "starting_restart": str(prior_restart),
                    "output_restart": str(
                        (Path(window.label) / stage_label / "restart.ncrst")
                    ),
                }
            )
            prior_restart = Path(window.label) / stage_label / "restart.ncrst"
        rows.append(
            {
                "label": window.label,
                "chart": window.chart,
                "biased_cv_names": list(window.bias.feature_names),
                "centers_A": list(window.bias.centers_A),
                "force_constants_kcal_mol_A2": list(window.bias.force_constants),
                "outer_extent_A": window.bias.outer_extent_A,
                "bias": {
                    "feature_names": list(window.bias.feature_names),
                    "centers_A": list(window.bias.centers_A),
                    "force_constants": list(window.bias.force_constants),
                    "outer_extent_A": window.bias.outer_extent_A,
                },
                "seed_source_window": window.seed_window,
                "seed_id": window.seed_id,
                "seed_source_sha256": window.seed_sha256,
                "seed_distance_in_biased_cv_space_A": window.seed_distance_A,
                "seed_offset_in_biased_CVs_A": list(window.seed_offset_A),
                "initial_bias_energy_kcal_mol": sum(
                    force * offset**2
                    for force, offset in zip(
                        window.bias.force_constants, window.seed_offset_A, strict=True
                    )
                ),
                "seed": seed_files[seed_key],
                "seed_evidence": seed_evidence[seed_key],
                "stages": stage_records,
            }
        )

    snapshot_hashes = {}
    for relative in SOURCE_SNAPSHOT_FILES:
        source = ROOT / relative
        destination = output / "source_snapshot" / relative
        snapshot_hashes[str(Path("source_snapshot") / relative)] = _copy_checked(
            source, destination
        )
    for relative in EVIDENCE_SOURCES:
        source = base / relative
        destination = output / "source_evidence" / relative
        snapshot_hashes[str(Path("source_evidence") / relative)] = _copy_checked(
            source, destination
        )
    copied_hashes.update(snapshot_hashes)

    path_3cv = equal_arclength_reference(base) if full_path else None
    physical_payload = {
        "topology_sha256": copied_hashes["inputs/system.prmtop"],
        "source_dimer_mdin_sha256": copied_hashes["inputs/source_dimer.mdin"],
        "source_seed_sha256": copied_hashes["inputs/source_seed.ncrst"],
        "parameters_3ob": parameters,
        "active_mask": MOBILE_MASK,
        "qm_mask": QM_MASK,
        "qm_charge": -2,
        "temperature_K": 300,
        "electronic_temperature_K": 0,
    }
    physical_hash = hashlib.sha256(
        json.dumps(physical_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    manifest = {
        "schema_version": 1,
        "campaign": "N-terminal relay full-path umbrella preparation" if full_path else "N-terminal relay local overlap/equilibration pilot",
        "native_calculations_launched": False,
        "stage_inputs_materialized": True,
        "launch_authorization": "not included; preparation-only user authorization",
        "parent_scientific_review_required": True,
        "temperature_K": 300,
        "production_data": True,
        "topology": "inputs/system.prmtop",
        "physical_hamiltonian_sha256": physical_hash,
        "path_3cv_A": path_3cv,
        "path_arclength_A": cumulative_arclength_A(base) if full_path else None,
        "full_path_chart_window_counts": {
            chart: sum(row["chart"] == chart for row in rows)
            for chart in ("activation", "delivery")
        } if full_path else None,
        "baseline_count_deviation": len(rows) - 55 if full_path else None,
        "full_path_max_axis_step_A": maximum_center_steps(windows) if full_path else None,
        "path_reference_images": 39 if full_path else None,
        "system": {
            "full_atoms": 104922,
            "active_atoms": 4165,
            "active_mask": MOBILE_MASK,
            "qm_atom_ids_mask": QM_MASK,
            "qm_real_atoms": 98,
            "qm_mm_link_atoms": 8,
            "qm_charge": -2,
            "qm_method": "DFTB3/3ob-3-1",
            "electronic_temperature_K": 0,
            "frozen_exterior": "Exact shared source coordinates; unchanged 4165-atom active selection.",
        },
        "atom_ids_1based": ATOM_IDS,
        "tracked_protons_1based": [13906, 13907, 13917, 19937],
        "proton_policy": "No identity walls; retain valid proton swaps and evaluate all-three CVs plus full-pool ownership diagnostics per frame.",
        "raw_features_per_frame": ["q_acyl", "q_act", "q_del"],
        "bias_convention": "Amber generalized distance; kcal/mol and Angstrom; E=k*(xi-center)^2 (no 1/2) with native linear tails at center +/-10 A.",
        "charts": ["activation=(q_acyl,q_act)", "delivery=(q_acyl,q_del)"],
        "pilot_window_count": len(rows) if not full_path else None,
        "pilot_window_ps": 4 if not full_path else None,
        "planned_production_blocks_per_window": 6 if full_path else 1,
        "planned_production_ps_per_window": 30 if full_path else 2,
        "windows": rows,
        "analysis_contract": {
            "cross_bias": "Calculate all named raw features for every frame; evaluate every window pair bias on every frame.",
            "primary": "MBAR on common physical Hamiltonian; shared relay samples must connect both charts.",
            "check": "Binned WHAM with bin-width sensitivity.",
            "patch_policy": "Joint common-overlap normalization; never add separately zeroed PMFs.",
            "no_pmf_promotion": "Full-path sampling requires all eight convergence gates before PMF interpretation." if full_path else "18-window, 4 ps/window output is an overlap/equilibration pilot only.",
            "es_search_bounds": [-0.10, 0.10],
            "acyl_search_bounds": [0.90, 1.10],
        },
        "runtime": {
            "sander": str(SANDER),
            "amberhome": str(AMBERHOME),
            "environment": {
                "AMBERHOME": str(AMBERHOME),
                "LD_LIBRARY_PATH": str(AMBERHOME / "lib"),
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
                "CUDA_VISIBLE_DEVICES": "",
            },
            "sander_host_specific_provenance_only": True,
            "parameter_sha256": parameters,
        },
        "origin_paths": {
            "topology": str(topology),
            "source_dimer_mdin": str(base / "model/inputs/dimer.mdin"),
            "source_restart": source_restart,
            "parameter_release": str(PARAMETER_RELEASE),
        },
        "external_runtime_sha256": external_hashes,
        "frozen_files_sha256": copied_hashes,
        "source_selection_note": "Accepted window restarts are nearest-seed candidates around the shared relay; root review is required. Future plateau-release evidence may motivate a separately versioned seed choice.",
    }
    prepared = output / "prepared.json"
    write_json(prepared, manifest)
    (output / "prepared.sha256").write_text(digest(prepared) + "\n")
    write_json(
        output / "preflight.json",
        {
            "prepared": True,
            "native_calculations_launched": False,
            "window_count": len(rows),
            "copied_seed_count": len(seed_files),
            "copied_file_count": len(copied_hashes),
            "external_runtime_file_count": len(external_hashes),
            "ready_for_root_review": True,
            "launch_allowed": False,
            "verification": "Call verify_prepared on this directory before any later runtime action.",
        },
    )
    return prepared


def prepare_full_path(base: Path, output: Path) -> Path:
    """Render the complete accepted-path schedule without launching MD."""
    return prepare_pilot(base, output, full_path=True)
