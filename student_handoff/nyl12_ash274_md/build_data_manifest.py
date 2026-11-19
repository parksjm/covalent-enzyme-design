"""Build a curated data allowlist from explicit existing local artifacts only."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import json
import logging
from pathlib import Path

from inventory_common import (
    BASE,
    ENTRIES,
    LEGACY,
    MISSING_OPTIONAL,
    OUT,
    RELEASE,
    ROOT,
    add,
    sha,
    write,
)
from inventory_datasets import dataset


def main() -> None:
    """Generate explicitallowlist and portable threeWTdataset artifactreceipts."""
    datasets = {}
    for identifier, root, topology, state, reps in (
        (
            "wt_asp272_ash274",
            BASE / "07_ash274_comparison",
            "ash274_solvated.prmtop",
            "ASP272/ASH274",
            ["ash274_replica_1", "ash274_replica_2"],
        ),
        (
            "wt_asp272_asp274",
            BASE,
            "md_solvated.prmtop",
            "ASP272/ASP274",
            ["neutral_replica_01", "neutral_replica_02"],
        ),
        (
            "wt_ash272_asp274",
            BASE / "06_ash272_comparison",
            "ash272_solvated.prmtop",
            "ASH272/ASP274",
            ["ash272_replica_1", "ash272_replica_2"],
        ),
    ):
        datasets[identifier] = dataset(identifier, root, topology, state, reps)
    write(OUT / "datasets/manifest.json", dict(schema_version=1, datasets=datasets))
    add(
        OUT / "datasets/manifest.json",
        "datasets/manifest.json",
        "portable_dataset_manifest",
    )
    manifest = json.loads((RELEASE / "manifest.json").read_text())
    for p in sorted(RELEASE.iterdir()):
        if p.is_file():
            if (
                p.name in manifest["artifact_sha256"]
                and sha(p) != manifest["artifact_sha256"][p.name]
            ):
                raise ValueError("Authoritative release artifact hash mismatch")
            add(p, f"parameters/{RELEASE.name}/{p.name}", "released_PA6_L4_RESP_GAFF2")
    for name in [
        "raw_model30.pdb",
        "complex_heavy_modular.pdb",
        "conversion.json",
        "whole_l4_source.mol2",
        "whole_l4_context_mapping.json",
    ]:
        add(
            BASE / "00_system" / name,
            f"original_inputs/model30/{name}",
            "original_WT_model_and_atom_mapping",
        )
    for p in [
        LEGACY / "Nyl12_refine13_no_tail.cif",
        LEGACY / "Nyl12_tetramer-symgen_no_tail.cif",
        LEGACY / "LICENSE",
    ]:
        add(
            p,
            f"original_inputs/templates/{p.name}",
            "original_structure_template",
            required=False,
        )
    y = (
        LEGACY
        / "output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/0_boltz/input/yaml_incremental_msa"
    )
    for p in sorted(y.glob("pa6_l4z_l2_nyl12_*.yaml")):
        add(
            p,
            f"original_inputs/historical_WT_mutant_references/{p.name}",
            "historical_mutation_inputs_not_new_AS Hsystems",
            required=False,
        )
    analysis = BASE / "08_asp_pair_md_comparison"
    for p in [
        analysis / "inputs/manifest.json",
        analysis / "job_retry/command.json",
        *sorted((analysis / "analysis_attempt2").glob("*.json")),
        *sorted((analysis / "analysis_attempt2").glob("*.csv")),
    ]:
        add(
            p, "analysis/existing_md_comparison/" + p.name, "completed_replica_analysis"
        )
    for p in [
        LEGACY / "pixi.toml",
        LEGACY / "pixi.lock",
        ROOT / "pyproject.toml",
        ROOT / "uv.lock",
    ]:
        add(
            p,
            "environment/"
            + ("historical_analysis_" + p.name if p.parent == LEGACY else p.name),
            "exact_environment_specification",
        )
    for name in [
        "run_flipped_model30_mm.py",
        "prepare_flipped_model30_ash272.py",
        "run_ash272_overnight.py",
        "run_model30_md_replica.py",
        "continue_model30_md.py",
        "compare_asp_pair_md.py",
    ]:
        add(
            ROOT / "scripts" / name,
            "preparation_source/scripts/" + name,
            "existing_project_recipe",
        )
    for name, source in [
        ("pilot_md_protocol.py", ROOT / "src/nyl12_covalent/pilot_md_protocol.py"),
        *[
            (name, BASE / "source_snapshot/src/nyl12_covalent" / name)
            for name in (
                "modular_pa6_md_setup.py",
                "modular_pa6_topology.py",
                "raw_pa6_modular.py",
                "modular_pa6_pdb.py",
            )
        ],
        (
            "active_site_geometry.py",
            ROOT / "src/nyl12_covalent/active_site_geometry.py",
        ),
        (
            "modular_es_snapshot_selection.py",
            ROOT / "src/nyl12_covalent/modular_es_snapshot_selection.py",
        ),
    ]:
        add(
            source,
            "preparation_source/src/nyl12_covalent/" + name,
            "retained_existing_project_recipe",
        )
    raw = []
    for row in json.loads((analysis / "analysis_attempt2/summary.json").read_text())[
        "replicas"
    ]:
        for trajectory in row["trajectory_sources"]:
            p = Path(trajectory["path"])
            raw.append(
                dict(
                    source=str(p),
                    destination=f"raw_trajectories/{row['state']}/{row['name']}/{p.name}",
                    profile="raw_reference_only",
                    required=False,
                    bytes=p.stat().st_size,
                    sha256=None,
                    category="excluded_bulk_MD_trajectory",
                    known_source_receipt_sha256=sha(
                        analysis / "analysis_attempt2/summary.json"
                    ),
                )
            )
    qm = BASE / "07_ash274_comparison/09_water_relay_acylation"
    for label, p in [
        ("campaign_README.md", qm / "README.md"),
        ("campaign_status.json", qm / "campaign_status.json"),
        ("acyl_endpoint_assessment.json", qm / "02_acyl/endpoint_assessment.json"),
        (
            "ES_endpoint_assessment.json",
            qm / "03_reverse_path/endpoint_assessment.json",
        ),
    ]:
        if p.is_file():
            add(p, "analysis/QM_summary/" + label, "QM_summary", required=False)
        else:
            MISSING_OPTIONAL.append(str(p))
    for label, pointer in [
        ("reverse_coordinate_path", qm / "progress/latest"),
        ("native_NEB", qm / "04_neb/progress/latest"),
    ]:
        d = pointer.resolve(strict=True)
        for name in [
            "progress.json",
            "images.csv",
            "images.json",
            "snapshot_energies.csv",
            "snapshot_coordinates.csv",
        ]:
            if (d / name).is_file():
                add(
                    d / name,
                    f"inspection/{label}/{name}",
                    "selected_QM_inspection",
                    profile="inspection",
                    required=False,
                )
        v = d / "vmd"
        for name in [
            "view_path.vmd",
            "path.parm7",
            "path.dcd",
            "path.csv",
            "manifest.json",
            "MANIFEST.sha256",
        ]:
            add(
                v / name,
                f"inspection/{label}/vmd/{name}",
                "selected_portable_QM_VMD_bundle",
                profile="inspection",
                required=False,
            )
    counts = Counter()
    size = Counter()
    for e in [*ENTRIES.values(), *raw]:
        counts[e["profile"]] += 1
        size[e["profile"]] += e["bytes"]
    report = dict(
        schema_version=1,
        created_utc=datetime.now(timezone.utc).isoformat(),
        entries=[*ENTRIES.values(), *raw],
        datasets=datasets,
        parameter_release=dict(
            id=RELEASE.name,
            source=str(RELEASE),
            repository="https://github.com/parksjm/polyamide-amber-ff.git",
            observed_repository_commit="91d180c730eee4d3757bf01c3f4bcb32b556d452",
            declared_release_revision=manifest["source_revision"],
            license_file_at_repository_or_release_root=False,
        ),
        summary=dict(entry_count_by_profile=dict(counts), bytes_by_profile=dict(size)),
        excluded=[
            "licensed executables/software installations",
            "Gaussian checkpoints/ESP training bulk",
            "fullrawMDtrajectories",
            "bulkQMMMevaluationtrees",
            "historicalmutants notnewvalidatedASHsystems",
        ],
        warnings=[
            "Historicalinputpaths are provenance; editexplicitly for rebuild. No sourcepathfallback.",
            "Baseline MDstart notforce-converged stationarypoint; existingMDcontactgate passes.",
            "CurrentFF repository/release root has noLICENSE file; no license terms invented.",
        ],
        missing_optional_sources=MISSING_OPTIONAL,
        source_data_modified=False,
        new_MD_or_QM_run=False,
    )
    write(OUT / "data_manifest.json", report)
    logging.basicConfig(level=logging.INFO)
    logging.info("Inventory summary: %s", json.dumps(report["summary"]))


if __name__ == "__main__":
    main()
