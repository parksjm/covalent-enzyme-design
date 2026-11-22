#!/usr/bin/env bash
# Run on the laptop. Download the curated inspection and VMD packages.
set -euo pipefail

usage() {
    echo "Usage: $0 [--dry-run] [USER@HOST [DESTINATION]]"
    echo "By default, pull from $(id -un)@violet into ./nyl12-inspection-20261021"
}

dry_run=0
if [[ ${1:-} == --dry-run ]]; then
    dry_run=1
    shift
fi

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    usage
    exit 0
fi
if [[ $# -gt 2 ]]; then
    usage >&2
    exit 2
fi
remote=${1:-$(id -un)@violet}
destination=${2:-./nyl12-inspection-20261021}
# Keep the SSH/rsync remote specification unambiguous; SSH aliases are supported.
if [[ ! $remote =~ ^[a-zA-Z0-9_][a-zA-Z0-9_.@-]*$ ]]; then
    echo "Use an SSH alias or USER@HOST (no spaces, colons, or options)." >&2
    exit 2
fi
if [[ -z $destination || $destination == -* ]]; then
    echo "Destination must be a directory path, not an option." >&2
    exit 2
fi
command -v rsync >/dev/null
command -v ssh >/dev/null
if command -v shasum >/dev/null; then
    checksum=(shasum -a 256 -c)
elif command -v sha256sum >/dev/null; then
    checksum=(sha256sum -c)
else
    echo "Install a SHA-256 checksum utility before downloading." >&2
    exit 2
fi

verify_nested_manifests() {
    local root=$1
    local manifest
    local found=0
    while IFS= read -r -d '' manifest; do
        found=1
        (cd "$(dirname "$manifest")" && "${checksum[@]}" "$(basename "$manifest")")
    done < <(find "$root" -name MANIFEST.sha256 -type f -print0)
    if [[ $found -ne 1 ]]; then
        echo "No checksum manifests found under $root" >&2
        return 1
    fi
}

source_dir=/data/jpc/covalent-enzyme-design/output_laptop_inspection_20261021
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$destination"
fi
# No --delete: files already on the laptop are never deleted.
# These exclusions are defensive; the workstation package is already curated.
rsync_args=(-rtpvz --safe-links --max-size=10m)
if [[ $dry_run -eq 1 ]]; then
    rsync_args+=(--dry-run)
fi
rsync "${rsync_args[@]}" \
    --exclude='*.chk' --exclude='*.fchk' --exclude='*.log' \
    --exclude='*.prmtop' --exclude='*.inpcrd' --exclude='*.rst7' \
    --exclude='*.nc' --exclude='*.dcd' --exclude='scratch/' \
    "$remote:$source_dir/" "$destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$destination" && "${checksum[@]}" MANIFEST.sha256)
fi

# This separately curated package contains only selected VMD inspection frames.
# Its DCD files are intentional; exclusions above apply to the older package.
vmd_source=/data/jpc/covalent-enzyme-design/output_vmd_inspection_20261024
vmd_destination="$destination/vmd_initial_final"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$vmd_destination"
fi
vmd_rsync_args=(-rtpvz --safe-links --max-size=100m)
if [[ $dry_run -eq 1 ]]; then
    vmd_rsync_args+=(--dry-run)
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.gz' --include='*.dcd' --include='*.vmd' \
    --include='*.tcl' --include='*.json' --include='*.md' \
    --include='*.png' --include='*.csv' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$vmd_source/" "$vmd_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$vmd_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Download verified. Read $vmd_destination/README.md for VMD instructions."
fi

# Retrospective control calibration: small tables, plots and reproducible sources.
calibration_source=/data/jpc/covalent-enzyme-design/analysis/nyl12_filter_calibration
calibration_destination="$destination/analysis/nyl12_filter_calibration"
if ssh "$remote" test -f "$calibration_source/report.json"; then
    if [[ $dry_run -eq 0 ]]; then mkdir -p "$calibration_destination"; fi
    rsync "${vmd_rsync_args[@]}" \
        --exclude='**/__pycache__/***' --exclude='**/.pytest_cache/***' \
        --exclude='/plots/initial_annotations/***' \
        --include='*/' --include='*.json' --include='*.csv' --include='*.md' \
        --include='*.png' --include='*.svg' --include='*.py' --include='*.sh' \
        --include='*.xlsx' --include='*.txt' \
        --include='*.toml' --include='*.lock' --include='*.log' --exclude='*' \
        "$remote:$calibration_source/" "$calibration_destination/"
    if [[ $dry_run -eq 0 ]]; then
        echo "WT/mutant calibration: $calibration_destination/README.md"
        echo "Calibration plots: $calibration_destination/plots/"
    fi
fi

# Exploratory sequence-design pilot: preserve the workstation directory layout.
pilot_source=/data/jpc/covalent-enzyme-design/design/nyl12_ligandmpnn_pilot
pilot_destination="$destination/design/nyl12_ligandmpnn_pilot"
if ssh "$remote" test -f "$pilot_source/accepted_parent/preparation_manifest.json"; then
    if [[ $dry_run -eq 0 ]]; then mkdir -p "$pilot_destination"; fi
    rsync "${vmd_rsync_args[@]}" \
        --exclude='**/__pycache__/***' --exclude='*.pt' \
        --exclude='/backbone_audit/***' --exclude='/prepared/***' \
        --exclude='/preparation_job/***' --exclude='/inputs/relaxed_pdb/***' \
        --exclude='/native_preflight/canonical_pool/***' \
        --exclude='/boltz_validation/source_snapshot/boltz/***' \
        --exclude='/boltz_validation/predictions/runs/***' \
        --exclude='/boltz_validation/predictions/WT/seed_01/***' \
        --exclude='/boltz_validation/checks/***' \
        --exclude='/boltz_validation/models/*/seed_*/relaxation/system/***' \
        --exclude='/boltz_validation/models/*/seed_*/relaxation/inputs/***' \
        --exclude='/boltz_validation/models/*/seed_*/relaxation/minimization/***' \
        --exclude='/boltz_validation/models/*/seed_*/preparation_history/***' \
        --include='*/' --include='*.json' --include='*.csv' --include='*.md' \
        --include='*.fasta' --include='*.fa' --include='*.pdb' --include='*.cif' \
        --include='*.py' --include='*.sh' --include='*.log' \
        --include='*.yaml' --include='*.txt' --include='*.jsonl' \
        --include='*.mol2' --include='*.lib' --include='*.frcmod' --exclude='*' \
        "$remote:$pilot_source/" "$pilot_destination/"
    if [[ $dry_run -eq 0 ]]; then
        echo "LigandMPNN pilot: $pilot_destination/README.md"
    fi
fi

# Current controlling plans, copied directly rather than from an older bundle.
planning_source=/data/jpc/covalent-enzyme-design/docs
planning_destination="$destination/docs"
if [[ $dry_run -eq 0 ]]; then mkdir -p "$planning_destination"; fi
rsync "${rsync_args[@]}" \
    --include='MASTER_PLAN_NYL12_COVALENT_REDESIGN.md' \
    --include='NYL12_MECHANISM_ROADMAP.md' \
    --include='plm_and_filter_calibration_20261121.md' \
    --include='nyl12_ligandmpnn_pilot.md' --include='CURRENT_JOBS.md' \
    --include='standalone_handoff_20261020.md' \
    --exclude='*' "$remote:$planning_source/" "$planning_destination/"

# Current QM/MM work: a small curated tree, organized by stage and window.
qmmm_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_asp272
qmmm_destination="$destination/qmmm/pa6_asp272"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$qmmm_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --exclude='/01_coordinate_drive_to_ti1/live_inspection/' \
    --exclude='/01_coordinate_drive_to_ti1/live_inspection_complete/' \
    --include='*/' --include='*.gz' --include='*.dcd' --include='*.vmd' \
    --include='*.tcl' --include='*.json' --include='*.md' --include='*.png' \
    --include='*.csv' --include='*.py' --include='*.mdin' --include='*.RST' \
    --include='*.cpptraj' --include='*.log' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$qmmm_source/" "$qmmm_destination/"
if [[ $dry_run -eq 0 ]]; then
    verify_nested_manifests "$qmmm_destination"
    echo "QM/MM: $qmmm_destination/plots/distances.png"
    echo "VMD: $qmmm_destination/vmd/view_path.vmd"
    echo "Active E-S optimization: $qmmm_destination/00_es_unrestrained/vmd/view_optimization.vmd"
    echo "Unrestrained TI1 validation: $qmmm_destination/02_ti1_unrestrained_validation/vmd_optimizer_trajectory/view_stage2.vmd"
    echo "Run scripts: $qmmm_destination/scripts/"
fi

# Completed frame-1603 reaction pathways. Each latest immutable bundle contains
# only stripped coordinates and its matching topology.
reaction_paths_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths
reaction_paths_destination="$destination/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths"
for path_name in 01_proton_transfer 02_asynchronous_es_to_ti1; do
    vmd_source="$reaction_paths_source/${path_name}/vmd_progress/latest"
    vmd_destination="$reaction_paths_destination/${path_name}/vmd_progress"
    if [[ $dry_run -eq 0 ]]; then
        mkdir -p "$vmd_destination"
    fi
    rsync "${vmd_rsync_args[@]}" \
        --include='*.parm7' --include='*.dcd' --include='*.vmd' \
        --include='*.json' --include='*.cpptraj' --include='*_export.log' \
        --include='MANIFEST.sha256' \
        --exclude='*' \
        "$remote:$vmd_source/" "$vmd_destination/"
    if [[ $dry_run -eq 0 ]]; then
        (cd "$vmd_destination" && "${checksum[@]}" MANIFEST.sha256)
        echo "Reaction-path VMD: $vmd_destination/reaction_path.vmd"
    fi
done

# Small potential-energy profiles for the direct proton-transfer and
# asynchronous E-S-to-TI1 coordinate drives.
for plot_name in 01_proton_transfer 02_asynchronous_es_to_ti1; do
    if [[ $plot_name == 01_proton_transfer ]]; then
        plot_source="$reaction_paths_source/01_proton_transfer/plots"
    else
        plot_source="$reaction_paths_source/02_asynchronous_es_to_ti1/plots"
    fi
    plot_destination="$reaction_paths_destination/${plot_name}/plots"
    if [[ $dry_run -eq 0 ]]; then
        mkdir -p "$plot_destination"
    fi
    rsync "${vmd_rsync_args[@]}" \
        --include='*.png' --include='*.json' --include='*.csv' \
        --include='MANIFEST.sha256' --exclude='*' \
        "$remote:$plot_source/" "$plot_destination/"
    if [[ $dry_run -eq 0 ]]; then
        (cd "$plot_destination" && "${checksum[@]}" MANIFEST.sha256)
        echo "Potential energy: $plot_destination/potential_energy.png"
    fi
done

# Direct proton-transfer umbrella diagnostics and the sequential, downsampled
# production-window movie. Preserve the workstation's campaign-relative layout.
umbrella_campaign_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/04_direct_pt_umbrella_initial
umbrella_campaign_destination="$destination/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/04_direct_pt_umbrella_initial"
overlap_source="$umbrella_campaign_source/analysis/03_production_2ps_overlap"
overlap_destination="$umbrella_campaign_destination/analysis/03_production_2ps_overlap"
mbar_source="$umbrella_campaign_source/analysis/04_mbar_pipeline_smoke_2ps"
mbar_destination="$umbrella_campaign_destination/analysis/04_mbar_pipeline_smoke_2ps"
direct_extension_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/07_direct_pt_umbrella_30ps_violet
direct_extension_destination="$destination/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/07_direct_pt_umbrella_30ps_violet"
mbar_7ps_source="$direct_extension_source/analysis/mbar_7ps_complete/result"
mbar_7ps_destination="$direct_extension_destination/analysis/mbar_7ps_complete/result"
umbrella_vmd_source="$umbrella_campaign_source/visualization/sequential_production_windows"
umbrella_vmd_destination="$umbrella_campaign_destination/visualization/sequential_production_windows"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$overlap_destination" "$mbar_destination" \
        "$mbar_7ps_destination" \
        "$umbrella_vmd_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.png' --include='*.json' --include='*.csv' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$overlap_source/" "$overlap_destination/"
rsync "${vmd_rsync_args[@]}" \
    --include='*.png' --include='*.json' --include='*.csv' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$mbar_source/" "$mbar_destination/"
rsync "${vmd_rsync_args[@]}" \
    --include='*.png' --include='*.json' --include='*.csv' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$mbar_7ps_source/" "$mbar_7ps_destination/"
umbrella_vmd_rsync_args=(-rtpvz --safe-links --max-size=500m)
if [[ $dry_run -eq 1 ]]; then
    umbrella_vmd_rsync_args+=(--dry-run)
fi
rsync "${umbrella_vmd_rsync_args[@]}" \
    --include='*.parm7' --include='*.dcd' --include='*.vmd' \
    --include='*.json' --include='*.cpptraj' --include='*_export.log' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$umbrella_vmd_source/" "$umbrella_vmd_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$overlap_destination" && "${checksum[@]}" MANIFEST.sha256)
    (cd "$mbar_destination" && "${checksum[@]}" MANIFEST.sha256)
    (cd "$mbar_7ps_destination" && "${checksum[@]}" MANIFEST.sha256)
    (cd "$umbrella_vmd_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Umbrella overlap: $overlap_destination/window_overlap.png"
    echo "Direct MBAR smoke test: $mbar_destination/provisional_pmf.png"
    echo "Direct interim 7-ps MBAR: $mbar_7ps_destination/pmf.png"
    echo "Umbrella VMD: $umbrella_vmd_destination/view_umbrella.vmd"
fi

# Immutable partial-block MBAR pipeline test for the asynchronous campaign.
async_mbar_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/06_asynchronous_umbrella_30ps_azalea/analysis/pipeline_pmf_block01_partial_20261109T2109EST/result
async_mbar_destination="$destination/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/06_asynchronous_umbrella_30ps_azalea/analysis/pipeline_pmf_block01_partial_20261109T2109EST/result"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$async_mbar_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.png' --include='*.json' --include='*.csv' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$async_mbar_source/" "$async_mbar_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$async_mbar_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Asynchronous MBAR smoke test: $async_mbar_destination/provisional_pmf.png"
fi

# Window-24 proton-identity and reactive-geometry audit underlying the PMF
# interpretation. Trajectory snapshots remain on violet; only small tables and
# the machine-readable summary are transferred.
proton_audit_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/06_asynchronous_umbrella_30ps_azalea/analysis/window24_proton_state_20261109T2123EST
proton_audit_destination="$destination/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/06_asynchronous_umbrella_30ps_azalea/analysis/window24_proton_state_20261109T2123EST"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$proton_audit_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.dat' --include='*.json' --include='MANIFEST.sha256' \
    --exclude='*' \
    "$remote:$proton_audit_source/" "$proton_audit_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$proton_audit_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Window-24 proton audit: $proton_audit_destination/proton_state_summary.json"
fi

# Small static E-S seed screen: canonical Tier A model and catalytic overlays.
static_es_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_static_es_seed
static_es_destination="$destination/qmmm/pa6_static_es_seed"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$static_es_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*/' --include='*.pdb' --include='*.json' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$static_es_source/" "$static_es_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$static_es_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Static E-S seed: $static_es_destination/pa6_l4z_l2_nyl12_wt_model_23_canonical.pdb"
    echo "Catalytic overlays: $static_es_destination/model_{23,40}_alignment/catalytic_site_overlay.pdb"
fi

# Small mechanistic literature and proton-donor reach audit; retain layout.
mechanism_evidence_source=/data/jpc/covalent-enzyme-design/analysis/nyl12_mechanism_evidence_20261118
mechanism_evidence_destination="$destination/analysis/nyl12_mechanism_evidence_20261118"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$mechanism_evidence_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.json' --include='*.md' --include='*.py' --include='*.csv' \
    --exclude='*' \
    "$remote:$mechanism_evidence_source/" "$mechanism_evidence_destination/"

# Bounded model-30 Thr activation/attack grid. Transfer protocol and named
# progress products only; full-system restarts and QM scratch stay on violet.
activation_grid_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid
activation_grid_destination="$destination/qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$activation_grid_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='/scripts/' --include='/scripts/*.py' \
    --include='/prepare_grid.py' --include='/run_row.py' \
    --include='/reference_checks/' --include='/reference_checks/*.py' \
    --include='/reference_checks/queue*.json' --include='/reference_checks/README.md' \
    --include='/reference_checks/results.json' --include='/reference_checks/selection.json' \
    --include='/reference_checks/gradient_analysis/' \
    --include='/reference_checks/gradient_analysis/acceptance.json' \
    --include='/reference_checks/gradient_analysis/results.json' \
    --include='/reference_checks/gradient_analysis/projections.csv' \
    --include='/reference_checks/gradient_analysis/directional_derivatives.png' \
    --include='/reference_checks/gradient_analysis/analyze.py' \
    --include='/reference_checks/gradient_analysis/plot.py' \
    --include='/reference_checks/gradient_analysis/verify_embedding.py' \
    --include='/reference_checks/gradient_analysis/test_analysis.py' \
    --include='/occupancy_callback_pilot/' \
    --include='/occupancy_callback_pilot/run_pilot.py' \
    --include='/occupancy_callback_pilot/resume_pilot.py' \
    --include='/occupancy_callback_pilot/runtime_support.py' \
    --include='/occupancy_callback_pilot/verify_callback.py' \
    --include='/occupancy_callback_pilot/monitor_progress.py' \
    --include='/occupancy_callback_pilot/export_support.py' \
    --include='/occupancy_callback_pilot/README.md' \
    --include='/occupancy_callback_pilot/release_endpoint.py' \
    --include='/occupancy_callback_pilot/monitor_release.py' \
    --include='/occupancy_callback_pilot/endpoint_audit/' \
    --include='/occupancy_callback_pilot/endpoint_audit/*.py' \
    --include='/occupancy_callback_pilot/endpoint_audit/*.json' \
    --include='/occupancy_callback_pilot/endpoint_audit/*.csv' \
    --include='/occupancy_callback_pilot/d_reference/' \
    --include='/occupancy_callback_pilot/d_reference/*.py' \
    --include='/occupancy_callback_pilot/d_reference/*.json' \
    --include='/occupancy_callback_pilot/d_reference/*.csv' \
    --include='/occupancy_callback_pilot/d_reference/gaussian_direct/' \
    --include='/occupancy_callback_pilot/d_reference/gaussian_direct/process.json' \
    --include='/occupancy_callback_pilot/d_reference/gfn2_xtb/' \
    --include='/occupancy_callback_pilot/d_reference/gfn2_xtb/process.json' \
    --include='/occupancy_callback_pilot/d_reference/gfn2_xtb/qm_link.xyz' \
    --include='/occupancy_callback_pilot/job_D_reference/' \
    --include='/occupancy_callback_pilot/job_D_release/' \
    --include='/occupancy_callback_pilot/job_D_release_progress/' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/*.json' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/accepted_iterations.jsonl' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/evaluate.mdin' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/accepted_iterations.csv' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/accepted_progress.png' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/accepted_progress.svg' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/progress.json' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/snapshot_stability.json' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/snapshot_energies.csv' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/snapshot_coordinates.csv' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/exported_frames.txt' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/exporter_source/' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/exporter_source/*.py' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/exporter_source/*.json' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/path.parm7' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/path.dcd' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/path.csv' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/energy.svg' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/distances.svg' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/angles.svg' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/view_path.vmd' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/manifest.json' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/MANIFEST.sha256' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/export.cpptraj' \
    --include='/occupancy_callback_pilot/pilot/point_D/unbiased_release/progress/vmd/cpptraj.log' \
    --include='/occupancy_callback_pilot/line_search_diagnosis/' \
    --include='/occupancy_callback_pilot/line_search_diagnosis/*.py' \
    --include='/occupancy_callback_pilot/line_search_diagnosis/*.json' \
    --include='/occupancy_callback_pilot/line_search_diagnosis/*.csv' \
    --include='/occupancy_callback_pilot/callback_proof/' \
    --include='/occupancy_callback_pilot/callback_proof/protocol.json' \
    --include='/occupancy_callback_pilot/callback_proof/evaluations.jsonl' \
    --include='/occupancy_callback_pilot/callback_proof/scientific_acceptance.json' \
    --include='/occupancy_callback_pilot/callback_proof/library_coverage.json' \
    --include='/occupancy_callback_pilot/callback_proof/verification.json' \
    --include='/occupancy_callback_pilot/job_pilot/' \
    --include='/occupancy_callback_pilot/job_callback_proof/' \
    --include='/occupancy_callback_pilot/job_progress/' \
    --include='/occupancy_callback_pilot/job_pilot_continuation/' \
    --include='/occupancy_callback_pilot/job_progress_continuation/' \
    --include='/occupancy_callback_pilot/job_*/command.json' \
    --include='/occupancy_callback_pilot/job_*/status.json' \
    --include='/occupancy_callback_pilot/job_*/run.sh' \
    --include='/occupancy_callback_pilot/pilot/' \
    --include='/occupancy_callback_pilot/pilot/protocol.json' \
    --include='/occupancy_callback_pilot/pilot/preflight.json' \
    --include='/occupancy_callback_pilot/pilot/completion.json' \
    --include='/occupancy_callback_pilot/pilot/continuation_001/' \
    --include='/occupancy_callback_pilot/pilot/continuation_001/protocol.json' \
    --include='/occupancy_callback_pilot/pilot/history/' \
    --include='/occupancy_callback_pilot/pilot/history/line_search_stop/' \
    --include='/occupancy_callback_pilot/pilot/history/line_search_stop/completion.json' \
    --include='/occupancy_callback_pilot/pilot/point_A/history/' \
    --include='/occupancy_callback_pilot/pilot/point_A/history/line_search_stop/' \
    --include='/occupancy_callback_pilot/pilot/point_A/history/line_search_stop/process.json' \
    --include='/occupancy_callback_pilot/pilot/point_A/history/line_search_stop/scientific_acceptance.json' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/stage.json' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/process.json' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/scientific_acceptance.json' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/process_failure.json' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/scientific_failure.json' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/accepted_iterations.jsonl' \
    --include='/occupancy_callback_pilot/pilot/point_[ABCD]/continuation_optimizer_events.jsonl' \
    --include='/occupancy_callback_pilot/pilot/progress/' \
    --include='/occupancy_callback_pilot/pilot/progress/accepted_iterations.csv' \
    --include='/occupancy_callback_pilot/pilot/progress/accepted_progress.png' \
    --include='/occupancy_callback_pilot/pilot/progress/accepted_progress.svg' \
    --include='/occupancy_callback_pilot/pilot/progress/progress.json' \
    --include='/occupancy_callback_pilot/pilot/progress/snapshot_stability.json' \
    --include='/occupancy_callback_pilot/pilot/progress/snapshot_energies.csv' \
    --include='/occupancy_callback_pilot/pilot/progress/snapshot_coordinates.csv' \
    --include='/occupancy_callback_pilot/pilot/progress/exported_frames.txt' \
    --include='/occupancy_callback_pilot/pilot/progress/exporter_source/' \
    --include='/occupancy_callback_pilot/pilot/progress/exporter_source/manifest.json' \
    --include='/occupancy_callback_pilot/pilot/progress/exporter_source/monitor_progress.py' \
    --include='/occupancy_callback_pilot/pilot/progress/exporter_source/export_support.py' \
    --include='/occupancy_callback_pilot/pilot/progress/exporter_source/export_water_path.py' \
    --include='/occupancy_callback_pilot/pilot/progress/exporter_source/water_path_vmd_support.py' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/path.parm7' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/path.dcd' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/path.csv' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/energy.svg' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/distances.svg' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/angles.svg' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/view_path.vmd' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/manifest.json' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/MANIFEST.sha256' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/export.cpptraj' \
    --include='/occupancy_callback_pilot/pilot/progress/vmd/cpptraj.log' \
    --include='/proton_occupancy_audit/' \
    --include='/proton_occupancy_audit/*.json' --include='/proton_occupancy_audit/*.csv' \
    --include='/proton_occupancy_audit/*.png' --include='/proton_occupancy_audit/*.py' \
    --include='/occupancy_coordinate_preflight/' \
    --include='/occupancy_coordinate_preflight/coordinate_definition.json' \
    --include='/occupancy_coordinate_preflight/analytic_implementation.json' \
    --include='/occupancy_coordinate_preflight/attempt_2/' \
    --include='/occupancy_coordinate_preflight/attempt_2/scientific_acceptance_corrected.json' \
    --include='/diagnostic_activation_release/' \
    --include='/diagnostic_activation_release/*.json' \
    --include='/diagnostic_activation_release/*.py' \
    --include='/diagnostic_activation_release/unrestrained.mdin' \
    --include='/diagnostic_activation_release/vmd/' \
    --include='/diagnostic_activation_release/vmd/*' \
    --include='/release_queue_status.json' --include='/queue_candidate_releases.py' \
    --include='/release_queue_support.py' --include='/release_execution.py' \
    --include='/protocol.json' --include='/grid_points.json' \
    --include='/grid_targets.csv' --include='/progress/' \
    --include='/progress/**/' --include='/progress/**/grid_heatmap.png' \
    --include='/progress/**/grid_points.csv' \
    --include='/progress/**/view_path.vmd' \
    --include='/progress/**/path.parm7' --include='/progress/**/path.dcd' \
    --include='/progress/**/energy.svg' --include='/progress/**/distances.svg' \
    --include='/progress/**/angles.svg' --include='/progress/**/path.csv' \
    --include='/progress/**/export.cpptraj' --include='/progress/**/cpptraj.log' \
    --include='/progress/**/manifest.json' \
    --include='/progress/**/MANIFEST.sha256' --exclude='*' \
    "$remote:$activation_grid_source/" "$activation_grid_destination/"
if [[ $dry_run -eq 0 ]]; then
    echo "Activation/attack grid protocol: $activation_grid_destination/protocol.json"
    echo "Activation/attack grid progress: $activation_grid_destination/progress/"
fi

# Relay acylation campaign: compact reports and latest accepted-state movie only.
relay_campaign_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation
relay_campaign_destination="$destination/qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$relay_campaign_destination"
fi
rsync "${vmd_rsync_args[@]}" --copy-links \
    --exclude='*.ncrst' --exclude='*.rst7' --exclude='*.nc' \
    --exclude='evaluations/' --exclude='frame_links/' --exclude='.building_*/' \
    --exclude='accepted_iterations.jsonl' --exclude='accepted_iterations.json' \
    --exclude='native_evidence.jsonl' --exclude='preparation_history/' \
    --exclude='construction_diagnostic*' \
    --exclude='*.so' --exclude='*.o' --exclude='*.exe' \
    --exclude='/progress/update_*/' \
    --exclude='/04_neb/progress/update_*/' \
    --exclude='/04_neb/progress/source_milestones/' \
    --exclude='/04_neb/history/' \
    --exclude='/04_neb/spacing_tests/**/progress/update_*/' \
    --exclude='/04_neb/spacing_tests/**/source_milestones/' \
    --exclude='/04_neb/spacing_tests/**/source_snapshot*/' \
    --exclude='/04_neb/spacing_tests/**/inputs/' \
    --exclude='/04_neb/spacing_tests/**/history/' \
    --include='/04_neb/spacing_tests/**/calculation/' \
    --include='/04_neb/spacing_tests/**/calculation/minimize.mdin' \
    --exclude='/04_neb/spacing_tests/**/calculation/***' \
    --include='/04_neb/calculation/' \
    --include='/04_neb/calculation/minimize.mdin' \
    --exclude='/04_neb/calculation/***' --exclude='*.xyz' --exclude='*.xyz.gz' \
    --exclude='*.mdout' --exclude='native_report.mdout' \
    --include='/04_neb/spacing_tests/**/progress/' \
    --include='/04_neb/spacing_tests/**/progress/latest/' \
    --include='/04_neb/spacing_tests/**/progress/latest/***' \
    --include='/04_neb/progress/' --include='/04_neb/progress/latest/' \
    --include='/04_neb/progress/latest/***' \
    --include='/04_neb/progress/progress.json' \
    --include='/04_neb/progress/monitor_status.json' \
    --include='/progress/' --include='/progress/latest/' \
    --include='/progress/latest/***' --include='/progress/progress.json' \
    --include='/progress/snapshot_stability.json' \
    --include='*/' --include='*.json' --include='*.jsonl' \
    --include='*.csv' --include='*.md' --include='*.py' --include='*.sh' \
    --include='*.mdin' --include='*.png' --include='*.svg' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$relay_campaign_source/" "$relay_campaign_destination/"
if [[ $dry_run -eq 0 ]]; then
    echo "Relay campaign accepted progress: $relay_campaign_destination/progress/latest/"
fi

# Final embedded-method comparison and DFT-corrected periodic-path diagnostic.
method_benchmark_source=/data/jpc/covalent-enzyme-design/qmmm/benchmarks/03_direct_pt_embedded_method_comparison/final_analysis
method_benchmark_destination="$destination/qmmm/benchmarks/03_direct_pt_embedded_method_comparison/final_analysis"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$method_benchmark_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.png' --include='*.json' --include='*.csv' \
    --include='*.md' --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$method_benchmark_source/" "$method_benchmark_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$method_benchmark_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Embedded-method benchmark: $method_benchmark_destination/dft_corrected_periodic_profile.png"
fi

# Download the acylenzyme/NEB artifacts already mirrored on violet.
# Workstation-to-workstation refreshes are separate from laptop downloads.
acyl_source="$reaction_paths_source/03_acylenzyme_neb"
acyl_destination="$reaction_paths_destination/03_acylenzyme_neb"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$acyl_destination" "$destination/docs"
fi
rsync "${vmd_rsync_args[@]}" --exclude='vmd_progress/' --exclude='clutter/' --include='*/' --include='latest' \
    --include='*.json' --include='*.csv' --include='*.svg' --include='*.png' \
    --include='*.dcd' --include='*.parm7' --include='*.vmd' \
    --include='*.cpptraj' --include='*_export.log' --include='MANIFEST.sha256' \
    --include='*.py' --include='*.sh' --include='*.mdin' --include='*.RST' \
    --exclude='*' "$remote:$acyl_source/" "$acyl_destination/"
rsync "${rsync_args[@]}" \
    "$remote:/data/jpc/covalent-enzyme-design/docs/acylenzyme_backward_path_20261110.md" \
    "$destination/docs/"
rsync "${rsync_args[@]}" \
    "$remote:/data/jpc/covalent-enzyme-design/docs/nyl12_active_site_conservation_20261110.md" \
    "$destination/docs/"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$destination/analysis/nyl12_active_site_conservation_20261110"
fi
rsync "${rsync_args[@]}" \
    "$remote:/data/jpc/covalent-enzyme-design/analysis/nyl12_active_site_conservation_20261110/" \
    "$destination/analysis/nyl12_active_site_conservation_20261110/"

# Download only the latest acylenzyme band, not every archived iteration.
if ssh "$remote" test -d "$acyl_source/vmd_progress/latest"; then
    if [[ $dry_run -eq 0 ]]; then mkdir -p "$acyl_destination/vmd_progress/latest"; fi
    rsync "${vmd_rsync_args[@]}" --exclude='source_restarts/' \
        "$remote:$acyl_source/vmd_progress/latest/" "$acyl_destination/vmd_progress/latest/"
    if [[ $dry_run -eq 0 ]]; then
        (cd "$acyl_destination/vmd_progress/latest" && "${checksum[@]}" MANIFEST.sha256)
    fi
else
    echo "Acylenzyme NEB VMD band is not available yet."
fi

# Small reproducible flipped-pose screening tables and mechanism assessment.
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$destination/analysis/pa6_binding_orientation_20261110" "$destination/docs"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.csv' --include='*.json' --exclude='*' \
    "$remote:/data/jpc/covalent-enzyme-design/analysis/pa6_binding_orientation_20261110/" \
    "$destination/analysis/pa6_binding_orientation_20261110/"
rsync "${vmd_rsync_args[@]}" \
    "$remote:/data/jpc/covalent-enzyme-design/docs/pa6_flipped_pose_screen_20261110.md" \
    "$remote:/data/jpc/covalent-enzyme-design/docs/nyl12_alternative_proton_routes_20261110.md" \
    "$destination/docs/"

# Intentional stop record for the retained direct-proton-transfer umbrella data.
rsync "${vmd_rsync_args[@]}" \
    "$remote:/data/jpc/covalent-enzyme-design/docs/direct_pt_umbrella_stop_20261110.md" \
    "$destination/docs/"

# Flipped model30: compact coordinates, run inputs and the complete latest VMD
# bundle. Materialize latest symlinks; skip immutable historical movie copies.
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$destination/qmmm/pa6_flipped_model30"
    echo "TS search VMD: $destination/qmmm/pa6_flipped_model30/05_neutral_reaction_path/06_dimer_ts/vmd_progress/latest/view_ts_optimization.vmd"
    echo "Converged TS VMD: $destination/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/vmd_progress/latest/view_ts_optimization.vmd"
fi
rsync "${vmd_rsync_args[@]}" --copy-links \
    --exclude='source_snapshot/' --exclude='clutter/' \
    --exclude='/11_leaving_group_relay_screen/coarse/' \
    --exclude='/11_leaving_group_relay_screen/coarse_nonflipped/' \
    --exclude='/11_leaving_group_relay_screen/refined/' \
    --exclude='/11_leaving_group_relay_screen/assisted_refined/' \
    --exclude='/11_leaving_group_relay_screen/pilot/' \
    --exclude='/11_leaving_group_relay_screen/selection_revision_02/' \
    --exclude='/11_leaving_group_relay_screen/source_snapshot_nonflipped/' \
    --exclude='/11_leaving_group_relay_screen/assisted_source_snapshot/' \
    --exclude='/11_leaving_group_relay_screen/assisted_review*/ranked_assisted_candidates.json' \
    --exclude='/10_wt_protonation_md/analysis_after_md/occupancy/*.csv' \
    --exclude='/10_wt_protonation_md/analysis_after_md/relay/new_four/*_paths.csv' \
    --include='/10_wt_protonation_md/analysis_after_md/logs/*.log' \
    --include='/10_wt_protonation_md/analysis_after_md/job/job.log' \
    --include='/10_wt_protonation_md/analysis_after_md/run_analysis.sh' \
    --exclude='/07_ash274_comparison/09_water_relay_acylation/' \
    --exclude='/07_ash274_comparison/07_water_to_asp272_ti1_like/history/' \
    --exclude='vmd_progress/update_*/' --include='*/' \
    --include='endpoint_complex.pdb' --include='raw_model30.pdb' \
    --include='*.json' --include='min.in' --include='measure_endpoint.py' \
    --include='*.mdin' --include='*.RST' --include='run.sh' \
    --include='*.cpptraj' --include='*_export.log' \
    --include='MANIFEST.sha256' \
    --include='*.vmd' --include='*.tcl' --include='*.dcd' \
    --include='stripped.parm7' --include='stripped.parm7.gz' \
    --include='trajectory_export.log' --include='topology_export.log' \
    --include='*.png' --include='*.csv' --include='*.md' \
    --exclude='*' \
    "$remote:/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/" \
    "$destination/qmmm/pa6_flipped_model30/"

# One compact 99-frame endpoint-minimization movie plus its small acceptance
# plot and analysis table. Preserve the workstation campaign-relative layout.
downhill_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K
downhill_destination="$destination/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K"
downhill_vmd_destination="$downhill_destination/vmd_progress/latest"
downhill_plots_destination="$downhill_destination/plots"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$downhill_vmd_destination" "$downhill_plots_destination"
fi
rsync "${vmd_rsync_args[@]}" --copy-links \
    --include='stripped.parm7' --include='trajectory.dcd' \
    --include='view_downhill_minimizations.vmd' --include='README.md' \
    --include='frame_manifest.json' --include='reactive_distances.csv' \
    --include='trajectory_export.log' --include='export_downhill_minimization_vmd.py' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$downhill_source/vmd_progress/latest/" "$downhill_vmd_destination/"
rsync "${vmd_rsync_args[@]}" \
    --include='endpoint_minimization_progress.png' \
    --include='endpoint_minimization_analysis.json' --exclude='*' \
    "$remote:$downhill_source/plots/" "$downhill_plots_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$downhill_vmd_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Endpoint minimization VMD: $downhill_vmd_destination/view_downhill_minimizations.vmd"
    echo "Endpoint minimization analysis: $downhill_plots_destination/endpoint_minimization_analysis.json"
    echo "Endpoint minimization plot: $downhill_plots_destination/endpoint_minimization_progress.png"
fi

# Current in-progress regular NEB band: a stable native history-matched frame.
# Keep its full-system source restarts so the nested export checksum verifies.
saved_neb_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/vmd_progress/saved_history_frame_000104
saved_neb_destination="$destination/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/vmd_progress/saved_history_frame_000104"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$saved_neb_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='source_restarts/' --include='source_restarts/*.ncrst' \
    --include='*.parm7' --include='*.dcd' --include='*.vmd' \
    --include='*.json' --include='*.cpptraj' --include='*_export.log' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$saved_neb_source/" "$saved_neb_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$saved_neb_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Current native NEB VMD (saved-history frame; cycle unverified): $saved_neb_destination/reaction_path.vmd"
fi

# Latest complete report from native NEB mdout. This energy-only profile is
# provisional and intentionally has no saved-coordinate-frame association.
neb_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/plots
neb_destination="$destination/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/plots"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$neb_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.py' --include='*.png' --include='*.svg' \
    --include='*.json' --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$neb_source/" "$neb_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$neb_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Latest provisional NEB energy profile: $neb_destination/latest_neb_energy_profile.png"
fi

ash_si_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/06_ash272_comparison/04_md_validation_20261113/results_attempt1/si_matched_panels
ash_si_destination="$destination/qmmm/pa6_flipped_model30/06_ash272_comparison/04_md_validation_20261113/results_attempt1/si_matched_panels"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$ash_si_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.png' --include='*.pdf' --exclude='*' \
    "$remote:$ash_si_source/" "$ash_si_destination/"

for ash_case in replica_1_tier_A replica_1_tier_B replica_2_tier_A replica_2_tier_B; do
    ash_vmd_source="/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/06_ash272_comparison/06_qmmm_optimizations/$ash_case/vmd_optimization"
    ash_vmd_destination="$destination/qmmm/pa6_flipped_model30/06_ash272_comparison/06_qmmm_optimizations/$ash_case/vmd_optimization"
    if [[ $dry_run -eq 0 ]]; then
        mkdir -p "$ash_vmd_destination"
    fi
    rsync "${vmd_rsync_args[@]}" \
        --include='*.dcd' --include='*.parm7' --include='*.vmd' \
        --include='*.cpptraj' --include='*.json' --include='*.log' \
        --include='MANIFEST.sha256' --exclude='*' \
        "$remote:$ash_vmd_source/" "$ash_vmd_destination/"
    if [[ $dry_run -eq 0 ]]; then
        (cd "$ash_vmd_destination" && "${checksum[@]}" MANIFEST.sha256)
    fi
done

# Completed, converged 24-image trace-seeded direct-route NEB movie and audits.
trace_neb_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/vmd_final_24img/final_24img
trace_neb_destination="$destination/qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/vmd_final_24img/final_24img"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$trace_neb_destination"
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.parm7' --include='*.dcd' --include='*.vmd' \
    --include='*.json' --include='*.csv' --include='*.cpptraj' \
    --include='*.log' --include='*.md' --include='*.png' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$trace_neb_source/" "$trace_neb_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$trace_neb_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Completed 24-image NEB VMD: $trace_neb_destination/reaction_path.vmd"
    echo "NEB profile: $trace_neb_destination/energy_profile.png"
    echo "Reactive distances: $trace_neb_destination/reactive_distances.csv"
fi

# NEB-seeded two-coordinate US: compact current diagnostics and window views.
trace_us_source="${trace_neb_source%/vmd_final_24img/final_24img}/umbrella_sampling"
trace_us_destination="${trace_neb_destination%/vmd_final_24img/final_24img}/umbrella_sampling"
if ssh "$remote" test -f "$trace_us_source/prepared.json"; then
    if [[ $dry_run -eq 0 ]]; then mkdir -p "$trace_us_destination"; fi
    rsync "${vmd_rsync_args[@]}" --copy-links \
        --exclude='/inputs/' \
        --exclude='/analysis_job/' --exclude='/preparation_history/' \
        --exclude='/vmd_dynamic/work/***' \
        --include='/vmd_dynamic/' --include='/vmd_dynamic/***' \
        --exclude='/native_boundary_probe*/' \
        --exclude='/execution_preflight/' --exclude='**/preparation_history/***' \
        --exclude='**/__pycache__/***' \
        --include='**/*source_snapshot/' --include='**/*source_snapshot/**/' \
        --include='**/*source_snapshot/**/*.py' \
        --include='**/*source_snapshot/source_sha256.json' \
        --exclude='**/*source_snapshot/***' \
        --exclude='**/incoming/***' --exclude='**/job/***' --exclude='**/sync_job/***' \
        --include='/analysis/' --include='/analysis/latest/' \
        --include='/analysis/latest/***' \
        --include='/analysis/pmf_latest/' --include='/analysis/pmf_latest/***' \
        --exclude='/analysis/*/' \
        --include='/window_*/progress/***' \
        --include='*/' --include='*.json' --include='*.csv' \
        --include='*.md' --include='*.mdin' --include='*.RST' \
        --include='*.png' --include='*.svg' --include='*.log' \
        --exclude='*' \
        "$remote:$trace_us_source/" "$trace_us_destination/"
    if [[ $dry_run -eq 0 ]]; then
        echo "NEB-seeded US diagnostics: $trace_us_destination/analysis/latest/"
        echo "NEB-seeded PMF: $trace_us_destination/analysis/pmf_latest/"
        echo "NEB-seeded US dynamic movie: $trace_us_destination/vmd_dynamic/view_umbrella.vmd"
    fi
fi

# Cumulative accepted frames of the ASH274 water-relay TI1-like diagnostic.
water_path_source=/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/vmd_progress
water_path_destination="$destination/qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/vmd_progress"
if ssh "$remote" test -f "$water_path_source/MANIFEST.sha256"; then
    if [[ $dry_run -eq 0 ]]; then mkdir -p "$water_path_destination"; fi
    rsync "${vmd_rsync_args[@]}" \
        --include='*.parm7' --include='*.dcd' --include='*.vmd' \
        --include='*.json' --include='*.csv' --include='*.svg' \
        --include='*.cpptraj' --include='*.log' \
        --include='MANIFEST.sha256' --exclude='*' \
        "$remote:$water_path_source/" "$water_path_destination/"
    if [[ $dry_run -eq 0 ]]; then
        (cd "$water_path_destination" && "${checksum[@]}" MANIFEST.sha256)
        echo "Water-relay VMD: $water_path_destination/view_path.vmd"
        echo "Water-relay energy: $water_path_destination/energy.svg"
    fi
fi
