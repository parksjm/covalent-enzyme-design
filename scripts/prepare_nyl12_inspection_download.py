"""Assemble a small, explicit inspection package without copying raw runs."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
from pathlib import Path

LOGGER = logging.getLogger(__name__)
VIEWER_SOURCE = Path(__file__).with_name("nyl12_inspection_view.pml")
MAX_FILE_BYTES = 10 * 1024 * 1024
MAX_TOTAL_BYTES = 25 * 1024 * 1024
MODEL_ALIASES = {
    "pa6_l4": {
        "source_ids": ["pa6_l2", "pa6_l4"],
        "substrate": "PA6 L4",
        "product": "L2",
    },
    "maxa_maxa": {
        "source_ids": ["poly_xama_maxa", "maxa_maxa"],
        "substrate": "MAXA–MAXA",
        "product": "MAXA",
    },
}


def selected_files(
    repo: Path,
    audit: Path,
    structures: Path,
    trajectories: Path,
    protein_source: Path | None = None,
    action_journal: Path | None = None,
    af3_audits: list[Path] | None = None,
) -> dict[str, Path]:
    """Select the documented small artifacts using explicit file inventories.

    Args:
        repo: Repository root containing retained calculation artifacts.
        audit: Completed optimized-conformer audit JSON.
        structures: Directory containing optimized XYZ and SDF exports.
        trajectories: Directory containing converted multi-frame XYZ files.
        protein_source: Directory with the two model PDBs and source metadata.
            Defaults to the retained legacy inspection snapshot.
        action_journal: Optional frozen journal to package instead of the legacy
            repository journal.
        af3_audits: Explicit result directories for diagnostic AF3 water models.

    Returns:
        Relative download names mapped to original source paths.
    """
    old = repo / "output_nyl12_resume_20261020"
    qm = old / "neutral_thr_nme_wb97xd/calculations"
    protein_source = protein_source or old / "nyl12_pymol_starting_models"
    files = {
        "README.md": repo / "docs/laptop_inspection_20261021.md",
        "qm/neutral_thr_nme/conformer_comparison.json": audit,
        "reports/ligand_confidence.json": old / "confidence_audit_v2/report.json",
        "pull_nyl12_inspection.sh": repo / "scripts/pull_nyl12_inspection.sh",
        "gaussian_log_to_xyz.py": repo / "scripts/gaussian_log_to_xyz.py",
        "provenance/agent_actions_at_export.jsonl": action_journal
        or repo / "provenance/agent_actions_20261021.jsonl",
        "provenance/cache_directory_migration.json": repo
        / "provenance/cache_directory_migration_20261021.json",
    }
    for directory in ("docs", "manuscript"):
        for path in sorted((repo / directory).glob("*.md")):
            files[f"{directory}/{path.name}"] = path
    source_files = {path.name: path for path in protein_source.iterdir()}
    for canonical, metadata in MODEL_ALIASES.items():
        source_ids = metadata["source_ids"]
        candidates = [f"{source_id}_ti1_starting.pdb" for source_id in source_ids]
        matches = [source_files[name] for name in candidates if name in source_files]
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous source models for {canonical} in {protein_source}: "
                f"{[path.name for path in matches]}"
            )
        source = matches[0] if matches else None
        if source is None:
            raise FileNotFoundError(
                f"No source model for {canonical} in {protein_source}; "
                f"expected one of {candidates}"
            )
        files[f"protein/{canonical}_ti1_starting.pdb"] = source

    metadata_suffixes = {
        "_ti1_audit.json": "ti1_audit.json",
        "_ti1_provenance.json": "ti1_provenance.json",
        "_residue_map.json": "residue_map.json",
    }
    for canonical, metadata in MODEL_ALIASES.items():
        for suffix, export_name in metadata_suffixes.items():
            source = next(
                (
                    source_files[f"{source_id}{suffix}"]
                    for source_id in metadata["source_ids"]
                    if f"{source_id}{suffix}" in source_files
                ),
                None,
            )
            if source is not None:
                files[f"protein/{canonical}_{export_name}"] = source
    files["protein/view_models.pml"] = VIEWER_SOURCE

    for canonical, metadata in MODEL_ALIASES.items():
        source_id = metadata["source_ids"][0]
        files[f"reports/{canonical}_relay_feasibility.json"] = (
            repo
            / f"output_nyl12_resume_20261021/relay_feasibility/{source_id}/report.json"
        )
        files[f"reports/{canonical}_active_site.json"] = (
            old / f"active_site_baseline/{source_id}.json"
        )
    for name in ("conformer_01", "conformer_02"):
        for filename in (
            "atom_map.json",
            "acceptance.json",
            "process.json",
            "calculation.json",
            "input.gjf",
        ):
            files[f"qm/neutral_thr_nme/{name}/{filename}"] = qm / name / filename
        for extension in ("xyz", "sdf"):
            files[f"qm/neutral_thr_nme/{name}/optimized.{extension}"] = (
                structures / f"{name}.{extension}"
            )
        for extension in ("xyz", "json"):
            files[f"qm/neutral_thr_nme/{name}/trajectory.{extension}"] = (
                trajectories / f"{name}.{extension}"
            )
    for result in af3_audits or []:
        label = result.parent.name
        prefix = f"diagnostics/af3_water/{label}"
        for name in (
            "report.json",
            "pa6_l4_alpha_beta_water_test.cif",
            "view_af3_water_test.pml",
        ):
            destination = f"{prefix}/{name}"
            if destination in files:
                raise ValueError(f"Duplicate AF3 audit label: {label}")
            files[destination] = result / name
    return files


def build_package(files: dict[str, Path], output: Path) -> dict[str, object]:
    """Copy an explicit size-limited inventory and write integrity manifests.

    Args:
        files: Relative destination names and source files.
        output: New package directory, which must not exist.

    Returns:
        Source paths, file sizes and SHA-256 digests of every copied artifact.

    Raises:
        ValueError: If a path is unsafe or the size limit would be exceeded.
        FileExistsError: If the output already exists.
    """
    if output.exists():
        raise FileExistsError(output)
    rows = {}
    for name, source in files.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts or "\n" in name:
            raise ValueError(f"Unsafe download path: {name}")
        if name in {"MANIFEST.sha256", "provenance.json"}:
            raise ValueError(f"Reserved download path: {name}")
        size = source.stat().st_size
        if not source.is_file() or size > MAX_FILE_BYTES:
            raise ValueError(f"Not a small regular file: {source}")
        rows[name] = {
            "source": str(source.resolve()),
            "bytes": size,
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        }
    if sum(row["bytes"] for row in rows.values()) > MAX_TOTAL_BYTES:
        raise ValueError("Download inventory exceeds 25 MiB")
    output.mkdir(parents=True)
    for name, source in files.items():
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if hashlib.sha256(target.read_bytes()).hexdigest() != rows[name]["sha256"]:
            raise RuntimeError(f"Source changed during packaging: {source}")
    report = {
        "files": rows,
        "model_aliases": MODEL_ALIASES,
        "payload_bytes": sum(r["bytes"] for r in rows.values()),
        "purpose": "Visual inspection; no fitted parameters or validated TI1 minimum",
    }
    provenance = output / "provenance.json"
    provenance.write_text(json.dumps(report, indent=2) + "\n")
    checksums = {name: row["sha256"] for name, row in rows.items()}
    checksums["provenance.json"] = hashlib.sha256(provenance.read_bytes()).hexdigest()
    (output / "MANIFEST.sha256").write_text(
        "".join(f"{digest}  {name}\n" for name, digest in sorted(checksums.items()))
    )
    return report


def main() -> None:
    """Build the explicit laptop package from accepted QM inspection exports."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--qm-audit", required=True, type=Path)
    parser.add_argument("--qm-structures", required=True, type=Path)
    parser.add_argument("--qm-trajectories", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--protein-source", type=Path)
    parser.add_argument("--action-journal", type=Path)
    parser.add_argument("--af3-audit", type=Path, action="append", default=[])
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    report = build_package(
        selected_files(
            args.repo,
            args.qm_audit,
            args.qm_structures,
            args.qm_trajectories,
            args.protein_source,
            args.action_journal,
            args.af3_audit,
        ),
        args.output,
    )
    LOGGER.info(
        "Packaged %d files (%d bytes)", len(report["files"]), report["payload_bytes"]
    )


if __name__ == "__main__":
    main()
