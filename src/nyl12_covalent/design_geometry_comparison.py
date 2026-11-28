"""Compare PA6 catalytic geometry across design variants and Boltz stages."""

from __future__ import annotations

import csv
import itertools
import json
import logging
from pathlib import Path
from typing import Any

from scripts import screen_pa6_binding_orientation as screen
from nyl12_covalent.native_thr_policy import geometry_eligible, policy_metadata

LOGGER = logging.getLogger(__name__)
SITE_NAMES = ("AE", "BF", "CG", "DH")
LIGAND_NAMES = ("C54", "O49", "N59", "C62")
ATOM_NAMES = ";".join(LIGAND_NAMES)
EXPECTED_LIGAND_SITES = {"I": "AE", "J": "BF", "K": "CG", "L": "DH"}


def _assign_ligands(
    ligands: list[str], d1_by_ligand: dict[str, dict[str, float]]
) -> dict[str, str]:
    """Assign ligands one-to-one to sites by minimum total Thr OG1–C54 distance."""
    if len(ligands) == 1:
        ligand = ligands[0]
        return {ligand: min(SITE_NAMES, key=lambda site: d1_by_ligand[ligand][site])}
    if len(ligands) != 4:
        raise ValueError(f"Expected one or four PA6 ligands, found {len(ligands)}")
    choices = (
        (sum(d1_by_ligand[ligand][site] for ligand, site in zip(ligands, sites)), sites)
        for sites in itertools.permutations(SITE_NAMES)
    )
    _, best_sites = min(choices)
    return dict(zip(ligands, best_sites, strict=True))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write rows with stable field order inferred from the first row."""
    if not rows:
        raise ValueError(f"Cannot write empty comparison table: {path.name}")
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def compare_models(
    entries: list[dict[str, Any]], reference_path: Path, output: Path
) -> dict[str, Any]:
    """Compare assigned and exhaustive PA6 geometry for Boltz design models.

    Args:
        entries: Model records containing ``variant``, ``seed``, ``stage``,
            absolute PDB ``path``, and explicit ``ligand_count`` (1 or 4).
        reference_path: Native Nyl12 catalytic reference mmCIF.
        output: New output directory for CSV tables and JSON summaries.

    Returns:
        Model-count and output metadata suitable for a campaign manifest.

    Raises:
        FileExistsError: If output already exists.
        FileNotFoundError: If an input structure or reference is missing.
        ValueError: If entry metadata, ligand inventory, or geometry is invalid.
    """
    if not entries:
        raise ValueError("At least one model entry is required")
    if not reference_path.is_file():
        raise FileNotFoundError(reference_path)
    reference = screen._cif_atoms(
        reference_path, "A", {residue for _, residue, _ in screen.SITE_RESIDUES}
    )
    all_rows: list[dict[str, Any]] = []
    occupied_rows: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []

    for entry in entries:
        required = {"variant", "seed", "stage", "path", "ligand_count"}
        missing = required - entry.keys()
        if missing:
            raise ValueError(f"Entry missing required keys: {sorted(missing)}")
        if entry["stage"] not in {"raw", "relaxed"}:
            raise ValueError(f"Unsupported stage: {entry['stage']}")
        expected_count = entry["ligand_count"]
        if expected_count not in {1, 4}:
            raise ValueError("ligand_count must explicitly be 1 or 4")
        path = Path(entry["path"])
        if not path.is_absolute():
            raise ValueError(f"Model path must be absolute: {path}")
        if not path.is_file():
            raise FileNotFoundError(path)
        atoms = screen.read_atoms(path)
        ligand_keys = sorted(
            {(chain, residue) for chain, residue, atom in atoms if atom == "C54"}
        )
        if len(ligand_keys) != expected_count:
            raise ValueError(
                f"{path}: expected {expected_count} PA6 ligands, found {len(ligand_keys)}"
            )
        for chain, residue in ligand_keys:
            absent = [
                name for name in LIGAND_NAMES if (chain, residue, name) not in atoms
            ]
            if absent:
                raise ValueError(f"{path}: ligand {chain}:{residue} missing {absent}")

        site_atoms = {
            site: screen.layout_site(atoms, index)
            for index, site in enumerate(SITE_NAMES)
        }
        pair_rows: list[dict[str, Any]] = []
        d1_by_ligand: dict[str, dict[str, float]] = {}
        for chain, residue in ligand_keys:
            ligand_id = f"{chain}:{residue}"
            ligand = {name: atoms[(chain, residue, name)][1] for name in LIGAND_NAMES}
            d1_by_ligand[ligand_id] = {}
            for site_name in SITE_NAMES:
                metrics = screen.measure(site_atoms[site_name], ligand, reference)
                metrics["selection_eligible"] = geometry_eligible(metrics["tier"])
                d1_by_ligand[ligand_id][site_name] = float(metrics["d1"])
                row = {
                    **entry,
                    "path": str(path),
                    "ligand": ligand_id,
                    "site": site_name,
                    "intended_pair": EXPECTED_LIGAND_SITES.get(chain) == site_name,
                    **metrics,
                }
                pair_rows.append(row)

        assignment = _assign_ligands(
            [f"{chain}:{residue}" for chain, residue in ligand_keys], d1_by_ligand
        )
        model_key = {key: site for key, site in assignment.items()}
        selected = [row for row in pair_rows if model_key[row["ligand"]] == row["site"]]
        for row in pair_rows:
            row["assigned_site"] = model_key[row["ligand"]]
            row["is_assigned"] = row["site"] == row["assigned_site"]
            row["assignment_policy"] = (
                "closest_site" if expected_count == 1 else "minimum_total_d1_one_to_one"
            )
        all_rows.extend(pair_rows)
        occupied_rows.extend(selected)
        count = len(selected)
        model_rows.append(
            {
                "variant": entry["variant"],
                "seed": entry["seed"],
                "stage": entry["stage"],
                "path": str(path),
                "ligand_count": count,
                "tier_a_count": sum(row["tier"] == "A" for row in selected),
                "tier_ab_count": sum(row["tier"] in {"A", "B"} for row in selected),
                "selection_eligible_count": sum(
                    row["selection_eligible"] for row in selected
                ),
                "joint_pass_count": sum(bool(row["joint_pass"]) for row in selected),
                "tier_a_fraction": sum(row["tier"] == "A" for row in selected) / count,
                "tier_ab_fraction": sum(row["tier"] in {"A", "B"} for row in selected)
                / count,
                "joint_pass_fraction": sum(bool(row["joint_pass"]) for row in selected)
                / count,
                "assignment": json.dumps(assignment, sort_keys=True),
            }
        )

    output.mkdir(parents=True, exist_ok=False)
    _write_csv(output / "all_sites.csv", all_rows)
    _write_csv(output / "occupied_sites.csv", occupied_rows)
    _write_csv(output / "model_summary.csv", model_rows)
    report = {
        "selection_policy": policy_metadata(),
        "model_count": len(model_rows),
        "occupied_site_count": sum(row["ligand_count"] for row in model_rows),
        "comparison_rows": len(all_rows),
        "reference_path": str(reference_path),
        "assignment_policy": {
            "one_ligand": "closest Thr233 OG1–PA6 C54 site",
            "four_ligands": "minimum-total-distance one-to-one permutation assignment",
            "intended_pairs": EXPECTED_LIGAND_SITES,
        },
        "replicate_unit": "model seed; sites within a model are nested and not independent replicates",
        "model_summary": model_rows,
    }
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info(
        "Compared %d models across %d ligand-site pairs", len(model_rows), len(all_rows)
    )
    return report
