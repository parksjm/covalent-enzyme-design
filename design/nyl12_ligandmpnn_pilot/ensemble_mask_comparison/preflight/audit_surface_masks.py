#!/usr/bin/env python3
"""Audit SASA-gated masks for the two accepted WT ensemble parents."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[4]

sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.design_surface_mask import (  # noqa: E402
    passes_all_copy_sasa,
    proposal_shell_positions,
    protein_molecule_from_pdb,
    read_bound_complex_heavy_coordinates,
    summarize_sasa,
)

LOGGER = logging.getLogger(__name__)
FOCUSED = {
    "alpha": [
        50,
        51,
        52,
        73,
        76,
        101,
        111,
        112,
        114,
        115,
        116,
        151,
        152,
        153,
        185,
        187,
        188,
    ],
    "beta": [39, 41],
}
PROTECTED = {"alpha": [113, 156, 186, *range(222, 233)], "beta": [1, 2, 3, 40, 42]}
CHAIN_GROUPS = {"alpha": "ABCD", "beta": "EFGH"}
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def sha256(path: Path) -> str:
    """Return SHA256 for an existing source or output file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_parent(path: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Calculate SASA and geometric mask evidence for one relaxed parent."""
    molecule = protein_molecule_from_pdb(path)
    sasa = summarize_sasa(molecule, expected_residues=1312)
    protein, ligands = read_bound_complex_heavy_coordinates(path)
    if {chain: len(residues) for chain, residues in protein.items()} != {
        **{chain: 232 for chain in "ABCD"},
        **{chain: 96 for chain in "EFGH"},
    }:
        raise ValueError(f"Unexpected chain lengths in {path}")
    sequence_by_family: dict[str, str] = {}
    for family, representative in (("alpha", "A"), ("beta", "E")):
        sequence = "".join(
            THREE_TO_ONE[protein[representative][n][0]]
            for n in sorted(protein[representative])
        )
        if any(
            "".join(THREE_TO_ONE[protein[chain][n][0]] for n in sorted(protein[chain]))
            != sequence
            for chain in CHAIN_GROUPS[family]
        ):
            raise ValueError(f"Equivalent {family} copies differ in {path}")
        sequence_by_family[family] = sequence
    sasas = {(row["chain"], row["residue_number"]): row for row in sasa["residues"]}
    protein_trees = {
        chain: cKDTree(
            np.concatenate([coords for _, coords in protein[chain].values()])
        )
        for chain in protein
    }
    ligand_tree = cKDTree(np.concatenate(list(ligands.values())))
    geometry: dict[str, dict[int, dict[str, float]]] = {}
    for family, group in CHAIN_GROUPS.items():
        for chain in group:
            for number, (_, coords) in protein[chain].items():
                ligand_distance = float(ligand_tree.query(coords)[0].min())
                other_chain_distance = min(
                    float(protein_trees[other].query(coords)[0].min())
                    for other in protein
                    if other != chain
                )
                geometry.setdefault(family, {}).setdefault(number, {})[chain] = {
                    "ligand_A": ligand_distance,
                    "other_chain_A": other_chain_distance,
                }
    detail = {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "sequence_by_family": sequence_by_family,
        "sasa": sasa,
        "geometry": geometry,
        "per_copy_rsa": {
            family: {
                str(number): {chain: sasas[(chain, number)]["rsa"] for chain in group}
                for number in sorted(protein[group[0]])
            }
            for family, group in CHAIN_GROUPS.items()
        },
    }
    return detail, protein, ligands


def run_audit(proposal_path: Path, output_dir: Path) -> dict[str, Any]:
    """Audit both bound parents and write surface exposure, masks, and hashes."""
    proposal = json.loads(proposal_path.read_text())
    details: dict[str, Any] = {}
    geometry: dict[str, Any] = {}
    rsa_values: dict[str, dict[int, list[float]]] = {
        family: {} for family in CHAIN_GROUPS
    }
    shell_candidates: dict[str, dict[str, set[int]]] = {
        family: {"6-12": set(), "12-16": set()} for family in CHAIN_GROUPS
    }
    parent_paths = [Path(row["path"]) for row in proposal["parents"]]
    for parent_row, path in zip(proposal["parents"], parent_paths, strict=True):
        detail, _, _ = audit_parent(path)
        parent_id = path.parents[1].name
        details[parent_id] = detail
        geometry[parent_id] = detail["geometry"]
        for family, group in CHAIN_GROUPS.items():
            for n in range(1, 233 if family == "alpha" else 97):
                vals = detail["per_copy_rsa"][family].get(str(n))
                if vals is None:
                    continue
                rsa_values[family].setdefault(n, []).extend(vals[c] for c in group)
            for shell in ("6-12", "12-16"):
                shell_candidates[family][shell] |= proposal_shell_positions(
                    proposal, parent_row["seed"], family, shell
                )
    mask: dict[str, Any] = {
        "focused_positions": FOCUSED,
        "hard_protected_positions": PROTECTED,
        "threshold_rsa_all_eight_copies": 0.25,
        "new_positions": {},
        "outer_shell_positions": {},
        "expanded_positions": {},
        "hard_zero_positions": {
            family: sorted(PROTECTED[family]) for family in PROTECTED
        },
        "zero_protection_violations": [],
        "candidate_decisions": {},
    }
    for family, candidate_shells in shell_candidates.items():
        decisions: dict[str, Any] = {}
        accepted: set[int] = set()
        accepted_outer: set[int] = set()
        for shell, candidates in candidate_shells.items():
            for number in sorted(candidates):
                key = f"{family}:{number}"
                values = rsa_values[family].get(number, [])
                nearest_interchain = [
                    geometry[parent][family][number][chain]["other_chain_A"]
                    for parent in details
                    for chain in CHAIN_GROUPS[family]
                ]
                shell_distances = [
                    geometry[parent][family][number][chain]["ligand_A"]
                    for parent in details
                    for chain in CHAIN_GROUPS[family]
                ]
                in_shell = (
                    any(6.0 < d <= 12.0 for d in shell_distances)
                    if shell == "6-12"
                    else any(12.0 < d <= 16.0 for d in shell_distances)
                )
                rejection_reasons = []
                if number in FOCUSED[family]:
                    rejection_reasons.append("already_in_focused_mask")
                if number in PROTECTED[family]:
                    rejection_reasons.append("hard_protected")
                if not passes_all_copy_sasa(values, 0.25):
                    rejection_reasons.append("rsa_below_0.25_in_at_least_one_copy")
                if any(distance <= 4.0 for distance in nearest_interchain):
                    rejection_reasons.append("interchain_heavy_contact_at_or_below_4_A")
                if not in_shell:
                    rejection_reasons.append("outside_requested_shell_in_both_parents")
                eligible = not rejection_reasons
                if eligible and shell == "6-12":
                    accepted.add(number)
                elif eligible:
                    accepted_outer.add(number)
                decisions[key] = {
                    "source_shell": shell,
                    "min_rsa_all_copies": round(min(values), 6),
                    "max_rsa_all_copies": round(max(values), 6),
                    "min_other_chain_heavy_A": round(min(nearest_interchain), 4),
                    "substrate_distance_range_A": [
                        round(min(shell_distances), 4),
                        round(max(shell_distances), 4),
                    ],
                    "eligible": eligible,
                    "rejection_reasons": rejection_reasons,
                }
        accepted -= accepted_outer
        expanded = set(FOCUSED[family]) | accepted | accepted_outer
        violation = expanded & set(PROTECTED[family])
        mask["candidate_decisions"][family] = decisions
        mask["new_positions"][family] = sorted(accepted)
        mask["outer_shell_positions"][family] = sorted(accepted_outer)
        mask["expanded_positions"][family] = sorted(expanded)
        mask["zero_protection_violations"].extend(
            f"{family}:{number}" for number in sorted(violation)
        )
    if mask["zero_protection_violations"]:
        raise ValueError("Final mask overlaps hard-protected positions")
    mask["counts"] = {
        family: {
            "focused": len(FOCUSED[family]),
            "6_12_A_additions": len(mask["new_positions"][family]),
            "12_16_A_additions": len(mask["outer_shell_positions"][family]),
            "expanded_per_copy": len(mask["expanded_positions"][family]),
            "hard_zero_positions": len(PROTECTED[family]),
        }
        for family in CHAIN_GROUPS
    }
    exposure = {
        "status": "complete",
        "method": "RDKit rdFreeSASA Lee-Richards / Protor; probe 1.4 A; protein ATOM records only",
        "reference": "Tien et al. 2013, PLOS ONE 8:e80635, Table 1 allowed maximum ASA",
        "radii_policy": "All atoms classified by Protor; nonpositive/unclassified radii are errors",
        "surface_eligibility": "RSA >= 0.25 in all four equivalent copies in both parents (8 values)",
        "interchain_policy": "Reject a new position if any protein heavy-atom contact to another chain is <=4.0 A in either parent",
        "parent_exposure": details,
        "mask_decisions": mask["candidate_decisions"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    exposure_path = output_dir / "surface_exposure.json"
    masks_path = output_dir / "final_masks.json"
    hashes_path = output_dir / "source_hashes.json"
    exposure_path.write_text(json.dumps(exposure, indent=2, sort_keys=True) + "\n")
    masks_path.write_text(json.dumps(mask, indent=2, sort_keys=True) + "\n")
    source_paths = [
        proposal_path,
        *parent_paths,
        Path(__file__).resolve(),
        ROOT / "src/nyl12_covalent/design_surface_mask.py",
    ]
    hashes = {str(path.resolve()): sha256(path) for path in source_paths}
    hashes.update(
        {str(path.resolve()): sha256(path) for path in (exposure_path, masks_path)}
    )
    hashes_path.write_text(json.dumps(hashes, indent=2, sort_keys=True) + "\n")
    LOGGER.info("Wrote exposure and final masks under %s", output_dir)
    return mask


def main() -> None:
    """Run the one-off native-parent audit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proposal", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    run_audit(args.proposal, args.output_dir)


if __name__ == "__main__":
    main()
