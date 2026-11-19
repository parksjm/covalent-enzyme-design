"""Extend the corrected 12 Å shell for the explicit QM water."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import parmed
from parmed.amber.mask import AmberMask
from scipy.spatial import cKDTree

from nyl12_covalent.ash272_qmmm import audit_qm_region, compact_residue_mask

RUN = Path(__file__).resolve().parent
BASE = RUN.parents[3]
SOURCE_TOPOLOGY = BASE / "05_neutral_reaction_path/inputs/system.prmtop"
TARGET_TOPOLOGY = BASE / "07_ash274_comparison/00_system/ash274_solvated.prmtop"
PRIOR_PREPARATION = (
    BASE
    / "06_ash272_comparison/06_qmmm_optimizations/replica_1_tier_A/preparation.json"
)
OLD_AUDIT = RUN.parent / "corrected_shell/shell_audit.json"
START = RUN.parent / "corrected_shell/optimized.ncrst"
RADIUS = 12.0


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    """Write a non-overwritable shell audit for the 102-atom QM region."""
    output = RUN / "shell_audit_water_qm.json"
    if output.exists():
        raise FileExistsError(output)
    source = parmed.load_file(str(SOURCE_TOPOLOGY))
    target = parmed.load_file(str(TARGET_TOPOLOGY), xyz=str(START))
    prior = json.loads(PRIOR_PREPARATION.read_text())
    old = json.loads(OLD_AUDIT.read_text())
    source_ids = [
        i
        for i, selected in enumerate(
            AmberMask(source, prior["qm_region"]["source_mask"]).Selection()
        )
        if selected
    ]
    water = source.residues[1575]
    if water.name != "WAT" or {atom.name for atom in water.atoms} != {"O", "H1", "H2"}:
        raise ValueError("Unexpected source WAT1576 inventory")
    source_ids.extend(atom.idx for atom in water.atoms)
    region = audit_qm_region(
        source,
        target,
        source_ids,
        {970: "ASH", 1066: "ASH", 1162: "ASH", 1258: "ASH"},
        reference_charge=-2,
        expected_link_h_count=8,
    )
    if (region.real_atom_count, region.formal_charge, region.electron_count) != (
        102,
        -1,
        378,
    ):
        raise ValueError("Unexpected QM water electron inventory")
    old_selected = {
        i
        for i, chosen in enumerate(
            AmberMask(target, old["new_active"]["mask"]).Selection()
        )
        if chosen
    }
    xyz = np.asarray(target.coordinates, dtype=float)
    distances = cKDTree(xyz[list(region.atom_indices)]).query(xyz, k=1)[0]
    near = set(np.flatnonzero(distances <= RADIUS).tolist())
    missing = near - old_selected
    added_residues = {target.atoms[i].residue.idx + 1 for i in missing}
    all_residues = {
        target.atoms[i].residue.idx + 1 for i in old_selected
    } | added_residues
    mask = compact_residue_mask(
        target,
        [
            atom.idx
            for residue in all_residues
            for atom in target.residues[residue - 1].atoms
        ],
    )
    new_selected = {
        i for i, chosen in enumerate(AmberMask(target, mask).Selection()) if chosen
    }
    if near - new_selected or set(region.atom_indices) - new_selected:
        raise ValueError("Expanded shell still excludes a QM-near atom")
    report = {
        "purpose": "12 Å complete-residue shell around 102-atom water-in-QM E-S region",
        "radius_angstrom": RADIUS,
        "qm": {
            "mask": region.mask,
            "charge": -1,
            "real_atoms": 102,
            "electrons": 378,
            "link_h": 8,
        },
        "old_active": {
            "mask": old["new_active"]["mask"],
            "atom_count": len(old_selected),
        },
        "new_active": {
            "mask": mask,
            "atom_count": len(new_selected),
            "residue_count": len(all_residues),
            "added_residues": sorted(added_residues),
            "added_atom_count": len(new_selected - old_selected),
            "near_qm_omissions_at_start": len(near - new_selected),
        },
        "source_sha256": {
            str(p.resolve()): sha256(p)
            for p in (
                Path(__file__),
                SOURCE_TOPOLOGY,
                TARGET_TOPOLOGY,
                PRIOR_PREPARATION,
                OLD_AUDIT,
                START,
            )
        },
    }
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
