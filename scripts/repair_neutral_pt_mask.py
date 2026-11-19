"""Compress an exact whole-residue movable set and retain a failed launch."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask

from prepare_model30_neutral_pt import compact_mask


def main() -> None:
    """Write an equivalent short mask and archive the initial parser failure."""
    root = Path(__file__).resolve().parents[1]
    run = root / "qmmm/pa6_flipped_model30/05_neutral_reaction_path"
    record = json.loads((run / "preparation.json").read_text())
    system = AmberParm(str(run / "inputs/system.prmtop"))
    active = set(record["active_atom_indices_one_based"])
    residues: list[int] = []
    for residue in system.residues:
        atoms = {atom.idx + 1 for atom in residue.atoms}
        if atoms & active:
            if not atoms <= active:
                raise ValueError("Partial residue requires a different encoding")
            residues.append(residue.idx + 1)
    mask = ":" + compact_mask(residues)[1:]
    actual = {
        index + 1
        for index, selected in enumerate(AmberMask(system, mask).Selection())
        if selected
    }
    if actual != active:
        raise ValueError("Compact mask differs from original movable atom set")
    source = run / "inputs/fixed_shell.mdin"
    target = run / "inputs/fixed_shell_compact_residues.mdin"
    if target.exists():
        raise FileExistsError(target)
    target.write_text(re.sub(r"active='[^']+'", f"active='{mask}'", source.read_text()))
    archive = run / "clutter/initial_mask_failure"
    archive.mkdir(parents=True, exist_ok=False)
    for name in ("job", "01_direct_pt"):
        (run / name).rename(archive / name)
    audit = {
        "reason": "Long explicit atom mask caused Amber stack underflow before optimization",
        "exact_same_atom_set": True,
        "movable_atoms": len(active),
        "compact_mask": mask,
        "compact_mask_length": len(mask),
        "corrected_input_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "archived_failure": str(archive),
        "root_verification": "Lead independently observed failure and approved compact equivalent residue mask; root logging attempt timed out.",
    }
    (run / "compact_mask_repair.json").write_text(json.dumps(audit, indent=2) + "\n")


if __name__ == "__main__":
    main()
