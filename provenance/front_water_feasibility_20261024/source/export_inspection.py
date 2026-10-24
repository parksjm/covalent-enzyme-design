"""Export diagnostic oxygen probes in the full fixed pocket, not a new seed."""

from __future__ import annotations
import json
import logging
from pathlib import Path
import parmed as pmd
from nyl12_covalent.local_repair_model import identity

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = Path(__file__).resolve().parent


def main() -> None:
    """Retain source geometry and label rejected oxygen positions for inspection."""
    source = (
        ROOT
        / "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation"
    )
    model = pmd.load_file(
        str(source / "model.prmtop"), xyz=str(source / "prepared.rst7")
    )
    for residue in model.residues:
        residue.chain, residue.number = identity(residue.idx)
    out = BASE / "inspection"
    out.mkdir(exist_ok=False)
    model.save(str(out / "full_pocket.pdb"), renumber=False)
    report = json.loads((BASE / "oxygen_feasibility.json").read_text())
    positions = [
        [-3.833305, 8.999226, -1.874442],
        report["best_clearance_positions"][0]["oxygen_A"],
    ]
    lines = []
    for i, (x, y, z) in enumerate(positions, 1):
        lines.append(
            f"HETATM{i:5d}  O   PRB X{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
        )
    (out / "rejected_oxygen_probes.pdb").write_text("".join(lines) + "END\n")
    (out / "view_front_search.pml").write_text("""load full_pocket.pdb, pocket
load rejected_oxygen_probes.pdb, probes
hide everything
show sticks, pocket within 5 of probes
hide everything, pocket and chain W
show spheres, probes
set sphere_scale, 0.3, probes
color red, probes
util.cnc pocket
label probes and resi 1, "User O: protein clash"
label probes and resi 2, "Best grid O: still clashes"
set valence, 0
zoom probes, 6
""")
    (out / "README.md").write_text(
        "Diagnostic full-pocket context. Both oxygen probes are rejected, not candidate waters. Probe1=user manual position; probe2=best clearance in tested front hemisphere. No water H are claimed for probes. Full-source water is retained as historical context. Original AmberAsp272 maps to chainE residue40; Thr233 maps E1. No protein/substrate coordinates were relaxed. Run @view_front_search.pml from this directory.\n"
    )
    logging.info("Exported full fixed context and two REJECTED oxygen probes")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
