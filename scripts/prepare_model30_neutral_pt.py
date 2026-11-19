"""Preserve the accepted neutral model30 shell for a bounded direct-PT test."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import re
import shutil

from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / (
    "qmmm/pa6_flipped_model30/04_md_screen_20261111/"
    "neutral_selected_replica02_frame00027_early_transient/"
    "autonomous_continuation/continued_mm_and_gated_qmmm/"
    "qmmm_shells/qmmm_dftb3_shell_12A/dftb3_mm"
)
OUTPUT = ROOT / "qmmm/pa6_flipped_model30/05_neutral_reaction_path"


def sha256(path: Path) -> str:
    """Return the exact file identity."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_mask(indices: list[int]) -> str:
    """Compress ordered one-based atom indices to an Amber mask."""
    groups: list[str] = []
    first = last = indices[0]
    for value in indices[1:]:
        if value == last + 1:
            last = value
        else:
            groups.append(str(first) if first == last else f"{first}-{last}")
            first = last = value
    groups.append(str(first) if first == last else f"{first}-{last}")
    return "@" + ",".join(groups)


def main() -> None:
    """Materialize the original movable shell and immutable launch sources."""
    preparation = json.loads((SOURCE / "preparation.json").read_text())
    command = preparation["command"]
    topology = Path(command[command.index("-p") + 1])
    original_start = Path(command[command.index("-c") + 1])
    endpoint = SOURCE / "optimized.ncrst"
    system = AmberParm(str(topology), xyz=str(original_start))
    indices = [
        index + 1
        for index, selected in enumerate(
            AmberMask(system, preparation["active_mask"]).Selection()
        )
        if selected
    ]
    qm = preparation["region_audit"]
    assert qm["real_atom_count"] == 98 and qm["link_H_count"] == 8
    assert qm["qmcharge"] == -2
    assert {atom["index"] for atom in qm["qm_atoms"]} <= set(indices)
    mask = compact_mask(indices)
    base, count = re.subn(
        r"active='[^']+'", f"active='{mask}'", (SOURCE / "minimize.mdin").read_text()
    )
    assert count == 1
    OUTPUT.mkdir(exist_ok=False)
    inputs = OUTPUT / "inputs"
    inputs.mkdir()
    for source, name in (
        (topology, "system.prmtop"),
        (endpoint, "accepted_es.ncrst"),
        (original_start, "original_shell_reference.ncrst"),
        (SOURCE / "minimize.mdin", "original_minimize.mdin"),
        (SOURCE / "preparation.json", "original_preparation.json"),
    ):
        shutil.copy2(source, inputs / name)
    (inputs / "fixed_shell.mdin").write_text(base)
    snapshot = OUTPUT / "source_snapshot"
    (snapshot / "scripts").mkdir(parents=True)
    shutil.copytree(
        ROOT / "src/nyl12_covalent",
        snapshot / "src/nyl12_covalent",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    for name in (
        "run_direct_thr233_proton_transfer.py",
        "export_direct_pt_progress.py",
        "plot_reaction_energies.py",
        "prepare_model30_neutral_pt.py",
    ):
        shutil.copy2(ROOT / "scripts" / name, snapshot / "scripts" / name)
    metadata = {
        "hypothesis": "Stepwise OG1-HG1 to neutral Thr929-N proton transfer",
        "coordinate": "d(OG1,HG1)-d(N,HG1)",
        "maximum_target_increment_A": 0.10,
        "final_target_A": 0.75,
        "only_reaction_restraint": "100 kcal/mol/A^2 difference coordinate",
        "release": "No reaction or positional restraints; fixed outer shell",
        "proton_inventory": qm["proton_inventory"],
        "qm_real_atoms": 98,
        "links": 8,
        "qmcharge": -2,
        "active_atom_indices_one_based": indices,
        "active_count": len(indices),
        "active_mask": mask,
        "original_active_expression": preparation["active_mask"],
        "shell_reference": str(original_start),
        "shell_reference_sha256": sha256(original_start),
        "fixed_shell_note": "Original spatial selection materialized once at original MM input; identical atoms move in all windows and release.",
        "acceptance": "Native convergence and H1/H2 N-bound; release additionally retains HG1 on N. No stationary-point order claim.",
        "bounded_decision": "Stop after direct PT/release; no automatic coupled fallback or umbrella sampling.",
        "sources": {
            str(path): sha256(path) for path in (topology, endpoint, original_start)
        },
        "snapshot_sha256": {
            str(path.relative_to(snapshot)): sha256(path)
            for path in snapshot.rglob("*")
            if path.is_file()
        },
    }
    (OUTPUT / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    logging.info("Prepared %s with %d fixed movable atoms", OUTPUT, len(indices))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
