"""Export the completed unbiased activation release and classify its endpoint."""

from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
import re
import subprocess
import sys

import parmed
import numpy as np
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from scripts.export_downhill_minimization_vmd import dcd_counts  # noqa: E402
from scripts.queue_modular_es_qmmm_after_mm import endpoint_geometry  # noqa: E402
from nyl12_covalent.proton_coordination import coordination_transfer  # noqa: E402

GRID = ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid"
CASE = GRID / "diagnostic_activation_release"
TOPOLOGY = GRID.parent / "00_system/ash274_solvated.prmtop"
LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Hash one input or output artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Write a portable trajectory with retained relay water and final acceptance."""
    logging.basicConfig(level=logging.INFO)
    folder = CASE / "vmd"
    folder.mkdir(exist_ok=False)
    script = folder / "export.cpptraj"
    script.write_text(
        f"parm {TOPOLOGY}\ntrajin {CASE / 'starting.ncrst'}\n"
        f"trajin {CASE / 'optimization.nc'}\ntrajin {CASE / 'optimized.ncrst'}\n"
        "strip !:1-1316,1576\n"
        f"parmwrite out {folder / 'path.parm7'}\ntrajout {folder / 'path.dcd'} dcd\nrun\nquit\n"
    )
    with (folder / "cpptraj.log").open("w") as log:
        subprocess.run(
            ["/data/jpc/code/ambertools26/venv/bin/cpptraj", "-i", str(script)],
            env=os.environ.copy()
            | {"LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib"},
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    frames, atoms = dcd_counts(folder / "path.dcd")
    old_manifest = json.loads((GRID / "progress/row_004/manifest.json").read_text())
    keep = old_manifest["stripped_to_original_atom_indices_0based"]
    if atoms != len(keep):
        raise ValueError("Stripped trajectory and atom mapping differ")
    mapping = {old + 1: new for new, old in enumerate(keep)}
    wo, wh, og = [mapping[i] for i in (20317, 20318, 13916)]
    text = (GRID / "progress/row_004/view_path.vmd").read_text()
    text = re.sub(
        r'puts "Loaded [^\n]*',
        f'puts "Loaded {frames} diagnostic frames: initial, saved optimizer frames, endpoint."',
        text,
    )
    text = text.replace("mol color ColorID 7", "mol color Element")
    text = text.replace("mol color ColorID 4", "mol color Element")
    text = text.replace(
        "mol representation DynamicBonds",
        "mol color Element\nmol representation DynamicBonds",
    )
    cleanup = (
        f"set sel [atomselect $molid all]\nset bonds [$sel getbonds]\n"
        f"foreach pair {{{{{wo} {wh}}} {{{og} {wh}}}}} {{\n"
        " lassign $pair a b\n foreach {i j} [list $a $b $b $a] {\n"
        "  set keep {}; foreach k [lindex $bonds $i] {if {$k != $j} {lappend keep $k}}\n"
        "  lset bonds $i $keep\n }\n}\n$sel setbonds $bonds\n$sel delete\n"
    )
    text = text.replace("mol delrep 0 $molid\n", "mol delrep 0 $molid\n" + cleanup)
    for center in (wo, og):
        text += (
            f"\nmol color Element\nmol representation DynamicBonds 1.4 0.1 1\n"
            f'mol selection "index {center} {wh}"\nmol addrep $molid\n'
            f"label add Bonds $molid/{center} $molid/{wh}\n"
        )
    text += "\ncolor Element H white\ncolor Element N blue\ncolor Element O red\n"
    (folder / "view_optimization.vmd").write_text(text)
    process = json.loads((CASE / "process.json").read_text())
    inventory = json.loads((CASE / "proton_inventory.json").read_text())
    topology = parmed.load_file(str(TOPOLOGY), xyz=str(CASE / "optimized.ncrst"))
    protocol = json.loads((GRID / "protocol.json").read_text())
    qm = AmberMask(topology, protocol["qm"]["mask"]).Selection()
    hydrogens = [a.idx for a in topology.atoms if qm[a.idx] and a.atomic_number == 1]
    heavy = np.array([a.idx for a in topology.atoms if a.atomic_number > 1])
    xyz = np.asarray(topology.coordinates)
    bound = {13905: [], 13916: []}
    for h in hydrogens:
        distances = np.linalg.norm(xyz[heavy] - xyz[h], axis=1)
        k = int(distances.argmin())
        owner = int(heavy[k]) + 1
        if owner in bound and distances[k] <= 1.25:
            bound[owner].append(h + 1)
    occupancy = coordination_transfer(xyz, 13915, 13904, hydrogens)
    energies = re.findall(
        r"EPtot\s*=\s*([-+\d.EeDd]+)", (CASE / "minimized.mdout").read_text()
    )
    baseline = json.loads((GRID / "grid_points.json").read_text())["points"][
        "q_001_attack_001"
    ]["physical_energy_kcal_mol"]
    identities = {
        a: inventory["endpoint"]["proton_nearest_heavy"][a]["nearest_heavy_atom_id"]
        for a in ("Thr233_HG1", "water1576_H1", "Thr233_N_H1", "Thr233_N_H2")
    }
    result = {
        "process_and_force_converged": process["exit_status"] == 0
        and process["force_converged"],
        "scc_converged": process["scc_converged"],
        "activation_retained": False,
        "endpoint_class": "E-S with exchanged Thr/water proton identities",
        "classification_criteria": "HG1 nearest water O; water H1 nearest OG1; original two N protons retained",
        "identities": identities,
        "covalent_real_QM_H_by_heavy_atom_id": bound,
        "smooth_CN_N": occupancy.acceptor,
        "smooth_CN_OG1": occupancy.donor,
        "smooth_occupancy_coordinate": occupancy.value,
        "relative_endpoint_energy_kcal_mol": float(energies[-1].replace("D", "E"))
        - baseline,
        "geometry": endpoint_geometry(TOPOLOGY, CASE / "optimized.ncrst"),
        "residue_number_note": "Use Amber residue idx+1. Raw proton_inventory residue.number fields follow ParmEd numbering.",
        "frame_count": frames,
        "stripped_atom_count": atoms,
        "source_sha256": {
            str(p): digest(p)
            for p in (
                Path(__file__),
                TOPOLOGY,
                CASE / "starting.ncrst",
                CASE / "optimized.ncrst",
                CASE / "optimization.nc",
            )
        },
    }
    expected = {
        "Thr233_HG1": 20317,
        "water1576_H1": 13916,
        "Thr233_N_H1": 13905,
        "Thr233_N_H2": 13905,
    }
    if identities != expected:
        raise ValueError("Endpoint differs from inspected diagnostic chemistry")
    if len(bound[13905]) != 2 or len(bound[13916]) != 1:
        raise ValueError("All-real-QM-H audit does not match neutral Thr-NH2/OH state")
    (CASE / "scientific_assessment.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    (folder / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    (folder / "MANIFEST.sha256").write_text(
        "".join(
            f"{digest(p)}  {p.name}\n" for p in sorted(folder.iterdir()) if p.is_file()
        )
    )
    LOGGER.info(
        "Exported %d frames, %d atoms; endpoint returned to E-S chemistry",
        frames,
        atoms,
    )


if __name__ == "__main__":
    main()
