"""Prepare and run one audited ASH274 DFTB3/MM E-S optimization."""

from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any

import numpy as np
import parmed
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[6]
BASE = ROOT / "qmmm/pa6_flipped_model30"
CASE = Path(__file__).resolve().parent
BASE_TOPOLOGY = BASE / "05_neutral_reaction_path/inputs/system.prmtop"
PRIOR_TOPOLOGY = BASE / "06_ash272_comparison/00_system/ash272_solvated.prmtop"
TOPOLOGY = BASE / "07_ash274_comparison/00_system/ash274_solvated.prmtop"
SYSTEM_ACCEPTANCE = BASE / "07_ash274_comparison/00_system/preparation_acceptance.json"
PRIOR_CASE = BASE / "06_ash272_comparison/06_qmmm_optimizations/replica_1_tier_A"
PRIOR_PREPARATION = PRIOR_CASE / "preparation.json"
PRIOR_MDIN = PRIOR_CASE / "minimize.mdin"
MM_CASE = BASE / "07_ash274_comparison/05_selected_mm/replica_1_tier_A"
MM_PREPARATION = MM_CASE / "preparation.json"
MM_STATUS = MM_CASE / "status.json"
MM_RESTART = MM_CASE / "optimized.rst7"
SOURCE_PREPARATION = BASE / "05_neutral_reaction_path/preparation.json"
SANDER = Path("/data/jpc/code/ambertools26/venv/bin/sander")
DLFIND = Path("/data/jpc/code/ambertools26/venv/lib/libdlfind.so")
SANDER_LIB = Path("/data/jpc/code/ambertools26/venv/lib/libsander.so")
PROTONATION = {970: "ASH", 1066: "ASH", 1162: "ASH", 1258: "ASH"}
SITE_PAIRS = (
    (929, 968, 970),
    (1025, 1064, 1066),
    (1121, 1160, 1162),
    (1217, 1256, 1258),
)
LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_indices(topology: Any, mask: str) -> list[int]:
    """Resolve an Amber mask to zero-based atom indices."""
    selection = AmberMask(topology, mask).Selection()
    indices = [index for index, selected in enumerate(selection) if selected]
    if not indices:
        raise ValueError(f"Amber mask selected no atoms: {mask}")
    return indices


def atom_key(atom: Any) -> tuple[int, str]:
    """Return one-based residue number and atom name for identity mapping."""
    return int(atom.residue.idx) + 1, str(atom.name).strip()


def cut_identities(topology: Any, atom_indices: list[int]) -> set[tuple[Any, Any]]:
    """Return stable atom-identity pairs for every QM/MM boundary bond."""
    selected = set(atom_indices)
    cuts = set()
    for bond in topology.bonds:
        first, second = bond.atom1.idx in selected, bond.atom2.idx in selected
        if first != second:
            cuts.add(tuple(sorted((atom_key(bond.atom1), atom_key(bond.atom2)))))
    return cuts


def validate_microstate(topology: Any, acceptance: dict[str, Any]) -> None:
    """Require four ASP272/ASH274 pairs and neutral Thr N termini."""
    if acceptance.get("protonated_residue") != 274 or not acceptance.get("neutral_Thr"):
        raise ValueError("ASH274 system acceptance does not match requested microstate")
    for thr, asp272, ash274 in SITE_PAIRS:
        thr_res, asp_res, ash_res = (
            topology.residues[index - 1] for index in (thr, asp272, ash274)
        )
        if (thr_res.name, asp_res.name, ash_res.name) != ("THR", "ASP", "ASH"):
            raise ValueError(
                f"Unexpected site residue states at {thr}/{asp272}/{ash274}"
            )
        thr_atoms = {atom.name: atom for atom in thr_res.atoms}
        if "H3" in thr_atoms or not {"N", "H1", "H2", "OG1"} <= thr_atoms.keys():
            raise ValueError(f"Neutral Thr N inventory mismatch at residue {thr}")
        for name, expected in (("N", -0.4844), ("H1", 0.1224), ("H2", 0.1224)):
            if abs(thr_atoms[name].charge - expected) > 1e-6:
                raise ValueError(f"Neutral Thr N charge mismatch at {thr}:{name}")
        ash_atoms = {atom.name: atom for atom in ash_res.atoms}
        if "HD2" not in ash_atoms or not any(
            bond.atom1 is ash_atoms["HD2"]
            and bond.atom2 is ash_atoms["OD2"]
            or bond.atom2 is ash_atoms["HD2"]
            and bond.atom1 is ash_atoms["OD2"]
            for bond in topology.bonds
        ):
            raise ValueError(f"ASH274 HD2-OD2 bond absent at {ash274}")
        if "HD2" in {atom.name for atom in asp_res.atoms}:
            raise ValueError(f"ASP272 is unexpectedly protonated at {asp272}")


def point(structure: Any, residue: int, atom_name: str) -> np.ndarray:
    """Return coordinates for one exact residue/atom identity."""
    atoms = [
        atom for atom in structure.residues[residue - 1].atoms if atom.name == atom_name
    ]
    if len(atoms) != 1:
        raise ValueError(f"Expected one atom {residue}:{atom_name}, found {len(atoms)}")
    return np.asarray(structure.coordinates[atoms[0].idx], dtype=float)


def geometry(structure: Any) -> dict[str, float]:
    """Measure E-S Tier geometry and scissile-bond lengths."""
    sys.path.insert(0, str(ROOT / "src"))
    from nyl12_covalent.active_site_geometry import (
        angle_degrees,
        projected_herschlag_phi,
    )

    og1 = point(structure, 929, "OG1")
    carbon = point(structure, 1314, "C")
    oxygen = point(structure, 1314, "O")
    nitrogen = point(structure, 1315, "N")
    adjacent = point(structure, 1314, "C5")
    return {
        "d1": float(np.linalg.norm(og1 - carbon)),
        "d2": float(np.linalg.norm(point(structure, 113, "N") - oxygen)),
        "d3": float(np.linalg.norm(point(structure, 186, "ND2") - oxygen)),
        "attack": angle_degrees(og1, carbon, oxygen),
        "projected_phi": projected_herschlag_phi(
            tuple(og1), tuple(oxygen), tuple(carbon), tuple(nitrogen), tuple(adjacent)
        ),
        "scissile_CO_A": float(np.linalg.norm(carbon - oxygen)),
        "scissile_CN_A": float(np.linalg.norm(carbon - nitrogen)),
    }


def _mdin(qm_mask: str, active_mask: str, charge: int) -> str:
    """Render the validated periodic PME/DL-FIND DFTB3/MM protocol."""
    return (
        "ASH274 ASP272/ASH274 neutral-Thr E-S Tier B DFTB3/MM optimization\n"
        "&cntrl\n nmropt=0, imin=1, ntmin=5, ncyc=10, dx0=0.02, drms=0.05,\n"
        " maxcyc=2000, ntb=1, igb=0, ntc=1, ntf=1, cut=10.0,\n"
        " ifqnt=1, jfastw=4, ntr=0, ntpr=25, ntwx=25, ntwr=100, ioutfm=1, ntxo=2,\n/\n"
        "&qmmm\n"
        f" qmmask='{qm_mask}', qmcharge={charge}, spin=1, qm_theory='DFTB3',\n"
        " qmmm_int=1, adjust_q=2, qm_ewald=1, qm_pme=1, qmcut=10.0,\n"
        " scfconv=1.0d-10, lnk_dis=1.09, lnk_atomic_no=1, lnk_method=1,\n/\n"
        "&dlfind\n"
        f" active='{active_mask}', optalg='LBFGS', crdrep='CART', tol=1.d-4,\n"
        " tole=1.d-5, maxcycle=2000, maxstep=0.5, verbosity=1,\n/\n"
    )


def prepare() -> dict[str, Any]:
    """Audit the target microstate and write immutable QM/MM inputs."""
    from nyl12_covalent.ash272_qmmm import audit_qm_region
    from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

    for path in (CASE / "preparation.json", CASE / "minimized.mdout"):
        if path.exists():
            raise FileExistsError(path)
    if (CASE / "starting.ncrst").exists() and sha256(CASE / "starting.ncrst") != sha256(
        MM_RESTART
    ):
        raise ValueError("Staged original restart hash differs from source")
    status = json.loads(MM_STATUS.read_text())
    mm_prep = json.loads(MM_PREPARATION.read_text())
    if status.get("exit_status") != 0 or not status["geometric_acceptance"]["accepted"]:
        raise ValueError("The ASH274 MM endpoint is not accepted")
    if mm_prep["geometric_acceptance"]["tier_after_mm"] != "B":
        raise ValueError("Expected the declared post-MM Tier B endpoint")
    if sha256(MM_RESTART) != mm_prep["output_sha256"]["optimized.rst7"]:
        raise ValueError("Selected ASH274 MM restart hash mismatch")

    source = parmed.load_file(str(BASE_TOPOLOGY))
    prior_top = parmed.load_file(str(PRIOR_TOPOLOGY))
    target_top = parmed.load_file(str(TOPOLOGY))
    baseline = json.loads(SOURCE_PREPARATION.read_text())
    prior = json.loads(PRIOR_PREPARATION.read_text())
    source_indices = selected_indices(source, prior["qm_region"]["source_mask"])
    link_count = int(prior["qm_region"]["link_h_count"])
    reference_charge = int(baseline["qmcharge"])
    prior_region = audit_qm_region(
        source,
        prior_top,
        source_indices,
        {968: "ASH", 1064: "ASH", 1160: "ASH", 1256: "ASH"},
        reference_charge=reference_charge,
        expected_link_h_count=link_count,
    )
    target_region = audit_qm_region(
        source,
        target_top,
        source_indices,
        PROTONATION,
        reference_charge=reference_charge,
        expected_link_h_count=link_count,
    )
    if (
        prior_region.formal_charge != prior["qm_region"]["charge"]
        or target_region.formal_charge != prior_region.formal_charge
        or target_region.real_atom_count != prior_region.real_atom_count
        or target_region.link_h_count != prior_region.link_h_count
        or target_region.electron_count != prior_region.electron_count
        or target_region.multiplicity != prior_region.multiplicity
    ):
        raise ValueError("Prior and target QM charge/link/electron inventories differ")
    prior_cuts = cut_identities(prior_top, list(prior_region.atom_indices))
    target_cuts = cut_identities(target_top, list(target_region.atom_indices))
    if prior_cuts != target_cuts:
        raise ValueError("QM/MM cut atom identities differ from validated ASH272 setup")
    prior_keys = {
        atom_key(prior_top.atoms[index]) for index in prior_region.atom_indices
    }
    target_keys = {
        atom_key(target_top.atoms[index]) for index in target_region.atom_indices
    }
    expected_swap = {(968, "HD2")}, {(970, "HD2")}
    if (prior_keys - target_keys, target_keys - prior_keys) != expected_swap:
        raise ValueError(
            "QM atom identity mapping differs beyond the declared proton swap"
        )

    audit = json.loads((CASE / "shell_audit.json").read_text())
    active_mask = audit["new_active"]["mask"]
    active_indices = selected_indices(target_top, active_mask)
    if len(active_indices) != int(audit["new_active"]["atom_count"]):
        raise ValueError("Corrected shell mask atom count differs from audit")
    if not set(target_region.atom_indices) <= set(active_indices):
        raise ValueError("Target QM region is not fully inside the 12A movable region")
    acceptance = json.loads(SYSTEM_ACCEPTANCE.read_text())
    validate_microstate(target_top, acceptance)
    if len(target_top.atoms) != 104934:
        raise ValueError("Unexpected full-solvent ASH274 topology atom count")
    if (
        abs(
            sum(atom.charge for atom in target_top.atoms)
            - acceptance["inherited_charge_residual_e"]
        )
        > 1e-5
    ):
        raise ValueError(
            "Target topology total charge differs from accepted ASH274 charge"
        )
    old_box, new_box = list(map(float, prior_top.box)), list(map(float, target_top.box))
    if any(abs(a - b) > 2e-5 for a, b in zip(old_box, new_box, strict=True)):
        raise ValueError("Target periodic box differs from validated ASH272 topology")

    initial_path = CASE / "starting.ncrst"
    if not initial_path.exists():
        shutil.copy2(MM_RESTART, initial_path)
    initial = parmed.load_file(str(TOPOLOGY), xyz=str(initial_path))
    if len(initial.atoms) != len(target_top.atoms):
        raise ValueError("The MM restart does not match the full solvated topology")
    from scipy.spatial import cKDTree

    coordinates = np.asarray(initial.coordinates, dtype=float)
    near = set(
        np.flatnonzero(
            cKDTree(coordinates[list(target_region.atom_indices)]).query(
                coordinates, k=1
            )[0]
            <= 12.0
        ).tolist()
    )
    if near - set(active_indices):
        raise ValueError(
            f"Corrected shell omits {len(near - set(active_indices))} atoms within 12 A at start"
        )
    initial_geometry = geometry(initial)
    if (
        classify_geometry(initial_geometry) != "B"
        or initial_geometry["projected_phi"] <= 0
    ):
        raise ValueError(
            "The starting restart does not pass Tier B positive-phi geometry"
        )

    mdin = CASE / "minimize.mdin"
    mdin.write_text(_mdin(target_region.mask, active_mask, target_region.formal_charge))
    sources = (
        BASE_TOPOLOGY,
        PRIOR_TOPOLOGY,
        TOPOLOGY,
        SYSTEM_ACCEPTANCE,
        PRIOR_PREPARATION,
        PRIOR_MDIN,
        SOURCE_PREPARATION,
        MM_PREPARATION,
        MM_STATUS,
        MM_RESTART,
        SANDER,
        DLFIND,
        SANDER_LIB,
        Path(__file__),
        CASE / "run.sh",
        initial_path,
        mdin,
    )
    record = {
        "case": CASE.name,
        "source_mm_case": str(MM_CASE.resolve()),
        "source_mm_restart_sha256": sha256(MM_RESTART),
        "microstate": {
            "ASP272": "ASP",
            "ASP274": "ASH(HD2-OD2)",
            "Thr233_N": "neutral NH2",
        },
        "system": {
            "atoms": len(target_top.atoms),
            "residues": len(target_top.residues),
            "total_charge_e": sum(atom.charge for atom in target_top.atoms),
            "full_solvent": True,
        },
        "qm_region": {
            "mask": target_region.mask,
            "charge": target_region.formal_charge,
            "real_atoms": target_region.real_atom_count,
            "link_h": target_region.link_h_count,
            "electrons": target_region.electron_count,
            "multiplicity": target_region.multiplicity,
            "link_cuts": [list(pair) for pair in sorted(target_cuts)],
        },
        "prior_comparison": {
            "charge_equal": True,
            "prior_charge": prior_region.formal_charge,
            "target_charge": target_region.formal_charge,
            "link_cut_identities_equal": True,
            "only_qm_atom_identity_change": {
                "prior_only": ["968:HD2"],
                "target_only": ["970:HD2"],
            },
        },
        "movable_region": {
            "radius_angstrom": 12.0,
            "mask": active_mask,
            "atom_count": len(active_indices),
            "residue_count": audit["new_active"]["residue_count"],
            "audit_sha256": sha256(CASE / "shell_audit.json"),
            "near_qm_atom_omissions_at_start": 0,
            "all_qm_atoms_movable": True,
        },
        "periodic_embedding": {
            "ntb": 1,
            "igb": 0,
            "qmmm_int": 1,
            "adjust_q": 2,
            "qm_ewald": 1,
            "qm_pme": 1,
            "qmcut_angstrom": 10.0,
            "prior_box": old_box,
            "target_box": new_box,
            "max_box_delta_angstrom": max(
                abs(a - b) for a, b in zip(old_box, new_box, strict=True)
            ),
        },
        "geometry_before_qmmm": initial_geometry,
        "method": "Amber Sander26 native DL-FIND / DFTB3-3OB-3-1 QM/MM",
        "optimization": {
            "algorithm": "LBFGS",
            "tol": 1e-4,
            "tole": 1e-5,
            "maxcycle": 2000,
            "restraints": "none",
        },
        "acceptance_gates": {
            "process": "Sander return code 0 and normal final results",
            "force_convergence": "DL-FIND reports 'Convergence reached' and 'Optimisation converged' before maxcycle=2000",
            "proton_retention": "E-site HD2 remains within 1.25 A of ASH274 OD2 and OD2 is its nearest listed catalytic acceptor",
            "neutral_thr_n": "Thr929 retains exactly H1/H2 at N-H <=1.25 A and no H3",
            "geometry": "existing Tier A/B classifier and positive projected phi; scissile C-N remains <=1.6 A",
        },
        "source_sha256": {str(path.resolve()): sha256(path) for path in sources},
        "status": "prepared; launch follows in tmux",
        "original_starting_sha256": sha256(CASE.parent / "starting.ncrst"),
        "starting_source": str((CASE.parent / "starting.ncrst").resolve()),
    }
    (CASE / "preparation.json").write_text(json.dumps(record, indent=2) + "\n")
    (CASE / "acceptance_criteria.json").write_text(
        json.dumps(record["acceptance_gates"], indent=2) + "\n"
    )
    return record


def assess() -> dict[str, Any]:
    """Apply the declared chemical, geometric, and force-convergence gates."""
    from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

    topology = parmed.load_file(str(TOPOLOGY))
    endpoint = parmed.load_file(str(TOPOLOGY), xyz=str(CASE / "optimized.ncrst"))
    validate_microstate(topology, json.loads(SYSTEM_ACCEPTANCE.read_text()))

    def distance(first: tuple[int, str], second: tuple[int, str]) -> float:
        """Measure a named interatomic distance in the endpoint."""
        return float(np.linalg.norm(point(endpoint, *first) - point(endpoint, *second)))

    proton_distances = {
        "HD2_to_ASH274_OD2_A": distance((970, "HD2"), (970, "OD2")),
        "HD2_to_Asp272_nearest_O_A": min(
            distance((970, "HD2"), (968, "OD1")), distance((970, "HD2"), (968, "OD2"))
        ),
        "HD2_to_Thr233_N_A": distance((970, "HD2"), (929, "N")),
        "HD2_to_Thr233_OG1_A": distance((970, "HD2"), (929, "OG1")),
        "HD2_to_substrate_N_A": distance((970, "HD2"), (1315, "N")),
        "Thr233_N_H1_A": distance((929, "N"), (929, "H1")),
        "Thr233_N_H2_A": distance((929, "N"), (929, "H2")),
    }
    geometry_metrics = geometry(endpoint)
    tier = classify_geometry(geometry_metrics)
    output = (CASE / "minimized.mdout").read_text(errors="replace")
    force = "Convergence reached" in output and "Optimisation converged" in output
    other_acceptors = [
        value
        for key, value in proton_distances.items()
        if key.startswith("HD2_to_") and key != "HD2_to_ASH274_OD2_A"
    ]
    gates = {
        "force_converged": force,
        "proton_retained_on_ash274": proton_distances["HD2_to_ASH274_OD2_A"] <= 1.25
        and proton_distances["HD2_to_ASH274_OD2_A"] < min(other_acceptors),
        "neutral_thr233_n_retained": proton_distances["Thr233_N_H1_A"] <= 1.25
        and proton_distances["Thr233_N_H2_A"] <= 1.25
        and "H3" not in {atom.name for atom in topology.residues[928].atoms},
        "tier_a_or_b_positive_phi": tier in {"A", "B"}
        and geometry_metrics["projected_phi"] > 0,
        "scissile_cn_intact": geometry_metrics["scissile_CN_A"] <= 1.6,
    }
    return {
        "accepted": all(gates.values()),
        "gates": gates,
        "tier": tier,
        "geometry": geometry_metrics,
        "proton_distances": proton_distances,
        "force_convergence_marker": force,
    }


def main() -> None:
    """Prepare, launch, and record one DFTB3/MM optimization and its gates."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if not (CASE / "preparation.json").exists():
        prepare()
    command = [
        str(SANDER),
        "-O",
        "-i",
        "minimize.mdin",
        "-p",
        str(TOPOLOGY),
        "-c",
        "starting.ncrst",
        "-o",
        "minimized.mdout",
        "-r",
        "optimized.ncrst",
        "-x",
        "optimization.nc",
        "-inf",
        "mdinfo",
    ]
    process = {
        "command": command,
        "cwd": str(CASE.resolve()),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
    }
    (CASE / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    environment = os.environ.copy() | {
        "AMBERHOME": "/data/jpc/code/ambertools26/venv",
        "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
        "OMP_NUM_THREADS": "1",
    }
    with (CASE / "engine.log").open("w") as log:
        result = subprocess.run(
            command,
            cwd=CASE,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    process.update(
        returncode=result.returncode,
        ended_utc=datetime.now(timezone.utc).isoformat(),
        status="complete" if result.returncode == 0 else "failed",
    )
    (CASE / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    if result.returncode != 0:
        (CASE / "status.json").write_text(
            json.dumps(
                {"exit_status": result.returncode, "scientific_acceptance": None},
                indent=2,
            )
            + "\n"
        )
        raise RuntimeError(f"Sander exited {result.returncode}")
    mdout = (CASE / "minimized.mdout").read_text(errors="replace")
    normal = (
        ("FINAL RESULTS" in mdout or "Run   done at" in mdout)
        and "FATAL" not in mdout
        and "NaN" not in mdout
    )
    acceptance = assess()
    acceptance["normal_engine_completion"] = normal
    acceptance["accepted"] = acceptance["accepted"] and normal
    if not normal:
        acceptance["gates"]["normal_engine_completion"] = False
    report = {
        "process": {
            "exit_status": result.returncode,
            "normal_engine_completion": normal,
        },
        "force_convergence": acceptance["gates"]["force_converged"],
        "chemistry_and_geometry": acceptance["gates"],
        "tier_a_b": {
            "tier": acceptance["tier"],
            "accepted": acceptance["gates"]["tier_a_or_b_positive_phi"],
        },
        "water_1576": {
            "included_in_active_shell": True,
            "initial_og1_distance_A": json.loads(
                (CASE / "shell_audit.json").read_text()
            )["restart_audits"]["original_starting"]["water_1576_OG1_distance_A"],
        },
        "case": CASE.name,
        "preparation_sha256": sha256(CASE / "preparation.json"),
        "acceptance": acceptance,
        "status": "accepted" if acceptance["accepted"] else "completed_but_gate_failed",
        "qm_mm_launched": True,
    }
    (CASE / "acceptance.json").write_text(json.dumps(report, indent=2) + "\n")
    artifacts = (
        "minimize.mdin",
        "starting.ncrst",
        "minimized.mdout",
        "optimized.ncrst",
        "optimization.nc",
        "engine.log",
        "mdinfo",
        "acceptance.json",
    )
    process["output_sha256"] = {
        name: sha256(CASE / name) for name in artifacts if (CASE / name).is_file()
    }
    (CASE / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    (CASE / "status.json").write_text(
        json.dumps(
            {
                "exit_status": 0,
                "scientific_acceptance": acceptance,
                "qm_mm_launched": True,
            },
            indent=2,
        )
        + "\n"
    )
    LOGGER.info(
        "DFTB3/MM completed; accepted=%s tier=%s",
        acceptance["accepted"],
        acceptance["tier"],
    )


if __name__ == "__main__":
    main()
