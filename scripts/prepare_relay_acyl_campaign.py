#!/usr/bin/env python3
"""Prepare one balanced acyl-enzyme hypothesis from the accepted E-S model."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import re
import shutil
import sys

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from nyl12_covalent.acylenzyme_seed import AcylAtomMap  # noqa: E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa: E402
from nyl12_covalent.relay_acyl_seed import (  # noqa: E402
    RelayAcylAtomMap,
)
from nyl12_covalent.relay_case import load_case, write_json  # noqa: E402
from nyl12_covalent.relay_seed_search import search_balanced_acyl_seed  # noqa: E402
from nyl12_covalent.relay_seed_refinement import refine_acyl_torsion_seed  # noqa: E402
from scripts.prepare_dimer_continuation import write_restart  # noqa: E402

LOGGER = logging.getLogger(__name__)
BASE = ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison"
CAMPAIGN = BASE / "09_water_relay_acylation"
SOURCE = BASE / "06_qmmm_optimization/replica_1_tier_A_seed/water_qm_reoptimization"
ATOM_MAP = {
    "thr_n": 13905,
    "thr_h1": 13906,
    "thr_h2": 13907,
    "thr_cb": 13910,
    "thr_og": 13916,
    "proton": 13917,
    "carbon": 19920,
    "oxygen": 19921,
    "leaving_n": 19933,
    "leaving_h": 19941,
    "water_o": 20317,
    "water_h1": 20318,
    "water_h2": 20319,
    "ash_oxygen": 14579,
    "ash_proton": 14580,
    "asp_o1": 14559,
    "asp_o2": 14560,
    "phe_n": 1638,
    "asn_nd2": 2732,
}


def prepare() -> Path:
    """Retain immutable sources and construct an explicitly unoptimized product.

    Returns:
        Path to the validated campaign case; no engine is launched.

    Raises:
        ValueError: If the accepted E-S state, masks or construction gates fail.
    """
    accepted = json.loads((SOURCE / "acceptance.json").read_text())
    if accepted["status"] != "accepted" or not accepted["acceptance"]["accepted"]:
        raise ValueError("Original water-QM E-S source is not accepted")
    folder = CAMPAIGN / "01_endpoint_inputs"
    folder.mkdir(exist_ok=False)
    topology_path = BASE / "00_system/ash274_solvated.prmtop"
    topology = parmed.load_file(str(topology_path))
    source = SOURCE / "optimized.ncrst"
    restart = Rst7.open(str(source))
    xyz = np.asarray(restart.coordinates).reshape(-1, 3).copy()
    text = (SOURCE / "minimize.mdin").read_text()
    qm_mask = re.findall(r"\bqmmask='([^']+)'", text)
    active_mask = re.findall(r"\bactive='([^']+)'", text)
    if len(qm_mask) != 1 or len(active_mask) != 1:
        raise ValueError("Expected exactly one source QM and active mask")
    qm = list(AmberMask(topology, qm_mask[0]).Selected())
    active = list(AmberMask(topology, active_mask[0]).Selected())
    prefix, separator, _ = text.partition("&dlfind")
    if not separator or "dftb_telec" in prefix.lower():
        raise ValueError("Unexpected source optimizer/electronic-temperature input")
    prefix, count = re.subn(r"\bmaxcyc=2000\b", "maxcyc=1", prefix)
    if count != 1:
        raise ValueError("Unexpected source maximum cycle count")
    prefix = prefix.replace("qm_theory='DFTB3',", "qm_theory='DFTB3', dftb_telec=0.0,")
    mdin = prefix + (
        "&dlfind\n optalg='LBFGS', crdrep='CART', dimer=.false.,\n"
        f" active='{active_mask[0]}', maxcycle=1, tol=1.d-4, tole=1.d-5,\n"
        " maxstep=0.1, verbosity=2,\n/\n"
    )
    shutil.copy2(source, folder / "accepted_es.ncrst")
    cfg = {
        "topology": str(topology_path),
        "source_restart": str(folder / "accepted_es.ncrst"),
        "mdin": mdin,
        "qm_atom_ids": [i + 1 for i in qm],
        "active_atom_ids": [i + 1 for i in active],
        "atom_map": ATOM_MAP,
        "qm_charge": -1,
        "amber_home": "/data/jpc/code/ambertools26/venv",
        "calibration": str(
            ROOT
            / (
                "qmmm/pa6_flipped_model30/05_neutral_reaction_path/"
                "07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113/"
                "native_gradient_0K_onecycle_20261113_attempt2/scientific_acceptance.json"
            )
        ),
        "max_atom_step_A": 0.05,
        "force_gate": 0.12,
        "memory": 7,
        "microstate": "All four beta N termini neutral; ASP272/ASH274; same proton inventory throughout",
        "movable_region": "Unchanged accepted complete-residue 12 A selection; exterior frozen; ntr=0",
        "independent_endpoint_test": True,
        "relay_ready_ES_demonstrated": False,
        "no_MM_product_parameters": True,
        "construction_heavy_movement_cap_A": 3.0,
        "construction_H_movement_cap_A": 5.0,
    }
    case_path = CAMPAIGN / "case.json"
    write_json(case_path, cfg)
    model = load_case(case_path)
    if int(model.active.sum()) != 4662:
        raise ValueError("Accepted 12 A active shell changed")
    acyl = AcylAtomMap(
        **{key: ATOM_MAP[key] for key in AcylAtomMap.__dataclass_fields__}
    )
    mapping = RelayAcylAtomMap(
        acyl,
        ATOM_MAP["water_o"],
        ATOM_MAP["water_h1"],
        ATOM_MAP["water_h2"],
        ATOM_MAP["ash_oxygen"],
        ATOM_MAP["ash_proton"],
        (ATOM_MAP["asp_o1"], ATOM_MAP["asp_o2"]),
    )
    base_product, base_evidence = search_balanced_acyl_seed(
        topology,
        xyz,
        mapping,
        active_atom_ids=tuple(cfg["active_atom_ids"]),
        qm_charge=-1,
    )
    product, evidence = refine_acyl_torsion_seed(
        topology,
        xyz,
        mapping,
        base_product,
        base_evidence,
        active_atom_ids=tuple(cfg["active_atom_ids"]),
        max_heavy_displacement_A=cfg["construction_heavy_movement_cap_A"],
        max_hydrogen_displacement_A=cfg["construction_H_movement_cap_A"],
    )
    write_json(folder / "search_summary.json", evidence["search"])
    if evidence["severe_contact_count"]:
        raise ValueError("Bounded seed construction retained severe nonbonded contacts")
    # This audit reports contacts rather than silently deleting/moving waters.
    write_json(folder / "construction.json", evidence)
    write_restart(folder / "acyl_guess.ncrst", product, restart)
    paths = [
        Path(__file__),
        SOURCE / "acceptance.json",
        SOURCE / "minimize.mdin",
        source,
        topology_path,
        case_path,
        ROOT / "src/nyl12_covalent/relay_acyl_seed.py",
        ROOT / "src/nyl12_covalent/acylenzyme_seed.py",
        ROOT / "src/nyl12_covalent/relay_seed_search.py",
        ROOT / "src/nyl12_covalent/relay_seed_refinement.py",
        ROOT / "src/nyl12_covalent/relay_case.py",
        ROOT / "src/nyl12_covalent/native_sander_cli.py",
        ROOT / "src/nyl12_covalent/proton_coordination.py",
        ROOT / "scripts/prepare_dimer_continuation.py",
    ]
    snapshot = folder / "source_snapshot"
    snapshot.mkdir()
    hashes = {}
    for index, path in enumerate(paths):
        hashes[str(path)] = source_sha256(path)
        copied = snapshot / f"{index:02d}_{path.name}"
        shutil.copy2(path, copied)
        if source_sha256(copied) != hashes[str(path)]:
            raise ValueError("Preparation source snapshot hash mismatch")
    write_json(
        folder / "manifest.json",
        {
            "input_sha256": hashes,
            "original_atoms": len(xyz),
            "qm_atoms": len(qm),
            "active_atoms": len(active),
            "link_atoms": 8,
            "seed_recipe": "285 leaving placements then 24 acyl-arm torsions; retain acute source ester angle",
            "seed_ester_angle_degrees": evidence["ester_angle_degrees"],
            "seed_is_optimized": False,
            "seed_severe_contact_count": evidence["severe_contact_count"],
            "source_is_already_MM_then_QMMM_optimized": True,
            "new_MM_minimization": False,
            "reason_no_product_MM": "All reacting atoms are QM; E-S MM connectivity is not a product force field",
            "relay_hypothesis_not_proven_by_endpoint_construction": True,
        },
    )
    LOGGER.info(
        "Prepared balanced product hypothesis; %d contacts require inspection",
        len(evidence["close_contacts"]),
    )
    return case_path


def main() -> None:
    """Prepare the fixed accepted-source endpoint test without launching an engine."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    LOGGER.info("Prepared case: %s", prepare())


if __name__ == "__main__":
    main()
