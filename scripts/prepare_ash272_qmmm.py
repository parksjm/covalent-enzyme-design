"""Prepare charge-audited, unrestrained ASH272 DFTB3/MM endpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np
import parmed
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
COMPARISON = ROOT / "qmmm/pa6_flipped_model30/06_ash272_comparison"
SOURCE_PREPARATION = (
    ROOT / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/preparation.json"
)
SOURCE_TOPOLOGY = (
    ROOT / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/inputs/system.prmtop"
)
ASH_TOPOLOGY = COMPARISON / "00_system/ash272_solvated.prmtop"
ASH_ACCEPTANCE = COMPARISON / "00_system/preparation_acceptance.json"
MM_STATUS = COMPARISON / "04_md_validation_20261113/mm_job_continuation/status.json"
SNAPSHOT_INDEX = (
    COMPARISON / "04_md_validation_20261113/results_attempt1/selected_snapshots.json"
)
MM_ROOT = COMPARISON / "05_selected_mm"
OUTPUT = COMPARISON / "06_qmmm_optimizations"
CASE_NAMES = (
    "replica_1_tier_A",
    "replica_1_tier_B",
    "replica_2_tier_A",
    "replica_2_tier_B",
)
QM_MASK = "@1621-1655,2727-2734,13905-13931,14555-14560,14574-14579,19914-19917,19925-19931,19937-19941"
PROTONATION_STATES = {968: "ASH", 1064: "ASH", 1160: "ASH", 1256: "ASH"}
LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _selected_indices(topology: Any, mask: str) -> list[int]:
    """Return selected zero-based atom indices from an Amber mask."""
    selection = AmberMask(topology, mask).Selection()
    indices = [index for index, selected in enumerate(selection) if selected]
    if not indices:
        raise ValueError(f"Amber mask selected no atoms: {mask}")
    return indices


def _endpoint_gate(
    topology_path: Path, restart: Path, acceptance: dict[str, Any]
) -> dict[str, Any]:
    """Apply existing ASH protonation and Tier A/B geometry gates."""
    source = COMPARISON / "04_md_validation_20261113/source"
    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(source))
    from ash272_geometry import build_indices, frame_geometry, validate_protonation
    from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

    topology = parmed.load_file(str(topology_path))
    protonation = validate_protonation(topology, acceptance)
    structure = parmed.load_file(str(topology_path), xyz=str(restart))
    coordinates = np.asarray(structure.coordinates, dtype=float)
    indices = build_indices(topology)
    metrics = frame_geometry(
        coordinates,
        indices,
        coordinates[indices["protein_ca"]],
        coordinates[indices["ligand_heavy"]],
    )
    tier = classify_geometry(
        {
            "d1": metrics["d1"],
            "d2": metrics["d2"],
            "d3": metrics["d3"],
            "attack": metrics["alpha_BD_deg"],
            "projected_phi": metrics["projected_phi_deg"],
        }
    )
    if tier is None or metrics["projected_phi_deg"] <= 0:
        raise ValueError(f"Endpoint failed Tier A/B or flipped-phi gate: {restart}")
    return {
        "accepted": True,
        "tier": tier,
        "projected_phi_deg": metrics["projected_phi_deg"],
        "d1_A": metrics["d1"],
        "d2_A": metrics["d2"],
        "d3_A": metrics["d3"],
        "attack_angle_deg": metrics["alpha_BD_deg"],
        "scissile_CO_A": metrics["scissile_CO_A"],
        "scissile_CN_A": metrics["scissile_CN_A"],
        "protonation": protonation,
        "gate": "existing modular Tier A/B and positive signed phi after MM minimization",
    }


def _input_hashes(path: Path) -> dict[str, str]:
    """Read and validate an endpoint's sha256sum manifest."""
    records = {}
    for line in path.read_text().splitlines():
        digest, separator, name = line.partition("  ")
        if not separator or len(digest) != 64:
            raise ValueError(f"Malformed endpoint hash record: {line!r}")
        records[str(Path(name).resolve())] = digest
    return records


def _mdin(qm_mask: str, active_mask: str, qm_charge: int) -> str:
    """Render a native Sander DFTB3/MM optimization with no restraints."""
    return (
        "ASH272 PA6 E-S charge-audited DFTB3/MM endpoint optimization\n"
        "&cntrl\n nmropt=0, imin=1, ntmin=5, ncyc=10, dx0=0.02, drms=0.05,\n"
        " maxcyc=2000, ntb=1, igb=0, ntc=1, ntf=1, cut=10.0,\n"
        " ifqnt=1, jfastw=4, ntr=0, ntpr=25, ntwx=25, ntwr=100, ioutfm=1, ntxo=2,\n/\n"
        "&qmmm\n"
        f" qmmask='{qm_mask}', qmcharge={qm_charge}, spin=1, qm_theory='DFTB3',\n"
        " qmmm_int=1, adjust_q=2, qm_ewald=1, qm_pme=1, qmcut=10.0,\n"
        " scfconv=1.0d-10, lnk_dis=1.09, lnk_atomic_no=1, lnk_method=1,\n/\n"
        "&dlfind\n"
        f" active='{active_mask}', optalg='LBFGS', crdrep='CART', tol=1.d-4,\n"
        " tole=1.d-5, maxcycle=2000, maxstep=0.5, verbosity=1,\n/\n"
    )


def prepare(output: Path = OUTPUT) -> dict[str, Any]:
    """Prepare four post-MM accepted inputs with protonation-derived QM charge."""
    from nyl12_covalent.ash272_qmmm import (
        audit_qm_region,
        compact_residue_mask,
        remap_active_indices,
    )

    if output.exists():
        raise FileExistsError(output)
    required = (
        SOURCE_PREPARATION,
        SOURCE_TOPOLOGY,
        ASH_TOPOLOGY,
        ASH_ACCEPTANCE,
        MM_STATUS,
        SNAPSHOT_INDEX,
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    old_prep = json.loads(SOURCE_PREPARATION.read_text())
    if json.loads(MM_STATUS.read_text()).get("exit_status") != 0:
        raise ValueError("The selected MM continuation did not complete successfully")
    selected = {
        f"{item['replica']}_tier_{item['tier']}": item
        for item in json.loads(SNAPSHOT_INDEX.read_text())
    }
    if set(selected) != set(CASE_NAMES):
        raise ValueError("Expected exactly the four declared selected MM cases")
    old_topology = parmed.load_file(str(SOURCE_TOPOLOGY))
    new_topology = parmed.load_file(str(ASH_TOPOLOGY))
    source_qm_indices = _selected_indices(old_topology, QM_MASK)
    source_region = audit_qm_region(
        old_topology,
        old_topology,
        source_qm_indices,
        {},
        reference_charge=int(old_prep["qmcharge"]),
        expected_link_h_count=int(old_prep["links"]),
    )
    region = audit_qm_region(
        old_topology,
        new_topology,
        source_qm_indices,
        PROTONATION_STATES,
        reference_charge=int(old_prep["qmcharge"]),
        expected_link_h_count=int(old_prep["links"]),
    )
    from nyl12_covalent.ash272_qmmm import ASP_STATE_CHARGE

    source_residues = {
        int(old_topology.atoms[index].residue.idx) + 1 for index in source_qm_indices
    }
    state_changes = {
        position: state
        for position, state in PROTONATION_STATES.items()
        if position in source_residues
        and str(old_topology.residues[position - 1].name).strip().upper() != state
    }
    expected_charge = source_region.formal_charge + sum(
        ASP_STATE_CHARGE[state]
        - ASP_STATE_CHARGE[
            str(old_topology.residues[position - 1].name).strip().upper()
        ]
        for position, state in state_changes.items()
    )
    if (
        region.real_atom_count != source_region.real_atom_count + len(state_changes)
        or region.link_h_count != source_region.link_h_count
        or region.formal_charge != expected_charge
        or region.nuclear_charge_with_links
        != source_region.nuclear_charge_with_links + len(state_changes)
        or region.electron_count != source_region.electron_count
        or region.multiplicity != source_region.multiplicity
    ):
        raise ValueError(
            "QM protonation, formal charge, link, or electron audit failed"
        )
    active_old = _selected_indices(old_topology, old_prep["active_mask"])
    active = remap_active_indices(
        old_topology, new_topology, active_old, PROTONATION_STATES
    )
    if not set(region.atom_indices) <= set(active):
        raise ValueError("QM atoms must lie inside the original movable 12 Å shell")
    active_mask = compact_residue_mask(new_topology, active)
    if _selected_indices(new_topology, active_mask) != sorted(active):
        raise ValueError("Compressed active mask changed atom membership")
    acceptance = json.loads(ASH_ACCEPTANCE.read_text())

    output.mkdir(parents=True)
    cases = []
    source_files = (*required, SOURCE_PREPARATION)
    for name in CASE_NAMES:
        item = selected[name]
        mm_case = MM_ROOT / name
        start = Path(item["restart"]["path"])
        if sha256(start) != item["restart"]["sha256"]:
            raise ValueError(f"Selected MD snapshot hash mismatch: {name}")
        manifest = _input_hashes(mm_case / "input_sha256.txt")
        if manifest.get(str(ASH_TOPOLOGY.resolve())) != sha256(ASH_TOPOLOGY):
            raise ValueError(f"MM topology provenance mismatch: {name}")
        if manifest.get(str(start.resolve())) != sha256(start):
            raise ValueError(f"MM starting snapshot provenance mismatch: {name}")
        endpoint, mdout = mm_case / "optimized.ncrst", mm_case / "minimized.mdout"
        if not endpoint.is_file() or "FINAL RESULTS" not in mdout.read_text(
            errors="replace"
        ):
            raise ValueError(f"MM optimization did not finish normally: {name}")
        geometry = _endpoint_gate(ASH_TOPOLOGY, endpoint, acceptance)
        case = output / name
        case.mkdir()
        starting = case / "starting.ncrst"
        shutil.copy2(endpoint, starting)
        mdin = case / "minimize.mdin"
        mdin.write_text(_mdin(region.mask, active_mask, region.formal_charge))
        sources = (
            *source_files,
            ASH_TOPOLOGY,
            endpoint,
            mdout,
            mm_case / "input_sha256.txt",
            starting,
            mdin,
        )
        record = {
            "case": name,
            "starting_mm_endpoint": str(endpoint.resolve()),
            "starting_endpoint_sha256": sha256(endpoint),
            "source_selected_frame": {
                "time_ns": item["time_ns"],
                "tier": item["tier"],
                "restart_sha256": item["restart"]["sha256"],
            },
            "mm_endpoint_geometry_acceptance": geometry,
            "qm_region": {
                "source_mask": QM_MASK,
                "remapped_mask": region.mask,
                "real_atom_count": region.real_atom_count,
                "link_h_count": region.link_h_count,
                "charge": region.formal_charge,
                "multiplicity": region.multiplicity,
                "nuclear_charge_with_links": region.nuclear_charge_with_links,
                "electron_count": region.electron_count,
                "protonation_states": PROTONATION_STATES,
                "active_atom_count": len(active),
                "active_mask": active_mask,
            },
            "method": "Amber Sander native DL-FIND / DFTB3-3OB-3-1 QM/MM",
            "restraints": "none; nmropt=0 and ntr=0",
            "force_convergence": "assessed from native Sander output after completion",
            "source_sha256": {str(path.resolve()): sha256(path) for path in sources},
            "status": "prepared; not launched",
        }
        (case / "preparation.json").write_text(json.dumps(record, indent=2) + "\n")
        (case / "run.sh").write_text(
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            ': "${SANDER_EXE:?set SANDER_EXE to the validated AmberTools binary}"\n'
            ': "${AMBERHOME:?set AMBERHOME to the validated AmberTools environment}"\n'
            ': "${LD_LIBRARY_PATH:?set the validated AmberTools library path}"\n'
            'exec "$SANDER_EXE" -O -i minimize.mdin -p '
            f'"{ASH_TOPOLOGY.resolve()}" -c starting.ncrst -o minimized.mdout '
            "-r optimized.ncrst -x optimization.nc -inf mdinfo\n"
        )
        (case / "run.sh").chmod(0o755)
        cases.append(record)
    report = {
        "qm_region": {
            "source_charge": int(old_prep["qmcharge"]),
            "declared_charge": region.formal_charge,
            "real_atom_count": region.real_atom_count,
            "link_h_count": region.link_h_count,
            "electron_count": region.electron_count,
            "multiplicity": region.multiplicity,
            "mask": region.mask,
            "protonation_states": PROTONATION_STATES,
            "state_logic": "formal ASP-to-ASH charge delta applies only inside the remapped QM region; catalytic Asp272 HD2 is selected",
        },
        "cases": cases,
        "status": "four inputs prepared; not launched",
    }
    (output / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info("Prepared %d ASH272 QM/MM endpoints in %s", len(cases), output)
    return report


def main() -> None:
    """Expose the preparation CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    prepare(args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
