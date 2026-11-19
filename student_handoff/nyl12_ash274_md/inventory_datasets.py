"""Existing WTdataset chemistry inventory and MD-onlyreceipt normalization."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from parmed.amber import AmberFormat, Rst7

from inventory_common import BASE, OUT, add, sha, write


def topology_audit(path: Path, microstate: str) -> tuple[dict, list[dict]]:
    """Measure real24protected identities, allfourneutralTHR andAsp/ASHcharges."""
    data = AmberFormat(str(path)).parm_data
    names, labels = data["ATOM_NAME"], data["RESIDUE_LABEL"]
    ptr = [p - 1 for p in data["RESIDUE_POINTER"]] + [data["POINTERS"][0]]
    charges = data["CHARGE"]
    records, states = [], []
    acid272, acid274 = microstate.split("/")
    for site, alpha, beta in zip(("AE", "BF", "CG", "DH"), "ABCD", "EFGH", strict=True):
        n = ord(alpha) - ord("A")
        for chain, seq, structure, amber, expected in (
            (alpha, 113, 113, 113 + 232 * n, "PHE"),
            (alpha, 156, 156, 156 + 232 * n, "LYS"),
            (alpha, 186, 186, 186 + 232 * n, "ASN"),
            (beta, 233, 1, 929 + 96 * n, "THR"),
            (beta, 272, 40, 968 + 96 * n, acid272[:3]),
            (beta, 274, 42, 970 + 96 * n, acid274[:3]),
        ):
            # Canonical microstate includes residue number suffix: ASP272/ASH274.
            expected = "ASH" if expected == "ASH" else expected
            if labels[amber - 1] != expected:
                raise ValueError(f"Protected residue mismatch {chain}/{seq}/{amber}")
            records.append(
                dict(
                    site=site,
                    chain=chain,
                    sequence_residue=seq,
                    structure_residue=structure,
                    amber_residue=amber,
                    resname=labels[amber - 1],
                )
            )
        for amber in (929 + 96 * n, 968 + 96 * n, 970 + 96 * n):
            first, last = ptr[amber - 1], ptr[amber]
            atoms = dict(zip(names[first:last], charges[first:last], strict=True))
            if amber == 929 + 96 * n:
                if (
                    "H3" in atoms
                    or not all(k in atoms for k in ("H1", "H2"))
                    or any(
                        abs(atoms[k] - q) > 1e-6
                        for k, q in (("N", -0.4844), ("H1", 0.1224), ("H2", 0.1224))
                    )
                ):
                    raise ValueError(
                        "Neutral beta N-terminal THR inventory/charges changed"
                    )
            else:
                expected_charge = 0.0 if labels[amber - 1] == "ASH" else -1.0
                if ("HD2" in atoms) != (expected_charge == 0.0) or abs(
                    sum(atoms.values()) - expected_charge
                ) > 1e-4:
                    raise ValueError("Asp/ASH charge/proton inventory changed")
            states.append(
                dict(
                    amber_residue=amber,
                    resname=labels[amber - 1],
                    charge_e=sum(atoms.values()),
                    atom_charges={
                        k: v
                        for k, v in atoms.items()
                        if k in ("N", "H1", "H2", "H3", "OD1", "OD2", "HD2")
                    },
                )
            )
    return dict(
        natom=data["POINTERS"][0],
        total_charge_e=sum(charges),
        four_neutral_beta_Thr=True,
        protected_site_charge_inventory=states,
    ), records


def dataset(
    identifier: str,
    root: Path,
    topology_name: str,
    microstate: str,
    replicas: list[str],
) -> dict:
    """Pin an existing unrestrained minimizedWTsystem, receipts andMDprotocols."""
    system = root / "00_system"
    topology = system / topology_name
    neutral = identifier == "wt_asp272_asp274"
    restart = (
        BASE / "01_mm_minimization/min.rst7"
        if neutral
        else system / "min_unrestrained.rst7"
    )
    prep = system / ("acceptance.json" if neutral else "preparation_acceptance.json")
    receipt = (
        BASE / "01_mm_minimization/md_start_acceptance.json"
        if neutral
        else system / "minimization_acceptance.json"
    )
    destination = f"datasets/{identifier}"
    top_art = add(topology, f"{destination}/topology.prmtop", "matched_WT_topology")
    rst_art = add(
        restart, f"{destination}/start.rst7", "unrestrained_minimized_WT_restart"
    )
    audit, protected = topology_audit(topology, microstate)
    rst = Rst7.open(str(restart))
    xyz = np.asarray(rst.coordinates).reshape(-1, 3)
    if (
        xyz.shape != (audit["natom"], 3)
        or not np.isfinite(xyz).all()
        or rst.box is None
        or not np.isfinite(rst.box).all()
    ):
        raise ValueError("Original minimized restart is unreadable/invalid")
    source_prep = add(
        prep,
        f"{destination}/evidence/source_preparation_acceptance.json",
        "source_acceptance",
    )
    source_min = add(
        receipt,
        f"{destination}/evidence/source_md_start_acceptance.json",
        "source_acceptance",
    )
    original = json.loads(receipt.read_text())
    if neutral:
        process_path = BASE / "01_mm_minimization/process_acceptance.json"
        process_complete = json.loads(process_path.read_text())["returncode"] == 0
        if not original["md_start_contact_gate_pass"]:
            raise ValueError("Baseline MDstart contactgate failed")
        add(
            process_path,
            f"{destination}/evidence/source_process_acceptance.json",
            "source_acceptance",
        )
    else:
        process_complete = original["process_complete"] is True
        if (
            original["restart_sha256"] != rst_art["sha256"]
            or original["topology_sha256"] != top_art["sha256"]
        ):
            raise ValueError("Original receipt hashes differ")
    if not process_complete:
        raise ValueError("Source minimization process not complete")
    local = OUT / "datasets" / identifier
    prep_wrapper = dict(
        accepted=True,
        scope="Existing topology/microstate inventory only; no new scientific acceptance",
        topology_sha256=top_art["sha256"],
        source_receipt=source_prep,
        source_receipt_absolute=str(prep),
        source_receipt_sha256=sha(prep),
        topology_audit=audit,
    )
    min_wrapper = dict(
        process_complete=True,
        restart_readable=True,
        coordinates_sha256=rst_art["sha256"],
        restart_sha256=rst_art["sha256"],
        topology_sha256=top_art["sha256"],
        restart_stage="unrestrained_minimized",
        force_convergence="not asserted; MD preparation only",
        stationary_point_claim=False,
        source_receipt=source_min,
        source_receipt_absolute=str(receipt),
        source_receipt_sha256=sha(receipt),
        source_assessment=original,
    )
    write(local / "prep_acceptance.json", prep_wrapper)
    write(local / "md_start_receipt.json", min_wrapper)
    pa = add(
        local / "prep_acceptance.json",
        f"{destination}/prep_acceptance.json",
        "normalized_preparation_receipt",
    )
    ma = add(
        local / "md_start_receipt.json",
        f"{destination}/md_start_receipt.json",
        "normalized_MD_start_receipt",
    )
    for name in [
        "input_complex_leap.pdb",
        "md_dry.pdb",
        "leap_final.in"
        if neutral
        else f"leap_ash{272 if 'ash272' in identifier else 274}_final.in",
    ]:
        add(system / name, f"{destination}/build/{name}", "system_rebuild_input")
    for name in [
        "preparation_manifest.json",
        "preparation_acceptance.json",
        "min_unrestrained.in",
    ]:
        if (system / name).exists():
            add(
                system / name,
                f"{destination}/build/{name}",
                "system_rebuild_provenance",
            )
    protocols = []
    for replica in replicas:
        d = root / "03_md" / replica
        protocol = json.loads((d / "protocol.json").read_text())
        files = []
        for p in sorted(d.iterdir()):
            if p.is_file() and (
                p.suffix == ".in"
                or p.name
                in ["protocol.json", "extension_100ns.json", "completion.json"]
                or p.name.endswith(".acceptance.json")
            ):
                files.append(
                    add(
                        p,
                        f"{destination}/source_md/{replica}/{p.name}",
                        "source_MD_protocol",
                    )
                )
        for name in ("04_equil_npt.rst7", "06_production_90ns.rst7"):
            add(
                d / name,
                f"{destination}/source_md/{replica}/{name}",
                "existing_equilibrated_or_production_restart",
                required=False,
            )
        protocols.append(
            dict(
                replica=replica,
                protocol={k: v for k, v in protocol.items() if k != "input_sha256"},
                files=files,
            )
        )
    return dict(
        microstate=microstate,
        topology=top_art,
        restart=rst_art,
        preparation_acceptance=pa,
        minimization_acceptance=ma,
        restart_stage="unrestrained_minimized",
        minimizedstart=rst_art,
        protected_residues=protected,
        source_md_protocols=protocols,
        scientifically_optional_control=identifier == "wt_ash272_asp274",
        requires_new_system_build_for_mutants=True,
    )
