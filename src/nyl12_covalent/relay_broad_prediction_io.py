"""Software-level acceptance for unguided eight-chain E-S predictions."""

from __future__ import annotations

import math
import shlex

from rdkit import Chem
from pathlib import Path

from .af3_postprocess import read_atom_sites
from .ligandmpnn_pilot_audit import THREE_TO_ONE


def audit_prediction(
    cif: Path, alpha: str, beta: str, atom_map: dict[str, str]
) -> dict:
    """Check mature-chain sequences, physical ligand atoms and finite coordinates.

    Args:
        cif: Sole native prediction's model mmCIF.
        alpha: Expected232-residue alpha sequence, repeatedA–D.
        beta: Expected96-residue beta sequence, repeatedE–H.
        atom_map: Actual predictor atomname to sourceMol2 name dictionary.

    Returns:
        Software/graph acceptance only; no contact or WT geometry gate.
    """
    sites = read_atom_sites(cif)
    if {s.chain for s in sites} != set("ABCDEFGHIJKL"):
        raise ValueError("Expectedexactly8matureproteinand4L4ligandchains")
    if not all(math.isfinite(v) for s in sites for v in s.xyz):
        raise ValueError("Nonfinitepredictedphysicalcoordinates")
    for chain in "ABCDEFGH":
        ca = sorted(
            (s for s in sites if s.chain == chain and s.name == "CA"),
            key=lambda s: s.sequence,
        )
        expected = alpha if chain in "ABCD" else beta
        if (
            len(ca) != len(expected)
            or "".join(THREE_TO_ONE[s.residue] for s in ca) != expected
        ):
            raise ValueError("Predictedmaturechainsequence/inventorymismatch:" + chain)
        if [s.sequence for s in ca] != list(range(1, len(expected) + 1)):
            raise ValueError("Predictedmaturechainnumberingunexpected:" + chain)
        molecule = Chem.MolFromSequence(expected)
        required = {}
        for atom in molecule.GetAtoms():
            info = atom.GetPDBResidueInfo()
            name = info.GetName().strip()
            if name != "OXT":
                required.setdefault(info.GetResidueNumber(), set()).add(name)
        for number, names in required.items():
            observed = {
                s.name for s in sites if s.chain == chain and s.sequence == number
            }
            if observed - {"OXT"} != names:
                raise ValueError("Predictedproteinheavyatominventorymismatch:" + chain)
    for chain in "IJKL":
        ligand = [s for s in sites if s.chain == chain]
        if len(ligand) != 33 or {s.name for s in ligand} != set(atom_map):
            raise ValueError("PredictedL4atomnametablechanged:" + chain)
        if any(s.element.upper() != atom_map[s.name][0] for s in ligand):
            raise ValueError("Predictedligandelementsdifferentfromsourcegraph")
        if len({s.sequence for s in ligand}) != 1:
            raise ValueError("ExpectedonesingleCCD/SMILESligandperchain")
    lines = cif.read_text().splitlines()
    key = "_struct_conn.conn_type_id"
    if key in lines:
        start = lines.index(key)
        while start and lines[start - 1].startswith("_struct_conn."):
            start -= 1
        columns = []
        index = start
        while index < len(lines) and lines[index].startswith("_struct_conn."):
            columns.append(lines[index])
            index += 1
        for line in lines[index:]:
            if line == "#" or line.startswith("loop_"):
                break
            if not line.strip():
                continue
            row = dict(zip(columns, shlex.split(line), strict=True))
            a, b = (
                row["_struct_conn.ptnr1_label_asym_id"],
                row["_struct_conn.ptnr2_label_asym_id"],
            )
            if row[key].lower().startswith("coval") and (
                (a in "IJKL") != (b in "IJKL")
            ):
                raise ValueError("Unexpectedcovalentproteinligandconnection")
    return dict(
        software_graph_acceptance=True,
        protein_chains=8,
        ligand_chains=4,
        heavy_ligand_atoms_per_chain=33,
        finite_coordinates=True,
        ligand_atom_maps={c: atom_map for c in "IJKL"},
        geometric_utility="Not evaluated; separate any-amide metrics; no WT-contact gate",
    )


def export_geometry(cif: Path, ligand_maps: dict, output: Path) -> dict:
    """Export lead-owned any-amide metrics separately from software acceptance.

    Args:
        cif: Accepted native predictor mmCIF with explicit mature-chain identities.
        ligand_maps: Per-I–L actual-model-name to source-Mol2-name dictionaries.
        output: Existing single-model directory for geometry JSON and CSV tables.

    Returns:
        Paths for all48 comparisons and the four one-to-one assigned records.
    """
    import csv
    import json

    from .relay_any_amide import evaluate_model, load_graph
    from .relay_broad_graph import CONFIG

    atoms = {
        (s.chain, s.sequence, s.name): (s.element, s.xyz) for s in read_atom_sites(cif)
    }
    report = evaluate_model(
        atoms, ligand_maps, load_graph(CONFIG / "whole_l4_source.mol2")
    )
    if len(report["all_amide_site_rows"]) != 48 or len(report["assigned"]) != 4:
        raise ValueError("Expected48amide/sitecomparisonsand4assignedpairs")
    report["software_acceptance_separate"] = True
    json_path = output / "geometry_any_amide.json"
    json_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    paths = {}
    for key, name in [
        ("all_amide_site_rows", "geometry_all_48.csv"),
        ("assigned", "geometry_assigned_4.csv"),
    ]:
        path = output / name
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(report[key][0]))
            writer.writeheader()
            writer.writerows(report[key])
        paths[key] = str(path)
    return dict(
        geometry_json=str(json_path),
        geometry_csv=paths,
        assigned_tier_ab_count=report["assigned_tier_ab_count"],
        interpretation="SeparateunrelaxedESgeometryproxies;notactivity;noWTcontactgate",
    )
