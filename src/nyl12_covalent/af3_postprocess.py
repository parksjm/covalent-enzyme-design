"""Audit and export atom-stable AF3 covalent poses for Amber and inspection."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

from .chemistry import Substrate, product_report


@dataclass(frozen=True)
class AtomSite:
    """One heavy-atom coordinate and its source mmCIF identity."""

    serial: int
    element: str
    name: str
    residue: str
    chain: str
    sequence: int
    xyz: tuple[float, float, float]
    confidence: float


def read_atom_sites(path: Path) -> list[AtomSite]:
    """Read the AF3 atom-site loop while retaining author chain and atom names.

    Args:
        path: AlphaFold 3 model mmCIF.

    Returns:
        All heavy atoms in the sole model, in original CIF order.
    """
    lines = path.read_text().splitlines()
    try:
        start = lines.index("_atom_site.group_PDB")
    except ValueError as exc:
        raise ValueError(f"No atom-site loop in {path}") from exc
    fields: list[str] = []
    index = start
    while index < len(lines) and lines[index].startswith("_atom_site."):
        fields.append(lines[index].removeprefix("_atom_site."))
        index += 1
    required = {
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id",
        "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
        "B_iso_or_equiv", "pdbx_PDB_model_num", "label_alt_id",
    }
    if not required.issubset(fields):
        raise ValueError(f"AF3 atom-site loop lacks {sorted(required - set(fields))}")
    sites: list[AtomSite] = []
    seen: set[tuple[str, int, str]] = set()
    while index < len(lines) and lines[index] != "#":
        words = lines[index].split()
        index += 1
        if not words:
            continue
        if len(words) != len(fields):
            raise ValueError(f"Unexpected AF3 atom-site record: {words}")
        row = dict(zip(fields, words, strict=True))
        if row["pdbx_PDB_model_num"] != "1" or row["label_alt_id"] not in {".", "?"}:
            raise ValueError("Only one AF3 model without alternate locations is supported")
        if row["type_symbol"].upper() in {"H", "D"}:
            continue
        site = AtomSite(
            serial=int(row["id"]), element=row["type_symbol"],
            name=row["label_atom_id"], residue=row["label_comp_id"],
            chain=row["auth_asym_id"], sequence=int(row["auth_seq_id"]),
            xyz=tuple(float(row[f"Cartn_{axis}"]) for axis in "xyz"),
            confidence=float(row["B_iso_or_equiv"]),
        )
        key = (site.chain, site.sequence, site.name)
        if key in seen:
            raise ValueError(f"Duplicate AF3 atom identity {key}")
        seen.add(key)
        sites.append(site)
    if not sites:
        raise ValueError(f"No AF3 atom records in {path}")
    return sites


def _distance(first: AtomSite, second: AtomSite) -> float:
    """Return a Cartesian distance in Å."""
    return math.dist(first.xyz, second.xyz)


def _pdb_line(site: AtomSite, serial: int) -> str:
    """Convert one exact AF3 atom identity to an inspectable PDB record."""
    record = "HETATM" if site.chain == "I" else "ATOM  "
    residue = "LIG" if site.chain == "I" else site.residue
    if len(site.name) > 4 or len(residue) > 3 or serial > 99999:
        raise ValueError("AF3 atom cannot be represented losslessly in PDB")
    sequence = 1 if site.chain == "I" else site.sequence
    x, y, z = site.xyz
    return (
        f"{record}{serial:5d} {site.name:4s} {residue:3s} "
        f"{site.chain}{sequence:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.0:6.2f}{site.confidence:6.2f}          {site.element:>2s}  "
    )


def canonicalize_af3_ti1(
    cif: Path, substrate: Substrate, output: Path
) -> dict[str, object]:
    """Export an AF3 TI1 pose and an explicit atom/geometry audit.

    Args:
        cif: Completed AF3 covalent model.
        substrate: Named, product-gated full substrate graph.
        output: New canonical export directory.

    Returns:
        Audit including source hash, key distances, atom counts, and clashes.
    """
    product_report(substrate)
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Canonical output already exists: {output}")
    sites = read_atom_sites(cif)
    lookup = {(site.chain, site.sequence, site.name): site for site in sites}
    ligand = [site for site in sites if site.chain == "I"]
    protein = [site for site in sites if site.chain != "I"]
    names = {atom.GetProp("name") for atom in substrate.mol.GetAtoms() if atom.GetAtomicNum() > 1}
    if {site.name for site in ligand} != names:
        raise ValueError("AF3 ligand atom names do not match the substrate graph")
    expected = {**dict.fromkeys("ABCD", 232), **dict.fromkeys("EFGH", 96)}
    counts = {
        chain: len({site.sequence for site in protein if site.chain == chain})
        for chain in expected
    }
    if counts != expected or {site.chain for site in protein} != set(expected):
        raise ValueError(f"AF3 protein chain inventory changed: {counts}")
    carbon = lookup[("I", 1, substrate.carbonyl_c)]
    oxygen = lookup[("I", 1, substrate.carbonyl_o)]
    nitrogen = lookup[("I", 1, substrate.leaving_n)]
    og1 = lookup[("E", 1, "OG1")]
    beta_n = lookup[("E", 1, "N")]
    nearest: list[dict[str, object]] = []
    for lig in ligand:
        for prot in protein:
            if lig is carbon and prot is og1:
                continue
            distance = _distance(lig, prot)
            if distance < 2.4:
                nearest.append({
                    "distance_A": round(distance, 3),
                    "ligand": lig.name,
                    "protein": f"{prot.chain}{prot.sequence}:{prot.name}",
                })
    nearest.sort(key=lambda item: item["distance_A"])
    bond_outliers = []
    for bond in substrate.mol.GetBonds():
        first, second = bond.GetBeginAtom(), bond.GetEndAtom()
        if first.GetAtomicNum() == 1 or second.GetAtomicNum() == 1:
            continue
        pair = (first.GetProp("name"), second.GetProp("name"))
        distance = _distance(lookup[("I", 1, pair[0])], lookup[("I", 1, pair[1])])
        if distance < 1.1 or distance > 1.9:
            bond_outliers.append({"atoms": pair, "distance_A": round(distance, 3)})
    audit: dict[str, object] = {
        "source_cif": str(cif.resolve()),
        "source_sha256": hashlib.sha256(cif.read_bytes()).hexdigest(),
        "substrate": substrate.key,
        "product_gate": product_report(substrate),
        "atom_count": len(sites),
        "protein_chain_residue_counts": counts,
        "ligand_heavy_atom_count": len(ligand),
        "covalent_connection": ["E1:OG1", f"I1:{carbon.name}"],
        "distances_A": {
            "OG1_carbon": _distance(og1, carbon),
            "carbon_oxyanion": _distance(carbon, oxygen),
            "carbon_leaving_N": _distance(carbon, nitrogen),
            "beta_N_OG1": _distance(beta_n, og1),
        },
        "protein_ligand_contacts_below_2p4_A": nearest[:30],
        "ligand_bond_outliers_below_1p1_or_above_1p9_A": bond_outliers,
        "chemistry_limit": "AF3 output has no hydrogen/protonation or relay water; QM/MM must repair the reactive geometry.",
    }
    output.mkdir(parents=True)
    for filename, selection in (
        ("protein.pdb", protein), ("ligand.pdb", ligand),
        ("complex_heavy.pdb", sites),
    ):
        lines: list[str] = []
        previous_chain = ""
        for site in selection:
            if previous_chain and previous_chain != site.chain and previous_chain != "I":
                lines.append("TER")
            lines.append(_pdb_line(site, len(lines) + 1))
            previous_chain = site.chain
        output.joinpath(filename).write_text("\n".join([*lines, "END"]) + "\n")
    atom_map = [
        {"source_cif_serial": site.serial, "chain": site.chain,
         "residue": site.sequence, "source_atom_name": site.name,
         "canonical_atom_name": site.name}
        for site in sites
    ]
    (output / "atom_map.json").write_text(json.dumps(atom_map, indent=2) + "\n")
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit
