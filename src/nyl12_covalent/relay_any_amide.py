"""Evaluate every PA6 amide with the established ES geometry conventions.

Tier thresholds reproduce analyze_nyl12_nma_four_site_geometry.py:tier.
They describe geometry, not catalytic activity. Input protein numbering is
explicitly alpha1-232/beta1-96; ligand maps come from actual predictor parsers.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import itertools
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from rdkit import Chem

from .active_site_geometry import angle_degrees, distance, projected_herschlag_phi

AtomKey = tuple[str, int, str]
Atoms = Mapping[AtomKey, tuple[str, Sequence[float]]]
TIERS = (
    ("A", 3.2, 3.8, 90.0, 120.0, 70.0, 110.0),
    ("B", 3.5, 4.0, 85.0, 125.0, 60.0, 125.0),
    ("C", 3.8, 4.2, 80.0, 130.0, 45.0, 140.0),
)


@dataclass(frozen=True)
class Amide:
    """Named coherent amide and its carbonyl-side carbon substituent."""

    carbon: str
    oxygen: str
    nitrogen: str
    adjacent_carbon: str


@dataclass(frozen=True)
class AmideGraph:
    """Heavy-element identities and all internal amides from one source graph."""

    elements: dict[str, str]
    amides: tuple[Amide, ...]
    source_sha256: str


def load_graph(path: Path) -> AmideGraph:
    """Read named amides from MOL2 connectivity, without perceiving spatial bonds.

    Args:
        path: Audited noncovalent PA6 graph; partial charges are not used.

    Returns:
        Heavy atoms and carbonyl-C/O/N/acyl-C groups, including every register.

    Raises:
        ValueError: The graph, atom names, or internal amides are invalid.
    """
    mol = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError(f"Invalid MOL2 graph: {path}")
    names = {a.GetIdx(): a.GetProp("_TriposAtomName") for a in mol.GetAtoms()}
    elements = {
        names[a.GetIdx()]: a.GetSymbol()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() != 1
    }
    if len(elements) != sum(a.GetAtomicNum() != 1 for a in mol.GetAtoms()):
        raise ValueError("Duplicate ligand heavy-atom names")
    amides = []
    for carbon in mol.GetAtoms():
        if carbon.GetAtomicNum() != 6:
            continue
        oxygen, nitrogen, adjacent = [], [], []
        for bond in carbon.GetBonds():
            other = bond.GetOtherAtom(carbon)
            key = (other.GetAtomicNum(), bond.GetBondTypeAsDouble())
            if key == (8, 2.0):
                oxygen.append(other.GetIdx())
            elif key == (7, 1.0):
                nitrogen.append(other.GetIdx())
            elif key == (6, 1.0):
                adjacent.append(other.GetIdx())
        if oxygen and nitrogen:
            if (len(oxygen), len(nitrogen), len(adjacent)) != (1, 1, 1):
                raise ValueError("Ambiguous carbonyl substituents")
            amides.append(
                Amide(
                    names[carbon.GetIdx()],
                    names[oxygen[0]],
                    names[nitrogen[0]],
                    names[adjacent[0]],
                )
            )
    if not amides:
        raise ValueError("No internal amides in source graph")
    return AmideGraph(
        elements, tuple(amides), hashlib.sha256(path.read_bytes()).hexdigest()
    )


def _point(atoms: Atoms, key: AtomKey, element: str) -> tuple[float, float, float]:
    if key not in atoms:
        raise ValueError(f"Missing mapped atom: {key}")
    actual, xyz = atoms[key]
    if actual.upper() != element.upper() or len(xyz) != 3:
        raise ValueError(f"Wrong element or coordinate shape: {key}")
    point = tuple(float(v) for v in xyz)
    if not all(math.isfinite(v) for v in point):
        raise ValueError(f"Nonfinite coordinate: {key}")
    return point


def _tier(row: dict[str, Any]) -> str:
    if not row["amide_geometry_valid"]:
        return "reject"
    for label, d1, oxy, amin, amax, pmin, pmax in TIERS:
        if (
            row["d1"] <= d1
            and row["d2"] <= oxy
            and row["d3"] <= oxy
            and amin <= row["a_bd"] <= amax
            and pmin <= abs(row["herschlag_phi"]) <= pmax
        ):
            return label
    return "reject"


def _measure(site: dict[str, tuple], ligand: dict[str, tuple], amide: Amide) -> dict:
    c, o, n, adj = (
        ligand[x]
        for x in (amide.carbon, amide.oxygen, amide.nitrogen, amide.adjacent_carbon)
    )
    og = site["og"]
    co, cn, cc = distance(c, o), distance(c, n), distance(c, adj)
    row = {
        "carbonyl_c": amide.carbon,
        "carbonyl_o": amide.oxygen,
        "leaving_n": amide.nitrogen,
        "acyl_adjacent_c": amide.adjacent_carbon,
        "r_c_o": co,
        "r_c_n": cn,
        "r_c_adjacent": cc,
        "amide_geometry_valid": 1.0 < co < 1.6 and 1.1 < cn < 1.7 and 1.2 < cc < 1.8,
        "d1": distance(og, c),
        "d2": distance(o, site["phe"]),
        "d3": distance(o, site["asn"]),
        "a_bd": angle_degrees(og, c, o),
        "herschlag_phi": projected_herschlag_phi(og, o, c, n, adj),
        "nterm_to_leaving_n": distance(site["nterm"], n),
        "og1_to_leaving_n": distance(og, n),
        "cb_og1_carbonyl_angle": angle_degrees(site["cb"], og, c),
    }
    for residue in (272, 274):
        row[f"asp{residue}_to_leaving_n"] = min(
            distance(site[f"asp{residue}_{i}"], n) for i in (1, 2)
        )
    row["asp_pair_min_oo"] = min(
        distance(site[f"asp272_{i}"], site[f"asp274_{j}"])
        for i in (1, 2)
        for j in (1, 2)
    )
    row["tier"] = _tier(row)
    return row


def evaluate_model(
    atoms: Atoms, ligand_maps: Mapping[str, Mapping[str, str]], graph: AmideGraph
) -> dict[str, Any]:
    """Evaluate every amide/site and assign four ligands without double counting.

    Args:
        atoms: Heavy atoms keyed by actual chain, beta-local residue, atom name.
        ligand_maps: I-L mappings from actual model atom names to MOL2 names.
        graph: Audited source ligand graph from :func:`load_graph`.

    Returns:
        All 48 local evaluations, four assigned records and explicit limitations.

    Raises:
        ValueError: Mapping, inventory, coordinates, or geometric planes fail.
    """
    if set(ligand_maps) != set("IJKL") or len(graph.amides) != 3:
        raise ValueError("Expected four L4 ligands and three internal amides each")
    sites = {}
    for alpha, beta in zip("ABCD", "EFGH", strict=True):
        specs = {
            "og": (beta, 1, "OG1", "O"),
            "nterm": (beta, 1, "N", "N"),
            "cb": (beta, 1, "CB", "C"),
            "phe": (alpha, 113, "N", "N"),
            "asn": (alpha, 186, "ND2", "N"),
        }
        for residue in (272, 274):
            for i in (1, 2):
                specs[f"asp{residue}_{i}"] = (beta, residue - 232, f"OD{i}", "O")
        sites[alpha + beta] = {k: _point(atoms, v[:3], v[3]) for k, v in specs.items()}
    rows = []
    for chain, mapping in ligand_maps.items():
        keys = [k for k in atoms if k[0] == chain]
        if (
            len({k[1] for k in keys}) != 1
            or len(keys) != len(graph.elements)
            or set(mapping) != {k[2] for k in keys}
            or len(set(mapping.values())) != len(mapping)
            or set(mapping.values()) != set(graph.elements)
        ):
            raise ValueError(f"Nonbijective or incomplete ligand atom map: {chain}")
        ligand = {
            mapping[k[2]]: _point(atoms, k, graph.elements[mapping[k[2]]]) for k in keys
        }
        for name, site in sites.items():
            for amide in graph.amides:
                rows.append(
                    {"ligand": chain, "site": name, **_measure(site, ligand, amide)}
                )
    ranks = {"A": 0, "B": 1, "C": 2, "reject": 3}
    best = {
        (lig, site): min(
            (r for r in rows if r["ligand"] == lig and r["site"] == site),
            key=lambda r: (ranks[r["tier"]], r["d1"], r["carbonyl_c"]),
        )
        for lig in "IJKL"
        for site in sites
    }
    choices = []
    for permutation in itertools.permutations(sites):
        picked = [
            best[lig, site] for lig, site in zip("IJKL", permutation, strict=True)
        ]
        key = (
            -sum(r["tier"] in {"A", "B"} for r in picked),
            sum(ranks[r["tier"]] for r in picked),
            sum(r["d1"] for r in picked),
            permutation,
        )
        choices.append((key, picked))
    _, assigned = min(choices, key=lambda item: item[0])
    return {
        "all_amide_site_rows": rows,
        "assigned": assigned,
        "assigned_tier_ab_count": sum(r["tier"] in {"A", "B"} for r in assigned),
        "assignment_policy": "max_AB_then_tier_then_total_d1_one_to_one",
        "source_graph_sha256": graph.source_sha256,
        "tier_thresholds": TIERS,
        "native_thr_hard_gate": False,
        "interpretation": "Unrelaxed ES geometry proxies; not activity or independent sites",
    }
