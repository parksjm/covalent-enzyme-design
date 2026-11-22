"""Screen static PA6 L4 central amides by signed phi and native Thr geometry."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LEGACY = Path("/data/jpc/enzyme-substrate-modeling")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(LEGACY / "scripts"))
from analyze_nyl12_nma_four_site_geometry import tier  # noqa: E402
from nyl12_covalent.active_site_geometry import (  # noqa: E402
    angle_degrees,
    distance,
    projected_herschlag_phi,
)
from nyl12_covalent.catalytic_alignment import (  # noqa: E402
    BACKBONE,
    FIT_RESIDUES,
    SITE_RESIDUES,
    _cif_atoms,
)
from nyl12_covalent.catalytic_alignment_geometry import (  # noqa: E402
    fit_coordinates,
    state_metrics,
)


def read_atoms(path: Path) -> dict:
    """Read unique heavy-atom PDB keys, rejecting ambiguous coordinates."""
    atoms = {}
    for line in path.read_text().splitlines():
        if line[:6] not in {"ATOM  ", "HETATM"}:
            continue
        element = line[76:78].strip().upper() or line[12:16].strip()[0]
        if element == "H":
            continue
        key = (line[21], int(line[22:26]), line[12:16].strip())
        if key in atoms or line[16] not in (" ", "A"):
            raise ValueError(f"Ambiguous atom {key}: {path}")
        xyz = np.array([float(line[i : i + 8]) for i in (30, 38, 46)])
        if not np.isfinite(xyz).all():
            raise ValueError(f"Nonfinite atom {key}: {path}")
        atoms[key] = (element, xyz)
    return atoms


def layout_site(atoms: dict, site: int) -> dict:
    """Map four symmetry-related sites to canonical biological residues."""
    result = {}
    if ("A", 929, "OG1") in atoms:
        layout = "flat_amber"
    elif ("E", 1, "OG1") in atoms:
        layout = "split"
    elif ("A", 233, "OG1") in atoms:
        layout = "canonical"
    else:
        raise ValueError("Unrecognized protein layout")
    for _, auth, _ in SITE_RESIDUES:
        if layout == "flat_amber":
            key_chain = "A"
            res = auth + 232 * site if auth <= 232 else 929 + 96 * site + auth - 233
        elif layout == "split":
            key_chain = "ABCD"[site] if auth <= 232 else "EFGH"[site]
            res = auth if auth <= 232 else auth - 232
        else:
            key_chain, res = "ABCD"[site], auth
        for (chain, number, name), value in atoms.items():
            if chain == key_chain and number == res:
                result[(auth, name)] = value
    return result


def measure(site: dict, ligand: dict, reference: dict) -> dict[str, Any]:
    """Measure signed orientation, catalytic geometry, and native Thr gates."""

    def pos(residue: int, atom: str) -> np.ndarray:
        return site[(residue, atom)][1]

    c, o, n, adj = [ligand[k] for k in ("C54", "O49", "N59", "C62")]
    if not (
        1.0 < distance(c, o) < 1.6
        and 1.1 < distance(c, n) < 1.7
        and 1.2 < distance(c, adj) < 1.8
    ):
        raise ValueError("Central amide connectivity inconsistent with atom mapping")
    og = pos(233, "OG1")
    row = {
        "d1": distance(og, c),
        "d2": distance(o, pos(113, "N")),
        "d3": distance(o, pos(186, "ND2")),
        "a_bd": angle_degrees(og, c, o),
        "herschlag_phi": projected_herschlag_phi(og, o, c, n, adj),
        "nterm_to_leaving_n": distance(pos(233, "N"), n),
        "og1_to_leaving_n": distance(og, n),
        "cb_og1_carbonyl_angle": angle_degrees(pos(233, "CB"), og, c),
    }
    for residue in (272, 274):
        row[f"asp{residue}_to_leaving_n"] = min(
            distance(pos(residue, a), n) for a in ("OD1", "OD2")
        )
    row["tier"] = tier(row)
    row["orientation"] = "positive_phi" if row["herschlag_phi"] > 0 else "negative_phi"
    keys = [(auth, a) for _, auth, _ in FIT_RESIDUES for a in BACKBONE]
    rotation, center, target, rmsd = fit_coordinates(
        np.array([site[k][1] for k in keys]), np.array([reference[k][1] for k in keys])
    )
    aligned = {k: (e, (v - center) @ rotation + target) for k, (e, v) in site.items()}
    native = state_metrics(reference, aligned)["thr233_held_out"]
    row.update(
        scaffold_rmsd=rmsd,
        og1_shift=native["heavy_atom_displacements_angstrom"]["OG1"],
        nterm_shift=native["heavy_atom_displacements_angstrom"]["N"],
        psi_delta=native["psi_delta_es_minus_reference_degrees"],
        chi1_delta=native["chi1_delta_es_minus_reference_degrees"],
    )
    row["native_pass"] = (
        row["og1_shift"] <= 0.75
        and row["nterm_shift"] <= 1.25
        and abs(row["psi_delta"]) <= 35
        and abs(row["chi1_delta"]) <= 20
    )
    row["joint_pass"] = row["native_pass"] and row["tier"] in ("A", "B")
    return row


def main() -> None:
    """Screen explicit manifest entries and preserve every input hash."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    entries = list(csv.DictReader(args.manifest.open()))
    reference_path = ROOT / "cif/n12c.cif"
    reference = _cif_atoms(reference_path, "A", {r for _, r, _ in SITE_RESIDUES})
    rows, sources = [], []
    for entry in entries:
        path = Path(entry["path"])
        atoms = read_atoms(path)
        sources.append(
            {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        )
        names = entry["atom_names"].split(";")
        ligands = sorted({(ch, res) for ch, res, atom in atoms if atom == names[0]})
        if not ligands:
            raise ValueError(f"No mapped PA6 central amide: {path}")
        for site_index in range(4):
            site = layout_site(atoms, site_index)
            for chain, res in ligands:
                ligand = {
                    a: atoms[(chain, res, source)][1]
                    for a, source in zip(
                        ("C54", "O49", "N59", "C62"), names, strict=True
                    )
                }
                row = measure(site, ligand, reference)
                rows.append(
                    {
                        **entry,
                        "site": "ABCD"[site_index],
                        "ligand": f"{chain}:{res}",
                        **row,
                    }
                )
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "all_sites.csv").open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for path in (
        reference_path,
        Path(__file__),
        ROOT / "src/nyl12_covalent/active_site_geometry.py",
        ROOT / "src/nyl12_covalent/catalytic_alignment_geometry.py",
        LEGACY / "scripts/analyze_nyl12_nma_four_site_geometry.py",
        args.manifest,
    ):
        sources.append(
            {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        )
    summary = {
        "models": len(entries),
        "site_ligand_pairs": len(rows),
        "input_hashes": sources,
        "interpretation": "Opposite signed catalytic phi range identifies flipped local pose with fixed atom mapping. Native gates unchanged. Raw and relaxed ensembles remain separate.",
        "joint_passes": [r for r in rows if r["joint_pass"]],
        "tier_ab": [r for r in rows if r["tier"] in ("A", "B")],
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    logging.info(
        "Screened %d models; %d joint site passes",
        len(entries),
        len(summary["joint_passes"]),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
