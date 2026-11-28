"""Compare unrelaxed predictions using protein-only physical subunit pairing."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.af3_postprocess import read_atom_sites
from nyl12_covalent.catalytic_alignment_geometry import fit_coordinates, state_metrics
from nyl12_covalent.relay_any_amide import TIERS
from nyl12_covalent.relay_broad_prediction_io import audit_prediction, export_geometry

ROOT = Path("/data/jpc/covalent-enzyme-design")
PILOT = ROOT / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot"
REFERENCE = ROOT / "analysis/native_thr_gate_calibration_20261123/inputs/n12c.cif"
LOG = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return the SHA256 of a retained input or output."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data: Any) -> None:
    """Write finite, readable authoritative JSON."""
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def atom_loop(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read the simple atom-site loop with exact source column identities."""
    columns, rows = [], []
    started = False
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            columns.append(line.strip().split(".")[1])
            started = True
        elif started and line.strip() == "#":
            break
        elif started and line.strip():
            row = line.split()
            if len(row) != len(columns):
                raise ValueError(f"Unsupported atom-site syntax: {path}")
            rows.append(row)
    if not rows:
        raise ValueError(f"No atom-site records: {path}")
    return columns, rows


def reference_atoms() -> dict[str, dict[tuple[int, str], tuple[str, np.ndarray]]]:
    """Read canonical precursor-numbered crystal proteins, retaining coordinates."""
    columns, rows = atom_loop(REFERENCE)
    result = {}
    for values in rows:
        d = dict(zip(columns, values, strict=True))
        if d["group_PDB"] != "ATOM" or d["label_alt_id"] not in (".", "?", "A"):
            continue
        chain = result.setdefault(d["auth_asym_id"], {})
        key = (int(d["auth_seq_id"]), d["label_atom_id"])
        chain[key] = (
            d["type_symbol"],
            np.array([float(d[f"Cartn_{a}"]) for a in "xyz"]),
        )
    return result


def pairing(atoms: dict, refs: dict) -> tuple[dict[str, str], dict]:
    """Assign four beta chains by unique maximum core CA interface contacts.

    Interface counts use alpha20-215 and beta10-85. Crystal fit diagnostics
    exclude Phe113, Asn186, Asp272/274.
    Thr233 is held out. No ligand atoms or catalytic geometry select partners.
    """
    scores, details = {}, {}
    for alpha, beta in itertools.product("ABCD", "EFGH"):
        fits = []
        for refname, ref in refs.items():
            keys = [(r, "CA") for r in range(20, 216) if r not in (113, 186)]
            keys += [(r + 232, "CA") for r in range(10, 86) if r not in (40, 42)]
            keys = [k for k in keys if k in ref]
            mobile = np.array(
                [
                    atoms[(alpha if r <= 232 else beta, r if r <= 232 else r - 232, a)][
                        1
                    ]
                    for r, a in keys
                ]
            )
            target = np.array([ref[k][1] for k in keys])
            rmsd = fit_coordinates(mobile, target)[3]
            fits.append((rmsd, refname, len(keys)))
        value, refname, count = min(fits)
        scores[alpha + beta] = value
        details[alpha + beta] = {
            "best_reference_chain": refname,
            "ca_count": count,
            "rmsd_A": value,
            "reference_fits": fits,
        }
    choices = sorted(
        (sum(scores[a + b] for a, b in zip("ABCD", p, strict=True)), "".join(p))
        for p in itertools.permutations("EFGH")
    )
    ca = {
        c: np.array(
            [
                x
                for (chain, r, a), (_, x) in atoms.items()
                if chain == c
                and a == "CA"
                and (20 <= r <= 215 if c in "ABCD" else 10 <= r <= 85)
            ]
        )
        for c in "ABCDEFGH"
    }
    contacts = {
        a + b: int(
            np.sum(np.linalg.norm(ca[a][:, None, :] - ca[b][None, :, :], axis=2) < 8.0)
        )
        for a, b in itertools.product("ABCD", "EFGH")
    }
    contact_choices = sorted(
        (
            (sum(contacts[a + b] for a, b in zip("ABCD", p, strict=True)), "".join(p))
            for p in itertools.permutations("EFGH")
        ),
        reverse=True,
    )
    best, perm = contact_choices[0]
    if best == contact_choices[1][0]:
        raise ValueError("Protein core-contact assignment is ambiguous")
    mapping = dict(zip("ABCD", perm, strict=True))
    report = {
        "physical_alpha_to_beta": mapping,
        "pair_fit_matrix": details,
        "all_24_assignments": choices,
        "best_sum_rmsd_A": choices[0][0],
        "core_CA_contacts_under8A": contacts,
        "all24_contact_assignments": contact_choices,
        "contact_assignment_gap": best - contact_choices[1][0],
        "runner_up_sum_rmsd_A": choices[1][0],
        "assignment_gap_A": choices[1][0] - choices[0][0],
        "selection": "maximum total protein core CA contacts below8A, one-to-one; unique optimum; no ligand coordinates. Crystal pair-fit scores retained as poor-fit diagnostic, not pairing selection.",
        "core": "alpha20-215 except113/186; beta10-85 except40/42; CA atoms present in crystal",
        "crystal_path": str(REFERENCE),
        "crystal_sha256": digest(REFERENCE),
    }
    return mapping, report


def canonical_cif(source: Path, output: Path, mapping: dict[str, str]) -> dict:
    """Write explicitly derived atom-site-only CIF with beta labels canonicalized."""
    rename = {
        actual: canonical
        for actual, canonical in zip(mapping.values(), "EFGH", strict=True)
    }
    columns, rows = atom_loop(source)
    for row in rows:
        for col in ("auth_asym_id", "label_asym_id"):
            idx = columns.index(col)
            row[idx] = rename.get(row[idx], row[idx])
    text = "data_DERIVED_PHYSICAL_PAIR_CANONICALIZATION\n#\nloop_\n"
    text += "\n".join("_atom_site." + col for col in columns) + "\n"
    text += "\n".join(" ".join(row) for row in rows) + "\n#\n"
    output.write_text(text)
    return {
        "actual_to_derived_beta_chain": rename,
        "source": str(source),
        "source_sha256": digest(source),
        "derived": str(output),
        "derived_sha256": digest(output),
        "changes": "Atom-site chain labels only; coordinates/atom ordering unchanged. Atom-site-only derived view; original full CIF immutable.",
    }


def failed_gates(row: dict, tier: tuple) -> list[str]:
    """List exact failed geometric inequalities for a requested tier."""
    _, d1, oxy, amin, amax, pmin, pmax = tier
    gates = {
        "amide_geometry_valid": row["amide_geometry_valid"],
        "d1": row["d1"] <= d1,
        "d2": row["d2"] <= oxy,
        "d3": row["d3"] <= oxy,
        "a_bd": amin <= row["a_bd"] <= amax,
        "abs_phi": pmin <= abs(row["herschlag_phi"]) <= pmax,
    }
    return [k for k, accepted in gates.items() if not accepted]


def analyze(source: Path, output: Path, predictor: str) -> list[dict]:
    """Audit one native model and export all48 plus four assigned geometries."""
    output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads(
        (PILOT / "prediction_preparation/ligand_graph_manifest.json").read_text()
    )
    maps = manifest["ligand_atom_maps"][predictor]
    sequences = json.loads((PILOT / "canonical_protection.json").read_text())[
        "source_sequences"
    ]
    write_json(
        output / "raw_software_acceptance.json",
        audit_prediction(source, sequences["alpha"], sequences["beta"], maps["I"]),
    )
    atoms = {
        (s.chain, s.sequence, s.name): (s.element, s.xyz)
        for s in read_atom_sites(source)
    }
    refs = reference_atoms()
    mapping, evidence = pairing(atoms, refs)
    write_json(output / "pairing_evidence.json", evidence)
    derived = output / "derived_pair_canonicalized.cif"
    receipt = canonical_cif(source, derived, mapping)
    write_json(output / "canonicalization_receipt.json", receipt)
    write_json(
        output / "derived_software_acceptance.json",
        audit_prediction(derived, sequences["alpha"], sequences["beta"], maps["I"]),
    )
    export_geometry(derived, maps, output)
    report = json.loads((output / "geometry_any_amide.json").read_text())
    rows = []
    for row in report["assigned"]:
        alpha = row["site"][0]
        beta = mapping[alpha]
        ref = refs[evidence["pair_fit_matrix"][alpha + beta]["best_reference_chain"]]
        mobile = {
            (r if c == alpha else r + 232, a): (e, np.array(x))
            for (c, r, a), (e, x) in atoms.items()
            if c in (alpha, beta)
        }
        keys = [(r, a) for r in (113, 186, 272, 274) for a in ("N", "CA", "C", "O")]
        rotation, mc, rc, rmsd = fit_coordinates(
            np.array([mobile[k][1] for k in keys]), np.array([ref[k][1] for k in keys])
        )
        aligned = {k: (e, (x - mc) @ rotation + rc) for k, (e, x) in mobile.items()}
        native = state_metrics(ref, aligned)
        write_json(
            output / f"native_thr_diagnostic_{alpha}{beta}.json",
            {"fit_rmsd_A": rmsd, **native},
        )
        rows.append(
            {
                "model": source.stem,
                "actual_physical_pair": alpha + beta,
                **row,
                "native_thr_og1_displacement_A": native["thr233_held_out"][
                    "heavy_atom_displacements_angstrom"
                ]["OG1"],
                **{f"failed_{t[0]}": ",".join(failed_gates(row, t)) for t in TIERS},
            }
        )
    write_json(output / "assigned_with_diagnostics.json", rows)
    LOG.info(
        "%s pairing=%s gap=%.3f A", source.name, mapping, evidence["assignment_gap_A"]
    )
    return rows


def main() -> None:
    """Run reproducible analysis for explicit native prediction CIF inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif", type=Path, action="append", required=True)
    parser.add_argument("--predictor", choices=("boltz", "af3"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for index, cif in enumerate(args.cif):
        rows.extend(
            analyze(
                cif.resolve(),
                args.output / f"{args.predictor}_{index + 1}_{cif.stem}_contact_pairs",
                args.predictor,
            )
        )
    write_json(args.output / f"{args.predictor}_assigned_all.json", rows)
    with (args.output / f"{args.predictor}_assigned_all.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    LOG.info(
        "Completed %s: %d assigned sites nested within %d models",
        args.predictor,
        len(rows),
        len(args.cif),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
