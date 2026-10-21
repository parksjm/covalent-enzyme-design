"""Audit ligand-specific confidence from the original Nyl12 AF3 outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import math
import shlex
import statistics
from pathlib import Path
from typing import Any

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from nyl12_covalent.chemistry import (  # noqa: E402
    build_poly_xama_maxa,
    product_report,
)

LOGGER = logging.getLogger(__name__)
CHAIN_ID = "I"
CONFIDENCE_BINS = (("lt50", 0.0, 50.0), ("50to70", 50.0, 70.0), ("ge70", 70.0, 101.0))


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file.

    Args:
        path: File to hash.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _single_match(directory: Path, suffix: str) -> Path:
    matches = sorted(directory.glob(f"*{suffix}"))
    if suffix == "_confidences.json":
        matches = [path for path in matches if "_summary_confidences.json" not in path.name]
    if len(matches) != 1:
        raise ValueError(f"Expected one *{suffix} in {directory}; found {len(matches)}")
    return matches[0]


def _cif_atom_rows(path: Path) -> list[dict[str, str]]:
    """Read atom_site rows from a CIF without adding a parser dependency."""
    lines = path.read_text(encoding="utf-8").splitlines()
    rows: list[dict[str, str]] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        headers: list[str] = []
        while index < len(lines) and lines[index].startswith("_atom_site."):
            headers.append(lines[index].strip())
            index += 1
        if not headers:
            continue
        while index < len(lines):
            line = lines[index].strip()
            if not line or line.startswith("#"):
                index += 1
                if line.startswith("#"):
                    break
                continue
            if line == "loop_" or line.startswith("_"):
                break
            tokens = shlex.split(line)
            if len(tokens) != len(headers):
                raise ValueError(f"Malformed atom_site row in {path}: {line}")
            rows.append(dict(zip(headers, tokens, strict=True)))
            index += 1
    if not rows:
        raise ValueError(f"No _atom_site rows found in {path}")
    return rows


def _distribution(values: list[float]) -> dict[str, Any]:
    """Summarize confidence values and common pLDDT bins."""
    if not values:
        raise ValueError("Cannot summarize an empty confidence set")
    bins = {
        label: {"count": sum(lower <= value < upper for value in values),
                "fraction": sum(lower <= value < upper for value in values) / len(values)}
        for label, lower, upper in CONFIDENCE_BINS
    }
    return {
        "count": len(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "minimum": min(values),
        "maximum": max(values),
        "population_sd": statistics.pstdev(values),
        "bins": bins,
    }


def _matrix_summary(matrix: list[list[float]], row_indices: list[int],
                    column_indices: list[int]) -> dict[str, float | int]:
    """Summarize a rectangular selection from a token-level matrix."""
    values = [matrix[row][column] for row in row_indices for column in column_indices]
    if not values:
        raise ValueError("PAE selection contains no ligand/protein token pairs")
    return {"count": len(values), "mean": statistics.fmean(values),
            "median": statistics.median(values), "minimum": min(values),
            "maximum": max(values)}


def _ligand_atoms(cif_rows: list[dict[str, str]], confidence: dict[str, Any],
                  chain_id: str) -> dict[str, float]:
    """Join chain atom names to AF3 per-atom pLDDT and validate CIF B values."""
    chain_ids = confidence["atom_chain_ids"]
    plddts = confidence["atom_plddts"]
    if len(chain_ids) != len(plddts) or len(chain_ids) != len(cif_rows):
        raise ValueError("AF3 atom confidence arrays do not align with CIF atom rows")
    result: dict[str, float] = {}
    for row, row_chain, raw_value in zip(cif_rows, chain_ids, plddts, strict=True):
        label_chain = row["_atom_site.label_asym_id"]
        if label_chain != row_chain:
            raise ValueError("AF3 atom chain IDs do not align with CIF atom rows")
        if row_chain != chain_id:
            continue
        name = row["_atom_site.label_atom_id"]
        if name in result:
            raise ValueError(f"Duplicate atom name in ligand chain {chain_id}: {name}")
        value = float(raw_value)
        cif_value = float(row["_atom_site.B_iso_or_equiv"])
        if (not math.isfinite(value) or not math.isfinite(cif_value)
                or not 0.0 <= value <= 100.0 or not 0.0 <= cif_value <= 100.0):
            raise ValueError(f"Invalid pLDDT value for ligand atom {name}")
        if abs(value - cif_value) > 0.02:
            raise ValueError(f"CIF B_iso and AF3 atom_plddts disagree for {name}")
        result[name] = value
    if not result:
        raise ValueError(f"No ligand atoms found on chain {chain_id}")
    return result


def _monad_report(atom_confidence: dict[str, float]) -> dict[str, Any]:
    """Map AF3 ligand atom confidence onto the chemistry builder's monads."""
    substrate = build_poly_xama_maxa()
    groups: dict[int, dict[str, Any]] = {}
    for atom in substrate.mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        position = atom.GetIntProp("monad_index")
        name = atom.GetProp("name")
        entry = groups.setdefault(position, {"kind": atom.GetProp("monad_kind"), "names": []})
        entry["names"].append(name)
    if set(atom_confidence) != {name for group in groups.values() for name in group["names"]}:
        raise ValueError("AF3 ligand heavy-atom names do not match the chemistry builder")
    monads: list[dict[str, Any]] = []
    for position, entry in sorted(groups.items()):
        values = [atom_confidence[name] for name in entry["names"]]
        monads.append({"position": position, "kind": entry["kind"],
                       "atom_names": entry["names"], "confidence": _distribution(values)})
    tail_names = [name for position in (1, max(groups)) for name in groups[position]["names"]]
    core_names = [name for position, group in groups.items() if position not in (1, max(groups))
                  for name in group["names"]]
    return {"monads": monads,
            "terminal_monads_1_and_8": _distribution([atom_confidence[name] for name in tail_names]),
            "internal_monads_2_through_7": _distribution([atom_confidence[name] for name in core_names])}


def audit_model(model_dir: Path, target: str) -> dict[str, Any]:
    """Audit one AF3 model's ligand confidence, PAE, scores, and provenance.

    Args:
        model_dir: Directory containing one model CIF and AF3 JSON/CSV outputs.
        target: Stable target label, either ``poly_xama_maxa`` or ``pa6_l2``.

    Returns:
        JSON-compatible audit record with ligand-specific confidence and hashes.

    Raises:
        ValueError: If required AF3 data are missing, inconsistent, or malformed.
    """
    if target not in {"poly_xama_maxa", "pa6_l2"}:
        raise ValueError(f"Unsupported target: {target}")
    model_dir = model_dir.resolve()
    cif = _single_match(model_dir, "_model.cif")
    confidence_file = _single_match(model_dir, "_confidences.json")
    summary_file = _single_match(model_dir, "_summary_confidences.json")
    ranking_file = _single_match(model_dir, "_ranking_scores.csv")
    confidence = json.loads(confidence_file.read_text(encoding="utf-8"))
    summary = json.loads(summary_file.read_text(encoding="utf-8"))
    rows = _cif_atom_rows(cif)
    ligand = _ligand_atoms(rows, confidence, CHAIN_ID)
    token_chains = confidence["token_chain_ids"]
    pae = confidence["pae"]
    if len(token_chains) != len(pae) or any(len(row) != len(token_chains) for row in pae):
        raise ValueError("AF3 PAE matrix does not match token_chain_ids")
    if any(not math.isfinite(float(value)) for row in pae for value in row):
        raise ValueError("AF3 PAE matrix contains non-finite values")
    ligand_tokens = [i for i, chain in enumerate(token_chains) if chain == CHAIN_ID]
    protein_tokens = [i for i, chain in enumerate(token_chains) if chain != CHAIN_ID]
    ids = sorted(set(token_chains))
    if len(summary["chain_pair_iptm"]) != len(ids):
        raise ValueError("Summary chain-pair arrays do not match AF3 chain IDs")
    chain_index = ids.index(CHAIN_ID)
    pair_iptm = summary["chain_pair_iptm"]
    pair_pae = summary["chain_pair_pae_min"]
    if (len(pair_iptm) != len(ids) or len(pair_pae) != len(ids)
            or any(len(row) != len(ids) for row in pair_iptm + pair_pae)):
        raise ValueError("Summary chain-pair arrays are not square")
    with ranking_file.open(newline="", encoding="utf-8") as handle:
        ranking_rows = list(csv.DictReader(handle))
    if len(ranking_rows) != 1:
        raise ValueError("Expected one AF3 ranking score row for the selected model")
    csv_ranking_score = float(ranking_rows[0]["ranking_score"])
    summary_ranking_score = float(summary["ranking_score"])
    if (csv_ranking_score != summary_ranking_score
            and round(csv_ranking_score, 2) != round(summary_ranking_score, 2)):
        raise ValueError("CSV and summary ranking scores disagree")

    record: dict[str, Any] = {
        "target": target,
        "source_model_directory": str(model_dir),
        "ligand_chain": CHAIN_ID,
        "ligand_atom_confidence": _distribution(list(ligand.values())),
        "ligand_atoms": ligand,
        "interchain_pae": {
            "ligand_to_protein": _matrix_summary(pae, ligand_tokens, protein_tokens),
            "protein_to_ligand": _matrix_summary(pae, protein_tokens, ligand_tokens),
            "by_protein_chain": {
                chain: {
                    "ligand_to_chain": _matrix_summary(
                        pae, ligand_tokens,
                        [i for i, token_chain in enumerate(token_chains) if token_chain == chain],
                    ),
                    "chain_to_ligand": _matrix_summary(
                        pae,
                        [i for i, token_chain in enumerate(token_chains) if token_chain == chain],
                        ligand_tokens,
                    ),
                }
                for chain in ids if chain != CHAIN_ID
            },
        },
        "summary_metrics": {
            "ranking_score": csv_ranking_score,
            "ranking_score_summary": summary_ranking_score,
            "ranking_score_precision_note": (
                "CSV and summary values are retained at source precision; acceptance "
                "requires exact equality or agreement after rounding to two decimals."
            ),
            "iptm": summary["iptm"], "ptm": summary["ptm"],
            "actifptm": summary["actifptm"],
            "chain_pair_iptm_ligand_vs_protein_mean": statistics.fmean(
                pair_iptm[i][chain_index] for i in range(len(ids)) if i != chain_index
            ),
            "chain_pair_pae_min_ligand_vs_protein_mean": statistics.fmean(
                pair_pae[i][chain_index] for i in range(len(ids)) if i != chain_index
            ),
        },
        "interpretation_limits": [
            "Ranking score and complex-level pTM/iPTM/actifPTM are not ligand atom pLDDT.",
            "Confidence scores do not establish chemical validity or product identity.",
        ],
    }
    if target == "poly_xama_maxa":
        substrate = build_poly_xama_maxa()
        reactive_names = (substrate.carbonyl_c, substrate.carbonyl_o, substrate.leaving_n)
        record["monad_confidence"] = _monad_report(ligand)
        record["chemical_identity"] = {
            "polymer_label": "poly(X,AMA)",
            "modeled_sequence": "MAXAMAXA",
            "scissile_display": "MAXA|MAXA",
            "product_report": product_report(substrate),
        }
    else:
        reactive_names = ("C54", "O49", "N59")
    record["reactive_atom_names"] = list(reactive_names)
    record["reactive_atoms"] = {name: ligand[name] for name in reactive_names}
    record["reactive_atom_confidence"] = _distribution([ligand[name] for name in reactive_names])

    source_root = Path(__file__).resolve().parents[1]
    input_paths = (cif, confidence_file, summary_file, ranking_file)
    record["sha256"] = {
        "inputs": {str(path): sha256_file(path) for path in input_paths},
        "audit_script": sha256_file(Path(__file__).resolve()),
        "chemistry_builder": sha256_file(source_root / "src/nyl12_covalent/chemistry.py"),
    }
    return record


def main() -> None:
    """Run the two-target AF3 ligand confidence audit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poly-model-dir", type=Path, required=True)
    parser.add_argument("--pa6-model-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {
        "schema_version": 1,
        "models": [audit_model(args.poly_model_dir, "poly_xama_maxa"),
                   audit_model(args.pa6_model_dir, "pa6_l2")],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    LOGGER.info("Wrote ligand confidence audit to %s", args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
