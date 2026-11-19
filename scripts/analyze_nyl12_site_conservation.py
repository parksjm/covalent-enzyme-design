"""Audit existing split-chain Boltz MSAs; no remote BLAST submission occurs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)
STANDARD = set("ACDEFGHIKLMNPQRSTVWYX-")
SITES = {"alpha": {113: 113, 156: 156, 186: 186}, "beta": {233: 1, 272: 40, 274: 42}}


def aligned_sequence(sequence: str) -> str:
    """Remove lowercase insertion characters without shifting match columns.

    Args:
        sequence: A3M-style aligned amino-acid sequence.

    Returns:
        Uppercase match-column sequence, retaining deletion gaps.

    Raises:
        ValueError: If the remaining sequence is empty or contains unknown codes.
    """
    result = re.sub("[a-z]", "", sequence)
    if not result or set(result) - STANDARD:
        raise ValueError("Invalid aligned amino-acid sequence")
    return result


def site_counts(sequences: list[str], position: int, residue: str) -> dict[str, Any]:
    """Count a site with explicit gap and ambiguous-residue denominators.

    Args:
        sequences: Validated equal-length aligned sequences.
        position: One-based match-column index.
        residue: Query amino acid at this column.

    Returns:
        Counts and query frequency among non-gap and standard-residue calls.

    Raises:
        ValueError: If no rows exist or the position is invalid.
    """
    if not sequences or not 1 <= position <= min(map(len, sequences)):
        raise ValueError("Empty alignment or invalid site")
    counts = Counter(sequence[position - 1] for sequence in sequences)
    total = len(sequences)
    nongap = total - counts["-"]
    standard = nongap - counts["X"]
    return {
        "counts": dict(sorted(counts.items())),
        "rows": total,
        "nongap": nongap,
        "gaps": counts["-"],
        "ambiguous_X": counts["X"],
        "gap_fraction": counts["-"] / total,
        "query_fraction_nongap": counts[residue] / nongap if nongap else None,
        "query_fraction_standard": counts[residue] / standard if standard else None,
    }


def analyze_msa(path: Path, chain: str, coverage: float) -> dict[str, Any]:
    """Summarize an existing MSA after exact deduplication and coverage filtering.

    Args:
        path: Existing Boltz CSV with key and sequence columns; query first.
        chain: Explicit alpha or beta numbering map.
        coverage: Required fraction of query columns occupied by non-gap calls.

    Returns:
        Source hash, query, filter accounting and per-site raw/filtered counts.

    Raises:
        ValueError: If query, rows, columns or filtering prerequisites fail.
    """
    if chain not in SITES or not 0 < coverage <= 1:
        raise ValueError("Invalid chain or coverage threshold")
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows or rows[0].get("key") != "0":
        raise ValueError("First Boltz CSV row must be query key0")
    sequences = [aligned_sequence(row["sequence"]) for row in rows]
    query = sequences[0]
    if "-" in query or "X" in query or {len(s) for s in sequences} != {len(query)}:
        raise ValueError("Query must be complete; match-column lengths must agree")
    unique = sorted(set(sequences))
    nonquery = [s for s in unique if s != query]
    accepted = [
        s for s in nonquery if sum(c != "-" for c in s) / len(query) >= coverage
    ]
    if not accepted:
        raise ValueError("No nonquery sequences pass coverage filter")
    sites = {}
    for full_position, column in SITES[chain].items():
        residue = query[column - 1]
        sites[str(full_position)] = {
            "match_column_1based": column,
            "query_residue": residue,
            "query_motif": query[max(0, column - 6) : column + 5],
            "raw": site_counts(sequences, column, residue),
            "unique_including_query": site_counts(unique, column, residue),
            "filtered_unique_nonquery": site_counts(accepted, column, residue),
        }
    result = {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "query": query,
        "query_length": len(query),
        "raw_rows": len(rows),
        "unique_rows_including_query": len(unique),
        "unique_nonquery": len(nonquery),
        "accepted_unique_nonquery": len(accepted),
        "coverage_threshold": coverage,
        "rejected_coverage": len(nonquery) - len(accepted),
        "sites": sites,
    }
    if chain == "beta":
        result["D272_D274_pair"] = {
            "DD": sum(s[39] == s[41] == "D" for s in accepted),
            "rows": len(accepted),
        }
    return result


def global_alignment(query: str, reference: str) -> tuple[str, str, dict[int, int]]:
    """Align explicit precursor sequences using reproducible linear-gap scores.

    Args:
        query: Ungapped Nyl12 concatenated modeled alpha and beta sequence.
        reference: Public NylC precursor sequence with native numbering.

    Returns:
        Aligned query/reference strings and one-based query-to-reference map.
        Scores are match2, mismatch-1, gap-2; ties prefer diagonal, up, left.
    """
    score = [[0] * (len(reference) + 1) for _ in range(len(query) + 1)]
    for i in range(len(query) + 1):
        score[i][0] = -2 * i
    for j in range(len(reference) + 1):
        score[0][j] = -2 * j
    for i in range(1, len(query) + 1):
        for j in range(1, len(reference) + 1):
            diagonal = score[i - 1][j - 1] + (
                2 if query[i - 1] == reference[j - 1] else -1
            )
            score[i][j] = max(diagonal, score[i - 1][j] - 2, score[i][j - 1] - 2)
    i, j = len(query), len(reference)
    left, right, mapping = [], [], {}
    while i or j:
        if (
            i
            and j
            and score[i][j]
            == score[i - 1][j - 1] + (2 if query[i - 1] == reference[j - 1] else -1)
        ):
            left.append(query[i - 1])
            right.append(reference[j - 1])
            mapping[i] = j
            i, j = i - 1, j - 1
        elif i and score[i][j] == score[i - 1][j] - 2:
            left.append(query[i - 1])
            right.append("-")
            i -= 1
        else:
            left.append("-")
            right.append(reference[j - 1])
            j -= 1
    return "".join(reversed(left)), "".join(reversed(right)), mapping


def main() -> None:
    """Write source-hashed analysis artifacts from explicitly supplied inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha-msa", type=Path, required=True)
    parser.add_argument("--beta-msa", type=Path, required=True)
    parser.add_argument("--nylc-fasta", type=Path, required=True)
    parser.add_argument("--nylc-origin-url", required=True)
    parser.add_argument("--coverage", type=float, default=0.8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    alpha = analyze_msa(args.alpha_msa, "alpha", args.coverage)
    beta = analyze_msa(args.beta_msa, "beta", args.coverage)
    fasta = args.nylc_fasta.read_text()
    reference = "".join(
        line.strip() for line in fasta.splitlines() if not line.startswith(">")
    )
    if (
        not fasta.startswith(">")
        or not reference
        or set(reference) - (STANDARD - {"-", "X"})
    ):
        raise ValueError("Expected one complete standard-residue public NylC FASTA")
    query = alpha["query"] + beta["query"]
    aligned_query, aligned_reference, mapping = global_alignment(query, reference)
    correspondence = {}
    for position in (113, 156, 186, 233, 272, 274):
        if position not in mapping:
            raise ValueError(f"Requested site {position} is not mapped")
        reference_position = mapping[position]
        correspondence[str(position)] = {
            "Nyl12": query[position - 1],
            "NylC_position": reference_position,
            "NylC": reference[reference_position - 1],
        }
    report = {
        "analysis": "existing Boltz MSA; no fresh BLAST or external sequence submission",
        "methods": {
            "lowercase_insertions": "removed",
            "deduplication": "exact aligned sequence",
            "query": "all exact query matches excluded",
            "coverage": "non-gap/query match columns",
            "phylogenetic_weighting": False,
            "alignment_scores": "match2 mismatch-1 gap-2",
            "limitations": "split-chain broad-family candidates; unweighted; identities/taxa unavailable; conservation does not establish catalytic roles",
        },
        "alpha": alpha,
        "beta": beta,
        "correspondence": correspondence,
        "reference": {
            "origin_url": args.nylc_origin_url,
            "fasta_header": fasta.splitlines()[0],
            "sha256": hashlib.sha256(args.nylc_fasta.read_bytes()).hexdigest(),
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "command_argv": sys.argv,
        "python_version": sys.version,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "conservation.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.output / "public_nylc_3axg.fasta").write_text(fasta)
    (args.output / "nyl12_nylc_numbering_alignment.fasta").write_text(
        f">Nyl12_modeled_alpha_plus_beta\n{aligned_query}\n>NylC_3AXG\n{aligned_reference}\n"
    )
    LOGGER.info("Saved conservation audit to %s", args.output)


if __name__ == "__main__":
    main()
