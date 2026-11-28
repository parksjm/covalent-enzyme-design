"""Split one homolog-origin precursor alignment into paired mature-chain MSAs."""

from __future__ import annotations

from collections.abc import Sequence
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYX-")


def read_precursor_a3m(path: Path) -> list[tuple[str, str]]:
    """Read multiline A3M records without inventing homolog relationships.

    Args:
        path: One full-precursor alignment, with its query as the first record.

    Returns:
        Ordered source header and alignment-string pairs. Optional single-entity
        ColabFold length metadata is validated; multi-entity metadata is rejected.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    metadata_length: int | None = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            fields = line[1:].split()
            if (
                header is not None
                or metadata_length is not None
                or len(fields) != 2
                or not fields[0].isdigit()
                or fields[1] != "1"
                or int(fields[0]) < 1
            ):
                raise ValueError(
                    "Only single-precursor ColabFold '#length 1' metadata is supported"
                )
            metadata_length = int(fields[0])
        elif line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header, parts = line[1:], []
            if not header:
                raise ValueError("A3M header is empty")
        elif header is None:
            raise ValueError("A3M sequence precedes its header")
        else:
            parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    if not records:
        raise ValueError("Full precursor A3M is empty")
    if metadata_length is not None and len(records[0][1]) != metadata_length:
        raise ValueError("ColabFold metadata length does not match its precursor query")
    return records


def split_precursor_records(
    records: Sequence[tuple[str, str]], *, alpha_query: str, beta_query: str
) -> dict[str, Any]:
    """Split query-aligned columns while preserving each single homolog origin.

    Lowercase insertion residues at the cleavage boundary join the beta row,
    consistent with native Boltz deletion counts attaching to the next column.
    A row duplicated in either mature part is excluded from both parts, matching
    native CSV deduplication without creating an orphaned or mispaired row.

    Args:
        records: Ordered records from one full-precursor A3M.
        alpha_query: Exact uppercase mature alpha sequence.
        beta_query: Exact uppercase mature beta sequence.

    Returns:
        Native alpha/beta CSV rows, stable source row keys and explicit exclusions.
    """
    query = alpha_query + beta_query
    if (
        not records
        or not alpha_query
        or not beta_query
        or any(c not in ALPHABET - {"-"} for c in query)
        or records[0][1] != query
    ):
        raise ValueError("First precursor row must exactly equal both mature queries")
    boundary = len(alpha_query)
    alpha_rows = [{"sequence": alpha_query, "key": -1}]
    beta_rows = [{"sequence": beta_query, "key": -1}]
    origins, excluded = [], []
    seen_alpha, seen_beta = {alpha_query}, {beta_query}
    for ordinal, (header, sequence) in enumerate(records):
        if not header or not sequence:
            raise ValueError("Every precursor record requires a header and sequence")
        alpha, beta, column = [], [], 0
        for letter in sequence:
            if letter.upper() not in ALPHABET:
                raise ValueError(f"Unsupported A3M character in source row {ordinal}")
            (alpha if column < boundary else beta).append(letter)
            if not letter.islower():
                column += 1
        if column != len(query):
            raise ValueError(f"Source row {ordinal} does not span the query columns")
        if ordinal == 0:
            continue
        a, b = "".join(alpha), "".join(beta)
        identity_a, identity_b = a.replace("-", "").upper(), b.replace("-", "").upper()
        reason = None
        if not identity_a or not identity_b:
            reason = "no_residues_in_one_mature_part"
        elif identity_a in seen_alpha or identity_b in seen_beta:
            reason = "native_duplicate_in_one_mature_part"
        if reason:
            excluded.append({"source_row": ordinal, "header": header, "reason": reason})
            continue
        seen_alpha.add(identity_a)
        seen_beta.add(identity_b)
        alpha_rows.append({"sequence": a, "key": ordinal})
        beta_rows.append({"sequence": b, "key": ordinal})
        origins.append({"key": ordinal, "source_row": ordinal, "header": header})
    if not origins:
        raise ValueError("No independently retained precursor homolog pairs")
    return {
        "alpha_rows": alpha_rows,
        "beta_rows": beta_rows,
        "row_origins": origins,
        "excluded_rows": excluded,
        "boundary_after_query_residue1": boundary,
        "query_length": len(query),
        "boundary_insertions": "beta; insertion counts precede its first aligned column",
        "source_row_count": len(records),
        "paired_homolog_row_count": len(origins),
        "query_pairing_key": -1,
        "homolog_keys": "original single-precursor source ordinal; never independent labels",
    }


def write_paired_precursor_msas(
    source: Path, output: Path, *, alpha_query: str, beta_query: str
) -> dict[str, Any]:
    """Write native CSVs and a reversible source-origin receipt to a new directory.

    Args:
        source: Validated single full-precursor A3M source.
        output: New isolated MSA output directory.
        alpha_query: Exact mature alpha query.
        beta_query: Exact mature beta query.

    Returns:
        Source and output hashes plus all retained/excluded homolog origins.
    """
    result = split_precursor_records(
        read_precursor_a3m(source), alpha_query=alpha_query, beta_query=beta_query
    )
    output.mkdir(parents=True, exist_ok=False)
    for name in ("alpha", "beta"):
        path = output / f"{name}_paired.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sequence", "key"])
            writer.writeheader()
            writer.writerows(result.pop(f"{name}_rows"))
    result["source"] = str(source.resolve())
    result["source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    result["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob("*.csv")
    }
    (output / "precursor_pairing_receipt.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result
