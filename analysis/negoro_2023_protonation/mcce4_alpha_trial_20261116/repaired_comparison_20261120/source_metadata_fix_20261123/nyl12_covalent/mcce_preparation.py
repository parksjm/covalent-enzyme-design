"""Prepare an unchanged protein PDB with MCCE-compatible atom-name fields."""

from __future__ import annotations

import hashlib
import json
import logging
import math
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)
CANONICAL_RESIDUES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)
HEAVY_ELEMENTS = frozenset({"C", "N", "O", "S"})
PDB_METADATA_RECORDS = frozenset(
    {b"HEADER", b"TITLE ", b"REMARK", b"SEQRES", b"CRYST1"}
)


def prepare_input(
    source: Path,
    output: Path,
    expected_sha256: str,
    expected_termini: dict[str, dict[str, int]],
) -> dict[str, Any]:
    """Create a new MCCE staging directory without changing protein atoms.

    The source digest and exact per-chain ``n``/``c`` residue numbers are
    required. Only atom-name columns 13–16 are normalized. The original PDB,
    including OXT, is copied byte-for-byte to the output directory.

    Args:
        source: Registered, HIS-normalized, protein-only PDB input.
        output: New directory for original, prepared and provenance files.
        expected_sha256: Reviewed SHA-256 digest of ``source``.
        expected_termini: Chain-to-``{"n": ..., "c": ...}`` residue numbers.

    Returns:
        A JSON-compatible normalization report.

    Raises:
        FileExistsError: If the output directory already exists.
        ValueError: If provenance, PDB structure, atom identities or termini
            differ from the explicit requirements.
    """
    if output.exists():
        raise FileExistsError(f"Output directory already exists: {output}")
    if not expected_sha256 or len(expected_sha256) != 64:
        raise ValueError("expected_sha256 must be a 64-character SHA-256 digest")
    if not expected_termini:
        raise ValueError("expected_termini must declare every protein chain")
    for chain, termini in expected_termini.items():
        if not chain or set(termini) != {"n", "c"}:
            raise ValueError("Each expected chain must declare exactly integer n and c")
        if any(type(value) is not int for value in termini.values()):
            raise ValueError("Expected n and c termini must be integer residue numbers")

    source_bytes = source.read_bytes()
    source_sha256 = hashlib.sha256(source_bytes).hexdigest()
    if source_sha256 != expected_sha256.lower():
        raise ValueError(
            f"Source SHA-256 mismatch: expected {expected_sha256}, got {source_sha256}"
        )

    prepared_bytes, changes, chains, atom_count, metadata_records = _normalize_pdb(
        source_bytes, expected_termini
    )
    output.mkdir(parents=True)
    (output / "original.pdb").write_bytes(source_bytes)
    prepared_path = output / "protein.pdb"
    prepared_path.write_bytes(prepared_bytes)
    report: dict[str, Any] = {
        "schema_version": 1,
        "source": str(source.resolve()),
        "source_sha256": source_sha256,
        "original_pdb": "original.pdb",
        "prepared_pdb": "protein.pdb",
        "prepared_sha256": hashlib.sha256(prepared_bytes).hexdigest(),
        "atom_count": atom_count,
        "chains": chains,
        "expected_termini": expected_termini,
        "changed_atom_names": changes,
        "coordinates_unchanged": True,
        "residue_names_unchanged": True,
        "atom_records_preserved": True,
        "non_atom_records_preserved": True,
        "oxt_records_preserved": True,
        "metadata_records_preserved": metadata_records,
    }
    (output / "normalization.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    LOGGER.info(
        "Prepared %d protein atoms; normalized %d atom-name fields in %s",
        atom_count,
        len(changes),
        output,
    )
    return report


def _normalize_pdb(
    source_bytes: bytes, expected_termini: dict[str, dict[str, int]]
) -> tuple[
    bytes,
    list[dict[str, Any]],
    dict[str, dict[str, int]],
    int,
    dict[str, dict[str, Any]],
]:
    lines = source_bytes.splitlines(keepends=True)
    output_lines: list[bytes] = []
    changes: list[dict[str, Any]] = []
    observed: dict[str, dict[str, int]] = {}
    identities: set[tuple[str, int, str]] = set()
    residue_atoms: set[tuple[str, int, str]] = set()
    ended_chains: set[str] = set()
    active_chain: str | None = None
    atom_count = 0
    metadata: dict[bytes, list[bytes]] = {}
    end_seen = False

    for line_number, raw_line in enumerate(lines, start=1):
        body = raw_line.rstrip(b"\r\n")
        ending = raw_line[len(body) :]
        record = body[:6]
        if end_seen and body.strip():
            raise ValueError(f"Nonblank record after END on line {line_number}")
        if record == b"ATOM  ":
            if len(body) < 78:
                raise ValueError(f"Truncated ATOM record on line {line_number}")
            if body[16:17] != b" ":
                raise ValueError(f"Alternate location on line {line_number}")
            if body[26:27] != b" ":
                raise ValueError(f"Insertion code on line {line_number}")
            atom_name = body[12:16].decode("ascii").strip()
            residue_name = body[17:20].decode("ascii").strip()
            chain = body[21:22].decode("ascii").strip()
            element = body[76:78].decode("ascii").strip().upper()
            if not atom_name or len(atom_name) > 4:
                raise ValueError(f"Invalid atom name on line {line_number}")
            if residue_name not in CANONICAL_RESIDUES:
                raise ValueError(
                    f"Noncanonical residue {residue_name!r} on line {line_number}"
                )
            if not chain:
                raise ValueError(f"Blank chain identifier on line {line_number}")
            if element not in HEAVY_ELEMENTS:
                raise ValueError(
                    f"Non-heavy or unsupported element {element!r} on line {line_number}"
                )
            try:
                residue_number = int(body[22:26])
                coordinates = [float(body[start : start + 8]) for start in (30, 38, 46)]
            except ValueError as error:
                raise ValueError(
                    f"Invalid residue number or coordinates on line {line_number}"
                ) from error
            if not all(math.isfinite(value) for value in coordinates):
                raise ValueError(f"Non-finite coordinates on line {line_number}")

            if chain != active_chain:
                if active_chain is not None:
                    ended_chains.add(active_chain)
                if chain in ended_chains:
                    raise ValueError(f"Chain {chain!r} occurs in multiple segments")
                active_chain = chain
            chain_terms = observed.setdefault(
                chain, {"n": residue_number, "c": residue_number}
            )
            chain_terms["c"] = residue_number
            residue_key = (chain, residue_number, body[17:20].decode("ascii"))
            identity = (chain, residue_number, atom_name)
            if identity in identities:
                raise ValueError(
                    f"Duplicate atom identity {chain}:{residue_number}:{atom_name}"
                )
            identities.add(identity)
            residue_atoms.add(residue_key)

            normalized = _canonical_atom_field(atom_name)
            if body[12:16] != normalized:
                changes.append(
                    {
                        "line": line_number,
                        "chain": chain,
                        "residue_number": residue_number,
                        "residue_name": residue_name,
                        "atom_name": atom_name,
                        "before_field": body[12:16].decode("ascii"),
                        "after_field": normalized.decode("ascii"),
                    }
                )
                body = body[:12] + normalized + body[16:]
            output_lines.append(body + ending)
            atom_count += 1
        elif record.strip() == b"TER":
            if active_chain is None:
                raise ValueError(f"TER without preceding chain on line {line_number}")
            ended_chains.add(active_chain)
            active_chain = None
            output_lines.append(raw_line)
        elif record.strip() == b"END":
            end_seen = True
            output_lines.append(raw_line)
        elif record in PDB_METADATA_RECORDS:
            metadata.setdefault(record, []).append(raw_line)
            output_lines.append(raw_line)
        elif not body.strip():
            output_lines.append(raw_line)
        else:
            raise ValueError(
                f"Unsupported PDB record {record.decode('ascii', 'replace')!r} "
                f"on line {line_number}"
            )

    if not end_seen:
        raise ValueError("PDB is missing its END record")
    if observed != expected_termini:
        raise ValueError(
            f"Observed termini {observed} do not match expected {expected_termini}"
        )
    if len(residue_atoms) == 0:
        raise ValueError("PDB contains no protein ATOM records")
    metadata_report = {
        record.decode("ascii").strip(): {
            "count": len(records),
            "sha256": hashlib.sha256(b"".join(records)).hexdigest(),
        }
        for record, records in sorted(metadata.items())
    }
    return b"".join(output_lines), changes, observed, atom_count, metadata_report


def _canonical_atom_field(atom_name: str) -> bytes:
    """Return the MCCE four-column atom-name representation."""
    if len(atom_name) == 4:
        return atom_name.encode("ascii")
    return f" {atom_name:<3}".encode("ascii")
