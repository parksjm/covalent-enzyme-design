"""Audit completed LigandMPNN focused-redesign pilot FASTA outputs."""

from __future__ import annotations

import json
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from .nyl12_sequence_protection import (
    MINIMUM_PROTECTED_ALPHA,
    MINIMUM_PROTECTED_BETA_LOCAL,
)

LOGGER = logging.getLogger(__name__)
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
THREE_TO_ONE = {
    three: one
    for three, one in zip(
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
        "A R N D C Q E G H I L K M F P S T W Y V".split(),
        strict=True,
    )
}
EXPECTED_CHAINS = {**dict.fromkeys("ABCD", 232), **dict.fromkeys("EFGH", 96)}
REQUIRED_JOBS = 4
SAMPLES_PER_JOB = 6
EXPECTED_SAMPLES = REQUIRED_JOBS * SAMPLES_PER_JOB


def _fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for line_number, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
            if not header:
                raise ValueError(f"Empty FASTA header in {path}:{line_number}")
        elif header is None:
            raise ValueError(f"Sequence precedes FASTA header in {path}:{line_number}")
        else:
            sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def _pdb_sequences(path: Path) -> dict[str, str]:
    residues: dict[str, dict[int, str]] = defaultdict(dict)
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if line[:6] != "ATOM  " or line[12:16].strip() != "CA":
            continue
        chain = line[21]
        if chain not in EXPECTED_CHAINS:
            raise ValueError(
                f"Unexpected protein chain {chain!r} in {path}:{line_number}"
            )
        try:
            position = int(line[22:26])
        except ValueError as exc:
            raise ValueError(
                f"Invalid CA residue number in {path}:{line_number}"
            ) from exc
        residue_name = line[17:20].strip()
        if residue_name not in THREE_TO_ONE:
            raise ValueError(
                f"Nonstandard CA residue {residue_name!r} in {path}:{line_number}"
            )
        if line[16].strip() or line[26].strip():
            raise ValueError(f"Alternate CA location in {path}:{line_number}")
        if position in residues[chain]:
            raise ValueError(f"Duplicate CA {chain}:{position} in {path}")
        residues[chain][position] = THREE_TO_ONE[residue_name]
    sequences: dict[str, str] = {}
    for chain, length in EXPECTED_CHAINS.items():
        positions = residues.get(chain, {})
        if sorted(positions) != list(range(1, length + 1)):
            raise ValueError(f"Input PDB {path} has incomplete {chain} CA sequence")
        sequences[chain] = "".join(
            positions[position] for position in range(1, length + 1)
        )
    return sequences


def _argument(argv: list[str], flag: str, job_id: str) -> str:
    try:
        index = argv.index(flag)
        return argv[index + 1]
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Job {job_id} is missing {flag} in argv") from exc


def _header_value(header: str, key: str, job_id: str) -> str:
    match = re.search(rf"(?:^|,\s*){re.escape(key)}=([^,\s]+)", header)
    if match is None:
        raise ValueError(f"Job {job_id} FASTA header is missing {key}")
    return match.group(1)


def _split_sequence(sequence: str, job_id: str, record_index: int) -> dict[str, str]:
    pieces = sequence.split(":")
    if len(pieces) != len(EXPECTED_CHAINS):
        raise ValueError(
            f"Job {job_id} record {record_index} must contain eight chains"
        )
    chains = dict(zip(EXPECTED_CHAINS, pieces, strict=True))
    for chain, expected_length in EXPECTED_CHAINS.items():
        value = chains[chain]
        if len(value) != expected_length:
            raise ValueError(
                f"Job {job_id} record {record_index} chain {chain} length "
                f"{len(value)}; expected {expected_length}"
            )
        invalid = sorted(set(value) - STANDARD_AA)
        if invalid:
            raise ValueError(
                f"Job {job_id} record {record_index} chain {chain} has "
                f"noncanonical amino acids {invalid}"
            )
    return chains


def _audit_ties_and_changes(
    chains: dict[str, str],
    native: dict[str, str],
    alpha_design: set[int],
    beta_design: set[int],
    protected_alpha: set[int],
    protected_beta: set[int],
    job_id: str,
    sample_id: int,
) -> tuple[list[dict[str, Any]], str, str]:
    protected_alpha = set(protected_alpha) | MINIMUM_PROTECTED_ALPHA
    protected_beta = set(protected_beta) | MINIMUM_PROTECTED_BETA_LOCAL
    alpha, beta = chains["A"], chains["E"]
    if any(chains[chain] != alpha for chain in "BCD"):
        raise ValueError(f"Job {job_id} sample {sample_id} violates alpha A-D ties")
    if any(chains[chain] != beta for chain in "FGH"):
        raise ValueError(f"Job {job_id} sample {sample_id} violates beta E-H ties")
    mutations: list[dict[str, Any]] = []
    for family, designed, protected, designed_positions, offset in (
        ("alpha", alpha, protected_alpha, alpha_design, 0),
        ("beta", beta, protected_beta, beta_design, 232),
    ):
        original = native["A" if family == "alpha" else "E"]
        for position, (before, after) in enumerate(
            zip(original, designed, strict=True), 1
        ):
            if position in protected and after != before:
                raise ValueError(
                    f"Job {job_id} sample {sample_id} changes protected {family} position {position}"
                )
            if position not in designed_positions and after != before:
                raise ValueError(
                    f"Job {job_id} sample {sample_id} changes untied/fixed "
                    f"{family} position {position}"
                )
            if after != before:
                mutations.append(
                    {
                        "family": family,
                        "physical_position": position,
                        "merged_position": position + offset,
                        "from": before,
                        "to": after,
                    }
                )
    return mutations, alpha, beta


def audit_job(manifest: dict[str, Any], job: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate one completed job FASTA and return its six audited designs.

    Args:
        manifest: Parsed preparation manifest containing the design mask and config.
        job: One frozen job entry from ``manifest["jobs"]``.

    Returns:
        Six sample records with parent model, settings, sequences, and mutations.

    Raises:
        FileNotFoundError: The job FASTA or input PDB is absent.
        ValueError: FASTA, native sequence, ties, or fixed-site checks fail.
    """
    mask = manifest["design_mask"]
    alpha_design = set(mask["alpha_positions"])
    beta_design = set(mask["beta_positions"])
    config = manifest["config"]
    protected_alpha = set(config["protected_alpha_positions"])
    protected_beta = set(config["protected_beta_positions"])
    job_id, argv = job["job_id"], job["argv"]
    out_folder = Path(_argument(argv, "--out_folder", job_id))
    input_pdb = Path(_argument(argv, "--pdb_path", job_id))
    temperature = float(_argument(argv, "--temperature", job_id))
    seed = int(_argument(argv, "--seed", job_id))
    batch_size = int(_argument(argv, "--batch_size", job_id))
    batches = int(_argument(argv, "--number_of_batches", job_id))
    if batch_size * batches != SAMPLES_PER_JOB:
        raise ValueError(
            f"Job {job_id} requests {batch_size * batches} samples; "
            f"expected {SAMPLES_PER_JOB}"
        )
    native = _pdb_sequences(input_pdb)
    fasta_path = out_folder / "seqs" / "ligandmpnn_input.fa"
    records = _fasta(fasta_path)
    if len(records) != SAMPLES_PER_JOB + 1:
        raise ValueError(
            f"Job {job_id} FASTA has {len(records) - 1} generated records; "
            f"expected {SAMPLES_PER_JOB}"
        )
    native_chains = _split_sequence(records[0][1], job_id, 0)
    if native_chains != native:
        raise ValueError(f"Job {job_id} native FASTA sequence disagrees with input PDB")
    if any(native[chain] != native["A"] for chain in "BCD"):
        raise ValueError(f"Input PDB {input_pdb} alpha chains are not identical")
    if any(native[chain] != native["E"] for chain in "FGH"):
        raise ValueError(f"Input PDB {input_pdb} beta chains are not identical")

    samples: list[dict[str, Any]] = []
    observed_sample_ids: list[int] = []
    for record_index, (header, sequence) in enumerate(records[1:], start=1):
        chains = _split_sequence(sequence, job_id, record_index)
        header_seed = int(_header_value(header, "seed", job_id))
        header_temperature = float(_header_value(header, "T", job_id))
        sample_id = int(_header_value(header, "id", job_id))
        if header_seed != seed or header_temperature != temperature:
            raise ValueError(f"Job {job_id} FASTA seed/temperature mismatches manifest")
        if not 1 <= sample_id <= SAMPLES_PER_JOB:
            raise ValueError(
                f"Job {job_id} has invalid generated sample id {sample_id}"
            )
        observed_sample_ids.append(sample_id)
        mutations, alpha, beta = _audit_ties_and_changes(
            chains,
            native,
            alpha_design,
            beta_design,
            protected_alpha,
            protected_beta,
            job_id,
            sample_id,
        )
        samples.append(
            {
                "sample_id": f"{job_id}:{sample_id}",
                "parent_model": job_id.rsplit("_t", 1)[0],
                "job_id": job_id,
                "temperature": temperature,
                "seed": seed,
                "header": header,
                "alpha_sequence": alpha,
                "beta_sequence": beta,
                "mutations": mutations,
                "mutation_count": len(mutations),
            }
        )
    if sorted(observed_sample_ids) != list(range(1, SAMPLES_PER_JOB + 1)):
        raise ValueError(f"Job {job_id} has missing or duplicate generated sample IDs")
    LOGGER.info("Audited %d LigandMPNN designs for job %s", len(samples), job_id)
    return samples


def audit_campaign(manifest_path: Path) -> dict[str, Any]:
    """Audit the frozen 24-sequence LigandMPNN pilot.

    Args:
        manifest_path: Preparation manifest with four frozen job argv vectors.

    Returns:
        Per-sample sequences, mutations, duplicate lineage, and summary.

    Raises:
        FileNotFoundError: A manifest, job FASTA, or input PDB is absent.
        ValueError: Job count, FASTA, native sequence, ties, or fixed-site checks fail.
    """
    manifest = json.loads(manifest_path.read_text())
    jobs = manifest["jobs"]
    if len(jobs) != REQUIRED_JOBS:
        raise ValueError(f"Expected {REQUIRED_JOBS} jobs; found {len(jobs)}")
    samples = [sample for job in jobs for sample in audit_job(manifest, job)]

    if len(samples) != EXPECTED_SAMPLES:
        raise ValueError(
            f"Expected {EXPECTED_SAMPLES} generated samples; found {len(samples)}"
        )
    lineage: dict[str, list[str]] = defaultdict(list)
    for sample in samples:
        sequence = sample["alpha_sequence"] + ":" + sample["beta_sequence"]
        lineage[sequence].append(sample["sample_id"])
    duplicate_groups = [
        {"sequence": sequence, "first_sample": ids[0], "duplicate_samples": ids[1:]}
        for sequence, ids in lineage.items()
        if len(ids) > 1
    ]
    mutation_counts = [sample["mutation_count"] for sample in samples]
    result = {
        "manifest": str(manifest_path.resolve()),
        "samples": samples,
        "duplicate_groups": duplicate_groups,
        "summary": {
            "job_count": len(jobs),
            "generated_sample_count": len(samples),
            "native_records_excluded": len(jobs),
            "unique_sequence_count": len(lineage),
            "duplicate_sample_count": len(samples) - len(lineage),
            "mutation_count_min": min(mutation_counts),
            "mutation_count_max": max(mutation_counts),
        },
    }
    LOGGER.info("Audited %d LigandMPNN designs across %d jobs", len(samples), len(jobs))
    return result
