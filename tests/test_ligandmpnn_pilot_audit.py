"""Focused validation tests for LigandMPNN pilot FASTA outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.nyl12_covalent.ligandmpnn_pilot_audit import audit_campaign, audit_job

THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def _native_sequences() -> dict[str, str]:
    alpha = list("A" * 232)
    for position, amino_acid in {113: "F", 186: "N", 229: "E", 232: "N"}.items():
        alpha[position - 1] = amino_acid
    beta = list("A" * 96)
    for position, amino_acid in {1: "T", 40: "D", 42: "D"}.items():
        beta[position - 1] = amino_acid
    return {
        **dict.fromkeys("ABCD", "".join(alpha)),
        **dict.fromkeys("EFGH", "".join(beta)),
    }


def _write_pdb(path: Path, sequences: dict[str, str]) -> None:
    lines = []
    serial = 1
    for chain, sequence in sequences.items():
        for position, amino_acid in enumerate(sequence, 1):
            lines.append(
                f"ATOM  {serial:5d}   CA {THREE[amino_acid]} {chain}{position:4d}    "
                "   0.000   0.000   0.000  1.00 80.00           C"
            )
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")


def _sequence(alpha: str, beta: str) -> str:
    return ":".join([alpha] * 4 + [beta] * 4)


def _campaign(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    native = _native_sequences()
    alpha, beta = native["A"], native["E"]
    alpha_design = alpha[:49] + "S" + alpha[50:]
    beta_design = beta[:9] + "V" + beta[10:]
    manifest_jobs = []
    fasta_paths: dict[str, Path] = {}
    for model in (40, 41):
        for temperature in (0.1, 0.2):
            job_id = f"model_{model}_t{temperature}"
            out_folder = tmp_path / "outputs" / job_id
            fasta_path = out_folder / "seqs" / "ligandmpnn_input.fa"
            fasta_path.parent.mkdir(parents=True)
            input_pdb = tmp_path / f"{job_id}.pdb"
            _write_pdb(input_pdb, native)
            seed = 20261100 + model + int(temperature * 10)
            native_header = f"ligandmpnn_input, T={temperature}, seed={seed}"
            records = [f">{native_header}\n{_sequence(alpha, beta)}"]
            records.extend(
                f">ligandmpnn_input, id={sample_id}, T={temperature}, seed={seed}\n"
                f"{_sequence(alpha_design, beta_design)}"
                for sample_id in range(1, 7)
            )
            fasta_path.write_text("\n".join(records) + "\n")
            manifest_jobs.append(
                {
                    "job_id": job_id,
                    "argv": [
                        "python",
                        "run.py",
                        "--out_folder",
                        str(out_folder),
                        "--pdb_path",
                        str(input_pdb),
                        "--temperature",
                        str(temperature),
                        "--seed",
                        str(seed),
                        "--batch_size",
                        "1",
                        "--number_of_batches",
                        "6",
                    ],
                }
            )
            fasta_paths[job_id] = fasta_path
    manifest = {
        "config": {
            "protected_alpha_positions": [113, 186, *range(222, 233)],
            "protected_beta_positions": [1, 2, 3, 40, 42],
        },
        "design_mask": {"alpha_positions": [50], "beta_positions": [10]},
        "jobs": manifest_jobs,
    }
    manifest_path = tmp_path / "preparation_manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    return manifest_path, fasta_paths


def test_audit_excludes_native_counts_samples_and_reports_duplicate_lineage(
    tmp_path: Path,
) -> None:
    manifest, _ = _campaign(tmp_path)
    manifest_data = json.loads(manifest.read_text())
    assert len(audit_job(manifest_data, manifest_data["jobs"][0])) == 6
    result = audit_campaign(manifest)
    assert result["summary"] == {
        "job_count": 4,
        "generated_sample_count": 24,
        "native_records_excluded": 4,
        "unique_sequence_count": 1,
        "duplicate_sample_count": 23,
        "mutation_count_min": 2,
        "mutation_count_max": 2,
    }
    assert len(result["samples"]) == 24
    assert result["samples"][0]["mutations"] == [
        {
            "family": "alpha",
            "physical_position": 50,
            "merged_position": 50,
            "from": "A",
            "to": "S",
        },
        {
            "family": "beta",
            "physical_position": 10,
            "merged_position": 242,
            "from": "A",
            "to": "V",
        },
    ]
    assert len(result["duplicate_groups"][0]["duplicate_samples"]) == 23
    assert result["samples"][0]["seed"] == 20261141


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ("untied", "untied/fixed alpha position 20"),
        ("tie", "violates alpha A-D ties"),
        ("protected", "changes protected alpha position 113"),
        ("noncanonical", "noncanonical amino acids"),
        ("count", "has 5 generated records"),
    ],
)
def test_audit_rejects_invalid_generated_fasta(
    tmp_path: Path,
    change: str,
    message: str,
) -> None:
    manifest, fasta_paths = _campaign(tmp_path)
    path = next(iter(fasta_paths.values()))
    records = path.read_text().splitlines()
    if change == "count":
        records = records[:-2]
    else:
        sequence_index = 3
        chains = records[sequence_index].split(":")
        if change == "untied":
            chains[0] = chains[0][:19] + "C" + chains[0][20:]
            chains[1:4] = [chains[0]] * 3
        elif change == "tie":
            chains[1] = chains[1][:49] + "C" + chains[1][50:]
        elif change == "protected":
            chains[0] = chains[0][:112] + "C" + chains[0][113:]
            chains[1:4] = [chains[0]] * 3
        elif change == "noncanonical":
            chains[0] = chains[0][:49] + "X" + chains[0][50:]
            chains[1:4] = [chains[0]] * 3
        records[sequence_index] = ":".join(chains)
    path.write_text("\n".join(records) + "\n")
    with pytest.raises(ValueError, match=message):
        audit_campaign(manifest)


def test_audit_rejects_fasta_native_sequence_mismatch(tmp_path: Path) -> None:
    manifest, fasta_paths = _campaign(tmp_path)
    path = next(iter(fasta_paths.values()))
    records = path.read_text().splitlines()
    chains = records[1].split(":")
    chains[0] = chains[0][:112] + "A" + chains[0][113:]
    records[1] = ":".join(chains)
    path.write_text("\n".join(records) + "\n")
    with pytest.raises(ValueError, match="native FASTA sequence disagrees"):
        audit_campaign(manifest)
