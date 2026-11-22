"""Prepare sequence-specific, template-conditioned Boltz input manifests."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

from nyl12_covalent.ligandmpnn_pilot_audit import THREE_TO_ONE, EXPECTED_CHAINS

LOGGER = logging.getLogger(__name__)
LIGAND_CHAINS = "IJKL"


def _native_sequences(path: Path) -> dict[str, str]:
    residues: dict[str, dict[int, str]] = {chain: {} for chain in EXPECTED_CHAINS}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if line[:6] != "ATOM  " or line[12:16].strip() != "CA":
            continue
        chain = line[21]
        if chain not in residues:
            raise ValueError(f"Unexpected CA chain {chain!r} at {path}:{line_number}")
        try:
            position = int(line[22:26])
        except ValueError as exc:
            raise ValueError(
                f"Invalid CA residue number at {path}:{line_number}"
            ) from exc
        residue = line[17:20].strip()
        if residue not in THREE_TO_ONE or position in residues[chain]:
            raise ValueError(f"Invalid or duplicate CA residue at {path}:{line_number}")
        residues[chain][position] = THREE_TO_ONE[residue]
    result: dict[str, str] = {}
    for chain, length in EXPECTED_CHAINS.items():
        if sorted(residues[chain]) != list(range(1, length + 1)):
            raise ValueError(f"Native PDB has incomplete {chain} CA sequence")
        result[chain] = "".join(residues[chain][i] for i in range(1, length + 1))
    if any(result[chain] != result["A"] for chain in "BCD"):
        raise ValueError("Native alpha chains A-D are not identical")
    if any(result[chain] != result["E"] for chain in "FGH"):
        raise ValueError("Native beta chains E-H are not identical")
    return result


def _read_msa(
    path: Path, native_sequence: str
) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if fields is None or not {"key", "sequence"}.issubset(fields):
        raise ValueError(f"MSA {path} must contain key and sequence columns")
    query_rows = [row for row in rows if row["key"] == "0"]
    if len(query_rows) != 1 or query_rows[0]["sequence"] != native_sequence:
        raise ValueError(f"MSA {path} must have one key=0 row matching native sequence")
    return fields, rows


def _write_msa(
    source: Path, destination: Path, native_sequence: str, sequence: str
) -> None:
    fields, rows = _read_msa(source, native_sequence)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            changed = dict(row)
            if row["key"] == "0":
                changed["sequence"] = sequence
            writer.writerow(changed)
    _, written_rows = _read_msa(destination, sequence)
    original_nonquery = [row for row in rows if row["key"] != "0"]
    written_nonquery = [row for row in written_rows if row["key"] != "0"]
    if original_nonquery != written_nonquery:
        raise ValueError(f"Non-query MSA rows changed while writing {destination}")


def _mutations(
    alpha: str,
    beta: str,
    native_alpha: str,
    native_beta: str,
    alpha_positions: set[int],
    beta_positions: set[int],
    protected_alpha: set[int],
    protected_beta: set[int],
) -> list[dict[str, Any]]:
    mutations: list[dict[str, Any]] = []
    for family, sequence, native, designed, protected, offset in (
        ("alpha", alpha, native_alpha, alpha_positions, protected_alpha, 0),
        ("beta", beta, native_beta, beta_positions, protected_beta, 232),
    ):
        if len(sequence) != len(native):
            raise ValueError(f"{family} sequence length differs from native")
        for position, (before, after) in enumerate(
            zip(native, sequence, strict=True), 1
        ):
            if before == after:
                continue
            if position in protected:
                raise ValueError(
                    f"Candidate changes protected {family} position {position}"
                )
            if position not in designed:
                raise ValueError(
                    f"Candidate changes untied/fixed {family} position {position}"
                )
            mutations.append(
                {
                    "family": family,
                    "physical_position": position,
                    "merged_position": position + offset,
                    "from": before,
                    "to": after,
                }
            )
    return mutations


def _prepare_yaml(
    template: dict[str, Any],
    alpha: str,
    beta: str,
    alpha_msa: Path,
    beta_msa: Path,
    ligand_count: int,
) -> dict[str, Any]:
    data = json.loads(json.dumps(template))
    proteins = [entry["protein"] for entry in data["sequences"] if "protein" in entry]
    ligands = [entry["ligand"] for entry in data["sequences"] if "ligand" in entry]
    if len(proteins) != 2 or len(ligands) != 1:
        raise ValueError("Template YAML must contain two proteins and one ligand")
    proteins[0]["sequence"], proteins[0]["msa"] = alpha, str(alpha_msa.resolve())
    proteins[1]["sequence"], proteins[1]["msa"] = beta, str(beta_msa.resolve())
    ligand_ids = list(LIGAND_CHAINS[:ligand_count])
    ligands[0]["id"] = ligand_ids[0] if ligand_count == 1 else ligand_ids
    if ligand_count == 1:
        data["properties"] = []
    else:
        contacts = [entry["contact"] for entry in data.get("constraints", [])]
        if len(contacts) != 2:
            raise ValueError("Four-ligand template requires the two original contacts")
        data["constraints"] = []
        for alpha_chain, beta_chain, ligand_chain in zip(
            "ABCD", "EFGH", LIGAND_CHAINS, strict=True
        ):
            for contact in contacts:
                mapped = beta_chain if contact["token1"][0] == "E" else alpha_chain
                data["constraints"].append(
                    {
                        "contact": {
                            "token1": [mapped, contact["token1"][1]],
                            "max_distance": contact["max_distance"],
                            "force": contact["force"],
                            "token2": [ligand_chain, contact["token2"][1]],
                        }
                    }
                )
        data.pop("properties", None)
    return data


def prepare_inputs(
    report_path: Path,
    native_pdb: Path,
    template_json: Path,
    output: Path,
    ligand_count: int,
    seeds: list[int],
) -> dict[str, Any]:
    """Prepare Boltz YAML/MSA inputs and a seeded job manifest.

    Args:
        report_path: Verified LigandMPNN campaign report with 24 samples.
        native_pdb: Accepted Model 23 eight-chain parent coordinates.
        template_json: JSON-formatted Boltz YAML emitted from the validated template.
        output: New, nonexistent directory for prepared inputs.
        ligand_count: Number of identical substrate copies, either one or four.
        seeds: Five distinct nonnegative inference seeds.

    Returns:
        Candidate and per-seed job manifest written to ``candidate_manifest.json``.

    Raises:
        FileNotFoundError: A required input or cached MSA is missing.
        ValueError: Sequence, schema, count, or MSA validation fails.
        FileExistsError: The output path already exists.
    """
    if ligand_count not in (1, 4):
        raise ValueError("ligand_count must be 1 or 4")
    if len(seeds) != 5 or len(set(seeds)) != 5 or any(seed < 0 for seed in seeds):
        raise ValueError("Exactly five distinct nonnegative seeds are required")
    if output.exists():
        raise FileExistsError(f"Output directory already exists: {output}")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    samples = report["samples"]
    if len(samples) != 24 or report["summary"]["unique_sequence_count"] != 8:
        raise ValueError("LigandMPNN report must contain 24 samples and eight uniques")
    manifest_path = Path(report["manifest"])
    source_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    alpha_positions = set(source_manifest["design_mask"]["alpha_positions"])
    beta_positions = set(source_manifest["design_mask"]["beta_positions"])
    protected_alpha = set(source_manifest["config"]["protected_alpha_positions"])
    protected_beta = set(source_manifest["config"]["protected_beta_positions"])
    native = _native_sequences(native_pdb)
    native_alpha, native_beta = native["A"], native["E"]
    template = json.loads(template_json.read_text(encoding="utf-8"))
    proteins = [
        entry["protein"] for entry in template["sequences"] if "protein" in entry
    ]
    if (
        len(proteins) != 2
        or proteins[0]["sequence"] != native_alpha
        or proteins[1]["sequence"] != native_beta
    ):
        raise ValueError("Template protein sequences do not match accepted native PDB")
    if "templates" not in template or len(template.get("constraints", [])) != 2:
        raise ValueError(
            "Template YAML must retain its existing template and two contacts"
        )
    alpha_msa_source, beta_msa_source = (
        Path(proteins[0]["msa"]),
        Path(proteins[1]["msa"]),
    )
    candidates: list[dict[str, Any]] = [
        {
            "id": "WT",
            "alpha_sequence": native_alpha,
            "beta_sequence": native_beta,
            "mutations": [],
            "mutation_count": 0,
            "lineage": [],
        }
    ]
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for sample in samples:
        key = (sample["alpha_sequence"], sample["beta_sequence"])
        grouped.setdefault(key, []).append(sample)
    if len(grouped) != 8:
        raise ValueError(
            f"Expected eight unique candidate sequences; found {len(grouped)}"
        )
    for index, ((alpha, beta), lineage) in enumerate(grouped.items(), 1):
        mutations = _mutations(
            alpha,
            beta,
            native_alpha,
            native_beta,
            alpha_positions,
            beta_positions,
            protected_alpha,
            protected_beta,
        )
        if any(sample["mutation_count"] != len(mutations) for sample in lineage):
            raise ValueError(
                "Reported mutation count disagrees with candidate sequence"
            )
        candidates.append(
            {
                "id": f"variant_{index:02d}",
                "alpha_sequence": alpha,
                "beta_sequence": beta,
                "mutations": mutations,
                "mutation_count": len(mutations),
                "lineage": [sample["sample_id"] for sample in lineage],
            }
        )
    output.mkdir(parents=True)
    jobs: list[dict[str, Any]] = []
    for candidate_index, candidate in enumerate(candidates):
        candidate_id = candidate["id"]
        msa_dir = output / "msas" / candidate_id
        alpha_msa, beta_msa = msa_dir / "alpha.csv", msa_dir / "beta.csv"
        _write_msa(
            alpha_msa_source, alpha_msa, native_alpha, candidate["alpha_sequence"]
        )
        _write_msa(beta_msa_source, beta_msa, native_beta, candidate["beta_sequence"])
        data = _prepare_yaml(
            template,
            candidate["alpha_sequence"],
            candidate["beta_sequence"],
            alpha_msa,
            beta_msa,
            ligand_count,
        )
        yaml_path = output / "inputs" / f"{candidate_id}.yaml"
        yaml_path.parent.mkdir(parents=True, exist_ok=True)
        # JSON syntax is a valid YAML document and needs no extra runtime package.
        yaml_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        for model_index, seed in enumerate(seeds):
            jobs.append(
                {
                    "variant": candidate_id,
                    "variant_index": candidate_index,
                    "ligand_ids": list(LIGAND_CHAINS[:ligand_count]),
                    "seed": seed,
                    "model_index": model_index,
                    "yaml_path": str(yaml_path.resolve()),
                    "expected_alpha_sequence": candidate["alpha_sequence"],
                    "expected_beta_sequence": candidate["beta_sequence"],
                }
            )
    result: dict[str, Any] = {
        "schema_version": 1,
        "report": str(report_path.resolve()),
        "native_pdb": str(native_pdb.resolve()),
        "template_yaml_json_source": str(template_json.resolve()),
        "ligand_count": ligand_count,
        "seeds": seeds,
        "candidates": candidates,
        "jobs": jobs,
        "calculation_launched": False,
    }
    (output / "candidate_manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info(
        "Prepared %d Boltz candidates and %d one-sample jobs in %s",
        len(candidates),
        len(jobs),
        output,
    )
    return result
