"""Tests for deterministic mask selection and matched Boltz preparation."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest
from nyl12_covalent.design_boltz_inputs import (
    _native_sequences,
    _prepare_yaml,
    _write_msa,
)
from nyl12_covalent.design_mask_candidates import (
    SEEDS,
    prepare_validation,
    select_candidates,
)

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
AA = "ACDEFGHIKLMNPQRSTVWY"


def _write_pdb(path: Path, alpha: str, beta: str) -> None:
    lines, serial = [], 1
    for chain in "ABCD":
        for position, aa in enumerate(alpha, 1):
            lines.append(
                f"ATOM  {serial:5d}   CA {THREE[aa]} {chain}{position:4d}    "
                "   0.000   0.000   0.000  1.00 80.00           C"
            )
            serial += 1
    for chain in "EFGH":
        for position, aa in enumerate(beta, 1):
            lines.append(
                f"ATOM  {serial:5d}   CA {THREE[aa]} {chain}{position:4d}    "
                "   0.000   0.000   0.000  1.00 80.00           C"
            )
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")


def _select_samples() -> list[dict[str, Any]]:
    rows = []
    for arm in ("focused", "broader"):
        for parent in ("p1", "p2"):
            for temperature in (0.1, 0.2):
                index = len(rows)
                sequence = (
                    "C" + "A" * 7
                    if parent == "p1" and temperature == 0.1
                    else AA[(index % len(AA))] * 8
                )
                rows.append(
                    {
                        "sample_id": f"{arm}-{parent}-{temperature}",
                        "arm": arm,
                        "parent": parent,
                        "temperature": temperature,
                        "alpha_sequence": sequence,
                        "beta_sequence": "A" * 4,
                    }
                )
    # Both arms choose this same first-stratum sequence, preserving both lineages.
    return rows


def _previous_campaign(root: Path) -> tuple[Path, str, str, Path, Path]:
    previous = root / "previous"
    (previous / "configuration").mkdir(parents=True)
    (previous / "prepared/inputs").mkdir(parents=True)
    (previous / "prepared/msas/WT").mkdir(parents=True)
    source_manifest = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "design/nyl12_ligandmpnn_pilot/boltz_validation/prepared/candidate_manifest.json"
        ).read_text()
    )
    alpha, beta = (
        _native_sequences(Path(source_manifest["native_pdb"]))["A"],
        _native_sequences(Path(source_manifest["native_pdb"]))["E"],
    )
    native = previous / "native.pdb"
    _write_pdb(native, alpha, beta)
    alpha_msa, beta_msa = (
        previous / "configuration/alpha.csv",
        previous / "configuration/beta.csv",
    )
    for path, sequence in ((alpha_msa, alpha), (beta_msa, beta)):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["key", "sequence", "note"])
            writer.writeheader()
            writer.writerow({"key": "0", "sequence": sequence, "note": "query"})
            writer.writerow(
                {"key": "1", "sequence": "A" * len(sequence), "note": "homolog"}
            )
    cif = previous / "configuration/n12c.cif"
    cif.write_text("fixture CIF\n")
    template = {
        "version": 1,
        "sequences": [
            {"protein": {"id": list("ABCD"), "sequence": alpha, "msa": str(alpha_msa)}},
            {"protein": {"id": list("EFGH"), "sequence": beta, "msa": str(beta_msa)}},
            {
                "ligand": {
                    "id": "I",
                    "smiles": "[NH3+]CCCCCC(NCCCCCC(NCCCCCC(NCCCCCC([O-])=O)=O)=O)=O",
                }
            },
        ],
        "constraints": [
            {
                "contact": {
                    "token1": ["E", 1],
                    "max_distance": 3.0,
                    "force": True,
                    "token2": ["I", "C54"],
                }
            },
            {
                "contact": {
                    "token1": ["A", 186],
                    "max_distance": 3.2,
                    "force": True,
                    "token2": ["I", "O49"],
                }
            },
        ],
        "properties": [{"affinity": {"binder": "I"}}],
        "templates": [
            {
                "cif": str(cif),
                "chain_id": list("AEBFCGDH"),
                "template_id": list("AABBCCDD"),
            }
        ],
    }
    template_path = previous / "configuration/template.json"
    template_path.write_text(json.dumps(template))
    old_alpha, old_beta = (
        previous / "prepared/msas/WT/alpha.csv",
        previous / "prepared/msas/WT/beta.csv",
    )
    _write_msa(alpha_msa, old_alpha, alpha, alpha)
    _write_msa(beta_msa, old_beta, beta, beta)
    old_yaml = previous / "prepared/inputs/WT.yaml"
    old_yaml.write_text(
        json.dumps(_prepare_yaml(template, alpha, beta, old_alpha, old_beta, 4))
    )
    release = previous / "parameters"
    release.mkdir()
    mol2, mapping = (
        previous / "configuration/whole_l4_source.mol2",
        previous / "configuration/whole_l4_context_mapping.json",
    )
    mol2.write_text("mol2 fixture\n")
    mapping.write_text("{}\n")
    relaxation = {
        "release": str(release),
        "mol2": str(mol2),
        "mapping": str(mapping),
        "reference_cif": str(cif),
    }
    (previous / "configuration/relaxation.json").write_text(json.dumps(relaxation))
    jobs, hosts = [], []
    for index, seed in enumerate(SEEDS):
        jobs.append(
            {
                "variant": "WT",
                "model_index": index,
                "seed": seed,
                "yaml_path": str(old_yaml),
                "expected_alpha_sequence": alpha,
                "expected_beta_sequence": beta,
            }
        )
        hosts.append(
            {
                "variant": "WT",
                "model_index": index,
                "seed": seed,
                "prediction_host": "violet" if index in (0, 2, 4) else "tulip",
            }
        )
        if index in (0, 2, 4):
            base = previous / "models/WT" / f"seed_{index + 1:02d}"
            (base / "raw").mkdir(parents=True)
            (base / "relaxation").mkdir()
            raw, relaxed = base / "raw/normalized.pdb", base / "relaxation/relaxed.pdb"
            raw.write_text(f"raw WT {seed}\n")
            relaxed.write_text(f"relaxed WT {seed}\n")
            (base / "request.json").write_text(
                json.dumps(
                    {
                        "variant": "WT",
                        "seed": seed,
                        "model_index": index,
                        "ligand_count": 4,
                        "sha256": hashlib.sha256(raw.read_bytes()).hexdigest(),
                    }
                )
            )
            (base / "raw/mapping.json").write_text(
                json.dumps(
                    {
                        "status": "accepted_graph_isomorphic_component_normalization",
                        "output_sha256": hashlib.sha256(raw.read_bytes()).hexdigest(),
                    }
                )
            )
            acceptance = {
                "status": "minimization_complete",
                "process_exit_code": 0,
                "heavy_export": {
                    "output_sha256": hashlib.sha256(relaxed.read_bytes()).hexdigest()
                },
            }
            (base / "relaxation/acceptance.json").write_text(json.dumps(acceptance))
    (previous / "prepared/candidate_manifest.json").write_text(
        json.dumps(
            {"native_pdb": str(native), "ligand_count": 4, "seeds": SEEDS, "jobs": jobs}
        )
    )
    (previous / "dispatch").mkdir(exist_ok=True)
    with (previous / "dispatch/prediction_host_assignments.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(hosts[0]))
        writer.writeheader()
        writer.writerows(hosts)
    (previous / "dispatch/runtime_versions.json").write_text(
        json.dumps(
            {
                "same_boltz_source": True,
                "checkpoint_sha256": "abc",
                "violet": {"python": "3.10"},
            }
        )
    )
    return previous, alpha, beta, template_path, cif


def _report(root: Path, previous: tuple[Path, str, str, Path, Path]) -> Path:
    source_root = (
        Path(__file__).resolve().parents[1]
        / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison"
    )
    source_report = json.loads(
        (source_root / "generation/results/report.json").read_text()
    )
    manifest_path = Path(source_report["manifest"])
    sample_rows = source_report["samples"]
    report_path = root / "generation_report.json"
    report_path.write_text(
        json.dumps(
            {
                "manifest": str(manifest_path),
                "source_sha256": source_report["source_sha256"],
                "samples": sample_rows,
            }
        )
    )
    (root / "acceptance.json").write_text(
        json.dumps(
            {
                "accepted": True,
                "report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
            }
        )
    )
    return report_path


def test_select_candidates_is_deterministic_diverse_and_deduplicates_arms() -> None:
    rows = _select_samples()
    first = select_candidates(rows)
    assert first == select_candidates(list(reversed(rows)))
    shared = next(row for row in first if row["alpha_sequence"] == "C" + "A" * 7)
    assert shared["arms"] == ["broader", "focused"]
    assert {lineage["arm"] for lineage in shared["selection_lineage"]} == set(
        shared["arms"]
    )


def test_select_candidates_rejects_missing_fields_and_duplicate_ids() -> None:
    rows = _select_samples()
    del rows[0]["temperature"]
    with pytest.raises(ValueError, match="missing fields"):
        select_candidates(rows)
    rows = _select_samples()
    rows[1]["sample_id"] = rows[0]["sample_id"]
    with pytest.raises(ValueError, match="Duplicate sample_id"):
        select_candidates(rows)


def test_prepare_validation_preserves_identity_controls_and_wt_reuse(
    tmp_path: Path,
) -> None:
    previous_context = _previous_campaign(tmp_path)
    previous = previous_context[0]
    report = _report(tmp_path, previous_context)
    result = prepare_validation(report, previous, tmp_path / "new_validation", SEEDS)
    assert result["ligand_count"] == 4
    assert result["raw_sample_count"] == 48
    assert result["new_prediction_job_count"] <= 52
    assert [row["model_index"] for row in result["reused_models"]] == [0, 2, 4]
    assert {row["seed"] for row in result["reused_models"]} == {
        SEEDS[0],
        SEEDS[2],
        SEEDS[4],
    }
    assert [row["id"] for row in result["candidates"][:3]] == [
        "WT",
        "control_K59D",
        "control_V111C",
    ]
    wt_jobs = [row for row in result["jobs"] if row["variant"] == "WT"]
    assert [(row["model_index"], row["seed"]) for row in wt_jobs] == [
        (1, SEEDS[1]),
        (3, SEEDS[3]),
    ]
    manifest = json.loads(
        (tmp_path / "new_validation/prepared/candidate_manifest.json").read_text()
    )
    assert manifest["jobs"] == result["jobs"]
    with (tmp_path / "new_validation/prepared/jobs.csv").open(newline="") as handle:
        assert len(list(csv.DictReader(handle))) == len(result["jobs"])
    prepared = json.loads(
        (tmp_path / "new_validation/prepared/inputs/WT.yaml").read_text()
    )
    assert prepared["sequences"][-1]["ligand"]["id"] == list("IJKL")
    assert len(prepared["constraints"]) == 8 and "properties" not in prepared
    with (tmp_path / "new_validation/prepared/msas/WT/alpha.csv").open(
        newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["sequence"] == prepared["sequences"][0]["protein"]["sequence"]
    assert rows[1]["note"] == "homolog" and rows[1]["sequence"] == "A" * 232
    assert (tmp_path / "new_validation/configuration/relaxation.json").is_file()
    assert result["runtime_requirement"]["host"] == "violet"


def test_prepare_validation_rejects_seed_or_reuse_identity_mismatch(
    tmp_path: Path,
) -> None:
    previous_context = _previous_campaign(tmp_path)
    previous = previous_context[0]
    report = _report(tmp_path, previous_context)
    with pytest.raises(ValueError, match="fixed five-seed"):
        prepare_validation(report, previous, tmp_path / "bad_seed", [1, 2, 3, 4, 5])
    request = previous / "models/WT/seed_01/request.json"
    value = json.loads(request.read_text())
    value["model_index"] = 4
    request.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="request does not match"):
        prepare_validation(report, previous, tmp_path / "bad_reuse", SEEDS)


def test_prepare_validation_rejects_wt_msa_content_mismatch(tmp_path: Path) -> None:
    previous_context = _previous_campaign(tmp_path)
    previous = previous_context[0]
    report = _report(tmp_path, previous_context)
    old_msa = previous / "prepared/msas/WT/alpha.csv"
    old_msa.write_text(old_msa.read_text().replace("homolog", "changed"))
    with pytest.raises(ValueError, match="scientific identity differs"):
        prepare_validation(report, previous, tmp_path / "bad_msa", SEEDS)


def test_prepare_validation_rejects_wt_reuse_on_wrong_host(tmp_path: Path) -> None:
    previous_context = _previous_campaign(tmp_path)
    previous = previous_context[0]
    report = _report(tmp_path, previous_context)
    assignment_path = previous / "dispatch/prediction_host_assignments.csv"
    with assignment_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["prediction_host"] = "tulip"
    with assignment_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="not predicted on Violet"):
        prepare_validation(report, previous, tmp_path / "bad_host", SEEDS)
