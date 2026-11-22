"""Tests for sequence-specific Boltz input preparation."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest
from nyl12_covalent.design_boltz_inputs import prepare_inputs

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


def _write_pdb(path: Path, alpha: str, beta: str) -> None:
    lines: list[str] = []
    serial = 1
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
    path.write_text("\n".join(lines) + "\nEND\n", encoding="ascii")


def _fixture(root: Path) -> tuple[Path, Path, Path, str, str]:
    motif = "ACDEFGHIKLMNPQRSTVWY"
    alpha = (motif * 12)[:232]
    beta = (motif * 5)[:96]
    native_pdb = root / "native.pdb"
    _write_pdb(native_pdb, alpha, beta)
    alpha_msa, beta_msa = root / "alpha.csv", root / "beta.csv"
    for path, sequence in ((alpha_msa, alpha), (beta_msa, beta)):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["key", "sequence", "note"])
            writer.writeheader()
            writer.writerow({"key": "0", "sequence": sequence, "note": "query"})
            writer.writerow(
                {"key": "1", "sequence": "A" * len(sequence), "note": "hit"}
            )
    template = root / "template.json"
    template.write_text(
        json.dumps(
            {
                "version": 1,
                "sequences": [
                    {
                        "protein": {
                            "id": list("ABCD"),
                            "sequence": alpha,
                            "msa": str(alpha_msa),
                        }
                    },
                    {
                        "protein": {
                            "id": list("EFGH"),
                            "sequence": beta,
                            "msa": str(beta_msa),
                        }
                    },
                    {"ligand": {"id": "I", "smiles": "CC(=O)N"}},
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
                        "cif": "template.cif",
                        "chain_id": list("AEBFCGDH"),
                        "template_id": list("AABBCCDD"),
                    }
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    source_manifest = {
        "design_mask": {
            "alpha_positions": [
                50,
                51,
                52,
                73,
                76,
                101,
                111,
                112,
                114,
                115,
                116,
                151,
                152,
                153,
                185,
                187,
                188,
            ],
            "beta_positions": [39, 41],
        },
        "config": {
            "protected_alpha_positions": [113, 156, 186, *range(222, 233)],
            "protected_beta_positions": [1, 2, 3, 40, 42],
        },
    }
    report_manifest = root / "preparation_manifest.json"
    report_manifest.write_text(json.dumps(source_manifest), encoding="utf-8")
    samples: list[dict[str, Any]] = []
    for index in range(8):
        # Each variant differs at a distinct editable position.
        position = [50, 51, 52, 73, 76, 101, 111, 112][index]
        replacement = next(
            aa for aa in "ACDEFGHIKLMNPQRSTVWY" if aa != alpha[position - 1]
        )
        candidate_alpha = alpha[: position - 1] + replacement + alpha[position:]
        samples.append(
            {
                "sample_id": f"sample-{index}",
                "alpha_sequence": candidate_alpha,
                "beta_sequence": beta,
                "mutation_count": 1,
                "mutations": [
                    {
                        "family": "alpha",
                        "physical_position": position,
                        "merged_position": position,
                        "from": alpha[position - 1],
                        "to": replacement,
                    }
                ],
            }
        )
    report_path = root / "report.json"
    report_path.write_text(
        json.dumps(
            {
                "manifest": str(report_manifest),
                "summary": {"unique_sequence_count": 8},
                "samples": samples * 3,
            }
        ),
        encoding="utf-8",
    )
    return report_path, native_pdb, template, alpha, beta


@pytest.mark.parametrize("ligand_count", [1, 4])
def test_prepare_inputs_writes_nine_conditions_and_45_seeded_jobs(
    tmp_path: Path,
    ligand_count: int,
) -> None:
    report, native, template, _, _ = _fixture(tmp_path)
    output = tmp_path / "prepared"
    result = prepare_inputs(
        report, native, template, output, ligand_count, [1, 2, 3, 4, 5]
    )
    assert [candidate["id"] for candidate in result["candidates"]] == [
        "WT",
        *(f"variant_{index:02d}" for index in range(1, 9)),
    ]
    assert len(result["jobs"]) == 45
    assert result["jobs"][0]["ligand_ids"] == list("IJKL")[:ligand_count]
    assert result["calculation_launched"] is False
    prepared = json.loads((output / "inputs/WT.yaml").read_text())
    assert "templates" in prepared
    assert prepared["sequences"][-1]["ligand"]["id"] == (
        "I" if ligand_count == 1 else list("IJKL")
    )
    if ligand_count == 4:
        assert "properties" not in prepared
        assert len(prepared["constraints"]) == 8
        assert [
            item["contact"]["token1"][0] for item in prepared["constraints"]
        ] == list("EAFBGCHD")
        assert [
            item["contact"]["token2"][0] for item in prepared["constraints"]
        ] == list("IIJJKKLL")
    with (output / "msas/variant_01/alpha.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["key"] == "0"
    assert rows[0]["sequence"] == result["candidates"][1]["alpha_sequence"]
    assert rows[1] == {"key": "1", "sequence": rows[1]["sequence"], "note": "hit"}


def test_prepare_inputs_rejects_protected_mutation(tmp_path: Path) -> None:
    report, native, template, alpha, _ = _fixture(tmp_path)
    data = json.loads(report.read_text())
    protected = 113
    replacement = next(
        aa for aa in "ACDEFGHIKLMNPQRSTVWY" if aa != alpha[protected - 1]
    )
    bad_alpha = data["samples"][0]["alpha_sequence"]
    bad_alpha = bad_alpha[: protected - 1] + replacement + bad_alpha[protected:]
    original_candidate = data["samples"][0]["alpha_sequence"]
    for sample in data["samples"]:
        if sample["alpha_sequence"] == original_candidate:
            sample["alpha_sequence"] = bad_alpha
            sample["mutation_count"] = 2
    report.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="protected alpha position 113"):
        prepare_inputs(report, native, template, tmp_path / "bad", 4, [1, 2, 3, 4, 5])


def test_prepare_inputs_rejects_existing_output_and_bad_seed_count(
    tmp_path: Path,
) -> None:
    report, native, template, _, _ = _fixture(tmp_path)
    with pytest.raises(ValueError, match="five distinct"):
        prepare_inputs(report, native, template, tmp_path / "out", 4, [1, 2])
    existing = tmp_path / "exists"
    existing.mkdir()
    with pytest.raises(FileExistsError):
        prepare_inputs(report, native, template, existing, 4, [1, 2, 3, 4, 5])
