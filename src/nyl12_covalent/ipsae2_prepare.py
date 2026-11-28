"""Prepare a bounded completed native Boltz campaign for optional ipSAE2."""

from __future__ import annotations

import importlib.util
import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from .ipsae2_scoring import file_sha256, require_qualification, verify_runtime
from .ipsae2_sidecar import read_json

LOGGER = logging.getLogger(__name__)


def _write(path: Path, value: Any) -> None:
    with path.open("x") as stream:
        stream.write(
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )


def physical_pairing(cif: Path, adapter: Path) -> dict[str, Any]:
    """Reuse the established protein-only core-contact assignment adapter."""
    specification = importlib.util.spec_from_file_location(
        "ipsae2_pairing_adapter", adapter
    )
    if specification is None or specification.loader is None:
        raise ValueError("Cannot load physical pairing adapter")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    atoms = {
        (site.chain, site.sequence, site.name): (site.element, site.xyz)
        for site in module.read_atom_sites(cif)
    }
    pairs, evidence = module.pairing(atoms, module.reference_atoms())
    return {
        "physical_pairs": pairs,
        "protomer_groups": [[alpha, beta] for alpha, beta in pairs.items()],
        "physical_pairs_verified": True,
        "evidence": evidence,
        "input_sha256": {
            str(path): file_sha256(path) for path in (cif, adapter, module.REFERENCE)
        },
    }


def prepare_completed_campaign(
    *,
    campaign: Path,
    output: Path,
    runtime_path: Path,
    qualification_path: Path,
    pairing_adapter: Path,
    boltz_python: Path,
    boltz_source_root: Path,
    root: Path,
) -> Path:
    """Create explicit per-case manifests, native mapping and scoring config.

    All six native receipts must already be completed. This wrapper does not
    launch predictions, score models, canonicalize coordinates or modify Amber.
    Run only under the parent agent's intent/outcome journal transaction.
    """
    runtime = read_json(runtime_path)
    verify_runtime(runtime)
    require_qualification(read_json(qualification_path), runtime)
    native_paths = sorted(
        (campaign / "prediction/cases").glob("*/native_prediction.json")
    )
    if len(native_paths) != 6:
        raise ValueError("Expected exactly six native prediction receipts")
    for path in native_paths:
        if read_json(path).get("status") != "native_prediction_completed":
            raise ValueError(f"Native prediction is not completed: {path}")
    output.mkdir(parents=True, exist_ok=False)
    cases = []
    for native_path in native_paths:
        native = read_json(native_path)
        case_id = native_path.parent.name
        if native["case"] != case_id:
            raise ValueError("Native receipt case identity mismatch")
        cif = Path(native["prediction_cif"])
        if file_sha256(cif) != native["prediction_cif_sha256"]:
            raise ValueError(f"Native CIF hash mismatch: {case_id}")
        prediction = cif.parent
        model_stem = cif.stem
        results = prediction.parent.parent
        artifacts = {
            "cif": cif,
            "pae": prediction / f"pae_{model_stem}.npz",
            "confidence": prediction / f"confidence_{model_stem}.json",
            "plddt": prediction / f"plddt_{model_stem}.npz",
            "structure": results / "processed/structures" / f"{case_id}.npz",
        }
        model = {
            "case_id": case_id,
            "candidate": native["candidate"],
            "seed": native["seed"],
            "native_receipt": {
                "path": str(native_path),
                "sha256": file_sha256(native_path),
            },
            "artifacts": {
                key: {"path": str(path), "sha256": file_sha256(path)}
                for key, path in artifacts.items()
            },
        }
        case_output = output / "prepared_cases" / case_id
        case_output.mkdir(parents=True)
        model_path = case_output / "model_manifest.json"
        _write(model_path, model)
        mapping_path = case_output / "native_mapping.json"
        command = [
            str(boltz_python),
            str(root / "scripts/verify_ipsae2_native_mapping.py"),
        ]
        for key, path in artifacts.items():
            command.extend([f"--{key}", str(path)])
        command.extend(
            [
                "--boltz-source-root",
                str(boltz_source_root),
                "--output",
                str(mapping_path),
            ]
        )
        with (case_output / "mapping.log").open("x") as log:
            process = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        _write(
            case_output / "mapping_process.json",
            {"command": command, "returncode": process.returncode},
        )
        if process.returncode:
            raise RuntimeError(
                f"Native mapping failed: {case_id}; retained {case_output}"
            )
        mapping = read_json(mapping_path)
        if (
            mapping["accepted"] is not True
            or mapping["artifacts"] != model["artifacts"]
        ):
            raise ValueError(f"Native mapping not bound to triplet: {case_id}")
        pairing = physical_pairing(cif, pairing_adapter)
        pairing_path = case_output / "physical_pairing.json"
        _write(pairing_path, pairing)
        cases.append(
            {
                "case_id": case_id,
                "native_receipt": str(native_path),
                "native_receipt_sha256": file_sha256(native_path),
                "model_manifest": str(model_path),
                "mapping_receipt": str(mapping_path),
                "protomer_groups": pairing["protomer_groups"],
                "physical_pairs_verified": True,
                "physical_pairing_receipt": str(pairing_path),
                "physical_pairing_sha256": file_sha256(pairing_path),
            }
        )
        LOGGER.info("Prepared native mapping and physical groups: %s", case_id)
    config = {
        "schema_version": 1,
        "runtime": str(runtime_path),
        "qualification": str(qualification_path),
        "scoring_root": str(output),
        "scope": "Supplementary native confidence only; no catalytic or PA6 score threshold",
        "journal_command": [
            str(root / ".venv/bin/python"),
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "ipsae2_sidecar",
        ],
        "cases": cases,
    }
    pins = [
        runtime_path,
        qualification_path,
        root / "src/nyl12_covalent/ipsae2_scoring.py",
        root / "src/nyl12_covalent/ipsae2_sidecar.py",
        root / "scripts/score_boltz_ipsae2.py",
    ]
    for case in cases:
        pins.extend(
            Path(case[key])
            for key in (
                "native_receipt",
                "model_manifest",
                "mapping_receipt",
                "physical_pairing_receipt",
            )
        )
    config["input_sha256"] = {str(path): file_sha256(path) for path in pins}
    config_path = output / "sidecar_config.json"
    _write(config_path, config)
    return config_path
