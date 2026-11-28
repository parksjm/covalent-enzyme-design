"""Run stock native CPU feature validation for two template-free paired WT inputs."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import torch
from pytorch_lightning import seed_everything
from boltz.data import const
from boltz.data.feature import featurizerv2 as native
from boltz.data.module.inferencev2 import (
    Boltz2InferenceDataModule,
    PredictionDataset,
    load_input,
)
from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
from boltz.data.types import Manifest
from boltz.main import predict, process_inputs


def sha(path: Path) -> str:
    """Hash an explicit source prerequisite."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(sequence: str) -> tuple[int, ...]:
    """Encode aligned A3M columns with the installed native residue alphabet."""
    return tuple(
        const.token_ids[const.prot_letter_to_token[c]]
        for c in sequence
        if not c.islower()
    )


def audit_case(case: dict, cache: Path) -> dict:
    """Validate actual native full dataset pairing, four ligands and no guidance."""
    torch.set_num_threads(2)
    msas = Path(case["msas"])
    alpha = list(csv.DictReader((msas / "alpha_paired.csv").open()))
    beta = list(csv.DictReader((msas / "beta_paired.csv").open()))
    assert [r["key"] for r in alpha] == [r["key"] for r in beta]
    expected = {
        int(a["key"]): (encode(a["sequence"]), encode(b["sequence"]))
        for a, b in zip(alpha, beta, strict=True)
    }
    assert len(expected) == len(alpha) > 1
    original_get = PredictionDataset.__getitem__
    entered = False

    def guarded(dataset: object, index: int) -> dict:
        nonlocal entered
        if entered:
            raise RuntimeError("Native dataset fallback attempted")
        entered = True
        try:
            return original_get(dataset, index)
        finally:
            entered = False

    PredictionDataset.__getitem__ = guarded
    defaults = {p.name: p.default for p in predict.params}
    seed_everything(case["seed"])
    out = Path(case["remote_output"]) / ("boltz_results_" + case["case"])
    if out.exists():
        raise FileExistsError("CPU audit requires fresh native processed outputs")
    source = Path(case["remote_input"])
    document = json.loads(source.read_text())
    assert set(document) == {"version", "sequences"}
    process_inputs(
        data=[source],
        out_dir=out,
        ccd_path=cache / "ccd.pkl",
        mol_dir=cache / "mols",
        msa_server_url=defaults["msa_server_url"],
        msa_pairing_strategy=defaults["msa_pairing_strategy"],
        boltz2=True,
    )
    target = out / "processed"
    manifest = Manifest.load(target / "manifest.json")
    assert len(manifest.records) == 1
    record = manifest.records[0]
    assert record.templates == []
    if record.inference_options is not None:
        assert not record.inference_options.pocket_constraints
        assert not record.inference_options.contact_constraints
    data = load_input(
        record,
        target / "structures",
        target / "msa",
        target / "constraints",
        target / "templates",
        target / "mols",
    )
    assert set(str(c["name"]) for c in data.structure.chains) == set("ABCDEFGHIJKL")
    assert not data.templates
    tokens = Boltz2Tokenizer().tokenize(data)
    names = {int(c["asym_id"]): str(c["name"]) for c in data.structure.chains}
    chains = {
        name: [
            int(t["token_idx"])
            for t in tokens.tokens
            if names[int(t["asym_id"])] == name
        ]
        for name in "ABCDEFGH"
    }
    actual_keys = {}
    for idx, msa in data.msa.items():
        name = names[int(idx)]
        if name in "ABCDEFGH":
            actual_keys[name] = [int(s["taxonomy"]) for s in msa.sequences]
    assert set(actual_keys) == set("ABCDEFGH")
    assert all(keys == actual_keys["A"] for keys in actual_keys.values())
    assert len(actual_keys["A"]) > 1 and actual_keys["A"][0] == -1
    assert actual_keys["A"] == [int(r["key"]) for r in alpha[: len(actual_keys["A"])]]
    module = Boltz2InferenceDataModule(
        manifest=manifest,
        target_dir=target / "structures",
        msa_dir=target / "msa",
        mol_dir=cache / "mols",
        num_workers=0,
        constraints_dir=target / "constraints",
        template_dir=target / "templates",
        extra_mols_dir=target / "mols",
    )
    features = module.predict_dataloader().dataset[0]
    dummy = native.load_dummy_templates_features(1, len(tokens.tokens))
    assert all(torch.equal(features[k], value) for k, value in dummy.items())
    paired = features["msa_paired"]
    sequences = features["msa"]
    protein = [idx for chain in chains.values() for idx in chain]
    source_keys = []
    for row in range(len(paired)):
        if not bool((paired[row, protein] == 1).all()):
            continue
        a = tuple(int(v) for v in sequences[row, chains["A"]])
        b = tuple(int(v) for v in sequences[row, chains["E"]])
        assert all(
            tuple(int(v) for v in sequences[row, chains[c]]) == a for c in "ABCD"
        )
        assert all(
            tuple(int(v) for v in sequences[row, chains[c]]) == b for c in "EFGH"
        )
        source_key = actual_keys["A"][row]
        assert expected[source_key] == (a, b)
        source_keys.append(source_key)
    assert len(source_keys) > 1 and any(key >= 1 for key in source_keys)
    assert len(source_keys) == len(set(source_keys))
    assert source_keys == actual_keys["A"][: len(source_keys)]
    report = {
        "case": case["case"],
        "input_sha256": sha(source),
        "record_sha256": sha(target / "records" / (case["case"] + ".json")),
        "record_templates_empty": True,
        "actual_native_template_features_equal_stock_dummy": True,
        "native_template_masks_and_coordinates_zero": True,
        "dummy_restype_categorical_one_hot": "native categorical placeholder; not coordinates",
        "native_paired_row_source_keys": source_keys,
        "native_nonquery_coherent_paired_rows": len(source_keys) - 1,
        "parsed_chain_keys": actual_keys,
        "total_precursor_pairs_in_CSV": len(alpha) - 1,
        "actual_native_CSV_limit_applied": len(actual_keys["A"]) < len(alpha),
        "full_native_dataset_path_passed": True,
    }
    PredictionDataset.__getitem__ = original_get
    report.update(
        native_featurizer_source_sha256=sha(Path(native.__file__)),
        native_featurizer_source_path=str(Path(native.__file__)),
        audit_source_sha256=sha(Path(__file__)),
        paired_source_receipt_sha256=sha(msas / "precursor_pairing_receipt.json"),
        status="accepted",
        stock_native_entry_only=True,
    )
    return report
