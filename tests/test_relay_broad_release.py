"""Test fixed four-template allocation and fail-closed generation release."""

import copy

import pytest

from nyl12_covalent.relay_broad_pilot import symmetry_manifest
from nyl12_covalent.relay_broad_release import (
    four_template_schedule,
    validate_generation_request,
    validate_idle_gpu,
    validate_original_lineage,
)


def test_allocation_requires_two_per_kind_and_exactly24() -> None:
    schedule = four_template_schedule(["o1", "o3"], ["r1", "r3"])
    assert [entry["raw_samples"] for entry in schedule] == [6] * 4
    assert sum(entry["raw_samples"] for entry in schedule[:2]) == 12
    for original, remodeled in (
        (["o1"], ["r1", "r3"]),
        (["o1", "o3"], []),
        (["x", "o3"], ["x", "r3"]),
    ):
        with pytest.raises(ValueError):
            four_template_schedule(original, remodeled)


def test_exact_mask_ties_and_method_required() -> None:
    mask = {"alpha": list(range(1, 173)), "beta": list(range(1, 58))}
    protection = {"mutable": mask, "counts": {"fixed_total": 99}}
    generation = {
        "mask": mask,
        **symmetry_manifest(mask),
        "temperature": 0.1,
        "native_identity_bias": 0,
        "max_raw_samples": 24,
    }
    validate_generation_request({"generation": generation}, protection)
    for field, replacement in (
        ("mask", {"alpha": [], "beta": []}),
        ("symmetry_weights", "1"),
        ("temperature", 0.2),
        ("max_raw_samples", 25),
    ):
        altered = copy.deepcopy(generation)
        altered[field] = replacement
        with pytest.raises(ValueError):
            validate_generation_request({"generation": altered}, protection)
    with pytest.raises(ValueError):
        validate_generation_request(
            {"generation": generation}, {"mutable": mask, "counts": {"fixed_total": 98}}
        )


def test_only_two_original5A_qualified_sources_release() -> None:
    parents = [
        {"parent_source_id": source, "noise_A": 5, "sequence_scope": "native_WT"}
        for source in ("medoid01", "medoid03")
    ]
    validate_original_lineage(parents)
    with pytest.raises(ValueError):
        validate_original_lineage(parents[:1])
    parents[1]["parent_source_id"] = "medoid01"
    with pytest.raises(ValueError):
        validate_original_lineage(parents)
    parents[1]["parent_source_id"] = "medoid03"
    parents[0]["noise_A"] = 2
    with pytest.raises(ValueError):
        validate_original_lineage(parents)


def test_gpu_process_memory_and_activity_hold_generation() -> None:
    validate_idle_gpu("38, 0\n", "")
    for rows, pids in (("38,0", "1290"), ("257,0", ""), ("38,1", "")):
        with pytest.raises(RuntimeError):
            validate_idle_gpu(rows, pids)


def test_native_template_preparation_uses_explicit_four_by_six(
    tmp_path, monkeypatch
) -> None:
    import importlib.util
    from pathlib import Path

    from nyl12_covalent.relay_broad_release import FOUR_TEMPLATE_PROTOCOL

    spec = importlib.util.spec_from_file_location(
        "pilot_release_test", Path("scripts/run_relay_broad_pilot.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = tmp_path / "parent.pdb"
    source.write_text("END\n")
    monkeypatch.setattr(module, "parser_audit", lambda *args: {"accepted": True})
    monkeypatch.setattr(
        module,
        "selected_remodels",
        lambda *args: [
            {"parent_id": parent, "pdb_path": str(source)} for parent in ("r1", "r3")
        ],
    )
    manifest = {
        "allocation_protocol": FOUR_TEMPLATE_PROTOCOL,
        "request": {
            "original_parent_ids": ["o1", "o3"],
            "original_parent_pdbs": {"o1": str(source), "o3": str(source)},
        },
        "accepted_manifests": [],
        "native": {
            "generation_engine": "run.py",
            "ES_template": {"id": "ES", "pdb": str(source)},
        },
    }
    out = tmp_path / "out"
    out.mkdir()
    templates = module.prepare_templates(manifest, out)
    assert len(templates) == 5
    assert [template["raw_samples"] for template in templates[:-1]] == [6, 6, 6, 6]
    assert templates[-1]["state"] == "ES"
    manifest["allocation_protocol"] = "unreviewed"
    with pytest.raises(ValueError, match="Unknown explicit"):
        module.prepare_templates(manifest, out)
