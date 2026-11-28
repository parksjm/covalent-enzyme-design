"""Focused checks for the corrected pilot release gate and provenance."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from nyl12_covalent.relay_ensemble_sequences import pool_samples  # noqa: E402
from run_corrected_pilot import (  # noqa: E402
    namespace_ids,
    verify_release,
    write_effective_manifest,
)


def test_namespace_ids_preserves_real_pool_lineage_references() -> None:
    samples = [
        {"sample_id": "parent01_1", "sequence": "AAA"},
        {"sample_id": "parent03_2", "sequence": "AAA"},
    ]
    candidates = pool_samples(samples)

    namespace_ids(samples, candidates)

    assert [sample["sample_id"] for sample in samples] == [
        "corr01_parent01_1",
        "corr01_parent03_2",
    ]
    assert candidates[0]["candidate_id"] == "relay_corr01_01"
    assert candidates[0]["sample_id"] == "corr01_parent01_1"
    assert candidates[0]["lineage"] == [
        "corr01_parent01_1",
        "corr01_parent03_2",
    ]


def test_release_requires_exact_manifest_and_wt_receipt(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}\n", encoding="utf-8")
    receipt = tmp_path / "wt.json"
    receipt.write_text('{"accepted": true}\n', encoding="utf-8")
    release = tmp_path / "release.json"
    release.write_text(
        json.dumps(
            {
                "approved": True,
                "accepted_WT_source_id": "WT_paired_seed_20261126",
                "force_converged": True,
                "positive_phi_A_B_sites": 4,
                "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                "allow_shared_host": True,
                "max_cpu_threads": 2,
                "gpu_uuid": "GPU-unit-test",
                "single_control_acceptance_path": str(receipt),
                "single_control_acceptance_sha256": hashlib.sha256(
                    receipt.read_bytes()
                ).hexdigest(),
            }
        ),
        encoding="utf-8",
    )

    assert verify_release(manifest, release)["gpu_uuid"] == "GPU-unit-test"


def test_effective_manifest_is_written_with_resolved_runtime_provenance(
    tmp_path: Path,
) -> None:
    wt_release = tmp_path / "wt_release.json"
    wt_release.write_text("{}\n", encoding="utf-8")
    resource_release = tmp_path / "resource_release.json"
    resource_release.write_text("{}\n", encoding="utf-8")
    manifest = {"native": {"environment": {"CUDA_VISIBLE_DEVICES": "GPU-unit-test"}}}
    output = tmp_path / "effective_manifest.json"

    receipt = write_effective_manifest(
        manifest,
        output,
        "source-manifest-sha",
        wt_release,
        resource_release,
        "GPU-unit-test",
        2,
    )

    assert json.loads(output.read_text(encoding="utf-8")) == manifest
    assert receipt["effective_cuda_visible_devices"] == "GPU-unit-test"
    assert (
        receipt["effective_manifest_sha256"]
        == hashlib.sha256(output.read_bytes()).hexdigest()
    )
    assert output.with_name("effective_manifest_acceptance.json").exists()
