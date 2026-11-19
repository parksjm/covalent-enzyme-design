#!/usr/bin/env python3
"""Publish a rendering-only relay update without rebuilding trajectory coordinates."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.export_relay_campaign_progress import (  # noqa: E402
    DEFAULT,
    atomic_json,
    digest,
    journal,
)
from nyl12_covalent.relay_progress import verify_mapping  # noqa: E402
from nyl12_covalent.relay_progress_vmd import (  # noqa: E402
    VMD_API_REFERENCES,
    substrate_indices,
    vmd_script,
)

LOGGER = logging.getLogger(__name__)


def verify(folder: Path) -> list[str]:
    """Verify every listed portable payload hash and return its relative names."""
    names = []
    for line in (folder / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if digest(folder / name) != expected:
            raise ValueError(f"Payload checksum mismatch: {folder / name}")
        names.append(name)
    return names


def refresh_view(campaign: Path, base: Path) -> Path:
    """Create a new immutable view update from verified unchanged coordinates.

    Args:
        campaign: Existing campaign with case.json and progress directory.
        base: Immutable completed update with root and VMD checksum manifests.

    Returns:
        Newly published immutable update directory.

    Raises:
        ValueError: If payload checksums or scientific atom mappings disagree.
    """
    verify(base)
    bundle_names = verify(base / "vmd")
    case_path = campaign / "case.json"
    case = json.loads(case_path.read_text())
    topology = verify_mapping(case)
    progress = campaign / "progress"
    numbers = [
        int(p.name.rsplit("_", 1)[1])
        for p in progress.iterdir()
        if p.name.startswith(("update_", ".building_"))
    ]
    number = max(numbers, default=0) + 1
    name = f"update_{number:06d}"
    pending = progress / f".building_{number:06d}"
    pending.mkdir()
    for path in base.iterdir():
        if path.is_file() and path.name not in {"MANIFEST.sha256", "progress.json"}:
            os.link(path, pending / path.name)
        elif path.is_dir() and path.name not in {"vmd", "source", "frame_links"}:
            raise ValueError(f"Unexpected base directory: {path}")
    bundle = pending / "vmd"
    bundle.mkdir()
    for path in (base / "vmd").iterdir():
        if not path.is_file():
            raise ValueError(f"Unexpected VMD payload: {path}")
        os.link(path, bundle / path.name)
    manifest = json.loads((bundle / "manifest.json").read_text())
    source_paths = [
        Path(path)
        for path in json.loads((base / "progress.json").read_text())["identity"][
            "source_sha256"
        ]
    ]
    source_paths.append(Path(__file__).resolve())
    source_paths = list(dict.fromkeys(source_paths))
    hashes = {str(p): digest(p) for p in source_paths}
    source = pending / "source"
    source.mkdir()
    for path in source_paths:
        shutil.copy2(path, source / path.name)
    atomic_json(
        source / "manifest.json",
        {"source_sha256": hashes, "base_update": str(base), "command": sys.argv},
    )
    substrate = substrate_indices(
        topology, case["atom_map"]["carbon"] - 1, case["atom_map"]["leaving_n"] - 1
    )
    manifest["view_correction"] = {
        "vmd_api_references": VMD_API_REFERENCES,
        "source_sha256": hashes,
        "base_update": str(base),
        "substrate_original_atom_ids_1based": sorted(i + 1 for i in substrate),
        "default_distance_atom_ids_1based": [
            [case["atom_map"][a], case["atom_map"][b]]
            for a, b in (
                ("thr_og", "carbon"),
                ("carbon", "leaving_n"),
                ("thr_og", "water_h1"),
                ("leaving_n", "water_h1"),
                ("oxygen", "phe_n"),
                ("oxygen", "asn_nd2"),
            )
        ],
        "camera": "Bounding box of rendered QM plus connected substrate only",
        "validation_limit": "VMD unavailable; topology, coordinates and headless Tcl checked",
    }
    # Replace paths rather than truncate files hardlinked to the old update.
    new_view = bundle / "view_path.vmd.pending"
    new_view.write_text(vmd_script(manifest, case, topology))
    new_view.replace(bundle / "view_path.vmd")
    atomic_json(bundle / "manifest.json", manifest)
    checksum = bundle / "MANIFEST.sha256.pending"
    checksum.write_text("".join(f"{digest(bundle / n)}  {n}\n" for n in bundle_names))
    checksum.replace(bundle / "MANIFEST.sha256")
    report = json.loads((base / "progress.json").read_text())
    report["identity"]["source_sha256"] = hashes
    report["remote_evidence_sha256"] = {
        "accepted_iterations.json": digest(pending / "accepted_iterations.json")
    }
    report["latest_update"] = name
    report["view_correction"] = {
        "base_update": str(base),
        "base_manifest_sha256": digest(base / "MANIFEST.sha256"),
        "coordinate_payload_reused": True,
    }
    atomic_json(pending / "progress.json", report)
    payload = sorted(
        p
        for p in pending.rglob("*")
        if p.is_file() and p != pending / "accepted_iterations.json"
    )
    (pending / "MANIFEST.sha256").write_text(
        "".join(f"{digest(p)}  {p.relative_to(pending)}\n" for p in payload)
    )
    verify(bundle)
    verify(pending)
    destination = progress / name
    pending.rename(destination)
    pointer = progress / f".latest_{os.getpid()}"
    pointer.symlink_to(name, target_is_directory=True)
    pointer.replace(progress / "latest")
    atomic_json(progress / "progress.json", report)
    return destination


def main() -> None:
    """Verify the base and publish one view-only update, journalling failures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=DEFAULT)
    parser.add_argument("--base", type=Path)
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    base = (
        args.base.resolve() if args.base else (campaign / "progress/latest").resolve()
    )
    journal("intent", f"Publish view-only relay correction from verified {base}")
    try:
        result = refresh_view(campaign, base)
    except Exception as error:
        journal(
            "outcome", f"View-only correction failed: {type(error).__name__}: {error}"
        )
        raise
    journal(
        "outcome", f"Published view-only correction {result}; coordinates unchanged"
    )
    LOGGER.info("Published %s", result)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
