"""Audit retained TS source, proton inventory and mature-local canonical exports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np

from nyl12_covalent.relay_broad_runtime import digest, export_local, write_new


def audit_sources(config: dict, destination: Path) -> dict:
    """Check source hashes, atom maps and beta renumbering without altering sources.

    Args:
        config: Explicit original-source audit records and template manifest.
        destination: New immutable audit receipt.

    Returns:
        Source mapping and retained proton inventory evidence.
    """
    template_path = Path(config["template_manifest"])
    template = json.loads(template_path.read_text())
    records = []
    for original in config["original_sources"]:
        entry = next(t for t in template["templates"] if t["name"] == original["id"])
        for field in ("heavy_pdb", "protonated_pdb", "restart"):
            item = entry[field]
            if digest(Path(item["path"])) != item["sha256"]:
                raise ValueError(f"Original source changed: {item['path']}")
        if (
            entry["heavy_pdb"]["residues"] != 1316
            or entry["heavy_pdb"]["ter_records"] != 9
        ):
            raise ValueError("TS parent must have 1316 residues and nine chain termini")
        if entry["native_proton_inventory"]["tracked_protons_1based"] != [
            13906,
            13907,
            13917,
            19937,
        ]:
            raise ValueError("TS source proton inventory changed")
        source = Path(entry["heavy_pdb"]["path"])
        mapping_path = Path(original["heavy_atom_map"])
        mapping = json.loads(mapping_path.read_text())
        numbering_path = Path(original["numbering_map"])
        numbering = json.loads(numbering_path.read_text())
        expected_numbering = {
            (chain, position, position - 232)
            for chain in "EFGH"
            for position in range(233, 329)
        }
        if {
            (a["chain"], a["source_residue"], a["input_residue"]) for a in numbering
        } != expected_numbering:
            raise ValueError("Source precursor-to-mature beta numbering map differs")
        reference_path = Path(original["reference_graph"])
        reference = json.loads(reference_path.read_text())
        if mapping != reference["atom_map"]:
            raise ValueError("Structural-to-Amber atom identity map changed")
        lines = [
            line
            for line in source.read_text().splitlines()
            if line[:6] in ("ATOM  ", "HETATM")
        ]
        if len(lines) != len(mapping) or len(
            {a["original_atom_id1based"] for a in mapping}
        ) != len(mapping):
            raise ValueError("Original structural-to-Amber map is not bijective")
        for line, atom in zip(lines, mapping, strict=True):
            if (line[21], int(line[22:26]), line[12:16].strip()) != (
                atom["chain"],
                atom["canonical_residue"],
                atom["atom_name"],
            ):
                raise ValueError(
                    "Original structural identities disagree with canonical atom map"
                )
        with tempfile.TemporaryDirectory(prefix="relay_ts_source_audit_") as directory:
            local = Path(directory) / "local.pdb"
            export_local(source, local)
            local_lines = [
                line
                for line in local.read_text().splitlines()
                if line[:6] in ("ATOM  ", "HETATM")
            ]
            keys = [
                (line[21], int(line[22:26]), line[12:16].strip())
                for line in local_lines
            ]
            expected_lines = [
                line
                for line in Path(reference["pdb"]).read_text().splitlines()
                if line[:6] in ("ATOM  ", "HETATM")
            ]
            expected_keys = [
                (line[21], int(line[22:26]), line[12:16].strip())
                for line in expected_lines
            ]
            if keys != expected_keys:
                raise ValueError(
                    "Mature-local export differs from native reference identities"
                )
            coords = np.array(
                [[float(line[i : i + 8]) for i in (30, 38, 46)] for line in local_lines]
            )
            if not np.array_equal(coords, np.asarray(reference["coords"])):
                raise ValueError("Canonical source export changed TS coordinates")
        records.append(
            {
                "id": original["id"],
                "accepted": True,
                "source": str(source),
                "source_sha256": digest(source),
                "reference_graph": str(reference_path),
                "reference_graph_sha256": digest(reference_path),
                "heavy_atom_map_sha256": digest(mapping_path),
                "numbering_map_sha256": digest(numbering_path),
                "atom_count": len(mapping),
                "residue_count": 1316,
                "structural_to_Amber_to_mature_local_bijective": True,
                "coordinates_exact": True,
                "native_proton_inventory": entry["native_proton_inventory"],
                "proton_owners": entry["proton_owners"],
                "proton_policy": template["proton_policy"],
            }
        )
    if len(records) != 2:
        raise ValueError("Require both original source audits")
    receipt = {
        "accepted": True,
        "template_manifest": str(template_path),
        "template_manifest_sha256": digest(template_path),
        "sources": records,
        "four_ligand_1328_residue_control_conversion_applied": False,
    }
    write_new(destination, receipt)
    return receipt


def main() -> None:
    """Write independent canonical-source audit evidence for the exact queue."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    audit_sources(json.loads(args.config.read_text()), args.destination)


if __name__ == "__main__":
    main()
