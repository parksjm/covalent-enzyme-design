"""Prepare provenance-stable template-conditioned Boltz pilot inputs."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from pathlib import Path

LOGGER = logging.getLogger(__name__)

CAMPAIGN_NAME = "template_conditioned_es"
SEQUENCE_NAME = "nyl12_pa6_l4_to_l2_template"
QUERY_CHAINS = ("A", "E", "B", "F", "C", "G", "D", "H")
# Boltz uses mmCIF label-asym identifiers A--D for the four precursor chains;
# their corresponding deposited auth-asym identifiers are A, F, H, and V.
TEMPLATE_CHAINS = ("A", "A", "B", "B", "C", "C", "D", "D")


def sha256_file(path: Path) -> str:
    """Return the SHA256 digest of a file.

    Args:
        path: Existing file to hash.

    Returns:
        Lowercase hexadecimal SHA256 digest.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_source_yaml(text: str) -> None:
    """Fail if the validated source no longer has the expected pilot chemistry.

    Args:
        text: Source Boltz YAML text.
    """
    if text.count("- ligand:") != 1:
        raise ValueError("source YAML must contain exactly one ligand")
    if text.count("- contact:") != 2:
        raise ValueError("source YAML must contain exactly two contacts")
    if "pa6" not in text.lower() and "[NH3+]CCCCCC" not in text:
        raise ValueError("source YAML does not contain the expected PA6 ligand")
    if "templates:" in text:
        raise ValueError("source YAML already contains a templates block")


def template_block(template_cif: Path) -> str:
    """Return the explicit split-chain template mapping block.

    Args:
        template_cif: Nyl12 crystallographic mmCIF path.

    Returns:
        YAML text accepted by the installed custom Boltz2 schema.
    """
    if not template_cif.is_file():
        raise FileNotFoundError(template_cif)
    query = ", ".join(QUERY_CHAINS)
    template = ", ".join(TEMPLATE_CHAINS)
    return (
        "templates:\n"
        f"- cif: {template_cif.resolve()}\n"
        f"  chain_id: [{query}]\n"
        f"  template_id: [{template}]\n"
    )


def build_template_yaml(source_yaml: Path, template_cif: Path) -> str:
    """Add the explicit n12c template mapping to a validated Boltz YAML.

    Args:
        source_yaml: Previously validated PA6 L4-to-L2 Boltz input.
        template_cif: Nyl12 crystallographic mmCIF.

    Returns:
        Complete template-conditioned YAML text.
    """
    text = source_yaml.read_text(encoding="utf-8")
    validate_source_yaml(text)
    return text.rstrip() + "\n" + template_block(template_cif)


def write_seed_manifest(path: Path, yaml_path: Path, seed_base: int) -> None:
    """Write a deterministic 12-seed pilot manifest.

    Args:
        path: Destination CSV path.
        yaml_path: Generated template-conditioned YAML.
        seed_base: First nonnegative random seed.
    """
    if seed_base < 0:
        raise ValueError("seed_base must be nonnegative")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "campaign",
                "sequence_name",
                "model_index",
                "seed",
                "yaml_path",
            ),
        )
        writer.writeheader()
        for model_index in range(12):
            writer.writerow(
                {
                    "campaign": CAMPAIGN_NAME,
                    "sequence_name": SEQUENCE_NAME,
                    "model_index": model_index,
                    "seed": seed_base + model_index,
                    "yaml_path": str(yaml_path.resolve()),
                }
            )


def prepare_pilot(
    source_yaml: Path,
    template_cif: Path,
    output_dir: Path,
    seed_base: int,
) -> dict[str, object]:
    """Write the pilot YAML, seed manifest, and provenance record.

    Args:
        source_yaml: Validated non-template Boltz input.
        template_cif: Nyl12 crystallographic mmCIF.
        output_dir: Destination directory for compact pilot inputs.
        seed_base: First of 12 deterministic seeds.

    Returns:
        JSON-serializable provenance record.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = output_dir / "nyl12_pa6_l4_to_l2_template.yaml"
    yaml_path.write_text(
        build_template_yaml(source_yaml, template_cif), encoding="utf-8"
    )
    manifest_path = output_dir / "seeds.csv"
    write_seed_manifest(manifest_path, yaml_path, seed_base)
    provenance: dict[str, object] = {
        "source_yaml": str(source_yaml.resolve()),
        "source_yaml_sha256": sha256_file(source_yaml),
        "template_cif": str(template_cif.resolve()),
        "template_cif_sha256": sha256_file(template_cif),
        "generated_yaml": str(yaml_path.resolve()),
        "generated_yaml_sha256": sha256_file(yaml_path),
        "seed_manifest": str(manifest_path.resolve()),
        "seed_manifest_sha256": sha256_file(manifest_path),
        "seed_count": 12,
        "seed_base": seed_base,
        "query_chains": list(QUERY_CHAINS),
        "template_chains": list(TEMPLATE_CHAINS),
        "calculation_launched": False,
    }
    provenance_path = output_dir / "provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info("Prepared template-conditioned pilot inputs in %s", output_dir)
    return provenance
