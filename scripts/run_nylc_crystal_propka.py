"""Screen NylC Asp306/Asp308 pKas from deposited mature crystal models."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import re
import subprocess

LOGGER = logging.getLogger(__name__)
SUMMARY = re.compile(r"^\s*ASP\s+(306|308)\s+([A-Z])\s+(-?\d+\.\d+)\s+")
TARGET_ATOMS = {"CG", "OD1", "OD2"}


def sha256(path: Path) -> str:
    """Hash one immutable deposited or generated file.

    Args:
        path: Existing file.

    Returns:
        SHA-256 hex digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_protein(source: Path, target: Path, mature: bool) -> list[str]:
    """Write protein-only PDB with the experimental cleavage state explicit.

    Args:
        source: Deposited PDB file.
        target: New protein-only PDB input.
        mature: Whether the Asn266/Thr267 peptide bond has been cleaved.

    Returns:
        Chain IDs with complete Asp306/Asp308 carboxylates.
    """
    inventory: dict[tuple[str, int], set[str]] = {}
    source_chains = sorted(
        {
            line[21]
            for line in source.read_text().splitlines()
            if line.startswith("ATOM")
        }
    )
    if len(source_chains) > 13 or " " in source_chains:
        raise ValueError("Require 1–13 identified chains for explicit beta-chain IDs")
    beta_map = {
        chain: chr(ord("A") + len(source_chains) + index)
        for index, chain in enumerate(source_chains)
    }
    beta_chains: set[str] = set()
    previous_segment: tuple[str, bool] | None = None
    with source.open() as input_handle, target.open("w") as output_handle:
        for line in input_handle:
            if line.startswith("MODEL"):
                raise ValueError(
                    f"Multi-model PDB requires explicit model selection: {source}"
                )
            if not line.startswith("ATOM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            chain = line[21]
            if chain == " ":
                raise ValueError("Deposited protein chain lacks ID")
            resid = int(line[22:26])
            beta = mature and resid >= 267
            segment = (chain, beta)
            if previous_segment is not None and segment != previous_segment:
                output_handle.write("TER   \n")
            previous_segment = segment
            output_chain = beta_map[chain] if beta else chain
            if beta and resid == 267 and line[12:16].strip() == "N":
                beta_chains.add(output_chain)
            if resid in (306, 308):
                if line[17:20] != "ASP":
                    raise ValueError(f"Expected ASP {chain}:{resid}")
                inventory.setdefault((output_chain, resid), set()).add(
                    line[12:16].strip()
                )
            output_handle.write(line[:21] + output_chain + line[22:])
        output_handle.write("TER   \n")
        output_handle.write("END\n")
    if mature and beta_chains != set(beta_map.values()):
        raise ValueError("Mature Thr267 N terminus absent in at least one chain")
    target_chains = beta_chains if mature else set(source_chains)
    for chain in target_chains:
        for resid in (306, 308):
            if not TARGET_ATOMS <= inventory.get((chain, resid), set()):
                raise ValueError(f"Incomplete ASP {chain}:{resid} in {source}")
    return sorted(target_chains)


def extract_pkas(path: Path, chains: list[str]) -> list[dict[str, object]]:
    """Extract every target residue from PROPKA's final summary table.

    Args:
        path: Raw PROPKA ``.pka`` output.
        chains: Deposited protein chain IDs.

    Returns:
        Chain and residue pKa records.
    """
    text = path.read_text()
    if text.count("SUMMARY OF THIS PREDICTION") != 1:
        raise ValueError(f"Missing PROPKA summary in {path}")
    values: dict[tuple[str, int], float] = {}
    for line in text.split("SUMMARY OF THIS PREDICTION", 1)[1].splitlines():
        match = SUMMARY.match(line)
        if match is None:
            continue
        key = match.group(2), int(match.group(1))
        if key in values:
            raise ValueError(f"Duplicated pKa group {key}")
        values[key] = float(match.group(3))
    expected = {(chain, resid) for chain in chains for resid in (306, 308)}
    if set(values) != expected:
        raise ValueError(f"PROPKA target coverage differs: {set(values) ^ expected}")
    return [
        {"chain": chain, "residue": resid, "pka": values[chain, resid]}
        for chain, resid in sorted(expected)
    ]


def main() -> None:
    """Prepare deposited protein inputs and run coupled-site PROPKA."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", action="append", type=Path, required=True)
    parser.add_argument("--propka", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.output.exists():
        raise FileExistsError(args.output)
    if len({p.stem for p in args.structure}) != len(args.structure):
        raise ValueError("Require distinct deposited mature NylC PDBs")
    version = subprocess.run(
        [str(args.propka), "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    args.output.mkdir(parents=True)
    all_rows = []
    manifest = []
    for source in args.structure:
        source = source.resolve(strict=True)
        pdb_id = source.name.split("_", 1)[0].upper()
        case_dir = args.output / pdb_id
        case_dir.mkdir()
        protein = case_dir / "protein.pdb"
        if pdb_id not in {"5XYG", "7YU0", "7YU1"}:
            raise ValueError(f"Cleavage status not declared for {pdb_id}")
        mature = pdb_id == "5XYG"
        chains = write_protein(source, protein, mature=mature)
        commands = []
        outputs = {}
        for state, options, filename in (
            ("default", [], "protein.pka"),
            (
                "coupled_alternative",
                ["--display-coupled-residues"],
                "protein_alt_state.pka",
            ),
        ):
            command = [str(args.propka), *options, protein.name]
            run = subprocess.run(command, cwd=case_dir, capture_output=True, text=True)
            (case_dir / f"propka_{state}.log").write_text(run.stdout + run.stderr)
            if run.returncode:
                raise RuntimeError(
                    f"PROPKA {state} failed for {pdb_id}: {run.returncode}"
                )
            raw = case_dir / filename
            if not raw.is_file():
                raise FileNotFoundError(raw)
            if mature:
                summary = raw.read_text().split("SUMMARY OF THIS PREDICTION", 1)[1]
                for chain in chains:
                    if not re.search(
                        rf"(?m)^\s*N\+\s+267\s+{chain}\s+-?\d+\.\d+", summary
                    ):
                        raise ValueError(
                            f"Mature N terminus Thr267 not assigned in {raw}: {chain}"
                        )
            rows = [
                {"structure": pdb_id, "state": state, **row}
                for row in extract_pkas(raw, chains)
            ]
            all_rows.extend(rows)
            outputs[state] = sha256(raw)
            commands.append(command)
        manifest.append(
            {
                "structure": pdb_id,
                "source": str(source),
                "source_sha256": sha256(source),
                "protein_input_sha256": sha256(protein),
                "raw_output_sha256": outputs,
                "chains": chains,
                "commands": commands,
                "propka_version": version,
                "protein_only": True,
                "mature": mature,
                "beta_chain_split": mature,
            }
        )
        LOGGER.info("%s: %d Asp pKa values across two states", pdb_id, 4 * len(chains))
    with (args.output / "asp306_asp308_pka.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("structure", "state", "chain", "residue", "pka")
        )
        writer.writeheader()
        writer.writerows(all_rows)
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
