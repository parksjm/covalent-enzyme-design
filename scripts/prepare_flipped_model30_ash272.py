"""Stage four-site ASH272 or ASH274 preparation without running Amber."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from math import dist
from pathlib import Path
import shutil

ASP272_SITES = (968, 1064, 1160, 1256)


def sha256(path: Path) -> str:
    """Return the SHA256 digest of a file."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def prepare(
    source: Path, destination: Path, oxygen: str, target_residue: int = 272
) -> None:
    """Stage source parameters and ASH inputs with explicit oxygen identity.

    Args:
        source: Immutable model30 00_system directory.
        destination: New preparation directory, which must not already exist.
        oxygen: Original OD1 or OD2 atom to protonate at all four sites.
        target_residue: Catalytic Asp residue to protonate (272 or 274).

    Raises:
        ValueError: Required residue mapping or source inputs are invalid.
        FileExistsError: The destination already exists.
    """
    if oxygen not in {"OD1", "OD2"}:
        raise ValueError("An explicit OD1 or OD2 choice is required")
    if target_residue not in {272, 274}:
        raise ValueError("Target residue must be Asp272 or Asp274")
    sites = tuple(site + (2 if target_residue == 274 else 0) for site in ASP272_SITES)
    source = source.resolve(strict=True)
    if destination.exists():
        raise FileExistsError(destination)
    pdb = source / "input_complex_leap.pdb"
    leap = source / "leap_final.in"
    lines = pdb.read_text().splitlines(keepends=True)
    required = set(ASP272_SITES) | {site + 2 for site in ASP272_SITES}
    atoms: dict[int, dict[str, tuple[float, float, float]]] = {}
    indices: dict[int, list[int]] = {}
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = int(line[22:26])
        if residue not in required:
            continue
        if line[17:20] != "ASP" or line[16:17].strip():
            raise ValueError(f"Expected unambiguous ASP at residue {residue}")
        name = line[12:16].strip()
        positions = atoms.setdefault(residue, {})
        if name in positions:
            raise ValueError(f"Duplicate atom {residue}:{name}")
        positions[name] = tuple(
            float(line[start : start + 8]) for start in (30, 38, 46)
        )
        indices.setdefault(residue, []).append(index)
    evidence: list[dict[str, object]] = []
    mapping = {"OD1": "OD2", "OD2": "OD1"} if oxygen == "OD1" else {}
    for site in sites:
        partner = site + (2 if target_residue == 272 else -2)
        for residue in (site, partner):
            if not {"OD1", "OD2"} <= atoms.get(residue, {}).keys():
                raise ValueError(f"Incomplete ASP oxygens at {residue}")
        distances = {
            f"{first}-{second}": dist(atoms[site][first], atoms[partner][second])
            for first in ("OD1", "OD2")
            for second in ("OD1", "OD2")
        }
        for index in indices[site]:
            line = lines[index]
            name = line[12:16].strip()
            if name.startswith("H"):
                raise ValueError(f"Expected heavy atoms only at ASH site {site}")
            lines[index] = (
                line[:12]
                + f"{mapping.get(name, name):>4}"
                + line[16:17]
                + "ASH"
                + line[20:]
            )
        evidence.append(
            {
                "asp272_amber": site if target_residue == 272 else partner,
                "asp274_amber": partner if target_residue == 272 else site,
                "protonated_residue": target_residue,
                "original_proton_oxygen": oxygen,
                "amber_proton_oxygen": "OD2",
                "oxygen_distances_angstrom": distances,
                "original_to_amber_names": mapping or {"OD1": "OD1", "OD2": "OD2"},
            }
        )
    prefix, marker, _ = leap.read_text().partition("solvateOct")
    if (
        not marker
        or "saveamberparm" not in prefix.lower()
        or "addions" in prefix.lower()
    ):
        raise ValueError("Expected dry topology save before solvation and ions")
    release = sorted(
        {
            path
            for pattern in ("*.lib", "*.off", "*.mol2", "*.frcmod", "leaprc*")
            for path in source.glob(pattern)
            if path.is_file()
        }
    )
    if not release:
        raise ValueError("Source force-field release absent")
    source_hashes = {path.name: sha256(path) for path in [pdb, leap, *release]}
    destination.mkdir(parents=True)
    for path in release:
        shutil.copy2(path, destination / path.name)
    shutil.copy2(leap, destination / "leap_final.source.in")
    (destination / pdb.name).write_text("".join(lines))
    (destination / "leap_dry.in").write_text(prefix + "charge complex\nquit\n")
    (destination / "leap_count.in").write_text(
        prefix + "solvateOct complex TIP3PBOX 12.0\n"
        "saveAmberParm complex count.prmtop count.inpcrd\nquit\n"
    )
    manifest = {
        "source_directory": str(source),
        "source_sha256": source_hashes,
        "prepared_sha256": {
            path.name: sha256(path) for path in destination.iterdir() if path.is_file()
        },
        "sites": evidence,
        "protonated_residue": target_residue,
        "expected_protein_charge_change_e": 4,
        "status": "staged only; dry charge and 150mM ions require recount",
        "acceptance": "not assessed; no calculations executed",
    }
    (destination / "preparation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    logging.getLogger(__name__).info(
        "Staged ASH%d inputs in %s", target_residue, destination
    )


def main() -> None:
    """Expose the staging-only command line interface."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--proton-oxygen", choices=("OD1", "OD2"), required=True)
    parser.add_argument("--target-residue", type=int, choices=(272, 274), default=272)
    args = parser.parse_args()
    prepare(args.source, args.destination, args.proton_oxygen, args.target_residue)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
