"""Reuse the established water-retaining trajectory exporter for pilot snapshots."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
GRID = HERE.parent
ROOT = GRID.parents[3]
WATER_EXPORT = GRID.parent / "07_water_to_asp272_ti1_like"
sys.path.insert(0, str(WATER_EXPORT))
from export_water_path import export  # noqa: E402


def digest(path: Path) -> str:
    """Return SHA-256 of a retained artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_table(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a rectangular inspection table."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def snapshot_exporter(output: Path) -> None:
    """Retain immutable exporter source identities before the first refresh."""
    paths = [
        Path(__file__),
        HERE / "monitor_progress.py",
        WATER_EXPORT / "export_water_path.py",
        WATER_EXPORT / "water_path_vmd_support.py",
    ]
    hashes = {str(p): digest(p) for p in paths}
    folder = output / "exporter_source"
    if folder.exists():
        recorded = json.loads((folder / "manifest.json").read_text())
        if recorded["source_sha256"] != hashes:
            raise ValueError("Exporter sources changed; use a new progress directory")
        return
    folder.mkdir()
    for path in paths:
        shutil.copy2(path, folder / path.name)
    (folder / "manifest.json").write_text(
        json.dumps(
            {
                "source_sha256": hashes,
                "command": sys.argv,
                "frame_definition": "Only accepted initial/every25/final restart snapshots; no trial steps",
            },
            indent=2,
        )
        + "\n"
    )


def snapshot_frames(
    pilot: Path, rows: list[dict[str, Any]]
) -> list[tuple[Path, dict[str, Any]]]:
    """Find complete snapshots corresponding to authoritative accepted records."""
    frames = []
    seen = set()
    for row in rows:
        if row["label"] == "accepted" and row["accepted_iteration"] % 25:
            continue
        name = f"accepted_{row['accepted_iteration']:05d}_{row['label']}.ncrst"
        path = pilot / row["stage"] / name
        # A live writer appends the ledger before it finishes the NetCDF file.
        # The caller supplies only files stable over two monitor checkpoints.
        if not path.is_file():
            continue
        identity = digest(path)
        if identity not in seen:
            frames.append((path, row))
            seen.add(identity)
    return frames


def correct_bundle(folder: Path, frame_records: list[dict[str, Any]]) -> None:
    """Use element colors and dynamic QM H bonds while retaining QM licorice."""
    manifest_path = folder / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    mapping = {
        old + 1: new
        for new, old in enumerate(manifest["stripped_to_original_atom_indices_0based"])
    }
    protocol = json.loads((GRID / "protocol.json").read_text())
    h = [
        mapping[a["amber_atom_id"]]
        for a in protocol["qm_atoms_ordered"]
        if a["element"].upper() == "H"
    ]
    if len(h) != 54 or manifest["retained_water_residue_1based"] != 1576:
        raise ValueError("Expected retained WAT1576 and all54 real QM H")
    qm = " ".join(map(str, manifest["qm_indices_stripped_0based"]))
    hydrogens = " ".join(map(str, h))
    cleanup = (
        f"set proton_indices {{{hydrogens}}}\n"
        "set sel [atomselect $molid all]\nset bonds [$sel getbonds]\n"
        "for {set i 0} {$i < [llength $bonds]} {incr i} {\n"
        " set keep {}\n"
        " if {[lsearch -exact $proton_indices $i] < 0} {\n"
        "  foreach j [lindex $bonds $i] {\n"
        "   if {[lsearch -exact $proton_indices $j] < 0} {lappend keep $j}\n"
        "  }\n }\n lset bonds $i $keep\n}\n"
        "$sel setbonds $bonds\n$sel delete\n"
    )
    path = folder / "view_path.vmd"
    text = path.read_text().replace(
        "mol delrep 0 $molid\n", "mol delrep 0 $molid\n" + cleanup
    )
    for color in (1, 4, 7):
        text = text.replace(f"mol color ColorID {color}", "mol color Element")
    text = text.replace(
        "Alternative Asp272-COOH / neutral-Thr-N path",
        "Occupancy pilot accepted snapshots; not dynamics",
    )
    text += (
        f"\nmol color Element\nmol representation DynamicBonds 1.4 0.1 1\n"
        f'mol selection "index {qm}"\nmol addrep $molid\n'
        "color Element H white\ncolor Element N blue\ncolor Element O red\n"
    )
    for first, second in [
        (13905, 13917),
        (13916, 19920),
        (13916, 20318),
        (13905, 20318),
        (19920, 19933),
    ]:
        text += f"label add Bonds $molid/{mapping[first]} $molid/{mapping[second]}\n"
    seen_labels: set[str] = set()
    lines = []
    for line in text.splitlines():
        if line.startswith("label add "):
            if line in seen_labels:
                continue
            seen_labels.add(line)
        lines.append(line)
    text = "\n".join(lines) + "\nlabel show Bonds all\nlabel show Angles all\n"
    path.write_text(text)
    manifest.update(
        {
            "state": "Occupancy-bias diagnostic accepted snapshots; not MD/MFEP/PMF",
            "accepted_frame_records": frame_records,
            "proton_rendering": "All54 real QM H static bonds removed; dynamic cutoff1.4A",
            "pilot_exporter_sha256": digest(Path(__file__)),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    checksum = folder / "MANIFEST.sha256"
    names = [line.split("  ", 1)[1] for line in checksum.read_text().splitlines()]
    checksum.write_text("".join(f"{digest(folder / name)}  {name}\n" for name in names))


def export_movie(output: Path, frames: list[tuple[Path, dict[str, Any]]]) -> None:
    """Reuse the validated stripped exporter with distinct frame-parent names."""
    if not frames:
        return
    signature = "\n".join(f"{path}:{digest(path)}" for path, _ in frames)
    marker = output / "exported_frames.txt"
    if marker.exists() and marker.read_text() == signature:
        return
    staging = output / "frame_links"
    staging.mkdir(exist_ok=True)
    links, energies, distances = [], [], []
    baseline = frames[0][1]["physical_energy_kcal_mol"]
    for path, row in frames:
        window = f"{row['stage']}_{path.stem}"
        folder = staging / window
        folder.mkdir(exist_ok=True)
        link = folder / "accepted.ncrst"
        if not link.exists():
            link.symlink_to(path.resolve())
        elif link.resolve() != path.resolve():
            raise ValueError("Frame-link collision")
        links.append(link)
        energies.append(
            {
                "window": window,
                "energy_kcal_mol": row["physical_energy_kcal_mol"] - baseline,
                "energy_bias_inclusive_kcal_mol": row["total_energy_kcal_mol"]
                - baseline,
            }
        )
        distances.append(
            {
                "window": window,
                "q_CN": row["q"],
                "CN_N": row["CN_N"],
                "CN_OG1": row["CN_OG1"],
                "OG1_C_A": row["r_A"],
            }
        )
    energy, distance = (
        output / "snapshot_energies.csv",
        output / "snapshot_coordinates.csv",
    )
    csv_table(energy, energies)
    csv_table(distance, distances)
    bundle = output / "vmd"
    export(
        GRID.parent / "00_system/ash274_solvated.prmtop",
        links,
        energy,
        distance,
        bundle,
    )
    correct_bundle(bundle, [row for _, row in frames])
    marker.write_text(signature)
