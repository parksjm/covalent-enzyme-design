"""Assemble verified stripped host exports into a portable umbrella-window movie."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import parmed

OUT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(OUT / "source_snapshot/src"))

from nyl12_covalent.neb_us_movie_vmd import render_movie_vmd  # noqa: E402
from nyl12_covalent.path_inspection import dcd_counts  # noqa: E402


def digest(path: Path) -> str:
    """Return a streaming SHA-256 digest without loading a trajectory into memory."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def verify_distances(
    campaign: Path, rows: list[dict], kept_ids: list[int], roles: dict[str, int]
) -> dict:
    """Compare every movie frame with its recorded production-frame distances."""
    mapping = {old: new + 1 for new, old in enumerate(kept_ids)}
    pairs = [
        ("r_c_n", "carbon", "leaving_n"),
        ("r_c_og1", "carbon", "og1"),
        ("r_og1_hg1", "og1", "hg1"),
        ("r_nleave_hg1", "leaving_n", "hg1"),
    ]
    lines = [f"parm {OUT / 'path.parm7'}", f"trajin {OUT / 'path.dcd'}"]
    for name, first, second in pairs:
        lines.append(
            f"distance {name} @{mapping[roles[first]]} "
            f"@{mapping[roles[second]]} out {OUT / 'verified_distances.dat'} noimage"
        )
    script = OUT / "verify_distances.cpptraj"
    script.write_text("\n".join([*lines, "run", "quit"]) + "\n")
    with (OUT / "verify_distances.log").open("w") as stream:
        subprocess.run(
            ["/data/jpc/code/ambertools26/venv/bin/cpptraj", "-i", str(script)],
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )
    cached = {}
    expected = []
    for row in rows:
        key = (row["window"], row["block"])
        if key not in cached:
            folder = campaign / key[0] / key[1] / "observables"
            metadata = json.loads((folder / "observables.json").read_text())
            if digest(folder / "observables.csv") != metadata["csv_sha256"]:
                raise ValueError(f"Original observable CSV changed: {key}")
            with (folder / "observables.csv").open() as stream:
                cached[key] = list(csv.DictReader(stream))
        original = cached[key][row["source_frame1based"] - 1]
        expected.append([float(original[name]) for name, _, _ in pairs])
    actual = np.loadtxt(OUT / "verified_distances.dat")
    np.testing.assert_allclose(actual[:, 1:], expected, atol=0.001, rtol=0)
    return {
        "all_frames_matched": True,
        "distance_columns": [p[0] for p in pairs],
        "tolerance_A": 0.001,
        "maximum_absolute_difference_A": float(
            np.max(np.abs(actual[:, 1:] - expected))
        ),
    }


def main() -> None:
    """Verify host topology identities and frame plans, then assemble the movie."""
    request = json.loads((OUT / "request.json").read_text())
    campaign = Path(request["campaign"])
    config = json.loads((campaign / "prepared.json").read_text())
    full_topology = Path(config["topology"])
    full = parmed.load_file(str(full_topology))
    host_parts = {}
    by_window = {}
    reference_ids = None
    qm_ids = None
    source_hashes = {str(full_topology): digest(full_topology)}
    for host, expected_windows in request["host_windows"].items():
        folder = OUT / "work" / host
        record = json.loads((folder / "manifest.json").read_text())
        if set(record["windows"]) != set(expected_windows):
            raise ValueError(f"Host window assignment mismatch: {host}")
        if record["topology_sha256"] != source_hashes[str(full_topology)]:
            raise ValueError(f"Original topology mismatch: {host}")
        ids = record["original_atom_ids_1based"]
        if reference_ids is None:
            reference_ids, qm_ids = ids, record["qm_atom_ids_1based"]
        if ids != reference_ids or record["qm_atom_ids_1based"] != qm_ids:
            raise ValueError(f"Original atom mapping mismatch: {host}")
        actual_top = parmed.load_file(str(folder / "part.parm7"))
        if len(actual_top.atoms) != len(ids):
            raise ValueError(f"Stripped topology count mismatch: {host}")
        for new, old_id in enumerate(ids):
            a, b = actual_top.atoms[new], full.atoms[old_id - 1]
            if (a.name, a.atomic_number, a.residue.name) != (
                b.name,
                b.atomic_number,
                b.residue.name,
            ):
                raise ValueError(f"Atom identity mismatch: {host}/{new}")
            np.testing.assert_allclose(a.charge, b.charge, atol=1e-7, rtol=0)
        if dcd_counts(folder / "part.dcd") != (record["frames"], len(ids)):
            raise ValueError(f"Part DCD header mismatch: {host}")
        host_parts[host] = folder
        for row in record["frame_plan"]:
            by_window.setdefault(row["window"], []).append((host, row))
        source_hashes[str(folder / "manifest.json")] = digest(folder / "manifest.json")
        source_hashes[str(folder / "part.dcd")] = digest(folder / "part.dcd")
        source_hashes[str(folder / "part.parm7")] = digest(folder / "part.parm7")
    if set(by_window) != set(request["windows"]):
        raise ValueError("Not all requested windows have exported frames")
    shutil.copy2(next(iter(host_parts.values())) / "part.parm7", OUT / "path.parm7")
    lines = [f"parm {OUT / 'path.parm7'}"]
    rows = []
    for label in request["windows"]:
        entries = by_window[label]
        owners = {host for host, row in entries}
        local_indices = [row["output_frame0based"] for host, row in entries]
        if (
            len(owners) != 1
            or local_indices != list(range(local_indices[0], local_indices[-1] + 1))
            or len(entries) != 6 * request["frames_per_block"]
        ):
            raise ValueError(f"Noncontiguous or incomplete host frames: {label}")
        host = next(iter(owners))
        lines.append(
            f"trajin {host_parts[host] / 'part.dcd'} "
            f"{local_indices[0] + 1} {local_indices[-1] + 1}"
        )
        for _, row in entries:
            rows.append({**row, "host": host, "output_frame0based": len(rows)})
    lines += [f"trajout {OUT / 'path.dcd'} dcd", "run", "quit"]
    script = OUT / "assemble.cpptraj"
    script.write_text("\n".join(lines) + "\n")
    with (OUT / "assemble.log").open("w") as stream:
        subprocess.run(
            ["/data/jpc/code/ambertools26/venv/bin/cpptraj", "-i", str(script)],
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )
    frames, atoms = dcd_counts(OUT / "path.dcd")
    if (frames, atoms) != (request["frames_expected"], len(reference_ids)):
        raise ValueError("Assembled DCD count mismatch")
    geometry_check = verify_distances(
        campaign, rows, reference_ids, request["roles_1based"]
    )
    (OUT / "view_umbrella.vmd").write_text(
        render_movie_vmd(full, reference_ids, qm_ids, request["roles_1based"], frames)
    )
    (OUT / "frame_manifest.json").write_text(json.dumps(rows, indent=2) + "\n")
    with (OUT / "frame_manifest.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "frames": frames,
        "atoms": atoms,
        "windows": request["windows"],
        "frames_per_window": 6 * request["frames_per_block"],
        "interpretation": request["interpretation"],
        "original_atom_ids_1based": reference_ids,
        "qm_atom_ids_1based": qm_ids,
        "roles_1based": request["roles_1based"],
        "geometry_verification": geometry_check,
        "source_sha256": source_hashes,
        "checks": [
            "all63_windows",
            "378_accepted_blocks",
            "1890_frames",
            "original_atom_mapping",
            "source_topology_charge_identity",
            "disjoint_host_assignments",
            "ordered_window_concatenation",
        ],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    names = [
        "path.parm7",
        "path.dcd",
        "view_umbrella.vmd",
        "manifest.json",
        "frame_manifest.json",
        "frame_manifest.csv",
        "request.json",
        "README.md",
    ]
    (OUT / "MANIFEST.sha256").write_text(
        "".join(f"{digest(OUT / name)}  {name}\n" for name in names)
    )
    logging.info(
        "Published %s frames, %s atoms, %s windows", frames, atoms, len(by_window)
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
