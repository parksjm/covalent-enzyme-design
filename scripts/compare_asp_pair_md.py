"""Compare matched model30 Asp-pair microstates using one Tier A definition."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import logging
import math
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import numpy as np
import parmed

from nyl12_covalent.active_site_geometry import projected_herschlag_phi
from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

LOGGER = logging.getLogger(__name__)
SITES = ((968, 970), (1064, 1066), (1160, 1162), (1256, 1258))
STATES = {
    "ASP272_ASP274": ("ASP", "ASP"),
    "ASH272_ASP274": ("ASH", "ASP"),
    "ASP272_ASH274": ("ASP", "ASH"),
}


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atom(topology: Any, residue: int, name: str) -> int:
    atoms = [
        atom.idx for atom in topology.residues[residue - 1].atoms if atom.name == name
    ]
    if len(atoms) != 1:
        raise ValueError(f"Expected one {residue}:{name}, found {len(atoms)}")
    return int(atoms[0])


def _angle(first: np.ndarray, vertex: np.ndarray, last: np.ndarray) -> float:
    left, right = first - vertex, last - vertex
    cosine = float(np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right)))
    if not math.isfinite(cosine):
        raise ValueError("Degenerate Bürgi–Dunitz angle")
    return math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))


def _indices(topology: Any, state: str) -> dict[str, Any]:
    if state not in STATES:
        raise ValueError(f"Undeclared Asp-pair state: {state}")
    site_indices = []
    for asp272, asp274 in SITES:
        actual = (
            topology.residues[asp272 - 1].name,
            topology.residues[asp274 - 1].name,
        )
        if actual != STATES[state]:
            raise ValueError(f"Microstate mismatch at {asp272}/{asp274}: {actual}")
        pair = {
            "asp272_o": [_atom(topology, asp272, name) for name in ("OD1", "OD2")],
            "asp274_o": [_atom(topology, asp274, name) for name in ("OD1", "OD2")],
        }
        if state != "ASP272_ASP274":
            protonated = asp272 if state == "ASH272_ASP274" else asp274
            pair["ash_hd2"] = _atom(topology, protonated, "HD2")
            pair["ash_od2"] = _atom(topology, protonated, "OD2")
        site_indices.append(pair)
    thr_names = {atom.name for atom in topology.residues[928].atoms}
    if not {"N", "H1", "H2", "OG1", "HG1"} <= thr_names or "H3" in thr_names:
        raise ValueError("Catalytic Thr N terminus is not the declared neutral state")
    return {
        "og1": _atom(topology, 929, "OG1"),
        "carbon": _atom(topology, 1314, "C"),
        "oxygen": _atom(topology, 1314, "O"),
        "amide_n": _atom(topology, 1315, "N"),
        "adjacent_c": _atom(topology, 1314, "C5"),
        "phe_n": _atom(topology, 113, "N"),
        "asn_nd2": _atom(topology, 186, "ND2"),
        "sites": site_indices,
    }


def _row(x: np.ndarray, indices: dict[str, Any]) -> dict[str, float | str]:
    def point(name: str) -> np.ndarray:
        return x[indices[name]]

    og1, carbon, oxygen = point("og1"), point("carbon"), point("oxygen")
    metrics = {
        "d1": float(np.linalg.norm(og1 - carbon)),
        "d2": float(np.linalg.norm(point("phe_n") - oxygen)),
        "d3": float(np.linalg.norm(point("asn_nd2") - oxygen)),
        "attack": _angle(og1, carbon, oxygen),
        "projected_phi": projected_herschlag_phi(
            tuple(og1),
            tuple(oxygen),
            tuple(carbon),
            tuple(point("amide_n")),
            tuple(point("adjacent_c")),
        ),
    }
    row: dict[str, float | str] = {
        **metrics,
        "tier": classify_geometry(metrics) or "none",
    }
    for ordinal, site in enumerate(indices["sites"], start=1):
        left, right = site["asp272_o"], site["asp274_o"]
        row[f"site{ordinal}_asp_oo_min_A"] = min(
            float(np.linalg.norm(x[a] - x[b])) for a in left for b in right
        )
        if "ash_hd2" in site:
            proton = x[site["ash_hd2"]]
            row[f"site{ordinal}_ash_oh_A"] = float(
                np.linalg.norm(proton - x[site["ash_od2"]])
            )
            partner = right if site["ash_od2"] in left else left
            row[f"site{ordinal}_ash_h_partner_o_min_A"] = min(
                float(np.linalg.norm(proton - x[atom])) for atom in partner
            )
    if any(
        not math.isfinite(value) for value in row.values() if isinstance(value, float)
    ):
        raise ValueError("Nonfinite frame geometry")
    return row


def analyze_case(case: dict[str, Any], output: Path) -> dict[str, Any]:
    """Screen one replica with the same atom mapping and Tier A classifier.

    Args:
        case: Declared state, topology and production trajectory paths.
        output: New CSV path for the per-frame audit.

    Returns:
        Replica-level and 10 ns-block Tier counts and source provenance.
    """
    state = str(case["state"])
    topology_path = Path(case["topology"]).resolve(strict=True)
    trajectories = [Path(path).resolve(strict=True) for path in case["trajectories"]]
    topology = parmed.load_file(str(topology_path))
    indices = _indices(topology, state)
    universe = mda.Universe(str(topology_path), [str(path) for path in trajectories])
    expected = int(case["expected_frames"])
    if len(universe.trajectory) != expected or expected % 1000:
        raise ValueError("Trajectory frame count/cadence differs from manifest")
    counts: Counter[str] = Counter()
    blocks: dict[int, Counter[str]] = {}
    with output.open("w", newline="") as handle:
        writer: csv.DictWriter[str] | None = None
        for frame, ts in enumerate(universe.trajectory, start=1):
            data = _row(np.asarray(ts.positions, dtype=float), indices)
            data = {"frame": frame, "time_ps": float(ts.time), **data}
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(data))
                writer.writeheader()
            writer.writerow(data)
            tier = str(data["tier"])
            counts[tier] += 1
            block = (frame - 1) // 1000 + 1
            blocks.setdefault(block, Counter())[tier] += 1
    return {
        "name": case["name"],
        "state": state,
        "frames": expected,
        "tier_counts": dict(counts),
        "tier_a_fraction": counts["A"] / expected,
        "block_tier_a_fraction": {
            str(block): values["A"] / 1000 for block, values in sorted(blocks.items())
        },
        "topology_sha256": _hash(topology_path),
        "trajectory_sources": [
            {"path": str(path), "size_bytes": path.stat().st_size}
            for path in trajectories
        ],
        "frame_csv": str(output),
    }


def main() -> None:
    """Run the declared cases and write one comparison summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    manifest = json.loads(args.manifest.read_text())
    cases = manifest["cases"]
    if len({case["name"] for case in cases}) != len(cases):
        raise ValueError("Duplicate replica names in comparison manifest")
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    summaries = []
    for case in cases:
        result = analyze_case(case, args.output / f"{case['name']}.csv")
        summaries.append(result)
        LOGGER.info("%s Tier A %.3f", result["name"], result["tier_a_fraction"])
    (args.output / "summary.json").write_text(
        json.dumps(
            {"manifest_sha256": _hash(args.manifest), "replicas": summaries},
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
