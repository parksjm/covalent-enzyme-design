"""Screen existing primary MD for observed N-alpha/water/leaving-N bridges."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np
import parmed
from parmed.amber import NetCDFRestart
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.catalytic_alignment import SITE_RESIDUES, _cif_atoms  # noqa:E402
from nyl12_covalent.nalpha_water_relay import (  # noqa:E402
    angle,
    candidate_rank,
    cell_matrix,
    close_waters,
    hbond,
    minimum_image,
    native_geometry,
    relay_class,
)

LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return streaming SHA-256 of an immutable source artifact."""
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def mapping(topology: Any, config: dict) -> dict[str, Any]:
    """Map named residues and topology-owned hydrogens without fixed atom indices."""
    result = {}
    for label, (resid, name) in config["atoms"].items():
        atoms = [a for a in topology.residues[resid - 1].atoms if a.name == name]
        if len(atoms) != 1:
            raise ValueError(f"Ambiguous mapping for {label}")
        result[label] = atoms[0].idx
    nalpha = topology.atoms[result["nalpha"]]
    result["nalpha_h"] = [a.idx for a in nalpha.bond_partners if a.atomic_number == 1]
    if len(result["nalpha_h"]) != 2:
        raise ValueError("Primary ES requires neutral catalytic NH2")
    for resid, expected in config["residue_names"].items():
        if topology.residues[int(resid) - 1].name != expected:
            raise ValueError(f"Unexpected residue {resid}")
    ligand = {a.idx for r in topology.residues[1312:1316] for a in r.atoms}
    result["ligand_bonds"] = [
        (b.atom1.idx, b.atom2.idx)
        for b in topology.bonds
        if b.atom1.idx in ligand and b.atom2.idx in ligand
    ]
    result["waters"] = []
    for residue in topology.residues:
        if residue.name != "WAT":
            continue
        oxygen = [a for a in residue.atoms if a.atomic_number == 8]
        if len(oxygen) != 1:
            raise ValueError("Unexpected water oxygen inventory")
        hydrogens = [a.idx for a in oxygen[0].bond_partners if a.atomic_number == 1]
        if len(hydrogens) != 2 or len(residue.atoms) != 3:
            raise ValueError("Water must retain two oxygen-owned hydrogens")
        result["waters"].append((residue.idx + 1, oxygen[0].idx, *hydrogens))
    return result


def scan_case(case: dict, config: dict, output: Path) -> tuple[dict, list[dict]]:
    """Screen all cached Tier A/B positive-phi frames in one declared replica."""
    topology = parmed.load_file(case["topology"])
    mapped = mapping(topology, {**config, "residue_names": case["residue_names"]})
    reference = _cif_atoms(
        Path(config["reference_cif"]), "A", {auth for _, auth, _ in SITE_RESIDUES}
    )
    with Path(case["geometry_csv"]).open(newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    if len(all_rows) != case["expected_frames"]:
        raise ValueError("Cached geometry frame count mismatch")
    eligible = {
        int(row["frame"]): row
        for row in all_rows
        if row["tier"] in ("A", "B") and float(row["projected_phi"]) > 0
    }
    counts: Counter = Counter(
        total_frames=len(all_rows), tier_ab_positive_phi=len(eligible)
    )
    waters = np.asarray(mapped["waters"], dtype=int)
    bonds = np.asarray(mapped["ligand_bonds"], dtype=int)
    rows, gates = [], []
    offset = 0
    for trajectory in case["trajectories"]:
        with netcdf_file(trajectory, mmap=True) as source:
            nframes = source.variables["coordinates"].shape[0]
            for local in range(nframes):
                frame = offset + local + 1
                if frame not in eligible:
                    continue
                xyz = np.array(source.variables["coordinates"][local], dtype=float)
                lengths = np.array(source.variables["cell_lengths"][local])
                angles = np.array(source.variables["cell_angles"][local])
                cell = cell_matrix(lengths, angles)
                time = float(source.variables["time"][local])
                cached = eligible[frame]
                if abs(time - float(cached["time_ps"])) > 1e-3:
                    raise ValueError("Cached frame time differs from trajectory")
                native = native_geometry(
                    xyz, topology, reference, cell, mapped["nalpha"]
                )
                bond_lengths = np.linalg.norm(
                    minimum_image(xyz[bonds[:, 0]] - xyz[bonds[:, 1]], cell), axis=1
                )
                limits = np.array(
                    [
                        1.3
                        if topology.atoms[a].atomic_number == 1
                        or topology.atoms[b].atomic_number == 1
                        else 1.85
                        for a, b in bonds
                    ]
                )
                intact = bool(np.all((bond_lengths > 0.65) & (bond_lengths <= limits)))
                counts["intact_pa6"] += int(intact)
                counts["native_pass"] += int(native["native_pass"])
                gates.append(
                    {
                        "frame": frame,
                        "time_ps": time,
                        "tier": cached["tier"],
                        "intact_pa6": intact,
                        **native,
                    }
                )
                if not intact or not native["native_pass"]:
                    continue
                counts["joint_eligible"] += 1
                near = close_waters(xyz[waters[:, 1]], xyz[mapped["nalpha"]], cell, 3.5)
                near2 = close_waters(
                    xyz[waters[near, 1]], xyz[mapped["nleave"]], cell, 3.5
                )
                matches = waters[near[near2]]
                counts["heavy_bridge_frames"] += int(bool(len(matches)))
                local_indices = [
                    value for value in mapped.values() if isinstance(value, int)
                ]
                local_points = minimum_image(
                    xyz[local_indices] - xyz[mapped["carbon"]], cell
                )
                local_xyz = dict(zip(local_indices, local_points, strict=True))
                for resid, oxygen, h1, h2 in matches:
                    outgoing = [
                        hbond(xyz, oxygen, h, mapped["nleave"], cell) for h in (h1, h2)
                    ]
                    best_i = min(
                        range(2),
                        key=lambda i: (not outgoing[i]["pass"], outgoing[i]["ha_A"]),
                    )
                    second = outgoing[best_i]
                    incoming = [
                        hbond(xyz, mapped["nalpha"], h, oxygen, cell)
                        for h in mapped["nalpha_h"]
                    ]
                    first = min(incoming, key=lambda r: (not r["pass"], r["ha_A"]))
                    kind = relay_class(bool(first["pass"]), bool(second["pass"]))
                    row = {
                        "case": case["name"],
                        "frame": frame,
                        "time_ps": time,
                        "trajectory": trajectory,
                        "local_frame": local + 1,
                        "tier": cached["tier"],
                        "relay_class": kind,
                        "water_resid": int(resid),
                        "water_o_index1": int(oxygen + 1),
                        "water_donor_h_index1": int((h1, h2)[best_i] + 1),
                        "nalpha_water_A": first["da_A"],
                        "water_nleave_A": second["da_A"],
                        "nalpha_h_to_water_A": first["ha_A"],
                        "nalpha_h_to_water_angle_deg": first["angle_deg"],
                        "water_h_to_nleave_A": second["ha_A"],
                        "water_h_to_nleave_angle_deg": second["angle_deg"],
                        **{
                            k: float(cached[k])
                            for k in ("d1", "d2", "d3", "attack", "projected_phi")
                        },
                        "cb_og1_carbonyl_angle_deg": angle(
                            *[local_xyz[mapped[k]] for k in ("cb", "og1", "carbon")]
                        ),
                        "nleave_scissile_c_c1_angle_deg": angle(
                            *[
                                local_xyz[mapped[k]]
                                for k in ("carbon", "nleave", "leaving_c1")
                            ]
                        ),
                        **native,
                    }
                    rows.append(row)
                    counts[kind + "_events"] += 1
                if frame % 500 == 0:
                    LOGGER.info("%s frame %d: %s", case["name"], frame, dict(counts))
            offset += nframes
    if offset != len(all_rows):
        raise ValueError("Trajectory count mismatch")
    with (output / f"{case['name']}_gates.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(gates[0])
            if gates
            else ["frame", "time_ps", "tier", "intact_pa6", "native_pass"],
        )
        writer.writeheader()
        writer.writerows(gates)
    return {
        "case": case["name"],
        "counts": dict(counts),
        "topology": case["topology"],
        "topology_sha256": digest(Path(case["topology"])),
    }, rows


def extract(row: dict, config: dict, output: Path, label: str) -> dict:
    """Write one exact source frame as Amber NetCDF restart with topology mapping."""
    case = next(c for c in config["cases"] if c["name"] == row["case"])
    topology = parmed.load_file(case["topology"])
    with netcdf_file(row["trajectory"], mmap=True) as source:
        i = row["local_frame"] - 1
        xyz = np.array(source.variables["coordinates"][i], dtype=float)
        box = np.concatenate(
            [np.array(source.variables[k][i]) for k in ("cell_lengths", "cell_angles")]
        )
    path = output / f"{label}.ncrst"
    restart = NetCDFRestart.open_new(
        str(path),
        len(xyz),
        True,
        False,
        title=f"Unmodified MD {row['case']} frame {row['frame']}",
    )
    restart.coordinates, restart.box = xyz, box
    restart.close()
    verify = NetCDFRestart.open_old(str(path))
    try:
        if not np.array_equal(xyz, np.asarray(verify.coordinates).reshape(-1, 3)):
            raise ValueError("Restart changed raw source coordinates")
    finally:
        verify.close()
    mapped = mapping(topology, {**config, "residue_names": case["residue_names"]})
    atoms = {
        key: {
            "index1": value + 1,
            "resid": topology.atoms[value].residue.idx + 1,
            "resname": topology.atoms[value].residue.name,
            "name": topology.atoms[value].name,
        }
        for key, value in mapped.items()
        if isinstance(value, int)
    }
    record = {
        **row,
        "restart": str(path),
        "restart_sha256": digest(path),
        "coordinates_float64_sha256": hashlib.sha256(
            xyz.astype("<f8").tobytes()
        ).hexdigest(),
        "topology": case["topology"],
        "topology_sha256": digest(Path(case["topology"])),
        "atom_mapping": atoms,
        "coordinates_unmodified": True,
        "relay_water_atoms_index1": [
            a.idx + 1 for a in topology.residues[row["water_resid"] - 1].atoms
        ],
    }
    (output / f"{label}.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    """Run the explicit bounded primary screen and extract ranked observed frames."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    args.output.mkdir(exist_ok=False)
    snapshot = args.output / "source_snapshot"
    snapshot.mkdir()
    for name in config["source_files"]:
        shutil.copy2(name, snapshot / Path(name).name)
    provenance = {str(path): digest(Path(path)) for path in config["hash_inputs"]}
    (args.output / "input_hashes.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    results, rows = [], []
    for case in config["cases"]:
        result, candidates = scan_case(case, config, args.output)
        results.append(result)
        rows.extend(candidates)
        (args.output / "screen.partial.json").write_text(
            json.dumps(results, indent=2) + "\n"
        )
        LOGGER.info("Completed %s: %s", case["name"], result["counts"])
    rows.sort(key=candidate_rank)
    with (args.output / "candidates.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]) if rows else ["case", "frame", "relay_class"],
        )
        writer.writeheader()
        writer.writerows(rows)
    selections = []
    for row in rows:
        if any(
            row["case"] == s["case"] and abs(row["frame"] - s["frame"]) < 100
            for s in selections
        ):
            continue
        label = "selected" if not selections else "reserve"
        selections.append(extract(row, config, args.output, label))
        if len(selections) == 2:
            break
    summary = {
        "cases": results,
        "candidate_events": len(rows),
        "selected_count": len(selections),
        "selection_classes": [s["relay_class"] for s in selections],
        "acceptance": {
            "tier": "cached A/B and positive signed phi",
            "native": config["native_gates"],
            "relay": "O within3.5A both N; water H...N<=2.5A and O-H...N>=120deg",
            "complete_relay_proven": False,
        },
        "limitations": "Observed ES geometry only; intact amide N proton affinity and product path untested. Neutral NH2 lacks future transferred HG1; missing incoming H bond does not reject leaving-leg candidate.",
        "input_hashes": provenance,
        "python_executable": sys.executable,
    }
    (args.output / "screen.json").write_text(json.dumps(summary, indent=2) + "\n")
    LOGGER.info(
        "Screen complete: %d water events, %d extracted frames",
        len(rows),
        len(selections),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
