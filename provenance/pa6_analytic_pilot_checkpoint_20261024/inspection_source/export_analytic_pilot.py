"""Export unchanged analytic-pilot coordinates with a fixed canonical pocket."""

from __future__ import annotations

import argparse
import gzip
import json
import logging
from pathlib import Path
import shutil
import sys

import numpy as np
from parmed.amber import AmberParm, NetCDFTraj

from pocket_export_io import record, render, zipped

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))
from nyl12_covalent.es_seed_display import _write_pdb  # noqa: E402


def _metrics(
    xyz: np.ndarray, roles: dict[str, int], hydroxyl_h: int
) -> dict[str, float]:
    pairs = {
        "d1_OG1_scissile_C": (roles["nucleophile"], roles["scissile_carbon"]),
        "d2_backbone_N_scissile_O": (
            roles["oxyanion_donor_backbone"],
            roles["scissile_oxygen"],
        ),
        "d3_sidechain_N_scissile_O": (
            roles["oxyanion_donor_sidechain"],
            roles["scissile_oxygen"],
        ),
        "scissile_C_N": (roles["scissile_carbon"], roles["leaving_nitrogen"]),
        "OG1_original_HG1": (roles["nucleophile"], hydroxyl_h),
        "terminal_N_original_HG1": (roles["thr_n_terminus"], hydroxyl_h),
    }
    return {
        name + "_A": float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))
        for name, (a, b) in pairs.items()
    }


def main() -> None:
    """Export verified PDB states and source hashes without changing calculation files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt-dir", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--single-frame-test-source", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    attempt, output = args.attempt_dir.resolve(), args.output.resolve()
    if args.single_frame_test_source is not None:
        starting = trajectory_path = endpoint = args.single_frame_test_source.resolve()
        label = "TEST_NOT_OPTIMIZATION"
    else:
        starting, trajectory_path, endpoint = (
            attempt / name
            for name in ("starting.ncrst", "optimization.nc", "optimized.ncrst")
        )
        label = "CONDITIONAL_ANALYTIC_ES_PILOT"
    tools = Path(__file__).resolve().parent
    paths = [
        args.topology.resolve(),
        args.inventory.resolve(),
        starting,
        trajectory_path,
        endpoint,
        Path(__file__),
        tools / "pocket_export_io.py",
        tools / "view_analytic_pilot.py",
        PROJECT / "src/nyl12_covalent/es_seed_display.py",
    ]
    sources = [record(path) for path in dict.fromkeys(paths)]
    output.mkdir(parents=True, exist_ok=False)
    (output / "source").mkdir()
    for path in paths[5:]:
        shutil.copy2(path, output / "source" / path.name)
    shutil.copy2(args.inventory, output / "source" / "inventory.json")
    (output / "source_manifest.json").write_text(json.dumps(sources, indent=2) + "\n")
    config = json.loads(args.inventory.read_text())
    variant = config["variants"][0]
    parm = AmberParm(str(args.topology), xyz=str(starting))
    initial = np.asarray(parm.coordinates).copy()
    if parm.get_coordinates().shape[0] != 1 or not np.isfinite(initial).all():
        raise ValueError("One finite starting coordinate frame required")
    qm_heavy = sorted(
        {
            atom["index"] - 1
            for group in variant["groups"]
            for atom in group["atoms"]
            if atom["element"] > 1
        }
    )
    distance2 = np.full(len(parm.atoms), np.inf)
    for index in qm_heavy:
        distance2 = np.minimum(
            distance2, np.sum((initial - initial[index]) ** 2, axis=1)
        )
    selected = {
        atom.residue.idx + 1
        for atom in parm.atoms
        if distance2[atom.idx] <= 36
        and (atom.residue.idx < 1313 or atom.residue.name == "WAT")
    }
    selected.add(1313)
    waters = sorted(
        number for number in selected if parm.residues[number - 1].name == "WAT"
    )
    for number in waters:
        atoms = parm.residues[number - 1].atoms
        if len(atoms) != 3 or sorted(atom.atomic_number for atom in atoms) != [1, 1, 8]:
            raise ValueError(f"Original water inventory is not O+2H: {number}")
    roles = variant["roles"]
    hydroxyl_h = [
        atom.idx + 1
        for atom in parm.atoms[roles["nucleophile"] - 1].bond_partners
        if atom.atomic_number == 1
    ]
    if len(hydroxyl_h) != 1:
        raise ValueError("Expected one original Thr hydroxyl hydrogen")
    start_text, mapping, start_error = render(parm, selected, waters, _write_pdb)
    (output / "pocket_initial.pdb").write_text(start_text)
    context, context_mapping, context_error = render(
        parm, set(range(1, 1314)), [], _write_pdb
    )
    zipped(output / "initial_fullprotein_ligand.pdb.gz", context)
    zipped(
        output / "atom_mapping.json.gz",
        json.dumps({"pocket": mapping, "context": context_mapping}, indent=2),
    )
    state_metadata = []
    model_texts = []
    reader = NetCDFTraj.open_old(str(trajectory_path))
    try:
        frames = np.asarray(reader.coordinates)
    finally:
        reader.close()
    if (
        frames.ndim != 3
        or frames.shape[1:] != initial.shape
        or not np.isfinite(frames).all()
    ):
        raise ValueError("Invalid actual NetCDF evaluation-state coordinates")
    if args.single_frame_test_source is not None and len(frames) != 1:
        raise ValueError("Single-frame test source contains multiple states")
    connectivity = ""
    for ordinal, xyz in enumerate(frames, 1):
        parm.coordinates = xyz
        text, observed_mapping, error = render(parm, selected, waters, _write_pdb)
        if observed_mapping != mapping:
            raise ValueError("Pocket atom/H mapping changed between states")
        atom_lines = [
            line
            for line in text.splitlines(True)
            if line.startswith(("ATOM  ", "HETATM"))
        ]
        if not connectivity:
            connectivity = "".join(
                line for line in text.splitlines(True) if line.startswith("CONECT")
            )
        model_texts.append(
            f"MODEL     {ordinal:4d}\n" + "".join(atom_lines) + "ENDMDL\n"
        )
        state_metadata.append(
            {
                "source_NetCDF_state": ordinal,
                "maximum_PDB_coordinate_component_error_A": error,
                "metrics": _metrics(xyz, roles, hydroxyl_h[0]),
            }
        )
    trajectory_text = (
        "REMARK OPTIMIZER EVALUATIONS; TRIAL/REPEATED STATES; NOT MD TIME\n"
        + "".join(model_texts)
        + connectivity
        + "END\n"
    )
    zipped(output / "optimization_evaluations.pdb.gz", trajectory_text)
    if (
        gzip.open(output / "optimization_evaluations.pdb.gz", "rt").read()
        != trajectory_text
    ):
        raise ValueError("Compressed evaluation trajectory failed round trip")
    end = AmberParm(str(args.topology), xyz=str(endpoint))
    if end.get_coordinates().shape[0] != 1:
        raise ValueError("One final endpoint frame required")
    endpoint_text, endpoint_mapping, endpoint_error = render(
        end, selected, waters, _write_pdb
    )
    if endpoint_mapping != mapping:
        raise ValueError("Endpoint changed static atom/H inventory")
    (output / "pocket_endpoint.pdb").write_text(endpoint_text)
    shutil.copy2(tools / "view_analytic_pilot.py", output / "view_analytic_pilot.py")
    report = {
        "label": label,
        "source_state_count": len(frames),
        "atom_count": len(mapping),
        "hydrogen_count": sum(row["atomic_number"] == 1 for row in mapping),
        "water_count": len(waters),
        "water_inventory": "Each retained original WAT has one O and two H",
        "static_topology_residues": sorted(selected),
        "qm_heavy_atom_indices": [i + 1 for i in qm_heavy],
        "selection": "Whole protein residues and whole WAT with any original atom within6A of any initial QM heavy atom; whole LIG always. Cartesian input coordinates; no coordinate edits or added H.",
        "display_only_exclusions": "All ions/non-WAT residues above topology1313 excluded from CHNOS renderer; raw calculations unchanged.",
        "states": state_metadata,
        "initial_metrics": _metrics(initial, roles, hydroxyl_h[0]),
        "endpoint_metrics": _metrics(np.asarray(end.coordinates), roles, hydroxyl_h[0]),
        "PDB_error_A": {
            "initial": start_error,
            "context": context_error,
            "endpoint": endpoint_error,
        },
        "all_original_selected_H_retained": True,
        "sources_unchanged": all(record(Path(row["path"])) == row for row in sources),
        "interpretation": "Actual NetCDF optimizer evaluation states may include trial and repeated coordinates; no MD time or accepted-step classification is assigned.",
    }
    if not report["sources_unchanged"]:
        raise ValueError("An export source changed during execution")
    (output / "export_report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    logging.info(
        "Exported %d original states, %d atoms, %d waters; %s",
        len(frames),
        len(mapping),
        len(waters),
        label,
    )


if __name__ == "__main__":
    main()
