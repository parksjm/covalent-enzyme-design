"""Export representative MD catalytic sites into the fixed n12c reference frame."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys
from typing import Any

import numpy as np
from parmed import load_file
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.catalytic_alignment import BACKBONE, FIT_RESIDUES, _pdb_line  # noqa: E402
from nyl12_covalent.catalytic_alignment_geometry import fit_coordinates  # noqa: E402
from nyl12_covalent.design_md_geometry import (  # noqa: E402
    frame_rows,
    image_near,
    index_mapping,
    native_reference,
)

LOGGER = logging.getLogger(__name__)
FOLLOW = ROOT / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup"
SITES = ("AE", "BF", "CG", "DH")


def sha(path: Path) -> str:
    """Return an existing artifact's content hash."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def align_site(
    reference: dict, protein: dict, ligand: dict
) -> tuple[dict, dict, dict[str, Any]]:
    """Rigidly align without using the held-out catalytic Thr in the fit.

    Args:
        reference: Canonical native heavy-atom map.
        protein: Mobile canonical atom map with elements and coordinates.
        ligand: Mobile named substrate coordinates.

    Returns:
        Aligned protein, ligand and reconstructible rigid transformation.
    """
    keys = [(auth, atom) for _, auth, _ in FIT_RESIDUES for atom in BACKBONE]
    rotation, center, target, rmsd = fit_coordinates(
        np.array([protein[key][1] for key in keys]),
        np.array([reference[key][1] for key in keys]),
    )
    aligned = {
        k: (e, (x - center) @ rotation + target) for k, (e, x) in protein.items()
    }
    lig = {k: (x - center) @ rotation + target for k, x in ligand.items()}
    return (
        aligned,
        lig,
        {
            "fit_residues": [113, 186, 272, 274],
            "fit_atoms": list(BACKBONE),
            "Thr233_excluded": True,
            "rotation_row_vector": rotation.tolist(),
            "mobile_center": center.tolist(),
            "reference_center": target.tolist(),
            "fit_RMSD_A": rmsd,
        },
    )


def write_pdb(path: Path, protein: dict, ligand: dict) -> None:
    """Save canonical residue numbering with no modified source coordinates.

    Args:
        path: Fresh output PDB.
        protein: Prealigned canonical protein map.
        ligand: Prealigned four-atom substrate core, or an empty dictionary.
    """
    names = {113: "PHE", 186: "ASN", 233: "THR", 234: "THR", 272: "ASP", 274: "ASP"}
    lines = ["REMARK Rigidly aligned to n12c chain A; Thr233 excluded from fit.\n"]
    for serial, ((resid, atom), (element, xyz)) in enumerate(
        sorted(protein.items()), 1
    ):
        lines.append(
            _pdb_line(serial, "ATOM", atom, names[resid], "A", resid, xyz, element)
        )
    lines.append("TER\n")
    for atom, xyz in ligand.items():
        lines.append(_pdb_line(len(lines), "HETATM", atom, "PA6", "I", 1, xyz, atom[0]))
    lines.append("END\n")
    with path.open("x") as stream:
        stream.writelines(lines)


def export(output: Path) -> dict[str, Any]:
    """Export fixed-rule medoid/closest frames and initial minimized controls.

    Args:
        output: Existing diagnostic bundle containing selection.json.

    Returns:
        Source-bound snapshot manifest; all original files remain unchanged.
    """
    selection = json.loads((output / "selection.json").read_text())
    reference_path = ROOT / "cif/n12c.cif"
    ref = native_reference(reference_path)
    write_pdb(output / "reference.pdb", ref, {})
    snapshots = []
    frame_dir = output / "coordinate_evidence"
    frame_dir.mkdir(exist_ok=False)
    for dataset in ("WT", "broader_mask07"):
        manifest_path = FOLLOW / "systems" / dataset / "datasets.json"
        entry = json.loads(manifest_path.read_text())["datasets"][dataset]
        topology = (manifest_path.parent / entry["topology"]["path"]).resolve(
            strict=True
        )
        restart = (manifest_path.parent / entry["restart"]["path"]).resolve(strict=True)
        if (
            sha(topology) != entry["topology"]["sha256"]
            or sha(restart) != entry["restart"]["sha256"]
        ):
            raise ValueError("Prepared topology/restart identity changed")
        source = Path(
            json.loads((FOLLOW / "selection.json").read_text())["source_selection"][
                dataset
            ]
        )
        parm = load_file(str(topology), xyz=str(restart))
        mapping = index_mapping(parm, source)
        requests = [("initial", {"site": "AE", "time_ps": None})]
        requests += [
            ("typical" if role == "AE_median_medoid" else "closest", row)
            for role, row in selection[dataset].items()
        ]
        for role, request in requests:
            site = SITES.index(request["site"])
            if role == "initial":
                xyz, box = np.asarray(parm.coordinates), np.asarray(parm.box)
                frame_hash = hashlib.sha256(xyz.tobytes()).hexdigest()
                source_file, frame_index = restart, None
            else:
                source_file = Path(request["trajectory_path"])
                frame_index = request["frame_index_zero_based"]
                with netcdf_file(source_file, "r", mmap=True) as nc:
                    if float(nc.variables["time"][frame_index]) != request["time_ps"]:
                        raise ValueError("Selected frame/time mismatch")
                    xyz = np.array(nc.variables["coordinates"][frame_index])
                    box = np.r_[
                        nc.variables["cell_lengths"][frame_index],
                        nc.variables["cell_angles"][frame_index],
                    ]
                frame_hash = hashlib.sha256(xyz.tobytes()).hexdigest()
            metrics = frame_rows(xyz, box, mapping, ref)[site]
            if role != "initial":
                for name, expected in request["metrics"].items():
                    if not np.isclose(metrics[name], expected, rtol=0, atol=1e-7):
                        raise ValueError(
                            f"Export frame does not reproduce stored {name}"
                        )
                if metrics["native_pass"] != request["native_pass"]:
                    raise ValueError("Export native gate differs from source")
            indices = mapping[site]["protein"]
            origin = xyz[indices[(233, "OG1")][0]]
            protein = {
                key: (element, pos)
                for (key, (_, element)), pos in zip(
                    indices.items(),
                    image_near(xyz[[v[0] for v in indices.values()]], origin, box),
                    strict=True,
                )
            }
            # Include genuine MD hydrogens of Thr233/234 without constructing any.
            for resid in (233, 234):
                amber_res = 929 + 96 * site + resid - 233
                for atom in parm.residues[amber_res - 1].atoms:
                    if atom.atomic_number == 1:
                        point = image_near(xyz[[atom.idx]], origin, box)[0]
                        protein[(resid, atom.name)] = ("H", point)
            ligand_ids = mapping[site]["ligand"]
            ligand = dict(
                zip(
                    ligand_ids,
                    image_near(xyz[list(ligand_ids.values())], origin, box),
                    strict=True,
                )
            )
            aligned, ligand_aligned, fit = align_site(ref, protein, ligand)
            identity = f"{dataset}_{role}"
            np.savez(
                frame_dir / f"{identity}.npz",
                keys=np.array([f"{r}:{a}" for r, a in protein]),
                coordinates=np.array([x for _, x in protein.values()]),
            )
            filename = f"{identity}.pdb"
            write_pdb(output / filename, aligned, ligand_aligned)
            snapshots.append(
                {
                    "id": identity,
                    "dataset": dataset,
                    "role": role,
                    "file": filename,
                    "label": f"{dataset} {role}",
                    "site": request["site"],
                    "time_ns": None if role == "initial" else request["time_ps"] / 1000,
                    "frame_index_zero_based": frame_index,
                    "source_file": str(source_file),
                    "coordinate_sha256": frame_hash,
                    "pdb_sha256": sha(output / filename),
                    "metrics": metrics,
                    "fit": fit,
                    "selection_rule": "Original minimized control"
                    if role == "initial"
                    else request["rule"],
                }
            )
            LOGGER.info(
                "Exported %s (%s, time %s)",
                identity,
                request["site"],
                request["time_ps"],
            )
    result = {
        "reference": "reference.pdb",
        "reference_source": str(reference_path),
        "reference_source_sha256": sha(reference_path),
        "selection_sha256": sha(output / "selection.json"),
        "export_script_sha256": sha(Path(__file__)),
        "scalar_report_sha256": sha(output / "scalar_gate_report.json"),
        "alignment": "Same 16 scaffold atoms for every snapshot; catalytic Thr excluded",
        "snapshots": snapshots,
    }
    (output / "snapshot_manifest.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    return result


def main() -> None:
    """Parse the output bundle and export the selected source frames."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    export(args.output.resolve(strict=True))


if __name__ == "__main__":
    main()
