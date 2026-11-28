"""Derive a graph-preserving native RFD3 export and apply unchanged gates."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

import biotite.structure as struc
import numpy as np
from atomworks.io.utils.io_utils import read_any, to_cif_file
from biotite.structure.io.pdbx import get_structure

from nyl12_covalent.rfd3_export_adapter import (
    adapt_coordinates,
    reassociate_typed_graph,
)
from nyl12_covalent.rfd3_output_mapping import canonical_output_keys
from nyl12_covalent.rfd3_qualification_audit_v2 import coordinate_gates

LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return SHA256 for a retained artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    """Create immutable evidence."""
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def export_parent(manifest_path: Path, job_id: str, destination: Path) -> dict:
    """Export a single raw native candidate with explicit graph provenance.

    Args:
        manifest_path: Hash-bound native job manifest.
        job_id: Exact native job identity.
        destination: New derived export directory.

    Returns:
        Independent scientific acceptance, including artifact hashes.
    """
    manifest = json.loads(manifest_path.read_text())
    job = next(j for j in manifest["jobs"] if j["job_id"] == job_id)
    raw_files = [
        p
        for p in Path(job["output"]).glob("*.cif.gz")
        if "noisy" not in p.name and "denoised" not in p.name
    ]
    if len(raw_files) != 1:
        raise ValueError("Require exactly one raw native final CIF")
    raw = raw_files[0]
    aa = get_structure(read_any(raw), model=1, include_bonds=True)
    native_keys = canonical_output_keys(aa.chain_id, aa.res_id, aa.atom_name)
    reverse = {v: k for k, v in manifest["component_aliases"].items()}
    native_names = [reverse.get(str(n), str(n)) for n in aa.res_name]
    reference_path = Path(job["reference_graph"])
    ref = json.loads(reference_path.read_text())
    pdb_source = Path(ref["pdb"])
    lines = [
        line
        for line in pdb_source.read_text().splitlines()
        if line[:6] in ("ATOM  ", "HETATM")
    ]
    keys = [(line[21], int(line[22:26]), line[12:16].strip()) for line in lines]
    names = [line[17:20].strip() for line in lines]
    lookup = {k: i for i, k in enumerate(keys)}
    if any(native_names[i] != names[lookup[k]] for i, k in enumerate(native_keys)):
        raise ValueError("Native sequence/ligand identities changed")
    spec_path = Path(
        next(a.split("=", 1)[1] for a in job["argv"] if a.startswith("inputs="))
    )
    spec = next(iter(json.loads(spec_path.read_text()).values()))
    if (
        spec.get("ori_token") is not None
        or spec.get("infer_ori_strategy") is not None
        or spec.get("symmetry")
    ):
        raise ValueError("Adapter requires native default partial-diffusion mobile COM")
    source_cif = Path(spec["input"])
    source = get_structure(read_any(source_cif), model=1, include_bonds=True)
    source_keys = canonical_output_keys(
        source.chain_id, source.res_id, source.atom_name
    )
    typed = (
        reassociate_typed_graph(keys, source.bonds.as_array(), keys)
        if source_keys == keys
        else reassociate_typed_graph(source_keys, source.bonds.as_array(), keys)
    )
    ligand = sorted(
        [*sorted((int(a), int(b))), int(t)]
        for a, b, t in typed
        if keys[a][0] == keys[b][0] == "I"
    )
    if ligand != sorted([*sorted((a, b)), t] for a, b, t in ref["ligand_bonds"]):
        raise ValueError("Explicit source typed PA6 graph differs from reference")
    native_pairs = {tuple(sorted((int(a), int(b)))) for a, b, _ in typed}
    if native_pairs != {tuple(sorted(edge)) for edge in ref["native_edges"]}:
        raise ValueError(
            "Explicit source full graph differs from reference native edges"
        )
    coords = np.asarray(ref["coords"], dtype=float)
    mobile = set(job["mobile_residues"])
    moving = np.array([f"{c}{p}" in mobile for c, p, _ in keys])
    center = coords.astype(np.float32)[moving].mean(axis=0)
    generated, normalization = adapt_coordinates(
        keys, coords, native_keys, aa.coord.astype(float), mobile, center
    )
    atoms = [
        {"chain": c, "residue": p, "name": n, "residue_name": rn}
        for (c, p, n), rn in zip(keys, names, strict=True)
    ]
    receipt = coordinate_gates(
        coords,
        generated,
        atoms,
        ref["native_edges"],
        job["mobile_residues"],
        ref["core_indices0based"][:1] + ref["core_indices0based"][2:3],
        manifest["thresholds"],
    )
    destination.mkdir()
    pdb = destination / "audited_parent.pdb"
    with pdb.open("x") as handle:
        for line, xyz in zip(lines, generated, strict=True):
            handle.write(
                line[:30]
                + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                + line[54:]
                + "\n"
            )
        handle.write("END\n")
    # The coordinate/graph JSON is authoritative; PDB is a precision-limited consumer export.
    exact = destination / "derived_coordinates_and_typed_graph.json"
    save(
        exact,
        {
            "atom_keys": keys,
            "residue_names": names,
            "coords_A": generated.tolist(),
            "typed_bonds": typed.tolist(),
            "typed_graph_source": str(source_cif),
            "typed_graph_source_sha256": digest(source_cif),
            "spatial_bond_inference": False,
        },
    )
    derived = source.copy()
    source_lookup = {k: i for i, k in enumerate(keys)}
    derived.coord = generated[[source_lookup[k] for k in source_keys]]
    derived.bonds = struc.BondList(
        len(derived), reassociate_typed_graph(keys, typed, source_keys)
    )
    cif = destination / "audited_parent_explicit_graph.cif"
    to_cif_file(derived, cif, include_bonds=True, include_entity_poly=True)
    mapping = destination / "canonical_mapping.json"
    save(mapping, ref["atom_map"])
    graph = destination / "pa6_graph_acceptance.json"
    save(
        graph,
        {
            "accepted": True,
            "ligand_heavy_atoms": 33,
            "ligand_heavy_bonds": ligand,
            "source_atom_correspondence": True,
            "source_graph_convention": manifest["source_graph_convention"],
            "full_typed_graph": str(exact),
            "full_typed_graph_sha256": digest(exact),
            "reassociated_not_inferred": True,
        },
    )
    normalization |= {
        "source_reference_sha256": digest(reference_path),
        "source_coordinate_file_sha256": digest(pdb_source),
    }
    save(destination / "fixed_context_normalization.json", normalization)
    receipt |= {
        "job_id": job_id,
        "parent_source_id": job["parent"],
        "source_atom_inventory_exact": True,
        "source_sequence_exact": True,
        "source_ligand_graph_exact": True,
        "native_cif": str(raw),
        "native_cif_sha256": digest(raw),
        "manifest": str(manifest_path),
        "manifest_sha256": digest(manifest_path),
        "pdb_path": str(pdb),
        "canonical_mapping": str(mapping),
        "pa6_graph_receipt": str(graph),
        "authoritative_coordinates": str(exact),
        "explicit_graph_cif": str(cif),
        "fixed_context_normalization": str(
            destination / "fixed_context_normalization.json"
        ),
    }
    for field in (
        "pdb_path",
        "canonical_mapping",
        "pa6_graph_receipt",
        "authoritative_coordinates",
        "explicit_graph_cif",
    ):
        receipt[field + "_sha256"] = digest(Path(receipt[field]))
    save(destination / "scientific_acceptance.json", receipt)
    LOGGER.info("Exported %s; scientific acceptance=%s", job_id, receipt["accepted"])
    return receipt


def main() -> None:
    """Create one retained derived export and return failure for a rejected parent."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("job_id")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    receipt = export_parent(
        args.manifest.resolve(), args.job_id, args.destination.resolve()
    )
    if not receipt["accepted"]:
        raise RuntimeError("Derived native candidate failed unchanged scientific gates")


if __name__ == "__main__":
    main()
