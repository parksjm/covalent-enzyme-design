#!/usr/bin/env python3
"""Export the four completed native controls through solved PA6 conventions."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import shutil

from nyl12_covalent.control_inspection import export_control, publish_subset
from nyl12_covalent.control_pymol import write_views
from nyl12_covalent.design_pa6_mapping import _pdb_atom_line

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = (
    ROOT
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/dgx_wt_control_qualification_20261127"
)
CONFIG = ROOT / "design/nyl12_ligandmpnn_pilot/boltz_validation/configuration"


def main() -> None:
    """Build a new immutable curated control inspection bundle."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    stage = args.destination.with_name("." + args.destination.name + "_staging")
    if args.destination.exists() or stage.exists():
        raise FileExistsError(args.destination)
    stage.mkdir(parents=True)
    config = json.loads((CONFIG / "relaxation.json").read_text())
    maps = json.loads(
        (BASE.parent / "prediction_preparation/ligand_graph_manifest.json").read_text()
    )["ligand_atom_maps"]
    summary = []
    names = []
    for directory in sorted((BASE / "comparison_20261127").glob("*_contact_pairs")):
        receipt = json.loads((directory / "canonicalization_receipt.json").read_text())
        evidence = json.loads((directory / "pairing_evidence.json").read_text())
        if evidence["physical_alpha_to_beta"] != dict(zip("ABCD", "EFGH", strict=True)):
            raise ValueError(
                "Explicit nonidentity physical pair map needs adapter review"
            )
        source = Path(receipt["source"])
        predictor = directory.name.split("_")[0]
        seed = "20261126" if "20261126" in source.stem else "20261127"
        name = f"{predictor}_seed_{seed}"
        names.append(name)
        result = export_control(
            source,
            maps[predictor],
            Path(config["mol2"]),
            Path(config["mapping"]),
            Path(config["release"]),
            stage / name,
        )
        # Only curated payloads are published: retain conversion but omit converter scratch copies.
        component = stage / name / "amber_sequential_input/components"
        # Preserve converter evidence outside the public tree; never delete it.
        retained = BASE / "inspection_conversion_evidence" / name
        retained.parent.mkdir(exist_ok=True)
        if retained.exists():
            raise FileExistsError(retained)
        component.rename(retained)
        for filename in (
            "pairing_evidence.json",
            "geometry_all_48.csv",
            "geometry_assigned_4.csv",
        ):
            shutil.copy2(directory / filename, stage / name / filename)
        summary.append(
            {
                "model": name,
                "roundtrip_coordinate_component_error_A": result[
                    "roundtrip_max_coordinate_component_error_A"
                ],
                "roundtrip_metric_errors": result["metric_max_errors"],
                "atoms": result["raw_atom_count"],
            }
        )
    if len(names) != 4:
        raise ValueError("Exactly four models required")
    reference = ROOT / "analysis/native_thr_gate_calibration_20261123/inputs/n12c.cif"
    # Reuse audited crystal atom-site parser retained with the comparison evidence.
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "control_analysis", BASE / "comparison_20261127/analyze_controls.py"
    )
    analysis = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analysis)
    cols, rows = analysis.atom_loop(reference)
    records = []
    fit = set()
    for values in rows:
        d = dict(zip(cols, values, strict=True))
        if d["group_PDB"] != "ATOM" or d["label_alt_id"] not in (".", "?", "A"):
            continue
        element = d["type_symbol"]
        name = d["label_atom_id"]
        chain = d["auth_asym_id"]
        number = int(d["auth_seq_id"])
        if element in ("H", "D"):
            continue
        records.append(
            _pdb_atom_line(
                "ATOM",
                len(records) + 1,
                name,
                d["label_comp_id"],
                chain,
                number,
                *[float(d[f"Cartn_{a}"]) for a in "xyz"],
                element,
            )
        )
        if (
            chain == "A"
            and name == "CA"
            and 20 <= number <= 215
            and number not in (113, 186)
        ):
            fit.add(number)
    (stage / "crystal_reference.pdb").write_text("\n".join(records + ["END"]) + "\n")
    shutil.copy2(reference, stage / "crystal_reference.cif")
    write_views(stage, names, sorted(fit))
    shutil.copy2(
        BASE.parent / "prediction_preparation/ligand_graph_manifest.json",
        stage / "ligand_graph_manifest.json",
    )
    shutil.copy2(CONFIG / "whole_l4_source.mol2", stage / "whole_l4_source.mol2")
    shutil.copy2(
        CONFIG / "whole_l4_context_mapping.json",
        stage / "whole_l4_context_mapping.json",
    )
    shutil.copy2(
        BASE / "comparison_20261127/comparison_report.md",
        stage / "geometry_comparison.md",
    )
    (stage / "roundtrip_acceptance.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    (stage / "README.md").write_text("""# Native WT controls for inspection

These four original unrelaxed predictions fail protein-pocket/ligand geometry validation. The exports preserve that geometry. They are not accepted MD starting states or parameterized topologies.

On macOS, open a terminal in this directory and run `pymol load_all_controls.pml`, or in PyMOL `cd /path/to/this/directory` then `@load_all_controls.pml`. Per-model `view_*.pml` files recall individual scenes. Each model is one object; toggle `crystal_reference` to inspect the crystal overlay. White background, valence0, cyan protein/orange ligand carbons; N blue/O red. A single rigid alpha-chain core CA fit to crystalA is applied in the viewer only, excluding catalyticPhe113/Asn186 and all Thr/ligand atoms; poor pocket geometry remains visible. Crystal-defined orientation puts Thr below Phe(left)/Asn(right). Raw files and canonical files are not transformed.

Each model contains `raw_prediction.cif` (byte-identical original), `canonical_precursor.pdb/.cif` (biological A-D1-232, E-H233-328; whole ligandsI-L and source-MOL2 names), `canonical_mature_local.pdb` (same proteins but betaE-H1-96, required by the existing converter), and `amber_sequential_input/complex_heavy_modular.pdb` (established global residue order, protein1-1312; four P6N/P6I/P6I/P6C ligand copies1313-1328). This last file is heavy-atom LEaP input only: no hydrogens/protonation choice, topology, charge or tLEAP acceptance is implied. No MD, relaxation or force field was applied.

Biological Thr233 = mature betaE1 = original Amber929. BetaF1/G1/H1 correspond to Amber1025/1121/1217. Alpha chains each have232 residues; beta chains96. Never pass biological precursor numbering as Amber topology indices. `mapping_and_acceptance.json` contains every native→mature-local→biological→Amber atom/residue identity, exact raw coordinates and checked maximum PDB rounding/geometry errors; `amber_sequential_input/conversion.json` is the existing132-ligand-atom converter receipt. The graph release is the project's validated modular P6N/P6I/P6I/P6C release, not a newly fitted parameter set.

After pmemd/sander, use the existing `nyl12_covalent.design_pa6_mapping.export_canonical_heavy(topology,restart,conversion_json,output)` on the matching dry1328-residue topology to recover mature-local canonical IDs and original source atom names. Keep the matching topology and conversion receipt; do not infer identities from engine-written chain letters or atom indices. `nyl12_covalent.control_inspection.restore_modular` also restores a sequential PDB with unchanged heavy-residue inventory. Biological precursor display beta numbering is mature-local+232; retain both numbering conventions before RFD3 input preparation. A solvated topology must first use the established matching dry topology/stripping workflow, not this dry-only export API directly. No engine restart was run in this export test; the real sequential-heavy PDB roundtrip is tested on all four controls.

Geometry scans all3 internal amides ×4 wholeL4 ligands ×4 sites=48 comparisons/model. Four assigned rows maximize tierA/B count, then tier, then totald1 under a one-to-one ligand/site assignment; they do not fix one scissile amide and are not necessarily the closest ligand for each site individually. `mapping_and_acceptance.json:nearest_d1_per_site` separately exports the unconstrained nearest-d1 diagnostic, which can reuse a ligand and is not a catalytic assignment. All controls remain rejected; signedphi/nativeThr are diagnostics.

MANIFEST.sha256 covers exactly this curated, immutable publication. Original source/conversion evidence remains on the workstation; no private/Hessian/column caches are included. Download preserves this logical project path and never deletes laptop files.
""")
    publish_subset(stage, args.destination)
    logging.info("Published %s", args.destination)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
