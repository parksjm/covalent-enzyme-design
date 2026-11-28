"""Publish actual-case Amber diagnostics with graph-derived amide registers."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

from rdkit import Chem

from .negative_amide_gate import evaluate_signed_amide_gate
from .relay_any_amide import load_graph


def _write(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def annotate_amide_registers(
    rows: Sequence[Mapping[str, Any]],
    *,
    whole_context_mapping: Path,
    source_mol2: Path,
) -> list[dict[str, Any]]:
    """Copy geometry rows and annotate their source-graph L4 cleavage register.

    Args:
        rows: Geometry rows with source carbonyl-C/O and leaving-N names.
        whole_context_mapping: Explicit index-to-monomer context JSON.
        source_mol2: Audited source molecule whose atom indices define the map.

    Returns:
        Copies with monomer_boundary and hypothetical_fragment_lengths fields.

    Raises:
        ValueError: Mapping coverage, elements, graph or unit boundary is invalid.
    """
    graph = load_graph(source_mol2)
    mol = Chem.MolFromMol2File(str(source_mol2), removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError("Invalid source MOL2")
    assignments = json.loads(whole_context_mapping.read_text())["blocks"][0][
        "assignments"
    ]
    indices = [item["atom_index_1based"] for item in assignments]
    if len(indices) != len(set(indices)):
        raise ValueError("Duplicate context atom indices")
    index_map = {item["atom_index_1based"]: item for item in assignments}
    source_atoms = {atom.GetIdx() + 1: atom for atom in mol.GetAtoms()}
    if set(index_map) != set(source_atoms):
        raise ValueError("Context mapping must cover every source atom exactly")
    units = {}
    for index, atom in source_atoms.items():
        item = index_map[index]
        unit = item["residue_index_1based"]
        if type(unit) is not int or not 1 <= unit <= 4:
            raise ValueError("Expected four integer monomer identities")
        if item["element"] != atom.GetSymbol():
            raise ValueError("Context/source element mismatch")
        units[atom.GetProp("_TriposAtomName")] = unit
    if set(units.values()) != {1, 2, 3, 4} or len(graph.amides) != 3:
        raise ValueError("Expected complete L4 with three internal amides")
    registers = {}
    for amide in graph.amides:
        c, o, n = amide.carbon, amide.oxygen, amide.nitrogen
        unit = units[c]
        if units[o] != unit or units[n] != unit + 1:
            raise ValueError("Carbonyl O must share C unit and leaving N follow C")
        registers[c, o, n] = {
            "monomer_boundary": [unit, unit + 1],
            "hypothetical_fragment_lengths": [unit, 4 - unit],
            "register_is_product_hypothesis": True,
        }
    if {tuple(value["monomer_boundary"]) for value in registers.values()} != {
        (1, 2),
        (2, 3),
        (3, 4),
    }:
        raise ValueError("Source graph lacks distinct L4 internal boundaries")
    result = []
    for row in rows:
        key = tuple(row[k] for k in ("carbonyl_c", "carbonyl_o", "leaving_n"))
        if key not in registers:
            raise ValueError(f"Geometry amide absent from audited graph: {key}")
        result.append({**row, **registers[key]})
    return result


def process_design_case(
    *,
    case: str,
    target_phi_sign: int,
    topology: Path,
    before_restart: Path,
    after_restart: Path,
    conversion: Path,
    output_dir: Path,
    legacy_processor: Path,
    whole_context_mapping: Path,
    source_mol2: Path,
) -> dict[str, Any]:
    """Export one actual design case and gate all amides before assignment.

    Args:
        case: Actual case identity, preserved unchanged in every publication.
        target_phi_sign: Required +1 or -1, qualified against a declared reference.
        topology: Explicit dry Amber topology.
        before_restart: Explicit dry native initial restart.
        after_restart: Explicit dry native final restart.
        conversion: Explicit native-to-canonical atom mapping JSON.
        output_dir: New output directory; existing outputs are never overwritten.
        legacy_processor: Retained Python source exposing process(Namespace).
        whole_context_mapping: Source graph index-to-monomer JSON.
        source_mol2: Audited complete L4 connectivity source.

    Returns:
        Actual-case identity, explicit inputs/hashes and before/after gate reports.

    Raises:
        ValueError: Case identity, graph/register, geometry or input changes fail.
        FileExistsError: The output directory already exists.
        ImportError: The retained processor cannot be imported.
    """
    if type(target_phi_sign) is not int or target_phi_sign not in (-1, 1):
        raise ValueError("target_phi_sign must be integer +1 or -1; bool is invalid")
    if not case or Path(case).name != case or case in {".", ".."}:
        raise ValueError("Expected one actual case identity, not a path")
    inputs = [
        topology,
        before_restart,
        after_restart,
        conversion,
        legacy_processor,
        whole_context_mapping,
        source_mol2,
    ]
    hashes = {
        str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in inputs
    }
    # Validate graph and context before the legacy exporter creates any artifacts.
    annotate_amide_registers(
        [], whole_context_mapping=whole_context_mapping, source_mol2=source_mol2
    )
    if output_dir.exists():
        raise FileExistsError(output_dir)
    spec = importlib.util.spec_from_file_location(
        "actual_case_retained_processor", legacy_processor
    )
    if spec is None or spec.loader is None:
        raise ImportError(legacy_processor)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output_dir.mkdir(parents=True)
    canonical = output_dir / "canonical_before_after"
    module.process(
        argparse.Namespace(
            case=case,
            topology=topology,
            before_restart=before_restart,
            after_restart=after_restart,
            conversion=conversion,
            output_dir=canonical,
        )
    )
    result = {
        "case": case,
        "target_phi_sign": target_phi_sign,
        "input_sha256": hashes,
        "canonical_evidence": str(canonical.resolve()),
        "snapshots": {},
    }
    for state in ("before", "after"):
        geometry = canonical / f"{case}_{state}" / "geometry_any_amide.json"
        report = json.loads(geometry.read_text())
        rows = annotate_amide_registers(
            report["all_amide_site_rows"],
            whole_context_mapping=whole_context_mapping,
            source_mol2=source_mol2,
        )
        gate = {
            "case": case,
            "snapshot": state,
            "all_amide_site_rows": rows,
            "native_geometry_sha256": hashlib.sha256(geometry.read_bytes()).hexdigest(),
            **evaluate_signed_amide_gate(rows, target_phi_sign=target_phi_sign),
            "signed_ab_counts": {
                "positive": sum(
                    r["tier"] in {"A", "B"} and r["herschlag_phi"] > 0 for r in rows
                ),
                "negative": sum(
                    r["tier"] in {"A", "B"} and r["herschlag_phi"] < 0 for r in rows
                ),
            },
        }
        _write(output_dir / f"pose_gate_{state}.json", gate)
        result["snapshots"][state] = gate
    if any(
        hashlib.sha256(path.read_bytes()).hexdigest() != hashes[str(path.resolve())]
        for path in inputs
    ):
        raise ValueError("Frozen postprocessing input changed during execution")
    _write(output_dir / "case_postprocessing.json", result)
    (output_dir / "README.md").write_text(
        f"# Amber design-case inspection: {case}\n\n"
        "Canonical before/after artifacts retain the native legacy evidence verbatim.\n"
        "All48 graph-audited amide/site/ligand rows are annotated by actual register; "
        f"Explicit target phi sign {target_phi_sign:+d}: TierA/B filtering precedes one-to-one assignment. "
        "Fragment lengths are hypotheses, not experimentally identified products.\n"
    )
    return result
