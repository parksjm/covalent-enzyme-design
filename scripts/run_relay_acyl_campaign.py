#!/usr/bin/env python3
"""Gate balanced product relaxation and a reverse asynchronous acylation path."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys
from typing import Any

import numpy as np
from parmed.amber import Rst7

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.reaction_bias import (  # noqa: E402
    CoordinationCV,
    CoordinationDifferenceCV,
    DistanceCV,
    HarmonicTarget,
    LinearCombinationCV,
    WeightedCV,
    evaluate_bias,
)
from nyl12_covalent.relay_case import RelayCase, load_case, state_record, write_json  # noqa: E402
from nyl12_covalent.relay_stage import run_stage  # noqa: E402

LOGGER = logging.getLogger(__name__)


def classify(model: RelayCase, xyz: np.ndarray) -> dict[str, Any]:
    """Apply explicit covalent-product and reactant geometric/inventory gates.

    Args:
        model: Explicit unchanged scientific case.
        xyz: Full-system coordinates in topology order, angstrom.

    Returns:
        Measured geometry, all-H audit and provisional chemistry acceptance.
    """
    state = state_record(model, xyz)
    owners: dict[int, list[int]] = {}
    ambiguous = []
    unassigned = []
    for h in state["proton_owners"]:
        if h["classification"] == "shared_ambiguous":
            ambiguous.append(h["hydrogen_1based"])
        if h["classification"] == "unassigned":
            unassigned.append(h["hydrogen_1based"])
        if h["nearest_distance_A"] <= 1.25:
            owners.setdefault(h["nearest_qm_heavy_1based"], []).append(
                h["hydrogen_1based"]
            )
    m = model.atom_map
    count = {
        key: len(owners.get(m[key], []))
        for key in (
            "thr_n",
            "thr_og",
            "leaving_n",
            "oxygen",
            "water_o",
            "ash_oxygen",
        )
    }
    asp_h = len(owners.get(m["asp_o1"], []) + owners.get(m["asp_o2"], []))
    acid_h = sum(len(owners.get(i, [])) for i in (14578, 14579))
    reactive = {
        m[key] for key in ("thr_n", "thr_og", "leaving_n", "oxygen", "water_o")
    } | {14578, 14579, m["asp_o1"], m["asp_o2"]}
    changed_spectators = []
    for index in model.qm:
        atom = model.topology.atoms[index]
        if atom.atomic_number == 1 or index + 1 in reactive:
            continue
        expected = sum(a.atomic_number == 1 for a in atom.bond_partners)
        if len(owners.get(index + 1, [])) != expected:
            changed_spectators.append(index + 1)
    spectators = (
        count["water_o"] == 2
        and acid_h == 1
        and asp_h == 0
        and not unassigned
        and not changed_spectators
    )
    inventory = (count["thr_n"], count["thr_og"], count["leaving_n"], count["oxygen"])
    product = (
        1.25 <= state["rOG1C"] <= 1.65
        and state["rCN"] >= 2.3
        and 1.15 <= state["rCO"] <= 1.32
        and inventory == (2, 0, 2, 0)
        and spectators
        and not ambiguous
    )
    reactant = (
        state["rOG1C"] >= 2.5
        and 1.2 <= state["rCN"] <= 1.6
        and 1.15 <= state["rCO"] <= 1.32
        and inventory == (2, 1, 1, 0)
        and spectators
        and not ambiguous
    )
    return {
        "product_geometry_retained": product,
        "reactant_geometry": reactant,
        "geometric_counts": count,
        "ash274_protons": acid_h,
        "asp272_protons": asp_h,
        "ambiguous_protons": ambiguous,
        "unassigned_protons": unassigned,
        "spectator_heavy_atoms_with_changed_H_counts": changed_spectators,
        "state": state,
        "geometry_is_not_hessian_stationarity_verification": True,
    }


def asynchronous_coordinate(model: RelayCase) -> LinearCombinationCV:
    """Build one scalar with orthogonal attack/proton/scission motion unrestrained."""
    m = {key: value - 1 for key, value in model.atom_map.items()}
    return LinearCombinationCV(
        (
            WeightedCV(DistanceCV(m["thr_og"], m["carbon"]), 1.0),
            WeightedCV(DistanceCV(m["carbon"], m["leaving_n"]), -1.0),
            WeightedCV(
                CoordinationDifferenceCV(
                    m["leaving_n"], m["thr_og"], tuple(model.hydrogens)
                ),
                1.0,
            ),
        ),
        output_unit="angstrom",
    )


def run(case_path: Path) -> None:
    """Run one endpoint test, then conditionally connect product back toward E-S.

    Args:
        case_path: Explicit campaign case and fixed source identities.
    """
    case_path = case_path.resolve()
    folder = case_path.parent
    model = load_case(case_path)
    source = folder / "01_endpoint_inputs/acyl_guess.ncrst"
    xyz = np.asarray(Rst7.open(str(source)).coordinates).reshape(-1, 3)
    m = {key: value - 1 for key, value in model.atom_map.items()}
    leave = CoordinationCV(m["leaving_n"], tuple(model.hydrogens))
    leaving_count = evaluate_bias(xyz, [HarmonicTarget(leave, 0, 0)]).values[0]
    terms = [
        HarmonicTarget(DistanceCV(m["thr_og"], m["carbon"]), 1.4, 200),
        HarmonicTarget(DistanceCV(m["carbon"], m["leaving_n"]), 3.0, 50),
        HarmonicTarget(leave, leaving_count, 100),
    ]
    write_json(
        folder / "protocol.json",
        {
            "hypothesis": "A balanced acyl-enzyme endpoint may persist without a stable TI1",
            "source_case": str(case_path),
            "primary_water_screen": "no eligible relay",
            "endpoint_test_independent_of_relay_screen": True,
            "seed_preoptimization_targets": [t.to_metadata() for t in terms],
            "seed_max_steps": 1500,
            "unbiased_max_steps": 2500,
            "reverse_path": {
                "launch_only_if_unbiased_product_force_converged_and_retained": True,
                "coordinate": HarmonicTarget(
                    asynchronous_coordinate(model), 0, 0
                ).to_metadata()["cv"],
                "formula": "r(OG1,C)-r(C,Nleave)+1A*(CN(OG1,H)-CN(Nleave,H)); all54realQMH",
                "windows": 36,
                "spring_kcal_mol_A2": 100,
                "max_steps_per_window": 1500,
                "Nalpha_occupancy_restrained": False,
                "independent_distances_or_proton_identities_restrained": False,
                "maximum_is_not_a_verified_TS_or_PMF": True,
            },
            "completion_and_chemistry_are_separate": True,
            "no_automatic_TI1_claim_or_TS_verification_claim": True,
        },
    )
    write_json(folder / "starting_assessment.json", classify(model, xyz))
    result = run_stage(
        case_path, source, folder / "02_acyl/01_seed_relaxation", terms, max_steps=1500
    )
    if not result["forceconv"]:
        return finish(folder, "stopped_seed_force_not_converged", result)
    release = run_stage(
        case_path,
        Path(result["endpath"]),
        folder / "02_acyl/02_unbiased_release",
        [],
        max_steps=2500,
    )
    endpoint = np.asarray(Rst7.open(release["endpath"]).coordinates).reshape(-1, 3)
    assessment = classify(model, endpoint)
    write_json(folder / "02_acyl/endpoint_assessment.json", assessment)
    if not release["forceconv"] or not assessment["product_geometry_retained"]:
        return finish(folder, "stopped_unbiased_product_not_accepted", release)
    cv = asynchronous_coordinate(model)

    def value(points: np.ndarray) -> float:
        return evaluate_bias(points, [HarmonicTarget(cv, 0, 0)]).values[0]

    start, end = value(endpoint), value(model.coordinates)
    path = folder / "03_reverse_path"
    path.mkdir()
    write_json(
        path / "path_protocol.json",
        {
            "direction": "acylenzyme_to_ES",
            "coordinate": HarmonicTarget(cv, 0, 0).to_metadata()["cv"],
            "start_target_A": start,
            "end_target_A": end,
            "target_values_A": np.linspace(start, end, 37).tolist(),
            "physical_reference": "First accepted acyl-enzyme endpoint; not an E-S-relative barrier",
        },
    )
    latest = Path(release["endpath"])
    for index, target in enumerate(np.linspace(start, end, 37)[1:], start=1):
        result = run_stage(
            case_path,
            latest,
            path / f"window_{index:03d}",
            [HarmonicTarget(cv, float(target), 100)],
            max_steps=1500,
        )
        if not result["forceconv"]:
            return finish(
                folder, f"stopped_path_window_{index:03d}_not_converged", result
            )
        latest = Path(result["endpath"])
    final = run_stage(
        case_path, latest, path / "unbiased_ES_release", [], max_steps=2500
    )
    geometry = np.asarray(Rst7.open(final["endpath"]).coordinates).reshape(-1, 3)
    audit = classify(model, geometry)
    write_json(path / "endpoint_assessment.json", audit)
    status = (
        "path_complete_requires_TS_verification"
        if final["forceconv"] and audit["reactant_geometry"]
        else "path_not_product_to_reactant_connected"
    )
    finish(folder, status, final)


def finish(folder: Path, status: str, result: dict[str, Any]) -> None:
    """Record bounded termination without conflating process and scientific pass."""
    write_json(
        folder / "campaign_status.json",
        {"finished": True, "status": status, "last_result": result},
    )
    LOGGER.info("Campaign finished: %s", status)


def main() -> None:
    """Run the declared persistent campaign and retain uncaught failure evidence."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    try:
        run(args.case)
    except BaseException as error:
        path = args.case.resolve().parent / "campaign_status.json"
        if not path.exists():
            write_json(
                path,
                {"finished": True, "status": "process_failure", "error": str(error)},
            )
        raise


if __name__ == "__main__":
    main()
