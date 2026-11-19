"""Verify real integrated preparation in a separate temporary campaign."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import sys

import numpy as np
from parmed.amber import Rst7

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts import prepare_relay_acyl_campaign as preparation  # noqa: E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa: E402
from nyl12_covalent.relay_case import write_json  # noqa: E402


def main() -> None:
    """Run preparation with only CAMPAIGN overridden and compare retained seeds."""
    output = Path(__file__).parent
    original = output.parent
    temporary = Path("/tmp/nyl12_acyl_recipe_reproduction_20261118")
    temporary.mkdir(exist_ok=False)
    preparation.CAMPAIGN = temporary
    preparation.prepare()
    fresh = temporary / "01_endpoint_inputs"
    old_restart, new_restart = (
        Rst7.open(str(folder / "acyl_guess.ncrst")) for folder in (original, fresh)
    )
    old_xyz = np.asarray(old_restart.coordinates).reshape(-1, 3)
    new_xyz = np.asarray(new_restart.coordinates).reshape(-1, 3)
    if old_xyz.shape != new_xyz.shape:
        raise ValueError("Reproduction atom count differs from running source seed")
    old = json.loads((original / "construction.json").read_text())
    new = json.loads((fresh / "construction.json").read_text())
    original_summary = old["search"]
    selected = dict(
        leaving=new["search"]["leaving"]["best_choice"],
        torsion=new["search"]["acyl_torsion"]["best_torsion_degrees"],
        axis=new["search"]["acyl_torsion"]["torsion_axis_1based"],
    )
    previous = dict(
        leaving=original_summary["previous_leaving_choice"],
        torsion=original_summary["best_torsion_degrees"],
        axis=original_summary["torsion_axis_1based"],
    )
    differences = np.linalg.norm(new_xyz - old_xyz, axis=1)
    manifest = json.loads((fresh / "manifest.json").read_text())
    report = dict(
        status="completed_no_engine_comparison",
        temporary_campaign=str(temporary),
        overridden_settings={"CAMPAIGN": str(temporary)},
        atom_count=len(old_xyz),
        original_seed_sha256=source_sha256(original / "acyl_guess.ncrst"),
        reproduced_seed_sha256=source_sha256(fresh / "acyl_guess.ncrst"),
        coordinate_arrays_exact=bool(np.array_equal(new_xyz, old_xyz)),
        periodic_box_exact=bool(np.array_equal(old_restart.box, new_restart.box)),
        max_atom_difference_A=float(differences.max(initial=0)),
        nonexact_atom_ids_1based=(np.flatnonzero(differences > 0) + 1).tolist(),
        coordinates_equal_within_1e_minus_12_A=bool(
            differences.max(initial=0) <= 1e-12
        ),
        candidate_choices_exact=selected == previous,
        selected_choices=selected,
        original_choices=previous,
        severe_contact_count=new["severe_contact_count"],
        original_source_sha256=source_sha256(original / "accepted_es.ncrst"),
        reproduced_source_sha256=source_sha256(fresh / "accepted_es.ncrst"),
        construction_sha256={
            "original": source_sha256(original / "construction.json"),
            "reproduced": source_sha256(fresh / "construction.json"),
        },
        preparation_source_sha256=manifest["input_sha256"],
        running_campaign_untouched=True,
        scientific_acceptance=False,
    )
    write_json(output / "comparison.json", report)
    logging.info(
        "Reproduction exact_coords=%s choice_exact=%s max_difference_A=%.12g severe=%d",
        report["coordinate_arrays_exact"],
        report["candidate_choices_exact"],
        report["max_atom_difference_A"],
        report["severe_contact_count"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
