"""Require the explicit water and its two atom restraints to survive AF3 features."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components
from alphafold3.data import featurisation
from alphafold3.model import features


def main() -> None:
    """Fail before inference if water or either restraint is missing."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.output.exists():
        raise FileExistsError(args.output)
    inputs = list(folding_input.load_fold_inputs_from_path(args.input))
    if len(inputs) != 1:
        raise ValueError("Expected exactly one fold input")
    fold = inputs[0]
    ccd = chemical_components.Ccd(user_ccd=fold.user_ccd)
    structure = fold.to_structure(ccd)
    water = (structure.chain_id == "J") & (structure.atom_name == "O")
    if np.count_nonzero(water) != 1:
        raise ValueError("Input structure lost the single water oxygen J1:O")
    batches = featurisation.featurise_input(fold, ccd, buckets=(512,))
    if len(batches) != 1:
        raise ValueError("Expected exactly one featurized seed")
    constraints = features.ContactConstraints.from_data_dict(batches[0])
    expected_count = len(fold.contact_constraints)
    if expected_count < 2 or int(np.count_nonzero(constraints.mask)) != expected_count:
        raise ValueError(
            "Every declared water contact restraint must map to retained atoms"
        )
    report = {
        "accepted": True,
        "input_water_oxygen_count": int(np.count_nonzero(water)),
        "expected_contact_count": expected_count,
        "mapped_contact_count": int(np.count_nonzero(constraints.mask)),
        "token_indices": np.asarray(constraints.token_indices)[
            np.asarray(constraints.mask, dtype=bool)
        ].tolist(),
        "min_distances_A": np.asarray(constraints.min_distances)[:2].tolist(),
        "max_distances_A": np.asarray(constraints.max_distances)[:2].tolist(),
        "limitations": "Input/features validation only; no hydrogen orientation or catalytic geometry acceptance",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    logging.info("Water and all declared contacts retained in AF3 features")


if __name__ == "__main__":
    main()
