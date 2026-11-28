"""Analyze one independent-state supplemental relay checkpoint, without publishing primary."""

import argparse
import json
from pathlib import Path

import numpy as np

from nyl12_covalent.relay_pooled_analysis import analyze
from nyl12_covalent.relay_pooled_sources import pooled_data


def main() -> None:
    """Use either validated cached61inputs or the complete94statecohort sources."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--horizon", type=int, required=True)
    parser.add_argument("--first61-inputs", type=Path)
    args = parser.parse_args()
    if args.first61_inputs:
        data = json.loads(args.first61_inputs.with_suffix(".json").read_text())
        arrays = np.load(args.first61_inputs.with_suffix(".npz"))
        data.update(q=arrays["q"], counts=arrays["counts"])
    else:
        data = pooled_data(args.cache, args.horizon, True)
        if data is None:
            raise ValueError("Exactcommon94statecheckpointnotready")
    analyze(data, args.output, workers=4)


if __name__ == "__main__":
    main()
