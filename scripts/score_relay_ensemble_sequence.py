#!/usr/bin/env python3
"""Inject a validated whole sequence into the unchanged native score entry point."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.ligandmpnn_pilot_audit import _pdb_sequences  # noqa: E402
from nyl12_covalent.relay_ensemble_sequences import audit_sequence  # noqa: E402


def main() -> None:
    """Score one external sequence while preserving native fixed atomic context."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    args = parser.parse_args()
    request = json.loads(args.request.read_text())
    engine = Path(request["score_engine"])
    sys.path.insert(0, str(engine.parent))
    import data_utils
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError("Authorized GPU runtime is unavailable; no CPU fallback")
    torch.set_num_threads(2)
    native = _pdb_sequences(Path(request["pdb"]))
    candidate = audit_sequence(request["sequence"], native, request["mask"])
    chains = dict(zip("ABCDEFGH", candidate["sequence"].split(":"), strict=True))
    parse_native = data_utils.parse_PDB

    def parse_candidate(*args: object, **kwargs: object) -> tuple:
        result = parse_native(*args, **kwargs)
        protein = result[0]
        original = protein["S"].clone()
        for i, (chain, position) in enumerate(
            zip(protein["chain_letters"], protein["R_idx"].tolist(), strict=True)
        ):
            protein["S"][i] = data_utils.alphabet.index(chains[chain][position - 1])
        if original.shape != protein["S"].shape:
            raise ValueError("Candidate sequence changed native shape")
        return result

    data_utils.parse_PDB = parse_candidate
    sys.argv = [str(engine), *request["argv"]]
    runpy.run_path(str(engine), run_name="__main__")


if __name__ == "__main__":
    main()
