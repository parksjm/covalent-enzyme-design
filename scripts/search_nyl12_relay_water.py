"""Search and save relay-water seeds for a hydrogenated Nyl12 TI1 model."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.relay_water import search_relay_water  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Run the deterministic water search and preserve accepted coordinates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--oxyanion", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    topology = AmberParm(str(args.topology), xyz=str(args.coordinates))
    coordinates = np.array(topology.coordinates, copy=True)
    report = search_relay_water(topology, coordinates, args.oxyanion)
    report["input_sha256"] = {
        str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (
            args.topology, args.coordinates, Path(__file__),
            ROOT / "src/nyl12_covalent/relay_water.py",
        )
    }
    (args.output / "water_search.json").write_text(json.dumps(report, indent=2) + "\n")
    accepted = [candidate for candidate in report["candidates"] if candidate["accepted"]]
    if not accepted:
        raise RuntimeError("No relay-water seed passed steric screening; see water_search.json")
    for index, candidate in enumerate(accepted, 1):
        moved = coordinates.copy()
        moved[report["water_atom_indices_zero_based"]] = candidate["water_xyz_A"]
        topology.coordinates = moved
        topology.save(str(args.output / f"candidate_{index:02d}.rst7"), format="rst7")
    LOGGER.info(
        "Saved %d distinct accepted seeds; %d oriented candidates passed; best shortfall %.3f A",
        len(accepted), report["accepted_oriented_candidate_count"],
        accepted[0]["maximum_contact_shortfall_A"],
    )


if __name__ == "__main__":
    main()
