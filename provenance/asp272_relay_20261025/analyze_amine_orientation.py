"""Audit terminal amine hydrogen bonds and donor alignment without moving atoms."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
from parmed.amber import Rst7

ROOT = Path("/data/jpc/covalent-enzyme-design")
RUN = ROOT / "output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2"
ATOMS = {
    "N": 13905,
    "H1": 13906,
    "H2": 13907,
    "CA": 13908,
    "OG1": 13916,
    "HG1": 13917,
    "Asp272_OD1": 14559,
    "Asp272_OD2": 14560,
    "Asp274_OD1": 14578,
    "Asp274_OD2": 14579,
}


def angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Return the a-b-c angle in degrees."""
    x, y = a - b, c - b
    return float(
        np.degrees(
            np.arccos(
                np.clip(np.dot(x, y) / np.linalg.norm(x) / np.linalg.norm(y), -1, 1)
            )
        )
    )


def main() -> None:
    """Write exact-coordinate contacts and a declared geometric lone-pair proxy."""
    records = {}
    for window in ("stage1_00", "stage1_13"):
        xyz = Rst7.open(str(RUN / window / "optimized.ncrst")).coordinates.reshape(
            (-1, 3)
        )
        a = {name: xyz[index - 1] for name, index in ATOMS.items()}
        contacts = []
        for hydrogen in ("H1", "H2"):
            for oxygen in ("Asp272_OD1", "Asp272_OD2", "Asp274_OD1", "Asp274_OD2"):
                contacts.append(
                    {
                        "donor": "N",
                        "hydrogen": hydrogen,
                        "acceptor": oxygen,
                        "H_O_A": float(np.linalg.norm(a[hydrogen] - a[oxygen])),
                        "N_H_O_degrees": angle(a["N"], a[hydrogen], a[oxygen]),
                    }
                )
        vectors = np.array([a[name] - a["N"] for name in ("CA", "H1", "H2")])
        proxy = -np.sum(vectors / np.linalg.norm(vectors, axis=1)[:, None], axis=0)
        target = a["HG1"] - a["N"]
        records[window] = {
            "contacts": contacts,
            "Asp272_OD1_HG1_N_angle_degrees": angle(a["Asp272_OD1"], a["HG1"], a["N"]),
            "Asp272_OD1_N_A": float(np.linalg.norm(a["Asp272_OD1"] - a["N"])),
            "lone_pair_proxy_to_incoming_HG1_degrees": float(
                np.degrees(
                    np.arccos(
                        np.clip(
                            np.dot(proxy, target)
                            / np.linalg.norm(proxy)
                            / np.linalg.norm(target),
                            -1,
                            1,
                        )
                    )
                )
            ),
            "proxy_definition": "Negative sum of unit N-CA/N-H1/N-H2 vectors; geometric diagnostic, not an electronic orbital calculation.",
        }
    (RUN / "amine_orientation.json").write_text(json.dumps(records, indent=2) + "\n")
    for name, record in records.items():
        logging.info(
            "%s donor angle %.1f deg, N...OD1 %.3f A, lone-pair proxy mismatch %.1f deg",
            name,
            record["Asp272_OD1_HG1_N_angle_degrees"],
            record["Asp272_OD1_N_A"],
            record["lone_pair_proxy_to_incoming_HG1_degrees"],
        )
        for contact in record["contacts"]:
            if contact["H_O_A"] < 2.6:
                logging.info(
                    "%s %s to %s: %.3f A, %.1f deg",
                    name,
                    contact["hydrogen"],
                    contact["acceptor"],
                    contact["H_O_A"],
                    contact["N_H_O_degrees"],
                )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
