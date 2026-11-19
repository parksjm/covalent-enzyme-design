"""Persist one bounded coordinate search; publish only a zero-severe-contact guess."""

from __future__ import annotations

import logging
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from nyl12_covalent.acylenzyme_seed import AcylAtomMap  # noqa: E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa: E402
from nyl12_covalent.relay_acyl_seed import RelayAcylAtomMap  # noqa: E402
from nyl12_covalent.relay_case import load_case, write_json  # noqa: E402
from nyl12_covalent.relay_seed_search import search_balanced_acyl_seed  # noqa: E402
from scripts.prepare_dimer_continuation import write_restart  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Retain immutable inputs and conditionally publish an unoptimized seed."""
    folder = Path(__file__).parent
    case = folder.parent / "case.json"
    destination = folder / "coordinate_seed_search"
    destination.mkdir(exist_ok=True)
    model = load_case(case)
    cfg, m = model.config, model.atom_map
    snapshots = destination / "source_snapshot"
    snapshots.mkdir()
    paths = [case, Path(cfg["source_restart"]), Path(cfg["topology"]), Path(__file__)]
    paths += [
        ROOT / "src/nyl12_covalent" / name
        for name in (
            "relay_seed_search.py",
            "relay_acyl_seed.py",
            "acylenzyme_seed.py",
            "relay_case.py",
        )
    ]
    hashes = {}
    for index, path in enumerate(paths):
        digest = source_sha256(path)
        copied = snapshots / f"{index:02d}_{path.name}"
        shutil.copy2(path, copied)
        if source_sha256(copied) != digest:
            raise ValueError("Source snapshot mismatch")
        hashes[str(path)] = digest
    acyl = AcylAtomMap(**{key: m[key] for key in AcylAtomMap.__dataclass_fields__})
    mapping = RelayAcylAtomMap(
        acyl,
        m["water_o"],
        m["water_h1"],
        m["water_h2"],
        m["ash_oxygen"],
        m["ash_proton"],
        (m["asp_o1"], m["asp_o2"]),
    )
    product, evidence = search_balanced_acyl_seed(
        model.topology,
        model.coordinates,
        mapping,
        active_atom_ids=tuple(cfg["active_atom_ids"]),
        qm_charge=cfg["qm_charge"],
    )
    summary = {key: value for key, value in evidence.items() if key != "construction"}
    summary["source_sha256"] = hashes
    summary["scientifically_accepted"] = False
    write_json(destination / "search_summary.json", summary)
    write_restart(destination / "best_diagnostic.ncrst", product, model.restart)
    if evidence["severe_contact_count"] == 0:
        for path in (folder / "construction.json", folder / "acyl_guess.ncrst"):
            if path.exists():
                raise FileExistsError(path)
        construction = dict(evidence["construction"])
        construction.update(
            source_sha256=hashes,
            search=summary,
            moved_indices_1based=evidence["moved_indices_1based"],
            close_contacts=evidence["close_contacts"],
            max_heavy_displacement_A=evidence["max_heavy_displacement_A"],
            max_hydrogen_displacement_A=evidence["max_hydrogen_displacement_A"],
        )
        ids = {key: value - 1 for key, value in m.items()}
        construction["distances_A"] = {
            key: float(np.linalg.norm(product[ids[a]] - product[ids[b]]))
            for key, a, b in (
                ("OG1_C", "thr_og", "carbon"),
                ("C_O", "carbon", "oxygen"),
                ("C_Nleaving", "carbon", "leaving_n"),
            )
        }
        write_json(folder / "construction.json", construction)
        write_restart(folder / "acyl_guess.ncrst", product, model.restart)
    write_json(
        destination / "scientific_acceptance.json",
        dict(
            accepted=False,
            optimized=False,
            coordinate_severe_contact_gate_passed=evidence["severe_contact_count"] == 0,
            reason="Coordinate hypothesis only; unbiased QM/MM endpoint acceptance remains pending",
        ),
    )
    LOGGER.info(
        "Best seed severe=%d contacts=%d heavy_move=%.3f H_move=%.3f choice=%s",
        evidence["severe_contact_count"],
        len(evidence["close_contacts"]),
        evidence["max_heavy_displacement_A"],
        evidence["max_hydrogen_displacement_A"],
        evidence["best_choice"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
