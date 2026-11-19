"""Coordinate-only broad-cap seed diagnostic; no optimization or acceptance."""

from __future__ import annotations

import logging
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.acylenzyme_seed import AcylAtomMap  # noqa: E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa: E402
from nyl12_covalent.relay_acyl_seed import RelayAcylAtomMap, build_balanced_acyl_seed  # noqa: E402
from nyl12_covalent.relay_case import load_case, write_json  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Retain source snapshots and report actual broad-cap movement/contact data."""
    folder = Path(__file__).parent
    case = folder.parent / "case.json"
    model = load_case(case)
    cfg, m = model.config, model.atom_map
    snapshot = folder / "construction_diagnostic_sources"
    snapshot.mkdir()
    paths = [
        case,
        Path(cfg["topology"]),
        Path(cfg["source_restart"]),
        Path(__file__),
        ROOT / "src/nyl12_covalent/relay_acyl_seed.py",
        ROOT / "src/nyl12_covalent/acylenzyme_seed.py",
        ROOT / "src/nyl12_covalent/relay_case.py",
    ]
    hashes = {}
    for index, path in enumerate(paths):
        digest = source_sha256(path)
        copied = snapshot / f"{index:02d}_{path.name}"
        shutil.copy2(path, copied)
        if source_sha256(copied) != digest:
            raise ValueError("Immutable diagnostic source snapshot mismatch")
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
    product, evidence = build_balanced_acyl_seed(
        model.topology,
        model.coordinates,
        mapping,
        active_atom_ids=tuple(cfg["active_atom_ids"]),
        qm_charge=cfg["qm_charge"],
        max_heavy_displacement_A=20.0,
        max_hydrogen_displacement_A=20.0,
    )
    lengths = np.linalg.norm(product - model.coordinates, axis=1)

    def atom_record(index: int) -> dict:
        """Return source identity and measured displacement of one atom."""
        atom = model.topology.atoms[index]
        return dict(
            atom_id_1based=index + 1,
            atom_name=atom.name,
            atomic_number=atom.atomic_number,
            residue_name=atom.residue.name,
            residue_id_1based=atom.residue.idx + 1,
            displacement_A=float(lengths[index]),
        )

    heavy = [
        i for i, atom in enumerate(model.topology.atoms) if atom.atomic_number != 1
    ]
    hydrogens = [
        i for i, atom in enumerate(model.topology.atoms) if atom.atomic_number == 1
    ]
    contacts = []
    for contact in evidence["close_contacts"]:
        first, second = (i - 1 for i in contact["atoms_1based"])
        initial = float(
            np.linalg.norm(model.coordinates[first] - model.coordinates[second])
        )
        contacts.append(
            dict(
                **contact,
                atoms=[atom_record(first), atom_record(second)],
                source_distance_A=initial,
                distance_change_A=contact["distance_A"] - initial,
            )
        )
    contacts.sort(key=lambda row: row["distance_A"])
    report = dict(
        status="diagnostic coordinate hypothesis only; no engine/acceptance",
        caps_A=dict(heavy=20.0, hydrogen=20.0),
        source_sha256=hashes,
        max_heavy=atom_record(max(heavy, key=lambda i: lengths[i])),
        max_hydrogen=atom_record(max(hydrogens, key=lambda i: lengths[i])),
        largest_movements=[atom_record(i) for i in np.argsort(lengths)[-25:][::-1]],
        contact_count=len(contacts),
        contacts=contacts,
        construction=evidence,
    )
    write_json(folder / "construction_diagnostic.json", report)
    write_json(
        folder / "construction_diagnostic_scientific_acceptance.json",
        dict(
            accepted=False,
            optimized=False,
            production_cap_change_authorized=False,
            reason="Broad-cap coordinate diagnostic only; root scientific review required",
        ),
    )
    LOGGER.info(
        "Diagnostic heavy max=%s H max=%s contacts=%d",
        report["max_heavy"],
        report["max_hydrogen"],
        len(contacts),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
