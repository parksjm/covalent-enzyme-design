"""One bounded acyl-arm torsion sweep; preserve best leaving placement."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import shutil
import sys

import numpy as np
from parmed.amber import Rst7

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from nyl12_covalent.acylenzyme_seed import AcylAtomMap, _component, _unit  # noqa: E402
from nyl12_covalent.native_sander_cli import source_sha256  # noqa: E402
from nyl12_covalent.relay_acyl_seed import (  # noqa: E402
    RelayAcylAtomMap,
    _source_graph,
    build_balanced_acyl_seed,
)  # noqa: E402
from nyl12_covalent.relay_case import load_case, write_json  # noqa: E402
from nyl12_covalent.relay_seed_search import _Contacts, _rotation  # noqa: E402
from scripts.prepare_dimer_continuation import write_restart  # noqa: E402


def main() -> None:
    """Retain all 24 scores and publish a seed only if severe contacts vanish."""
    folder = Path(__file__).parent
    out = folder / "acyl_torsion_search"
    out.mkdir(exist_ok=True)
    model = load_case(folder.parent / "case.json")
    m, cfg = model.atom_map, model.config
    previous = json.loads(
        (folder / "coordinate_seed_search/search_summary.json").read_text()
    )
    start = folder / "coordinate_seed_search/best_diagnostic.ncrst"
    xyz = np.asarray(Rst7.open(str(start)).coordinates).reshape(-1, 3)
    mapping = RelayAcylAtomMap(
        AcylAtomMap(**{k: m[k] for k in AcylAtomMap.__dataclass_fields__}),
        m["water_o"],
        m["water_h1"],
        m["water_h2"],
        m["ash_oxygen"],
        m["ash_proton"],
        (m["asp_o1"], m["asp_o2"]),
    )
    _, construction = build_balanced_acyl_seed(
        model.topology,
        model.coordinates,
        mapping,
        active_atom_ids=tuple(cfg["active_atom_ids"]),
        qm_charge=cfg["qm_charge"],
        max_heavy_displacement_A=20.0,
        max_hydrogen_displacement_A=20.0,
        idealize_ester_angle=False,
    )
    graph, ids, _ = _source_graph(model.topology, model.coordinates, mapping)
    c, o, n = ids["carbon"], ids["oxygen"], ids["leaving_n"]
    anchor = next(i for i in graph[c] if i not in (o, n))
    graph[c].remove(anchor)
    graph[anchor].remove(c)
    arm = np.array(sorted(_component(graph, anchor)), dtype=int)
    if c in arm or o in arm or n in arm:
        raise ValueError("Acyl torsion cut failed to isolate acyl arm")
    graph[c].add(anchor)
    graph[anchor].add(c)
    for operation, pairs in construction["bond_edits_1based"].items():
        for a, b in pairs:
            a, b = a - 1, b - 1
            if operation == "remove":
                graph[a].remove(b)
                graph[b].remove(a)
            else:
                graph[a].add(b)
                graph[b].add(a)
    moving = np.array(
        sorted(
            set(
                np.flatnonzero(
                    np.linalg.norm(xyz - model.coordinates, axis=1) > 1e-10
                ).tolist()
            )
            | set(arm.tolist())
        )
    )
    if not set((moving + 1).tolist()) <= set(cfg["active_atom_ids"]):
        raise ValueError("Torsion arm outside original active region")
    scorer = _Contacts(model.topology, model.coordinates, moving, graph)
    axis = _unit(xyz[anchor] - xyz[c])
    best_score, best, choice, clashes = None, None, None, None
    rows = []
    for angle in range(-180, 180, 15):
        candidate = xyz.copy()
        candidate[arm] = (xyz[arm] - xyz[c]) @ _rotation(axis, angle).T + xyz[c]
        score, contacts = scorer.assess(candidate)
        rows.append(dict(torsion_degrees=angle, score=list(score)))
        if best_score is None or score < best_score:
            best_score, best, choice, clashes = score, candidate, angle, contacts
    snapshots = out / "source_snapshot"
    snapshots.mkdir()
    hashes = {}
    paths = [
        Path(__file__),
        folder.parent / "case.json",
        start,
        Path(cfg["topology"]),
        Path(cfg["source_restart"]),
    ]
    paths += [
        ROOT / "src/nyl12_covalent" / name
        for name in (
            "relay_seed_search.py",
            "relay_acyl_seed.py",
            "acylenzyme_seed.py",
            "relay_case.py",
        )
    ]
    for index, path in enumerate(paths):
        digest = source_sha256(path)
        target = snapshots / f"{index:02d}_{path.name}"
        shutil.copy2(path, target)
        if source_sha256(target) != digest:
            raise ValueError("Torsion input snapshot mismatch")
        hashes[str(path)] = digest
    displacement = np.linalg.norm(best - model.coordinates, axis=1)
    heavy = np.array([a.atomic_number != 1 for a in model.topology.atoms])
    summary = dict(
        status="coordinate hypothesis only; no engine or acceptance",
        candidate_count=24,
        torsion_axis_1based=[c + 1, anchor + 1],
        torsion_arm_1based=(arm + 1).tolist(),
        best_torsion_degrees=choice,
        severe_contact_count=int(best_score[0]),
        close_contacts=clashes,
        candidates=rows,
        previous_leaving_choice=previous["best_choice"],
        source_sha256=hashes,
        max_heavy_displacement_A=float(displacement[heavy].max(initial=0)),
        max_hydrogen_displacement_A=float(displacement[~heavy].max(initial=0)),
    )
    write_json(out / "search_summary.json", summary)
    write_restart(out / "best_diagnostic.ncrst", best, model.restart)
    if best_score[0] == 0:
        for path in (folder / "construction.json", folder / "acyl_guess.ncrst"):
            if path.exists():
                raise FileExistsError(path)
        construction.update(
            search=summary,
            source_sha256=hashes,
            close_contacts=clashes,
            moved_indices_1based=(np.flatnonzero(displacement > 1e-10) + 1).tolist(),
            max_heavy_displacement_A=summary["max_heavy_displacement_A"],
            max_hydrogen_displacement_A=summary["max_hydrogen_displacement_A"],
        )
        construction["distances_A"] = {
            key: float(np.linalg.norm(best[a] - best[b]))
            for key, a, b in (
                ("OG1_C", ids["thr_og"], c),
                ("C_O", c, o),
                ("C_Nleaving", c, n),
            )
        }
        write_json(folder / "construction.json", construction)
        write_restart(folder / "acyl_guess.ncrst", best, model.restart)
    write_json(
        out / "scientific_acceptance.json",
        dict(
            accepted=False,
            optimized=False,
            severe_coordinate_gate_passed=best_score[0] == 0,
        ),
    )
    logging.info(
        "Torsion best=%s severe=%d contacts=%d heavy_move=%.4f H_move=%.4f",
        choice,
        best_score[0],
        len(clashes),
        summary["max_heavy_displacement_A"],
        summary["max_hydrogen_displacement_A"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
