"""Prepare three unchanged-inventory N-terminal relay release seeds."""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations, product
from pathlib import Path
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from netCDF4 import Dataset  # noqa: E402
from parmed.amber import AmberParm, Rst7  # noqa: E402
from parmed.amber.mask import AmberMask  # noqa: E402

from nyl12_covalent.relay_release import aligned_rmsd, geometry, minimization_input  # noqa: E402
from nyl12_covalent.stationary_native import NativeState, digest, write_json  # noqa: E402
from scripts.prepare_dimer_continuation import write_restart  # noqa: E402

SEEDS = (
    ("ti1_window032", 32, 4, 65),
    ("ti1_window033", 33, 5, 46),
    ("early_window028", 28, 4, 98),
)


def prepare(base: Path, output: Path) -> None:
    """Extract reviewed source frames without modifying proton inventories.

    Args:
        base: Existing umbrella_sampling directory.
        output: New single bounded campaign directory.
    """
    reference = (
        base.parent
        / "ts_image08_attempt1/validation_loop_20261121/verified_frame_validation_attempt05"
    )
    state = NativeState(reference)
    if digest(state.topology) != digest(base / "inputs/system.prmtop"):
        raise ValueError("US/native topology mismatch")
    parm = AmberParm(str(state.topology))
    source_qm = re.findall(r"qmmask='([^']+)'", state.mdin_source)[0]
    qm = np.asarray(AmberMask(parm, source_qm).Selection(), dtype=bool)
    if sorted({a.atomic_number for a, yes in zip(parm.atoms, qm) if yes}) != [
        1,
        6,
        7,
        8,
    ]:
        raise ValueError("Unexpected QM element inventory")
    output.mkdir(parents=True, exist_ok=False)
    snapshot = output / "source_snapshot"
    snapshot.mkdir()
    sources = [
        ROOT / "scripts/prepare_neb_us_relay_releases.py",
        ROOT / "scripts/run_neb_us_relay_release.py",
        ROOT / "scripts/prepare_dimer_continuation.py",
    ]
    sources += [
        ROOT / "src/nyl12_covalent" / name
        for name in (
            "relay_release.py",
            "relay_release_run.py",
            "stationary_native.py",
            "native_sander_cli.py",
            "dftb_scc_validation.py",
        )
    ]
    for source in sources:
        target = snapshot / source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    parameters = output / "parameters_3ob"
    parameters.mkdir()
    for a, b in product(("H", "C", "N", "O"), repeat=2):
        source = (
            Path("/data/jpc/code/ambertools26/venv/dat/slko/3ob-3-1") / f"{a}-{b}.skf"
        )
        shutil.copy2(source, parameters / source.name)
    manifests, coordinates, local_mm = [], {}, {}
    links = sum(bool(qm[b.atom1.idx]) != bool(qm[b.atom2.idx]) for b in parm.bonds)
    qm_electrons = (
        sum(a.atomic_number for a, yes in zip(parm.atoms, qm) if yes) + links + 2
    )
    if qm_electrons % 2:
        raise ValueError("Odd QM electron count")
    for name, window, block, frame in SEEDS:
        source = base / f"window_{window:03d}/03_production_block_{block:02d}_5ps"
        frozen = (
            base
            / "analysis/pmf_refresh_20261122/frozen_inputs"
            / source.relative_to(base)
        )
        with Dataset(source / "trajectory.nc") as ds:
            xyz = np.asarray(ds["coordinates"][frame - 1], dtype=float)
            time_ps = float(ds["time"][frame - 1])
            box = np.r_[ds["cell_lengths"][frame - 1], ds["cell_angles"][frame - 1]]
        if len(xyz) != len(state.xyz):
            raise ValueError("Source frame atom count mismatch")
        audit = geometry(xyz)
        if audit["proton_bonds"] != {
            "H1": ["N"],
            "H2": ["N"],
            "HG": ["N"],
            "HN": ["NL"],
        }:
            raise ValueError("Reviewed proton state differs from source coordinates")
        observable = json.loads((source / "observables/observables.json").read_text())
        trajectory_hash = digest(source / "trajectory.nc")
        if (
            observable["source_sha256"][str(source / "trajectory.nc")]
            != trajectory_hash
        ):
            raise ValueError("Trajectory differs from accepted observable source")
        if digest(frozen / "observables/observables.csv") != digest(
            source / "observables/observables.csv"
        ):
            raise ValueError("Source observables differ from frozen accepted release")
        seed = output / name
        inputs = seed / "inputs"
        inputs.mkdir(parents=True)
        for filename in ("system.prmtop", "dimer.mdin", "scientific_acceptance.json"):
            shutil.copy2(reference / "inputs" / filename, inputs / filename)
        ref = Rst7.open(str(state.restart_path))
        ref.box = box
        write_restart(inputs / "optimized.ncrst", xyz, ref)
        for filename in ("acceptance.json", "run.mdin", "input_identity.json"):
            shutil.copy2(source / filename, inputs / f"source_{filename}")
        shutil.copy2(
            frozen / "observables/observables.csv", inputs / "source_observables.csv"
        )
        (seed / "minimize.mdin").write_text(minimization_input(state.mdin_source))
        coordinates[name] = xyz[qm]
        local_mm[name] = xyz[state.active & ~qm]
        identity = {
            "seed": name,
            "source_trajectory": str(source / "trajectory.nc"),
            "source_trajectory_sha256": trajectory_hash,
            "window": window,
            "block": block,
            "frame_1based": frame,
            "time_ps": time_ps,
            "full_atom_count": len(xyz),
            "topology_sha256": digest(inputs / "system.prmtop"),
            "coordinate_sha256": digest(inputs / "optimized.ncrst"),
            "active_atom_count": int(state.active.sum()),
            "active_atom_ids_sha256": hashlib.sha256(
                state.indices.astype("<i8").tobytes()
            ).hexdigest(),
            "qm_atom_count": int(qm.sum()),
            "qm_charge": -2,
            "topology_charge_e": float(sum(a.charge for a in parm.atoms)),
            "qm_mm_link_count": links,
            "qm_electrons_with_links": qm_electrons,
            "qm_electron_parity_even": qm_electrons % 2 == 0,
            "qm_atom_ids_sha256": hashlib.sha256(
                np.flatnonzero(qm).astype("<i8").tobytes()
            ).hexdigest(),
            "initial_geometry": audit,
        }
        write_json(seed / "source_identity.json", identity)
        manifests.append(identity)
    rmsd = {
        f"{a}__{b}": aligned_rmsd(coordinates[a], coordinates[b])
        for a, b in combinations(coordinates, 2)
    }
    mm_rmsd = {
        f"{a}__{b}": aligned_rmsd(local_mm[a], local_mm[b])
        for a, b in combinations(local_mm, 2)
    }
    if any(value < 0.05 and mm_rmsd[key] < 0.1 for key, value in rmsd.items()):
        raise ValueError("Insufficient geometric distinction")
    write_json(
        output / "prepared.json",
        {
            "seeds": manifests,
            "qm_aligned_rmsd_A": rmsd,
            "local_mm_aligned_rmsd_A": mm_rmsd,
            "hypothesis": "N-terminal-protonated sampled structures may retain a TI1 basin after bias removal",
            "boundary": "Unrestrained within unchanged 4165-atom 12A movable list; all other atoms frozen; no reaction-coordinate/position potentials",
            "acceptance": "Native termination, SCC, max force<=1e-5 Eh/bohr, energy decrease, fixed shell and proton geometry recorded independently; no Hessian/minimum claim",
            "source_sha256": {
                str(p.relative_to(snapshot)): digest(p)
                for p in snapshot.rglob("*")
                if p.is_file()
            },
            "parameter_sha256": {p.name: digest(p) for p in parameters.iterdir()},
            "selection_scope": "Locally preserved accepted payloads windows028-034; 17 metadata-only blocks excluded explicitly, not a population estimate",
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.base.resolve(), args.output.resolve())
