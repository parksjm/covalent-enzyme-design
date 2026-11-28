"""Prepare exactly six reviewed boundary histories without launching MD."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask
from scipy.io import netcdf_file

from nyl12_covalent.neb_umbrella_campaign import digest
from nyl12_covalent.neb_umbrella_observables import (
    frame_observables,
    topology_proton_pool,
    _distances,
)
from nyl12_covalent.relay_activation_support import prepare_reviewed_targets
from nyl12_covalent.relay_targeted import verify_target
from nyl12_covalent.relay_umbrella_diagnostics import ROLES

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K"
)
WORK = Path(__file__).resolve().parent
OUTPUT = BASE / "activation_support_boundary_20261128"


def validate_six(targets: list[dict[str, Any]]) -> None:
    """Reject any departure from three approved centers and two distinct starts."""
    expected = {
        (x, start)
        for x in (-1.0, -0.8, -0.6)
        for start in ("original_es", "dominant_es")
    }
    actual = {
        (float(row["centers_A"][1]), row["source_lineage"]["start_kind"])
        for row in targets
    }
    if len(targets) != 6 or actual != expected:
        raise ValueError("Exactly three paired ES bias locations required")
    for row in targets:
        if (
            row["centers_A"][0] != -2.0
            or row["feature_names"] != ["q_acyl", "q_act"]
            or row["k"] != 25
        ):
            raise ValueError("Unapproved bias or hidden q_del restraint")


def main() -> None:
    """Freeze selected-frame provenance, audited restarts and immutable inputs."""
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    parent = json.loads((BASE / "prepared.json").read_text())
    support = json.loads(
        (BASE / "activation_support_20261127/prepared.json").read_text()
    )
    old_hashes = {
        str(f): digest(f)
        for f in [
            BASE / "prepared.json",
            BASE / "activation_support_20261127/prepared.json",
            BASE
            / "analysis/main100_continuation_20261127/current_ext100/checkpoint_032_main70ps_supplement43/sources.json",
        ]
    }
    topology = BASE / "inputs/system.prmtop"
    top = parmed.load_file(str(topology))
    qm = [
        i + 1 for i in AmberMask(top, parent["system"]["qm_atom_ids_mask"]).Selected()
    ]
    active = np.array(
        list(AmberMask(top, parent["system"]["active_mask"]).Selected()), dtype=int
    )
    assert (
        len(top.atoms) == 104922
        and len(qm) == 98
        and len(active) == 4165
        and parent["system"]["qm_charge"] == -2
    )
    pool, initial_owners, owners = topology_proton_pool(topology, qm)
    stage = (
        BASE / "activation_support_20261127/activation_support_3_03/04_production_5ps"
    )
    trajectory = stage / "trajectory.nc"
    stat = trajectory.stat()
    with netcdf_file(str(trajectory), "r", mmap=True) as stream:
        xyz = np.array(stream.variables["coordinates"][51], dtype=float, copy=True)
        box = np.r_[
            np.array(stream.variables["cell_lengths"][51], copy=True),
            np.array(stream.variables["cell_angles"][51], copy=True),
        ]
        time = float(stream.variables["time"][51])
    assert xyz.shape == (104922, 3) and np.isfinite(xyz).all()
    restart = Rst7(
        natom=len(xyz),
        title="Exact native frame52; fresh velocities required",
        time=time,
    )
    restart.coordinates = xyz
    restart.box = box
    rare = WORK / "dominant_frame52.rst7"
    restart.write(str(rare))
    rare_read = Rst7.open(str(rare))
    assert not rare_read.hasvels
    np.testing.assert_allclose(
        rare_read.coordinates.reshape(-1, 3), xyz, atol=5.1e-8, rtol=0
    )
    ordinary_row = next(
        row for row in support["windows"] if row["label"] == "activation_support_3_01"
    )
    ordinary = Path(ordinary_row["source_restart"])
    assert digest(ordinary) == ordinary_row["source_restart_sha256"]
    source_receipts = [
        Path(ordinary_row["source_lineage"][key]) for key in ["acceptance", "process"]
    ]
    assert json.loads(source_receipts[0].read_text())["accepted"]
    assert json.loads(source_receipts[1].read_text())["exit_status"] == 0
    frozen = np.ones(len(xyz), dtype=bool)
    frozen[active] = False
    qc = {}
    for name, source in [("original_es", ordinary), ("dominant_es", rare)]:
        value = Rst7.open(str(source))
        coordinates = np.asarray(value.coordinates).reshape(-1, 3)
        obs = frame_observables(
            coordinates, np.asarray(value.box), ROLES, pool, owners, initial_owners
        )
        obs["r_nalpha_hg1"] = float(
            _distances(coordinates, [13905], [13917], value.box)[0, 0]
        )
        assert [
            obs[key] for key in ["hg1_owner", "nalpha_h1_owner", "nalpha_h2_owner"]
        ] == [13916, 13905, 13905]
        assert obs["r_c_n"] < 1.6 and obs["r_c_og1"] > 2.5
        qc[name] = {"source": str(source), "sha256": digest(source), "observables": obs}
        if name == "dominant_es":
            distances = np.array(
                [
                    obs[key]
                    for key in [
                        "r_og1_hg1",
                        "r_nleave_hg1",
                        "r_nalpha_hg1",
                        "r_c_og1",
                        "r_c_n",
                    ]
                ]
            )
            np.testing.assert_allclose(
                distances, [0.9521, 4.1486, 1.7658, 3.3707, 1.3690], atol=6e-5, rtol=0
            )
        else:
            deviation = float(np.max(np.abs(coordinates[frozen] - xyz[frozen])))
            assert deviation < 1e-3
            qc["frozen_shell_between_starts_max_abs_A"] = deviation
    assert (
        trajectory.stat().st_size == stat.st_size
        and trajectory.stat().st_mtime_ns == stat.st_mtime_ns
    )
    qc["selected_native_frame"] = {
        "trajectory": str(trajectory),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "frame_one_based": 52,
        "time_ps": time,
        "coordinates_float64_sha256": hashlib.sha256(xyz.tobytes()).hexdigest(),
        "cell_sha256": hashlib.sha256(box.tobytes()).hexdigest(),
        "whole_trajectory_hash": "not scanned; selected-coordinate payload hash plus native accepted-stage receipts",
        "velocity_policy": "Fresh300K velocities in first1ps stage; no trajectory velocities inherited.",
    }
    qc["unchanged_manifests_sha256"] = old_hashes
    (WORK / "seed_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
    plan = dict(support["targets_plan"])
    plan.update(parent_prepared=str(BASE / "prepared.json"), output=str(OUTPUT))
    manifest = {
        str(Path(__file__).resolve()): digest(Path(__file__).resolve()),
        str(WORK / "seed_qc.json"): digest(WORK / "seed_qc.json"),
    }
    for f in [
        topology,
        ordinary,
        rare,
        *source_receipts,
        stage / "acceptance.json",
        stage / "process.json",
        stage / "raw_features.csv",
    ]:
        manifest[str(f)] = digest(f)
    for f in (ROOT / "src/nyl12_covalent").glob("*.py"):
        manifest[str(f)] = digest(f)
    plan["source_manifest"] = manifest
    used_seeds = {spec["seed"] for row in support["windows"] for spec in row["stages"]}
    used_seeds.update(
        spec["seed"] for row in parent["windows"] for spec in row["stages"]
    )
    targets = []
    for i, center in enumerate([-1.0, -0.8, -0.6]):
        for j, (kind, path) in enumerate(
            [("original_es", ordinary), ("dominant_es", rare)]
        ):
            seed = 2030120000 + 20 * (2 * i + j)
            assert not used_seeds.intersection(range(seed, seed + 8))
            targets.append(
                {
                    "label": f"activation_boundary_{i + 1:02d}_{kind}",
                    "source_restart": str(path),
                    "source_restart_sha256": digest(path),
                    "initial_ig": seed,
                    "feature_names": ["q_acyl", "q_act"],
                    "centers_A": [-2.0, center],
                    "k": 25,
                    "source_lineage": {
                        "start_kind": kind,
                        "original_lineage": ordinary_row["source_lineage"]
                        if j == 0
                        else qc["selected_native_frame"],
                    },
                    "validation": {
                        "proton_owners_expected": [13916, 13905, 13905],
                        "seed_qc_sha256": digest(WORK / "seed_qc.json"),
                    },
                }
            )
    plan["targets"] = targets
    plan["validation"] = {
        "status": "local_preparation_only_no_launch",
        "physical_hamiltonian_unchanged": True,
        "source_qm_charge": -2,
        "source_system": parent["system"],
        "require_fresh_resource_approval": True,
    }
    path = WORK / "boundary_targets_plan.json"
    path.write_text(json.dumps(plan, indent=2) + "\n")
    prepared = prepare_reviewed_targets(path, OUTPUT.name, validate_six)
    cfg = verify_target(prepared, digest(prepared))
    assert len(cfg["windows"]) == 6 and all(
        len(row["stages"]) == 8 for row in cfg["windows"]
    )
    for filename, expected in old_hashes.items():
        assert digest(Path(filename)) == expected
    receipt = {
        "status": "prepared_not_launched",
        "prepared": str(prepared),
        "sha256": digest(prepared),
        "six_histories": 6,
        "stages_per_history": 8,
        "equilibration_ps": 5,
        "production_ps": 30,
        "old_manifests_unchanged": old_hashes,
        "seeds_unique_against_original57_and33": True,
        "release_gate": "PI Geranium shared-host overnight reservation plus root final QC; no runtime authorized by this preparation",
    }
    (WORK / "preparation_acceptance.json").write_text(
        json.dumps(receipt, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
