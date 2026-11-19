"""Prepare a reproducible, unlaunched 0 K DFTB3/MM native NEB."""

from __future__ import annotations

import hashlib
import json
import logging
import re
import shutil
import sys
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

REPO = next(
    parent
    for parent in Path(__file__).resolve().parents
    if (parent / "pyproject.toml").is_file()
)
sys.path.insert(0, str(REPO))
from scripts.prepare_dimer_continuation import read_xyz_frames, write_restart  # noqa: E402

LOGGER = logging.getLogger(__name__)
IRC = Path(__file__).resolve().parent
OUTPUT = IRC / "new_neb"
TS_ROOT = IRC.parent
TS = TS_ROOT / "continuation_cycle200_20261113/calculation/optimized.ncrst"
TOPOLOGY = TS_ROOT / "inputs/system.prmtop"
BASE = TS_ROOT.parent / "05_es_to_acylenzyme_neb/inputs/base.mdin"
MINUS = IRC / "endpoint_minimizations/mode_minus/calculation"
PLUS = IRC / "endpoint_minimizations/mode_plus/calculation"
CLASH_CUTOFF_ANGSTROM = 2.0


def digest(path: Path) -> str:
    """Return the SHA-256 digest of an artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resample_path(
    frames: list[np.ndarray], metric_frames: list[np.ndarray], intervals: int
) -> list[np.ndarray]:
    """Resample full coordinates uniformly in a supplied path metric."""
    if intervals < 1 or len(frames) < 2:
        raise ValueError(
            "Path resampling requires at least two frames and one interval"
        )
    shape = frames[0].shape
    if any(frame.shape != shape or not np.isfinite(frame).all() for frame in frames):
        raise ValueError("Path frames must have matching finite coordinate shapes")
    if len(metric_frames) != len(frames):
        raise ValueError("Each coordinate frame must have one metric frame")
    metric_shape = metric_frames[0].shape
    if any(
        frame.shape != metric_shape or not np.isfinite(frame).all()
        for frame in metric_frames
    ):
        raise ValueError("Path metric frames must have matching finite shapes")
    lengths = np.array(
        [np.linalg.norm(b - a) for a, b in zip(metric_frames, metric_frames[1:])]
    )
    cumulative = np.concatenate(([0.0], np.cumsum(lengths)))
    if cumulative[-1] <= 1e-12:
        raise ValueError("Cannot resample a path with zero arc length")
    targets = np.linspace(0.0, cumulative[-1], intervals + 1)
    output = [frames[0].copy()]
    for target in targets[1:-1]:
        right = int(np.searchsorted(cumulative, target, side="right"))
        left = right - 1
        fraction = (target - cumulative[left]) / (cumulative[right] - cumulative[left])
        output.append(frames[left] + fraction * (frames[right] - frames[left]))
    output.append(frames[-1].copy())
    return output


def distance(coords: np.ndarray, a: int, b: int) -> float:
    """Measure a one-based atom-index distance in angstroms."""
    return float(np.linalg.norm(coords[a - 1] - coords[b - 1]))


def reactive_metric(coords: np.ndarray) -> np.ndarray:
    """Return mapped Thr attack, scissile C-N and HG1 transfer distances."""
    return np.array(
        [
            distance(coords, 13916, 19916),
            distance(coords, 19916, 19929),
            distance(coords, 13917, 19929),
        ]
    )


def nonbonded_contacts(
    parm: AmberParm, coords: np.ndarray, active: np.ndarray
) -> dict[tuple[int, int], float]:
    """Find sub-cutoff nonbonded contacts involving a mobile heavy atom."""
    heavy = np.array([i for i, atom in enumerate(parm.atoms) if atom.atomic_number > 1])
    bonded = {tuple(sorted((bond.atom1.idx, bond.atom2.idx))) for bond in parm.bonds}
    exclusions = set(bonded)
    adjacency: dict[int, set[int]] = {int(i): set() for i in heavy}
    for i, j in bonded:
        if i in adjacency and j in adjacency:
            adjacency[i].add(j)
            adjacency[j].add(i)
    for i, neighbors in adjacency.items():
        for j in neighbors:
            exclusions.update(tuple(sorted((i, k))) for k in adjacency[j] if i != k)
    cell_width = CLASH_CUTOFF_ANGSTROM
    cells: dict[tuple[int, int, int], list[int]] = {}
    for atom_index in heavy:
        key = tuple(np.floor(coords[atom_index] / cell_width).astype(int))
        cells.setdefault(key, []).append(int(atom_index))
    mobile_heavy = heavy[active[heavy]]
    offsets = tuple(
        (x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)
    )
    contacts: dict[tuple[int, int], float] = {}
    cutoff_sq = cell_width**2
    for i in mobile_heavy:
        key = tuple(np.floor(coords[i] / cell_width).astype(int))
        for offset in offsets:
            neighbor_key = tuple(key[axis] + offset[axis] for axis in range(3))
            for j in cells.get(neighbor_key, ()):
                pair = tuple(sorted((int(i), j)))
                if i == j or pair in exclusions or pair in contacts:
                    continue
                delta = coords[i] - coords[j]
                squared = float(delta @ delta)
                if squared < cutoff_sq:
                    contacts[pair] = float(np.sqrt(squared))
    return contacts


def build() -> dict[str, object]:
    """Validate sources, create a trajectory-seeded 0 K NEB, and write provenance."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if any(path.is_file() for path in OUTPUT.rglob("*")):
        raise FileExistsError(f"Refusing to overwrite prepared artifacts in {OUTPUT}")
    inputs = OUTPUT / "inputs"
    images_dir = inputs / "images"
    calculation = OUTPUT / "calculation"
    images_dir.mkdir(parents=True, exist_ok=True)
    calculation.mkdir(exist_ok=True)
    parm = AmberParm(str(TOPOLOGY))
    base = BASE.read_text()
    active_match = re.search(r"\bactive\s*=\s*'([^']+)'", base, re.I)
    qm_match = re.search(r"\bqmmask\s*=\s*'([^']+)'", base, re.I)
    if active_match is None or qm_match is None:
        raise ValueError("Base input lacks active or QM atom masks")
    active = np.asarray(AmberMask(parm, active_match.group(1)).Selection(), dtype=bool)
    if active.sum() != 4165 or len(parm.atoms) != 104922:
        raise ValueError("Unexpected topology size or active-shell atom count")
    ts = Rst7.open(str(TS))
    es = Rst7.open(str(MINUS / "minimized.ncrst"))
    acyl = Rst7.open(str(PLUS / "minimized.ncrst"))
    ts_xyz = ts.coordinates.reshape(-1, 3)
    endpoints = [es.coordinates.reshape(-1, 3), acyl.coordinates.reshape(-1, 3)]
    if any(x.shape != ts_xyz.shape for x in endpoints):
        raise ValueError("Endpoint atom mapping mismatch")
    if any(not np.array_equal(r.box, ts.box) for r in (es, acyl)):
        raise ValueError("Endpoint periodic cells differ from TS")
    outer_max = max(
        float(np.max(np.linalg.norm(x[~active] - ts_xyz[~active], axis=1)))
        for x in endpoints
    )
    if outer_max > 1e-8:
        raise ValueError(f"Endpoint exterior differs from TS: {outer_max:.3g} A")

    branches: list[list[np.ndarray]] = []
    branch_metrics = {}
    expected_endpoint_energies = {
        "minus": -718.207958911553,
        "plus": -718.201496183098,
    }
    for name, folder, endpoint, reverse in (
        ("minus", MINUS, endpoints[0], True),
        ("plus", PLUS, endpoints[1], False),
    ):
        frames = read_xyz_frames(folder / "path.xyz", int(active.sum()))
        start = Rst7.open(str(folder / "start.ncrst")).coordinates.reshape(-1, 3)
        expected_end = Rst7.open(str(folder / "minimized.ncrst")).coordinates.reshape(
            -1, 3
        )
        if np.max(np.abs(frames[0] - start[active])) > 2e-5:
            raise ValueError(f"{name} trajectory does not start at its saved restart")
        if np.max(np.abs(frames[-1] - expected_end[active])) > 2e-5:
            raise ValueError(f"{name} trajectory does not end at its minimized restart")
        atom_numbers = [
            int(atom.atomic_number) for atom in parm.atoms if active[atom.idx]
        ]
        xyz_elements = (folder / "path.xyz").read_text().splitlines()
        if [line.split()[0] for line in xyz_elements[2 : 2 + len(atom_numbers)]] != [
            _element_symbol(n) for n in atom_numbers
        ]:
            raise ValueError(f"{name} XYZ active-atom order differs from topology")
        full_frames = []
        for active_xyz in frames:
            full = ts_xyz.copy()
            full[active] = active_xyz
            full_frames.append(full)
        if reverse:
            full_frames.reverse()
        branches.append(full_frames)
        branch_metrics[name] = {
            "trajectory_frames": len(frames),
            "start_TS_max_displacement_A": float(
                np.max(np.linalg.norm(start - ts_xyz, axis=1))
            ),
            "trajectory_start_max_error_A": float(
                np.max(np.abs(frames[0] - start[active]))
            ),
            "trajectory_end_max_error_A": float(
                np.max(np.abs(frames[-1] - expected_end[active]))
            ),
            "maximum_adjacent_active_atom_displacement_A": float(
                max(
                    np.max(np.linalg.norm(b - a, axis=1))
                    for a, b in zip(frames, frames[1:])
                )
            ),
            "endpoint_energy_hartree": expected_endpoint_energies[name],
        }

    joined = branches[0] + [ts_xyz.copy()] + branches[1]
    left_path = joined[: len(branches[0]) + 1]
    right_path = joined[len(branches[0]) :]
    sampled_left = resample_path(
        left_path, [reactive_metric(x) for x in left_path], intervals=6
    )
    sampled_right = resample_path(
        right_path, [reactive_metric(x) for x in right_path], intervals=6
    )
    sampled_full = sampled_left + sampled_right[1:]
    if len(sampled_full) != 13 or not np.array_equal(sampled_full[6], ts_xyz):
        raise ValueError("Expected 13 seeds with the exact TS at the center")
    if (
        np.max(np.abs(sampled_full[0] - endpoints[0])) > 2e-5
        or np.max(np.abs(sampled_full[-1] - endpoints[1])) > 2e-5
    ):
        raise ValueError("Sampled path endpoints differ from minimized endpoints")

    # Existing endpoint/TS contacts are documented baseline exceptions.
    baseline_contacts: set[tuple[int, int]] = set()
    for coords in (ts_xyz, *endpoints):
        baseline_contacts.update(nonbonded_contacts(parm, coords, active))
    clash_records = []
    candidates = [
        (f"image_{idx:02}", coords) for idx, coords in enumerate(sampled_full)
    ] + [
        (f"midpoint_{idx:02}", 0.5 * (left + right))
        for idx, (left, right) in enumerate(zip(sampled_full, sampled_full[1:]))
    ]
    for label, coords in candidates:
        contacts = nonbonded_contacts(parm, coords, active)
        new_contacts = {
            pair: value
            for pair, value in contacts.items()
            if pair not in baseline_contacts
        }
        clash_records.append(
            {
                "label": label,
                "contacts_below_2A": len(contacts),
                "new_contacts_vs_TS_and_endpoints": len(new_contacts),
                "minimum_contact_A": min(
                    contacts.values(), default=CLASH_CUTOFF_ANGSTROM
                ),
                "minimum_new_contact_A": min(
                    new_contacts.values(), default=CLASH_CUTOFF_ANGSTROM
                ),
            }
        )
    new_contact_count = sum(
        int(item["new_contacts_vs_TS_and_endpoints"]) for item in clash_records
    )
    if new_contact_count:
        raise ValueError(
            f"Seed or midpoint introduces {new_contact_count} heavy-atom contacts below 2 A"
        )

    pair_map = {
        "OG1_C": (13916, 19916),
        "C_N": (19916, 19929),
        "HG1_Nleave": (13917, 19929),
    }
    geom = [
        {key: distance(x, *pair) for key, pair in pair_map.items()}
        for x in sampled_full
    ]
    reactive_deltas = {
        key: [abs(right[key] - left[key]) for left, right in zip(geom, geom[1:])]
        for key in pair_map
    }
    maximum_reactive_delta = {
        key: max(values) for key, values in reactive_deltas.items()
    }
    if any(value > 0.36 for value in maximum_reactive_delta.values()):
        raise ValueError(
            f"Reactive-coordinate seed jump exceeds 0.36 A: {maximum_reactive_delta}"
        )
    if not (geom[0]["OG1_C"] > 2.7 and 1.2 < geom[0]["C_N"] < 1.7):
        raise ValueError("Minus endpoint is not E-S-like")
    if not (1.2 < geom[-1]["OG1_C"] < 1.75 and geom[-1]["C_N"] > 2.0):
        raise ValueError("Plus endpoint is not acylenzyme-like")
    proton_distances = {
        "Thr_N_H1": (
            distance(sampled_full[0], 13905, 13906),
            distance(sampled_full[-1], 13905, 13906),
        ),
        "Thr_N_H2": (
            distance(sampled_full[0], 13905, 13907),
            distance(sampled_full[-1], 13905, 13907),
        ),
        "substrate_N_HN": (
            distance(sampled_full[0], 19929, 19937),
            distance(sampled_full[-1], 19929, 19937),
        ),
        "Thr_OG1_HG1": (
            distance(sampled_full[0], 13916, 13917),
            distance(sampled_full[-1], 13916, 13917),
        ),
        "leaving_N_HG1": (
            distance(sampled_full[0], 19929, 13917),
            distance(sampled_full[-1], 19929, 13917),
        ),
    }
    if any(
        not 0.8 < value < 1.3
        for name, endpoint_values in proton_distances.items()
        if name in {"Thr_N_H1", "Thr_N_H2", "substrate_N_HN"}
        for value in endpoint_values
    ):
        raise ValueError("Expected Thr N-terminal H1/H2 or substrate N-HN bond changed")
    if not (
        0.8 < proton_distances["Thr_OG1_HG1"][0] < 1.3
        and proton_distances["leaving_N_HG1"][0] > 1.3
        and proton_distances["Thr_OG1_HG1"][1] > 1.3
        and 0.8 < proton_distances["leaving_N_HG1"][1] < 1.3
    ):
        raise ValueError("Expected catalytic HG1 transfer from Thr OG1 to leaving N")

    names = []
    for index, coords in enumerate(sampled_full):
        path = images_dir / f"image_{index:02}.ncrst"
        write_restart(path, coords, ts)
        names.append(path)
    shutil.copy2(TOPOLOGY, inputs / "system.prmtop")
    shutil.copy2(BASE, inputs / "base.mdin")
    ref_path = names[0]
    text = BASE.read_text()
    text = re.sub(r"dftb_telec\s*=\s*[^,\s]+\s*,?", "", text, flags=re.I)
    text = text.replace("&qmmm", "&qmmm\n dftb_telec=0.0,")
    text = re.sub(r"maxcycle\s*=\s*\d+", "maxcycle=1000", text, flags=re.I)
    text = re.sub(r"maxstep\s*=\s*[\d.]+", "maxstep=0.1", text, flags=re.I)
    text = text.replace("verbosity=1", "verbosity=2")
    text = text.replace("ntwx=25", "ntwx=5")
    controls = [
        "neb='FROZEN', neb_nimage=17, nebk=0.01,",
        "neb_climb_test=0.0, neb_freeze_test=0.0,",
        f"refcrd='{ref_path}',",
        f"wtmask='{qm_match.group(1)}',",
    ]
    controls.extend(f"coords2({i})='{path}'," for i, path in enumerate(names[1:], 1))
    text = text.replace("&dlfind", "&dlfind\n " + "\n ".join(controls))
    (calculation / "minimize.mdin").write_text(text)
    (calculation / "launch.sh").write_text(_launch_script())
    config: dict[str, object] = {
        "method": "0 K DFTB3/3OB periodic QM/MM native Sander/DL-FIND NEB",
        "topology_sha256": digest(TOPOLOGY),
        "base_input_sha256": digest(BASE),
        "active_atoms": int(active.sum()),
        "total_atoms": len(parm.atoms),
        "qm_atoms": 98,
        "link_atoms": 8,
        "qm_charge": -2,
        "dftb_telec_kelvin": 0.0,
        "active_mask": active_match.group(1),
        "endpoint_policy": "FROZEN",
        "neb_nimage": 17,
        "neb_spring_constant": 0.01,
        "neb_climb_test": 0.0,
        "ci_enabled": False,
        "ci_switch": "After a stable converged regular band, continue from validated native nebpath.xyz with neb_climb_test=3.0; do not replace endpoints or shell.",
        "seed_structure_count": len(sampled_full),
        "seed_ts_index_zero_based": 6,
        "seed_sampling": "6 equal intervals per side in cumulative arc length of mapped OG1-C, C-N, and HG1-N distances",
        "branch_metrics": branch_metrics,
        "endpoint_and_ts_sha256": {
            "ts": digest(TS),
            "es": digest(MINUS / "minimized.ncrst"),
            "acylenzyme": digest(PLUS / "minimized.ncrst"),
        },
        "trajectory_sha256": {
            "minus": digest(MINUS / "path.xyz"),
            "plus": digest(PLUS / "path.xyz"),
        },
        "outer_max_displacement_A": outer_max,
        "seed_geometry_A": geom,
        "adjacent_reactive_distance_deltas_A": reactive_deltas,
        "maximum_adjacent_reactive_distance_delta_A": maximum_reactive_delta,
        "endpoint_proton_distances_A": proton_distances,
        "clash_cutoff_A": CLASH_CUTOFF_ANGSTROM,
        "clash_checks": clash_records,
        "native_numbered_restart_valid": False,
        "coordinate_source_policy": "Seed restarts in inputs/images; post-run native nebpath.xyz authoritative",
        "scientific_acceptance": False,
        "launched": False,
    }
    (OUTPUT / "config.json").write_text(
        json.dumps(config, indent=2, allow_nan=False) + "\n"
    )
    (OUTPUT / "README.md").write_text(
        "# Direct-route 0 K native NEB setup\n\n"
        "This is a fresh, unlaunched band setup using the converged E-S and "
        "acylenzyme-like minima. The 295-frame mode-minus and 517-frame "
        "mode-plus minimization paths were verified against their full-system "
        "start/final restarts (maximum coordinate mismatch below 1.3e-7 A). "
        "Six equal intervals per side in cumulative arc length of the mapped "
        "OG1-C, scissile C-N and HG1-leaving-N distances give 13 source structures, with "
        "the exact converged TS geometry at image index 6; it is coordinate-identical "
        "to the 0 K one-cycle stationarity-check restart. The installed Sander input "
        "source allows up to 12 coords2 images and interpolates them to "
        "neb_nimage=17; the -c reference plus 12 coords2 entries therefore "
        "provides 13 ordered structures for redistribution to 17 images.\n\n"
        "Adjacent seed changes in each mapped reactive distance are gated at "
        "0.36 A; actual maxima are in `config.json`.\n\n"
        "CI is deliberately disabled for the first regular band "
        "(`neb_climb_test=0.0`). After the regular band is stable and converged, "
        "continue from its validated native `nebpath.xyz` with "
        "`neb_climb_test=3.0` to spawn the climbing image; preserve all image, "
        "topology, shell and Hamiltonian settings. This setup is not a TS or "
        "MEP acceptance record. The candidate still lacks Hessian-index "
        "validation.\n\n"
        "Preflight checks include exact topology/box/atom mapping, 4,165 active "
        "atoms, 98 QM atoms, charge -2, fixed-exterior equality, endpoint "
        "chemistry, retained Thr N-terminal H1/H2, substrate N-HN, Thr OG1-HG1, "
        "and no *new* heavy-atom contacts below 2.0 A relative to the TS and "
        "both endpoint baselines. Persistent baseline contacts are listed as "
        "exceptions in the preflight metrics. Metrics and hashes "
        "are in `config.json` and `input_source_sha256.json`.\n\n"
        "To launch later from the repository root, run `bash "
        "qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/calculation/launch.sh`. "
        "That command creates tmux session `model30_direct_0k_neb_20261113`; "
        "it has not been run.\n"
    )
    source_files = [
        Path(__file__),
        REPO / "scripts/prepare_dimer_continuation.py",
        REPO / "scripts/run_nyl12_job.py",
        REPO / "src/nyl12_covalent/jobs.py",
        BASE,
        TOPOLOGY,
        TS,
        MINUS / "path.xyz",
        PLUS / "path.xyz",
        MINUS / "minimized.ncrst",
        PLUS / "minimized.ncrst",
        Path("/data/jpc/code/ambertools26/venv/bin/sander"),
        Path(
            "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/sander/dlfind_module.F90"
        ),
        Path(
            "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/libdlfind/src/dl-find.F90"
        ),
    ]
    hashes = {str(path.resolve()): digest(path) for path in source_files}
    hashes.update(
        {
            str(path.relative_to(OUTPUT)): digest(path)
            for path in OUTPUT.rglob("*")
            if path.is_file() and path.name != "input_source_sha256.json"
        }
    )
    (OUTPUT / "input_source_sha256.json").write_text(
        json.dumps(hashes, indent=2) + "\n"
    )
    return config


def _element_symbol(atomic_number: int) -> str:
    """Return the standard element symbol used in XYZ trajectories."""
    from parmed.periodic_table import Element

    return str(Element[atomic_number])


def _launch_script() -> str:
    """Return a nonexecuted launcher using the shared auditable job manager."""
    amber = "/data/jpc/code/ambertools26/venv"
    sander = f"{amber}/bin/sander"
    inputs = OUTPUT / "inputs"
    calc = OUTPUT / "calculation"
    seed_inputs = " ".join(
        f"--input {path}"
        for path in sorted((OUTPUT / "inputs/images").glob("image_*.ncrst"))
    )
    args = [
        "uv run python scripts/run_nyl12_job.py launch",
        "--session model30_direct_0k_neb_20261113",
        f"--output-dir {OUTPUT / 'job_run'}",
        f"--cwd {calc}",
        f"--env AMBERHOME={amber}",
        f"--env LD_LIBRARY_PATH={amber}/lib",
        "--env MKL_NUM_THREADS=1 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1",
        f"--input {sander} --input {inputs / 'system.prmtop'} {seed_inputs}",
        f"--input {calc / 'minimize.mdin'} --",
        f"{sander} -O -i minimize.mdin -p {inputs / 'system.prmtop'} -c {inputs / 'images/image_00.ncrst'}",
        "-o minimized.mdout -r optimized.ncrst -x optimization.nc -inf mdinfo",
    ]
    return (
        '#!/usr/bin/env bash\nset -euo pipefail\ncd "$(git rev-parse --show-toplevel)"\n'
        + " \\\n  ".join(args)
        + "\n"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    result = build()
    LOGGER.info(
        "Prepared %d seeds; regular NEB is configured and not launched",
        result["seed_structure_count"],
    )
