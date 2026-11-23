# NEB-US N-terminal relay releases, 2026-11-23

The three approved native DFTB3/MM optimizations started on violet at
2026-11-23 22:21 UTC (17:21 EST). They test whether sampled N-terminal-protonated
structures retain a TI1-like basin after removing the umbrella potentials.
These are geometry optimizations, not a PMF, timed reaction trajectory, or
validated stationary-point classification. Raw seed energies must not be read
as cross-seed barriers: each seed contains a different finite-temperature MM
configuration.

Campaign, relative to the repository:
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/umbrella_sampling/nterminal_relay_releases_20261123/`.

| Seed | Accepted source window/block/frame (1-based) | Source time, ps | C–OG1 / C–N / C–O, Å |
|---|---|---:|---|
| `ti1_window032` | 032 / 04 / 65 | 24.25 | 1.558 / 1.471 / 1.307 |
| `ti1_window033` | 033 / 05 / 46 | 28.30 | 1.556 / 1.478 / 1.317 |
| `early_window028` | 028 / 04 / 98 | 25.90 | 1.897 / 1.495 / 1.275 |

All selected source coordinates have actual HG1–Nalpha bonds (1.028, 1.045,
1.068 Å), H1/H2 still bonded to Nalpha, and H19937 bonded to the leaving N.
OG1 and carbonyl O have no proton among these four tracked hydrogens. Selection
used actual bond distances and all six carbon-centered angles, not nearest-owner
labels. No atom/proton was added, removed, renamed, restrained, or reassigned.
A later H1/H2 donation to the leaving N is physically allowed and will be
reported from geometry. The two TI1 seeds have carbon-centered angle ranges
98.6–120.0° and 101.0–114.3°. Their aligned QM/local-MM RMSDs are 0.792/1.353 Å;
all pairwise QM RMSDs are 0.787–0.802 Å. Local MM means movable atoms outside QM.

The bounded scan inspected locally retained accepted trajectory payloads for
windows 028–034; 17 metadata-only block payloads were explicitly excluded.
It found 2207 HG1–Nalpha-bonded frames, including 1622 with both C–OG1 and C–N
below 1.7 Å. This local subset is not a population estimate. Full extracted
distances, inspection helper and exclusions are in `selection_evidence/`.
The frozen whole-campaign nearest-owner count (3891/37800) also is a raw biased
count, not an unbiased basin population.

Every seed retains the original 104922-atom topology, QM charge −2, 98 real QM
atoms and eight link hydrogens (368 QM electrons including links; even).
The inherited topology partial-charge sum is −0.004004953 e; this small
noninteger residual was preserved. Topology SHA-256 is
`448560852cfe65e86e6ef094f2bfc9551029dd1b3cf9847b17fb1a89116a8f63`.
The exact original 4165-atom 12 Å movable list matches both US and native inputs;
its zero-based little-endian int64 list hash is
`2cdbe5a2d210eb529fbf72d3fa71fc752c88f5f2769260cbeadb46b397a5a6b9`.
QM list hash uses the same convention:
`831552968095827287aad4c43711c513f173eb3ab2e79b0b99b50bffc6a2cb0e`.

“Unrestrained” means no reaction-coordinate or positional potentials
(`nmropt=0`, `ntr=0`) within that unchanged movable boundary; all other atoms
remain frozen at each seed's original coordinates. The original QM partition,
periodic embedding, DFTB3/3ob and zero electronic temperature are retained.
Native Cartesian LBFGS has maxcycle=10000, tol=1e-5 Eh/bohr, tole=1e-7.
AMBERHOME is `/data/jpc/code/ambertools26/venv`; all 16 H/C/N/O ordered-pair SK
files are copied/hashed. Each job has one CPU thread, no GPU allocation.

Tmux sessions are `nyl12_relay_ti1_window032_20261123`,
`nyl12_relay_ti1_window033_20261123`, and
`nyl12_relay_early_window028_20261123`. The project launcher preserves exact
commands, environment overrides, timestamps, hashes and process status in each
`job/`. `prepared.json`, `execution_manifest.json`, `source_identity.json`,
`inputs/`, `parameters_3ob/` and `source_snapshot/` retain provenance and inputs.

The first independent native gradients passed SCC and evaluated-coordinate
identity gates (maximum serialization error ≤5e-8 Å), with no unresolved SCC,
fatal markers or temperature changes. Starting maximum forces are
0.10428/0.09035/0.09858 Eh/bohr for the table order. They are thermal snapshots,
not converged structures. Initial energies are −711.9176456, −712.0976002 and
−711.9701318 Eh; compare only each seed with its own later optimization energy.

The runner independently records process exit, native convergence, SCC history,
maximum final force, within-seed energy decrease, outer-shell displacement,
active displacement, carbon geometry and H1/H2/HG1/H19937 bond alternatives.
Numerical acceptance requires max force ≤1e-5 Eh/bohr and unchanged frozen atoms.
Chemical retention remains a separate geometry assessment. No Hessian is run;
positive stationary-point order is not established by force convergence alone.
Original `optimization.nc`, `path.xyz`, and `path_force.xyz` are retained.

Each seed has a canonical `vmd_progress_verified/` bundle: `view_relay.vmd`,
`optimization.dcd`, `stripped.parm7`, complete original/display atom mapping,
source snapshot, and independently checked float32 coordinates/periodic boxes.
Frame zero is the exact sampled seed; later frames are saved optimizer
iterations, not MD time. Regenerate with
`scripts/export_neb_us_relay_release_vmd.py --seed SEED --output NEW_EXPORT`.
The single `nyl12_relay_movies_20261123` monitor refreshes canonical bundles every 300 seconds and once after termination, with at most 48 frames (one seed plus 47 sampled saved iterations). Its separate immutable source is in `monitor_source_snapshot/`. It holds an exclusive lock on `.vmd_export.lock` across all-three-bundle publication; the laptop pull takes the matching shared lock. Verified files publish into the stable canonical directories; temporary builders stay outside the campaign.
A failed initial export remains in `ti1_window032/vmd_initial_progress/`;
cpptraj's prohibition on parmstrip after trajectory setup was corrected by
separate topology export. The subsequent exports passed all checks. The first automatic monitor exited on a /tmp-to-/data cross-filesystem rename before publication; `monitor_job/` and its source are retained. The corrected monitor uses a checked copy to a same-directory pending file followed by atomic replacement, under the same exclusive lock. Its continuation is recorded in `monitor_continuation_job/` and `monitor_publication_fix_manifest.json`.

Checks: 18 focused and artifact-backed integration tests passed; combined custom-source coverage is 91% (geometry/input 100%, native runner 92%, preparation 86%, exporter 91%, monitor 90%). Ruff lint and format checks passed for all seven new source/test files. Preparation reproducibility, DCD coordinates/boxes, canonical publication, monitor lifetime/failure logging and locking were exercised; retained reports are in `checks/`. No dependency was added. Original artifacts and unrelated files remain untouched; no scientific job was stopped. Final scientific acceptance remains pending.

The corrected monitor's first all-three publication was verified at 38/38/39 frames; all coordinate/box checks passed, the generated dynamic-bond viewer Tcl parsed/executed with VMD command stubs, no pending files remained, and a nonblocking shared flock succeeded between refreshes. Native VMD GUI rendering was not available. The failed initial temporary builder was archived with a verified manifest hash under `clutter/failed_movie_builder_before_portable_publication/`; its /tmp location is absent. The system refused `gio trash` on its internal /tmp mount, so no unsafe deletion was substituted. Exclude campaign `clutter/` from routine laptop pulls.
