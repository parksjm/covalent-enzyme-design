# Model30 MD checkpoint — 2026-11-10

Two neutral replicas are in 10 ns production on violet. Two charged replicas
started 100 ps NVT heating on azalea GPUs 0/1 at 22:29:53/54 EST, seeds
831110/831127. All use the matched protocol in `03_md/*/protocol.json` beneath
`qmmm/pa6_flipped_model30`: heating, 400 ps staged NPT, then 10 ns unrestrained
NPT. Completion and scientific acceptance remain separate.

Charged preparation preserves neutral-minimized heavy coordinates and validates
every atom/type/charge of all four beta N termini against standard ff19SB NTHR.
The dry charge increases by 4.004: standard templates correct the legacy neutral
residual of -0.001 per terminus. Final solvated charge is -0.000004953 e.
Re-solvation uses TIP3P, 12 A buffer and 150 mM salt (123 Na+, 119 Cl-).

Charged unrestrained MM minimization attempt2 completed 10000 cycles, exit 0,
RMS gradient 0.10032 and GMAX 5.7733; it did not meet drms=0.001. Its projected
phi is 118.756 degrees, Tier B, native=false, scaffold RMSD 0.342 A,
OG1 shift 2.081 A and N-terminal shift 2.586 A. The charged MD is an explicitly
authorized diagnostic protonation comparison, **not an accepted catalytic
seed**. Exact measurements and input hashes are in
`02_charged_system/01_mm_minimization_attempt2/endpoint_geometry.json`.

Azalea uses the existing tested executable
`/data/jpc/code/pmemd26_blackwell_build_20261005/pmemd26_src/build_sysgcc_cuda128_v2/src/pmemd/src/pmemd.cuda_SPFP`
and its previously validated library environment. The first charged min attempt
with the default incompatible GPU executable failed before minimization and is
preserved. No software was installed or rebuilt.

The independent NEB job ended after 437 cycles at its energy evaluation limit;
peak 86.939 kcal/mol at image 8, not converged. It was left stopped.

## Overnight extension authorized

Charged replicas each have a persistent 90 ns checkpoint continuation queued
after successful acceptance of the current 10 ns segment, for 100 ns total.
Verified both `extension_100ns.json` states are
`queued_waiting_for_accepted_10ns` and tmux sessions
`model30_charged90_r1_20261110` / `model30_charged90_r2_20261110` are active.
Continuation seeds are 881110 / 881127. The exact production protocol is kept;
coordinates and velocities come from the 10 ns restart without reheating.
Measured Blackwell rates are 301.83 / 298.79 ns/day: approximately 8 hours per
100 ns, with completion around 06:30–06:35 EST on November 11 if rates hold.
Queue provenance is `03_md/charged_100ns_queue.json` on azalea.

## Overnight extension authorized — 22:46 EST

Both neutral replicas retain their current 10 ns trajectories and have queued
90 ns checkpoint continuations, for **100 ns total production each**. Confirmed
watchers are `nyl12_model30_neutral01_to100ns_20261110` and
`nyl12_model30_neutral02_to100ns_20261110`; live queue states are each replica
folder's `extension_100ns.json`. They require accepted 10 ns completion before
using the exact hashed `05_production_10ns.rst7`; new outputs are
`06_production_90ns.*` in the same replica directories. Current MD is unchanged.

The continuation retains 2 fs SHAKE-H, Langevin 300 K, Monte Carlo NPT and
`ntr=0,nmropt=0`; only segment length (45 million steps) and independent
Langevin seeds (781110/781127) change. Source snapshots and queue manifests are
`03_md/neutral_extension_source` and `03_md/neutral_100ns_extension_manifest.json`;
job records are beneath each existing `03_md/jobs/neutral_replica_0*/extension_to_100ns`.
Measured speeds were 280.89/275.96 ns/day: approximately 8.2/8.3 hours remained
when queued, implying completion around 07:00 EST on November 11 if speed holds.
Four focused protocol/continuation tests and Ruff passed; coverage is unmeasured.
