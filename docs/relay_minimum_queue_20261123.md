# Conditional relay local-minimum verification, 2026-11-23

PI-authorized queue: leave all three original unbiased releases running normally.
Only force-accepted **ti1_window033** triggers curvature verification; other jobs
need not finish first. No further path work is dispatched by this controller.

Campaign: existing `umbrella_sampling/nterminal_relay_releases_20261123` under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113`.
One new directory is `ti1_window033/minimum_verification/`, with stages:

- `00_queue/`: source snapshot, exact command/input hashes and tmux launcher job.
- `01_frozen/`: force-accepted restart, original topology/MDIN/calibration,
  precision serialization library, original build receipt and SHA-256 manifest.
- `02_gradient_repeats/`: original serialization plus two exact precise repeats.
- `03_step_probes/`: full-length reactive C, Thr930 HB13925 and first active MM
  xyz columns at h=0.001 and 0.0005 Angstrom (9 columns; 36 gradients).
- `04_full_hessian/`: unchanged `full_hessian(..., expected_negative=0)` using
  all 12495 Cartesian DOF, h=0.0005 Angstrom, fixed exterior.

Original native Cartesian LBFGS jobs retain maxcycle/maxcyc=10000,
force tolerance 1e-5 Eh/bohr and energy tolerance 1e-7. Existing sessions,
movie monitor, trajectories and originals are untouched. Budget exhaustion or
terminal nonconvergence requires separate reviewed restart coordination; this
queue creates no automatic optimizer continuation.

Controller session: `nyl12_relay_minimum_w033_20261123`; polling interval120 s,
explicit **168 h total wall budget**, including an owned verifier subprocess.
At that deadline only the controller-owned verification process group can be
stopped; original optimizer processes are never signaled. The existing movie
monitor continues to record the other releases.

The geometry screen requires C–OG1<=1.8, C–N<=1.7 and C–O>=1.28 Angstrom.
Within0.001 Angstrom of a boundary, stop for manual review. Retain angles,
all four reactive proton alternatives and the complete QM-hydrogen inventory;
no N-terminal NH3+ or particular H1/H2/HG assignment is required. Passing this
screen alone does not establish TI1 or a local minimum.

Prerequisites: exact original4165-atom12-Angstrom active list,98 real QM plus
8 links, QM charge−2,0K electronic temperature, original topology/parameter
hashes; AMBERHOME=/data/jpc/code/ambertools26/venv. All16 installed3ob SKFs must
match the campaign snapshot. Use the existing analysis Python environment.
Copy the validated serialization-only library from
`ts_image08_attempt1/validation_loop_20261121/verified_frame_validation_attempt05/library_build/libdlfind.so`
and the parent `build_receipt.json`, verifying its recorded SHA-256. No library
rebuild, dependency, parameter, proton inventory or QM partition changes.

Numerical gates: exact precise energy/gradient repeats, old/new energy parity
<=1e-9 Eh, gradient parity within original rounding bound0.50001e-7/0.529177²
Eh/Angstrom, force<=1e-5 Eh/bohr; every probe's relative column norm difference
<=1e-4. The step criterion is a declared engineering screen, not a rigorous
operator error bound. Fail/inconclusive stops retain all evidence and never
loosen criteria. Full matrix requires zero robust negatives, zero unresolved
signs and eigenvector residuals<=1e-8. Preserve literal inertia and caveat:
**2*skew is a numerical heuristic, not a rigorous truncation/error bound**.
A pass supports numerical local curvature within that scope, not an analytic
proof, product assignment or sampling convergence.

At most8 nice single-CPU workers, all BLAS/OpenMP pools explicitly1. Initial
Violet preflight passed48 logical CPUs/load4.4/~485GiB available RAM/~877GiB
free disk; resource headroom is checked again at dispatch. Source stage hashes
and the active append-only action journal retain intent, outcome and rollback.

Validation:38 focused pytest regressions pass; coverage86% across the three new
modules (contracts94%, run99%, CLI67%); Ruff check/format pass. Native gradients
are mocked in tests. Live read-only preflight verifies all40 prerequisite
artifacts and exact probe indices; no preconvergence native gradient or Hessian
is launched. The full repository suite is not repeated for this isolated queue.
