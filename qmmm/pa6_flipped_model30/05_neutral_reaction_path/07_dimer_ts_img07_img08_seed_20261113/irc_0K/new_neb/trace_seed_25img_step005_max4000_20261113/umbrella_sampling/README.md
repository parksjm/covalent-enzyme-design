# Direct acylation NEB-seeded umbrella sampling

Authorized 2026-11-19. Read `status.json`, `workers/*/status.json` and native
stage `process.json`/`acceptance.json` for execution status. A submitted worker
is not evidence of completed MD. Initial filename-buffer failures before any
MD step are retained in `preparation_history/`; no failed samples enter MBAR.

## Scientific identity

The parent regular NEB converged at cycle 477 with a physical profile maximum
of 35.8089536393 kcal/mol (image 8). It connects E–S to acylenzyme and has not
provided a verified first-order saddle. This is the historical neutral-Thr-N,
ASP272/ASP274 system: 104,922 atoms, 98 real QM atoms, charge −2, singlet,
DFTB3/3ob, electronic temperature 0 K. Use `prepared.json` for exact masks,
source/checkpoint/executable hashes and one-based atom identities.

The exact original 12 Å movable region contains 4,165 atoms; the other
100,757 atoms remain frozen. There are no global positional restraints.
The resulting free energies are conditional on this exterior and declared
proton inventory; they are not an unrestricted full-protein PMF. The native
belly GROUP selection reproduces this boundary without the 256-character
mask truncation. A native probe verifies zero frozen-atom displacement and
thermostat RNDF=12,495 (3 × 4,165), not full-system degrees of freedom.

## Independent chemical biases

- `q_acyl = d(C,N_leaving) − d(C,OG1)`.
- `q_H = d(OG1,HG1) − d(N_leaving,HG1)`.

Each receives a separate native Amber RK=100 kcal/mol/Å² term. Amber uses
`E = k*(q−center)^2`, without a factor of one half; the native linear tails
outside ±10 Å are included in cross-window bias evaluation. NEB springs are
removed. No Cartesian image restraints or N–H holding restraints are used.
The chemical coordinates fluctuate independently about their two centers.

There are 63 windows, retaining all 24 original image centers and inserting
CV centers until adjacent 2D center distance is at most 0.10 Å. Coordinates
are never interpolated: every seed is the nearest original, validated full
NEB restart. `neb_image_centers.csv`, `seed_cv_audit.json` and `prepared.json`
retain the source/target offset, original gaps and all identities. This is
sampling around a path-shaped ribbon, not guaranteed complete 2D coverage.

## Execution and analysis

Eight single-thread CPU workers share the queue and preserve existing GPU MD,
Gaussian and NEB jobs. Each window runs 1 ps heating (25→300 K), 5 ps
300 K equilibration, then six 5 ps production blocks: at least 30 ps/window.
Time step is 0.5 fs, NVT/Langevin gamma=2 ps⁻¹, no SHAKE, including QM bonds.
Heating acceptance has an explicit 25–350 K final-temperature gate because
finite thermostat response lags the ramp; equilibration/production use
250–350 K. Gate acceptance is not proof of thermal equilibration. Production
analysis must separately inspect temperature and stationarity after equilibration.
All stage Langevin seeds are distinct and pinned in the prepared manifest.

Retain all finite electronically converged samples, including unexpected
proton transfers. Numerical/electronic failures stop the affected window;
chemical branches are flagged rather than censored. Full-precision NetCDF
coordinates generate CVs and exact cross-window biases; rounded DUMPAVE
records are diagnostics only. Track all three Thr-N hydrogens, transferring
HG1, and a permutation-invariant coordination-number progress variable over
all nine exchangeable real-QM N/O hydrogens; monitor Thr-N coordination too.
The N/O owner search is scoped to the declared QM atoms and initial parents,
not every MM water in the box. Link hydrogens are not exchangeable particles.

`analysis/latest/` publishes CPU-only MBAR coverage, correlation/ESS, overlap
and interim free-energy surfaces after completed blocks. It analyzes both
(q_acyl,q_H) and the proton-permutation-invariant observable. It suppresses
a global surface while the sampled overlap graph is disconnected, preserves
unsampled bins, and makes no automatic scalar barrier/convergence claim.
At least 30 ps is a minimum allocation; extensions depend on connected overlap,
independent effective samples, chemical branches and block/replica stability.

## Finding and reproducing files

- Native inputs, logs, restart and full trajectory: `window_*/<stage>/`.
- Completed-block compact structural views: `window_*/progress/block_*/view.vmd`.
  These contain stage endpoints and are explicitly **not full MD movies**.
- Small plots and numerical diagnostics: `analysis/latest/`.
- CPU queue status and retained tmux job records: `workers/`.
- Frozen worker source/environment/input identities: `source_snapshot/`,
  `prepared.json` and per-job manifests. Read source snapshots for this run,
  not a later checkout of similarly named tools.
- Stable entry point: repository `scripts/run_neb_umbrella.py`; analysis:
  `scripts/analyze_neb_umbrella.py`.

The laptop `scripts/pull_nyl12_inspection.sh` preserves this same tree and
pulls small current diagnostics and complete compact views. Raw full-system
trajectories, restarts and preparation/history directories remain on violet.
No new production dependencies were installed for this campaign.
