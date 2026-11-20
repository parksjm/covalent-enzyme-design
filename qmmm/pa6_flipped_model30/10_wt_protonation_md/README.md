# Matched WT protonation MD comparison

Authorized 2026-11-19. Active tmux controller: `wt_three_protonation_20261119`.
This is a fresh matched comparison, not a continuation of the earlier replicas.

| State directory | All four Thr233 N termini | Asp272 | Asp274 | Na / Cl |
|---|---|---|---|---|
| nh2_asp272_asp274 | NH2 | ASP | ASP | 114 / 106 |
| nh2_asp272_ash274 | NH2 | ASP | ASH, HD2 on original OD2 | 110 / 106 |
| nh3_asp272_asp274 | NH3+ | ASP | ASP | 110 / 106 |

Each state uses the same parent heavy-atom coordinates, 28,245 water molecules
with identical coordinates/order and truncated-octahedron box. Both protonated
conditions remove the same four remote Na ions, selected using periodic minimum
distances. All systems have 104,922 atoms and charge magnitude below 1e-5 e.
Equal atom counts do not imply equal atom indices: extra protons occur before
the PA6 atoms, while removed ions occur later. Use each state's explicit mapping.

ff19SB/TIP3P and the pinned PA6 modular RESP/GAFF2 release are retained.
Stock NTHR is used for NH3+. Both NH2 states use the earlier custom neutral
recipe with its -0.001 e/residue bookkeeping residual corrected explicitly by
changing N from -0.4844 to -0.4834 e. Other inherited neutral-Thr types/terms
are unchanged. This charge normalization is not a new neutral-residue QM
parameterization or an independent validation of its bonded/LJ parameters.
`source_manifest.json`, preparation receipts and `final_inventory_audit.json`
record source hashes, proton inventories and exact shared-coordinate checks.
The failed ASCII-restart serialization attempt was retained; production inputs
use NetCDF precision so the exact coordinate gate passes without relaxation.

The frozen controller first runs fully unrestrained pmemd.cuda minimization
for all three states (10,000-cycle ceiling, target RMS 0.01). Process, finite
coordinates/proton bonds and MD-preparation RMS <=0.1 are required; meeting the
stricter force target is reported separately. No positional/contact restraints
are present in this minimization. Accepted starts then enter the existing
300 K schedule: 100 ps heating, 100 ps density relaxation, 100 ps weak
equilibration, 200 ps unrestrained equilibration, then 10+90 ns unrestrained
production. Initial three preparatory stages retain the documented 5/2/0.5
kcal/mol/A^2 heavy-solute restraints. Production ntr=0, nmropt=0.

Four independent 100 ns replicas per state are queued: 12 runs, 1.2 us total,
with 72 globally distinct stage seeds. Two persistent GPU workers each run
one replica at a time; a failed replica stops its worker explicitly. Queues
use a frozen 12-file protocol snapshot to avoid changes from later code edits.
Read `orchestration_job/status.json`, per-state minimization acceptance,
`queue/` and `md/<state>/<replica>/` for process and scientific records.
Do not rerun `execute` or overwrite this attempt to resume a failure.

Preparation/launch scripts are versioned in the parent repository. Classical
MD keeps each declared protonation topology fixed; it cannot demonstrate
proton transfer. Comparison of Tier A occupancy must also retain native-Thr
geometry, substrate binding, site identity and replica/block uncertainty.
No new umbrella sampling was launched as part of this WT campaign.
