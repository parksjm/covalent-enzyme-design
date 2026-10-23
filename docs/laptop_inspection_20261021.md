# Nyl12 laptop inspection package — 2026-10-21

## Completed TI1 reference results — 2026-10-23

The latest export adds `qm/ti1_bonded_reference_completed/`, with four
case directories containing final 37-atom XYZ/SDF coordinates, all20 explicit
hydrogens, atom maps and per-case assessments. After running your usual pull
script, change to a case directory and run `@view_final.pml` in PyMOL.
These constrained capped models intentionally contain no relay water.

All four geometry/frequency/conditional-stability checks passed, but the R
cap-pair bond-length sensitivity failed the predeclared0.03A limit. Final
coordinates are inspectable references, **not validated native parameters or
accepted enzyme TI1 models**. Compare R/−120 against R/−60, then the S pair;
inspect the short backbone N–H/oxyanion contacts. Source fchk files retain
full coordinate precision and are excluded from the laptop transfer.

Read `docs/ti1_qm_results_20261023.md` for the numerical results and
`docs/af3_water_placement_decision_20261023.md` for the prior water-test status
and why improved AF3 guidance is not a current prerequisite. The separately
labeled interim trajectories below remain as historical snapshots.

## TI1 Gaussian reference progress — 2026-10-22

The refreshed package includes `qm/ti1_bonded_reference_progress/`: four small
37-atom capped-anion references with **20 explicit hydrogens each**. These are
unconverged Gaussian coordinate snapshots captured at about 19:45 UTC, not
accepted minima or validated force-field parameters. Relay water is deliberately
absent from this local bonded-fitting surrogate; the enzyme model retains it.

Run `./pull_nyl12_inspection.sh` again. In PyMOL, change to a case directory,
for example `qm/ti1_bonded_reference_progress/attack_S_cap_phi_minus60/`, then
run `@view_progress.pml`. The viewer labels Thr OG1, the tetrahedral C*, the
oxyanion and leaving-group N. Use the movie controls for all saved tables;
VMD can read `trajectory.xyz`. `latest_geometry.xyz`/`.sdf` contain the latest
table, which may be an unevaluated trial. SDF bonds represent the intended graph.

Inspect the S/-60 case against R/-60 first. Both S references show short
oxyanion contacts to the capped Thr backbone N–H (about 1.77–1.93 angstrom at
capture), potentially the gas-phase H bonding that motivated torsion constraints.
Frozen phi/psi/chi1 do not by themselves exclude such contacts. No proton
transfer is detected in the latest snapshots, and the reactive bonds remain
within the declared interim screens. One earlier S/-120 trial shortened C*–O
to 1.225 angstrom and failed that screen; it is retained as frame 4.

Each case includes `trajectory_metrics.csv`, `geometry_checks.json`,
`atom_map.json` and the exact `input.gjf`. The four original live Gaussian logs
are under `/data/jpc/covalent-enzyme-design/output_nyl12_ti1_parameters_20261022/pa6_local_bonded/<case>/optimization_01/optimization.log`.
The snapshot audit is under `output_nyl12_parameter_validation_20261022/progress_audit_01/`.
The raw logs/checkpoints are excluded from the laptop transfer.

## Latest local-repair diagnostic

Viewer correction on 2026-10-22: both water hydrogens were present in all 100
trajectory states but hidden by the original display commands. The corrected
viewers show water sticks explicitly. A semicolon in a comment caused PyMOL to
parse the remaining comment text as Python; the comment is corrected. Actual
headless PyMOL 3.1.0 checks now verify O/H1/H2, both O–H bonds and visible H atoms
in all 100 states. Reloading the viewer replaces its own object instead of
appending another 100 states. Run the pull script again, then rerun
`@view_optimization.pml` (or `@view_local_repair.pml`) from its directory.

The workstation `optimized.pdb` was originally exported under the descriptive
name `pa6_l4_dftb3_diagnostic.pdb`. The refreshed package includes **both names**
with identical contents in `diagnostics/local_repair/pa6_l4_alpha_beta/`.

For the actual optimization snapshots, run `@view_optimization.pml` from that
directory in PyMOL, then use the frame slider or `mplay`. The file
`optimization_trajectory.pdb.gz` contains100 full-system frames of5054 atoms.
It comes from a separate replay with only coordinate-output settings changed;
its final coordinates match the original run exactly. Frames preserve actual
Sander minimizer snapshots, including possible trial/line-search or repeated
points; frame numbers are not MD time or guaranteed accepted iterations.
Source hashes and frame counts are in `trajectory_report.json`. No snapshots
were reconstructed by interpolation. Future MM and QM/MM stages now save
coordinate output explicitly. PyMOL accepts compressed PDB input; see its
[load documentation](https://pymolwiki.org/Cmd_load). For VMD, decompress a copy
first with `gzip -dk optimization_trajectory.pdb.gz` and load the resulting PDB.

Run the existing `./pull_nyl12_inspection.sh` again to retrieve the refreshed
package. In PyMOL, change to `diagnostics/local_repair/pa6_l4_alpha_beta/` and run
`@view_local_repair.pml`. The DFTB3/MM output is
`pa6_l4_dftb3_diagnostic.pdb`; the hydrogenated AF3 start is alongside it.
Water sterics improved, but the OG1–substrate bond stretched to 2.195 A and the
directional hydrogen-bond checks failed. This is a **failed TI1 diagnostic**,
not an accepted optimized intermediate. Acceptance/process reports accompany
the coordinates. Both AF3 water tests are in `diagnostics/af3_water/`.

This package contains small structures, validation reports, documentation,
manuscript drafts and provenance. Gaussian checkpoint/scratch files, raw logs,
force-field topologies and simulation trajectories are excluded. The included
XYZ files are small coordinate histories of the two Gaussian calculations.

## Download from your laptop

Copy the script from violet to your laptop, then run it from the directory where
you want the download:

```bash
scp -p violet:/data/jpc/covalent-enzyme-design/scripts/pull_nyl12_inspection.sh .
./pull_nyl12_inspection.sh
```

It connects as your current laptop username to `violet` and downloads into
`./nyl12-inspection-20261021/`. `./pull_nyl12_inspection.sh --dry-run` previews
the transfer. Remote and destination arguments remain optional overrides.

Only the curated `output_laptop_inspection_20261021/` directory is transferred.
The script uses the laptop's existing `rsync`, SSH and `shasum` or `sha256sum`.
It verifies the downloaded files against `MANIFEST.sha256` and does not delete
local files. `provenance.json` maps every download to its original workstation
source and SHA-256 digest. The export is a dated snapshot; subsequent scientific
results require a new explicitly prepared snapshot.

## Optimized QM structures

Inspect `qm/neutral_thr_nme/conformer_01/optimized.sdf` and the equivalent
`conformer_02/optimized.sdf`. These contain **final optimized coordinates**;
the original Gaussian preparation SDF files contain starting coordinates and
are deliberately not substituted for the optimized exports. XYZ alternatives,
atom maps, Gaussian inputs, acceptance records and the conformer comparison
report are alongside them.

These are isolated **neutral free-NH2 L-Thr-NME models**, not enzyme-bound TI1
clusters. Geometry/frequencies used wB97XD/def2TZVP, with no additional
dispersion correction. Separate HF/6-31G(d) calculations supplied the ESP.
Both retain 21 atoms, neutral charge, CA S / CB R stereochemistry and 57
positive frequencies. No RESP charges have yet been fitted or transferred.
Conformer 02 lies 4.6017 kcal/mol above conformer 01 in electronic DFT energy;
this is not a thermal free-energy or population estimate.

From the package directory, use these PyMOL commands:

```pml
load qm/neutral_thr_nme/conformer_01/optimized.sdf, thr_01
load qm/neutral_thr_nme/conformer_02/optimized.sdf, thr_02
show sticks, thr_01 or thr_02
util.cnc
disable thr_02
zoom thr_01
```

Check the free amino group has two H atoms, Thr OH retains its proton, the
N-methylamide cap remains intact, stereochemistry is correct, and there is no
unexpected proton transfer or distorted covalent geometry. Compare backbone
and side-chain conformations and intramolecular hydrogen bonds. A gas-phase
minimum alone does not establish the preferred enzyme conformation or parameter
transferability. Use SDF for explicit bond connectivity: XYZ has no bond orders.

## Gaussian coordinate histories

The package includes `gaussian_log_to_xyz.py`, which requires Python 3.10 or
newer and uses only its standard library. For another log, run:

```bash
python3 gaussian_log_to_xyz.py calculation.log trajectory.xyz
```

Automatic selection takes the Standard-orientation stream when present;
otherwise it takes Input orientation (needed for the current NoSymm jobs).
Use `--orientation input` or `--orientation standard` to select explicitly.
The output preserves repeated frames: frequency and ESP sections can repeat
coordinates, so frame numbers are **not guaranteed to equal optimization
iteration numbers**. These are coordinate snapshots, not a physical-time MD
trajectory or normal-mode animation. The JSON sidecar records source hashes,
frame counts and selection details. Malformed tables fail; `--allow-incomplete`
explicitly permits skipping an unfinished final table in a running log.

The current histories are already converted (64 frames for conformer 01;
24 for conformer 02). For example, in PyMOL:

```pml
load qm/neutral_thr_nme/conformer_01/trajectory.xyz, thr_path, discrete=0
show sticks, thr_path
util.cnc
mplay
```

Alternatively, load `optimized.sdf` first and use
`load_traj qm/neutral_thr_nme/conformer_01/trajectory.xyz, thr_01, state=1, format=xyz`
if your PyMOL build exposes XYZ through its trajectory plugin. In VMD's Tk
Console, use:

```tcl
mol new qm/neutral_thr_nme/conformer_01/trajectory.xyz type xyz waitfor all
```

Use the viewer's frame controls to examine the optimization history. Viewer
bond perception is a display convention; it is not chemical validation.
References: [PyMOL load](https://www.pymolwiki.org/Load),
[PyMOL load_traj](https://www.pymolwiki.org/index.php/Load_traj),
[VMD mol commands](https://www.ks.uiuc.edu/Research/vmd/vmd-1.8.2/ug/node123.html).

## Protein structures

Open `protein/pa6_l4_ti1_starting.pdb` and
`protein/maxa_maxa_ti1_starting.pdb`. PA6 L4 is the substrate; L2 is its
hydrolysis product. The other substrate is MAXA|MAXA. Historical identifiers
remain in original source paths and provenance records only.
To use the saved view, change PyMOL's
working directory to the package's `protein` directory and run
`@view_models.pml`. It retains the requested commands:

```pml
set valence, 0
show sticks, resi 113+186 and chain A
util.cnc
```

These are **unrefined, hydrogenated covalent TI1 starting models**. Neither is
an accepted optimized TI1 minimum. Their relay-water steric screen failed,
and current classical parameters remain provisional. Inspect:

- Thr **E1** (full sequence Thr233; Amber residue 929): N has three H atoms;
  OG1 has no H and is covalently bonded to the substrate carbon.
- Phe **A113** backbone N and Asn **A186** ND2 around the substrate oxyanion.
- Ligand **I1**: PA6 C54/O49/N59 or poly C033/O034/N035; the substrate C–N
  bond remains intact in TI1. Relay water is **J1**.
- Crowding near Asn186 and the relay water, substrate fit and uncertain tails.
  Confidence is particularly
  poor in the tails, so their displayed placement should not be treated as fixed.

`reports/` contains original active-site metrics and ligand-confidence evidence.
No noncovalent ES, true TS1 or acyl-enzyme validation model is included yet.

Water repair and all exports are scripted; see
`docs/relay_water_repair_20261021.md`. The intended TI1 bridge is
Thr233 N–H…water O–H…Thr233 OG1. A replacement structure requires both
directional hydrogen-bond and full-environment clash acceptance. Renaming alone
does not establish that a water position has been repaired.

The pull script preserves local files. Previously downloaded files bearing the
historical names may remain beside the corrected names; the regenerated viewer
loads only the canonical filenames above.


### Parallel preparation checkpoint, 2026-10-22

The `reports/ti1_parameter_preparation/` directory contains graph-only native
R/S atom mappings, proposed RESP constraints, derivative checks and the clean
Amber source-unit verification. These are preparation reports, not fitted
parameters or optimized native structures. Updated handoff, protocol and
Methods remain under `docs/` and `manuscript/`. Existing inspectable models
and trajectories are unchanged by this report-only update.


### Force-field parameterization Methods draft

After running `./pull_nyl12_inspection.sh`, open
`manuscript/METHODS_DRAFT.md`, section "Force-field parameterization".
Its bibliography is `manuscript/REFERENCES.md`; the implementation/status
record is `manuscript/SOURCE_LEDGER.md`. The draft distinguishes completed
reference preparation from planned GAFF2/RESP fitting and validation.


### Backbone, symmetry and access strategy

The current design recommendation and limitations are in
`docs/redesign_backbone_symmetry_and_access_20261022.md`. It compares the
LigandMPNN-first baseline, optional fixed-length local RFD3 and a separate
future indel branch, and records the D2-like scaffold audit and loop-access
requirements. This is a strategy record, not a completed design campaign.

## Native TI1 starting seed — 2026-10-23

After the usual zero-argument pull, inspect
`qm/native_ti1_seed/native_seed.sdf` or run its `view_seed.pml` from that
directory. It has all60 explicit hydrogens, including two on water, and
geometry-derived R at C54. This is an unoptimized103-atom starting cluster;
read preflight.json and coordinate_provenance.json. It is not the final
Gaussian model and not a validated force field. The frozen-scaffold Gaussian
pilot is running on violet; no final native trajectory is available yet.
The protocol is included at docs/native_ti1_reference_20261023.md.

## Current truncated TI1 pilot — 2026-10-23

Run the usual zero-argument pull, then inspect
`qm/truncated_ti1_seed/native_seed.sdf` or run `view_seed.pml` in that
directory. This 42-atom unoptimized starting cluster contains 25 explicit H,
including both water H, native Thr NH3+, the oxyanion and C54 R. The retained
atoms match the previous seed; two alkyl boundaries are methyl capped.
Read atom_map.json, coordinate_provenance.json, preflight.json and input.gjf.
The 103-atom Gaussian pilot was stopped at the user's request; its seed is
historical. The 42-atom optimization is active on 16 physical cores. No final
optimized 42-atom structure or trajectory is available at this checkpoint.

## Current water-free TI1 starting candidate — 2026-10-23

Run the usual zero-argument pull, then inspect
`qm/waterfree_ti1_seed/native_seed.sdf` or run its `view_seed.pml`.
This is the unoptimized 39-atom, water-free candidate (23 explicit H, net
charge 0, Thr NH3+, oxyanion, C54 R), not an accepted final structure.
Input, atom map, coordinate provenance and preflight are included. The prior
42-atom Gaussian job failed imposing constraints; its starting model is
historical. The current job uses explicit GIC constraints on 16 physical cores.
No final water-free optimized structure or trajectory exists at this checkpoint.

## Water-free optimization movie — 2026-10-23

The current snapshot contains14frames,39atoms and23explicitH per frame.
After the usual pull, change PyMOL directory to
`nyl12-inspection-20261021/qm/waterfree_ti1_optimization`, then run
`@view_optimization.pml` and `mplay`. VMD can load optimization_trajectory.xyz.
These are in-progress geometry snapshots, including possible trial/repeated
points. PDB bonds retain the input graph for display; consult the metrics.
At frame14 OG1-C54=2.02549A and C54-O49=1.21644A, outside TI1 windows.
This is not an accepted TI1 reference; parameter fitting remains blocked.
The immutable log snapshot was captured at20:00:43UTC(16:00:43EDT).
