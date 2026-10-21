# Nyl12 laptop inspection package — 2026-10-21

This package contains small structures, validation reports, documentation,
manuscript drafts and provenance. Gaussian checkpoint/scratch files, raw logs,
force-field topologies and simulation trajectories are excluded. The included
XYZ files are small coordinate histories of the two Gaussian calculations.

## Download from your laptop

Copy `scripts/pull_nyl12_inspection.sh` from the workstation to your laptop,
then run it in a terminal. Replace `jpc@violet` with your usual SSH destination
or alias if needed:

```bash
bash pull_nyl12_inspection.sh --dry-run jpc@violet ./nyl12-inspection-20261021
bash pull_nyl12_inspection.sh jpc@violet ./nyl12-inspection-20261021
```

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

Open `protein/pa6_l2_ti1_starting.pdb` and
`protein/poly_xama_maxa_ti1_starting.pdb`. To use the saved view, change PyMOL's
working directory to the package's `protein` directory and run
`@view_models.pml`. It retains the requested commands:

```pml
show sticks, resi 113 +186 and chain A
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
  The retained `poly_xama_maxa` filename is a stable historical identifier;
  its audited chemical sequence is **MAXA|MAXA**. Confidence is particularly
  poor in the tails, so their displayed placement should not be treated as fixed.

`reports/` contains original active-site metrics and ligand-confidence evidence.
No noncovalent ES, true TS1 or acyl-enzyme validation model is included yet.
