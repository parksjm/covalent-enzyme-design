# Software, environments and reproduction

## Required for new MD

Use Linux with a CUDA-capable NVIDIA GPU and a compatible driver. Install an
Amber distribution providing **pmemd.cuda**, plus AmberTools **tleap,pdb4amber,
cpptraj** and the standard ff19SB/TIP3P/GAFF2 data. The actual source MD headers
identify Amber26 PMEMD Release26; the recorded cpptraj is V7.6.2. The export
contains generated scientific inputs, not scientific-software executables.
Obtain an appropriate distribution and review its terms at
https://ambermd.org/GetAmber.php and https://ambermd.org/Installation.php .

Install Git,rsync,SSH and tmux through the Linux distribution. Install uv from
https://docs.astral.sh/uv/getting-started/installation/ . Do not reuse the PI's
absolute environment paths or modify machine mounts/drivers to reproduce them.
The source filenames retain those paths as provenance, not installation orders.

## Python environment

`analysis_environment/pyproject.toml` exports the established analysis package
versions: Python3.12,NumPy1.26.4,ParmEd4.3.0,MDAnalysis2.10.0,
matplotlib3.10.8,SciPy1.16.3,pandas2.3.3,Biopython1.86. RDKit2025.3.6
matches the repository's declared graph-preparation dependency. The original
analysis environment JSON and Pixi manifest/lock are also retained for audit.
No program installation was performed when building the handoff.

From the extracted bundle root:

```bash
uv sync --project analysis_environment --locked
source /path/to/amber/amber.sh
export AMBERHOME=/path/to/amber
export PMEMD_CUDA=/path/to/amber/bin/pmemd.cuda
uv run --project analysis_environment python tools/run_prepared_md.py plan \
  --repo "$PWD/project_code" --config tools/config_wt_asp272_ash274.json
```

The archive includes a small `project_code/` source snapshot. Alternatively,
clone https://github.com/parksjm/covalent-enzyme-design and give `--repo` that
checkout. The exact file hashes in the package, not an unpinned branch tip,
define the delivered code. A live checkout can contain later changes.

Plan mode validates supplied systems and writes frozen plans; it launches no MD.
Review GPU IDs,seeds and output paths. The actual launch is a separate command:

```bash
uv run --project analysis_environment python tools/run_prepared_md.py launch \
  --repo "$PWD/project_code" --config tools/config_wt_asp272_ash274.json
```

Use one active replica per GPU. The two-replica starter configs form wave1;
run the supplied primary/baseline wave2 configs only after wave1 has finished
to obtain four independent replicas without oversubscribing the same GPUs.
Use the same allocation strategy for each condition. Never replace a missing
GPU engine with CPU execution or another force field silently.

A short hardware-stability pilot is appropriate after installation; successful
installation tests do not replace study-specific topology/charge/geometry checks.
Conversely, this handoff does not require rerunning an installation-wide suite
on every student run. GPU trajectories need not reproduce coordinates bitwise
across architectures; match methods and compare independent ensembles.

## Prior analysis and figures

Cached six-replica geometry tables and summaries are available in the core
archive without downloading the approximately75.55 GB raw trajectories. The
portable raw-analysis manifest is relative to the extracted root; its raw
trajectory paths are intentionally absent until separately obtained. The
recorded command in `analysis/compare_analysis_command.json` runs the exact
site-AE comparison implementation and must be run from the bundle root.

The existing ASH272 matched SI panel,generator and supporting tables are
included with their original labels; the package does not invent an ASH274
SI panel that was not located. Prior WT/mutant figures are historical style
references. First reproduce supplied counts/plots, then generate consistently
styled six-metric time-series/histogram panels for new trajectories, retaining
signed phi. The study plan explicitly requires a four-site-aware analysis
extension before design ranking; the existing AE script is not that extension.

Full raw trajectories can be requested using the source paths and receipt
metadata in the manifest. They are optional transfers, not hidden missing files
in the archive checksum list. Imaging,atom names,residue mapping,frame cadence
and scissile-amide convention must be checked before comparing results.

## Visualization and optional later tools

VMD is available from https://www.ks.uiuc.edu/Research/vmd/ . Extract complete
VMD bundles and source their `view_path.vmd`; keep their topology,trajectory,
Tcl files and checksum manifest together. These QM/MM movies are reaction-path
candidates, not verified transition-state trajectories or PMFs.

LigandMPNN,Boltz/AF3 and QM/MM engines are not required for the initial MD
comparison. Add them only for the later explicitly scoped design/mechanism arm,
using pinned releases,approved models and the protected/symmetric residue map.
The full force-field parameterization workflow remains in the separate
polyamide repository; no new charge fitting is required to rerun PA6 L4 MD.
