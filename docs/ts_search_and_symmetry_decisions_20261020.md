# TS search, AF3 restraints and symmetry decisions — 2026-10-20

## Scope and status

This records user instructions, inspected implementation evidence and proposed
methods. No TS1 search, Hessian calculation or optimizer benchmark was run for
this assessment. No dependency was added. Native DL-FIND is acceptable to the
user if freely available; newer alternatives were also requested for comparison.

## Native Sander capability: verified in the installed source

Installation: `/data/jpc/code/ambertools26/venv`.
Source root: `/data/jpc/code/ambertools26/build_assets/src/ambertools26_src`.

- `AmberTools/src/sander/sander.F90`, `ntmincase`, case 5 calls `run_dlfind`.
  The entry point is `imin=1, ntmin=5`, with a separate `&dlfind` namelist.
- `AmberTools/src/sander/dlfind_module.F90` exposes `optalg='PRFO'`,
  `hessini='TWOPOINT'`, `hessupd='BOFILL'`, `maxupdate`, `delta`, `tol`,
  coordinate representations and an explicit movable-atom `active` mask.
  PRFO maps to library `iopt=10` and is marked as a TS search.
- The Hessian callback calls the same Sander energy/gradient evaluator at
  positive and negative Cartesian displacements, forming central differences
  and symmetrizing the Hessian. This includes the chosen QM/MM Hamiltonian;
  it does not require an analytic QM/MM Hessian.
- `optalg='HESS'` exports a fresh Cartesian Hessian to NetCDF; `hess_nseg`
  and `hess_iseg` support separate column segments. These are capability
  observations, not evidence that a segmented workflow has been tested here.
- The installed `libdlfind/src/dlf_formstep.F90` diagonalizes the Hessian for
  P-RFO and initializes `follow=0`, selecting the lowest mode. Arbitrary
  user-selected overlap tracking is not exposed by this inspected interface.
  The Sander interface also sets `imicroiter=0`: do not promise automatic
  outer-MM microiterations merely because general DL-FIND supports them.
- `ldd` with `LD_LIBRARY_PATH=/data/jpc/code/ambertools26/venv/lib` resolves
  `libdlfind.so`. The bare shell lacked this library path; retain the explicit
  environment used by the existing job runner.

A minimal *control fragment*, not a complete or launch-ready scientific input:

```text
&cntrl
  imin=1, ntmin=5,
/
&dlfind
  optalg='PRFO', hessini='TWOPOINT', hessupd='BOFILL',
  maxupdate=5,
/
```

The actual QM/MM settings, coordinates, active mask, restraints, tolerances,
coordinate representation and displacement-convergence study must be supplied
explicitly before use. In particular, do not accept the default all-atom active
mask for the approximately 20,000-atom model without evaluating Hessian cost.

## Saddle search and acceptance protocol — proposed

1. Repair and validate TI1; obtain ES/path guesses by reversing bond formation
   and proton transfer. A TI1 minimum is a hypothesis until verified. Sample
   coupled OG1–C approach and proton coordinates; allow additional stationary
   points if attack and proton transfer are stepwise.
2. Use a driven path or NEB to obtain a starting geometry near TS1. A converged
   climbing-image NEB can locate a saddle without a full Hessian; lack of a
   Hessian alone does not make it mathematically incapable of doing so. Our
   acceptance nevertheless requires independent curvature and connectivity
   verification, rather than identifying the highest sampled image as TS1.
3. Refine with P-RFO using the same DFTB3/3ob QM/MM energy/gradient model.
   Remove reaction-coordinate driving restraints. Define all movable atoms,
   including essential MM environment atoms, and the fixed outer environment.
   Recompute initial/periodic Hessians and inspect the unstable mode to avoid
   converging a tail rotation instead of the intended chemistry.
4. Recalculate a fresh Hessian at the converged point. Require small gradients
   for every unconstrained coordinate and exactly one significant negative
   eigenvalue after appropriate mass weighting and coordinate projection.
   Test displacement and SCF-tolerance sensitivity; tiny numerical negative
   modes cannot simply be ignored without evidence. Retain the full spectrum,
   gradient norms, displacement sizes, units and animated reaction mode.
5. Follow both signs of the unstable mode, preferably with a mass-weighted
   steepest-descent/IRC implementation, then minimize endpoints and audit their
   chemistry. If only displaced downhill minimizations are used, report them
   as endpoint tests, not as a formal IRC. The inspected Sander namelist does
   not establish an available IRC driver. Document any separate path driver.
6. A fixed-environment active-space Hessian establishes a first-order saddle
   only in that declared subspace. A QM-only Hessian does not certify mobile
   MM coordinates. Expand the movable environment and assess sensitivity.
   Frozen environments do not automatically have six free rigid-body modes;
   project only actual invariant/constraint directions. A stationary point on
   a restrained potential is not automatically one on the unrestrained PES.
7. Benchmark the reaction chemistry against higher-level calculations before
   quantitative barriers. A DFT single point on a DFTB saddle does not make it
   a DFT saddle. Electronic-energy paths are not free-energy rate predictions.

A central-difference Hessian costs approximately 6N gradient evaluations for
N movable atoms. The 109/181-atom QM regions alone imply about 654/1,086 calls;
including a mobile MM shell increases these counts. This motivates evaluating
partial-Hessian iterative methods, without omitting final curvature checks.

## Free availability and newer alternatives — checked 2026-10-20

DL-FIND is freely available under GNU LGPL, with no purchase or license key.
Its license is not the separate ChemShell distribution license. The installed
library includes `COPYING.LESSER` version 3 and LGPL notices in its source.
The [official project page](https://www.chemshell.org/dl-find/) identifies a
maintainer and newer GPR work; a 2009 methods paper does not establish that the
project is abandoned. The separately offered GPR snapshot is explicitly dated
2021 and labeled unsupported development code; it is not verified in our Amber
build. Sander's source comments identify the integrated lineage as Py-ChemShell
21.0.1, which must not be confused with every current upstream capability.

| Option | Relevant capability | Fit to this project |
|---|---|---|
| Native Sander/DL-FIND | P-RFO, numerical Hessian, Bofill; already linked | Reference implementation with fewest integration changes |
| Sella | Iterative Hessian diagonalization, internal-coordinate restricted-step P-RFO, constraints | Strongest newer candidate for reducing gradient calls; ASE/Sander QM/MM interface must be verified |
| pysisyphus | Restricted-step P-RFO, mode selection, numerical Hessians and IRC workflows | Attractive reaction-path orchestration; no direct Sander QM/MM integration was verified |
| geomeTRIC | Internal-coordinate TS optimization, initial/final numerical Hessians, Bofill | Another modern option through ASE/custom engine; not a demonstrated improvement for this enzyme |

Sella has a 2022 molecular benchmark, GPL-family free licensing (LGPL-3.0 in
its repository), and ongoing releases. Pysisyphus is GPL-3.0. The currently
inspected geomeTRIC license is BSD-derived with an additional prohibition on
using its code/binaries in AI training datasets; do not label it unmodified
BSD-3-Clause. All license statements concern inspected sources, not guarantees
about every packaged version. No alternative has been installed or adopted.

Recommendation: establish native DL-FIND as the reproducible baseline, then
compare Sella on identical small reactive models and representative QM/MM TS
guesses. Measure gradient calls, wall time, convergence frequency, mode identity
and endpoint connectivity under identical Hamiltonian, masks and tolerances.
Published small-molecule benchmarks do not establish superiority for Nyl12.
ASE has an Amber calculator, but DFTB3 QM/MM force correctness, sign/units,
coordinate mapping, link atoms, electrostatic embedding and frozen atoms need
explicit interface tests before using it with Sella. Adding Sella/ASE would
require the project's normal production-dependency confirmation.

References: [DL-FIND](https://doi.org/10.1021/jp9028968),
[Sella](https://doi.org/10.1021/acs.jctc.2c00395),
[Sella repository/license](https://github.com/zadorlab/sella),
[pysisyphus](https://doi.org/10.1002/qua.26390),
[geomeTRIC TS documentation](https://geometric.readthedocs.io/en/latest/transition.html),
[ASE Amber interface](https://docs.ase-lib.org/ase/calculators/amber.html).

## AF3 contact restraints and four equivalent active sites

User supplied and implementation inspected:
`/data/jpc/code/alphafold3-contact-restraints`, HEAD
`ef5dc91d65b4845f63e95ea03316f1613c6dfb5e` (clean worktree when inspected).
The user reports successful tests and production runs. Existing source/tests
and historical campaign notes corroborate the capability. It is available as
an explicit guided-prediction arm; it must not be rediscovered or reimplemented.

Production evidence was located at
`/data/jpc/enzyme-substrate-modeling/output_nyl12_af3_contact_pa6_l4z_l2_four_site_static_20261009/`.
RUN_NOTES and manifests identify the same fork commit, WT plus eight variants,
ten seeds each and four ligands per model. Pre/post-relax SUMMARY.json files
contain 90 models and 360 sites, with 250/193 Tier A/B sites respectively.
Original inference guidance diagnostics pass for all 90 models. Post-relax
reports retain those inference diagnostics; they are not a new measurement
that every imposed contact remained satisfied after relaxation.

Historical analyzer:
`/data/jpc/enzyme-substrate-modeling/scripts/analyze_nyl12_af3_contact_four_site.py`.
Shared site definitions are imported from
`analyze_unrestrained_polyamide_four_site_geometry.py`:

| Site | Alpha | Beta | Historical four-ligand campaign chain |
|---|---|---|---|
| AE | A | E | I |
| BF | B | F | J |
| CG | C | G | K |
| DH | D | H | L |

The historical analyzer assumes the four-ligand assignment above and PA6
reactive atom names C12/O4/N3/C11. It expects 12 guidance diagnostics (three
contacts per site). These are not the C54/O49/N59/C62 names in the current PA6
TI1 workflow; graph mapping is required before reusing selections.

For a single unrestrained ligand, chain I can occupy any of the four equivalent
sites. Score every ligand–site pairing using an identical predeclared geometric
criterion, retain all site results and the chosen correspondence, and allow a
productive pose at any equivalent site. Do not require ligand I at AE or select
the best atom/register independently at each site. For multiple ligands, use
an explicit nonduplicating assignment and retain ambiguities. Report local
symmetry-aligned core geometry, all four site outcomes and model-level any-site
success; do not let tail RMSD or arbitrary chain labels cause rejection.
Sites within a model are nested observations, not four independent replicates.
Use matched seed/model-level statistics for WT and designs. Existing single-AE
TI1 metric code is not yet a symmetry-aware ES evaluator.

Keep guided and unguided AF3 arms separate in manifests. Restraint satisfaction
documents compliance with imposed geometry, not spontaneous binding or
independent validation. The preferred independent ES arm remains unrestrained.
