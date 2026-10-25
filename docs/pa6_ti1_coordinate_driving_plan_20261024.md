# Next focused task: enzyme-supported TI1 candidate

The user authorized preparing coordinate-driving inputs on 2026-10-24.
This is a proposed construction pilot, not a completed reaction pathway.
The restrained ES reference came from WT replica 1, full frame 651 at
32.55 ns production time (legacy CSV label 32.50 ns). It was the eligible frame closest to the midpoint
of the predefined 20–45 ns window, not an end-of-run snapshot. Provenance:
`output_nyl12_pa6_es_audit_20261024/reuse_audit_attempt1/reuse_audit.json`.
The full-solvent frame was extracted and imaged around the substrate without
discarding waters, ions or hydrogens. Its original topology/protonation was
retained. A 25-cycle restrained DFTB3/MM minimization was followed by a
107-cycle continuation to the 0.05 kcal/mol/Å RMS target (final 0.042903),
using the same original MM positional reference. No late-MD endpoint or
AF3 model was substituted for this seed.

Use the converged restrained full-system DFTB3/MM ES endpoint and its unchanged
QM region. Generate a TI1 candidate through short constrained optimization
steps, using two explicitly recorded coordinates rather than relying on a
single combination that can conceal an incorrect proton or bond state:

- Attack: distance Thr233 OG1–substrate C54.
- Direct proton-transfer hypothesis: distance OG1–HG1 minus distance
  N-terminal Thr233 N–HG1, following the same original hydroxyl hydrogen.

Couple the schedules, relax the remaining coordinates at each step, and retain
the original MM positional reference. Exact schedules, force constants and
convergence criteria must be declared before launch. Do not assign an energy
barrier from restraint-biased energies or label a path maximum a verified TS.

Monitor substrate C54–O49 and C54–N59, tetrahedral geometry, Thr proton
connectivity, both Asp residues, oxyanion donor contacts, d1/d2/d3, a_BD and
herschlag_phi. The candidate target is covalent tetrahedral geometry with
Thr NH3+ and an oxyanion; that chemical assignment requires coordinate review.
The direct proton route is a hypothesis, not an established Nyl12 mechanism.
Do not impose a water relay. Unexpected transfer toward an Asp should prompt
chemical review rather than stronger restraints to force the chosen route.

At a candidate endpoint, reduce and release the reaction-coordinate restraints
while retaining the original MM positional reference. Persistence or collapse
tests local stability in this model; neither establishes the complete mechanism.
Only then consider selected hybrid-DFT checks and native TI1 parameterization.
No driving calculation has been launched at this checkpoint.

The user subsequently authorized one sequential forward/reverse overnight
construction pilot. Proposed input specification:
`configs/pa6_ti1_coordinate_driving_20261024.json`. Use 18 forward targets and
17 reverse targets, with reverse execution conditional on a geometrically
plausible forward TI1 endpoint. Stop at the first failed or unconverged window.
Save full-system optimizer coordinates every ten output intervals and exact
endpoints; retain the original MM restraint reference for both directions.

AF3 is not required to construct this WT TI1 candidate. It remains useful for
redesigned-sequence and alternative-pose models. Short QM/MM MD does not
require classical TI1 bonded parameters. Holding QM coordinates fixed while
evaluating QM/MM forces on mobile MM atoms is another possible relaxation
stage. Simply freezing atoms in an ordinary MM simulation does not remove
the need for appropriate TI1 electrostatics and nonbonded interactions.
Equilibration and parameter fitting are deferred from tonight's path pilot.

## Subsequent TS1 refinement

Refine a chemically appropriate guess near the local ES-to-TI1 path maximum
with the installed Sander/DL-FIND interface: P-RFO with a declared movable
active region, a numerical initial Hessian from analytic QM/MM gradients and
Bofill updates. Remove driving restraints for the saddle search; distinguish
any retained environmental restraints from an unrestrained PES. Do not form a
full Hessian for all 103251 atoms. Existing implementation evidence and
limitations: `docs/ts_search_and_symmetry_decisions_20261020.md`.

Require a fresh active-space Hessian with one significant chemically relevant
negative mode and downhill connectivity toward ES and TI1, plus adequate
stationarity of the declared movable coordinates. A restricted active-space
saddle is not proof of a full-system saddle. Refine each local elementary-step
barrier independently; there is no need to wait for the complete catalytic
cycle or assume a single highest point covers all transition states.

Reference: [DL-FIND transition-state algorithms, official ChemShell tutorial](https://www.chemshell.org/wp-content/static_files/py-chemshell/tutorial/build/html/transition_states.html).

Save actual optimizer snapshots and exact initial/final coordinates for every
window. Provide stripped VMD trajectories with a fixed nearby-water selection,
explicit source/frame mapping, retained hydrogens, and QM atoms in Licorice.
