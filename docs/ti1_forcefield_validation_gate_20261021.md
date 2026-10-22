# TI1 force-field validation prerequisite — 2026-10-21

The user requires the following order before further AF3 refinement:

1. Generate and validate the force field for the intended covalent TI1 state.
2. Relax the hydrogenated AF3 enzyme/substrate/water model classically with that
   validated parameter set and documented restraints.
3. Optimize the classically relaxed structure with DFTB3/MM.

This decision supersedes using the present provisional reactive topology for
general classical TI1 relaxation. No further AF3 refinement calculation should
start until a parameter-validation record supports the exact topology and state.
Preparation-only exports and read-only inventory audits remain appropriate.

## What exists and what is missing

The current topology declares the TI1 bond graph and proton inventory, but its
reactive bonded/charge terms have not been validated. The OG1–C54 topology bond
has k=390.8 kcal/mol/A^2 and req=1.3584 A; those numbers are not a fitted TI1
reference. The source includes explicitly labeled junction analogs and
placeholder terms. Arbitrary charge balancing is not a RESP fit. Completed
Gaussian calculations concern neutral free-NH2 Thr-NME, not the covalent
NH3+/oxyanion TI1 unit; they do not validate this parameterization.

Required work is to define the coupled TI1 chemical graph, protonation, atom
mapping, stereochemistry and representative capped model(s); obtain suitable
QM reference geometries/frequencies and ESP; fit constrained RESP charges and
derive/check the needed reactive-center bonded terms. Use the previously
approved wB97XD/def2-TZVP geometry protocol without extra dispersion and separate
HF/6-31G(d) ESP unless an explicit scientific decision changes them. A free
analogue need not support an enzyme-stabilized TI1 minimum: document any
environment, constraints or restrained reference, and do not label it a free
minimum. Do not silently change protonation or fit through chemical collapse.

Validate parameterized structures against independent QM geometry/curvature and
relevant conformational scans, including oxyanion/amide geometry and the new
Thr–substrate junction. Check charge totals and cap-transfer constraints, missing
or guessed atom types/parameters, stereochemistry, connectivity and torsions.
Define numerical acceptance tolerances before fitting and record reference
inputs/methods/hashes, fit residuals, held-out checks, software versions and the
final topology digest. A classical covalent bond remaining intact is expected
from its bonded topology and is not evidence that TI1 is a QM minimum.

The water-only TIP3P preconditioning already performed fixed every nonwater
atom. It cannot validate the provisional TI1 parameters. Existing ff19SB and
TIP3P environment parameters can be retained; compatibility and transfer of the
custom reactive region still require the checks above. Full-octamer neutral
nonreacting beta-NH2 termini require their separate pending parameter validation.

## Interpretation of the failed QM/MM diagnostic

Both OG1 and C54 were QM atoms. Installed Sander removes their classical bond
term, so DFTB3 governs their interaction. The bond elongated from1.335 to2.195 A
and C54–O49 shortened to1.270 A. This is consistent with movement toward a
reactant-like structure from a suboptimal initial geometry, as the user
suggested. It does not establish a validated ES minimum, rule out a true TI1,
or establish the accuracy of DFTB3 in this QM/MM model.

Even after classical validation, test QM/MM refinement for unintended chemical
changes. If temporary state-preserving restraints are needed to prepare TI1,
declare them explicitly and assess their release. A restrained geometry must
not be presented as an unrestrained stationary point.
