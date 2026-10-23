# TI1 held-out validation plan — 2026-10-22

This is a prospective validation plan, not a launched scan or accepted force
field. Complete the geometry, conditional-curvature and cap-sensitivity gates
in `ti1_parameter_protocol_20261022.md` before producing reference-specific
inputs. Keep both attack configurations separate until their comparison is
resolved. The actual native NH3+/oxyanion model requires its own charge and
transfer validation.

## Independent finite-displacement checks

First test the standard ff19SB/GAFF2 candidate as described in the parameter
protocol. If targeted bond/angle refitting is needed, derive candidate constants
from the accepted reference Hessian, not from the energy points below. Reserve
all listed finite-displacement energies for validation. If a point later informs a refit, label it training and generate
new independent validation points; do not continue calling it held out.

Use the same atom order, electronic state, DFT method/grid and reference
conformer as the accepted optimization. Preserve phi, psi and chi1 to the
existing 0.10-degree tolerance at every point. Generate each distortion from
the same reference, rather than accumulating coordinate edits. No point is
eligible if it changes stereochemistry, proton inventory or intended topology.

| Coordinate family | Initial check | Expanded check if initial passes |
| --- | --- | --- |
| Four bonds from C*: OG1, O−, N, acyl C | ±0.02 A, 8 points | ±0.04 A, 8 additional points |
| Six angles between those four C* neighbors | ±3 degrees, 12 points | ±6 degrees, 12 additional points |
| Thr CB–OG1–C* junction angle | ±3 degrees, 2 points | ±6 degrees, 2 additional points |

The initial panel has 22 distorted points plus the reference for one accepted
conformer. These offsets are project choices declared before observing their
energies; they are not literature-certified bounds. This panel specifically
checks local reactive-center transfer. It does not cover every bond, angle,
proper torsion or improper in the final force field.

For these local checks, evaluate DFT single-point energies on the exact
scripted coordinates without relaxing other atoms. Evaluate the complete
surrogate MM energy at those same coordinates, including its separately fitted
surrogate charge model, nonbonded terms, exclusions and Amber 1–4 scaling.
Subtract each method's own energy of the same reference geometry. Do not compare
an isolated fitted bond term to the total DFT energy, or a relaxed MM point to
an unrelaxed QM point. Do not use native-state charges on the capped anion.

Enforce the existing <=1 kcal/mol held-out relative-energy RMSE criterion in
the tested local basin, and also retain signed errors, maximum absolute error
and errors by coordinate family so an aggregate cannot conceal one bad term.
Finite displacements are a local energy check, not a reaction path or a barrier
calculation. If geometry/proton screening fails, retain the point as out of the
declared TI1 basin; do not fit it into a harmonic TI1 term or silently drop it
from the report. Local minima optimized with the candidate FF must separately
pass the existing bond/angle geometry tolerances.

## Torsions and remaining coverage

Retain standard ff19SB peptide-backbone terms. The three frozen reference
dihedrals do not supply independent data for refitting those terms. Explicitly
inventory all new atom types and every bond, angle, proper torsion and improper
touched by the covalent junction or tetrahedral/oxyanion atom-type changes;
missing terms must fail rather than inherit guessed defaults.

Two new-junction torsions requiring an explicit coverage decision are
CA–CB–OG1–C* and CB–OG1–C*–N. Their exact origin IDs and moving fragments must
be taken from the accepted atom map. Local energy checks cannot establish full
rotational barriers. Before fitting a torsion, declare its scan range, relaxed
coordinates, sampling and distinct held-out offsets in a separate manifest.
Retain the reference's three frozen torsions, monitor all relevant proton and
bond distances, and stop interpretation as a TI1 torsion scan if the graph or
protonation changes. No torsion-scan jobs have been launched at this checkpoint.

## Minimum records and acceptance

Each generated point must record source geometry and atom-map hashes, parent
conformer, coordinate definition, displacement, generation-code hash, intended
training/validation role, charge/multiplicity, frozen torsions, exact engine
input and geometry-screen result. Each calculation uses a fresh scripted tmux
attempt with explicit resources, source snapshots and process records.

Keep parameter-estimation success, engine completion, local validation,
cap-sensitivity acceptance and native-state transfer acceptance distinct.
Neither manufactured angle/bond checks nor the 22-point local panel alone
qualifies the complete TI1 force field for protein refinement.
