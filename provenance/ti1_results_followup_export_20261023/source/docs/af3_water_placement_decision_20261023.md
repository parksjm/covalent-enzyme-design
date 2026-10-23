# AF3 relay-water placement: status and priority — 2026-10-23

## Decision

Continue TI1 parameter development and physical local refinement using the
existing AF3 models. Improving the user's contact-restraint implementation
remains useful parallel work, but is not a prerequisite for the current
parameterization stage. Leave `/data/jpc/code/alphafold3-contact-restraints`
unchanged here. No new AF3 production/test run is requested by this decision.

This is a route to test, not a claim that physical refinement has already
produced an acceptable water bridge or that the AF3 implementation is correct
in every respect. The user reports limited progress on the separate AF3 code
project; the decision does not depend on that project being finished first.

## What the existing tests established

The minimal system contained one alpha chain, one beta chain, covalently bound
PA6 L4 and one water oxygen. Actual AF3 feature preflight retained water and
both Thr233 N/OG1-to-water contacts. Inference completed in both attempts.

| Attempt | Observed geometry | Interpretation |
|---|---|---|
| 01: two contact restraints | N–Ow 3.051 A; OG1–Ow 2.610 A; water–substrate C62 1.673 A | Requested 2.6–3.2 A distances met, but severe overlap; not an acceptable water placement |
| 02: extra exclusion restraints | N–Ow 4.847 A; OG1–Ow 3.172 A; residual clashes | Lost the N contact; not an accepted improvement |

The 28 extra exclusions in attempt02 are not adopted as the refinement
protocol. These two small tests do not establish either reliable AF3 water
placement or the general impossibility of obtaining it. An oxygen-only
prediction does not establish water-H orientation or catalytic protonation.

Independent source records:

- `output_nyl12_af3_water_audit_20261021/pa6_l4_alpha_beta/attempt_01/result/report.json`
- `output_nyl12_af3_water_audit_20261021/pa6_l4_alpha_beta/attempt_02/result/report.json`

## What physical preparation established

Explicit hydrogens were constructed while retaining the original AF3 oxygen.
Water uses WAT/O,H1,H2 and the audited TIP3P charges, types and Lennard–Jones
parameters. Names alone do not assign parameters. The two water hydrogens were
present in the diagnostic trajectory; corrected PyMOL viewers now display
them and have been checked in actual PyMOL.

User-approved water-only Amber/TIP3P minimization moved the oxygen 1.408215 A
while nonwater displacement was 0.0 A. It exhausted its cycle budget above
the convergence target, so it was coordinate preconditioning, not an accepted
minimum. The later DFTB3/MM endpoint passed the water steric screen but failed
directional bridge geometry and TI1 connectivity geometry: Thr OG1–C54
lengthened to 2.195122 A, and the optimization gradient target was not met.
Later refinement stages were withheld. No successful water repair is claimed.

Both reactive atoms were QM atoms, so their interaction followed the DFTB3
Hamiltonian, not the provisional classical topology's bond term. This failure
does not demonstrate that better AF3 placement is required, prove a stable ES
minimum, or rule out a genuine enzyme-stabilized TI1. The poor starting state,
reference/parameter quality and QM/MM method/setup remain relevant concerns.

Retained evidence: `docs/local_repair_checkpoint_20261021.md` and
`output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/`.

## Next experimental sequence in the computational workflow

1. Complete geometry, cap-artifact, electrostatic, bonded/torsional and
   native-state transfer validation for the intended TI1 state. A completed
   capped-cluster Gaussian calculation is only a prerequisite, not a validated
   force field. Cap sensitivity must not be hidden by averaging references.
2. Assemble the exact accepted topology and add/check all hydrogens. Perform
   scripted restrained classical relaxation with the water and defined local
   environment mobile, preserving source models and explicit restraint records.
3. Assess stereochemistry, oxyanion contacts, proton inventory, clashes and
   the directional N–H...Ow–H...OG1 bridge. Reduce/release the bridge restraints
   in a separate assessment; satisfying imposed bounds is not independent proof.
4. Attempt DFTB3/MM only from a chemically acceptable prepared state. Retain
   trajectories and test TI1 preservation and convergence again. Declared
   temporary state-preserving restraints cannot certify an unrestrained minimum.
5. Revisit AF3 development if validated physical preparation still cannot give
   a usable starting state, or if repeated design-scale preparation makes the
   initial oxygen placement a demonstrated bottleneck.

Reference parameter development may proceed now. Protein relaxation remains
gated by `docs/ti1_forcefield_validation_gate_20261021.md`. Changes to AF3
guidance should later be tested with retained actual features, seed-matched
predictions and independent geometry/steric assessment, rather than success
being inferred from the input JSON alone.

All actions remain in `output_agent_actions/20261021.jsonl`; calculation inputs,
source snapshots, process status and scientific acceptance are separate records.
