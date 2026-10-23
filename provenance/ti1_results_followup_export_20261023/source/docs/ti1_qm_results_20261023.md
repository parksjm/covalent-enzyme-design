# Completed TI1 Gaussian references — 2026-10-23

## Findings and limits

All four 37-atom capped-anion references completed wB97XD/def2TZVP
optimization and separate frequency/Hessian calculations normally. They retain
the intended graph, stereochemistry and proton inventory, with no detected
oxyanion proton transfer. Frozen phi/psi/chi1 remain within the predeclared
0.10-degree tolerance. Each raw frequency report contains 105 positive modes.

The gradient/Hessian assessment also passes at both numerical derivative
steps: all 102 internal directions of the constrained Lagrangian Hessian have
positive curvature above the predeclared tolerance. This establishes local
stability **under the three frozen torsions**, not unconstrained enzyme TI1
stability, a transition state, or a validated force field.

| Attack / cap phi | OG1–C* (A) | C*–O− (A) | C*–N (A) | C*–acyl C (A) | Lowest raw frequency (cm−1) | Smallest constrained eigenvalue (Eh/bohr²) |
|---|---:|---:|---:|---:|---:|---:|
| R / −120 | 1.50616 | 1.30676 | 1.49657 | 1.53617 | 22.0614 | 4.4377e−5 |
| R / −60 | 1.55943 | 1.27510 | 1.49749 | 1.54614 | 25.2106 | 1.26825e−4 |
| S / −120 | 1.52256 | 1.29553 | 1.48034 | 1.54692 | 21.6190 | 1.26378e−4 |
| S / −60 | 1.50858 | 1.30295 | 1.47287 | 1.54720 | 27.4340 | 1.55771e−4 |

These are two attack configurations and two declared cap orientations, not
four independent biological samples. The model is
ACE-THR(-O-C*(O−)(CH3)(NHCH3))-NME, C10H20N3O4−, charge−1, singlet,
with 20 explicit H. It lacks the native catalytic NH3+, full PA6 L4 and relay
water. No native charge or parameter transfer follows from this table.

The four source calculation directories are
`output_nyl12_ti1_parameters_20261022/pa6_local_bonded/attack_<R_or_S>_cap_phi_<minus120_or_minus60>/`.
Within each, retain `optimization_01/scientific_acceptance.json`,
`frequency_01/scientific_acceptance.json` and
`conditional_curvature_01/assessment.json`, as well as their immutable inputs,
checkpoints, logs and process/source identities. The third record is required;
positive raw frequencies alone were not the acceptance rule.

## Cap sensitivity is a failed transfer prerequisite

Compare the two cap orientations within each attack configuration. For R,
OG1–C* changes by 0.05327 A and C*–O− by 0.03166 A. Both exceed the predeclared
0.03 A bond-length sensitivity limit. S passes that length screen, with a
maximum change of about 0.0140 A. This is not grounds to select S merely
because it passes or to average away the R failure.

The capped Thr backbone N–H has H...O− contacts of 1.6630 A (R/−120),
1.7218 A (S/−120) and 1.5817 A (S/−60); R/−60 has no comparably short
N–H...oxyanion contact. Frozen backbone torsions did not prevent these
intramolecular contacts. They are consistent with a gas-phase cap/backbone
bias, but do not isolate its cause or prove any one conformer is native-like.
The cap sensitivity gate remains failed regardless of a plausible-looking
trajectory or successful constrained optimization.

As specified in the prospective protocol, excessive cap sensitivity requires
a larger reference before native transfer. The already-planned native
Thr-NH3+/full-PA6-L4/oxyanion/NME/water graph has 103 atoms and total charge0;
it needs separately mapped and accepted 3D geometries. Its native N terminus
has no phi, and the source backbone must not be assigned an invented native
phi. Additional environment or constraints, if needed, require explicit
definitions; no solvent, protonation or reaction-state substitution is made
by this checkpoint. A lower small-cluster gas-phase energy is not an enzyme
attack-face selection rule.

## Work started after these results

- Independently audit completed geometry/energy records and export final
  XYZ/SDF coordinates; final exports replace neither original data nor the
  separately labeled interim trajectories.
- Quantify raw-Hessian bond-constant sensitivity with the installed, previously
  convention-tested AmberTools extractor. These are diagnostic estimates,
  not emitted parameter overrides. Keep each cap/attack result separate.
  Report both symmetric percent differences and the conservative max/min
  force-constant ratio limit1.20; do not hide denominator ambiguity.
  The first executable diagnostic stopped during AmberTools import because
  AMBERHOME was missing; a fresh explicit-environment attempt then reached
  R/−60 OG1–C* and rejected its extracted constant. Neither attempt emitted
  parameters. Convention checks on manufactured Hessians do not guarantee
  usable constants for these molecular constrained references; retain this
  additional method limitation rather than force a positive replacement.
- Prepare separate HF/6-31G(d) MK ESP calculations on these accepted geometries,
  with source checkpoint copies and separate output checkpoints. Their scope
  is capped-surrogate diagnostics and later surrogate baseline testing,
  **not** fitting/transferring native charges or certifying the TI1 force field.
- A strict standalone ESP checker was added because the older checker assumes
  a combined optimization/frequency/ESP job. The new checker tests unchanged
  geometry, electronic state/method, complete finite espgen records and atomic
  center agreement. Process completion and scientific acceptance stay separate.

The full protein TI1 relaxation gate remains closed. Resolve reference-model
sensitivity, prepare/validate native charge-fit conformers, perform constrained
RESP and baseline/targeted bonded validation, then test classical AF3-model
relaxation followed by DFTB3/MM. Held-out scan inputs are not generated by this
checkpoint. No reaction barrier, improved activity or accepted enzyme TI1
minimum is claimed.

## Related decisions and provenance

`docs/af3_water_placement_decision_20261023.md` consolidates the existing AF3
and physical water-preparation results. Better AF3 contact restraints can
remain parallel work; validated physical preparation is the current route to
test, not an already completed repair. The user's AF3 fork is unchanged.

Ranking remains defined in `docs/design_ranking_and_iteration_20261022.md`:
primary independent geometry recovery refers to unguided noncovalent ES;
TI1/TS have separate state-specific validation. Interface_delta is optional
pending demonstrated value. No design campaign has launched.

All new diagnostic calculations use new dated directories, source hashes,
explicit environments and tmux. Active journal:
`output_agent_actions/20261021.jsonl`. Command-approval timeouts and bwrap
failures are execution failures, not scientific results, and do not establish
that a proposed job started. Inspect retained process records before resuming.

## Completed follow-up diagnostics — 2026-10-23

All four standalone HF/6-31G(d) MK calculations completed and passed
the separate source/geometry/state/grid acceptance checks. Grid counts:

- `attack_R_cap_phi_minus120`: 8606 points.
- `attack_R_cap_phi_minus60`: 8843 points.
- `attack_S_cap_phi_minus120`: 8712 points.
- `attack_S_cap_phi_minus60`: 8897 points.

No RESP charges were fitted or transferred. The raw-Hessian probe for
R/−60 OG1–C* returned a finite signed Amber-form K of −5.88965426
kcal/mol/A², with zero imaginary component. The tuple convention
was checked; no absolute-value replacement or parameter emission
was made. A negative local block projection is not a failed
full-molecule constrained-curvature test. The cap sensitivity and
parameter-transfer gates remain closed; larger/native reference
geometry and baseline/targeted validation are the next scientific
requirements. No protein relaxation or AF3 fork change occurred.

Records: `output_nyl12_parameter_validation_20261023/surrogate_esp_acceptance_01/` and `bond_sensitivity_diagnostic_04/`.
