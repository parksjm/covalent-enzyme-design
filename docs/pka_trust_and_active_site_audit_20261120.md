# pKa trust and active-site audit — 2026-11-20

## Decision

The completed Nyl12 MCCE pilot is **not qualified for quantitative catalytic
protonation assignments**. Preserve its outputs as a diagnostic result.
Successful native execution and correct terminal identification do not establish
pKa accuracy. Neither this audit nor PROPKA changes any MD/QM proton inventory.

## What passed, and what remains uncertain

The corrected MCCE input retains all 9,928 original heavy atoms and eight true
N/C terminal pairs without moving heavy atoms. Its starting-coordinate SHA-256
is identical to the newly HIS-normalized PROPKA input:
`afbca64d1e9277bc54ec13b9fa1fd7be7ff746c91237dcb54912b507ebd301c2`.
Both default and coupled PROPKA outputs use this input. Historical HIE-named
PROPKA results are a different preparation and are not substituted here.

The four MD active sites have the same sequence, but different instantaneous
geometries. Relative to alpha A/beta E, CA RMSDs for B/F, C/G and D/H are
1.514, 1.490 and 1.306 Å; corresponding reference-pocket heavy-atom RMSDs are
0.951, 1.072 and 1.094 Å. The pocket comprises 176 matched atoms within 8 Å
of E:Thr1 N/OG1 or E:Asp40/Asp42 CG after joint alpha/beta CA alignment.
Thus identical chain pKas are not required, but the observed shifts still
require validation.

MCCE used protein-only coordinates, omitted substrate/waters/ions, one
dielectric pair (4/80), and level-1 conformers. The latter explicitly disables
heavy-atom rotations and heavy relaxation cycles. Most target transitions are
not sampled within pH 0–14. More Monte Carlo steps alone cannot restore absent
conformations, missing chemical context or an unresolved titration midpoint.

## Matched-coordinate core-site comparison

| Group | PROPKA default, E/F/G/H | PROPKA coupled alternative, E/F/G/H | Native MCCE labels, E/F/G/H |
| --- | --- | --- | --- |
| Thr233 N terminus | 9.70 / 9.31 / 7.06 / 9.12 | Same | 6.513 / <0 / <0 / <0 |
| Asp272 | 1.02 / 1.39 / 1.99 / 1.36 | 2.24 / 3.88 / 3.52 / 2.97 | <0 / 1.326 / <0 / <0 |
| Asp274 | 2.95 / 4.37 / 3.94 / 3.03 | 1.72 / 1.88 / 2.41 / 1.42 | <0 / <0 / <0 / <0 |

The native NTR-E midpoint agrees with its occupancy curve (approximately 6.53).
The F-chain Asp272 fit is poorly described by a single sigmoid: fitted pKa
1.326, slope 0.176 and 1000×chi² 43.935, versus a direct interpolated 50%
crossing near 0.74. The crossing is descriptive, not a replacement accepted pKa.
Most other core curves are saturated; their bounds do not resolve exact values.
Small nonmonotonic changes and coupled-site alternative assignments are retained.

Native `mc_out` contains six-run target occupancy statistics. For example,
NTR-E positive-state occupancy is 0.297 ± 0.040 at pH 7 and 0.059 ± 0.015 at
pH 8 (sample SD across runs). These are not pKa confidence intervals and do
not measure structural/model error.

## A concrete output inconsistency and parameter concern

For His257 in beta E, native `pK.out` reports >14, although its summed positive
conformer occupancy is only 0–0.040 over the entire pH grid. That label cannot
be interpreted as a valid high-pKa bound. Its curve is unresolved/nonmonotonic;
the native single-sigmoid fit/output must be checked against occupancy data.
The F-chain His257 fit is also not a clean midpoint. Do not propagate those
fit labels into catalytic microstate claims.

The pinned upstream MCCE source also contains inconsistent proline charge
templates: `param/pro.ftpl` explicitly makes the PRO backbone uncharged, while
`param04/pro.tpl` gives N/CA/C/O charges −0.10/+0.10/+0.55/−0.55.
The observed zeros match the FTPL source, so the previous sparse-charge audit
established provenance, not electrostatic validity. This discrepancy is not
yet proven to cause the anomalous pKas. No parameters were patched.
Sources: [pinned FTPL](https://github.com/GunnerLab/MCCE4-Alpha/blob/08066fd303647a1287296bfae702bf4291e3baab/param/pro.ftpl)
and [bundled legacy template](https://github.com/GunnerLab/MCCE4-Alpha/blob/08066fd303647a1287296bfae702bf4291e3baab/param04/pro.tpl).

## Other nearby titratable residues

The screen includes titratable residues in any chain whose heavy atoms approach
within 6 Å of the union of Thr233/Asp272/Asp274, with explicit checks of the
paired alpha-chain Lys156, Phe113 and Asn186. This geometric neighborhood is
not evidence of a connected proton-transfer pathway.

| Nyl12 residue | PROPKA default range across copies | Interpretation limited to this static screen |
| --- | ---: | --- |
| Lys156, alpha chains A–D | 8.46–9.91 | Depressed relative to model 10.5; a neutral population near pH 8 is worth testing, but MCCE instead retains positive Lys at pH 7–8. No neutral-Lys assignment is established |
| His257, beta-local 25 | 4.85–5.22 | Mostly neutral near pH 8 in this predictor; tautomer and hydrogen-bond orientation may matter. MCCE fit labels are unreliable here |
| Tyr267, beta-local 35 | 14.00–14.57 | Elevated relative to model 10.0; both methods retain neutral phenol near pH 8 in this preparation |
| Asp45, alpha chains A–D | −1.97 to −0.06 | Strong downward shift, favoring anionic Asp in the static predictor; a possible electrostatic-network residue, not an established acid donor |
| Thr233 N terminus | 7.06–9.70 | Conformation-dependent terminal state; the MCCE discrepancy prevents consensus |

Phe113 and Asn186 are not ordinary titratable groups in these pKa models.
Their oxyanion-hole hydrogen bonding remains important independently of a
side-chain pKa screen. PA6 omission and single-snapshot bias limit every row.
The Lys189 neutral-amine proposal for NylC must not be assumed to establish
neutral Lys156 in Nyl12.

## NylC status and repair

Historical MCCE steps 1–4 completed only for mature 5XYG. That result was
rejected because missing-coordinate fragments created artificial termini.
Uncleaved 7YU0/7YU1 inputs were staged and held; there are no precursor MCCE
results. PROPKA was run on all three deposited sequences, without backmutation.
Even the corrected mature PROPKA input later failed the stricter full-terminal
inventory check, so those crystal comparisons retain preparation caveats.

The PI now authorized PDBFixer installation and 5XYG repair. PDBFixer 1.12.0
(source commit `94cfa4c0ca551cdc5f13320f9a658efd59f2b881`) and OpenMM 8.6.1
are installed in the isolated uv project `tools/pka-preparation`.
The repair must split mature alpha 1–266 from beta 267–355 **before** filling
SEQRES gaps, preserve the deposited sequence/observed heavy coordinates and
record every modeled residue. It must never create a peptide bond across
autocatalytic cleavage. Rebuilt missing loops/tails are modeled, not observed
X-ray coordinates. Read the repair's status/acceptance before using its output.

## Evidence and next gates

Canonical analysis directory:
`analysis/nyl12_mcce_step1_preflight_20261116/canonical_atom_names_repair_20261120/titration_audit_20261120/`.
It contains source-hashed matched comparisons, full charge-state curves,
nearby-site data, NylC history and an explicit trust verdict. The reusable
plot script is `scripts/plot_mcce_titration_audit.py`; raw MCCE files remain
unchanged in the parent `prepared/` directory.

Finish and audit 5XYG repair; then repeat matched protein-only predictors on
the repaired input with correct physical termini. Before promoting MCCE
numbers, resolve the charge-template discrepancy against an authoritative
parameter set and validate fitting/occupancy consistency. Subsequent comparisons
should vary one factor at a time: conformer sampling, dielectric and structural
ensemble. A single outlying predictor result does not authorize another
protonation-state simulation campaign.

## Native 5XYG reconstruction checkpoint

PDBFixer generated the mature split candidate but failed the acceptance gate:
new alpha-chain A9–A10 C–N was 1.967 Å. The candidate is not accepted for pKa
prediction. The original X-ray atoms and failed attempt are preserved.
One controlled geometry continuation using the installed PDBFixer soft potential
was authorized: optimize newly built atoms only, freeze all original heavy
coordinates, tighten minimizer tolerance and retain the original scientific
gates. This is structural reconstruction, not Amber production relaxation or
an accepted pKa run. Check its status/acceptance under
`analysis/negoro_2023_protonation/pdbfixer_repair_20261120/` before proceeding.

### Controlled reconstruction continuation launched

At 21:41 EST Nov20, `nylc_5xyg_softmin_continue_20261121` started the
new-atoms-only refinement from the retained failed candidate. Work directory:
`analysis/negoro_2023_protonation/pdbfixer_repair_20261120/repair_run_20261121/`.
Read `softmin_continuation.log` and `softmin_continuation_status.json`.
The separate `_softmin.pdb` requires its validation receipt before use.
No valid repaired structure or new NylC pKa result is claimed at launch.

## Three-structure rerun authorized — 2026-11-20 follow-up

The PI requested repair/rerun of mature 5XYG and the uncleaved 7YU0/7YU1
controls now. The earlier batch was held after the mature preparation failure;
it never ran precursor MCCE. Their deposited mutations remain unchanged, so
this is a structure/maturation/sequence comparison, not an isolated cleavage
free-energy difference.

Technical diagnosis found no insertion atom-copy bug: all 357 modeled identities
match and bad peptide bonds have their native equilibrium lengths. PDBFixer's
soft potential has weak angular restraints, bounded repulsion and a rescue
screen excluding same-residue contacts. It therefore did not yield acceptable
peptide/ring geometry here. A generic Amber19SB/GBn2 refinement now uses full
bonded/torsion/LJ interactions, freezes observed heavy atoms, and moves new atoms
and temporary hydrogens. Building hydrogens at pH8 is explicitly recorded;
they are removed before pKa prediction, and this is not a catalytic protonation
assignment. Mature and precursor structures have explicit distinct terminal
contracts; 267 is internal in both precursors.

Independent tmux controllers are being prepared for each structure, so failure
of one cannot silently prevent the others. They may proceed to native MCCE
only after coordinate, modeled-clash, bond, sequence, OXT and exact-terminal
checks pass. The pinned unchanged MCCE parameter set is retained for this
**diagnostic comparison**; completion will not promote pKas to quantitative
acceptance while the parameter/curve concerns remain unresolved. Read native
per-structure statuses in
`analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/`.
