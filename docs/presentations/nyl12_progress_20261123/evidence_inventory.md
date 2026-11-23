# Nyl12 progress deck evidence inventory

Date: 2026-11-23. This inventory is a source guide for the deck author. The
JSON file is authoritative for exact absolute paths, values, interpretation and
caveats. Existing source artifacts were not copied or modified.

## Suggested slide assets

- **PES:** `vmd_final_24img/final_24img/energy_profile.png` and
  `reactive_distances.png` in the 2026-11-13 NEB export. Label this as a PES;
  the 35.80895 kcal/mol NEB rise is not a validated barrier. The separate TS
  candidate is 46.41135 kcal/mol above independently relaxed E-S, meets the
  gradient criterion, but has no verified first-order mode or IRC.
- **PMF:** `compact_export/primary_q_acyl_pmf.png` and `mbar_overlap.png` from
  the 2026-11-22 PMF refresh. The apparent maximum is 47.0475 kcal/mol versus
  46.9414 at the same E-S density bin. It has 63 completed windows, 30 ps per
  window, target ESS 7.29 and three drift flags. This PMF is not converged.
- **Matched WT MD:** `analysis_after_md/plots/tier_occupancy.png` and
  `occupancy/occupancy.png`. Four matched 100 ns replicas per state. Tier A OR B occupancy is 27.77%
  (NH2/Asp272/Asp274), 5.18% (NH2/Asp272/ASH274), and 28.99%
  (NH3+/Asp272/Asp274); native occupancies are 0.0175%, 0.245%, and 0%.
  These are structural occupancies, not activity or pKa results.
- **Design:** `validation/analysis/plots/tier_native_comparison.png` and
  `geometry_before_after.png`. Counts (broader 12/80, focused 8/80, WT 1/20)
  are exploratory structural counts from selected sequences/five seeds, not
  activity predictions.
- **Parameterization:** current PA6 L4 manuscript figures include modular
  architecture, eight bonded profiles and 500 ps MD stability, plus a recent
  independent charge comparison at
  `/data/jpc/polyamide-amber-ff/manuscript/independent_charge_comparison_20261123/independent_charge_comparison.png`.
  Current release scope is zwitterionic PA6 L4 only. The compact-conformer ESP
  holdout fails; arbitrary lengths, bulk polymer and covalent catalytic states
  remain unvalidated.
- **VMD:** `umbrella_sampling/vmd_dynamic/view_umbrella.vmd` and `path.dcd` are
  available, but this inventory found no GUI-rendered screenshot. A placeholder
  should say so plainly if an image is not rendered. The concatenated biased
  windows are not a continuous reaction trajectory.

## pKa readout for the deck

The latest MCCE outputs for mature 5XYG and uncleaved controls 7YU0/7YU1
finished natively, but the qualification receipts explicitly say
`native_completed_pka_unqualified`, `assigned_states=false`, and
`sampling_convergence=false`. They used six Monte Carlo runs × 2,000 iterations,
protein-only inputs, level-1 conformers and dielectric 4/80; a known PRO charge
template inconsistency remains unresolved.

| Structure | Current MCCE target labels | How to describe them |
| --- | --- | --- |
| Mature 5XYG | Thr267 NTR C/D >14; Asp306 C/D <0; Asp308 C/D <0 | Bounds outside pH 0–14, not exact pKas. Diagnostic, no states assigned. |
| Uncleaved 7YU0 | Asp306 A <0, B 2.771; Asp308 A 10.085, B <0 | Finite fits are poor: 1000χ² 47.926 and 58.469. |
| Uncleaved 7YU1 | Asp306 A/B <0; Asp308 A 1.195, B 7.883 | Asp308 fits are poor: 1000χ² 64.212 and 993.765. |

Corrected PROPKA mature 5XYG estimated N+267 at 6.35–6.55 (default) or
4.59–4.60 (coupled alternative), Asp306 at 2.19–2.21 and Asp308 at 5.49–6.49
(default) or 7.24–8.45 (coupled). The older precursor PROPKA comparison is not
matched to the repaired MCCE structures. The precursor 7YU0/7YU1 carry a
Cys267 mutation and do not have a mature Thr267 N terminus. Present both methods
as unresolved, structure-sensitive evidence; do not pick a protonation state.

## Eyring table caveat

No validated activation free energy is available here. Do not use the NEB PES rise or the unconverged PMF maximum as ΔG‡ for predicted Eyring rates; if shown for teaching, label the calculation as illustrative and not a rate prediction.

## Existing PA6 force-field deck

`/data/jpc/polyamide-amber-ff/output/pa6_l4_validation_20261029/presentation_checkpoints/PA6_L4_VALIDATION_DECISION_20261031_ACTIVE_MD_CHECKPOINT.pptx`

The October 31 checkpoint contains useful historical parameterization pipeline,
charge and bonded-scan visuals. Its old withheld-production disposition is
superseded by the current scoped PA6 L4 release README and manifest. Current
release: standard two-stage multiconformer RESP/GAFF2, three improper
overrides, eight bonded-profile gates and three accepted 500 ps L4 replicas;
limited to zwitterionic L4.
