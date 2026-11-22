# 45-model Boltz/Amber validation comparison

This report summarizes the retained raw and relaxed structures for eight design variants and WT (5 seeds each). Each model contains four assigned PA6 sites. A seed is the model unit; the four sites are nested observations.

## Per-variant counts

Counts are site observations out of 20 per variant and stage; ‘models with any joint’ is out of 5. ‘Force-converged’ applies to relaxed outputs only.

| Variant | Stage | Tier A | Tier A/B | Native pass | Joint | Models with any joint | Force conv. models | Conv. sites | Conv. joint sites |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| WT | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| WT | relaxed | 2/20 | 5/20 | 9/20 | 3/20 | 2/5 | 5/5 | 20/20 | 3/20 |
| variant_01 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_01 | relaxed | 2/20 | 5/20 | 12/20 | 4/20 | 3/5 | 5/5 | 20/20 | 4/20 |
| variant_02 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_02 | relaxed | 1/20 | 2/20 | 6/20 | 1/20 | 1/5 | 5/5 | 20/20 | 1/20 |
| variant_03 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_03 | relaxed | 11/20 | 15/20 | 8/20 | 7/20 | 4/5 | 5/5 | 20/20 | 7/20 |
| variant_04 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_04 | relaxed | 4/20 | 12/20 | 14/20 | 7/20 | 4/5 | 5/5 | 20/20 | 7/20 |
| variant_05 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_05 | relaxed | 4/20 | 8/20 | 13/20 | 3/20 | 2/5 | 4/5 | 16/20 | 1/20 |
| variant_06 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_06 | relaxed | 4/20 | 10/20 | 18/20 | 8/20 | 3/5 | 5/5 | 20/20 | 8/20 |
| variant_07 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_07 | relaxed | 9/20 | 15/20 | 9/20 | 7/20 | 2/5 | 5/5 | 20/20 | 7/20 |
| variant_08 | raw | 0/20 | 0/20 | 20/20 | 0/20 | 0/5 | — | — | — |
| variant_08 | relaxed | 1/20 | 4/20 | 12/20 | 0/20 | 0/5 | 5/5 | 20/20 | 0/20 |

## Findings and limits

- Raw assigned sites: 0/180 Tier A/B, 180/180 native pass, and 0/180 joint. Relaxed: 76/180 Tier A/B, 101/180 native pass, 40/180 joint; 21/45 models have any joint site.
- Force convergence: 44/45 models (176/180 sites). Of 40 joint sites overall, 38/176 sites in converged models jointly pass. The nonconverged model is variant_05 seed 2026112204 (final RMS gradient 0.1455 kcal mol⁻¹ Å⁻¹); it has two joint sites retained in all-model diagnostics, outside the force-converged set.
- Assigned ligand-to-site pairings matched the intended I:AE, J:BF, K:CG, L:DH mapping for all models at both stages; no swaps were observed.
- See `comparison_summary.json` for raw/relaxed d1, d2, d3, Burgi-Dunitz angle and signed phi distributions; individual Tier B/native gate-failure counts; joint-failure component counts; source hashes; and up to three unique force-converged joint-passing models with same-seed WT path/gradient controls.
- The A/B and native gates are inherited unchanged from `scripts/screen_pa6_binding_orientation.py` and its Tier function source; no thresholds were tuned here. Missing/nonfinite values fail the audit.
- Historical 900-model PA6 calibration records are in `analysis/nyl12_filter_calibration/`; WT, K59D and V111C context is summarized in the JSON. Cohorts and protocols differ, so no pooling or recalibration was performed.
- These are static ES geometries after prediction/minimization, not catalytic rates, transition states, or barriers. Model seeds are not independent experimental replicates.

## Provenance

Exact SHA-256 digests for source tables, status/manifest, runtime, and geometry implementation are in `comparison_summary.json`. Boltz 2.2.1 was split across violet and tulip; Amber relaxation ran on geranium. Host assignment counts by variant are retained in the JSON. The campaign records that same-number seeds do not imply numerically identical cross-build realizations, so seed matching is a design label rather than proof of paired stochastic draws. The campaign completed 45 normalized and 45 relaxed models; retained tables contain 90 model-stage rows, 360 assigned sites, 1,440 ligand-site comparisons, and 45 relaxation acceptance records.
