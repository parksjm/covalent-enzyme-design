# DIAGNOSTIC ONLY — NOT FOR PUBLICATION

Curiosity-only omission of one accepted observation; production MBAR is unchanged.
Full refit: all native bias terms, 81,999 observations, 100 states; decrement only state 85 from 100 to 99.
Fixed original weights: delete the same observation and renormalize; no MBAR refit.
Both use the original 0.025-wide s bins and zero at s=0.0125. No new uncertainty or accepted barrier.

| Estimate | F at s=0.2625 | ESS there | Maximum frame share | Whole-profile maximum |
|---|---:|---:|---:|---:|
| baseline | 2.899793870 | 1.066235 | 96.83889% | 26.280034298 at s=0.6875 |
| full_refit | 4.959096384 | 9.768901 | 26.96487% | 26.280000779 at s=0.6875 |
| fixed_weights | 4.959081723 | 9.771817 | 26.95855% | 26.280034298 at s=0.6875 |

Energies are kcal/mol. Full authoritative evidence: comparison.json and omitted_observation.json.
Fixed-weight analytical local shift: 2.059287852 kcal/mol.
Full-refit maximum common-bin change: 2.059302514 kcal/mol; outside affected bin: 0.000629840 kcal/mol.
Gauge residual: 2.13e-14 kBT; overlap components: 1.
New dominant observation in affected bin: 78214 (activation_grid33/activation_support_3_04).

This diagnoses sensitivity to one high-weight accepted observation. Its removal is not justified as data cleaning.
Neither estimator establishes sampling convergence, a physical transition state or a validated activation barrier.
