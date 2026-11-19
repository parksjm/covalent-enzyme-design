# Static pKa screens of Nyl12 and mature/precursor NylC

Status: exploratory, 2026-11-16. PROPKA 3.5.1 was already available in
`/data/jpc/enzyme-substrate-modeling/.pixi/envs/analysis`; no dependency was
installed or added to this repository. The reproducible tmux job and hashes
are in `qmmm/pa6_flipped_model30/09_propka_static_tier_a/job/` and
`screen_attempt1/manifest.json`. Raw `.pka` files and the complete four-site
table are in `screen_attempt1/`; `catalytic_pka.csv` is the compact table.

## Inputs and method

Three 300 K model30 MD frames were selected as the temporal median among
Tier A frames in one replica per Asp-pair microstate, using the identical
Tier A classifier in `08_asp_pair_md_comparison/analysis_attempt2/`. The
frames were ASP272/ASP274 replica 2 frame 284, ASH272/ASP274 replica 1
frame 4868, and ASP272/ASH274 replica 1 frame 3492. Each has positive
signed Herschlag phi (97.3°, 87.3°, 102.1°); OG1–C distances are
3.137, 3.131, and 3.184 Å, respectively. These are catalytically oriented
**noncovalent** snapshots, not reaction intermediates.

All eight mature protein chains were retained, with explicit chain breaks;
beta chains E–H have Thr233 as residue 1, Asp272 as residue 40, and Asp274
as residue 42. Protein hydrogens and the nonstandard PA6 ligand were omitted
for this first PROPKA screen. ASH heavy-atom records were written as ASP so
the same titratable-group model was applied to all three geometries. Thus the
three calculations compare **conformations**, not the free energies or
occupancies of specified protonation microstates. PROPKA's predicted pKa
cannot verify the input classical proton assignment.

## Results

The substrate-facing beta site E predictions are:

| MD starting microstate | Asp272 pKa | Asp274 pKa | Thr233 N-terminal pKa |
| --- | ---: | ---: | ---: |
| ASP272/ASP274 | 1.02 | 3.44 | 9.70 |
| ASH272/ASP274 | 0.88 | 6.95 | 8.93 |
| ASP272/ASH274 | 2.37 | 5.84 | 9.62 |

The other three beta sites vary substantially within the same snapshots.
Across the three structures and four beta sites, Asp272 spans 0.88–6.67,
Asp274 spans 0.83–7.98, and Thr233 N terminus spans 6.53–9.87. The
complete 36-value table is the authoritative result. All three processes
exited successfully and every target group appeared in its final summary
table. Several catalytic Asp entries carry PROPKA's coupling asterisk in
the detailed output. The `--display-coupled-residues` option was not used,
so alternative coupled-site assignments have not yet been inspected.

These predictions provide **no consistent support for specifically
protonated Asp272** at the substrate-facing E site. Some Asp274 estimates
are elevated, while other sites favor Asp272; the result depends strongly
on conformation. Most N-terminal estimates exceed pH 8 and would favor
NH3+ in a simple isolated-site titration model, but terminal coupling,
the missing substrate, and an uncertain assay pH at 75 °C limit that
interpretation. No force-field protonation state should be changed solely
on this screen. A ligand-aware calculation or a method that explicitly
samples coupled microstates would be needed for stronger evidence.

## Relation to Negoro et al.

Negoro et al. (*FEBS J.* **290**, 3400–3421, 2023,
[doi:10.1111/febs.16755](https://doi.org/10.1111/febs.16755)) proposed
Asp306-COOH as a general acid in **precursor autocleavage**. Their support
was a drawn proton-transfer route, Asp306 proximity to the precursor
scissile N (3.84 Å to Cys267 amide N and 2.66 Å to Thr267 amide N in a
structural superposition), abolition of maturation by D306A, and
hydroxylamine sensitivity consistent with a thioester intermediate. They
did not directly measure Asp306 protonation or a microscopic pKa. In
**mature substrate hydrolysis**, they proposed an Asp308–Asp306–Thr267
relay by analogy with BapA; that does not independently establish the
Asp306 microstate. Our audit of all their Table 1 deposited structures
found closest Asp306/Asp308 O–O distances of 3.58–4.17 Å, without a
short direct Asp–Asp hydrogen bond. See `docs/negoro_2023_mechanism_audit.md`.

## Consolidated documentation and corrected crystal comparison (2026-11-19)

The manuscript Methods and Results now consolidate both the Nyl12 screen
and the corrected mature/uncleaved NylC comparison, including tables.
NylC 7YU0/7YU1 retain their deposited maturation-blocking mutations; they
were not back-mutated to WT. Cleavage effects are therefore confounded by
sequence and conformation.

Use mature 5XYG values only from
`analysis/negoro_2023_protonation/crystal_propka_20261116/attempt3/results/asp306_asp308_pka.csv`.
The TER-correction audit, including N-terminal predictions, is
`analysis/negoro_2023_protonation/terminal_pka_20261116/before_after.json`.
Precursor 7YU0/7YU1 values are in the attempt2 table; **exclude its superseded
5XYG rows**. The student `evidence/pka/` copies provide that filtered
comparison and full Nyl12 table with source hashes.

Mature 5XYG: Asp306 2.19–2.21; Asp308 5.49–6.49 default or 7.24–8.45
coupled alternative; Thr267 N terminus 6.35–6.55 default or 4.59–4.60
alternative. Precursor predictions are chain/coupling dependent. These are
static protein-only predictions, not protonation measurements.

The mature NylC MCCE4 result was rejected because internal missing-coordinate
boundaries generated artificial termini. Its apparent Asp <0/N-terminal >14
bounds must not be treated as accepted evidence. Complete Nyl12 MCCE step 1
also failed (zero NTR/CTR groups; 2,840 diagnostics); no full Nyl12 titration
was run. Rejection receipts are preserved under
`analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/` and
`analysis/nyl12_mcce_step1_preflight_20261116/neutral_replica_02_frame284_attempt2/`.
