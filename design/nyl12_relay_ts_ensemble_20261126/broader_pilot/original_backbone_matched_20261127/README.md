# Matched original-backbone Boltz/Amber cohort

Authorized and queued 27 November 2026. Candidates relay_prot101_06 and relay_prot101_03 are ranks 1/2 under the original-only local TS likelihood diagnostic (0.611583/0.636214). They contain 92/99 canonical substitutions and differ at 41/328 precursor positions. All101 protected identities and four-copy sequence ties passed the selection audit. Likelihood is not a predicted catalytic rate.

Two exact-query local ColabFold searches run in separate tmux sessions on Violet GPUs0/1 with two CPU threads each. The validated WT precursor MSA is reused, but both WT predictions are fresh under the same Violet stock Boltz runtime as the designs. Mature alpha/beta homologs are paired using the common precursor origin and query-aligned232|233 split.

Six predictions: WT plus both designs, seeds20261126/20261127, eight mature protein chains and four unguided PA6L4 substrates. No protein template, contact restraint or fixed substrate register. Native CPU feature acceptance checks actual pairing and stock dummy templates. A separate real WT native CPU diagnostic already passed.

`jobs/matched_pipeline01/` owns session `nyl12_original_matched_20261127`. It waits for both search jobs, then uses two balanced three-model GPU queues, followed by two matched Amber queues. No interactive request is needed to start relaxation. Exact commands, frozen source manifests, cache identities and process receipts are retained.

Canonical chain/atom maps and PA6 connectivity are checked before relaxation. Severe interchain heavy-atom contacts below1Å reject pathological assemblies before Amber. Remaining models use matched unrestrained ff19SB/TIP3P/150mM/modularRESP-GAFF2 minimization, at most10000cycles, RMS0.05kcalmol^-1Å^-1; ThrN-termini neutral and Asp272/Asp274 anionic. Force convergence and chemical/geometric acceptance are separate.

All48 site/ligand/amide combinations are evaluated per snapshot. Positive signed phi and TierA/B gates precede one-to-one assignment. Native-Thr similarity is diagnostic. Post-Amber exports restore canonical identities from actual topology maps and retain sequential Amber numbering separately. Canonical PyMOL-ready structures and before/after metrics are saved per case. No fixed amide or residue929 assignment is assumed.

`controller/summary.json` and `controller/comparison.csv` are generated after all six cases; detailed all-four-site/all-amide evidence is in each case postprocessing directory. Failure receipts remain in the census. QM/MM and umbrella-sampling jobs were left unchanged.
