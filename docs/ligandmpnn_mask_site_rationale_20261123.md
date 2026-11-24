# LigandMPNN mask rationale and flipped-pose follow-up, 2026-11-23

All positions below use the 328-residue full precursor: alpha 1–232 and beta 233–328. Beta local 39 and 41 are precursor 271 and 273. The focused mask is the previously specified local binding-pocket design set, retained unchanged for this comparison; it is a design hypothesis, not a claim that each residue is catalytic. Its historical position-level prior includes the modeled 7 Å pocket described in `docs/nyl12_ligandmpnn_pilot.md`, not a newly applied uniform 7 Å cutoff.

- **Focused 19 eligible positions:** 50, 51, 52, 73, 76, 101, 111, 112, 114, 115, 116, 151, 152, 153, 185, 187, 188, 271, 273.
- **Broader 29 eligible positions:** 28, 37, 38, 50, 51, 52, 73, 76, 92, 96, 97, 101, 111, 112, 114, 115, 116, 151, 152, 153, 160, 161, 162, 164, 185, 187, 188, 271, 273.
- The ten added alpha positions are 160, 161, 162, 164 from the 6–12 Å shell, and 28, 37, 38, 92, 96, 97 from the 12–16 Å shell. A candidate qualified for its shell if at least one copy in either accepted WT parent had substrate distance in the respective interval (6 < d ≤ 12 or 12 < d ≤ 16 Å). Added sites required relative solvent accessibility ≥0.25 in every equivalent copy of both parents and no protein heavy-atom contact to another chain ≤4 Å in either parent. This screen applies to additions; the original focused sites were retained without retroactive application of these expansion gates.

The two parents were accepted force-converged WT seed 02 and seed 03 mature tetramers, each with at least one occupied site passing the existing geometry and native-Thr gates. Explicit symmetry ties linked A–D alpha and E–H beta equivalents with equal 0.25 weights; chain lengths and sequences were audited. Hard protections are alpha 113, 156, 186, 222–232 and beta local 1–3, 40, 42 (precursor 233–235, 272, 274). These cover catalytic/maturation context; no claim is made that every protected position is independently proven catalytic. The 2 masks × 2 parents × 2 temperatures × 6 samples yielded 48 raw samples with paired parent/temperature RNG seeds; eligible positions need not mutate.

The selected **mask_03** has S50A, V51A, S153A, A187T, V188I, L271F (beta local L39F); **mask_07** has S50A, V51A, A92K, N97D, S153T, A187T, V188T. Exact sequence/chain verification and hashes are in `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/sequence_comparison.json`. These substitutions preserve the stated hard protections.

With a consistent scissile-amide mapping, opposite signed projected Herschlag phi values identify opposite local orientations; in this campaign negative phi is the flipped orientation. The PI-authorized flipped-only follow-up tests a possible Thr233 N-terminal proton relay to the leaving amide nitrogen; it does not imply that this orientation places the two nitrogens closer. The saved relaxed structures instead show the following Nα–leaving-N distances (Å), by assigned site:

| Model | AE | BF | CG | DH |
| --- | ---: | ---: | ---: | ---: |
| WT seed 04 | 4.59 | 4.87 | 4.75 | 4.83 |
| mask_03 seed 03 | 4.45 | 5.00 | 4.30 | 5.19 |
| mask_07 seed 04 | 4.89 | 4.30 | 4.83 | 4.58 |
| V111C seed 02 | 5.31 | 5.45 | 5.26 | 5.33 |

Across the five WT endpoints, negative-phi sites (n=16) average 4.65 Å (4.24–5.47 Å), while positive-phi sites (n=4) average 3.61 Å (3.51–3.67 Å). These are descriptive, nested structural observations, not independent activity measurements. Direct transfer from these flipped starting geometries would require rearrangement; water mediation is another hypothesis. Neither shorter positive-phi distances nor pose retention proves a mechanism or exclusive binding. All historical orientations remain recorded. The current campaign does not add positive-pose MD and uses no orientation-holding production restraints. Source: `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/analysis/occupied_sites.csv`.

The follow-up plans two independent 100 ns replicas each for mask_07 seed 04, mask_03 seed 03, V111C seed 02 and WT seed 04, with four PA6 L4 ligands and matched neutral Thr233/ASP272/ASP274 conditions. All selected starting sites have negative phi. WT seed 04 has four Tier A sites and zero native-Thr passes; its initial conformation and the selected designs differ, so comparisons are conditional on starting structures, not activity measurements. Queue receipts under `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/queues/` are authoritative for execution status.

Provenance: `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/protocol.json`, `preflight/final_masks.json`, `preflight/surface_exposure.json`, `preflight/parent_and_mask_proposal.json`, `preflight/symmetry_protection_acceptance.json`, `preflight/audit_surface_masks.py`, `generation/results/acceptance.json`, and `validation/prepared/candidate_manifest.json` (paths after `protocol.json` are relative to the ensemble comparison directory). See also `validation/analysis/occupied_sites.csv` and `md_followup/selection.json` in that directory.
