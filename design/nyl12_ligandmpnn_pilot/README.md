# Nyl12 LigandMPNN pilot

Completed successfully: **24 samples, eight unique sequences**, with four to six
unique precursor-position substitutions each. All sequence/tie/protection and
tensor/FASTA checks passed; outer job exit status is zero.

This is the PI-approved exploratory E–S sequence-design pilot from WT **Model23**.
The original Model40/41 pair failed native-Thr geometry preflight; its artifacts
remain evidence and were not used for inference. Model23's earlier QM/MM
geometry deterioration remains a separate negative result.

The reproducible run is `accepted_parent/run/run_approved_pilot.sh`; execution is
launched in tmux on one low-priority CPU on violet, with GPUs hidden. Read
`accepted_parent/run/job/status.json` for process status and
`accepted_parent/run/results/status.json` for sequence-generation/audit progress.

Four jobs generate six samples each: sampling temperatures 0.1 and 0.2, with two
independent random seeds per temperature. Temperatures here are dimensionless
sequence-sampling controls. No automatic resampling follows duplicates.

Useful results after completion:

- `accepted_parent/run/results/unique_sequences.fasta`: deduplicated eight-chain
  sequences, ordered alpha A/B/C/D then beta E/F/G/H; colon separators.
- `accepted_parent/run/results/all_samples.fasta`: all 24 raw samples.
- `accepted_parent/run/results/samples.csv`: precursor-position mutations, settings
  and alpha/beta sequences; beta physical residue numbers have offset +232.
- `accepted_parent/run/results/report.json`: duplicate lineage and full sample audit.
- `accepted_parent/run/results/scientific_acceptance.json`: fixed/protected sites,
  chain lengths and explicit alpha/beta copy ties.
- `accepted_parent/run/results/tensor_fasta_acceptance.json`: independent native
  parser and generated tensor/FASTA consistency checks.
- `accepted_parent/backbones/pa6_l4z_l2_nyl12_wt_model_23/ligandmpnn_input.pdb`:
  input tetramer and PA6 L4, physical alpha A–D / beta E–H / substrate I.
- `accepted_parent/run/preparation_manifest.json`, `launch_preflight.json` and
  `input_source_manifest.json`: exact commands, settings and input/source hashes.

The 19 designable precursor positions are 50,51,52,73,76,101,111,112,114,115,116,
151,152,153,185,187,188,271,273, tied across all four copies. Catalytic and
maturation protections remain fixed. Only one site contains substrate; sequence
symmetry does not impose identical conformations on unoccupied sites.

Generated PDB outputs are **not relaxed candidate structures**: side-chain
packing was disabled. These are candidate sequences, without validated folding,
binding, catalytic geometry or activity. No activity ordering is inferred from
noisy crude-lysate MS fold signals. Next comes human review and a separately
approved prediction/relaxation campaign with WT controls.

The existing laptop pull script copies this directory with the same layout,
excluding model weights, raw Torch statistics and rejected preparation trees.
Original data on violet are preserved. Scientific QM/MM, US and WT MD jobs run
independently and are not modified by this pilot.

The first startup stopped before sampling because bundled OpenFold source was
missing from the initial snapshot. Those failed receipts remain under
`accepted_parent/job/` and `accepted_parent/results/`. The active run uses
`accepted_parent/run/` and copied, hashed installed support files; no new
dependency was installed.

Verification: the executor's six focused tests pass (82% coverage); the FASTA
validator's seven tests pass (85% coverage). Ruff checks pass for those modules.
The laptop script passes Bash syntax and its complete violet-only dry-run test.
Preparation/runtime artifact scripts pass syntax/fatal-lint checks and actual
input/parser gates; their code coverage was not measured. Native process
completion, sequence integrity and biological performance are separate gates.

## Generated variants

Counts are generator sampling frequencies, not activity rankings. Each mutation
is introduced in all four equivalent copies. Numbering is the 328-aa precursor.

| Substitutions | Samples / 24 |
|---|---:|
| S50N, V111T, A187S, V188T | 10 |
| S50A, V111T, A187S, V188T | 7 |
| S50N, V111T, T116I, A187S, V188T | 2 |
| S50A, V111T, T116I, A187S, V188T | 1 |
| S50A, V51T, V111T, T116I, A187S, V188T | 1 |
| S50A, V111T, Y115F, T116V, A187S, V188T | 1 |
| S50A, V51T, V111T, T116L, A187S, V188T | 1 |
| S50N, V111T, I112L, A187S, V188T | 1 |

Repeated V111T/A187S/V188T substitutions reflect this single-parent, fixed-
backbone sampling context. They do not establish functional benefit. Preserve
all eight unique variants and WT for the next human-reviewed structural test.
