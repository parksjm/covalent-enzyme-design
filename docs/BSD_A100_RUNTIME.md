# DGX access, resources and bounded scientific work

Verified 27 November 2026. This is a shared workstation; visible hardware is
not an exclusive allocation. The PI authorizes **at most two GPUs** for this
project and requires the existing QM/MM and umbrella-sampling jobs to continue.

## Access and paths

The laptop alias is `ssh bsd-a100`, user `jpc`, authenticating with the laptop's
`~/.ssh/id_ed25519`. The PI subsequently enrolled Violet's existing public key,
allowing management from this project with:

```bash
ssh -i /home/jpc/.ssh/id_ed25519_aster \
  -o BatchMode=yes -o StrictHostKeyChecking=yes jpc@bsd-a100
```

No private key was transferred. The PI installed Codex CLI 0.172.0 and connected
the Desktop App. A Desktop connection alone does not expose cross-host command
execution to this Violet shell; the separate SSH route was necessary.

Use `/home/jpc/code` for existing software and `/home/jpc/projects` for work.
There is no `/data`. Both home paths are on the root filesystem, which had
only 16.7 GiB free at audit; the PI is cleaning home storage. `/raid` has about
3.37 TiB free but its top level is not writable by `jpc`; do not assume an
authorized writable directory there. Reuse existing weights, stage small
inputs and enforce a 10 GiB free-space floor for control preparation. Do not
delete existing data or copy/download additional weight libraries blindly.

## Verified resources

| Resource | Observed state |
| --- | --- |
| GPUs | Eight NVIDIA A100-SXM4, 40 GiB each |
| Other-user GPU work | GPUs 0 and 1 occupied at about 99% utilization |
| Proposed project GPUs | 2 and 3; idle at audit and compatibility check |
| CPU | Two AMD EPYC 7742 sockets; 128 physical cores, 256 hardware threads |
| RAM | About 1 TiB total, 954 GiB available at audit |
| Scheduler | No `sinfo`, `squeue`, `srun`, `sbatch` or `qstat` on inspected PATH |
| Containers | Docker installed; this account lacks socket access |
| Basic tools | `uv`, Python, tmux, rsync and Git present |

Use physical UUIDs rather than relying solely on ordinal device names:

- GPU 2: `GPU-a752cdba-475a-6afe-14c5-4a58ea101a56`
- GPU 3: `GPU-37c49915-4468-ef72-3b24-e0f9bcfcd9ee`

Refresh process occupancy immediately before launch. Restrict each independent
prediction to one UUID with `CUDA_VISIBLE_DEVICES`; serialize architectures
on each GPU so at most two predictions execute concurrently. Bound CPU thread
counts and affinity; the project is not reserving all 256 hardware threads.

## Existing runtimes and qualification

**AF3 compiler correction, 27 November:** the first two full WT inferences
failed because the launch environment selected `/usr/bin/ptxas` (CUDA 12.0,
PTX 8.0) while JAX emitted PTX 8.4. Using the AF3 interpreter by absolute path
does not activate its executable search path. Prepend
`/home/jpc/code/alphafold3/.conda/envs/af3/bin` to `PATH`; this selects the
already installed CUDA 12.6 assembler. Explicit PTX 8.4 assembly for `sm_80`
and a JAX CUDA kernel passed after that correction. Two unchanged-input,
unchanged-model retries are running as `nyl12_af3_wt_20261126_ptxas126` and
`nyl12_af3_wt_20261127_ptxas126`, with new `_ptxas126_retry01` outputs.
Original failures remain preserved. Both full retries finished successfully
at 12:25:55 EST and passed native output atom/sequence inventory. Their
catalytic geometry nevertheless failed, as described below. Compiler/runtime
acceptance is separate from structural or catalytic acceptance.

Boltz uses `/home/jpc/code/miniforge3/envs/boltz-env/bin/python` and its `boltz`
entry point: Boltz 2.2.1, Torch 2.9.0+cu128, Python 3.11.14. Existing weights and
CCD cache are in `/home/jpc/code/boltz`; use the explicit `boltz2_conf.ckpt`
checkpoint and the existing `mols` directory. Actual CUDA matrix multiplication
passed on GPUs 2 and 3 with Ampere capability 8.0. This qualifies that small
kernel, not full-model memory use or a scientific prediction.

AF3 is installed in `/home/jpc/code/alphafold3/.conda/envs/af3`, version 3.0.1,
with JAX/JAXlib 0.4.34 and existing weights in
`/home/jpc/code/alphafold3/models`. The installed native `folding_input` module
imports successfully and supports `userCCDPath`. Do not force
`PYTHONPATH=/home/jpc/code/alphafold3/src`: that source tree lacks the compiled
`alphafold3.cpp` extension and caused a diagnostic import failure. Use the
installed package and verify the exact entry point/flags before inference.

Streaming SHA256 comparison confirms that both existing DGX checkpoints match
Violet byte-for-byte: Boltz conformation checkpoint
`090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`
and AF3 compressed weights
`333320e3cb2532b07541ee165e554a3bc45a12f0faab2f374a0e37520185a685`.
Retain the runtime/source checks as well; weight identity alone does not
certify complete numerical equivalence between software installations.

Amber/Sander/pmemd are not on the inspected PATH and are **not qualified on
this host**. Do not infer compatibility from CPU architecture or GPU presence.
The current DFTB3/Sander umbrellas are CPU jobs; an A100 does not accelerate
their present implementation. Read the Amber runtime policy before any future
Amber GPU launch. No new dependency or model-weight installation was performed.

## Priority and scientific gates

**Current calibration stage, 27 November:** the query-only controls below are
complete and canonicalized. Additional matched-seed 20261126 controls use the
existing 3,137-row alpha and 4,913-row beta protein MSAs, preserving unpaired
homologs (no invented pairing); all four PA6 molecules remain unguided. Both
native MSA controls completed at13:29 EST, with successful canonical coordinate
and metric roundtrips but zero assigned Tier A/B sites. This triggers at most
one prepared protein-only crystal-template follow-up per predictor after fresh
GPU/free-space guards. An inference/identity failure stops the conditional
workflow. Template preparation preserves observed coordinates, missing-residue
gaps and the H228 ALA/query PRO substitution. AF3 receives eight individual
mature-chain templates, so no interchain assembly restraint is implied; Boltz
uses a nonforced joint protein template. No ligand/TS geometry enters either.
Use the same GPU 2/3 UUIDs and native environments above, particularly AF3's
CUDA 12.6 executable PATH. Both MSA and template source/receipt workspaces are
under `wt_msa_controls_20261127` and `wt_protein_template_controls_20261127`
beside the original `host_qualification_20261127` workspace. No new models,
packages or model server were installed. Old statuses below are historical.

Historical execution checkpoint at 12:11 EST, 27 November: both WT Boltz controls finished
with exit zero and no failed examples, approximately 5 min 18 s each including
startup. Their native output inventories passed the existing acceptance API:
eight mature protein chains and four 33-heavy-atom ligand chains per model.
The two AF3 jobs were released by the GPU-idle helpers after Boltz exited and
remain active on the same physical GPU 2/3 slots. Sessions are
`nyl12_af3_wt_seed_20261126` and `nyl12_af3_wt_seed_20261127`; individual native
logs/status and the finished Boltz receipts are under the workspace's `jobs/`.
No scientific convergence, catalytic-geometry or activity claim follows from
these process/software checks.

1. Prepare and qualify the existing **WT controls**, seeds 20261126/20261127,
   with Boltz and, after native input validation, AF3. These consume four WT
   control slots within the existing eight-prediction comparison budget.
   They are independent of redesigned-parent acceptance and are not a new
   mutant campaign. Preserve the mature eight-chain assembly, four zwitterionic
   PA6 L4 molecules, query-only MSAs, and absence of contacts/templates/TS input.
2. Release redesigned variants only after the existing backbone/protection/
   ligand-graph gates pass. RFD3 native exit zero did not establish scientific
   acceptance; the two current parents remain unaccepted. No gate is relaxed
   merely to occupy an available GPU.
3. Reserve later GPU capacity for the qualified symmetry-aware LigandMPNN and
   unguided prediction workflow, then a very small validated-MD shortlist.
   Additional long MD or model-server work is not the current priority.

The bounded control workspace is
`/home/jpc/projects/nyl12_relay_redesign/host_qualification_20261127`.
Before releasing the remaining four redesigned-sequence predictions, bind the
qualified WT inputs/runtime/output receipts into the comparison. Do not run
the untouched eight-job controller in addition to these four WT controls;
count all completed controls against the agreed total and review comparability.
Input parsing, successful inference, structural inventory and catalytic geometry
are separate acceptance records. Host qualifications do not automatically make
older and newer software outputs a matched scientific comparison; retain
versions/checkpoints/inputs and review comparability before pooling results.

Before exporting four-site geometry, establish the physical alpha–beta partner
assignment in the unrestrained predictions. Symmetry-equivalent chain labels
can permute: the current `relay_any_amide` implementation assumes A/E, B/F,
C/G and D/H. Verify this against the predicted assembly and crystal-core
alignment, or make an explicit validated derived relabeling; retain the raw
outputs. Otherwise a correct assembly could be scored with mismatched alpha
and beta partners. This is a postprocessing prerequisite, not evidence that
the new models have passed or failed catalytic geometry.

For the two completed Boltz controls, protein-only contact assignments uniquely
confirm A/E, B/F, C/G and D/H; crystal-core residue identities also match.
All eight assigned sites fail even Tier C; their OG1–carbonyl-C distances are
10.16–28.26 Å. This does not measure enzyme activity. The unrelaxed query-only
eight-chain WT prediction protocol has not recovered a catalytically competent
E–S baseline. See the complete [control comparison](../design/nyl12_relay_ts_ensemble_20261126/broader_pilot/dgx_wt_control_qualification_20261127/comparison_20261127/comparison_report.md).

The two corrected AF3 controls also completed and passed output inventory;
all eight AF3 sites reject even Tier C. Across both predictors, all 16 assigned
sites and all 192 register/ligand/site comparisons reject. The contact-based
physical pairs are AE/BF/CG/DH in all four models. Protein-only oxyanion donor
separations are severely non-native, independently of substrate placement.
The two DGX project GPU tasks have now ended naturally. Do not automatically
launch additional variants or treat these WT models as a validated activity
filter. The runtime is usable; protein-pocket modeling still needs calibration.

The read-only audits and CUDA receipt are under
[`provenance/bsd_a100_20261127`](../provenance/bsd_a100_20261127/).
The separate design readiness report is
[`dgx_gpu_readiness_20261127.md`](../provenance/dgx_gpu_readiness_20261127.md).
