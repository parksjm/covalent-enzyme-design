# Amber GPU runtime: mandatory workstation preflight

Use the **execution host's validated binary**, not an executable copied from
another workstation. This distinction has already been resolved repeatedly.
No installation, GPU-driver update, reboot or change of scientific method is
needed for the known Azalea failure described below.

| Host | GPU architecture | Validated executable |
|---|---|---|
| Violet | RTX 6000 Ada, 48 GB | `/data/jpc/code/pmemd26/bin/pmemd.cuda_SPFP` (`pmemd.cuda` resolves here) |
| Azalea | RTX PRO 6000 Blackwell Max-Q, 96 GB | `/data/jpc/code/pmemd26_blackwell_build_20261005/pmemd26_src/build_sysgcc_cuda128_v2/src/pmemd/src/pmemd.cuda_SPFP` |

Validated SHA-256 values:

- Violet Ada: `7ed3f439d09d6beff94210bb3f356585e6918cd3f500dcc58046067ac8fd9602`
- Azalea Blackwell CUDA 12.8:
  `f4b0b32eda07d1e6b4fbe597e8c3bd4d39cc492a2b01c19c13db7ebe92ff52b1`

The generic `/data/jpc/code/pmemd26/bin/pmemd.cuda_SPFP` **also exists on
Azalea** and has the Violet hash. It is the wrong build for Blackwell. File
existence and a clean `ldd` result do not test GPU kernels. Its known error is:

```text
cudaMemcpyToSymbol: SetSim copy to cSim failed no kernel image is available for execution on the device
```

## Required preflight and explicit environment

1. Read the current resume note and this guide. Check user sessions, GPU/CPU
   processes and queued work, not only the instantaneous utilization.
2. Pin the exact host-specific executable above and verify its hash. A changed
   hash needs new provenance/compatibility assessment, not silent acceptance.
3. Use `AMBERHOME=/data/jpc/code/ambertools26/venv`. Preserve the exact successful
   library environment in the job manifest. Check for missing shared libraries.
4. Before a substantial queue, confirm that its first native stage actually
   advances finite MD/minimization steps on the intended GPU. A submitted tmux
   session is not a passed native compatibility test. Preserve any failure.
5. Hash the executable and inputs; launch only from an explicit, frozen plan.
   Do not silently switch to a CPU engine or another force field/protonation.

For the corrected November 23 Azalea design-MD launch, the explicit library
environment is:

```bash
export AMBERHOME=/data/jpc/code/ambertools26/venv
export LD_LIBRARY_PATH=/data/jpc/code/cuda128_findcuda_compat/lib64:/data/jpc/code/pmemd26/lib:/data/jpc/code/ambertools26/venv/lib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
# Set CUDA_VISIBLE_DEVICES to the explicitly assigned physical GPU, 0 or 1.
```

The Blackwell build also has recorded runtime search paths. A historical
successful frame1603 command did not require a recorded library override.
Neither `CUDA_FORCE_PTX_JIT` nor electronic smearing is a remedy for selecting
the wrong executable. No such workaround is introduced here.

Verified Azalea orchestration Python is
`/data/jpc/code/ambertools26/venv/bin/python` (Python 3.14, NumPy 2.5.3,
ParmEd 4.3.1 on this date). Violet uses the pinned project uv environment for
the new design-MD controller. Per-host runtime records distinguish these
orchestration environments; the Amber Hamiltonian and MD protocol are unchanged.

## Provenance and recurrence

The accepted historical native command and executable hash are recorded in
`qmmm/pa6_resp_gaff2_wt/09_mm_qmmm_selected_candidates_20261106/frame_1603_attempt2/provenance/attempt_01/command.json`.
The existing explanation is also in
`docs/qmmm_resp_gaff2_restart_20261104.md`, “Three-seed unrestrained MM outcome.”

On November 23, the new design follow-up incorrectly selected the generic
Azalea binary again. Two children failed before any MD timestep. This was a
preparation error, not a chemistry or force-field failure. Failed native outputs
and original plans are preserved under the remote campaign's
`clutter/azalea_generic_cuda_startup/`; failed controller receipts remain in
`queues/azalea_gpu*/`. Corrected plans explicitly name the Blackwell executable
and use `configs/azalea_blackwell_*.json` and
`configs/queue_azalea_gpu*_blackwell.json`. Native advancement must be checked
before declaring the corrected runs active.

Campaign root on both hosts:
`/data/jpc/covalent-enzyme-design/design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/`.
The persistent instruction is also in root `AGENTS.md` so future agents encounter
this rule before launching jobs.
The design follow-up's `launch_host.py` additionally checks the exact host's
executable hash and frozen config selection before preparing or launching a job.
Its tests reject both a changed executable and a config naming the other build.
The corrected native runs have now passed heating, density and weak-restraint
stages; GPU utilization and native stage-acceptance records confirm advancement.
