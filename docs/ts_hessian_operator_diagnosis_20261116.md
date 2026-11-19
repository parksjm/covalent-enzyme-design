# Direct-route TS Hessian operator diagnosis

Status: numerical diagnostic completed 2026-11-16; **first-order saddle
not yet verified**. The 0 K DFTB3/MM dimer-refined direct-route structure
converged at cycle 99. Its prior block-LOBPCG probe estimated one
negative and one positive low eigenvalue but left residual norms
0.0505 and 0.0608 Eh/(Å² amu), far above the 10⁻⁴ target.

The bounded diagnostic used the same full 4,165-atom movable region,
topology, optimized restart, calibrated native Sander gradient and zero
electronic smearing. Five central Hessian-vector products required ten
native gradient evaluations. The chemical dimer mode and an orthogonal
seeded random direction were tested at 0.05 and 0.10 Å√amu; the chemical
mode product was also repeated at 0.05. The exact command, input hashes,
environment and logs are in
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/ts_image08_attempt1/hess_job_20261116/`;
the scalar result is in its sibling `hess_op_20261116/diagnosis.json`.

| Check | Result |
| --- | ---: |
| Repeated chemical-mode HVP relative difference | 0 |
| Chemical-mode HVP relative change, 0.05 vs 0.10 Å√amu | 0.93% |
| Random-mode HVP relative change, 0.05 vs 0.10 Å√amu | 0.085% |
| Bilinear symmetry relative difference, v·Hw vs w·Hv | 0.37% |
| Chemical-mode Rayleigh quotient | −0.07862 Eh/(Å² amu) |
| Orthogonal random-mode Rayleigh quotient | +0.39618 Eh/(Å² amu) |

These bounded checks do not show gross nondeterminism, a force-sign
reversal, or large step-size sensitivity. They do not establish that the
full Hessian has exactly one negative eigenvalue: a random positive
direction cannot exclude other unstable modes. The stalled ARPACK and
block-LOBPCG searches may reflect ill conditioning or many soft
protein/solvent modes in the 12 Å active region. Do not spend hundreds
more gradients on the same unpreconditioned search unchanged. Next,
audit feasible preconditioning or a chemically focused mode subspace,
then verify any proposed second-lowest mode in the **full** movable
region. Two-sided downhill connectivity from a near-identical older
candidate is supportive but also does not establish Hessian index one.

Three setup attempts preceding the successful run failed before usable
gradients: missing repository import path, missing Amber shared-library
path, and Fortran input-filename truncation at an absolute path over
255 characters. Their job logs are retained in the same TS directory;
the successful run used the prior validated Amber environment and a
shorter path without changing scientific inputs.

## Bounded warm-start continuation (2026-11-16)

The existing block spectrum contained **partial Ritz vectors only**; no
converged vectors were available to reuse. Added a warm-start option that
loads those two saved vectors and orthonormalizes them for the next block
LOBPCG solve. The full 4,165-atom movable region, fixed exterior, TS restart,
topology, dimer mode, 0 K calibrated native Sander gradient, and 0.05
Å√amu displacement were unchanged. The pilot was capped at 40 HVPs and
stopped after 38 HVPs (76 native gradients).

The two Ritz estimates were −0.230595 and +0.004731 Eh/(Å² amu), with
residuals 4.132×10⁻⁴ and 2.197×10⁻². Compared with the prior 68-HVP block
probe, the residuals improved by factors of 122 and 2.77, respectively;
the residuals still exceed the 10⁻⁴ acceptance threshold. This is partial
progress, not evidence of Hessian index one or TS acceptance. The positive
estimate moved closer to zero, so its identity and stability require further
full-region verification before interpretation.

The process record, manifest, spectrum, and separate scientific acceptance
record are under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/ts_image08_attempt1/mode_probe_block_warmstart_shortpath_20261116/`.
The short `/tmp/ts_hess_work` symlink kept Fortran filenames short while
writing the run artifacts into this existing TS directory. Two earlier
launch attempts failed before Hessian products (output-directory collision;
then Fortran path truncation); both attempts and logs are retained separately.

## Wider-block follow-up (2026-11-16)

The prior two-vector iteration history showed the first residual approaching
4×10⁻⁴ while the second oscillated between roughly 0.022 and 0.039. That is
consistent with unresolved mixing among low modes, but does not prove a
cluster because no additional Ritz values were available. A six-vector
lowest-algebraic LOBPCG search was therefore run on the unchanged full active
region, seeded by the dimer direction, both saved partial Ritz vectors, and
three deterministic independent random directions. At the same 0.05 Å√amu
step, the 60-HVP cap was reached before Ritz pairs and their full-Hessian
residuals could be returned. The process wrapper exited 0 after writing an
explicit `hvp_budget_exhausted` result; this is not scientific convergence.

The exact job, manifest, evaluation artifacts, spectrum status, and separate
acceptance record are under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/ts_image08_attempt1/`:
job metadata is in `mode_probe_block6_warmstart_20261116/`, and scientific
results are in `block6_warmstart_20261116/results/`. This bounded follow-up
does not distinguish soft-mode crowding from insufficient block-iteration
and residual budget. It produced no signed low-mode estimates, so the number
of negative modes and numerical index one remain unestablished. Do not infer
that this candidate is a first-order saddle from either process completion or
the six-vector budget-limited run.

## Final bounded six-vector run (2026-11-16)

Added an incremental orthonormal cache of completed `(v, Hv)` pairs. On a
hard HVP interruption, the probe can now diagonalize the projected symmetric
Hessian and report lifted Ritz vectors with full-space residuals from cached
HVP images. Synthetic tests verify exact Ritz recovery and residual
calculation. This fallback was not exercised by the final native run: its
20-iteration limit stopped the solve after 144 HVPs (288 gradients), below
the 180-HVP cap. The run was fresh, seeded with the dimer direction, the two
prior partial Ritz vectors, and three deterministic random vectors; SciPy
LOBPCG exposes no continuation state, so this was not a restart.

The six signed Ritz estimates were −0.230524, +0.002467, +0.003150,
+0.003267, +0.003528, and +0.003748 Eh/(Å² amu). Their residuals were
0.007598, 0.012981, 0.013607, 0.014536, 0.014987, and 0.015606. All fail
the 10⁻⁴ gate. The five positive estimates lie within a narrow interval, but
their residuals are larger than the spacings between them, so neither their
individual signs nor the number of negative modes is resolved. Residuals for
the lowest mode initially fell below 6×10⁻⁴, then increased to 7.6×10⁻³;
the other five ended between 1.3×10⁻² and 1.6×10⁻². This supports unresolved
soft-mode mixing/poor subspace conditioning but does not establish the
spectral cause. Numerical index one remains unverified.

The full 23×6 residual history and separate scientific acceptance record are
in
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/ts_image08_attempt1/block6_hvp180_20261116/results/`;
job status and launch provenance are in the sibling
`mode_probe_block6_hvp180_20261116/` directory.
