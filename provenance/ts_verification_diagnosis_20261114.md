# Dimer TS and direct NEB diagnosis (2026-11-14)

## Scope and artifacts

This is a read-only diagnosis of the PA6 model30 neutral acyl-enzyme path. No new
optimization, Hessian calculation, or scientific job was launched.

- Dimer candidate: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113/calculation/optimized.ncrst`.
- Direct route: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/`.
- Low-mode probes: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/mode_probe_attempt{2..6}/spectrum.json`.
- Force calibration: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113/native_gradient_0K_onecycle_20261113_attempt2/scientific_acceptance.json`.

The dimer candidate and the direct-route seed's central structure have the same
mapped reaction-center geometry (Å): Thr OG1–scissile C 2.085, C–Nleave 1.538,
and Thr HG1–Nleave 1.218. This is a coupled, asynchronous addition/proton
transfer structure, with the C–N bond still substantially intact; it is not a
late acyl-enzyme product geometry. The trace's `nebinfo` relative profile peaks
at about 0.0571 Eh near path length 8.14. Its NEB setup uses the same topology,
QM charge, active shell, and endpoint policy as the dimer candidate (4,165
active atoms, 98 QM atoms, QM charge -2). The profile is consistent with
continuation from the dimer-derived path, but the peak alone does not establish
a stationary point or a validated barrier. The trace README explicitly leaves
runtime geometry acceptance outstanding; the Sander exit status is not that
acceptance.

## Why attempts 2–6 returned no low modes

All five `spectrum.json` files contain no eigenpairs and mark the run
`incomplete_or_partial` because the finite-difference Hessian-vector-product
budget was exhausted:

| Attempt | FD step (Å√amu) | HVP cap reached | Result |
| --- | ---: | ---: | --- |
| 2 | 0.01 | 60 | process failed at the HVP budget |
| 3 | 0.01 | 200 | budget exhausted |
| 4 | 0.1 | 200 | budget exhausted |
| 5 | 0.1 | 600 | budget exhausted |
| 6 | 0.1 | 800 | budget exhausted |

Attempts 3–6 exited their Python wrappers with status 0 because partial
diagnostics were written; they did not converge eigenpairs. These were ARPACK
Lanczos `eigsh(which='SA', k=2)` calculations in 12,495 mass-weighted
coordinates. Increasing the HVP cap and ARPACK iteration/subspace settings
through attempt 6 did not cure the issue. The failure is therefore not
explained by a single small step alone; the request for the two globally lowest
modes over a large, soft active shell is the primary convergence bottleneck.

The precision evidence does not support a gross native-force sign or unit
mistake. A native 0 K energy matches the known 300 K endpoint evaluation within
5.8e-8 kcal/mol. The force writer sign and scaling were traced to the local
AmberTools sources. At 0.1 Å the independent energy and force curvature checks
agree within 2.03%. At 0.01 Å, however, the directional energy derivative is
only 0.028 kcal mol-1 Å-1 and the acceptance record flags mdout print precision
as limiting component-wise agreement.

A bounded additional ±0.05 Å check was run at the same saved candidate and
mode, using two one-cycle native Sander evaluations only. The energies were
-718.1345066 Eh (+0.05 Å) and -718.1346112 Eh (-0.05 Å), with projected raw
path forces +0.0050344202 and -0.0067279712. The energy curvature was
-0.42322 Eh/Å²; applying the calibrated 3.5710676726 Eh/Å per raw path-force
unit gave -0.42015 Eh/Å², a 0.73% difference. At ±0.1 Å, the persisted force
curvature is -0.42590 Eh/Å²; the energy curvature recomputed from rounded mdout
energies is -0.41334 Eh/Å², consistent with the documented print-resolution
limitation. These tests argue against force precision as the main cause of
the eigensolver failures along this mode. Finite-difference cancellation may
still affect other low-frequency directions. The repeated HVP budget
exhaustion, even at 0.1 Å, points to the global-lowest-mode ARPACK problem over
the large, soft active shell as the main observed failure. Missing
eigensolutions must not be interpreted as positive or negative curvature
evidence.

## Executable verification algorithm

1. Keep the exact 12 Å model, proton inventory, -2 QM charge, 98-atom QM set,
   4,165-atom active shell, topology, and endpoints fixed. First verify the
   direct NEB band has met its recorded force criterion; independently audit
   every image for mapped OG1–C, C–Nleave, HG1–Nleave, key C–O geometry,
   proton inventory, and new contacts. Save an explicit scientific-acceptance
   record separate from the process exit status.
2. Select the highest-energy regular-band image and optimize it with a
   mass-weighted dimer using the existing DFTB3/MM native Sander path and the
   full active shell. Preserve coordinates, mode, force history, input hashes,
   and logs. Require stable chemistry and a small active-shell gradient before
   calling it a candidate stationary point.
3. Before any Hessian-vector products, repeat native force evaluations at the
   center and at ±0.05 and ±0.1 Å along the normalized dimer mode. Check force
   reproducibility, the energy/force sign relation, and curvature agreement;
   use the coarser step if the fine-step directional energy change is at the
   printed-energy resolution.
4. Confirm negative curvature along the dimer mode with a central force
   difference, then run short downhill trajectories in both signs. Minimize
   both endpoints without restraints and require them to recover the specified
   ES and acyl-enzyme states with the same proton transfer and product
   connectivity. This checks reaction connectivity but does not alone prove
   the absence of additional negative modes.
5. For first-order-saddle evidence, use a matrix-free, mass-weighted block
   Lanczos calculation on the **entire** 4,165-atom active shell. Seed the
   block with the dimer mode and independent low-frequency vectors, use a
   preconditioner, require eigenpair residual convergence and repeat the
   lowest eigenpairs at two force-validated FD steps. Set a staged HVP budget
   (for example 100 then 250 products); stop on nonconvergence and report it
   rather than scaling blindly. Accept the saddle only with one reproducible
   negative mode, nonnegative remaining converged low modes, small gradient,
   and successful two-sided endpoint recovery.

The current attempts do not yet satisfy steps 1, 4, or 5. Do not label the
dimer candidate or NEB maximum a validated transition state until those gates
are met.
