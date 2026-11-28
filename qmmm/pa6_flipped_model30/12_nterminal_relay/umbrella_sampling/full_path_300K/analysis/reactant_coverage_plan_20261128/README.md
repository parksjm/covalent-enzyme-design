# Reactant coverage proposal, 28 November 2026

**Execution update:** the PI approved Geranium and all six histories launched
at 15:56 UTC on 28 November in tmux `nyl12_es_boundary6_20261128`.
The accepted-stage return monitor and separate 106-history MBAR checkpoint
monitor are running. See `preparation/release_receipt_20261128.json` and
`../boundary_pool_20261128/release_receipt_20261128.json`. All native input
hashes and original proposal evidence are retained. The first 1 ps initializes
fresh velocities at 300 K; it is followed by 4 ps continued equilibration and
30 ps production, with no heating ramp.

**Original pre-approval proposal:**

Proposal only; no MD launched and no scientific inputs changed. `plan.json` is authoritative.

The 33-state grid reaches the dominant ES frame, but samples it as an acyl-coordinate tail: q_acyl=-2.0018 is 2.31 harmonic sigmas beyond the closest center -1.75. Of 82,000 retained frames, only two lie within ±0.15 Å of it in all three CVs; both belong to activation_support_3_03. There are 32 within the same two biased CVs. Thus two-coordinate coverage does not establish unrestrained q_del mixing.

At that frame the MBAR denominator responsibilities are 56.25%, 25.33%, 16.28% for support3_03/3_02/3_04; these are state-classification responsibilities, not equilibrium state populations. All33 together account for99.90%. The dominant frame's96.84% target-bin weight is a separate statistic.

A new q_acyl=-2.0 row at q_act=-1,-.8,-.6, each with two geometric starts and fresh distinct velocities, is a bounded ES-coverage test. It preserves q_del unrestrained and therefore cannot guarantee mixing; judge paired-start q_del distributions and block stability. Keep all original33/100 inputs and six new histories separate (106 total).

Geranium's12 prior native5ps receipts give median2.676h per5ps: six concurrent5ps-equil+30ps-production histories project to18.73h wall and112.45CPUh, excluding staging/QC. First10/15ps production checkpoints project to8.03/10.71h. Violet33 have10ps accepted production at pinned baseline; remaining20ps project to27.69h from66 native5ps receipts, median6.923h. These are throughput projections, not guaranteed finish times.

Geranium had no observed scientific CPU/GPU processes and jpc's at queue was empty; sxo/jta have sessions and tmux, so a>few-hour reservation remains unverified. No host-free or launch-permission claim is made.

Checks: exact source MBAR logweights reproduced within1e-5; vectorized native bias energies checked against validated evaluator; Ruff passed. Coverage percentage for this analysis script was not measured; no production code changed.
