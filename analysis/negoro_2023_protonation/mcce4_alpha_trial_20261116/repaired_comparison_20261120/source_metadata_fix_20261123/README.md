# NylC MCCE reruns, 2026-11-23

The fixed normalizer preserves an explicit allowlist of standard PDB metadata
and retains strict atom, residue, alternate-location, chain/terminus and
coordinate validation. Actual accepted inputs required no atom-name changes:
their prepared PDB bytes and SHA256 digests equal the accepted upstream files.
No structures were rebuilt again, no precursor mutations were back-mutated,
and no protonation was reassigned. Unknown PDB records remain errors.

Each case retains previous_current.json, profile.json, preflight/, job/ and
result/ under rerun_20261123/. Profiles fix exact source/native executable
hashes and use the installed MCCE4-Alpha/DelPhi backend unchanged. The source
snapshot is immutable; subsequent job state lives in per-case status files.
Read/run the retained job/run.sh with its declared environment, not a guessed
Python interpreter. Three independent tmux sessions run on Violet.

Scientific acceptance remains separate from process completion. Native steps
1–4, complete conformer electrostatics, finite normalized occupancies and
expected termini must pass before interpreting curves. The existing FTPL
qualification, endpoint-saturated curves, coupled/multisite behavior and chain
asymmetry still limit quantitative pKa claims. No state is assigned automatically.

Mature5XYG has Asp306/Asp308 and Thr267 NTR on chains C/D. Precursor7YU0/7YU1
use continuous chains A/B through355; there is no NTR at267. Existing crystal
PROPKA predictions are reference context; a method comparison requires PROPKA
on the exact same newly prepared structures. The old Nyl12 plotting script
hardcodes different targets/chains and must not be reused unchanged.
