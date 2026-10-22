# Local water repair checkpoint — 2026-10-21, evening

## Terminal outcome at 22:37 Eastern

The first 100-cycle DFTB3/MM stage completed, then the chemistry gate stopped
the workflow. No calculation remains active. Water's maximum contact shortfall
is 0.029899 A (passes 0.10 A), but OG1–C54 is 2.195122 A (fails 1.2–1.8 A),
both directional H-bond screens fail, and RMS 0.63057 exceeds 0.5. Chirality
was retained. Later pocket/release stages did not run. This is a failed TI1
repair with useful steric cleanup, not an accepted minimum. The diagnostic PDB
and acceptance JSON are in `calculation/01_hydrogens_water/` of attempt05.
Review state-preserving restraints and the DFTB3/starting-geometry assumptions
before another experiment; do not repeat this attempt unchanged.
Running-state descriptions below record the earlier checkpoint.

The user authorized roughly positioned AF3 water -> explicit hydrogens -> local
restrained minimization. They subsequently rejected a geometric oxygen shift
and approved preoptimization with the existing Amber/TIP3P parameters.

## Current physical calculation

Tmux: `nyl12_local_preconditioned_20261021` on violet.
Artifacts: `output_nyl12_local_repair_20261021/pa 6_l 4_alpha_beta/attempt_05_tip3p_preconditioned/`.
Check `job/status.json`, `job/job.log`, `calculation/*/process.json` and
`acceptance.json` before deciding what completed. Absence of a session is not
scientific success. All inputs, source snapshots, engine-data snapshots, exact
commands and hashes are retained. The wrapper continues without the laptop.

The source is the two-restraint AF3 alpha A/beta E/PA6 L4 I/water J prediction.
The exact AF3 oxygen is preserved when adding hydrogens. Extracted source
topology residues 1–232,929–1024,1313–1314 become 330 residues,5054 atoms.
All 2516 AF3 heavy atoms mapped; no terminal heavy atoms needed reconstruction.
Retained Thr 233 is NH3+, OG1 covalently bonds substrate C54, substrate C–N remains
intact, and water is neutral H2O. Nonwater hydrogens are transferred using local
bonded frames; water H orientations are constructed explicitly.

Water remains WAT/O,H1,H2 with OW/HW types. Charges are-0.834/+0.417/+0.417;
O Rmin/2=1.7683 A and epsilon=0.152 kcal/mol. Names do not by themselves assign
parameters: these values were checked and copied from the original topology.

`00_tip3p_water_only`: only the three water atoms move, with the full enzyme and
substrate fixed. GB8, two N/OG1–Ow bounds 2.6–3.2 A with rk 2=rk 3=10,200 cycles,
drms 0.001. A temporary water-first topology satisfies Sander's GB belly rule;
the explicit atom permutation and inverse mapping are retained and atom types,
charges, LJ parameters and GB radii/screens are checked. **Measured nonwater
displacement 0.0 A.** Water O moves 1.408215 A by energy minimization. Final RMS
0.042551 exceeds 0.001, so this is coordinate preconditioning, not a converged
minimum or classical TI1 energy validation.

DFTB3/MM then runs three separately logged stages, using the original prepared
positional reference:100 cycles restraining nonwater heavy atoms at 5;300 cycles
with local residues free and outer heavy atoms restrained at 5;300 cycles with
bridge restraints released. Bridge force constants 10,5,0; drms 0.5,0.05,0.05.
Whole residues within 5 A of Thr 233/ligand/water and beta residue 2 are local
(29 residues here, including local backbones). This permits adjustment rather
than deliberate pocket widening; inspect resulting displacements on completion.
QM mask`:233,329,330|(:234&!@C,O)`,charge 0,109 atoms,364 electrons with link H;
one CA234–C234 boundary. Water is QM, with topology LJ retained for its MM
environment interactions. Numerical controls:dftb_maxiter 200,itrmax 2000,
scfconv 1e-8; electronic temperature unchanged. The first 10 optimization steps
had converged SCC evaluations after the classical preconditioning.

An explicit SCC inaccurate-force warning stops the stage. Inventory, reactive
bond lengths, proton bonds, stereochemical signed volumes, water contacts,
bridge geometry and gradient convergence are assessed separately from exit
status. No stationary point, full-octamer compatibility or parameter validation
is claimed. Local MM reactive parameters remain provisional.

## Retained earlier attempts

- `preparation_01`: preparation-only topology/mapping preflight passed.
- `attempt_01`: shifted oxygen by scripted 1.12 A; stopped at user correction.
- `attempt_02_original_oxygen`: no oxygen shift; SCC forces unconverged; stopped.
- `attempt_03_tip3p_premin`: rejected GB belly mobile-atom order before minimization.
- `attempt_04_tip3p_water_first`: ParmEd AmberParm copy/LJ-table bookkeeping error
  before engine launch. Fixed by converting to parameterized generic Structure
  before reordering; a regression verifies mapping and fixed coordinates.

Tests run in tmux under `verification_*/`; source/geometry checks and mocked
engine tests are distinct from real calculations. Optional real-topology tests
use `NYL12_REPAIR_FIXTURE` pointing to the retained preparation job directory.
All original failures and subsequent corrected attempts remain available.

Active append-only action journal:`output_agent_actions/20261021.jsonl`.
There were approval-review timeouts before command execution and the existing
bwrap RTM_NEWADDR failure; neither is a successful action or safety rejection.

Next: inspect the running job and each stage acceptance; export actual relaxed
coordinates and per-residue displacements; review pocket clashes/protonation
and released bridge before transferring any result to the full octamer. Refresh
the laptop inspection bundle with accepted/diagnostic labels as appropriate.
