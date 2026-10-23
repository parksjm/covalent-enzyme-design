# Viewer correction and TI1 parameter audit — 2026-10-22

## Water coordinates and PyMOL display

The trajectory contains O/H1/H2 in every one of its 100 frames. Actual headless
PyMOL 3.1.0 loaded all 5054 atoms per frame and identified both water O–H bonds.
The old viewer showed a sphere only for water O and omitted water sticks,
leaving its hydrogens invisible. A semicolon in a PML comment caused the text
after it (`failed TI1 diagnostic.`) to be parsed as Python and raise SyntaxError.

Both corrected viewers explicitly show water sticks and avoid the semicolon.
Reloading replaces only their own named objects, preventing state accumulation.
Actual PyMOL checks load each viewer twice: movie 100 states, final model 1 state,
two visible water H, no parser errors, and correct water inventory/bonds in every
state. No atoms were added and no coordinates changed. Source/viewer hashes,
console logs and detailed per-state reports are under
`output_nyl12_viewer_fix_20261022/{before_movie,after_movie,after_final}/`.
The repeatable audit is `scripts/check_nyl12_pymol_view.py`, run with the existing
`/data/jpc/code/miniforge3/envs/nyl12-pymol/bin/python`. Standard-library trace
measured 92% line coverage in the actual PyMOL checks. No dependencies added.

Run `./pull_nyl12_inspection.sh` again and then `@view_optimization.pml` or
`@view_local_repair.pml` from the downloaded diagnostic directory. The original
trajectory remains a failed TI1 diagnostic; the display correction does not
change its scientific status.

## Current reactive parameter inventory

Read-only script: `scripts/audit_nyl12_ti1_parameter_inventory.py`.
Output: `output_nyl12_parameter_validation_20261022/inventory_01/inventory.json`.
Source is the retained attempt05 dimer topology; no parameters were changed.

| Atom/term | Present value | Consequence |
|---|---|---|
| C54 type | `c` | Installed GAFF defines this as sp2 carbonyl C, not tetrahedral TI1 C. |
| O49 type/partial charge | `o`, -1.465499 | Charge was shifted for bookkeeping, not obtained from coupled-state RESP. |
| C54–O49 equilibrium bond | 1.2183 A, k 637.7 | Retains a carbonyl-like target inappropriate as a validated tetrahedral C–O target. |
| OG1–C54 equilibrium bond | 1.3584 A, k 390.8 | Inherited provisional term, not fitted TI1 evidence. |
| Several angles at C54 | Approximately 123 degrees | Retain carbonyl-like geometry rather than a validated TI1 reference. |
| Junction source terms | GAFF analogs and a zero torsional placeholder | Derivation and independent validation remain unfinished. |

The inventory lists 97 atoms across reacting Thr/ligand/water, 5 relevant bonds,
7 angles and 3 junction torsions, with topology SHA-256. Its readiness flags are
explicitly false. Actual run coverage of the inventory script is 89%.
Installed type definitions are in
`/data/jpc/code/ambertools26/venv/dat/leap/parm/gaff.dat`.

The next task is to construct the coupled TI1 reference, obtain suitable QM
data, fit charges and bonded terms, and validate before classical AF3 relaxation.
See `ti1_forcefield_validation_gate_20261021.md`. Neutral free-NH2 Thr-NME QM
calculations cannot certify the NH3+/oxyanion covalent adduct. Initial reference
scope was proposed as peptide-capped Thr–full PA6 L4 plus relay water. The user's
subsequent substrate-terminus clarification is recorded in
`substrate_termini_decision_20261022.md`; it does not authorize capping PA6 or
select between the one- and two-residue protein fragments. Do not infer that
an isolated adduct must be a free minimum: specify reference stabilization and
any constraints explicitly before fitting. No new refinement has been launched.

The user is improving the AF3 restraint fork independently. Leave it unchanged
and wait for their update before new AF3 predictions. All actions, including
delegated source audit and failures, remain in `output_agent_actions/20261021.jsonl`;
the filename is historical and records contain their actual UTC timestamps.
