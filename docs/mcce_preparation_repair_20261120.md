# MCCE preparation repair and complementary pKa pilot

The PI requested complementary pKa methods and repair of MCCE preparation on
November 20, 2026. The complete Nyl12 Tier A model was used before attempting
incomplete NylC crystal models. PROPKA and MCCE therefore use the same starting
heavy-atom coordinates and physical-chain identities; no protonation state is
assigned from either prediction automatically.

## Corrected diagnosis

The previous Nyl12 writer right-justified PDB atom names (`   N`, `  CA`).
MCCE's native parser compares the four-character names literally with ` N  `,
` C  ` and template names such as ` CA `. This caused zero detected termini
and thousands of atom/torsion lookup errors. The detector's `ddvv < 3.0` test
uses **squared distance**, corresponding to approximately 1.732 Å, not 3 Å.
No MCCE source, geometric cutoff or charge model was patched.

The new preparation changes only columns 13–16. All 9,928 heavy atoms,
coordinates, residue/chain numbering, and eight OXT atoms are retained. Native
step1 completed with **8 NTR and 8 CTR**, no spurious termini, zero error
messages, and zero heavy-atom displacement. The earlier OXT removal was an
unsuccessful workaround; valid terminal oxygens are now retained.

Entry points:

- `scripts/prepare_mcce_input.py`: source-hashed atom-name normalization.
- `scripts/audit_mcce_preparation.py`: exact terminal, atom/coordinate and
  native-process checks, failing on missing/extra groups or native errors.
- Focused validation: 19 tests passed, 86–87% module coverage; Ruff passed.

All source/input hashes, raw output and structural acceptance are beneath
`analysis/nyl12_mcce_step1_preflight_20261116/canonical_atom_names_repair_20261120/`.
Previous failed inputs/results remain unchanged and must not be interpreted as
accepted pKas.

## Native warning audit

The conformer stage emitted 227 sparse-charge/default notices. Every one was
checked against the pinned official FTPL files and matched an **explicit zero**;
no numerical parameters were substituted. The complete lookup/source-hash audit
is `charge_defaults_audit.json`. A native generated zero-amplitude ammonium
H–N–CA–C torsion is retained as an explicit model assumption for review. These
are not fitted Amber/RESP parameters and do not imply that methods share an
identical electrostatic model.

## Complementary pilot

After structural acceptance, one protein-only Nyl12 model proceeds through
MCCE4-Alpha commit `08066fd303647a1287296bfae702bf4291e3baab`, level 1 conformers,
DelPhi dielectric 4/80, 0.15 M salt and pH 0–14 in unit increments. Native Monte
Carlo uses six runs, 2,000 iterations per native protocol, 298.15 K, extended
titration sampling and explicit seed 20261120. Complete native configuration
files are retained; completion alone does not demonstrate sampling convergence.

An initial electrostatics launch could not find the already installed DelPhi
because the new runner omitted its bin directory from PATH. That failed attempt
was stopped and preserved in `history/step3_missing_path/`. Successful conformer
generation was retained. The corrected PBE/titration continuation uses an exact
solver PATH, four independent PBE workers with one thread each and nice 10 to
prioritize existing US work. No solver was replaced or dependency installed.

Active tmux session: `nyl12_mcce_pbe_titration_20261120`. Read
`pbe_retry_job/status.json` and `pilot_status.json`. The continuation fails on
native errors, failed PBE conformers or incomplete energy/occupancy tables.
Time budgets are recorded in `pbe_retry_profile.json`; source/scripts and prior
native outputs remain available. A persistent source/binary copy is under the
private mechanism-agent runtime rather than relying on `/tmp` cleanup behavior.

The script produces `mcce_pka_summary.json` and `mcce_pka_sites.csv` when complete.
Out-of-grid pKas remain bounds, not substituted point estimates. Results remain
a preliminary method/ensemble/dielectric-sensitive comparison; no MD/QM state
or reaction path is altered. NylC missing-coordinate/terminal chemistry needs
its own validated preparation and is not repaired by inventing free termini.

## Native pilot completion — 2026-11-20, 21:02 EST

The corrected PBE/titration continuation finished with exit status zero.
All 2,527 head3 entries have native completion flag t, corresponding opp
files exist, no current failed-conformer report is present, and native
step3/step4 completion messages and all four summary artifact hashes match.
Archived missing-DelPhi-PATH diagnostics are earlier failed-attempt evidence.

The beta-site output (canonical beta1=Thr233, beta40=Asp272, beta42=Asp274):

| Site | Chain E | Chain F | Chain G | Chain H |
| --- | ---: | ---: | ---: | ---: |
| Thr233 N terminus | 6.513 | <0 | <0 | <0 |
| Asp272 | <0 | 1.326 | <0 | <0 |
| Asp274 | <0 | <0 | <0 | <0 |

Values below zero are native bounds on the pH 0–14 grid, not exact pKas of
zero. This protein-only, single-preparation pilot does not favor elevated
Asp pKas, but the strong chain variation and unqualified sampling/dielectric/
input-ensemble sensitivity prevent a definitive catalytic assignment.
Native conformer-occupancy stdev maxima over the six MC runs at each pH
range 0.027–0.071; they do not establish pKa convergence. Compare matched
PROPKA results and inspect target titration curves/microstate coupling before
changing any MD/QM protonation. No additional simulation was launched.

Canonical machine-readable results: `mcce_pka_summary.json` and
`mcce_pka_sites.csv` in this note's experiment directory. Native occupancy
and fit diagnostics remain `prepared/fort.38`, `pK.out` and `mc_out`.
