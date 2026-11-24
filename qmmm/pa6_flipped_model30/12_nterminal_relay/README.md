# N-terminal relay: conditional paths from window033

The PI approved this bounded continuation on November 23. The persistent tmux
session is `nyl12_w033_relay_after_minimum_20261123`. It currently waits for
window033's complete minimum-verification result; preparing or viewing this
directory does not mean that a reaction path has started.

The existing window033 Hessian and all design/WT MD continue unchanged. The
watcher requires a successful terminal full-Hessian process receipt, matching
scientific acceptance and summary, zero resolved negative modes, zero unresolved
signs, the original 12,495-coordinate dimension, and unchanged frozen inputs.
Failure or ambiguity stops dispatch; no QM-only result substitutes for this gate.
The recorded numerical qualifications still apply to a passing result.

## Two branches from one accepted geometry and fixed environment

| Directory | Intended chemistry | Asynchronous composite coordinate |
|---|---|---|
| `01_activation_reverse` | Return from the tetrahedral-like state toward E–S: N-terminal proton returns to OG1 as C–OG1 opens | d(Nalpha,H) − d(OG1,H) + d(C54,OG1) |
| `02_acylation_forward` | Proton delivery from the N terminus to the leaving amide N, with C54–N59 cleavage and carbonyl recovery | d(Nalpha,H) − d(N59,H) + d(C54,N59) |

Both branches use exactly window033's complete coordinates and fixed exterior,
104,922-atom topology, 4,165 movable atoms, 98 real QM atoms, eight link atoms,
QM charge −2 and zero electronic temperature. The QM region contains no water;
this does not test a water-mediated route. DFTB3/3ob parameters, electrostatics,
link treatment and source masks are pinned and checked before native execution.
Earlier window028/032 releases have different fixed outer coordinates and are
not used as interchangeable endpoint energy references.

HG1 (atom13917) is the initial biased proton because its source geometry is
favorable: Nalpha–HG1 1.050 Å, HG1–N59 1.975 Å and Nalpha–HG1–N59 179.4°.
All three original N-terminal protons and the pre-existing substrate N–H are
monitored. Endpoint classification allows a different N-terminal proton to
transfer. No proton-identity wall, fixed OG1–C54 bond length or positional
restraint is added. The sole harmonic bias acts on the stated composite;
its component distances are free to change asynchronously.

Targets increase in at most0.15 Å steps to4.2 Å, stopping earlier at the defined
product-side geometry. The force constant is the established Amber100
kcal/mol/Å² convention. These are driving parameters, not mechanistic predictions.
Each native biased window has at most4,000 optimizer cycles and the declared
1e-4 Eh/bohr gradient /1e-7 Eh energy controls. Incomplete convergence, fatal SCC,
changed exterior, or unintended non-relay proton displacement stops that branch
with evidence preserved. Distance changes above0.35 Å are flagged in the plots.

At a candidate endpoint the entire reaction bias is removed. Up to10,000 native
minimizer cycles and an independent1e-5 Eh/bohr force check are required, together
with retained proton inventory and expected heavy-atom connectivity. A native
process exit does not establish a chemically accepted endpoint.

## Conditional TS searches and interpretation

After an accepted unbiased endpoint, a strict interior maximum of the
bias-removed potential-energy profile seeds a native DL-FIND dimer search.
Its initial direction follows adjacent accepted path geometries on the original
movable atoms. There is no automatic TS launch for a monotonically rising path,
an endpoint maximum or a rejected endpoint. Searches use up to5,000 cycles,
the existing0.05 maximum-step control and zero electronic temperature, followed
by an independent force check. These searches produce **TS candidates**;
first-order curvature classification and downhill connectivity remain required
before a verified TS or barrier is claimed. No new umbrella campaign is queued.

## Inspect progress

- `minimum_gate_watcher/wait.json`: actual wait/dispatch/failure state.
- `00_conditional_queue/job/`: exact persistent command, environment and logs.
- Each branch's `branch_status.json` and `progress.json`: native/scientific state.
- Each branch's `vmd_progress/view_path.vmd`: stripped, mapped endpoint trajectory.
- Each branch's `plots/energies.png` and `plots/distances.png`: physical and biased
  PES curves, tracked proton/bond distances and discontinuity flags.
- `01_coordinate_drive/window_0NN/`: successive immutable native inputs/results.
- `02_unrestrained_endpoint/` and `03_ts_search/`: eligible follow-up calculations.

Before dispatch, the two movies contain only the accepted window033 starting
frame. Every completed window adds its exact restart and refreshes plots/movie.
Real cpptraj export, frame coordinates, atom mapping and periodic-cell encoding
were checked before queuing. Invalid analysis-only preflight output is retained
in `clutter/`; original coordinates and scientific model were never changed.

Reproducible entrypoints are `scripts/prepare_window033_relay_path.py`,
`scripts/watch_window033_relay.py` and `scripts/run_window033_relay_path.py` in
the repository; `source_snapshot/` retains the imported source closure actually
used by this campaign. Direct invocation also enforces the minimum gate and
source hashes. Laptop pulls retain the same `qmmm/pa6_flipped_model30/12_nterminal_relay/`
layout and include the compact viewers/plots/receipts, not raw Hessian columns.
