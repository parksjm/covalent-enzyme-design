# TI1 attempt and local-repair protocol — 2026-10-20

## Scientific scope

The state remains the first tetrahedral intermediate: E-chain Thr1 (full
Thr233; Amber residue 929) is NH3+, OG1 is bonded to substrate carbon, the
substrate C–N bond remains intact, and the relay is neutral H2O. F/G/H beta
termini remain NH2. Product registers remain PA6 L2 and poly(X,AMA) MAXA.
This work changes execution and screening, not the Hamiltonian or charge model.

## Attempt interface

Use `scripts/run_nyl12_ti1_qmmm.py` through `scripts/run_nyl12_job.py` with
absolute paths. Required arguments are `--topology`, `--coordinates`,
`--reference`, `--output`, `--amber-bin`, and `--target`. A positive
`--max-cycles` controls the budget. `--output` must not exist. A build followed
by `--run-qmmm` is intentionally rejected before creating a model: a built
water seed has not passed screening.

The audited mask is `:929,1313,1314|(:930&!@C,O)`, formal charge zero, with
one CA930–C930 cut. The force-field partial sum need not equal that formal
charge. DFTB3, GB8, SCF tolerance 1e-8, RMS gradient target 0.05, and the
existing backbone restraint mask/weight (5.0) are retained. Keeping these
settings does not certify their physical accuracy or the embedding boundary.

Each attempt retains copies of topology, starting coordinates, the explicit
restraint reference, calculation sources, and uv environment specification.
It records input/source/engine hashes and installed Python package versions.
DFTB data are hashed and copied before an engine launch. The tmux wrapper
records its own terminal status separately from `process.json` (engine status)
and `acceptance.json` (scientific screens). Preflight rejection has no engine
start/exit code and must never be described as a failed DFTB optimization.

For continuation, pass the exact prior `ti1_dftb3.rst7` as `--coordinates`
and its directory as `--checkpoint-attempt`. The prior attempt must have a
terminal process record and a matching recorded restart digest. Topology,
restraint reference, target, method specification, Sander executable, and
DFTB parameter data must match. The cycle budget may change. A vanished job
without a terminal record remains unresolved and is not automatically resumed.
Historical restarts without these records are not silently adopted.

The preflight checks coordinate shape/finiteness, the TI1 and nonreacting
beta-N proton inventories, reactive connectivity, reactive bond lengths,
water proton bonds, water contacts against the full system, and distances
from water O to beta N/OG1. Water uses the existing element-specific geometric
limits with at most 0.10 A shortfall; N/OG1–water O must be 2.5–3.5 A. Reactive
carbon bonds are screened at 1.20–1.80 A, beta N–H at 0.80–1.25 A, and water
O–H at 0.80–1.20 A. These deliberately broad ranges are rejection screens,
not equilibrium geometry targets. Postflight also requires finite final Amber
`FINAL RESULTS` RMS at or below 0.05 and no recognized SCF/engine failure.

These checks do not prove stereochemistry, bond orders, full-pocket absence
of clashes, full parameter quality, or an unrestrained stationary point.
Product identity is inherited from the previously product-gated topology,
not re-derived from a PDB or a matching atom count. Successful restrained
minimization must still receive a complete contact/geometry review and a
restraint-release study before mechanistic interpretation.

## Local geometry repair strategy

The existing rigid-pocket water search failed for both targets. Its best
candidates include Asn186 side-chain and backbone heavy-atom crowding;
rotating only beta-N hydrogens cannot remove those obstacles. Do not relax
contact thresholds merely to get the current input into Sander.

1. Preserve each failed water seed and its input/source hashes. Use the
   hydrogenated canonical inspection models to identify the offending atoms
   and their bond graph before choosing flexible coordinates.
2. First explore Asn186 side-chain chi1/chi2 changes and rotation of the beta
   NH3+ hydrogen group around its N–CA axis, rebuilding hydrogens consistently.
   Preserve bond lengths/angles, Thr stereochemistry, substrate register, and
   the TI1 proton inventory. Re-run the existing water placement/orientation
   search for each explicit conformer. Record all rotamers and failures.
3. Screen every moved side chain against its full environment, excluding only
   topology-defined bonded/angle partners, and re-audit the entire ligand
   pocket. A water improvement that introduces another severe clash fails.
   Compare candidates to the same original heavy-atom coordinates and report
   per-residue displacement; a low geometric clash score is not an energy.
4. If backbone O still excludes all candidates, side-chain sampling is
   insufficient. Define a separate bounded backbone-repair experiment with
   explicit displaced residues and restrained coordinates. This requires a
   documented restraint specification and separate provenance; it is not a
   checkpoint continuation of the fixed-reference problem. Do not use the
   provisional classical reactive parameters to rank such structures.
5. Only screened seeds enter the unchanged DFTB3 protocol. Retain process,
   convergence, contact, protonation, and stereochemical assessments separately.
   Name outcomes restrained starting models or restrained minima according
   to their evidence. Subsequent restraint release and method/parameter checks
   are prerequisites for stationary-point or mechanistic claims.

No rotamer/backbone-repair implementation or successful repair is claimed by
this protocol. Gaussian/RESP neutral-N template fitting, the separate NH2 ES
state, and full bonded/charge validation remain subsequent work. The Gaussian
installation preflight is not a fitted-charge result.
