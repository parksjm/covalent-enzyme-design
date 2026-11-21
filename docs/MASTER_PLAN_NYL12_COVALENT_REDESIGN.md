# Master plan: Nyl12 PA6 mechanism and covalent-enzyme redesign

Status: authoritative working plan
Version: 0.2 (2026-11-20; optional design-resource update)
Primary system: mature tetrameric Nyl12 with PA6 L4 in an L4-to-two-L2 register
Secondary system: poly(MAXA), after the PA6 workflow passes

## Purpose and success definition

This document is the controlling scientific plan. Detailed methods, provenance,
failed attempts, and chronological status remain in the linked records below;
they are not duplicated here. A handoff records *where work stopped*. This plan
records *why work is being done, what evidence is required, and when to stop or
pivot*.

The current bounded mechanism execution sequence and manuscript-scope decision
are in [NYL12_MECHANISM_ROADMAP.md](NYL12_MECHANISM_ROADMAP.md). Its completed
state is tracked in the handoff and run manifests, not inferred from this plan.

The program has three linked objectives:

1. Determine a reproducible, chemically meaningful DFTB3/MM free-energy path
   for PA6 hydrolysis by wild-type Nyl12, including both acylation and
   deacylation.
2. Refine important minima and transition structures with embedded DFT/MM and
   verify each transition structure as a true first-order saddle connected to
   the intended adjacent states.
3. Use ensembles of the validated states for symmetry-aware LigandMPNN design,
   followed by Boltz and/or AlphaFold 3 prediction, physical relaxation,
   geometric filtering, multistate scoring, and experimental selection.

Success is not defined as obtaining a visually plausible path, one low-energy
model, or high prediction confidence. Success requires meeting the stage gates
below with preserved negative results. Experimental improvement remains the
ultimate test and cannot be guaranteed computationally.

## Locked facts, scope, and unresolved hypotheses

### Locked facts

- Nyl12 is treated as the mature symmetric alpha/beta tetramer. Catalytic
  residue identities, canonical chain mapping, and the N-terminal Thr233 state
  must remain explicit in every manifest.
- PA6 L4-to-two-L2 cleavage is the chosen modeling hypothesis, not an
  experimentally established precursor/product assignment. Published Nyl12
  products are primarily linear L3; L2 and L4 are also observed. Record the
  modeled substrate graph, termini, register and products explicitly.
- The primary soluble-oligomer model uses the PA6 L4 zwitterion (NH3+/COO-;
  net charge zero), as expected for solvent-exposed amino and carboxyl termini
  near neutral pH. This is a pH- and environment-dependent modeling hypothesis,
  not a claim that termini buried in bulk or crystalline nylon are always
  charged. A neutral NH2/COOH sensitivity state requires separate charges and
  must not be mixed into the zwitterion RESP fit.
- `cif/n12c.cif` is the PI's symmetry-generated tetramer derivative of
  crystallographic 9Z0B. Possible minor edits and additive removal are not
  fully documented; preserve its identity separately from the deposited model.
- Previous WT simulations provide three stable 100-ns PA6 L4 binding replicas
  using AM1-BCC/GAFF-family parameters. They are starting evidence, not proof
  that a new QM/MM Hamiltonian is equilibrated or that chemistry is correct.
- The current Model 23 whole-system DFTB3/MM optimization is a negative-control
  trajectory: its catalytic geometry and Thr233 orientation progressively
  deteriorate. It must not become the production E-S reference merely because
  it eventually converges.
- Contact restraints are generation aids. Passing restraints is not evidence of
  spontaneous binding. Final noncovalent validation must release them and allow
  occupation of any of the four symmetry-equivalent active sites.
- Substrate hydrolysis must not be conflated with Nyl12 autocatalytic cleavage.
- The modular zwitterionic PA6 L4 RESP/GAFF2 release dated 2026-10-31 is
  MD-ready for its declared P6N–P6I–P6I–P6C chemistry; it does not certify
  arbitrary lengths, bulk nylon, TI1 or acylenzyme parameters. Existing
  AM1-BCC results remain a distinct parameter family.

### Mechanistic hypotheses to test, not assume

- Direct transfer from Thr233 OG1 to the N-terminal amine during attack, as in
  the hASNaseIII precedent.
- Asp272-assisted proton transfer, with or without intervening water.
- Stepwise versus concerted nucleophilic attack and proton transfer.
- The identities and roles of Asp272 and the second nearby Asp in substrate
  hydrolysis. Their roles in maturation cannot be transferred automatically.
- Whether ordered water is mechanistically required. The previously proposed
  static OG1/water/N-terminal pseudo-ring is geometrically incompatible with
  the fixed Nyl12 pocket and is not a default assumption.

## Required chemical states

The reference mechanistic study covers:

1. open/apo enzyme and noncovalent E-S ensemble;
2. reactive E-S microstate(s), including declared protonation states;
3. acylation transition state(s), allowing a proton-transfer intermediate if
   it is a genuine minimum;
4. tetrahedral intermediate TI1;
5. acyl enzyme plus the first product;
6. hydrolytic-water/deacylation reactant state;
7. deacylation transition state(s) and TI2 if present;
8. product complex.

The first design round may begin after validated E-S, acylation TS ensemble,
and TI1 ensembles exist. Acyl-enzyme and deacylation states are required before
claiming a balanced catalytic-cycle design. Apo/open structures are required to
avoid designing only an inaccessible closed state.

## Stage gates and measurable acceptance criteria

Metrics are gates unless explicitly labeled diagnostic. Thresholds may be
revised only through a recorded decision with the evidence that motivated it;
they must not be relaxed after seeing a favored candidate.

### G0 — identity, provenance, and chemistry integrity

Required before any production calculation:

- atom graph, atom mapping, total charge, proton inventory, termini, cleavage
  register, products, residue mapping, chain symmetry, topology, coordinates,
  software versions, inputs, environment, and hashes recorded in a manifest;
- no missing/duplicate atoms, nonfinite coordinates, severe unresolved clashes,
  or silent parameter/method substitutions;
- reaction and restraint atoms resolved by stable chemical identities, not only
  transient serial numbers;
- scripted rerun in tmux for persistent calculations, with stdout/stderr,
  restart files, trajectory, process acceptance, and scientific acceptance kept
  separately.

Failure: stop. Do not repair chemistry by adding a fallback or changing charge
without a new declared input state.

### G1 — template-conditioned noncovalent E-S ensemble

Use `n12c.cif` as a protein template. Start with a bounded pilot rather than a
production campaign: at least 12 independent Boltz seeds and, after its current
implementation is validated, at least 12 contact-guided AF3 seeds. Preserve the
mature tetramer and one declared PA6 L4-to-L2 register. Score every equivalent
active site.

Established PA6 near-attack gates:

- d1 (Thr233 OG1 to substrate carbonyl C) <= 3.2 Angstrom;
- d2 and d3 <= 3.8 Angstrom;
- Burgi-Dunitz angle 90–120 degrees;
- absolute projected Herschlag angle 70–110 degrees.

Native-Thr orientation pilot gates after crystal alignment:

- Thr233 heavy-atom RMSD <= 1.0 Angstrom;
- OG1 displacement <= 0.75 Angstrom;
- N displacement <= 1.25 Angstrom;
- absolute chi1 change <= 20 degrees;
- absolute psi change <= 35 degrees.

The CB–OG1–substrate-carbonyl-C angle is reported for every model. The current
105–120-degree interval is **diagnostic**, not a hard gate: hASNaseIII deposited
structures give 135.21 degrees for TS1 and 116.95 degrees for TS2, and do not
define a universal Nyl12 optimum.

Preparation is staged: solvent/hydrogens, local pocket, then short independent
relaxation or MD. Generation contacts and ligand restraints must be released
before acceptance. Promote only if at least three independent models retain
Tier A geometry and native-like Thr233 after release; report the full denominator
and distributions. A single best model is insufficient.

Failure: review template mapping, substrate register/protonation, and relaxation
protocol once. If neither Boltz nor validated AF3 produces reproducible passing
models, run the same minimal pilot on NylB as a lower-complexity control before
deciding whether Nyl12 is presently tractable.

### G2 — QM/MM region and Hamiltonian validation

Define a fixed atom/proton inventory for each chemical path. The initial QM
region includes the reacting PA6 amide, Thr233 backbone and side chain, Phe113
backbone donor, Asn186 side chain, both nearby Asp groups, and only waters
selected by a declared structural rule. Audit every link atom, boundary charge,
exclusion, and water Lennard-Jones treatment. Compare at least two movable-MM
shell radii; do not immobilize every non-QM atom.

Benchmark DFTB3/3OB against omegaB97X-D/def2-SVP embedded or representative
cluster calculations for small, matched models of proton transfer and attack.
Use larger/diffuse-basis single points for charge-separated anionic states where
needed. Preliminary acceptance targets are:

- same qualitative proton localization, bonding state, and state ordering;
- key reactive distances within 0.20 Angstrom and angles within 15 degrees;
- relative stationary-state energies within 5 kcal/mol and barrier differences
  within 5 kcal/mol for the small benchmark set.

These tolerances are deliberately wider than chemical accuracy and are a gate
against qualitative failure, not a validation claim. If DFTB3 changes the
mechanism, loses required minima, or fails these limits systematically, do not
force a production DFTB path. Record a Hamiltonian pivot decision and test a
more suitable supported method/engine.

### G3 — DFTB3/MM paths and free energies

Coordinate driving is for seed generation and diagnosis, not a final barrier.
Use coupled collective variables for bond formation/breaking and proton motion;
include hidden solvent or side-chain reorganization coordinates if sampling
shows they are slow. Prefer string/path optimization followed by umbrella
sampling or an equivalently validated free-energy method.

Path acceptance requires:

- intended atom/proton inventory and chemical connectivity throughout;
- no unexplained adjacent-window jump greater than 0.20 Angstrom in a driven
  reactive distance; any larger jump is investigated rather than smoothed;
- neighboring umbrella overlap and no unsampled gap along each essential CV;
- forward/reverse or independent-seed agreement in state identity and path;
- PMF profiles stable to additional sampling and block analysis to within
  1.0 kcal/mol over the reported barrier region;
- reported bias convention, unbiased PMF, estimator, uncertainty, hysteresis,
  and conditional protonation-state population;
- both acylation and deacylation assessed before naming the catalytic bottleneck.

A PMF maximum is an ensemble bottleneck, not automatically a first-order saddle.

### G4 — DFT/MM stationary-point refinement

Optimize selected minima and path maxima with embedded DFT/MM using
omegaB97X-D/def2-SVP initially, followed by basis/region sensitivity and
larger-basis single points. Use DFT for potential-energy refinement, not routine
free-energy sampling.

A transition state is accepted only if:

- the unrestrained gradient converges under the declared movable region;
- a fresh curvature calculation has exactly one significant imaginary mode;
- the mode is localized on the intended bond/proton rearrangement rather than a
  tail, solvent, or boundary motion;
- downhill following in both mode directions, preferably IRC, reaches the
  intended adjacent minima;
- the result is stable to reasonable Hessian displacement, movable-shell, QM
  region, and basis checks.

Native Sander/DL-FIND P-RFO is the baseline. A DFTB saddle or DFT single point
on one is not a DFT saddle. A fixed-environment result is a saddle only in the
declared subspace.

### G5 — state ensembles and force fields

Generate a consistent PA6 L4 RESP/GAFF2 model if production classical sampling
is used for the new workflow; compare its conformational and binding behavior
with the historical AM1-BCC family rather than silently replacing it. Use the
agreed HF/6-31G* ESP/RESP charge convention, with modern DFT geometries and
declared conformer selection.

State-specific TI1/acyl-enzyme bonded parameters are needed only for classical
sampling of those states. They are not a prerequisite when the reacting core
remains QM or frozen while the MM environment is relaxed. Every force field
must pass isolated-molecule, bonded-geometry, charge, energy, and short-dynamics
tests before enzyme use.

### G6 — multistate symmetric design

Begin with fixed-WT-length, focused, symmetry-tied LigandMPNN design. Backbone
remodeling remains optional. When justified by backbone-placement gates,
prioritize evaluating RFD4-Proteina over RFD3 for a local same-length comparison;
retain existing RFD3 inputs/checkpoints for reproducibility. Consider indels
only in a separately validated mapping pipeline. RFD4 co-designs sequence and
structure, so subsequent LigandMPNN is an optional controlled comparison.

The PI's optional resource list also includes standalone tmol, ProLIF/PLIP,
the existing PA6 RESP/GAFF2 release and AF3 contact-restraint implementation.
See [the source-checked inventory](redesign_tool_inventory_20261120.md) for
roles, local presence, scope limitations and qualification gates. Experts may
recommend a relevant subset; none is a mandatory dependency or new campaign.

Design against ensembles representing open/apo access, E-S, the refined
acylation TS region, TI1 if it is a genuine minimum, and later
acyl-enzyme/deacylation states. Preserve
tetrameric sequence symmetry while allowing instantaneous conformational
asymmetry. Rank by gates followed by a diverse Pareto shortlist, not by a single
weighted score.

Candidate assessment includes the following validity checks and observables.
Chemical identity and structural integrity are hard gates; quantitative
confidence/geometry/occupancy cutoffs require the control calibration below:

- correct mature assembly, chain mapping, catalytic residue identities, and
  substrate graph;
- acceptable fold/assembly confidence calibrated against WT controls;
- recovery of catalytic geometry after physical relaxation;
- no severe effective strain or clashes;
- compatible geometry across the selected catalytic states;
- unrestrained noncovalent E-S recovery at any equivalent site, without contact
  restraints or a covalent bond;
- open-state/access compatibility and no obvious irreversible gate closure;
- reproducibility over independent seeds and no dependence on one lucky model.

AF3/Boltz confidence, sequence likelihood/convergence, geometry, effective
strain, and optional interface-delta scores measure different things. None may
rescue failure of a chemistry or structural-integrity gate. Interface-delta
remains exploratory until calibrated.

**Known-active calibration comes before enlarging the design campaign.** The
PI reports that experimentally active single mutants were rejected by the
current model/MD filters. First quantify known-positive retention, individual
gate failures and replica/pose sensitivity on matched WT/mutant controls.
The short L4 substrate is a tractable local proxy for a likely longer polymer;
chain direction and proximity of its physical termini are unresolved. Do not
require both charged ends as a biological binding assumption. Later interior-
segment versus near-N/near-C comparisons need explicit cap chemistry, mapping
and parameter/charge qualification; current systems remain unchanged.
Chemical identity, valid preparation and catastrophic clashes remain hard
checks; unvalidated geometric bins, Tier A occupancy and model confidence
must not eliminate variants as if they were measured activity. Keep the current
filter version frozen for that assessment; do not widen cutoffs to rescue
preferred candidates or restart PA6 parameterization solely because occupancy
is low. The bounded plan and optional PLM roles are in
[`plm_and_filter_calibration_20261121.md`](plm_and_filter_calibration_20261121.md).
ESMC sequence scores are secondary priors; ESMFold2 is a candidate independent
ligand-aware structure predictor. ESM-IF1 lacks LigandMPNN's explicit ligand
context by default. No new PLM installation or campaign is implied.

### G7 — experimental shortlist

Select a sequence-diverse, mechanism-diverse set rather than the top values of
one score. Include WT and appropriate negative/positive controls, record every
filter denominator, and predeclare assay-detectable products. Computational
success means a traceable shortlist with calibrated risks; improved activity is
established only experimentally.

## Data, provenance, and organization

- Keep this file as the authoritative plan. Keep
  `docs/standalone_handoff_20261020.md` as the current operational checkpoint.
- Store method details and literature synthesis in topic documents; link them
  here rather than copying evolving claims into multiple files.
- Use `qmmm/<system>/<NN_stage_name>/` for reaction work. Each stage contains
  `source/`, `inputs/`, `outputs/`, `trajectories/`, `plots/`, `vmd/`, and
  `provenance/` as needed. Do not create a parallel generic `jobs/` hierarchy.
- Use descriptive names; retain hashes in manifests rather than directory names.
- Keep immutable raw inputs and outputs. Derived plots and stripped trajectories
  must name their source hashes. Move obsolete convenience artifacts to
  `clutter/` for review; do not delete scientific evidence during cleanup.
- Version source, compact inputs, manifests, documentation, and tests in Git.
  Do not commit large trajectories, restart files, or proprietary/unpublished
  structures unless explicitly approved.
- Log every agent/tool action in the append-only journal before and after
  mutations. Record failed commands and permission mode. Persistent calculations
  use the scientific tmux launcher contract with environment and hashes.
- Every campaign has two conclusions: process acceptance (did it run as
  specified?) and scientific acceptance (did it pass the declared gate?).

## Immediate execution order

1. Let the current Model 23 optimization finish without consuming a GPU; archive
   it as a negative-control basin and do not promote it on convergence alone.
2. Finish the template-conditioned E-S preflight and launch the bounded Boltz
   pilot using `n12c.cif`, one declared PA6 L4-to-L2 register, canonical mapping,
   and the established contact set. Keep AF3 as a second arm after the user's
   restraint implementation is validated.
3. Run staged relaxation, release generation contacts, and apply G1 to all
   seeds. Promote at least three independent structures or stop/pivot as stated.
4. On promoted seeds, validate QM-region and movable-shell choices and run the
   DFTB3-versus-DFT small-model benchmark before production free-energy sampling.
5. Generate multidimensional path seeds for the direct and Asp-assisted
   hypotheses, then advance only mechanisms that pass G2/G3.

## Decision register

| Date | Decision | Reason |
|---|---|---|
| 2026-10-27 | Keep Model 23 running as a negative control, not the production E-S structure | At >6700 cycles it is only Tier B and shows large native-relative Thr233 displacement/psi drift. |
| 2026-10-27 | Build a template-conditioned ensemble before considering a Nyl12-to-NylB pivot | One failed seed/global optimization does not establish enzyme intractability. |
| 2026-10-27 | Treat CB–OG1–C as diagnostic pending calibration | The proposed 105–120-degree range is chemically plausible but not a validated universal gate. |
| 2026-10-27 | Use fixed-length symmetry-tied LigandMPNN before RFD3 | It minimizes mapping and maturation risk while providing a clear baseline. |
| 2026-11-20 | Prefer evaluating RFD4 over RFD3 for optional later backbone remodeling; retain the LigandMPNN baseline | PI-selected resource update; Nyl12 symmetry, catalytic conditioning, mapping and runtime remain to qualify. |
| 2026-10-27 | Map split Boltz alpha/beta query pairs to n12c mmCIF label-asym precursor chains A--D | The installed parser exposes label-asym IDs A--D; author IDs A/F/H/V are not valid template identifiers in this interface. |

## Detailed records incorporated by reference

- `docs/standalone_handoff_20261020.md` — live checkpoint and current jobs.
- `docs/template_conditioned_es_plan_20261027.md` — E-S pilot rationale.
- `docs/nyl12_mechanistic_project_scope_20261024.md` — mechanistic scope and
  relation to redesign and NylB.
- `docs/pa6_dftb_validation_pilot_20261024.md` — QM-region and DFTB validation.
- `docs/pa6_ti1_coordinate_driving_plan_20261024.md` — path preparation history.
- `docs/reaction_path_bias_provenance_20261024.md` — bias/energy conventions.
- `docs/ts_search_and_symmetry_decisions_20261020.md` — DL-FIND, saddle
  verification, AF3 restraints, and equivalent-site handling.
- `docs/design_ranking_and_iteration_20261022.md` — multistate ranking and
  LigandMPNN symmetry.
- `docs/redesign_backbone_symmetry_and_access_20261022.md` — RFD3 decision,
  chain-length contract, and substrate access.
- `docs/hasnaseiii_repository_audit_20261026.md` — relevant hASNaseIII evidence
  and limits of transfer to Nyl12.
- `docs/qm_method_and_substrate_checkpoint_20261020.md` — QM method and substrate
  parameterization decisions.
- `docs/action_logging_and_rollback.md` — audit and rollback requirements.
