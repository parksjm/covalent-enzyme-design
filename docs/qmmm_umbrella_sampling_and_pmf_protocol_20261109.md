# Nyl12 QM/MM umbrella-sampling and PMF protocol

## Matched pathway comparisons and publication acceptance (2026-11-26)

This section governs comparisons of the direct proton-transfer mechanism and
the N-terminal proton-relay mechanism. Apply one versioned preparation,
sampling, postprocessing and acceptance protocol to both. Earlier PMFs remain
historical controls until their compatibility is audited; differences between
unmatched calculations are not evidence for a mechanistic preference.

The initial source-hashed audit is
`provenance/pmf_pathway_comparison_audit_20261126.json`. Recorded native
production inputs match in topology, QM mask/charge, DFTB3, electronic
temperature, embedding, link treatment, 300 K Langevin bath, timestep and
movable-mask settings. SCC tolerance is 1e-12 in the relay input and 1e-10
in the legacy direct-seeded input; the former also requests tight population
convergence. The umbrella CV combinations and strengths differ (Amber
rk = 25 versus rk = 100 kcal mol⁻¹ Å⁻²). Correct reweighting removes the biases; equal rk is not a
mathematical requirement, but convergence and channel checks are mandatory.
The follow-up physical-model audit verified identical active/QM atom selections
and boxes; frozen-exterior positions differ by at most 6.32 × 10⁻⁶ Å, consistent
with ASCII/NetCDF coordinate precision. Current 3ob parameter hashes match the
relay record, but the historical direct campaign did not pin their hashes;
direct-era parameter identity remains unverified. Representative native stage
controls were compared, rather than certifying every historical input. These
qualifications remain explicit. No new MD was launched by this audit.

**Terminology.** Burn-in is the initial part of a continuous trajectory omitted
from equilibrium statistical estimation while its distribution adjusts to the
starting structure and applied bias. It is an analysis exclusion, not deletion
of native data, a repeated restart procedure or proof that subsequent frames
are equilibrated. In the current relay calculation, 1 ps heating and 5 ps
dedicated equilibration were excluded first; the additional 5 ps production
discard leaves 20 ps analyzed out of a 25 ps production trajectory. Restarted
5 ps production files are concatenated before applying that single discard.

**Controls that must match:** enzyme/substrate identity and atom mapping,
proton inventory and modeled protonation microstate, total/QM charges,
force-field and QM method/parameters, QM partition/link treatment, physical
temperature, embedding/electrostatics, boundary/movable-region treatment,
integration settings, bath protocol and physical endpoint-basin definitions.
Changing any of these defines a separate controlled comparison. The 300 K and
348.15 K branches must not be pooled. Microstate-population and conformational
preorganization penalties are separate from conditional chemical barriers.

**Route-dependent choices that must be explicit:** geometrical reference path,
active coordinate combinations, umbrella centers and necessary bridge windows.
Use consistent physical coordinate definitions, units, cross-bias conventions
and bias-strength policy. A scientifically different path need not have the
same numerical window centers or window count. Validate energies and coordinate
mapping independently. Analyze actual proton-transfer behavior; a direct-seeded
trajectory sampling the N-terminal relay is not a pure direct-mechanism sample.

**Common current analysis specification:**

- Exclude heating/dedicated equilibration, concatenate accepted continuous
  production segments per window, then discard an additional fixed 5 ps once.
- Fit MBAR to every subsequent frame using actual per-window counts, identical
  bias conventions, a documented common physical endpoint reference and
  consistent density/bin conventions. Do not assign IID analytical errors to
  correlated frames or remove configurations based on desired chemistry.
- Use 200 within-window contiguous 5 ps block-bootstrap refits and pointwise
  95% intervals, retaining missing-reference/support and disconnected draws.
  Evaluate 2.5 ps blocks, production discards of 0/2.5/10 ps and temporal halves
  in both pathways. These block sizes/discards remain provisional.
- Compare matched sampled horizons as a diagnostic and final independently
  qualified profiles as the scientific result. Different windows may need
  different sampling lengths to meet the same convergence criteria. Equal
  runtime does not establish equal precision or adequate equilibration.
- If evidence requires changing burn-in, block length, binning, reference or
  another estimator setting, version the revision and reanalyze every pathway
  in the comparison. Do not tune one pathway to obtain a desired barrier.

**Publication gate:** require connected and adequate local overlap in physical
CV space, resolved/bracketed E–S and product basins, temporal and orthogonal
stability, correlation-aware uncertainty and independent-initialization checks.
Agreement of PMFs from extended sampling to approximately 1 kcal/mol in the
barrier region is a prospective practical target, not a universal convergence
theorem. A narrow bootstrap interval cannot establish these other requirements.
Use committor/downhill connectivity checks before identifying a transition-state
ensemble. Preserve method/model uncertainty separately from sampling error.

**Direct-route reanalysis, 26 November:** the archived 47.05 kcal/mol direct
profile is a marginal density along q_acyl, whereas the relay profile is a
three-CV path density. Matching only their discard and error estimator is
insufficient. Reconstruct q_act and q_del from the native direct trajectories,
verify q_H = q_act + q_del and the original cross-bias arrays, and apply the
same equal-Angstrom three-CV projection with extended terminal tangents.
Retain each route's original physical reference path: the direct reference
consists of its 24 source-hashed NEB images. Normalized s labels are therefore
route-specific and are not identical molecular states across panels. Use
the same 0.025 bin grid and fixed s = 0.0125 E–S reference; extend outer edges
on that grid when necessary to retain terminal fluctuations, without changing
interior bins, clipping or chemical censoring. Analyze the matched 25 ps
horizon and the full available 30 ps direct horizon with the shared estimator.
Use actual native biases, not substituted relay biases. Report direct-seeded
N-terminal proton visits explicitly. These projections alone cannot yield a
mechanism-specific rate ratio or remove unequilibrated orthogonal sampling.

**Current priority:** finish the already authorized 60 ps relay sampling;
inspect slow leaving-N–HG1 behavior in windows 025/026 and resolve endpoint
bracketing before adding windows or launching independent-start tests. Audit
the direct campaign against these controls before commissioning a matched
rerun or presenting a direct-versus-relay free-energy difference. No new
mechanistic branch or altered Hamiltonian is authorized by this protocol note.

Retain the exact inputs, environments/lockfiles, random seeds, source hashes,
accepted trajectory lineage, analysis requests, all failed/superseded attempts,
uncertainty replicates and plotting source. The code checkpoint and curated
evidence are versioned in Git; native trajectories and full source records
remain on the workstations. A publication archive must provide the required
inputs and analysis data, with its persistent identifier recorded in the
manuscript; a Git commit alone does not reproduce unarchived trajectory data.

## Execution priority: full-path coverage first (2026-11-25)

The PI has prioritized completing the current acylation US/MBAR profile.
Accepted unbiased E–S and acylenzyme endpoints and their source-bound driven
connection permit preparation/sampling while PES stationary-point verification
continues independently. Do not make full-Hessian TS certification or a stable
TI1 well a prerequisite for umbrella dynamics. Native CV/bias correctness,
chemical mapping, source integrity and sampling acceptance remain mandatory.

Use activation `(q_acyl, q_act)` and delivery `(q_acyl, q_delivery)` charts over
the whole path, including both released endpoints. Densify centers to at most
0.15 Å in each biased axis, including the coarse activation proton-transfer
bracket. The 18-window near-relay grid is an input/compatibility fixture, not
the full campaign. First obtain one accepted 5 ps production block in every
window, then extend least-sampled windows in retained 5 ps blocks toward at
least30 ps/window and the eight convergence gates below. Analyze incomplete
coverage diagnostically; do not invent a continuous PMF across disconnected
components. Use all three proton/acyl features for cross-bias evaluation and
one shared MBAR normalization across the chart families.

## Current decision: N-terminal-relay connectivity before a new PMF (2026-11-24)

The PI authorized further structural evaluation and E–S → relay → acylenzyme
connectivity work, with the requirement for a credible, rigorously obtained
barrier. A value near an expected kinetic range is not an acceptance test.
The historical 35.81 kcal/mol PES path and approximately 47.05 kcal/mol apparent
PMF are retained negative results; this authorization does not restart that
campaign. The old campaign status sections below remain historical evidence.

**Connectivity gate.** Join the reverse-activation branch in reverse order to
the forward-acylation branch at the identical window033 restart. Verify every
accepted predecessor/restart hash, atom and proton inventory, topology, QM/MM
Hamiltonian, frozen exterior and energy reference. Require the released E–S
and acylenzyme endpoints to pass unbiased force/chemical retention checks.
Inspect all four tracked protons and the remaining QM N/O hydrogen pool, the
C–OG1/C–Nleave/C–O distances, oxyanion-hole contacts and discontinuities. A
reproducible geometric connection does not establish a minimum-energy path or
TS. Window033 full-minimum status stays unresolved; an obligatory TI1 minimum
is not a prerequisite for proceeding.

**Coordinate decision, after observed paths.** The candidate heavy-atom
coordinate is q_acyl = r(C,Nleave) − r(C,OG1). Resolve activation and delivery
separately when they exhibit distinct slow modes. Existing distance-difference
coordinates, q_act = r(OG1,H) − r(Nalpha,H) and
q_delivery = r(Nalpha,H) − r(Nleave,H), can support two local 2D patches with
q_acyl, provided proton identity remains adequate. Their joint distributions
must overlap in a common relay region before considering patch stitching.
Separate one-dimensional PMFs cannot be added using independently chosen zeros.
If the relay is not a basin, or the local charts exhibit hysteresis, use a
single connected multidimensional/path-coordinate description of the full
reaction rather than inventing two elementary rate constants.

Exchangeable hydrogens are chemically interchangeable. Reuse the implemented
smooth all-QM-N/O hydrogen-pool coordination diagnostics: CN(Nalpha), CN(OG1)
and CN(Nleave), plus CN(Nleave)−CN(OG1). These identify hidden proton relays
and label swaps; a label swap is not automatically an invalid configuration.
These quantities are presently analysis observables. Do not claim that Amber's
distance-restraint interface biases nonlinear coordination functions. Any new
coordination-number bias requires an explicitly verified energy/analytic-force
implementation and correct MBAR cross-bias evaluation before sampling. Do not
add arbitrary proton-identity walls or fix C–OG1 at 1.55 Å to obtain a result.

**Sampling and interpretation gate.** Retain all eight numerical convergence
gates below: overlap/connectivity, correlation-aware effective sample size,
block/independent-start stability, confidence intervals, estimator agreement,
and stationary orthogonal chemical/environment distributions. Thirty ps is an
allocation, not proof of convergence. Start with a bounded overlap/equilibration
pilot after connectivity; extend sampling based on these diagnostics. Examine
Nalpha coordination, both Asp oxygens, oxyanion-hole contacts, Thr orientation,
substrate register and relevant waters as possible slow coordinates. Plan a
small staged committor test near any candidate PMF saddle, with unbiased
shooting, independently drawn thermal velocities, explicit basin definitions
and unresolved trajectories retained. A PMF maximum alone does not prove that
the chosen collective variables identify a transition-state ensemble.

Before a quantitative mechanism claim, compare sparse matched embedded DFT/MM
energies/forces for new E–S, relay/TS-region and product structures; use exactly
the same partition, electrostatics and exterior. A free-energy correction needs
overlap of the reference/target energy-difference distributions; arbitrary
energy offsets are not corrections. Report method sensitivity if reweighting
is not supported. A true PES TS additionally needs first-order curvature and
downhill connectivity. Report protonation-conditional barriers explicitly;
include microstate-population costs before claiming an overall assay barrier.
The same requirement applies to conformational preorganization: a PMF started
from a rare relay-ready E–S pose is conditional on that pose unless the sampled
reactant ensemble connects to the dominant bound conformations. Monitor
Nalpha–Nleave distance/orientation and substrate register as candidate slow
coordinates. A population penalty −RT ln(p_ready) can only be used with an
adequately sampled, consistently defined bound-state probability; representative
snapshots or Tier occupancy alone do not supply that probability.

**Temperature/state provenance.** Keep current 300 K comparisons unchanged.
The authorized 348.15 K WT branch is prepared, not launched, for full-system
thermal adaptation/review. The 130,815-atom four-substrate design MD system and
104,922-atom one-substrate legacy QM/MM system are different systems. Any new
assay-temperature mechanistic state requires explicit preparation/mapping and
bath equilibration; it cannot be obtained by heating only the movable shell.
The pH-adjustment temperature of the 75 °C Tris assay is unknown. Neither the
assay's conversion nor heating a simulation determines a unique activation
free energy. A high but converged barrier remains a valid negative result and
must trigger an explicit mechanism/model decision rather than relaxed gates.

## Status and interpretation of the initial campaigns

The 20-window direct-proton-transfer and 25-window asynchronous E--S-to-TI1
campaigns use 1 ps heating, 1 ps biased equilibration and 2 ps biased
production per window. These short campaigns are engine, stability and overlap
pilots. They are not production PMFs and must not be used to report a barrier.

## Estimator decision

Use MBAR as the primary estimator and binned WHAM as an independent
cross-check. MBAR uses all cross-window reduced bias energies without choosing
histogram bins and provides an overlap matrix, effective sample information and
uncertainty estimates. WHAM remains useful because it is conventional for
enzyme QM/MM umbrella sampling and offers an implementation-independent check,
but its result must be tested against bin width and convergence tolerance.

Do not interpret agreement between estimators as proof of convergence. Sampling
duration is adaptive: 30 ps/window is the present literature-aligned allocation,
not the definition of convergence. Analyze cumulative retained blocks after
automated equilibration detection and require the following preregistered project
gates in *two successive* 5 ps additions before PMF promotion:

1. continuous neighboring-window overlap with no unsampled gap and adjacent
   histogram-intersection overlap of at least 0.20;
2. integrated-autocorrelation/statistical-inefficiency analysis for every
   window and decorrelation or an equivalent correlated-data treatment, with
   at least 100 effective uncorrelated samples per window for the biased
   coordinate after discarded equilibration;
3. a connected MBAR overlap matrix with no isolated state or bottleneck that
   disappears under modest resampling;
4. changes of no more than 0.5 kcal/mol in the barrier and basin free-energy
   differences and no more than 0.10 A in their coordinate locations between
   successive cumulative analyses;
5. a contiguous-block-bootstrap 95% confidence-interval half-width of no more
   than 1.0 kcal/mol for the barrier and basin free-energy differences;
6. MBAR and bin-width-sensitive WHAM profiles that agree within 1.0 kcal/mol
   over the thermodynamically relevant coordinate range;
7. independent-start or forward/reverse estimates whose barrier and basin
   differences agree within 1.0 kcal/mol or their joint 95% uncertainty; and
8. stationary component-distance and important orthogonal-coordinate
   distributions without unresolved drift, hysteresis or hidden multimodality.

These numerical values are declared Nyl12 project acceptance thresholds, not
universal constants. Failure of an overlap gate calls for new/intermediate
windows; failure of effective-sample-size or uncertainty/stability gates calls
for longer sampling; failure caused by hidden multimodality calls for a better
reaction coordinate rather than merely longer one-dimensional trajectories.

The asynchronous scalar coordinate

`xi = d(OG1-HG1) - d(N-HG1) - d(OG1-Ccarbonyl)`

can map distinct combinations of proton transfer and nucleophilic attack onto
the same value. It is also blind to proton identity: the 2026-11-09 extension
satisfied the HG1-based coordinate by transferring a different ammonium proton
(Thr929 H1) back to OG1. If conditional component-distance distributions are
multimodal, show hysteresis or permit proton exchange, do not solve that problem
merely by extending the 1D run. Use independent proton-transfer and attack
coordinates with proton-ownership diagnostics, followed by a targeted 2D
umbrella patch or a path/string coordinate in the affected region.

## Sampling plan

### NEB-seeded campaign: coverage-first scheduling (2026-11-20)

The PI authorized reprioritizing the existing 63-window, 32-worker campaign:
finish each active 5 ps native block, then run the queued windows through
their unchanged heating, equilibration and first 5 ps production block.
Extensions are held until all 63 windows have accepted, diagnostic-ready
first blocks. Thereafter, each worker extends its least-sampled assigned
window one block at a time until every window reaches at least 30 ps.
Host ownership, the two CV biases, temperature, proton inventory, Hamiltonian,
movable region, native stage inputs and Langevin seeds are unchanged.

The historical runner has no cooperative pause/resume interface. Its queue
controllers are replaced using verified process identities. Only the Python
window driver is suspended while its child Sander finishes normally; Sander
is never signaled. A separate controller records the native exit status,
preserves the original process metadata, validates the completed checkpoint
using the original acceptance gates, and starts the next unattempted stage.
Completed stages are never rerun or overwritten. Native coordinates,
velocities, trajectories and existing source/input snapshots are retained.
The scheduling-only source snapshot, process identities, actions and release
receipts remain under the existing campaign's
`dispatch/coverage_first_20261120/` directory.

Complete first-block coverage permits an earlier overlap assessment; it does
not establish connected overlap, a converged PMF or a viable mechanism.
The existing MBAR connectivity and sampling-convergence gates still apply.

### Earlier direct-proton-transfer and asynchronous campaigns

Retain the current pilots. After their overlap and structural diagnostics pass:

1. continue each accepted window from its restart for 5 ps of additional
   biased equilibration;
2. accumulate production in separately retained 5 ps blocks;
3. analyze every cumulative 5 ps increment rather than only selected nominal
   durations;
4. treat 30 ps/window as the initial publication-quality allocation, excluding
   heating and discarded equilibration, while applying the gates above; and
5. extend selected or all windows beyond 30 ps when overlap, effective sample
   size, uncertainty, PMF stability or independent-start tests fail.

With one campaign using violet, the initial measured rate was approximately
64 minutes per simulated picosecond per window. When all 45 direct and
asynchronous windows shared violet, the observed rate slowed to approximately
1.4--1.9 wall-hours per simulated picosecond. Consequently, one 5 ps block is
expected to require approximately 5.3 hours on a dedicated host or 7--10 hours
under the present shared load. Production remains segmented into retained
blocks for restart safety and convergence analysis, but sampling continues to
at least 30 ps/window.

The direct-PT extension is assigned to violet and the asynchronous extension
to azalea. Each host uses its native AmberTools 26 installation; architecture-
dependent executable hashes are recorded rather than treated as an error. The
topology and scientific input hashes must match their declared sources.

## Literature basis

- Pan et al., *J. Chem. Theory Comput.* **2023**, 19, 8234--8246,
  DOI: [10.1021/acs.jctc.3c00973](https://doi.org/10.1021/acs.jctc.3c00973).
  Their QM/MM protocol used 80 windows, Hamiltonian replica exchange, 50 ps per
  window and the final 30 ps for MBAR/TI analysis.
- Paz et al., *J. Comput. Chem.* **2025**, DOI:
  [10.1002/jcc.70109](https://doi.org/10.1002/jcc.70109). Their DFTB3/MM tests
  used replicated profiles and 10--40 ps production per window depending on
  model size, analyzed by WHAM.
- Maschio et al., *Chem. Sci.* **2024**, DOI:
  [10.1039/D4SC02908A](https://doi.org/10.1039/D4SC02908A). Their SCC-DFTB
  enzyme protocol used 2 ps per window only with ten independent repeats, or
  20 ps aggregate sampling per window; it does not support a single 2 ps
  production trajectory as a converged PMF.
- Zheng et al., *Commun. Chem.* **2026**, DOI:
  [10.1038/s42004-026-01888-w](https://doi.org/10.1038/s42004-026-01888-w).
  This polymer-hydrolase QM/MM study used 30 ps per window, analyzed the final
  15 ps with WHAM, bootstrapped the PMF and performed convergence and committor
  checks.
- Andjelkovic et al., *J. Chem. Inf. Model.* **2023**, 63, 5676--5688,
  DOI: [10.1021/acs.jcim.3c00900](https://doi.org/10.1021/acs.jcim.3c00900).
  Their adaptive-string protocol used neighbor Hamiltonian exchange and then
  path-coordinate umbrella sampling with WHAM. Their reported 2 ps criterion
  concerns string convergence, not evidence that 2 ps/window is sufficient
  production sampling for this Nyl12 PMF.
- Chodera, *J. Chem. Theory Comput.* **2016**, 12, 1799--1805,
  DOI: [10.1021/acs.jctc.5b00784](https://doi.org/10.1021/acs.jctc.5b00784).
  The automated equilibration criterion selects the retained region that
  maximizes the number of effectively uncorrelated samples.
- Feddersen and Biggin, *J. Chem. Theory Comput.* **2026**, 22, 5299--5310,
  DOI: [10.1021/acs.jctc.6c00246](https://doi.org/10.1021/acs.jctc.6c00246).
  Their workflow benchmarks emphasize direct convergence/error analysis and
  warn that apparently smooth PMFs can remain undersampled.

## Analysis and visualization

Use `pymbar==4.0.3`, pinned through `uv`, for MBAR analysis rather than a
project-local reimplementation. Treat PMF stability, effective sample sizes,
autocorrelation and adjacent-window overlap as scientific acceptance criteria;
the presence of a reconstructed PMF alone is insufficient.

`scripts/analyze_umbrella_overlap.py` accepts either one production-stage
directory or a stage glob such as `*production*`. It concatenates only stages
whose `process.json` records a successful exit. This permits cumulative overlap
analysis across retained 5 ps blocks without merging or modifying source
trajectories.

Use `scripts/export_umbrella_vmd.py` to create an immutable, stripped VMD
bundle. It samples frames evenly from each completed window in campaign order,
then concatenates the windows in reaction-coordinate order. The output includes
source hashes, a frame manifest, stripped topology, DCD and VMD script. The VMD
display follows the coordinate-driving convention: the QM region is licorice,
only MM substrate carbons are orange, reactive bonds are dynamic, and the
oxyanion-hole distances and Burgi--Dunitz angle are displayed. A sequential
umbrella trajectory illustrates each window's sampled ensemble; it is not a
continuous dynamical reaction movie.
