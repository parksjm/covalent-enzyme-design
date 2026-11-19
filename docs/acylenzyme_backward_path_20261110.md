# Backward acylation pathway and mechanistic scheme corrections

## Decision and scope

On 2026-11-10 the user authorized constructing an acylenzyme endpoint and
connecting it to E–S with DFTB3/MM NEB. The previous biased TI1-like seed reverted
toward E–S during unbiased DFT/MM release. This rejects that seed as a stable
minimum under that setup; it does not establish that Nyl12 lacks TI1.

DFTB3/MM is the exploratory Hamiltonian because its Amber implementation is
already operational. The benchmark does not validate its barrier heights.
GFN2-xTB is not promoted merely because it matches some attack-region energies:
its proton-transfer errors remain material. DFT is reserved for selected checks.

## Corrections to the supplied catalytic scheme

The attached scheme concerns mature-enzyme substrate hydrolysis, not autocleavage.

1. In the lower-left product complex, neutral carboxylic acid plus Thr–OH and
   Thr–NH3+ requires an additional proton relative to the drawn cycle. Donation
   from Thr–NH3+ to the departing Thr oxygen regenerates Thr–OH and Thr–NH2.
   Draw that neutral state, or explicitly show subsequent proton transfer to
   give carboxylate/Thr–NH3+. Balance solvent proton uptake/release explicitly.
2. Label the deacylation tetrahedral species TI2. Label both TI1 and TI2 as
   proposed species until their stability is demonstrated for this system.
3. Distinguish activation, addition, leaving-group protonation and C–N cleavage;
   a single arrow need not imply they are concerted, nor that each has a minimum.
4. Water is a reactant for deacylation. A relay during acylation is a separate
   hypothesis, not required by the drawing. Failure of the previously tested
   pocket placement does not exclude all water-assisted mechanisms.

Andjelkovic et al. provide a related N-terminal threonine hydrolase precedent,
not a demonstrated Nyl12 mechanism: [JCIM 2023, 63, 5676–5688](https://doi.org/10.1021/acs.jcim.3c00900).
An enzyme QM/MM NEB example is [PMC4388629](https://pmc.ncbi.nlm.nih.gov/articles/PMC4388629/).

## Endpoint construction and acceptance

Retain the original topology, atom order, names, substrate fragments and total
charge. In the acylenzyme seed, form the Thr OG1–substrate carbonyl-C ester,
cleave the substrate C–N bond, and transfer HG1 to the departing amine. Thr's
terminal N retains H1/H2 and is neutral. Retain the amine fragment in the pocket;
do not delete it or change atom count between path endpoints. All changing
bonds are internal to the existing QM region.

Coordinate editing is only a reproducible starting guess. Relax the constructed
endpoint with explicitly recorded temporary chemical restraints, then remove
them and test whether the acylenzyme connectivity survives. No classical
acylenzyme force field is implied by preserving names: classical MD would need
new connectivity, charges and bonded parameters. QM/MM uses the QM Hamiltonian
for those changing bonds.

Use the existing periodic DFTB3/3OB embedding and a consistently defined 12 Å
movable environment; the outer environment remains fixed. This is not a fully
mobile enzyme. Preserve the same fixed coordinates and active-atom membership
across endpoints and images. Do not add an artificial relay water.

## NEB and provisional minima

Use native Amber DL-FIND NEB, without installing another implementation.
Seed the band through existing PT/attack structures and the constructed product,
with consistent atom mapping. A fixed-endpoint band requires accepted endpoints;
if an endpoint has not converged, use freely optimizing endpoints and explicitly
report this limitation, or stop on chemically invalid endpoint release.

The installed native interface evaluates images serially. This is distinct from
concurrent single-CPU endpoint/minimum-release jobs. Record physical potential
energies, not artificial NEB spring energies, in the profile. NEB is a PES method,
not a PMF calculation.

Interior low-energy images of a complete iteration are provisional candidates.
Snapshot them with their iteration and hashes before running independent
reaction-unrestrained optimizations. Retain the 12 Å shell and frozen exterior.
Do not call them intermediates until chemical identity and optimization
convergence pass. A maximum on an unconverged band is not an activation barrier;
eventual TS identification needs saddle-order and connectivity validation.

## Geometry and catalytic-network checks requested during setup

Measurements from the original corrected `attack_09` TI1-like seed (2026-11-10):

| Pair | Distance (Å) |
|---|---:|
| Thr233 terminal N – leaving-group N | 4.937 |
| HG1 – leaving-group N | 4.093 |
| Asp272 nearest carboxylate O – Thr233 N | 2.605 |
| Asp274 nearest carboxylate O – Thr233 N | 2.743 |
| Asp272 nearest carboxylate O – leaving-group N | 5.988 |
| Asp274 nearest carboxylate O – leaving-group N | 6.915 |
| Thr233 OG1 – leaving-group N | 2.448 |

Direct terminal-N to leaving-N proton donation is not preorganized in this
snapshot. It would require conformational change or a relay; an endpoint with
transferred HG1 does not establish how that transfer occurs. The NEB seed uses
the existing terminal-N PT/attack branch and therefore has an initialization
bias. A direct hydroxyl-to-leaving-N route or an Asp/water-assisted route cannot
be excluded merely because this initialized band fails.

The closest water oxygen to Thr N is Amber atom 49306 (residue 11239), at
3.141 Å. No water oxygen is within 5 Å of the leaving nitrogen in this snapshot.
Thus no ready-made one-water bridge is demonstrated here. Existing waters are
MM waters and cannot exchange protons under the current Hamiltonian. A justified
relay-water test must promote the same explicitly identified water into QM in
both endpoints and all images. Do not place an arbitrary water to force a route.

The 2023 NylC study discusses the homologous Asp308–Asp306–Thr267 network for
substrate hydrolysis as well as autocleavage. Its distinction matters: failure
of maturation alone does not isolate a residue's catalytic role in the mature
enzyme. See [NylC crystallographic and mutational study, FEBS Journal, DOI 10.1111/febs.16755](https://doi.org/10.1111/febs.16755).
Both Nyl12 Asp residues are already represented in the QM subsystem. Their
likely activation/electrostatic role is not equivalent to demonstrated direct
acid donation to the distant leaving group. A deprotonated Asp must first
acquire a proton to act as its donor.

The common shell selects 4,620 movable atoms. The fixed exterior is copied from
E–S into every seed (largest reset among old seed frames: 0.01338 Å), and DL-FIND
`refcrd` fixes membership across jobs. A 12 Å shell is an exploratory choice:
check larger-shell sensitivity on promising stationary-point candidates before
accepting structural or energetic conclusions. These releases have no reaction
or positional bias, but are not globally unconstrained optimizations.

## Launched workflow

Host: azalea. Session: `nyl12_acylenzyme_neb_20261110`.
Directory under the frame1603 campaign:
`reaction_paths/03_acylenzyme_neb/`.

1. `01_endpoint_preparation`: DFTB3/MM relaxation with six declared local
   chemical distance restraints, no positional restraints. Initial seed OG1–C
   1.65 Å, C–N 2.30 Å; temporary targets favor ester formation and separation.
   This construction is not a reaction pathway.
2. `02_endpoint_release`: remove all six construction restraints; measure
   ester, carbonyl, C–N and all tracked N–H distances and native convergence.
   Stop the queue if the intended acylenzyme chemistry is lost.
3. `03_neb`: freely optimizing endpoints, 16 regular images plus climbing slot;
   no chemical restraints. Endpoint numerical convergence is recorded separately;
   no unconverged endpoint is frozen and mislabeled a stationary point.

`nyl12_acyl_neb_monitor_20261110` waits for NEB and then updates `plots/`.
Source snapshots and hashes, configs, original topology hashes and process logs
are retained. The prepared model preserves names/indices; it is not a classical
acylenzyme parameter set. Existing source files and old pathways are untouched.

## First endpoint result and running NEB

The reaction-unrestrained endpoint release converged and retained the acylenzyme
on 2026-11-10 (azalea). Distances: OG1–C1.32384Å, C–O1.25576Å,
C–leaving-N4.11101Å, leaving-N–HG1 1.01974Å. ThrN–H1/H2 remain
1.01814/1.02675Å. No Hessian was calculated, so this is a converged local
optimization with retained chemistry, not a proven zero-negative-mode minimum.

Native NEB is now running with zero reaction restraint energy. In addition to
`nyl12_acyl_neb_monitor_20261110` for profiles,
`nyl12_acyl_candidates_20261110` captures coherent full bands and publishes
`vmd_progress/latest/reaction_path.vmd`. Band/energy association is checked using
unbuffered Sander save ordering and stable cross-file signatures. Full source
restarts are retained on azalea for reproducibility; laptop exports omit them.

After cycle5, a dip present in two distinct observed reports is eligible for
an independent unbiased optimization. At most three distinct images are queued
per campaign to bound exploratory work. These are heuristic candidate-selection
rules, not convergence criteria. Each trial keeps the same QM/MM method, fixed
outer environment, reference shell and atom mapping. It records native
convergence and final chemistry without assuming that the state is acylenzyme.

Conservation analysis is in `docs/nyl12_active_site_conservation_20261110.md`.
The analysis uses existing Boltz MSAs, not a newly submitted BLAST search.
The updated laptop pull script refreshes small azalea artifacts through violet
and includes the scheme/procedure and conservation documentation.

## Native NEB restart export defect and repair

A live reactive-distance audit on2026-11-10 found zero distances in the numbered
NetCDF restart exports. The installed Amber `write_custom_rst_file` calls
`write_nc_restart(...,x,x,...)` instead of supplying the coordinate offset
`x(lcrd)`. Internal NEB energies and native active-atom XYZ coordinates are not
written through that defective routine. Ordinary endpoint minimization restarts
use a different routine and were independently measured successfully.

The first VMD exports and three candidate trials (images7,6,9) are REJECTED;
those trials exited with failure, not accepted chemical outcomes. Artifacts and
old monitoring state were retained in run-local
`clutter/native_restart_export_defect/`. No scientific conclusion should use them.

The corrected monitor reconstructs all103254 coordinates from the fixed E–S
reference plus the4620 mapped native XYZ coordinates. It verifies element/order,
nonreactive backbone bonds, report/save ordering and byte-consistent capture.
It retains original XYZ and report files. Independent image10 evaluation from
cycle67 reproduced the native energy exactly at printed precision:
-706.7081487Eh, difference0.000000kcal/mol (acceptance threshold0.01kcal/mol).
This validates coordinate/energy correspondence, not the QM Hamiltonian.
NEB was briefly paused/resumed for that validation snapshot; its inputs and
optimization were otherwise unchanged. No Amber rebuild was performed.

The active candidate monitor is now `nyl12_acyl_candidates_xyz_20261110`, using
`source_xyz_verified/`. VMD updates are every10 cycles, through the same
`vmd_progress/latest/reaction_path.vmd` entry. The periodic energy monitor is
unchanged. The pull script excludes clutter and full source restarts and
checksums the current VMD bundle. Corrected helper/monitor tests:27passed;
helper statement coverage94%.

## Current proton-transfer interpretation (provisional)

The inspected XYZ band shows HG1 moving from OG1 to Thr terminal N in an early
portion, returning to OG1, then transferring to leaving N as acylation occurs.
This is a candidate route involving proton donation from OG1 during addition/
cleavage, not demonstrated direct donation across the long terminal-N/amide-N
distance. The detour through terminal N may be an initialization artifact that
NEB can remove. The path is not converged and the provisional peak remains
about40kcal/mol; it is not an established viable catalytic barrier or mechanism.
No QM water relay is represented in this run.

The user accepted the proposed Lys156/Gly273-backbone QM-region sensitivity
assessment. It is a separate matched-region comparison; the running NEB keeps
its original QM definition. Exact fragment boundaries/charge are being audited
before launching comparisons.

## Follow-up: protonation, alternate substrate orientation, and NEB trend

The direct OG1-HG1-to-leaving-N acylation proposal does not logically require
neutral Thr233 terminal NH2 as the proton acceptor. Whether NH3+ supports this
route must be tested as a distinct protonation microstate with consistent atom
count and total/QM charge. Compare barriers relative to each microstate's E-S,
not absolute energies across systems differing by a proton. The current neutral
model cannot establish the preferred N-terminal protonation state.

Opposite signed Herschlag-phi poses merit a matched structural screen before
another full pathway calculation. Preserve the scissile bond mapping, assess
OG1 attack, donor-H-leaving-N geometry and oxyanion-hole contacts under the same
relaxation/QM protocol. With consistent scissile-amide mapping, opposite catalytic phi signs identify
the flipped local binding poses; other metrics assess catalytic suitability. Test orientation and protonation separately to keep interpretation
manageable; no additional NEB campaign was launched for this assessment.

At this assessment 136 NEB cycles were complete, without convergence. The
latest sampled band/movie was cycle133; its peak was 68.755 kcal/mol above its
E-S endpoint (image10). Earlier sampled peaks were 35.352 at cycle43, 54.384 at
113, 82.044 at123, and 67.563 at130. This supersedes the earlier approximately
40 kcal/mol snapshot: the peak is fluctuating, not steadily decreasing. It is
an unconverged band maximum, not an established activation barrier. VMD exports
update approximately every10 cycles, not every iteration; each current movie
contains16 images along one band. The laptop must pull again to see updates.

## Four-center transition-state precedent (literature check 2026-11-10)

Movilla, Marti and Roca, *Unrevealing the Proteolytic Activity of RgpB
Gingipain from Computational Simulations* (JCIM 2021),
https://doi.org/10.1021/acs.jcim.1c00666, explicitly evaluated a four-membered
O-to-N proton-transfer transition structure. Its reported free-energy cost was
30.9 kcal/mol relative to I1, or 46.5 relative to R1. This is an unfavorable
alternative in a cysteine protease, not evidence for an efficient Nyl12 route.

A PET-hydrolysis DFT model study (JCIM 2021, 61, 764–776),
https://doi.org/10.1021/acs.jcim.0c00797, describes direct proton transfer from
the attacking hydroxyl group to the ester leaving oxygen through a four-membered
cyclic transition structure (RP1), alongside a water-assisted six-membered
alternative (RP2). Primary manuscript: https://par.nsf.gov/servlets/purl/10304064.
This is an ester and model-system analogy, not an amide/Thr acylation demonstration.

These precedents establish that such structures can be investigated and located;
they do not establish a low barrier for our proposed OG1-HG1-N-C route. Reported
free energies above must not be equated with our unconverged NEB potential energies.
Latest plot refreshed here: cycle204, regular-band peak image9, 37.787 kcal/mol
above same-band E-S; no calculation settings changed.

## Mechanism decision checkpoint after cycle263

Latest exported regular-band maximum: 56.427 kcal/mol at image9, cycle263;
native convergence check cycle266 failed all criteria. No climbing image in
export. The rise from37.787 at204 is not evidence of converging barrier height.
Continue current run provisionally as a baseline, but first diagnose localized
forces/image spacing and proton identities near images8–11 before extending or
interpreting it. Do not enable climbing blindly on a poorly relaxed band.

Prioritize screening existing static models for genuinely reversed amide poses
(atom-mapped orientation, not phi sign alone), catalytic filters, native Thr
conformation, and terminal-N/leaving-N proximity. Check whether signed-phi
filters themselves excluded the alternative. Shortlisted poses should receive
unrestrained MM minimization before committing to MD replicas and QM/MM work.
No passing flipped candidate is established by the current records check.

A protonated N-terminus is a separate matched microstate test; direct donation
from NH3+ to leaving N would leave neutral NH2 while OG1 still needs activation.
Keep a full proton ledger, including Asp/water roles. Compare barriers from each
microstate's own E-S, accounting for protonation populations if comparing rates.
Retain Asp272/274 in QM; assess spatially supported relay routes, not arbitrary
proton chains. Lys156/Gly273 expansion remains a controlled sensitivity proposal.

The PET metalloenzyme analogy is peripheral, not support for a Nyl12 mechanism.
Movilla's high-energy four-center route is a warning, not a viable catalytic
precedent. No new campaign or input changes launched at this decision checkpoint.

User clarification: retain OG1–HG1 as leaving-N donor in the NH3+ microstate.
Completed flipped-pose screen: see docs/pa6_flipped_pose_screen_20261110.md.
