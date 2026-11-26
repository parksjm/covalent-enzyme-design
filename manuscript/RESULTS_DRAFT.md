## Reassessment of direct-route TS connectivity (2026-11-23)

Inspection of the installed Amber DFTB implementation showed that 70-cycle
inner SCC/Broyden trial failures can recover within the overall SCF budget.
Our previous log gate had rejected these recovered trials. Corrected,
source-validated history parsing identified five/nine recovered trials in the
two saved downhill branches, without unrecovered errors or temperature changes.
Fresh 0 K native gradients independently accepted both endpoints, with maximum
forces 9.38e-6 and 9.51e-6 Eh/bohr. The minus branch is E-S (OG1-C54=3.123 A,
OG1-HG1=0.981 A, C54-N59=1.360 A); the plus branch is acylenzyme
(OG1-C54=1.329 A, HG1-N59=1.019 A, C54-N59=3.062 A).
The candidate-to-relaxed-E-S potential-energy difference is 46.41 kcal/mol,
and the acylenzyme is 2.71 kcal/mol above this E-S minimum under the same
frozen-state Hamiltonian. These values are not a validated activation free
energy or evidence for a viable catalytic route.

Three independent evaluations at the candidate produced identical gradients;
its maximum force 9.24e-5 Eh/bohr meets the original 1e-4 criterion. First-order
saddle certification remains incomplete: the symmetrized Hessian has one
negative eigenvalue, but the conservative uncertainty screen leaves soft-mode
signs unresolved. The endpoint results establish diagnostic mode-directed
connectivity, not a formal IRC. Original data and prior false-negative gates
are preserved alongside the corrected analysis.

## Matched sequence-mask structural pilot (2026-11-22)

A matched pilot compared four selected sequences from a 19-position mask and
four from a nested 29-position mask, alongside WT, K59D and V111C. Each
sequence was predicted with five seeds and four PA6 L4 ligands per tetramer.
All 55 model-seed endpoints reached the Amber minimization force tolerance;
52 were newly computed and three qualified WT controls were reused. Two
previous WT seeds were regenerated to match the Violet prediction runtime.

| Cohort | Tier A/B sites | Native-Thr sites | Joint-pass sites | Model seeds with ≥1 joint pass |
| --- | --- | --- | --- | --- |
| WT | 10/20 | 5/20 | 1/20 | 1/5 |
| Focused 19-position mask | 24/80 | 21/80 | 8/80 | 4/20 |
| Broader 29-position mask | 20/80 | 43/80 | 12/80 | 7/20 |
| K59D control | 10/20 | 6/20 | 0/20 | 0/5 |
| V111C control | 10/20 | 16/20 | 7/20 | 3/5 |

Joint pass denotes Tier A/B catalytic geometry together with the native-Thr233
conformation gate. The broader arm retained native-like Thr conformations more
often but did not improve Tier A/B-only occupancy uniformly. Retaining only
negative signed Herschlag phi yielded 11/80, 3/80 and 0/20 joint-pass sites
for broader, focused and WT, respectively. Positive-phi models remain separate
local binding-pose hypotheses. Broader mask_07 provided six joint-pass sites
across three of five seeds, all negative phi; focused mask_03 provided three
across two seeds. The four active sites within a model and the five seeds
within a sequence are nested observations, not independent activity assays.
The small pilot establishes neither statistical significance nor improved
chemical barriers or polymer turnover. K59D's zero joint passes despite
activity in noisy lysate/MS assays underscores the limited sensitivity of
these structural gates as activity classifiers. No activity-fitted cutoffs
were introduced. Source tables, acceptance records and scripts are under
`design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/analysis/`.

In canonical full-precursor numbering, mask_03 carries S50A, V51A, S153A,
A187T, V188I and L271F (beta local L39F); mask_07 carries S50A, V51A,
A92K, N97D, S153T, A187T and V188T. The shared exact substitutions are
S50A, V51A and A187T. All declared hard-protected sites remain unchanged;
candidate, generation and relaxed-chain sequence checks are documented in
`design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/sequence_comparison.json`.
For the flipped-only follow-up WT seed 04 has four negative-phi Tier A sites,
but none passes the native-Thr gate. This does not invalidate WT as a control;
it limits direct attribution of differences to sequence. Positive-phi results
above remain separate observations, not evidence of a forbidden binding mode.
Across five WT relaxed endpoints, Nα–leaving-N distances averaged 4.65 Å
(4.24–5.47 Å; 16 negative-phi sites) versus 3.61 Å (3.51–3.67 Å; four
positive-phi sites). Sites are nested within models. These data do not support
an automatic proximity advantage for the flipped pose; direct transfer would
require rearrangement, with water mediation another hypothesis. The shorter
positive-phi distances alone do not establish a relay mechanism either.

## Exploratory E–S sequence-design pilot (2026-11-22)

WT and eight designed sequences were modeled with four bound PA6 L4 molecules
and five prediction seeds each (45 tetramers). All minimization processes
completed; 44 achieved the common force tolerance. Raw predictions yielded
no Tier A/B sites. After relaxation, 76/180 assigned sites met Tier A/B, and
40/180 also passed the native Thr233 conformation gate. Excluding the single
unconverged model gave 38/176 joint-passing sites across 20 model seeds. WT
had three joint sites across two of five models; designed variants 03/04 each
had seven across four models, and variant 06 had eight across three models.
All four sites are nested within a model; the small ensemble and mixed original
prediction-host builds preclude a validated activity ranking. These E–S
geometry results motivate a controlled broader-mask comparison and matched
known-active controls, not a claim of transition-state stabilization or
improved catalysis. All models, failed gates and the unconverged case were retained.

# Results — working draft

## Leaving-group relay candidates — exploratory screen, 2026-11-21

Screening sixteen 100 ns WT trajectories and neighboring saved frames identified
one unique assisted configuration with all directional hydrogen-bond legs and
Tier B catalytic geometry: matched NH2/ASP272/ASP274 replica 03, frame 9386,
with Lys156–water29282–water4517–leaving-amide-N connectivity. Its d1/d2/d3
were 3.287/3.815/2.865 Å, attack angle 86.69° and signed phi +74.06°. However,
Thr OG1 and Nalpha displacements from the reference were 2.396 and 2.958 Å,
outside the current native-conformation criteria. This is a structural donor
network, not evidence that Lys acts as the catalytic acid.

An ASH272 replica 1, frame 7316 nearpass contained an OG1–water9397–leaving-N
network with Tier A geometry and terminal H–N distance 2.50964 Å, only 0.00964 Å
above the screening cutoff. It also failed native Thr displacement/psi criteria.
No frame among the twenty bounded native-audit selections passed all native
criteria. Three unchanged diagnostic restarts were retained; no accepted new
QM/MM seed or mechanism was assigned. These results do not exclude other
frames, water reorganization, or productive conformations outside the sampled
ensembles. Known direct-route controls were not counted as new assisted relays.

## Preliminary NEB-seeded umbrella sampling (2026-11-21)

**Closeout, 2026-11-22.** All 63 windows completed six accepted 5 ps
production blocks (30 ps/window; 1.89 ns aggregate). Using the same native
two-coordinate bias model, histogram grid and fixed E–S density-reference bin,
the updated marginal PMF maximum was 47.05 kcal/mol (previously 46.94).
All 50 within-window 5 ps block-bootstrap refits retained connected overlap;
the pointwise interval at the maximum bin was 45.41–47.87 kcal/mol. The
unbiased-target effective sample size remained only 7.29, and windows 028,
030 and 032 retained block-drift flags. This interval does not account for
unresolved slow-state sampling or Hamiltonian error. We therefore do not
interpret the profile maximum as a converged activation free energy. Further
sampling of this direct-route hypothesis was discontinued because of its
persistently unfavorable preliminary profile and low expected value relative
to testing alternative pathways. All trajectories, restarts and analyses were
retained; this decision does not establish that alternative mechanisms are
inaccessible or that the enzyme cannot catalyze PA6 hydrolysis.

An exploratory 300 K DFTB3/MM umbrella campaign seeded from the converged
regular direct-acylation NEB contained 1,190 ps of accepted production across
all 63 windows (10–30 ps/window) at this checkpoint. Connected sampled
overlap permitted a preliminary MBAR marginal along q_acyl = d(C,N_leaving)
− d(C,OG1). Its maximum was 46.50 kcal/mol relative to a fixed E–S
free-energy-density bin; the pointwise 95% within-window block-bootstrap
interval at that fixed maximum bin was 44.54–47.56 kcal/mol. This is not a
validated activation free energy: pooled target ESS was 28.28, three windows
showed block drift, and independent-start/hysteresis and slow-branch sampling
were not established. The 0 K source NEB maximum was 35.80895 kcal/mol;
these different reference definitions do not isolate an entropic correction.

HG1 was nearest Thr233 N-terminal N in 2,491 of 23,800 raw biased frames,
principally windows 028–034, rather than nearest the intended leaving N.
Those counts are not equilibrium populations. The alternative branch closely
tracked both distance-difference umbrellas: for example, window 031 had
mean q_H = 0.4417 Å versus target 0.4568 Å while mean HG1–leaving-N was
1.9161 Å. Thus the present CVs do not uniquely resolve the proton-transfer
route; merely increasing their force constants would not remove configurations
with the same CV values. All frames were retained. Sampling was ongoing at that checkpoint and
was closed after 30 ps/window, as reported above. The high projected maximum does not presently support a viable
direct-acylation mechanism. Frozen figures, tables, limitations and hashes are
in [the preliminary record](../docs/neb_us_preliminary_findings_20261121.md).


## Exploration of acylation pathways

We considered a direct concerted acylation route in which Thr233 OG1 attacks
the scissile substrate carbonyl carbon while HG1 transfers to the substrate
amide nitrogen, with cleavage of the scissile C–N bond leading to the
acylenzyme. This route involves a proposed four-membered arrangement of
OG1, HG1, the substrate nitrogen, and the acyl carbon. Coordinate driving
and an E–S-to-acylenzyme NEB calculation were used to explore this hypothesis.
Here, concerted denotes a coupled reaction hypothesis and does not require
synchronous changes in the individual bond distances. These exploratory
calculations have not yet established a validated first-order saddle point.

The original calculations treated Asp272 and Asp274 as deprotonated. The
triad-mediated mechanism and protonation hypothesis proposed for NylC [R25]
motivated a subsequent comparison with Asp272 protonated, Asp274 deprotonated,
and the Thr233 N-terminus neutral. The initial trajectory comparison is
reported below; reaction-path results for this state remain pending.

**Exploratory result, 2026-11-13.** Two 100 ns ASH272/neutral-N replicas
yielded 3,886 and 1,724 Tier A frames among 10,000 frames each. The prior
ASP272/neutral-N replicas yielded 37 and 238 under the existing geometric
screen. The larger ASH272 occupancies suggest greater near-reactive pose
persistence in these starts, but different initial ensembles and only two
replicas per state preclude a causal protonation or rate conclusion. A
one-water Thr233 N–Asp274 OD2 bridge occurred in 11.45% of sampled frames of
ASH replica 1, and an OG1–ASH272 OD1 bridge in 4.0% of ASH replica 2;
occupancy does not prove proton transfer. Four independently selected ASH272
E–S snapshots converged under unrestrained DFTB3/MM optimization and passed
the geometric Tier A filter. The replica-1 Tier A start gave the shortest
OG1–substrate-carbonyl C distance (2.778 Å) with oxyanion-hole distances
2.961 and 2.981 Å, αBD 103.84°, signed Herschlag φ +85.72°, and
CB–OG1–C 113.45°. However, in every optimized endpoint the original
Asp272 HD2 moved to Asp274 OD2 (H–O 1.02–1.05 Å), so none preserves the
nominal ASH272 microstate. This is an observed endpoint proton transfer, not
evidence by itself that Asp274 is protonated during catalysis. The native
scissile amide C–N and Thr233 OG1–HG1 bonds remained intact.

An earlier 16-image direct-route DFTB3/MM NEB stopped at its 500-cycle limit
without force convergence and had a ~70.6 kcal/mol image outlier. A later
25-image trace-seeded, smaller-step calculation converged by the native
optimizer at cycle 477. Among its 24 regular images, the maximum was
35.80895 kcal/mol relative to image 1 (image 8). Reactive-coordinate and
endpoint-identity checks passed, whereas a complete nonbonded-clash screen
and first-order saddle verification remain outstanding. The ~46.36 kcal/mol
direct dimer candidate was referenced to a separately minimized E–S
structure and is not the same observable or necessarily the same geometry;
these numbers cannot be interpreted as a measured barrier reduction. A
bounded ±0.05 Å finite-difference probe found energy and force directional
curvatures within 0.73%, but the global low-mode eigensolver failed to
converge eigenpairs within 800 Hessian-vector products. No transition state
is yet verified.

## Protonation-state interpretation and testable hypotheses

The ASH272 MD replicas had many more geometrically Tier A frames than the
earlier ASP272 replicas, but those trajectories imposed the ASH272 protonation
through their fixed-charge topology. The comparison also began from different
initial ensembles. Thus the occupancy difference supports a **testable
pose-stabilization effect of the chosen microstate**, not an experimental
assignment of the catalytic proton. The NylC crystal-panel Asp–Asp oxygen
distances below likewise do not establish Asp306–COOH. In all four subsequent
unrestrained Nyl12 QM/MM E–S optimizations, the proton initially on Asp272
localized on Asp274. This suggests that, for those four structures and that
Hamiltonian, Asp274–COOH is favored over the prepared Asp272–COOH local
minimum. It does not determine equilibrium populations in solution or during
turnover, where solvent, pH, temperature, and protein motion matter.

Three bounded hypotheses follow. First, a proton on either member of the
Asp272/Asp274 pair may improve binding geometry relative to two carboxylates,
independent of which Asp carries it. This can be tested with matched 300 K MD
replicas initialized from the same E–S pose for ASP/ASP, ASH272/ASP274, and
ASP272/ASH274, with identical sampling and Tier filters. Second, the proton
may exchange within the pair as the active site reorganizes; QM/MM
proton-transfer scans or enhanced sampling on matched snapshots can test
that without assuming the original ASH272 assignment. Third, the apparent
QM/MM transfer could reflect DFTB3, embedding, or a specific starting pose;
sparse higher-level QM/MM single-point/gradient comparisons on the *same
geometries* and independent E–S seeds can test its robustness. A classical
ASP272/ASH274 comparison needs a verified topology, proton placement, total
charge, and counterion inventory.

## Bounded water-assisted activation scan (2026-11-17)

Starting from the force-converged Tier A ASP272/ASH274 E–S model with a
neutral Thr233 N terminus, we tested a water-assisted activation branch using
a 102-real-atom DFTB3/MM subsystem containing WAT1576. The 25 restrained
optimization windows coupled a composite two-proton relay coordinate to an
independently restrained OG1–substrate-carbonyl distance. Each window reached
the engine force-convergence criteria. At the endpoint, Thr HG1 was bound to
the relay water (O–HG1 1.015 Å), and the water's original H2 was bound to
Asp272 (OD1–H2 0.985 Å). Asp274 retained its proton; the Thr N terminus
remained neutral. Thus this branch protonated Asp272 rather than the Thr
N terminus and left both Asp side chains neutral.

The endpoint did not satisfy the predefined tetrahedral-geometry gate:
OG1–carbonyl C was 1.852 Å versus a target of 1.650 Å, with carbonyl C–O
1.300 Å and scissile C–N 1.416 Å. The physical QM/MM potential energy rose
monotonically to +54.02 kcal/mol relative to the accepted E–S seed; including
the remaining 4.13 kcal/mol restraint energy gives +58.15 kcal/mol. These
values describe a restrained PES scan, not an activation free energy or a
converged saddle-point barrier. The conditional unbiased endpoint release
was not launched because the geometric gate failed. The result disfavors
this sampled branch but does not establish that TI1 is absent under another
protonation pattern, reaction-coordinate timing, protein conformation, or
Hamiltonian. In particular, it does not directly test the proposed
Thr-NH3+/oxyanion tetrahedral state.

Across the calculations examined so far we have not established a stable TI1
or a verified first-order acylation saddle. Earlier TI1-like releases reverted
toward E–S under their respective setups; those negative results and this
incomplete scan are not equivalent tests of the same chemical state. A
tetrahedral geometry may also lie on a concerted addition–elimination pathway
without forming a persistent minimum. The next discriminating test is a
bounded two-coordinate activation/attack map on one accepted E–S model,
followed by unbiased release and sparse matched-embedding method checks.

Evidence: `qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/`
(`window_025/acceptance.json`, `protocol.json`, `energy_ledger.csv`). No new
PMF, validated TS, or intermediate-stability claim follows from this scan.

## Independent activation/attack grid: completed diagnostic (2026-11-18)

All 25 reported restrained grid points completed. No point passed the
declared tetrahedral-intermediate geometry gate. The terminal attack points
in the four rows retaining the intended proton inventory had OG1–C distances
of 1.887, 1.886, 1.860 and 1.783 Å, with physical potential energies
+28.34, +32.71, +36.62 and +54.58 kcal/mol relative to the common E–S seed.
Attack remained uphill in each row and lagged its nominal 1.50 Å target.
This incomplete coverage of actual reaction coordinates prevents a conclusion
that TI1 is absent from the full potential-energy surface.

The highest nominal q_N row followed an atom-label exchange branch. At its
first grid point, Thr HG1 was bound to water O (0.986 Å), while water H1
was bound to Thr OG1 (0.980 Å); the Thr N terminus remained neutral.
Thus the −4.54 kcal/mol energy of this point describes a reorganized
reactant geometry, not a proton-activated intermediate. The signed distance
difference of a single labelled hydrogen did not distinguish exchange with
solvent from chemical activation. Subsequent reaction-coordinate design must
monitor proton occupancy of both donor and acceptor, allowing exchange of
chemically equivalent protons. Sparse selected-state method comparisons are
still in progress; no new barrier or unbiased intermediate is established.

Evidence: `08_activation_attack_grid/grid_points.json` and each point's
`acceptance.json`, under the model30 ASH274 comparison tree. Hydrogen
destinations were checked against the nearest full-system heavy atoms for
the first and last points in row 5, rather than inferred from q_N alone.

## Unbiased activation release and proton-occupancy audit (2026-11-18)

We audited all 25 grid points using the 54 real QM hydrogen atoms, excluding
fictitious link atoms, and retained each hydrogen's nearest heavy-atom
partners and declared ambiguity thresholds. The audit confirmed that the
highest labelled-hydrogen coordinate row exchanged Thr HG1 with water H1;
it did not produce the intended persistent proton-activated species.

A separate unbiased DFTB3/MM minimization from q_004_attack_001, whose HG1
was initially nearer Thr N, reached force/SCC convergence. The same proton
exchange occurred: HG1 bound water O at 0.985 Å and water H1 bound OG1 at
0.980 Å. A full real-QM hydrogen audit found two protons on Thr N and one
on OG1, consistent with E–S chemistry. OG1–substrate C was 3.205 Å. This
endpoint was 4.7757 kcal/mol below the original E–S reference and retained
Tier B geometry (d2/d3 3.103/2.968 Å; attack angle 117.39°; signed phi
+102.52°). These results identify a relaxed alternative reactant endpoint,
not a stable proton-activated intermediate or TI1. Its lower potential energy
does not establish global conformational or sampling convergence.

The smooth, label-independent N–H minus OG1–H coordination coordinate
reproduced the independently audited values at every point. A native
Sander/PMD integration test nevertheless applied none of its intended
0.20 kcal/mol occupancy bias during the tested DL-FIND minimization. Native
CV input files alone are therefore insufficient for the next optimization
campaign. The analytic coordinate implementation was tested independently;
no occupancy-driven reaction path was launched using the inactive hook.

Evidence: `08_activation_attack_grid/proton_occupancy_audit/`,
`diagnostic_activation_release/scientific_assessment.json`, and
`occupancy_coordinate_preflight/attempt_2/scientific_acceptance_corrected.json`.

## Saved-gradient diagnosis and verified occupancy callback (2026-11-18)

Saved DFT and GFN2-xTB gradients showed local resistance to nucleophilic
attack in all three tested geometries. The table reports derivatives with
respect to increasing OG1–C distance while holding q_occ fixed; a negative
value means shortening the distance is locally uphill.

| Geometry | ωB97X-D/def2-SVP | GFN2-xTB |
| --- | ---: | ---: |
| E–S | −10.49 | −9.51 |
| Partially activated (q_004_attack_001) | −17.78 | −19.41 |
| Attack (q_003_attack_005) | −59.59 | −55.26 |

Units are kcal/mol/Å for the matched fixed-point-charge **embedded-QM**
component. The corresponding q_occ slopes (kcal/mol per coordination unit)
were +1.95/−0.51, −0.98/+4.12 and +9.98/+5.64 for DFT/xTB, respectively.
Thus xTB reproduced attack resistance but disagreed in sign on two small
proton-transfer derivatives. These three local comparisons do not validate
xTB globally or establish the absence of TI1. Native DFTB3 total QM/MM
derivatives are retained separately; their raw magnitudes are not directly
comparable to this embedded-QM-only table. No new QM calculations were needed.

The replacement analytic callback passed a real optimizer integration proof:
the prescribed initial 0.20 kcal/mol bias and its nonzero analytic force were
present, and input geometry, energy units and force composition passed their
explicit gates. Four native calls supported two deliberately bounded optimizer
iterations; this was a coupling test, not evidence of molecular convergence.
A four-point occupancy/attack pilot was then launched, with process, force
and chemical-state acceptance recorded separately. It has not yet established
a new intermediate, reaction barrier or PMF.

Evidence: `08_activation_attack_grid/reference_checks/gradient_analysis/`
and `occupancy_callback_pilot/callback_proof/scientific_acceptance.json`.

The first pilot point stopped at 158 accepted optimizer iterations before
force convergence (maximum component 0.277869 versus 0.12 kcal/mol/Å).
The final trial displacement was only 4.77×10⁻⁹ Å, yet the total energy changed
by −0.0010040 kcal/mol rather than the +3.10×10⁻⁸ kcal/mol predicted from
the analytic gradients. One ordinary MM pair, GLU754 OE2–GLN757 O (topology
atoms 11272/11330), crossed the 10 Å cutoff. Its calculated screened-Coulomb
plus Lennard–Jones removal energy was −0.0009792 kcal/mol; the difference
from the observed jump was below native energy print resolution. Predicted
and observed equal/opposite gradient jumps likewise agreed within force print
resolution. Real-QM cutoff membership was unchanged. This numerical
discontinuity explains failure of the smooth-energy line search; it does not
indicate a chemical barrier or electronic convergence failure. A gradient-based
DL-FIND continuation was launched without changing the potential or force gate.
Evidence is retained in `occupancy_callback_pilot/line_search_diagnosis/`.
With the gradient-based continuation, A met the unchanged force criterion
after three additional steps (maximum component 0.109355 kcal/mol/Å), and
the next attack point started. The force-converged biased A geometry retained
three N-terminal protons, no covalent proton on OG1, intact water and the
ASH274 proton. This verifies the numerical repair, not an unbiased activated
minimum or a stable tetrahedral intermediate.

## Completed occupancy-controlled attack pilot (2026-11-18)

All four points met the declared maximum biased-gradient criterion. The
intended net proton inventory persisted: Thr N carried three protons, OG1
none, water two and ASH274 one. The only real-QM hydrogen ownership change
relative to original E–S was transfer of HG1 from OG1 to Thr N. The smooth
coordinate therefore avoided the previous proton-label exchange artifact.

| Point | OG1–C (Å) | Substrate C–O (Å) | Substrate C–N (Å) | Physical ΔE vs E–S (kcal/mol) |
| --- | ---: | ---: | ---: | ---: |
| A | 3.1403 | 1.2696 | 1.3508 | 17.052 |
| B | 2.6623 | 1.2689 | 1.3551 | 20.633 |
| C | 2.2266 | 1.2699 | 1.3734 | 28.416 |
| D | 1.8537 | 1.2919 | 1.4090 | 42.763 |

At D, the carbonyl's out-of-original-neighbor-plane height increased to
0.2783 Å from 0.0098 Å in E–S, with angle sum 348.57° and Bürgi–Dunitz
angle 108.24°. These descriptors show attack-like pyramidalization, but D
did not meet the established attack ≤1.75 Å and C–O ≥1.30 Å TI1 screens.
Its bias energy was 6.630 kcal/mol. These restrained diagnostic points do
not establish a potential-energy barrier, free-energy barrier or stable TI1.

The substrate oxygen lost the Phe113 backbone N contact (3.161→4.110 Å),
while the Asn186 ND2 contact remained short (2.837 Å versus 2.915 Å in E–S).
Loss of one oxyanion-hole contact is a testable contributor to unfavorable
attack geometry/energy; it does not by itself identify the cause or quantify
its energetic contribution. One unbiased D release and a matched sparse
electronic-gradient comparison were initiated to discriminate persistence
and local method sensitivity. Their completed outcomes are given below.

Evidence: `occupancy_callback_pilot/endpoint_audit/` and each final point's
`scientific_acceptance.json` and accepted-iteration ledger.

## Completed D release and limits of a TI1-centered interpretation

Removing both driving biases from D yielded a force-converged E–S geometry
after 249 accepted steps/250 native evaluations. OG1–C returned to 3.211 Å,
C–O to 1.264 Å and C–N to 1.356 Å. Thr N retained only H1/H2, original HG1
bound water O, and water H1 bound OG1. Water remained H2O and ASH274 retained
its proton. This is an explicit negative persistence result for the selected
attack geometry; it does not establish that all TI1 basins or productive
hydrolysis paths are absent. D itself had not met the covalent/tetrahedral
geometry screens before release, so its return is also compatible with release
from the reactant side of a later reaction barrier.

The completed finite-embedding electronic checks at the exact D geometry
showed that both ωB97X-D/def2-SVP and GFN2-xTB oppose additional attack at
fixed q_occ: dE/dr_OG1-C = −31.80 and −26.81 kcal/mol/Å, respectively.
Their q_occ derivatives were +8.55 and +7.81 kcal/mol per coordination unit.
These matched electronic-component trends do not support attributing the
current difficulty solely to DFTB3, but do not validate its total periodic
QM/MM surface. Shortening the Phe113 contact at fixed q_occ and attack also
has a positive local electronic derivative (+10.91/+17.05 kcal/mol/Å of
shortening); isolated restoration of this contact is not demonstrated to
rescue the geometry. These are directional derivatives, not hydrogen-bond
energies or barriers, and exclude some mechanical contributions.

The diagnostic pilot varied attack while favoring retained protonation of
Thr N. It did not independently explore cleavage of the substrate C–N bond
and proton delivery to its nitrogen. In the usual N-terminal base/acid cycle,
net transfer of a proton onward from NH3+ to the leaving group would reduce
Thr-N proton occupancy; a fixed high q_occ penalizes that branch unless a
different proton simultaneously replenishes it. Thus the pilot tests one
activation/attack branch, rather than the complete acylation mechanism.
The return to E–S does not distinguish a stepwise route with an unsampled TI1
from a coupled acylation route without a stable TI1. Literature precedent
[R48] supports treating these as alternatives, rather than imposing a minimum
as a prerequisite for productive chemistry.

Primary artifacts: `pilot/point_D/unbiased_release/release_assessment.json`,
`process.json`, and `d_reference/{results.json,gradient_projections.json}` under
the occupancy callback pilot. Full stationary-point classification by Hessian
has not been performed; no new catalytic barrier is claimed.

A separate coordinate audit clarifies the leaving-group problem. The distance
between Thr Nalpha and the substrate leaving N is 5.80 Å in original E–S,
5.23 Å at D and 5.88 Å after release. The closest Nalpha-bound proton remains
4.49 Å from the leaving N at D. Direct Nalpha proton donation is therefore
not arranged in these selected structures. At D the QM water is closer
(O–N 3.35 Å; nearest owned H–N 2.60 Å), whereas the ASH274 proton is 8.19 Å
away. These distances motivate testing a reorganized proton relay and nitrogen
orientation; they do not establish that such a relay is energetically viable.
The audit is in `analysis/nyl12_mechanism_evidence_20261118/acid_reach.json`.
It reinforces that attack geometry alone does not specify a complete chemical
pathway for amide cleavage.

## Existing-MD screen for leaving-group proton delivery (2026-11-18)

A bounded search of the six existing 100 ns, neutral-Thr-NH2 E–S replicas
at 10 ps/frame found no preorganized one-water Nalpha→water→leaving-amide-N
bridge under the unchanged structural gates. Of 60,000 frames, 18,958 passed
Tier A/B geometry with positive signed Herschlag φ, but only 51 also passed
the native-Thr and intact-PA6 gates. The primary ASP272/ASH274 replicas
contributed eight and zero such frames; the ASP272/ASP274 controls contributed
zero and 31, and the ASH272/ASP274 controls 12 and zero. All 51 lacked even
the required ≤3.5 Å water-O contact to both nitrogens; therefore none reached
the leaving-leg or strict two-leg hydrogen-bond classes. The six-case table
and exact unchanged thresholds are in the executed Methods section below.

Native displacement, rather than Tier geometry alone, was the main filter.
In primary replica 1, OG1 displacement exceeded 0.75 Å in 5,851/5,865
Tier A/B positive-φ frames and Nalpha displacement exceeded 1.25 Å in
5,747/5,865 (overlapping failures); their median displacements were 2.03
and 2.35 Å. In the eight fully eligible primary frames, Nalpha–Nleave
was 4.06–4.56 Å. The water minimizing the larger of its two N–O distances
was WAT22725 in frame 93: Nalpha–O 7.09 Å, O–Nleave 4.99 Å, and its best
owned H–Nleave contact 5.14 Å at 76.1°. No starting restart was selected.

This negative result concerns preorganization in the sampled fixed-charge
E–S ensembles for mature-substrate hydrolysis. It does not exclude transient
water reorganization, multiple-water networks or other conformations, and
does not test enzyme maturation/autocleavage. Neither ensemble-averaged
acidity nor a catalytic free-energy barrier is inferred. Authoritative
counts and primary nearpasses are in `00_screen/bounded_screen_acceptance.json`
and `00_screen/diagnostics_attempt1/diagnostics.json` under
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/`.

## NEB numerical-spacing diagnostics — provisional (2026-11-19)

The ASP272/ASH274 acyl-path NEB showed repeated native kink repairs after
optimizer step capping. Reset-associated physical maxima increased
35.08→82.74 kcal/mol (cycles 244→245) and 32.83→75.15 (269→270), interrupting
local relaxation. Thus a low instantaneous maximum does not establish sustained
progress. Four matched-source numerical trials now compare uniform 24/48-image
and adaptive 24-image L-BFGS bands, with a short 24-image CG diagnostic. All four
movies have published; native and comparison monitors remain active. Their
early profiles are initialization diagnostics: neither improved convergence
nor a validated barrier, transition state or kinetic mechanism is established.
Force reduction, resets, coherent chemical continuity and actual work remain
the comparison criteria. The retained source is cycle 268; exact protocols
and seed/startup exceptions are documented in Methods and the spacing-test
README. These numerical seed reviews do not certify chemical stationarity.

## Barrier scale and kinetic interpretation

The assay used 75 °C (348.15 K) for 6 h, whereas the present MD used about
300 K. With transition-state theory, `k = κ(k_B T/h)exp(−ΔG‡/RT)`; setting
the transmission factor κ to one gives illustrative single-step rates at
348.15 K of about 0.11 s⁻¹ for ΔG‡ = 22 kcal/mol, 0.0062 s⁻¹ for
24 kcal/mol, 0.00035 s⁻¹ for 26 kcal/mol, and 1.1×10⁻⁶ s⁻¹ for
30 kcal/mol. Roughly 22–26 kcal/mol is therefore an initial plausibility
range for a bound chemical step on the assay timescale, **not an acceptance
cutoff**. Crude-lysate conversion of polymer powder also depends on active
enzyme concentration, substrate accessibility and binding, competing steps,
and product release. A NEB potential-energy maximum is not ΔG‡ and cannot
be inserted into this rate expression. A converged PMF at assay temperature,
with an appropriate reactant-state reference, would offer a more relevant
comparison. This calculation has not yet been done.

## Asp–Asp geometry across the NylC crystal panel

Across the nine Table 1 entries of Negoro et al. (2023), we measured 31
deposited sites, including mature enzymes, uncleaved precursors and a
chemically cleaved Cys267 variant. The smallest Asp306/Asp308 oxygen-pair
distance at each site ranged from 3.5794 Å (5XYO chain B) to 4.1736 Å
(7YU2 chain A). None was below 3.5 Å. Thus the deposited panel does not
provide a short Asp–Asp oxygen contact that independently establishes a
proton-shared carboxyl pair. These distances do not determine protonation,
exclude substrate-induced rearrangement, or invalidate an acid–base mechanism.

| PDB | Enzyme/variant | State | Deposited sites | Site-minimum O–O range (Å) |
| --- | --- | --- | ---: | ---: |
| 5XYG | NylCp2 wild type | Mature | 2 | 4.0985–4.1535 |
| 3AXG | NylCA wild type | Mature | 15 | 3.8092–4.0669 |
| 5XYO | NylCp2 D122G | Mature | 2 | 3.5794–3.8284 |
| 5XYT | NylCp2 H130Y | Mature | 2 | 3.8701–4.0311 |
| 5Y0L | NylCp2 D122G/H130Y | Mature | 2 | 3.8618–3.9588 |
| 5Y0M | NylCp2 D122G/H130Y/D36A/E263Q | Mature | 2 | 3.9110–3.9163 |
| 7YU0 | NylCp2 H130Y/N266A/T267A | Precursor | 2 | 4.1220–4.1417 |
| 7YU1 | NylCp2 D122G/H130Y/T267C | Precursor | 2 | 3.9752–4.0344 |
| 7YU2 | NylCp2 D122G/H130Y/T267C | Hydroxylamine-cleaved | 2 | 4.0995–4.1736 |

The 3AXG count refers to 15 deposited chains with complete target oxygen
pairs; additional nominal or symmetry-related sites were not imputed.
The mature Fig. 10 reference, 5Y0M, gave minima of 3.9163 and 3.9110 Å
in its two independent sites. Consequently, the mechanistic proposal that
Asp306 becomes protonated should be treated as circumstantial support for
testing explicit microstates, rather than a protonation assignment measured
from this crystal panel. NylC Tyr146 maps to Nyl12 Phe113, while Lys189,
Asn219, Thr267, Asp306 and Asp308 map to Lys156, Asn186, Thr233, Asp272
and Asp274, respectively. The Tyr phenolic OH mechanism therefore cannot be
transferred unchanged to the actual Nyl12 sequence.

## Static pKa sensitivity: cleavage state and Asp-pair hypotheses

The corrected NylC predictions are summarized below. Ranges span resolved
chains, not statistical confidence intervals. Values are model predictions;
none constitutes an experimentally determined protonation assignment.

| NylC input | PROPKA output | Asp306 pKa | Asp308 pKa | Thr267 free N-terminal pKa |
|---|---|---:|---:|---:|
| Mature 5XYG, corrected termini | Default | 2.19–2.21 | 5.49–6.49 | 6.35–6.55 |
| Mature 5XYG, corrected termini | Coupled alternative | 2.19–2.21 | 7.24–8.45 | 4.59–4.60 |
| Uncleaved 7YU0 | Default | 4.92–6.34 | 4.18–6.37 | Not a free terminus |
| Uncleaved 7YU0 | Coupled alternative | 4.99–6.22 | 5.07–5.53 | Not a free terminus |
| Uncleaved 7YU1 | Default | 6.36–6.39 | 4.24–6.53 | Not a free terminus |
| Uncleaved 7YU1 | Coupled alternative | 5.01–6.39 | 5.59–6.53 | Not a free terminus |

The precursor structures were not back-mutated. Differences therefore combine
cleavage state, mutations and structural environment; they cannot be attributed
to maturation alone. The corrected mature predictions do not favor specifically
protonated Asp306. The coupled alternative instead permits an elevated Asp308
pKa. Its correspondence to Nyl12 Asp274 motivates a competing microstate test,
not a transferred protonation assignment. Likewise, a relatively low mature
NylC terminal pKa is consistent with testing a neutral N-terminal base but
does not establish the corresponding Nyl12 state.

For the substrate-facing Nyl12 beta site E, the original default predictions
were:

| MD starting geometry | Asp272 pKa | Asp274 pKa | Thr233 free N-terminal pKa |
|---|---:|---:|---:|
| ASP272/ASP274 replica 2, frame 284 | 1.02 | 3.44 | 9.70 |
| ASH272/ASP274 replica 1, frame 4868 | 0.88 | 6.95 | 8.93 |
| ASP272/ASH274 replica 1, frame 3492 | 2.37 | 5.84 | 9.62 |

Across all four sites in these three snapshots, Asp272 spans 0.88–6.67,
Asp274 0.83–7.98 and the Thr233 amino terminus 6.53–9.87. Thus neither
specifically protonated Asp272 nor a universally neutral Nyl12 Thr233 amino
terminus is established. The site-E terminal estimates tend toward NH3+ near
pH 8 in an isolated-site interpretation, contrasting with mature NylC. Missing
substrate, conformational sampling, group coupling, crystal mutations and the
uncertain assay-temperature pH limit that comparison. No fixed protonation
assignment was changed automatically on the basis of these screens.

The MCCE results cannot resolve the uncertainty: mature NylC's apparent
Asp pKas below 0 and terminal pKas above 14 were scientifically rejected
because of artificial internal termini, while complete Nyl12 failed input
preflight before titration. Treating these as a valid independent prediction
disagreeing with PROPKA would be misleading. Current ASH274 MD comparisons
therefore test sensitivity to an explicit hypothesis; increased geometric
occupancy does not establish that microstate's equilibrium population.


## Diagnostic MCCE–PROPKA comparison and nearby-site flags (2026-11-20)

The corrected complete-Nyl12 MCCE pilot finished but did not establish reliable
quantitative catalytic pKas. Most target curves were saturated over pH 0–14.
A single-sigmoid Asp272 fit disagreed with its observed midpoint, and the
native His257-E >14 label contradicted an almost entirely neutral-state curve.
The pinned FTPL proline backbone charges also differed from the bundled legacy
template; this parameter inconsistency has not been shown to explain the
anomaly. These findings limit interpretation independently of native exit status.

Matched PROPKA predicted Thr233 N-terminal pKas 7.06–9.70, Asp272 1.02–1.99,
and Asp274 2.95–4.37 for the HIS-normalized snapshot. Coupled alternatives and
all four site values are retained separately. The protein copies have different
MD conformations, so exact pKa symmetry is not expected. Nearby Lys156 had
PROPKA values 8.46–9.91 versus its model 10.5; His257 4.85–5.22, Tyr267
14.00–14.57 and Asp45 −1.97 to −0.06 were additional shifted groups.
These are hypotheses about electrostatics/proton availability, not established
catalytic roles or a basis for automatic protonation changes. Substrate removal,
the single structural snapshot and method disagreement remain limitations.


## Matched-design MD: two independent replicas (2026-11-24)

All eight WT, focused mask_03, broader mask_07 and V111C production runs
completed 100 ns and passed native completion/final-temperature gates. Each
was analyzed at 1,000 times (100 ps spacing) across four fixed substrate/site
assignments, giving 4,000 site-frame observations per replica. The two velocity
replicas are independent units; sites and time-correlated frames are nested.
Native-Thr similarity is reported separately and does not censor these values.

| Sequence | Replica 1 Tier A / A-or-B (%) | Replica 2 Tier A / A-or-B (%) |
|---|---:|---:|
| WT | 10.3 / 45.7 | 12.6 / 55.8 |
| Broader mask_07 | 30.5 / 69.5 | 17.4 / 46.4 |
| Focused mask_03 | 4.8 / 23.3 | 10.9 / 38.7 |
| V111C | 5.7 / 33.5 | 4.7 / 35.6 |

Mask_07 had more Tier A frames than WT in both replicas, but its broader
A-or-B advantage did not reproduce in replica 2. Conditional Tier A among
A-or-B observations was 43.9/37.5% for mask_07 and 22.5/22.6% for WT.
Focused mask_03 and the known-active V111C control retained these geometries
less often than WT. These results support mask_07 as a structural follow-up,
not a prediction of improved catalysis or justification for experimental
spending by themselves. Only two replicas per sequence were analyzed, starting
poses differ, and none of these nonreactive simulations measures a reaction
barrier. Source-bound tables, first/second-half summaries, native receipts and
plot are retained under `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/`
`md_followup/analysis/matched_replica_comparison.{json,csv,md,png}`.


## Native-Thr gate sensitivity and structural interpretation (2026-11-24)

The rare joint native-Thr/Tier A/B passes are dominated by simultaneous
OG1/Nalpha displacement and psi deviations, rather than chi1 rotamer changes.
Among Tier A/B frames, WT versus mask_07 failure rates are 99.89/98.67% for
OG1, 99.84/99.03% for Nalpha, 83.82/90.21% for psi and 5.47/3.71% for chi1.
Modest diagnostic relaxation (1.0 Å OG1, 1.5 Å N, 45-degree psi, 30-degree chi1)
still admits only 2/1829 and 14/2779 frames. The original static WT starts are
borderline OG1 failures (0.757–0.880 Å versus 0.75 Å); the MD distributions
show larger real shifts. Independent alpha- and beta-scaffold fits do not
remove the discrepancy in a bounded WT sample. The original torsion
definitions/circular differences were verified; N-terminal phi is not used.

No numerical thresholds changed. Following the diagnostic audit and PI
approval, native resemblance is retained as a continuous structural diagnostic
and does not reject or rank candidates; it is not a validated activity classifier.
The crystal and ligand-bound MD need not have identical ensembles, and known-active
control rejection already limits hard-filter interpretation. Matched-site
representative and rare closest-native overlays, initial controls, scalar
sensitivity and portable PyMOL scenes are in
`design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/analysis/native_thr_md_gate_audit_20261124/`.


## Initial retained-product acylenzyme hydration screen (2026-11-24)

The initial 6 ps reconnaissance (1 ps heating and 5 ps 300 K NVT) was screened
at every saved frame, including its exact starting endpoint: 121 configurations
and 28,245 water molecules. The closest observed water-O–acyl-C distance was
5.662 Å for WAT23286; its Ow–C–Ocarbonyl angle at that closest approach was
37.24°. No saved configuration met the exploratory combined near-attack
annotation (Ow–C≤3.5 Å and angle 95–125°). The final nearest distance was 5.786 Å,
compared with 6.246 Å initially. The acylenzyme remained covalent at the accepted
endpoint (C–OG1=1.310 Å; C–leaving-N=3.250 Å).

This is a limited hydration observation with the amine product retained and a
fixed exterior, not evidence that catalytic water cannot enter or that
hydrolysis is unfavorable. Water remains MM and therefore does not carry out
QM bond changes in this pilot. Heating and 300 K trajectory fractions are
reported separately. The PI authorized an additional 20 ps in four checkpointed
5 ps blocks to examine further approach before selecting any expanded-QM state.
No deacylation barrier or sampling-convergence claim is made.

The first accepted 5 ps continuation brought the analyzed trajectory to 11 ps
(221 saved configurations including heating and the exact initial endpoint).
WAT23286 reached a closest Ow–acyl-C distance of 5.422 Å at 9.7 ps, with an
Ow–C–Ocarbonyl angle of 47.59°. No water met the exploratory combined
near-attack annotation. The accepted 11 ps endpoint retained the acyl bond
(C–OG1=1.373 Å). These are interim hydration observations; sampling continues.


## N-terminal-relay acylation: retained endpoints and TS candidate (2026-11-25)

The reverse-activation and forward-acylation branches reached E–S-like and
acylenzyme-like states, respectively. Removing their reaction-coordinate
biases retained the respective connectivities in the same 12 Å movable model
with its original frozen exterior. Independent maximum gradients were
4.79e−6 and 9.88e−6 Eh/bohr. E–S C–OG1/C–Nleave distances were 2.963/1.380 Å,
whereas acylenzyme distances were 1.327/3.019 Å. The proton HG1 resided on OG1
in E–S and on the leaving nitrogen in acylenzyme, while Nalpha retained two
hydrogens at both endpoints. The shared window 033 source remained minimum-
unverified; its role here was a common exploratory relay-like seed.

Relative to the regenerated E–S endpoint, the combined driven physical PES
maximum was 16.31 kcal/mol, and the acylenzyme endpoint was 4.15 kcal/mol lower.
The forward dimer search converged to a candidate 16.312 kcal/mol above E–S
using independent same-Hamiltonian energies (−718.1008556 versus
−718.1268512 Eh). Its independent maximum gradient was 3.40e−5 Eh/bohr.
First-order curvature and downhill connectivity remain unverified. These PES
values are not activation free energies. The subsequent sampling status is
reported below.
The conformational reference/fast-QM accuracy and sampling qualifications in
the Methods still apply. The reverse candidate search was stopped after repeated SCC failures/resets,
large energy/curvature excursions and an E–S-like saved geometry. Its
SCC-contaminated curvature is not used for TS qualification.

A complementary extraction of the eight matched classical-MD CSVs found no
Tier A/B site-frame observation with Nalpha–Nleave ≤3.5 Å. Medians conditional
on Tier A/B were 5.99–6.34 Å. These descriptive distances reinforce that the
original attack/oxyanion-hole geometry filters do not establish proton-relay
preorganization; they are not new hard selection thresholds or evidence that
induced-fit/water-mediated pathways are impossible.


## Origin and limits of the retained TI1-like candidate (2026-11-25)

The retained window 033 candidate was selected from an existing biased QM/MM
MD trajectory (window 033, block 05, frame 46; recorded time 28.30 ps), then
released from the umbrella potential and optimized with the unchanged DFTB3
Hamiltonian, QM partition/proton inventory and 12 Å movable-region protocol.
C–OG1 lengthened from 1.556 to 1.740 Å, C–Nleave remained 1.490 Å and C–O
remained elongated at 1.309 Å; HG1 stayed on Nalpha at 1.050 Å. The independent
maximum gradient was 7.94e−6 Eh/bohr. Full-minimum curvature qualification was
unresolved, and the geometry-only classifier returned `unassigned`. We therefore
report a retained partially tetrahedral/addition-like candidate, rather than a
verified stable TI1 minimum or a finite-temperature intermediate lifetime.

A matched release from window 032 lost the addition contact (C–OG1=2.661 Å,
C–Nleave=1.371 Å, C–O=1.262 Å), although HG1 remained on Nalpha at 1.079 Å.
Its geometry-only `reactant_like` label describes loss of addition character,
not return to the neutral-Thr protonation microstate. Similar initial addition
distances and the same model/proton inventory led to different outcomes.
The results support conformational/environmental sensitivity; they do not
identify a unique stabilizing contact or demonstrate that an altered force
field, Asp protonation state or QM water caused retention. Reactive-state-
conditioned QM/MM sampling supplied a useful alternative starting arrangement,
but it did not itself establish an unbiased intermediate population.

### Bias-free release of coordinate-driving plateau images (2026-11-26)

Ten additional plateau-image releases (images 11–19 and 21) reached native
force convergence and passed independent maximum-gradient checks, all below
9.1×10⁻⁶ Eh/bohr. Image 20 reused the retained window 033 structure. Reaction-
coordinate biases were removed; the common 12 Å movable region and frozen
exterior were retained. The releases converged into two groups of addition-
like geometries rather than returning to the noncovalent E–S geometry.

| Starting images | C–OG1 (Å) | C–Nleave (Å) | Carbonyl C–O (Å) | Relative potential energy (kcal/mol) |
| --- | --- | --- | --- | --- |
| 11–15 | 1.808–1.811 | 1.476–1.477 | 1.300 | 0–0.012 |
| 16–19, 21 | 1.739–1.743 | 1.490 | 1.309 | 1.386–1.392 |

HG1 remained on Nalpha at approximately 1.050 Å in both groups. Energies
are relative to the lowest released endpoint (image 12) under the common
Hamiltonian and fixed exterior; they are neither free energies nor estimates
of basin populations. The elongated C–OG1 contacts support partially
tetrahedral/addition-like character but do not establish a conventional fully
formed covalent TI1. The geometry classifier remains `unassigned`.
No endpoint in this release series has been qualified as a local minimum.
Window033's full-Hessian verification remains inconclusive because of
unresolved eigenvalue signs; its subsequent selected-column repair failed
the numerical step-stability gate. Native acceptance and endpoint structures
are retained in `qmmm/pa6_flipped_model30/12_nterminal_relay/03_plateau_releases/`.

### Relay umbrella-sampling progress (2026-11-26)

All 57 windows completed 20 ps of production each (1.14 ns aggregate) by the
26 November observation. The first completed MBAR checkpoint used only the
first 5 ps/window: its highest sampled bin was 20.61 kcal/mol above the E–S
reference, with unsampled projection bins and insufficient effective sampling.
No global activation free energy is assigned. Later accepted native blocks
were temporarily excluded by a derived trajectory-time validation failure.
Analysis-only recovery regenerates their compact records with corrected time
serialization handling, preserving native inputs, trajectories and scientific
acceptance gates. The recovered cumulative 10–20 ps profiles filled the
earlier internal projection gap. At 15 and 20 ps, however, the sampled
reactant minimum lacked occupied neighboring bins, so the predeclared
reference-resolution gate withheld a global barrier. Local overlap and
effective-sampling gates also remained unmet. These estimates remain
preliminary; neither apparent shoulders nor a favorable maximum establish
a converged activation free energy.

## Why the N-terminal relay is a plausible lower-energy route (2026-11-25)

The relay offers a chemically plausible division of acid/base catalysis.
The neutral Thr233 alpha-amine can accept the hydroxyl proton, increasing
OG1 nucleophilicity, and the resulting alpha-ammonium group can subsequently
protonate the leaving nitrogen as addition disrupts amide conjugation and
C–N cleavage becomes favorable. This is consistent with the substrate-
hydrolysis mechanism calculated for human asparaginase III (R41), and with
an N-terminal proton donor in the cysteine Ntn-hydrolase study R54. These are
analogies, not evidence that Nyl12 shares every intermediate or rate-limiting
step. A neutral amide nitrogen has a resonance-delocalized lone pair and is a
poorer proton acceptor than an ordinary amine; proton delivery can become
more favorable after carbonyl addition changes that electronic structure.

Unlike direct OG1-to-leaving-N transfer during attack, a relay need not satisfy
a compressed OG1–H–Nleave–C four-center arrangement while preserving a useful
carbonyl-attack angle. The retained window 033 candidate provides direct
geometric evidence for delivery preorganization: Nalpha–HG1=1.050 Å,
HG1–Nleave=1.975 Å and Nalpha–HG1–Nleave=179.44 degrees. The released
acylenzyme instead has HG1–Nleave=1.026 Å and C–Nleave=3.019 Å. These distances
are consistent with the mapped relay connectivity; they do not determine
electronic bond orders, transmission coefficients or a TS ensemble. The
source-hashed measurements are in
`qmmm/pa6_flipped_model30/12_nterminal_relay/proton_relay_geometry_20261125.json`.

Formation of an OG1-anion/alpha-ammonium pair still has an energetic cost.
Local electrostatics and hydrogen bonding may offset that cost, while the
Phe113 backbone NH/Asn186 side-chain NH donors can stabilize developing
carbonyl-oxygen charge. Our present calculations have not isolated their
individual energetic contributions or established that a discrete ion-pair
minimum is necessary. Earlier direct paths and the present relay also have
different conformational histories and E–S references: their apparent PES
maxima cannot be subtracted to assign a relay-specific catalytic free-energy
advantage. The ongoing jointly normalized US/MBAR campaign is the next test;
for identical endpoints, a different route can lower the activation barrier
but cannot change the overall reaction free energy.
Hamiltonian accuracy, orthogonal sampling and conformational reference remain
independent qualifications. No additional chemical restraint or target barrier
was introduced to support this interpretation.
# PMF statistical reporting checkpoint, 26 November 2026

The preliminary N-terminal-relay PMF at 300 K contains 20 ps of accepted
production sampling in each of 57 umbrella windows (1.14 ns aggregate), with
3,150 configurations retained by the existing blockwise equilibration and
decorrelation procedure. Pointwise statistical uncertainty was assessed using
200 within-window resamples of complete 5 ps production blocks, retaining the
original blockwise selections and refitting MBAR for each resample. Histogram
edges and the sampled E–S reference bin were held fixed. All 43 occupied bins
were supported in all resamples. The largest sampled PMF value was
22.15 kcal mol⁻¹ relative to the fixed reference bin, with pointwise 95%
percentile limits of 20.64–23.33 kcal mol⁻¹ at that coordinate. These limits
describe a fixed-bin free-energy difference; they are not confidence limits
for a selected activation barrier. Basin-reference and sampling-convergence
criteria remain unresolved.

The visibly irregular profile coincides with uneven effective support.
For example, the bin centered at normalized path coordinate s=0.3125 contains
81 selected configurations but a bin-local weight effective sample size of
approximately 1.02; at s=0.1125, only five configurations are selected, with
weight effective sample size approximately 2.56. Statistical intervals are
conditional on the sampled configurations and do not exclude slower protein,
solvent or proton-transfer relaxation. The previous direct-seeded campaign
used 30 ps per window, 7,871 selected configurations and a different primary
coordinate and histogram grid; its smoother appearance alone does not imply
better convergence or permit a quantitative comparison of activation barriers.

In a selection sensitivity check, the same equilibration/correlation estimator
was applied once to each continuous 20 ps window rather than independently to
four successive 5 ps segments. The median statistical inefficiency increased
from 3.08 to 6.73 (maximum 11.73 versus 33.23); 2,840 configurations were
selected compared with 3,150 in the archived analysis. Thus the short-segment
procedure may miss correlations spanning segment boundaries, even though it
discards more configurations as putative equilibration. This diagnostic did
not replace the archived PMF or its selections; whole-window reanalysis and
time stability are required before quantitative barrier acceptance.

## Continuous-window sensitivity of the preliminary relay PMF

At the same 20 ps/window checkpoint, we concatenated the four continuous
production segments of each umbrella window before estimating equilibration
and correlation times. This procedure selected 2,840 configurations, compared
with 3,150 from segmentwise analysis. The largest sampled PMF value, located at
s=0.6875 under both analyses, increased from 22.15 to 25.32 kcal mol⁻¹ relative
to the same fixed E–S plotting reference. Using matched PyMBAR asymptotic
covariance estimates, the pointwise 95% intervals at that coordinate were
20.92–23.38 and 23.92–26.72 kcal mol⁻¹, respectively. The largest absolute
difference between profiles over their common occupied bins was
4.52 kcal mol⁻¹. Since both estimates use the same trajectory data, their
intervals do not constitute independent-replicate significance tests.

The continuous-window estimate used fewer configurations but had a larger
global target-weight effective sample size (29.29 versus 6.73). This does not
establish convergence: 46 of 57 windows retained fewer than 100 selected
configurations, and eight retained less than 5 ps after equilibration removal.
Two hundred within-window bootstrap refits of the globally decorrelated frames
gave pointwise percentile limits of 24.12–26.76 kcal mol⁻¹ at s=0.6875.
These resampling units differ from the archived five-ps block bootstrap;
matched analytical intervals were therefore used for the main comparison.
The sensitivity to selection policy and unresolved endpoint references
preclude a validated activation-free-energy claim at this checkpoint.

## Audit of extreme automatic production-data truncation

Time-series inspection of the eight latest automatic cutoffs showed that
extreme discards were not uniformly attributable to a sustained initial
transient. In activation window 002, the proton-coordination progress metric
had five-ps means of −0.0434, −0.0477, −0.0441 and −0.0461. Nonetheless, the
default lag-summation approximation selected a 17 ps discard; using full lag
summation moved the maximum cutoff to 6.85 ps. Window023 showed a substantial
initial change around5–6 ps and a later excursion; its 18.2 ps cutoff isolated
the final quieter suffix. Window024 exhibited a late change in the biased
acylation coordinate, so neither universal early truncation nor the automatic
latest-cutoff rule was established as a validated choice.

The latest-of-thirteen individual optima, followed by maximum-g thinning of
their common suffix, did not maximize the final common sample count. A
diagnostic 5 ps discard retained 18 configurations in window 023 versus 6 under
its automatic choice, and 18 versus 9 in window 002. These larger counts do
not establish equilibration and were not substituted into the reported PMF.
HG1 ownership and the two original N-terminal N–H bonds remained intact
throughout these eight trajectories, indicating that the cutoff-driving
fluctuations were not transfers of those protons between bonded owners.
The observed sensitivity requires further temporal and burn-in checks before
quantitative free-energy barriers can be accepted.

These observations qualify interpretation of the automatic discards. Chodera's
criterion maximizes an estimated effective sample count after truncation [R55];
it does not establish a unique physical equilibration time. Clark, Cole and
Michel showed that noisy correlation estimates can create spurious late
optima and that their tested truncation heuristics did not reliably identify
insufficient sampling [R56]. Their benchmark used synthetic series modeled
on protein–ligand free-energy data, providing methodological context rather
than system-specific validation for Nyl12 QM/MM trajectories. Our additional
latest-of-thirteen aggregation is a separate source of sensitivity identified
by the present audit. Neither the 18.2 ps automatic choice nor a universal
5 ps replacement is therefore established as a validated policy. The current
confidence intervals describe conditional statistical uncertainty and do not
resolve incomplete conformational sampling or uncertainty in burn-in selection
[R57].

## Direct-route umbrella sampling reanalyzed with the common protocol

The existing direct-NEB-seeded data were reanalyzed with the same fixed 5 ps
production discard, all-frame MBAR estimator, three-CV projection algorithm,
histogram grid, fixed E–S reference and contiguous-block bootstrap procedure as
the relay data. The 25 ps comparison retained 25,200 direct and 22,800 relay frames;
the full 30 ps direct analysis retained 31,500 frames. Bias energies were evaluated
using each campaign's actual native restraints. All 378 direct production
segments passed source, timing, distance and cross-bias parity checks.

| Sampling campaign | Production/window | Sampled maximum, kcal/mol | Pointwise 95% interval at maximum |
|---|---:|---:|---|
| Direct-NEB-seeded | 25 ps | 50.920 | Unresolved (194/200 draws support bin) |
| Direct-NEB-seeded | 30 ps | 50.947 | Unresolved (195/200 draws support bin) |
| N-terminal-relay-seeded | 25 ps | 24.807 | 23.938–25.553 |

The direct maximum occurs at s=0.4125 and contains only 15 frames, including no
additional visits in the final 5 ps. Missing-bin draws were retained as unsupported
and no finite interval was assigned. This is not evidence of a converged
50.9 kcal/mol activation barrier. At 25 ps, varying the initial discard from 0 to
10 ps changes the direct maximum from 49.229 to 52.670 kcal/mol; early/late halves
peak at 50.333 and 48.540 kcal/mol. The corresponding relay temporal maxima are
23.008 and 26.318 kcal/mol. Connected MBAR networks and relatively narrow bands
elsewhere do not remove this sampling sensitivity.

Direct-seeded sampling also visits the N-terminal proton channel: HG1 is nearest
to Nα in 2,763/25,200 retained frames at 25 ps and 3,463/31,500 at 30 ps. These are
biased-ensemble nearest-owner counts, not equilibrium populations or covalent
assignments. No such frames were censored. The comparison therefore does not
isolate a pure direct mechanism or establish a mechanism-specific rate ratio.

The historical 47.05 kcal/mol direct profile was a q_acyl marginal, whereas this
revision uses the same three-CV projection convention as the relay. Its changed
height must not be attributed to new MD sampling alone. Each path retains its
own reference geometry, so identical s labels do not mean identical molecular
states. The physical audit found identical topology bytes, QM/active selections,
charge and boxes; exterior positions agree within 6.32 × 10⁻⁶ Å. SCC tolerances differ
(10⁻¹⁰ direct, 10⁻¹² relay with tight population convergence). Current 3ob hashes
match the relay record, but historical direct-era hashes were not recorded.
These limitations preclude calling this a fully qualified matched kinetic
comparison. The ongoing relay sampling and queued 60 ps continuation were left
unchanged. Figures, numerical reports and source identities are in the direct
campaign's `plots/` and `analysis/fixed_burnin_20261126/` directories.

## Fixed-burn-in reanalysis of the 25 ps relay checkpoint

All 57 windows completed 25 ps of production. Native heating and equilibration
were excluded, and an additional fixed 5 ps production discard retained 20 ps
and 400 consecutive configurations per window (22,800 frames total). This
source-bound reanalysis used identical physical Hamiltonian, umbrella biases,
projection, histogram edges and E–S plotting reference. No configurations
were excluded by proton ownership or catalytic geometry.

The highest sampled PMF value was **24.807 kcal/mol at s = 0.6875**, relative
to the fixed sampled E–S reference bin. The pointwise 95% percentile interval
at that coordinate was **23.938–25.553 kcal/mol**, from 200 full MBAR refits
resampling contiguous 5 ps chunks within each window. This is a conditional
interval at the reported bin, not a confidence interval for a selected global
activation barrier. Each window provided only four such chunks. All 200
bootstrap fits remained connected; 43 of 45 populated bins had support in
every draw. Intervals were not plotted for the two incompletely supported
terminal bins. Resampling 2.5 ps chunks gave 24.151–25.449 kcal/mol at the same
coordinate, but agreement between these block sizes does not exclude slower
relaxation.

| Production burn-in | Highest sampled PMF, kcal/mol |
|---|---:|
| 0 ps | 24.458 |
| 2.5 ps | 24.546 |
| 5 ps, provisional primary analysis | 24.807 |
| 10 ps | 25.495 |

The temporal comparison was more concerning than the initial-discard
sensitivity: the 5–15 ps and 15–25 ps profiles reached **23.008** and
**26.318 kcal/mol**, respectively. The late profile differed from the combined
profile by up to 4.302 kcal/mol over common populated bins. Thus the relatively
narrow conditional bootstrap band did not establish temporal convergence.
The largest estimated g × frame interval was 5.326 ps in activation window026
for the HG1–leaving-N distance (5.192 ps in window025), further cautioning
against interpreting 5 ps chunks as fully independent. These estimates are
themselves uncertain on short records [R56–R58].

The global target-weight Kish ESS was 1.72 before accounting for temporal
correlation; the corresponding bin-local weight ESS values at the displayed
reference and maximum were 114.2 and 167.8. These measure different aspects of
reweighting concentration and are not independent-sample counts. No converged
activation free energy was accepted. The older automatically thinned estimate
had an unresolved basin-reference gate; its gate labels were not transferred
to the all-frame estimate. In the new histogram, the candidate E–S and acyl
minima lie inside the declared search intervals at s = −0.0125 and 1.0375,
respectively, with occupied neighboring bins. However, the acyl minimum has
only 27 frames and bin-local Kish ESS 1.24, and the next bin has one frame.
Numerical occupancy therefore does not establish a well-sampled product basin.
Its early/late minimum also shifts in position and free energy. Further sampling
and basin/initialization checks remain necessary; neither a low sampled
maximum nor narrow conditional intervals justify adjusting the model to
match a desired rate.

The existing US trajectories and queued extension to 60 ps/window were left
running. The new analysis publisher applies this fixed protocol at subsequent
complete 5 ps checkpoints. Evidence is in
`full_path_300K/analysis/fixed_burnin_20261126/fixed_25ps_analysis.json`, with
PMF, burn-in/temporal and block-size sensitivity figures in the campaign's
existing `plots/` directory. Methods and bibliography provide the statistical
and recent QM/MM precedents [R55–R60].
