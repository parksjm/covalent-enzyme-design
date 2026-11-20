# Nyl12–PA6 mechanism: decision roadmap

Status: working mechanism subproject plan, updated 2026-11-19. The overarching
`MASTER_PLAN_NYL12_COVALENT_REDESIGN.md` remains controlling for redesign;
this document is its **single forward-looking mechanism execution plan**.
The append-only action journal, run manifests and manuscript source ledger hold
execution history. Earlier dated mechanism plans are provenance, not queues.

## Active queue and next gates

### Launched — matched WT protonation ensembles (2026-11-19)

The user authorized all three WT conditions. The common-box campaign in
`qmmm/pa6_flipped_model30/10_wt_protonation_md/` now runs unrestrained GPU
MM before four independent 100 ns replicas per state, sequentially allocated
to two free Violet GPUs. All three actual preparation audits passed; full
MD acceptance and occupancy results remain pending. Read its README and
frozen orchestration status. These simulations test fixed-microstate
structural ensembles; they do not allow proton exchange or establish pKas.
Proton movement from cationic Thr N to either Asp would reach a proton
assignment already considered, although a different protein/water ensemble
can affect access to the reactive geometry. Do not equate this with finding
a new stoichiometric intermediate.

### Proposed — umbrella sampling from the converged direct NEB

The 35.80895 kcal/mol historical band is force-converged at cycle 477,
with ASP272/ASP274, neutral Thr N, 98 real QM atoms and charge -2. It
connects E–S to acylenzyme by direct OG1-proton delivery to leaving N.
A regular-band maximum is not a verified TS, but the images can seed US.
The current ASH274+water band is a different model/seed, not an acid-only
comparison; its original HG1 has already exchanged onto water.

NEB's Cartesian path spans the declared movable QM/MM atoms; it does not
force an umbrella coordinate choice. A native-Sander proposal uses independent
chemical CVs: q_acyl=d(C,N_leaving)-d(C,OG1), and
q_H=d(OG1,Hgamma)-d(N_leaving,Hgamma). Each window receives separate
harmonic terms about the image-derived centers. Monitor all four component
distances and exchangeable proton identities so degeneracy or competing
proton transfer is not hidden by a summed coordinate. Seed a local pilot
from the existing band, then add windows/orthogonal coverage based on overlap
and hysteresis; do not interpret an unsampled narrow ribbon as a complete
two-dimensional free-energy surface.

Alternatively, a path CV can express progress s and distance z from the
reference path using these reactive distances. This needs a reviewed native
implementation or explicit dependency approval; no PLUMED installation is
implied. Avoid restraining every QM Cartesian coordinate to each image.
NEB springs are removed. All applied umbrella/orthogonal biases must be
represented in MBAR reduced potentials; an extra confining wall that cannot
be unweighted defines a conditional PMF and must be reported. Preserve the
reference Hamiltonian and boundary conditions for a controlled PES/PMF
comparison. Production is at least 30 ps/window, with extension based on
overlap, autocorrelation/ESS and block/replica stability, not elapsed time alone.
No new US production was launched by this proposal. A lower free-energy
barrier is not assumed, and a stable TI1 is not required for concerted acylation.

Method context: Branduardi, Gervasio and Parrinello, JCP 126, 054103 (2007),
DOI 10.1063/1.2432340; Lonsdale, Harvey and Mulholland, Chem. Soc. Rev. (2012),
DOI 10.1039/C2CS15297E. These support CV/bias methodology, not a Nyl12 mechanism.


### Cationic Thr233 comparison and NEB checkpoint — 2026-11-19

Cationic Thr233 N termini with ASP272/ASP274 should be retained as a competing
WT hypothesis, not assumed superior or substituted across all mutants. Existing
PROPKA site-E estimates (8.93–9.70) motivate this comparison but omit PA6 and
are conformation/coupling sensitive. NH3+ cannot accept HG1 as a neutral
amine base would; its reactive pathway must involve prior deprotonation or
another base/relay. A minority neutral population could still carry catalytic
flux. No new microstate or MD campaign is authorized by documenting this option.

Two charged-N WT Model30 replicas already completed 100 ns with the same
PA6 RESP/GAFF2 release and ff19SB/TIP3P. Their Tier A counts are 487/10,000
and 601/10,000 (4.87%, 6.01%), versus baseline 37 and 238 and ASH274 1,677
and 0. The charged joint Tier/native screen yielded only 0 and 1 frame, with
no temporally supported seed; its one diagnostic frame subsequently failed
post-MM geometry. No charged QM/MM followed. Thus simple Tier A persistence
is more consistent between those two charged replicas, but improved catalytic
competence is unestablished. Reuse the cached diagnostics before new runs.

A bounded WT comparison would include neutral-N/ASP274 baseline, neutral-N/
ASH274 and cationic-N/ASP274, keeping Asp272 deprotonated. The last two have
the same dry proton count and formal charge, making proton location a useful
controlled variable. Existing boxes are not identical (118,829 versus 104,934
atoms; Na/Cl 123/119 versus 110/106), and dry partial charges differ by 0.004 e
from the neutral-Thr recipe's -0.001 e per copy. Preserve these differences;
do not call the historical simulations an identical-solvent proton-relocation
experiment. Validate terminal charge bookkeeping before new matched builds.
Original cached evidence: `docs/md_model30_charged_screen_20261111.md` and
`qmmm/pa6_flipped_model30/04_md_screen_20261111/`.

At 2026-11-20 02:12:36 UTC (Nov 19 local), the extended CG NEB reached cycle
260: physical profile maximum 43.642 kcal/mol, maximum force 8.093 times
target, RMS force 0.182 times target, and zero geometric image resets. It
exactly replayed the 79-cycle pilot before advancing. Uniform24 L-BFGS was
at cycle 426 / 51.405 kcal/mol / 59.842 times max-force target; uniform48 at
219 / 47.187 / 23.572; adaptive24 at 428 / 44.553 / 44.484. Original NEB
was at 734 / 49.799 / 71.076. All five engines were verified running.
None passes the maximum-force criterion; these are provisional PES maxima,
not verified barriers or PMFs. No calculation settings changed for this audit.

Actual native geometric image-reset counts were 50 for original, 5 for
uniform24, 11 for uniform48, 9 for adaptive24 and 0 for either CG. The old
monitor's `reset_count` counts BFGS memory resets and must not be interpreted
as this geometric-reset count. The explicit corrected audit and frozen plots
are under the existing `04_neb/plots/status_refresh_20261120_021148/`.


### Active — long CG identical-start replay (2026-11-19)

The user authorized a much larger CG budget, including starting over if needed.
The finished 79-cycle pilot has no usable optimizer checkpoint. The supported
same-start replay is running in `04_uniform_24_cg/long_run/`, under the same
`04_neb/spacing_tests/` campaign. Session: `nyl12_neb_cg_long_20261119`.
Only iteration ceilings changed: DL-FIND maxcycle 250→32,000 and Amber maxcyc
100,000→300,000. The resulting 224,000-evaluation budget is 128× the pilot.
All physical settings, thirteen control files and tolerances are unchanged;
input equivalence was checked by inverse replacement of those exact tokens.
The pilot and all other NEB jobs remain unchanged.

The first native profile reproduces all 24 pilot energies and Max/RMS gradient
values exactly. An independent tmux checker compares through cycle 79 and
records whether the replay subsequently advances beyond it. Full overlap is
pending; allow about 65 minutes to retrace the pilot, not to converge.
Current plots and movies publish under `long_run/plots/` and
`long_run/progress/latest/vmd/`; recursive laptop pull filters include them.
Read `long_run/README.md` and `replay_verification/last.json` before reporting
this as an exact optimizer restart or a scientifically converged path.


### Running — four parallel NEB numerical-spacing trials (2026-11-19)

**Early progress checkpoint:** uniform24/adaptive24 L-BFGS reached cycle80,
dense48 reached40, and the short CG diagnostic reached78. At approximately
comparable energy/gradient work, all four new bands have zero automatic image
resets; the original band has31. Profile maxima are39.8,55.1,42.0 and46.4
kcal/mol respectively; these are provisional physical-image maxima, not
validated barriers. Maximum-gradient/target ratios are113,81,64 and21; all
four RMS-gradient criteria pass, but no band has converged. CG currently has
the smallest maximum gradient despite a higher profile maximum; this is
promising numerical evidence, not a basis to rank physical mechanisms.
Its short1,750-evaluation budget is nearly exhausted by design. Continue the
planned matched-work checks before selecting a production continuation.
Snapshot: `04_neb/spacing_tests/plots/early_comparison_status_20261119.json`.


The user authorized four violet trials under the canonical
`09_water_relay_acylation/04_neb/spacing_tests/` directory. The retained original
NEB continues as context; it is not a matched restart control.

| Trial | Regular images | 13-control seed | Optimizer and scope |
|---|---:|---|---|
| `01_uniform_24_lbfgs` | 24 | Uniform QM arc length | Primary numerical baseline |
| `02_uniform_48_lbfgs` | 48 | Same uniform controls | Primary image-density comparison |
| `03_adaptive_24_lbfgs` | 24 | Double weight near source images 9–10 | Primary seed-detail comparison |
| `04_uniform_24_cg` | 24 | Same uniform controls | Short curvature-history diagnostic |

All four reuse the immutable cycle-268 band (`04_neb/progress/update_000075`),
periodic DFTB3/3ob at electronic temperature 0 K, QM102/54 real H/eight link
atoms/charge −1, ASP272/ASH274 and neutral Thr Nalpha. Frozen endpoints,
exterior, force tolerances and the complete-residue 12 Å/4,662-atom movable
region are unchanged. No chemical or positional restraint was added. Common
`maxstep=0.02` bohr limits each Cartesian step component. Native initialization
still uses piecewise Cartesian linear interpolation: these trials change its
control distribution or image count, not its interpolation formula.

The installed native code repairs bends **above 90° after the optimizer step
cap**, interpolating neighboring images; reducing `maxstep` does not cap those
repairs. In the original band, reset-associated maxima jumped
35.08→82.74 kcal/mol (cycles 244→245) and 32.83→75.15 (269→270). This is a
concrete numerical diagnosis, not evidence that a lower transient maximum
represents progress. L-BFGS remains primary; CG only tests the role of its
curvature history. CG has a 1,750-evaluation budget (`maxcycle=250` in this
adapter), roughly 79 full 24-image cycles; exhaustion is not convergence.
The three L-BFGS trials retain 112,000 evaluations each.

The 48-image initialization has one explicitly reviewed shared reacting H:
H20318 in image 15 is 1.521 Å from leaving N and 1.590 Å from OG1. The strict
retention gate remains false; `inputs/uniform/shared_proton_review.json` pins
this exact unconverged seed. All other images and spectator checks pass; the
exception neither accepts a stationary state nor adds a chemical restraint.
Per-trial `startup_precision_review.json` accepts native atom-map agreement at
1e-6 Å; original stricter 1e-7 Å failures remain preserved. Startup mapping
acceptance is separate from force convergence and chemical acceptance.

All four portable movies have published under `<trial>/progress/latest/vmd/`.
Native engines and individual coherent exporters persist in recorded tmux
sessions; comparison session `nyl12_neb_spacing_comparison_20261119` writes
`spacing_tests/plots/comparison.png` and source-hashed JSON/CSV/checkpoints.
Compare Max/RMS gradient-to-target ratios, reset counts/timing, coherent
chemical-coordinate gaps and elapsed work at 50/100/250 cycles where budgets
allow. With frozen endpoints, 24/48-image bands require 22/46 interior native
evaluations per full cycle; compare actual `Calculating the energy` counts,
not cycle counts or maxima alone. No validated barrier or improved convergence
has yet been established. Curvature and downhill connectivity remain later gates.

Authoritative current details: `04_neb/spacing_tests/README.md`, per-trial
configs/reviews/job records, and the latest section of
`docs/standalone_handoff_20261020.md`. Earlier queue entries below retain the
execution history; they do not supersede this active comparison.

### Running context — original ASP272/ASH274 NEB (2026-11-19)

The user authorized NEB. Current canonical stage is
`09_water_relay_acylation/04_neb/` within the ASH274 campaign. Native session
`nyl12_ash274_acyl_neb_20261119` and coherent-band export session
`nyl12_ash274_neb_progress_20261119` are launched. The24regular-image band is
seeded by13controls from38optimized stage endpoints plus12window025 trace
checkpoints. Exact unbiased endpoints/Hlabels/native placement were checked;
no controls dropped. Earlier draft seed-order error was rejected before engine
launch and retained under preparation_history. QM102/−1/54realH/8links, zero
Kelvin electronic temperature and4,662mobile atoms in12Å remain unchanged.
Maximum16000cycles implies112000energy evaluations in this installed adapter;
no climbing image yet. There are no extra reaction/positional restraints.

Watch native process/scientific records and `04_neb/progress/latest/`, not
unusable native numbered restarts. Curvature/downhill verification remains
separate. VMD's PRO1312display error is fixed in `progress/latest/`; original
coordinates unchanged. Laptop pull includes both latest bundles and NEB mdin.


### Overnight outcome — product retained, reverse scan complete (2026-11-19)

The unbiased acyl-enzyme candidate force-converged and retained the ester and
balanced proton inventory. All36reversewindows and a final unbiased E–S
release force-converged. The product is +11.06 kcal/mol relative to returned
E–S; window020 is the sampled physical maximum at +50.89 kcal/mol (bias
0.000108kcal/mol). Neither curvature nor a continuous MEP is verified.
Window24→25 drops18.59kcal/mol with a proton-occupancy change; inspect/refine
that interval before treating the profile as a continuous reaction path.
Thr Nalpha remains NH2 and the tracked transferring H is original waterH20318;
originalHG1 remains on water, so no Nalpha-water relay is established.
**Next bounded work:** retain the product endpoint, inspect/refine that jump,
and test window020 as a saddle candidate with curvature and downhill checks.
No endpoint/reverse-scan jobs remain active; the original NEB and the four
parallel numerical trials above are separate running work. Source/summary:
`09_water_relay_acylation/overnight_summary_20261119.json` relative to ASH274;
final inspection under the same campaign's `progress/latest/`.


### Authorized execution — balanced acyl endpoint before another path

The user authorized the revised sequence on 2026-11-18. The bounded screen of
all six existing 100 ns replicas is complete: 60,000 correlated frames include
18,958 Tier A/B frames with positive signed phi, but only 51 also satisfy the
unchanged native-Thr gates. None contains the specified preorganized one-water
Nalpha–water–leaving-N bridge. Counts by microstate are in
`09_water_relay_acylation/00_screen/bounded_screen_counts.csv`, relative to the
ASH274 campaign. This excludes a ready bridge in that screened subset; it does
not exclude water reorganization, multiple waters or other conformations.
No MD snapshot or alternative microstate was silently promoted.

The independent endpoint test uses the **already accepted water-QM E–S source**,
not a claimed positive result from that screen. It preserves QM102/54 real H,
charge −1, eight link H, the original WAT1576 oxygen, periodic DFTB3/3ob at
zero electronic temperature, and the complete-residue 12 Å/4,662-atom movable
selection. There are no added positional restraints; the exterior remains
frozen as in the accepted calculation. Product construction conserves every
atom/proton and the ASP272/ASH274 inventory. The source E–S MM topology is
**not** used as a classical force field for the acyl-enzyme: reacting atoms are
QM throughout. Repeating E–S MM minimization is unnecessary for this accepted
source.

The first ideal planar-ester construction failed contact checks: imposing a
120° ester angle rotates an entire polymer arm and creates severe overlaps.
The source approach angle OG1–C–acyl-Cadj is 83.97°. The bounded 285-position leaving-fragment search plus a 24-value acyl torsion
sweep yielded zero severe contacts, with maximum heavy/H moves 2.172/2.366 Å.
The selected +15° torsion about bonded C19920–C19919 retains the acute source-like
83.97° ester angle for subsequent optimization. This is a coordinate hypothesis,
not an optimized endpoint or a reaction path. TIP3P's topological H–H constraint is omitted only
from the builder's chemical graph; the topology itself is unchanged. Failed
construction artifacts and the contact diagnostic are retained.

Once a seed passes the contact/inventory checks, the scripted tmux sequence is:

1. Product-seed QM/MM preoptimization with explicit temporary harmonic terms
   for ester formation, amide separation and leaving-amine proton occupancy.
2. Remove every driving term and optimize for up to 2,500 accepted steps.
   Require the unchanged force gate and balanced acyl-enzyme/amine retention;
   report every real-QM-H owner and all spectator proton inventories.
3. Only after that gate, run 36 reverse windows toward E–S with the single
   coordinate `r(OG1,C) − r(C,Nleave) + 1 Å*[CN(OG1,H) − CN(Nleave,H)]`.
   All 54 real QM H contribute; attack, cleavage and proton motion can vary
   independently orthogonal to this scalar. Nalpha occupancy is unrestrained.
   No stable TI1 is required. This scalar can still be chemically degenerate,
   so a coordinate value or energy maximum alone never establishes success.
4. Release the reverse endpoint without bias, audit connectivity and inventory,
   then nominate a chemically valid saddle for separate curvature/downhill
   verification. Production PMFs and redesign remain downstream gates.

Entry points are `scripts/prepare_relay_acyl_campaign.py` and
`scripts/run_relay_acyl_campaign.py`; one logical campaign is
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/`.
The campaign `protocol.json`, stage records and `campaign_status.json` will
distinguish running, force convergence, chemical retention and process failure.
The live accepted-state movie/plots use `progress/latest/`, and the laptop pull
includes compact inspection artifacts. Verify launcher records before calling
a queued/prepared calculation running. The native workflow is launched in tmux `nyl12_acyl_endpoint_20261118`;
the accepted-state monitor is `nyl12_acyl_progress_20261118`. Inspect
`02_acyl/job/` and stage scientific records for current progress.

### Controlling revision — productive hydrolysis, not a required TI1 minimum

The user's mechanistic review request motivates this revision (2026-11-18).
**Success is a chemically balanced, kinetically credible product-connected
pathway and verified relevant transition state(s). A stable TI1 is optional.**
Literature supports both stepwise Ntn chemistry (Andjelkovic 2023, R41) and
amide acylation without a metastable tetrahedral intermediate in the calculated
path (Palermo 2015 FAAH, R48). “Concerted addition” can still yield a TI, so
identify the actual coupled bond changes (Nutho 2019, R46). Negoro's mature
NylC Fig. 10 is a BapA-guided proposal; precursor autocleavage evidence must
remain separate. Nitrogen inversion/proton access are explicit design concerns
(R44,R45), not merely oxyanion geometry.

Both outstanding calculations completed normally. The unbiased D release
force-converged back to E–S after 249 steps: OG1–C 3.211 Å, Thr-NH2/OG1H,
with Thr/water proton-label exchange. DFT and xTB embedded-QM gradients at D
both oppose additional attack while retaining the specified proton occupancy.
This is a negative result for this configuration/branch, not proof that no TI1
or catalytically viable pathway exists. D was still only partially added before
release. No new barrier, minimum or broadly validated Hamiltonian follows.

The missing chemical test is **leaving-group protonation plus C–N scission**.
Holding q_occ near 2.10 favors Thr-NH3+ and penalizes the normal net donation
of its extra proton onward unless another proton replenishes it. The pilot
therefore was not a complete test of an N-terminal general-base/general-acid
cycle. Coordinate audit also finds Nalpha–leaving-N 5.80 Å in E–S and 5.23 Å
at D; direct donation is not arranged. At D, QM-water O–leaving-N is 3.35 Å
and its nearest H is 2.60 Å. A reorganized water relay is a concrete hypothesis,
not an established path. Asp272 oxygens and the ASH274 proton are too distant
for direct donation in these snapshots. Nyl12 Phe113 also lacks the phenolic
OH of the corresponding proposed NylC residue Tyr146; its backbone NH remains.
Do not silently add Tyr chemistry or neutralize Lys156 to obtain a desired path.

**Bounded next sequence (proposed; no new jobs launched by this review):**

1. Screen existing equilibrated structures for a credible complete proton
   network to the **leaving amide N**, retaining useful carbonyl/oxyanion
   geometry. Prioritize an Nalpha–water–leaving-N relay if geometry supports it;
   inspect Asp-assisted alternatives only when they supply a reachable relay.
   This is distinct from requiring the earlier fixed-pocket OG1–Nalpha ring.
   Keep at most one primary route and one justified control.
2. Prepare and relax a chemically balanced acyl-enzyme/amine endpoint in the
   same atom mapping and explicitly declared microstate. Thr covalent catalysis
   is the best family-guided hypothesis, not a consequence of net hydrolysis
   stoichiometry; direct activated-water attack remains a lower-priority
   alternative if the covalent route is not supported.
3. Connect E–S to that endpoint with independently evolving attack, C–N
   cleavage and proton-transfer coordinates. Permit the catalytic Nalpha and
   water occupancies to change; do not lock in a TI1 protonation pattern or
   demand a minimum at a prescribed window. Allow nitrogen pyramidalization,
   inversion and polymer reorganization. Classify any TI only if it emerges.
4. Verify the chemically relevant saddle and downhill endpoints, then assess
   deacylation before advancing multistate design or production PMFs. A
   persistent TI1 is not required for LigandMPNN/redesign; a credible TS,
   reactant and product/intermediate ensemble is required.

Do not reopen another TI1-only grid or repeat Hamiltonian swaps on the same
geometrically incomplete branch. Percent polymer conversion in crude lysate
cannot determine kcat or a unique chemical barrier: adsorption/accessibility,
productive-state population, enzyme amount and turnover steps are confounded.
All completed negative results remain evidence. See manuscript R44–R48 and
`analysis/nyl12_mechanism_evidence_20261118/` for literature/coordinate provenance.

### Previous decision — complete pilot, test D release and electronic response

All four points are force-converged. The proton-label-invariant coordinate
retains the intended N3/OG0 state without water proton exchange; nevertheless
attack reaches only 1.854 Å and costs +42.763 kcal/mol relative to original
E–S. Carbonyl C–O is 1.292 Å, so D fails the established ≤1.75 Å attack and
≥1.30 Å tetrahedral C–O screens. There is no accepted TI1 or validated barrier.

The endpoint audit also identifies loss of one oxyanion contact: substrate O
to Phe113 N grows 3.161→4.110 Å, while Asn186 ND2–O is 2.837 Å. This motivates
testing both geometry persistence and the local electronic response rather
than increasing the bias or starting another scan. Evidence:
`occupancy_callback_pilot/endpoint_audit/{audit.json,summary.csv,supplement.json}`.

The user-authorized next tests are launched:

1. Remove **both** driving biases from D; optimize in the unchanged 12 Å
   active region with native gradient L-BFGS and the same force gate. Session
   `occupancy_D_unbiased_release`; output `pilot/point_D/unbiased_release/`.
   Audit all proton ownership and original bonds on completion. A returned
   E–S endpoint is a negative result for this starting geometry/microstate,
   not proof that the enzyme lacks a viable hydrolysis mechanism.
2. One exact-D ωB97X-D/def2-SVP/GFN2-xTB force comparison in the existing
   finite embedding, under `d_reference/` (session `occupancy_D_method_reference`).
   Inspect projected attack/proton-transfer and oxyanion-contact derivatives.
   Do not compare embedded-QM-only slopes directly to native total QM/MM
   slopes or promote a Hamiltonian on a single result.

**Next decision:** review both results before another path/sampling campaign.
TI1-like force convergence would nominate a stationary-point candidate;
minimum/TS verification requires separate curvature/connectivity evidence.

### Current repair — native gradient-based continuation, 2026-11-18

The line-search failure is explained quantitatively by ordinary MM pair
11272/11330 crossing the 10 Å nonbonded cutoff. The resulting ≈0.001 kcal/mol
energy jump exceeds the predicted smooth energy change by about 32,000-fold;
pair energy and force predictions match within native print precision. This
is not an SCC failure, a changed QM membership or evidence against activation.
`occupancy_callback_pilot/line_search_diagnosis/` retains reproducible evidence.

The user-authorized fix resumes from the final accepted A snapshot using the
installed DL-FIND gradient-based L-BFGS routine, with 0.05 Å maximum per-atom
steps and explicit history resets on nonpositive curvature. No dependency,
potential, protonation, spring, QM/active selection or convergence criterion
changed. Session: `occupancy_attack_pilot_continuation`; script:
`occupancy_callback_pilot/resume_pilot.py`. The original failed terminal records
are preserved in named history folders; native calls and accepted frames append
in place. Points B–D still require A force convergence. The same diagnostic
scope applies: no automatic unbiased release or claim of TI1/TS/PMF.

Repair verification succeeded: point A reached maximum biased gradient
0.109355 kcal/mol/Å in three additional steps (iteration 161), retaining
three Thr-N protons, deprotonated OG1 and intact QM water/ASH274. Point B is
now active. The accepted-step plot/movie monitor was restarted against the
same progress folder. This establishes successful numerical continuation,
not stability of the biased activated geometry after release.

### Latest execution checkpoint — point A incomplete, 2026-11-18

The four-point pilot stopped at A after 158 optimizer iterations/284 native
evaluations, with L-BFGS-B status 2 (`ABNORMAL` line-search termination).
Native/SCC checks passed, but maximum active biased gradient 0.277869
kcal/mol/Å exceeded the 0.12 criterion. Points B–D remain unstarted.
`pilot/completion.json` is `incomplete`; launcher exit 0 is process completion,
not scientific acceptance. Retain this attempt. Diagnose the final line-search
energy/gradient behavior before a recorded continuation from its final accepted
snapshot; do not loosen the tolerance or change chemistry silently.

### Current decision — verified callback and four-point pilot, 2026-11-18

The analytic occupancy bias now reaches the actual optimizer through
`ActiveBiasCallback`, using native Sander single-geometry energies and analytic
gradients. The two-iteration functional proof passed all seven gates: a known
0.20 kcal/mol energy addition, nonzero added force, correct gradient composition,
native energy conversion, optimizer execution, and requested geometry mapping.
It used four native evaluations; its iteration-limit exit was intentional and
does not establish stationary-point convergence. Evidence:
`08_activation_attack_grid/occupancy_callback_pilot/callback_proof/`.
Native NFE/PMD remains rejected for this DL-FIND minimization backend.

The authorized four-point diagnostic is launched in tmux
`occupancy_attack_pilot`, under the same grid's `occupancy_callback_pilot/pilot/`.
Starting at q_004_attack_001, it targets q_occ=2.10 independently of OG1–C
targets 3.1251784, 2.60, 2.10 and 1.60 Å (points A–D, serial). Harmonic
potentials are 0.5*kq*(q−q0)^2 + 0.5*kr*(r−r0)^2, with kq=100
kcal/mol per coordination-unit squared and kr=200 kcal/mol/Å². These smooth
occupancies are not formal protonation labels. Both component counts, all
54 real QM H owners, original QM bonds and attack distance must be inspected.
QM102/charge −1, 4,662 movable atoms, periodic electrostatics and 0 K electronic
temperature are retained. There are no additional positional restraints or
bounds. Each point permits 1,500 L-BFGS-B iterations/3,000 function calls and
must reach maximum active biased-gradient component ≤0.12 kcal/mol/Å before
continuation. Process failure or an unconverged point stops the pilot explicitly.
This is a bounded diagnostic, not a converged reaction path, PMF or barrier.
Inspection products are `occupancy_callback_pilot/pilot/progress/accepted_progress.png`
and `progress/vmd/view_path.vmd`, refreshed from accepted steps and selected
accepted restart snapshots. Plot ΔE is relative to the first pilot physical
energy, not the original E–S reference. The laptop pull preserves these paths
and excludes full native-evaluation artifacts.

Saved reference-gradient analysis is complete, without new QM evaluations:
`reference_checks/gradient_analysis/`. At the three sampled states, DFT and
GFN2-xTB embedded-QM gradients both oppose shortening OG1–C at fixed q_occ;
the strongest opposing slopes are about 60 and 55 kcal/mol/Å. Some small
proton-transfer slopes have opposite signs. GFN2-xTB therefore is not promoted
as globally validated. Native DFTB3 *total* QM/MM gradients are recorded
separately and cannot be compared directly to embedded-QM-only gradients.

**Next gate:** inspect the pilot's actual occupancy/connectivity and force
acceptance. Only a chemically suitable force-converged candidate warrants a
separately recorded unbiased release or further matched method checks. Do not
restart the labelled-H grid or umbrella sampling from a nominal coordinate
value alone. Earlier checkpoints below remain execution history.

### Previous checkpoint — occupancy audit and unbiased release, 2026-11-18

The authorized audit of all 25 points is complete: 54 real QM hydrogens per
point (1,350 assignments; fictitious link H excluded), with nearest-heavy
ownership, ambiguity thresholds and smooth coordination counts. Rows 1–3
keep HG1 nearest OG1, row 4 places it nearest Thr N, and row 5 exchanges it
with a water proton. No detected hydrogen escaped the real QM heavy-atom
region. Threshold-based assignments do not establish electronic bond orders
or protonation free energies. Evidence is in `08_activation_attack_grid/proton_occupancy_audit/`.

The separately authorized one-candidate diagnostic release completed with
normal engine/SCC/force convergence in
`08_activation_attack_grid/diagnostic_activation_release/`. It started from
q_004_attack_001 and removed both reaction restraints while retaining the
QM102/−1 setup and 4662-atom movable shell. It did **not** retain activation:
HG1 moved to water, water H1 moved to OG1, and Thr N retained only its original
two hydrogens. The endpoint is therefore E–S chemistry, 4.7757 kcal/mol lower
in native potential energy than the original E–S reference. Its geometry is
Tier B (d1/d2/d3 3.205/3.103/2.968 Å; attack angle 117.39°; signed phi
+102.52°). Do not promote it automatically over the original Tier A seed or
call it a stable activated intermediate. The helper's ligand_rmsd=0 field
is not a measured displacement from the original E–S reference. A 13-frame
stripped movie retaining WAT1576 is in `diagnostic_activation_release/vmd/`.

The proton-label-invariant definition is
`q_occ = CN(Thr N,H) − CN(Thr OG1,H)`, where each count is
`sum_H 1/(1+(r/1.25 Å)^6)` over the 54 real QM hydrogens. Both separate
counts must also be monitored: their difference alone is not a complete
chemical-state assignment. The reusable analytic value/gradient implementation
is `src/nyl12_covalent/proton_coordination.py`; 12 focused tests passed with
100% line coverage, including hydrogen permutation invariance and a custom
derivative check on synthetic coordinates. It agrees with the independent
audit for all 25 points to 1.2×10⁻¹⁵. These smooth counts are not integer
valences, formal charges or pKas.

**Backend preflight failed:** the installed Sander/DL-FIND minimization did
not initialize native NFE/PMD occupancy restraints. A deliberately nonzero
0.20 kcal/mol bias gave zero energy contribution, with no PMD/CV output.
See `occupancy_coordinate_preflight/attempt_2/scientific_acceptance_corrected.json`;
the original 0.40 expectation is retained with its half-k correction. No
Amber installation was changed. Before another occupancy-driven path, wire
the tested analytic bias into a supported energy/gradient callback and prove
that it contributes to the optimizer's energy and forces, or explicitly repair
the native initialization with equivalent verification. Merely adding native
`&pmd` inputs to current DL-FIND minimizations is invalid. No new driven path
or sampling campaign has launched.

The three-state finite-embedding reference calculations also finished. Relative
energies (E–S zero) for the partially activated/attack structures are
18.56/39.23 kcal/mol with DFTB3, 13.84/33.45 with ωB97X-D/def2-SVP, and
19.95/27.20 with GFN2-xTB. These energies alone do not validate a fast
Hamiltonian or locate a minimum; reaction-direction gradient interpretation
remains pending. Preserve the two current priorities: correct coordinate
integration, and interpretation of the saved reference gradients. Do not
restart the labelled-H grid or umbrella sampling.

### Execution result — completed grid, 2026-11-18

All 25 reported points and all five row sweeps completed. No point passed the
declared TI1 geometry gate. For rows 1–4, the last attack points reached
OG1–C = 1.887, 1.886, 1.860 and 1.783 Å, with physical relative energies
+28.34, +32.71, +36.62 and +54.58 kcal/mol. Energy increased along attack
in every row; these are constrained PES observations, not barriers or proof
that no TI1 basin exists. The 1.50 Å attack target was not achieved.

The largest-q row is a **proton-label exchange branch**: Thr HG1 binds to
WAT1576 O (0.986 Å in its first grid point), while water H1 binds Thr OG1
(0.980 Å). Thr remains an alcohol and water remains water; Thr N is not
protonated by this exchange. Its −4.54 kcal/mol first point is a rearranged
reactant geometry, not a favorable activated intermediate. The labelled-H
coordinate therefore did not sample complete intended N-terminal activation
over its full nominal range. Any follow-up must distinguish chemical proton
occupancy from atom labels, e.g. through smooth O–H/N–H coordination variables;
do not merely extend this degenerate coordinate or increase its spring.

The E–S three-method reference finished normally, but selected-state exports
stopped because Amber truncates its Gaussian executable field to 20
characters. A real Sander export preflight verified the repair: a short local
wrapper plus the required Gaussian template. Session
`ash274_reference_continuation` now runs only the two unfinished activated/
attack references concurrently, reusing completed E–S and DFTB3 results.
Original failures remain under each state's
`extern_export_failed_pathlength/`. No TI1-like state was available for the
third selection. No unbiased release has run; the reference-review gate is
still pending. The next decision is to inspect the completed matched relative
energies/gradients and then select whether either actual candidate warrants
release. Do not launch a replacement pathway campaign before that decision.

### Controlling decision — after completed water-relay scan, 2026-11-17

**Current state: no active calculation in this path; no validated TI1 or
first-order acylation TS.** The 25-window water-to-Asp272 scan completed with
force-converged windows and continuous monitored reactive distances. Its final
physical potential energy is +54.0235 kcal/mol relative to its own E–S seed
(bias-inclusive +58.1547; bias +4.1312). OG1–C reached 1.852 Å despite a
1.650 Å restraint target. The endpoint is `endpoint_negative_or_incomplete`;
the conditional unbiased release did not run. Both intended proton transfers
occurred: Thr HG1 moved to water and water H2 moved to Asp272. ASH274 retained
HD2, leaving **both Asp residues protonated and Thr N neutral**. This particular
charge redistribution is energetically unfavorable along the sampled scan.
It does not establish the absence of a tetrahedral basin elsewhere on the PES.
Evidence: `qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/`
(`protocol.json`, `window_025/acceptance.json`, and `energy_ledger.csv`).

**Immediate objective:** discriminate an inadequate prescribed path from
absence of a persistent tetrahedral intermediate in this microstate. A stable
TI1 is a hypothesis, not a prerequisite imposed on the mechanism. Experimental
hydrolysis establishes overall chemistry but does not identify a microscopic
route, intermediate lifetime, or the correctness of the modeled binding pose.

1. **One bounded alternative on the same accepted E–S seed.** Prepare a coarse
   5 × 5 relaxed PES map with independent coordinates
   `q_N = d(Thr929 OG1,HG1) − d(Thr929 N,HG1)` and
   `r_attack = d(Thr929 OG1,substrate carbonyl C)`.
   Suggested initial ranges are the measured E–S q_N (~−1.78 Å) to +1.05 Å,
   and r_attack 3.10 to 1.50 Å. Freeze the exact grid and restraints in a
   reviewed protocol before execution. This tests direct N-terminal activation
   and different timings of attack without tying the two targets together.
   Retain the same 102 QM atoms, −1 total QM charge, water, ASH274/ASP272
   starting microstate, topology, and 12 Å movable shell. Treat any unexpected
   proton relocation as a distinct chemical branch requiring review. Preserve
   full proton accounting; water and the Asp pair may reorganize. Use the
   grid to find candidate valleys, not a barrier or a proof of a minimum.
2. **Validate only informative points.** At most four representative E–S,
   activated, attack, and tetrahedral-like geometries receive matched-embedding
   DFTB3/GFN2-xTB/ωB97X-D checks, reusing previous results where identical.
   Inspect reaction-direction forces as well as relative energies. Keep the
   point-charge membership and boundary treatment consistent in cross-method
   comparisons. The existing benchmark found DFTB3's relative TI1-like error
   +7.42 kcal/mol versus DFT and xTB's direct-PT error +8.02; neither is a
   validated universal replacement. Stop fast-method production if the new
   reaction-region comparison disagrees qualitatively with the reference.
3. **Release candidates, then connect endpoints.** Release reaction restraints
   on at most two plausible candidate basins. Require force convergence and
   retained chemistry. If retained, refine the valley with a string/NEB seeded
   from the map and then verify its saddle(s). If candidates collapse, one
   bounded restrained QM/MM equilibration followed by independent releases may
   test environmental reorganization; declare sampling and restraints in
   advance. A constrained basin alone is not a stable TI1.
4. **Explicit decision if TI1 remains absent.** Seek an E–S → acylenzyme path
   allowing attack, leaving-group protonation, and C–N cleavage to couple;
   do not insist that it pause in TI1 or reuse only the strained direct
   four-center route. Use the observed proton network to select the relay.
   If the grid/reference checks still leave no credible route, make one
   targeted QM-region/protonation sensitivity test (e.g. the documented
   Lys156/Gly273 expansion), or obtain external QM/MM review of the setup.
   Avoid another broad collection of paths with untested assumptions.

**Work limits:** this is the next proposed calculation sequence; updating this
document does not mean these calculations have launched. Keep the old direct
route and its unresolved Hessian index as a comparator, defer further costly
full-shell eigensolver attempts and umbrella sampling until a credible path is
selected, and preserve all prior negative results. Do not shrink the movable
shell. Refine the grid only where a valley or discontinuity gives a specific
reason. A second starting conformation is a robustness check after the first
mechanism passes, not another simultaneous campaign. Success means connected,
chemically assigned endpoints, stationary gradients, a reproducible index-one
TS and downhill connectivity, followed by selected DFT checks. PMFs follow
that gate. No computational plan guarantees those results for the present
model; each step must resolve a named uncertainty.

The historical updates below remain evidence. This dated decision supersedes
their outstanding next-run instructions where they conflict.

**Execution status — 2026-11-18:** The bounded activation/attack grid is
**RUNNING** as five independent tmux row workers. Its 25 reported points use
135 bounded optimization stages. The same QM102/−1 model, 4662-atom movable
shell and periodic Hamiltonian are retained. The E–S ωB97X-D/def2-SVP force
reference queue is now running with up to two eight-core DFT workers; it uses
a separately declared nonperiodic full-PC sensitivity model. Gaussian runs
directly on the exported 110-atom QM/link geometry and 104,824 adjusted charges
to avoid its field-query centre limit. The earlier Sander/EXTERN exit-0 result
failed Gaussian's termination/gradient gates and is rejected. The progress
monitor refreshes a measured-coordinate energy plot and per-row VMD bundles.
The candidate-release watcher requires explicit comparison review before
running at most two releases. Launcher records, frozen inputs and logs are under
`08_activation_attack_grid/job_row_N/` and `reference_checks/`; running status
does not establish force convergence or a retained intermediate.

**Update — corrected-shell ASH274 E–S accepted, 2026-11-17.** The matched
ASP272/ASH274 seed now has a force-converged, accepted Tier A DFTB3/MM
endpoint in `qmmm/pa6_flipped_model30/07_ash274_comparison/06_qmmm_optimization/replica_1_tier_A_seed/corrected_shell/`.
It retains ASH274 HD2 and neutral Thr233 N, with d1/d2/d3 =
3.07665/3.06706/2.94462 Å, attack angle 116.807°, and signed phi +97.925°.
The corrected movable shell contains 4,650 atoms; water 1576 is now movable
and its oxygen shifted 0.1465 Å. This corrects the inherited shell that froze
111 residues within 12 Å; the earlier accepted optimization was an incomplete
relaxation and is superseded. The final audit
`corrected_shell/water_1576_full_hbond_audit_final.json` finds
Thr233 OG1–HG1···water O at 167.9° and H···O 1.811 Å, and water
O–H2···Asp272 OD1 at 173.4° and H···OD1 1.533 Å. It does not establish a
bridge to Thr233 N or substrate leaving N. This local geometry is not evidence
of water occupancy or reactivity. A one-water-in-QM (102 atoms, charge −1)
reoptimization completed in tmux `ash274_water_qm_reopt` in the sibling
`water_qm_reoptimization/`. It started from the corrected-shell converged
endpoint, preserved the −1 QM charge and eight link cuts, and expanded the
12 Å movable region to 4,662 atoms to include four boundary waters. Sander
exited normally and DL-FIND converged at cycle 248. The endpoint retained
Tier A geometry (d1/d2/d3 = 3.100/3.063/2.959 Å; attack angle 116.44°;
signed phi +98.75°) and the declared Thr, water and Asp274 covalent proton
inventory. Thr HG1···water O is 1.932 Å at 165.5°, and water H2···Asp272
OD1 is 1.730 Å at 172.7°; see `scientific_water_assessment.json`. There is
still no bridge to the substrate leaving N (water O···N = 4.107 Å). This is
a stable *local geometry*, not evidence of water occupancy, proton transfer,
or a viable barrier. Next test one bounded water-assisted Thr activation and
OG1 attack route toward TI1, with explicit proton/geometry gates and a
separate unrestrained endpoint release. Compare with the direct route only
after identifying chemically continuous states; do not broaden the campaign.

1. **Asp-pair microstate comparison:** all six matched 100 ns, 300 K
   replicas completed. Tier A counts per 10,000 frames were 37/238 for
   ASP272/ASP274, 3886/1724 for ASH272/ASP274, and 1677/0 for
   ASP272/ASH274. Site-1 median closest Asp-pair O–O distances were
   4.89/4.90, 2.56/2.56, and 2.73/2.65 Å, respectively. The large
   between-replica variation, including zero Tier A frames in one ASH274
   replica, prevents ranking protonation microstates from frame counts
   alone. Next inspect ligand register, the 10 ns block series, and
   proton/water geometry before choosing at most one follow-up snapshot
   per microstate. Fixed-charge MD cannot exchange the Asp proton.
   A three-frame, protein-only PROPKA 3.5.1 screen is recorded in
   `docs/propka_asp_pair_pilot_20261116.md`; it does not select a
   protonation microstate. Assess ligand/water sensitivity before choosing
   the next QM/MM E–S model.
   The matched 10 ns block and water-oxygen analysis is in
   `docs/asp_pair_md_water_screen_20261116.md`. Water proximity is common
   in one ASH274 replica but the inspected candidate lacks a complete
   hydrogen-bond relay.
   That two-leg check has now completed: only 1, 6 and 32 Tier A frames
   in the three screened replicas pass both hydrogen bonds to Thr N,
   with at most three consecutive 10 ps samples for one water. None
   passes both legs to the leaving-group N while Tier A. Treat an
   ordered water relay as unproven and screen the remaining substrate
   register/protonation gates before choosing QM/MM seeds.
   Compare Asp306/Asp308 pKas from mature NylC 5XYG with **uncleaved**
   precursor controls 7YU0 and 7YU1, explicitly representing the new
   Thr267 N terminus only in the mature structure. These precursors carry
   maturation-blocking substitutions; differences cannot be assigned solely
   to cleavage. The first comparison completed with process exit 0, but a
   malformed `TER` record prevented recognition of mature Thr267's free
   N terminus. An atom-identical corrected input restores N+267 and removes
   its warning. Corrected PROPKA predictions for mature 5XYG are N+267
   6.35–6.55 (default) or 4.59–4.60 (coupled alternative), Asp306
   2.19–2.21, and Asp308 5.49–6.49 (default) or 7.24–8.45 (coupled
   alternative). The earlier mature pKas must not be used. The corrected
   audit is `analysis/negoro_2023_protonation/terminal_pka_20261116/before_after.json`;
   the corrected scripted run passed end to end under
   `analysis/negoro_2023_protonation/crystal_propka_20261116/attempt3/`.
   the precursor comparison is under
   `analysis/negoro_2023_protonation/crystal_propka_20261116/attempt2/`.
   Precursor predictions vary by chain and coupled-state treatment. These
   are structure-dependent predictions, not measured protonation assignments.
   Compare standard and coupled-site PROPKA predictions, then use an
   isolated MCCE4 trial on a small, matched structure set if installation and
   input validation succeed. Neither method establishes protonation by itself.
   The isolated mature-5XYG MCCE4/DelPhi trial finished, but its pKa bounds
   (both Asp <0; Thr267 N terminus >14) are **scientifically rejected**:
   MCCE generated 13 charged termini at internal missing-coordinate
   boundaries. See `analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/scientific_acceptance.json`.
   Do not use these bounds to assign an Asp state or run the mutant precursor
   controls with the same input preparation. Next preflight one complete
   Nyl12 Tier A protein model with MCCE step 1 only and require exactly its
   true mature-chain termini; launch a full pKa calculation only if that
   structural gate passes. Independently compare one matched Tier A
   ASP272/ASH274 E–S seed with the previous ASP/ASP and ASH272/ASP274
   simulations using identical QM/MM atom mapping, charge and relaxation
   gates. This is a stability test, not an absolute pKa determination.
   The complete Nyl12 Tier A MCCE4 step-1 preflight also **failed**:
   it retained all eight mature protein chains but assigned zero NTR/CTR
   groups and emitted 2,840 `Error!` diagnostics. No full pKa run was
   launched. See `analysis/nyl12_mcce_step1_preflight_20261116/neutral_replica_02_frame284_attempt2/scientific_acceptance.json`.
   Stop MCCE4 troubleshooting for this mechanism campaign unless a tested
   structural-input fix becomes available. Use the matched microstate
   stability comparison and, if needed, a targeted free-energy calculation
   to estimate proton occupancy; static PROPKA remains a sensitivity screen.
   Defer constant-pH MD until these cheaper screens identify a specific,
   consequential Asp/Thr-N microstate ambiguity and a validated protocol can
   handle the tetramer, termini and nonstandard PA6 substrate.
2. **Direct-route saddle:** 0 K DFTB3/MM dimer refinement from the
   converged NEB maximum finished at cycle 99 without SCC failure and with
   native force convergence. Its geometry matches the earlier downhill-tested
   candidate within 0.006 Å RMSD over the movable region. The bounded
   block-LOBPCG probe finished, but its approximate eigenvalues
   (−0.227, +0.033 Eh/(Å² amu)) have residual norms 0.051/0.061,
   far above 10⁻⁴. The Hessian index is unresolved. Require exactly one
   reproducible negative mode and two-sided endpoint connectivity before
   a true-TS claim. Do not report the
   35.81 kcal/mol NEB image maximum as a validated activation barrier.
   A bounded operator check now shows exact repeatability, <1% step
   sensitivity and 0.37% bilinear asymmetry; see
   `docs/ts_hessian_operator_diagnosis_20261116.md`. Gross force noise
   is unlikely to explain the stalled solver. Audit preconditioning or
   focused mode seeds before another full-region eigenvalue search. The next
   bounded test should recycle approximate Ritz vectors, improve conditioning
   without shrinking the 12 Å movable region, and compare residual-versus-
   gradient-evaluation curves against the previous probe. Accept a first-order
   saddle only after a reproducible single negative mode with converged
   residual and two-sided downhill connectivity from this same TS structure.
   A first warm-start pilot reused the two saved **partial** Ritz vectors and
   completed 38 HVPs/76 native gradients. Residuals improved to
   4.13×10⁻⁴ (negative candidate) and 2.20×10⁻² (next candidate), still
   above the 10⁻⁴ gate. Thus the Hessian index remains unverified; a
   better-conditioned follow-up is needed before a true-TS claim. The pilot,
   failures and acceptance record are under the direct-route TS tree's
   `mode_probe_block_warmstart_shortpath_20261116/` directory.
   A six-vector lowest-algebraic pilot then exhausted its 60-HVP budget
   before returning Ritz pairs; process exit 0 does not imply scientific
   acceptance. Its inputs and failure record are under
   `mode_probe_block6_warmstart_20261116/`. Require the eigensolver to save
   intermediate Ritz/residual checkpoints at its budget cap, then run one
   budgeted wider-block continuation. If the second low mode still fails
   to converge, diagnose a preconditioner or alternative eigensolver before
   spending more gradients; do not call the current structure a verified TS.
   The final six-vector run returned one negative and five near-zero positive
   *estimates*, but all full-Hessian residuals were 0.0076–0.0156 versus
   the 10⁻⁴ gate; index one is unresolved. Test numerical restart/damping or
   a positive damped preconditioner on synthetic clustered modes first,
   then allow one bounded native pilot if it actually improves convergence.
   A scalar Hessian shift cannot untangle eigenvectors because it preserves
   their spacings; adding physical restraints would change the PES and
   cannot certify the unrestricted saddle. Prioritize the triad-assisted
   mechanism if this high-barrier direct route remains difficult to verify.
   A 12-gradient precision gate on the unchanged system found exact repeated
   HVPs, but 0.025 vs 0.0125 Å√amu products differed by 8.44×10⁻⁴ and
   1.79×10⁻³ Eh/(Å² amu), and bilinear asymmetry worsened at the smaller
   step. The nominal 10⁻⁴ residual target is not numerically supported by
   this central-difference operator. Do not relax the gate arbitrarily;
   test a higher-order finite-difference HVP or analytic/mixed Hessian method
   on a bounded set of directions and calibrate the tolerance before another
   full eigensolver run. Evidence is in the direct-route TS tree's
   `hvp_precision_gate_20261116/attempt2/`.

Full commands, hashes and exact paths are at the top of
`docs/standalone_handoff_20261020.md`. Do not open a new mechanism campaign
until these gates are assessed.

## Question and boundaries

Determine a chemically credible, reproducible route for WT Nyl12 cleavage of
PA6 L4 at the L2-producing bond and obtain stationary E–S, intermediate,
transition-state and acylenzyme structures useful for multi-state redesign.
Keep the tetramer and explicit solvent. Use the flipped, positive signed
Herschlag-phi pose as the first test, while retaining the opposite pose as a
bounded alternative. Treat substrate hydrolysis separately from NylC/Nyl12
autocleavage. Do not expand to PA66, poly(MAXA), mutants or a full polymer
binding free energy in this mechanism campaign without a new decision record.

The working chemical comparison keeps the Thr233 N terminus neutral and tests
ASP272/ASP274, ASH272/ASP274, and ASP272/ASH274. The Negoro NylC protonation
proposal is an inference, not established Nyl12 protonation. Test a direct four-membered
OG1-HG1-to-leaving-N route against one triad/water-assisted route. Add another
route only when geometry, conserved residues, or a failed chemical gate gives
specific evidence. The experimental assay at 75 °C and uncertain Tris pH at
that temperature cannot be equated to the 300 K MD elementary rate.
The initial neutral Thr233 N terminus is also a hypothesis. If it is NH3+,
it lacks a lone pair and cannot accept Thr HG1 as a general base in that
microstate. Then test a carboxylate- or water-assisted Thr activation, or
transient N-terminal deprotonation, with an explicit proton/charge inventory.
An NH3+ group might instead donate a proton to the departing substrate N
during tetrahedral-state collapse; favorable electrostatics alone do not
establish this catalytic role. Do not carry a neutral-N reaction coordinate
unchanged into an NH3+ calculation.

## Current evidence and what it does not prove

* Two ASH272/neutral-N 100 ns replicas completed. The geometry screen found
  3,886 and 1,724 Tier A frames versus 37 and 238 in the earlier ASP272
  replicas. The three campaigns share a static model30 protein/substrate pose
  and 300 K protocol, but solvent placements and random seeds differ. This
  is a promising association, not a causal protonation or pKa result.
* The first direct-route NEB stopped unconverged at 500 cycles with a
  ~70.6 kcal/mol outlier. A later 24-regular-image, trace-seeded NEB converged
  at cycle 477 and has a 35.81 kcal/mol image maximum relative to its first
  image. Its endpoint identities and local image continuity pass; a full
  first-order saddle is unverified. The older dimer point is ~46.4 kcal/mol
  above a separately minimized E–S reference, so these are not comparable
  activation barriers or PMFs.
* All four ASH272 charge-audited DFTB3/MM E–S optimizations converged with
  99 real QM atoms, eight links and charge −1. All four are geometric Tier A,
  but each transfers the prepared Asp272 proton to Asp274. Treat this as a
  microstate observation under the tested Hamiltonian, not proof of in-vivo
  protonation. Do not label their endpoints ASH272-protonated E–S minima.

## Execution gates, in order

1. **Resolve the E–S microstate comparison.** The six matched 100 ns
   replicas, Tier A/10 ns block analysis, and initial pKa and water screens
   are complete. Confirm the PA6 register and Thr orientation in candidate
   frames, test whether including PA6 or alternative pKa methods change
   the Asp/Thr interpretation, and select at most one accepted
   representative per microstate for follow-up QM/MM. Do not multiply
   reaction paths by six or assign protonation from static PROPKA alone.
2. **Validate electronic setup.** Parameterize the QM/MM input by declared
   protonation and charge; fail on mismatch. Keep the predeclared 12 Å movable
   region initially. Check charge, link cuts, periodic embedding, explicit
   waters, and MM-only versus QM/MM endpoint geometry. Stop a candidate if it
   loses the PA6 register or catalytic geometry during unbiased relaxation.
3. **Choose two mechanisms, then screen cheaply.** From a chemically accepted
   E–S microstate, build chemically explicit acylenzyme and possible tetrahedral
   states; equilibrate the surrounding MM degrees of freedom if an endpoint
   collapses, with the QM state restrained and all restraints recorded. First
   test triad/water-assisted proton transfer and attack. Retain the direct
   route as a comparator until its TS and energy are verified. Use restrained driving
   only to construct seeds; it is not an energy barrier.
4. **Converge a path before TS refinement.** Optimize both endpoints without
   reaction-coordinate restraints. Run regular NEB/string with smooth,
   geometry-checked images and enough cycles to meet predeclared force and
   displacement tolerances; only then use climbing-image refinement. Keep
   image-wise energies, forming/breaking distances, proton inventory and VMD
   frames updated. If the band stalls, diagnose its highest-force image and
   missing coordinate, then change **one** element of the path setup. Do not
   equate a 500-cycle stop or shoulder with convergence or a new minimum.
5. **Verify stationary points.** Refine each candidate TS with DL-FIND or an
   appropriate mode-following method at the screening Hamiltonian. Require
   stationary gradient, exactly one relevant negative Hessian mode, and
   downhill/IRC connectivity to the intended minima. The 4,165-atom numerical
   Hessian eigensolver is currently inefficient; diagnose precision and
   eigensolver strategy before merely raising its evaluation cap. Perform
   sparse hybrid-DFT/MM checks of key geometries and relative energies. If a
   fast Hamiltonian changes state ordering or proton identity relative to
   DFT, do not promote its PMF to mechanistic evidence.
6. **Extend to deacylation only after acylation is credible.** Add water attack,
   the second tetrahedral state, product-forming transition state(s) and final
   products. For redesign, use a consistent ensemble of E–S, acylation TS,
   acylenzyme and deacylation TS when available; an early design pilot may use
   validated acylation states but must be labeled incomplete-cycle design.

## Free-energy decision

Umbrella sampling is a **second-phase mechanistic result**, not a substitute
for locating a viable path. Launch a small 2D pilot only after a fast method
passes the chemistry/DFT checks and the relevant proton and bond coordinates
are identified. Check neighboring-window overlap, equilibration, forward/
reverse hysteresis, block stability and replica variation; extend windows
based on these diagnostics rather than a fixed 30 ps rule. If these gates pass,
compute an acylation PMF, then consider deacylation. A PES TS and a PMF
maximum answer different questions and should be reported separately.

## Manuscript scope and structure

The mechanism paper should include the **released PA6 parameter set and its
validation as a concise Methods foundation**, with versioned citation to
`/data/jpc/polyamide-amber-ff`; detailed modular residue derivation, broad
polyamide coverage and parameter benchmarks belong in that separate force-
field paper/SI. Include WT MD for E–S pose/protonation selection and sampling
limits. Include mutant MD only if it answers a predeclared mechanistic test;
otherwise leave redesign validation for a subsequent paper. The central
Results should be competing QM/MM paths, stationary-point validation, the
role of protonation/water, and DFT sensitivity. Include umbrella sampling in
the main paper if the convergence gates pass and it changes or strengthens the
mechanistic conclusion; otherwise report a bounded pilot as a limitation, not
a headline free-energy barrier. A full chemical-cycle claim requires both
acylation and deacylation, not only TS1.

Maintain `manuscript/INTRODUCTION_DRAFT.md` for prior NylC/Nyl12 work,
experimental motivation and competing hypotheses; `METHODS_DRAFT.md` for
executed protocols (main-text concise description plus SI detail); and
`RESULTS_DRAFT.md` for claims linked to exact run manifests and negative
results. Keep references in `REFERENCES.md`, claims/artifacts in
`SOURCE_LEDGER.md`, and decisions here. Move only mature, accepted findings
from the run ledger into manuscript prose. Do not paste long run logs into the
manuscript.

## Student MD comparison arm — 2026-11-19

A transferable handoff is in `student_handoff/nyl12_ash274_md/START_HERE.md`.
The immediate hypothesis is a matched WT protonation-state effect, followed by
rebuilding the prior eight mutants with the same PA6 RESP/GAFF2 release.
ASH274 has 8.385% pooled Tier A occupancy versus 1.375% for ASP/ASP, but its
replicas are 16.77% and 0%; joint native-Thr gates do not show a corresponding
advantage. This is a testable hypothesis, not an established catalytic benefit.
Use independent 100 ns replicas, matched 300 K protocols and replica-level
uncertainty. Preserve the historical site-AE comparison and implement validated
four-site tracking before ranking LigandMPNN designs. The old mutant panel
uses legacy AM1-BCC/GAFF substrate parameters and cannot isolate protonation.
No new MD was launched while preparing the handoff. Three WT microstates are
portable prepared starts; mutant coordinate references still require rebuilding.
This MD subproject does not replace the current mechanism/TS priorities.

## Operating rules

### Current bounded water-relay diagnostic (2026-11-17)

The accepted ASH274/ASP272 E–S endpoint with WAT1576 in the 102-atom QM
region is the sole starting structure. A 25-window DFTB3/MM scan drives the
composite Thr929 OG1–HG1 / water / Asp272 OD1 proton relay and OG1 attack on
the substrate carbonyl. The remaining MM atoms in the established 12 Å shell
are movable; no global positional restraints are added. Each accepted window
updates a stripped VMD trajectory and energy/distance plots. The runner stops
on failed convergence, proton identity, or geometry continuity. Only an
accepted endpoint triggers an unbiased minimization. This route would yield
Asp272-COOH and neutral Thr N: call it an **alternative TI1-like state**, not
canonical TI1, whose required Thr beta N is NH3+. Neither a biased scan nor
an unconverged endpoint establishes a stable intermediate or transition state.

The first execution stopped at window 5 on a false continuity rejection of
a distant ASH274–water separation; its reactive geometry and proton inventory
passed. Revision 5 preserves that original decision, checks only the seven
declared reactive distances for continuity, and resumes in place from window
6 after validating all completed output hashes. The path remains a diagnostic
until the endpoint and unbiased release pass their separate gates.


One active reference system, at most two competing routes, and one next
scientific decision at a time. Each queued run has a one-page hypothesis,
atom mapping, charge/proton inventory, method, restraints, input hashes,
acceptance gates and stop condition. Script and run in tmux; check only
meaningful milestones. Archive failed attempts without renaming them as
successful. Use a short weekly dashboard: completed, running, blocked,
next decision, and the single figure or structure the PI should inspect.
Change this roadmap only when a dated decision record cites new data and the
specific gate it changes. Avoid new calculations merely because hardware is
idle.

Selected methodological context: Negoro et al., *FEBS J.* 2023,
DOI 10.1111/febs.16755; Kästner and Sherwood, *J. Chem. Phys.* 2008,
DOI 10.1063/1.2815812. The manuscript reference ledger controls final
bibliography and claim verification.
