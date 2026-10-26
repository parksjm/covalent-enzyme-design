# Asp272 two-stage HG1 shuttle pilot

User explicitly authorized this construction hypothesis on 2026-10-25.
Active run: `output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/`.
tmux: `nyl12_asp272_relay2_20261025`. Read status.json and per-window records
for progress; initial native optimizer evaluations were verified at launch.

Stage 1 uses 17 coupled targets: shorten OG1–C54 from its initial distance
to 1.55 Å while increasing d(OG1,HG1)−d(Asp272 OD1,HG1) from its measured
initial value to +0.8 Å. OD1 is selected because it accepts the observed
hydroxyl hydrogen bond. All other QM coordinates remain free.

Stage 2 starts only after a converged, geometrically plausible covalent
intermediate with HG1 on Asp272. It holds the attack target at 1.55 Å and
uses 17 targets for d(Asp272 OD1,HG1)−d(terminal N,HG1), ending at +0.8 Å.
The original HG1 identity and both original terminal N–H bonds are tracked.
The intended final state has HG1 on terminal N, an oxyanion, and restored
Asp272 carboxylate. These assignments require final coordinate inspection.

The 103251-atom solvated system, 98 original QM atoms plus eight links,
charge −2/singlet, DFTB3/3ob-3-1 Hamiltonian and original MM restraint
reference/constant are retained. No water is added and no classical TI1
parameters are introduced. Start from the original ES endpoint, not the
failed direct-N forward scan or the stalled strict-ES endpoint.

Use native Sander `ntmin=5` with DL-FIND Cartesian L-BFGS and analytic
gradients. The maximum gradient-component tolerance is 1e-4 Hartree/Bohr
(approximately 0.119 kcal/mol/Å), together with the native RMS-gradient,
step and energy-change checks. This replaces the inadequate global-RMS-only
acceptance and avoids the prior conjugate-gradient line search. All atoms
remain active. The native optimizer is being tested here; convergence is
not guaranteed. A failed optimizer stops the sequence without a fallback.

Ordinary harmonic NMR restraints use rk2=rk3=100 kcal/mol/Å² in Amber's
kΔ² convention. Each window must report native convergence, finite coordinates,
the unchanged box and both achieved coordinates within 0.15 Å of target.
Chemical gates check the covalent bond, C–O/C–N geometry and proton identity
before advancing between stages. No reverse scan, MD or TS search is queued.

Keep ordinary DISANG separate from the optional `&dlfind disang` hard-
constraint converter; that converter is not used. NMR reads its weight and
redirection records first, so `&dlfind` follows those records in the input.
Attempt 1 failed before gradients because of the original namelist ordering
and remains preserved. Attempt 2 reached L-BFGS energy/gradient evaluations.

ntwx=5 records DL-FIND evaluation callbacks, including trial evaluations;
these are not necessarily accepted optimizer iterations or physical times.
Completed windows also export the exact endpoint and a stripped DCD with
the fixed protein/substrate/three-water selection. Original full-system
trajectories remain unchanged. The construction scan is not itself a
validated reaction barrier. Keep restraint energies explicit in later plots.

Specification: configs/asp272_relay_pilot_20261025.json. Source, input hashes,
environments, logs and process/scientific records are retained per attempt.
Ruff passed and the first real native optimizer evaluations were checked;
test coverage was not measured for this one-off experimental controller.

## Asp274 and runtime assessment (2026-10-25)

In the original ES, Asp274 OD2 is 3.535 Å from terminal N; its H2 contact
is long (H···O 2.752 Å, angle 133.9°). Asp274 is 5.493 Å from OG1.
Its closest carboxylate O–O separation from Asp272 is 5.854 Å. Thus this
starting conformation does not show a direct Asp–Asp hydrogen bond or a
well-positioned Asp274 acceptor for HG1. A working hypothesis is electrostatic
and conformational modulation of the terminal amine/active-site network;
a specific proton-transfer role or pKa shift has not been established.
The NylC substrate-hydrolysis proposal assigns an Asp308–Asp306–Thr267 network,
including terminal N (paper numbering; Nyl12 sequence mapping still required):
https://febs.onlinelibrary.wiley.com/doi/10.1111/febs.16755 .

At 17:03 EDT, the first native window was around 148 energy evaluations after
approximately five minutes, still unconverged. A planning estimate of 2–5 h
for phase 1 assumes roughly 200–500 evaluations per window at about 2 s per
evaluation. This is provisional, not a completion prediction: convergence and
target tracking remain untested at a completed window, and either can stop
the run. Subsequent warm starts may reduce the time; difficult windows may
increase it. Do not present this estimate as measured full-phase throughput.

## Energy convention and VMD progress snapshot (2026-10-25)

The optimizer minimizes E_QMMM + V_position + V_drive. Preserve and report
that total for diagnostics. The underlying restrained-environment potential at
the same coordinates is E_total - V_drive; subtracting the artificial driving
spring does not remove physical bond distortion, electrostatics, or protein
strain. Remaining MM positional restraint energy is still included and must be
identified. Neither curve is currently a converged barrier or free energy.

An explicit primary example of restrained scans subtracting the spring energy
is Likhitsup et al., *Apparent non-statistical binding in a ditopic receptor for
guanosine* (2009), DOI https://doi.org/10.1039/b812969j, computational methods
(author manuscript https://wrap.warwick.ac.uk/id/eprint/27963/23/WRAP_Marsh_OBCb812969J_postprint.pdf).
This is a torsional scan, not an enzymatic QM/MM study. Developer documentation
also explicitly separates bias during minimization from reported scan energy:
https://t-0hmura.github.io/pdb2reaction/scan.html (accessed 2026-10-25).
Thus keep both curves, labeled, rather than interpret the artificial spring
penalty as part of a physical reaction barrier.

A reproducible tmux export now provides 239 actual optimization evaluation
frames across the first 10 completed Asp272-relay windows (stage1_00–09),
19976 atoms, 57.3 MB. Every exact endpoint is retained; evaluation frames are
sampled with a declared stride and include optimizer trial points, not physical
time. No interpolation or alignment. Source script/hashes/frame mappings:
provenance/asp272_relay_movie_20261025/build_movie.py and
output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/vmd_progress_snapshot_20261025/.
Curated laptop files: output_vmd_inspection_20261024/asp272_relay_progress.vmd,
matching .dcd and asp272_relay_frames.json. Existing no-argument rsync script
includes these files; checksum verification passed for every package entry.
QM Licorice, orange non-QM substrate carbons, cartoon off by default.
VMD itself is absent on violet, so no graphical VMD test is claimed.

## Tracking stop, continuation policy and automatic inspection

Window stage1_10 converged natively and retained valid coordinates/box, but
OG1–C=2.300403 Å missed its 2.146639 Å harmonic target by 0.153764 Å.
The proton difference was 0.007425 Å vs target 0.125425 Å (miss 0.117999 Å).
Thus an arbitrary 0.150 Å tracking gate stopped a converged construction window;
this is not evidence of an inaccessible TS or failed mechanism. User authorized
continuation with a 0.150 Å warning and 0.200 Å hard tracking limit, retaining
all optimizer convergence and stage-end chemical gates. No force constant,
Hamiltonian, reference, target schedule or existing geometry was changed.

User explicitly requests staying in the same run directory. Continue there
from stage1_10/optimized.ncrst, skipping all existing windows and proceeding to
stage1_11. Preserve status_before_continuation_20261025.json, the original
acceptance records and inputs.json. continuation_tracking_policy.json records
window10 reassessment, policy and source hashes; its original 0.15 Å acceptance
remains false as historical evidence. New records use the declared 0.20 Å field.
Controller: provenance/asp272_relay_20261025/run_relay_continuation.py.

Every future reaction-path controller must either call publish_inspection after
each finalized window and at termination, as this continuation does, or use
scripts/run_nyl12_path_inspection.py --watch to supervise its prepared run.
The latter accepts a controller command after --. Use the tmux job launcher
for either approach. Inspection errors fail visibly; do not silently skip them.
Do not modify frozen running controllers just to add inspection.

Stable workstation viewer/plot paths are:
output_vmd_inspection_20261024/reaction_path.vmd,
reaction_path.dcd, reaction_distances.png, reaction_distances.csv and
reaction_distances.json. Re-run the ordinary no-argument laptop pull and reload
the viewer to see new frames. A loaded VMD molecule does not automatically
reread an updated DCD. The script rebuilds the small endpoint DCD with one exact
frame per finalized window, publishes files atomically and updates package
checksums last. Detailed evaluation-frame movies remain available separately.
No coordinate interpolation, fitting or alignment is applied.

Sixteen explicitly mapped distances include the forming OG1–C bond, substrate
C–O and C–N, original HG1 distances to OG1/Asp272 OD1/OD2/Thr N, original
N–H bonds, donor–acceptor distances and contacts to both nearby Asp residues.
Explicit atom names/residue IDs are checked against topology. CSV/JSON and
plots distinguish rejected endpoints, and stage boundaries are marked.
A configurable 0.25 Å adjacent-window change is a review flag, not a scientific
acceptance limit or proof of discontinuity. Smooth endpoints do not prove
there is no hidden conformational jump. At window10 no monitored jump exceeds
0.25 Å. Examine geometry and consider finer windows/alternative coordinates
before adding constraints; never add them automatically.

Immutable internal inspection updates retain input hashes and frame provenance;
latest selects the newest successful update. The scientific run stays in the
same directory. Existing ff19SB/GAFF1/AM1-BCC/TIP3P model is unchanged.

### Expert enzyme QM/MM references replacing the earlier scan example

- van der Kamp, Chaudret and Mulholland (2013), *QM/MM modelling of ketosteroid
  isomerase reactivity indicates that active site closure is integral to
  catalysis*, FEBS J 280, 3120–3131. DOI https://doi.org/10.1111/febs.12158.
  Restrained adiabatic profiles and higher-level single points; continuity is
  investigated by forward/backward scans. This does not provide an explicit
  sentence about subtracting the driving spring.
- Lonsdale et al. (2014), *A Multiscale Approach to Modelling Drug Metabolism
  by Membrane-Bound Cytochrome P450 Enzymes*, PLoS Comput Biol 10:e1003714.
  DOI https://doi.org/10.1371/journal.pcbi.1003714. Enzyme QM/MM minimizations
  use harmonic coordinate restraints (1000 kcal/mol/Å², 0.1 Å increments).
  Again, do not claim an explicit subtraction instruction from these methods.
- Lonsdale, Harvey and Mulholland (2012), *A practical guide to modelling
  enzyme-catalysed reactions*, Chem Soc Rev, DOI https://doi.org/10.1039/C2CS15297E.
  Defines adiabatic mapping and cautions that local-minimum changes can cause
  discontinuities. This is a methodological review, not primary research.
- Explicit implementation evidence: CHARMM documentation, “Create a reaction
  coordinate for QM/MM”: minimize with RESDistance, write the minimized frame,
  then SKIP RESD and ENERGY to recompute without the driving restraint before
  writing the energy profile. https://charmm-gui.org/charmmdoc/cons.html
  (lines 558–587, accessed 2026-10-25).

Conclusion: report optimizer total and artificial-bias-removed potential
separately. The latter still contains actual molecular distortion and, in our
model, declared MM positional restraints. Neither is a free-energy barrier.

## Human-readable QM/MM layout

User-facing current results now live at `qmmm/pa6_asp272/` on violet and at
`nyl12-inspection-20261021/qmmm/pa6_asp272/` after the updated laptop pull.
Open plots/distances.png, vmd/view_path.vmd or scripts/run_path.py. The script
is the exact executed continuation source; initial_driver.py is historical.

The two elementary-step hypotheses are numbered01_og1_to_asp272 and
02_asp272_to_nterm. Within each, window_00, window_01, ... preserve scan order;
the original stage1_XX/stage2_XX names remain in the unchanged computational
run. window_order.json maps every movie frame and readable directory to its
original window. Each readable window includes inputs, targets and the original
geometry/acceptance record; large raw engine outputs and full solvent restarts
are not duplicated. Each window is a scan point, not a separate elementary
chemical reaction.

The curated package is generated by scripts/publish_nyl12_qmmm.py using
configs/qmmm_inspection_layout_20261025.json. A persistent tmux publisher
follows finalized inspection updates and refreshes this layout. Calculation
files and the active controller are not moved or modified. The root qmmm guide,
builder and configuration are versioned; generated payloads are excluded from
Git. The updated pull script transfers this curated tree and verifies checksums.

Verification at this checkpoint: final full `uv run pytest -q` in tmux passed
681 tests with13 skips. Focused helper coverage: distances100%, path inspection90%,
download publisher100%, structured QM/MM package98%. Ruff, shell syntax and
package SHA-256 checks passed. The first resumed window(stage1_11) converged
with a0.170179 Å tracking miss, correctly warning rather than stopping; stage1_12
started. The movie and structured package automatically increased from11 to12
frames. VMD GUI validation remains for the laptop; VMD is absent on violet.

## Reactive distance labels in the VMD viewer

The preferred `qmmm/pa6_asp272/vmd/view_path.vmd` viewer now creates native VMD
bond-distance labels in Å, checking atom index/name identity before labeling.
They update with animation and select the appropriate three bonds per stage:

- Stage1: OG1–HG1, HG1–Asp272 OD1, OG1–substrate C54.
- Stage2: Asp272 OD1–HG1, HG1–terminal N, OG1–substrate C54.

The actual stripped zero-based atom indices are declared in
configs/qmmm_inspection_layout_20261025.json; they are not inferred by nearest
atom or display position. scripts/nyl12_reactive_distances.tcl installs native
labels and switches visibility on frame changes. Reversed-order bonds share one
label, and unrelated labels are preserved. Existing Licorice and carbon colors
are unchanged. Re-run the usual pull and reload view_path.vmd to use the labels.
Only the inspection publisher was restarted as nyl12_qmmm_labels_20261025;
the scientific calculation and frozen controller were not changed.

API reference: https://www.ks.uiuc.edu/Research/vmd/current/ug/node133.html
(accessed2026-10-25). VMD is absent on violet; Tcl tests mock the VMD commands
to verify atom checks, stage switching and label ownership, rather than claim
a graphical VMD test.

## Energy plot in the preferred inspection layout

`qmmm/pa6_asp272/plots/energies.png` now accompanies distances.png and updates
from the same finalized endpoint sequence. Numeric values and native source
hashes are in data/energies.csv/.json; scripts/plot_reaction_energies.py and
scripts/path_energies.py expose the implementation. Both total energy and
explicit driving-bias-removed potential energy are shown, relative to the first
window in kcal/mol. The implemented extraction explicitly checks the existing
100 kcal/mol/Å² springs, targets and harmonic-region bounds; it rejects changed
restraint parameters rather than silently applying the old formula. The MM
positional restraints remain included in both curves.

The current plot has14 endpoints: windows00–12 accepted by their documented
tracking policy and converged window13 rejected by the0.20 Å gate. Native
DL-FIND Etot is used, not the coarse legacy NSTEP value. This rising incomplete
scan is not a validated barrier. Only derived inspection publication changed.

## Warning-only continuation and mechanistic diagnosis

User questioned repeated tracking-limit stops and requested progress. Window13
converged but missed its target by0.204562 Å, barely exceeding the previous
0.200 Å gate. The current in-place continuation treats tracking as a diagnostic
warning above0.15 Å with no tracking hard stop. Native optimizer convergence,
finite coordinates/fixed box and intended chemical endpoint gates remain.
No spring strength, targets, Hamiltonian or coordinate reference changed.
Controller: provenance/asp272_relay_20261025/run_relay_warning_only.py, frozen
in the run/source directory. Session: nyl12_phase1_warning_only_20261025.
It resumes window13 without rerunning it, then runs windows14–16 of phase1.
It holds phase2 for donor/amine orientation review, even if phase1 chemistry
passes, rather than treating poor proton-transfer alignment as resolved.
Original acceptance/status/policy records are preserved alongside the active
warning_only_tracking_policy.json. Three targeted continuation tests pass;
full-suite check:701 passed,13 skipped. Orchestration coverage is not quantified.

A scripted exact-coordinate review at window13 confirms N–H1···Asp272 OD1
H···O2.150 Å /149.6°, N–H2···Asp274 OD2 H···O2.278 Å /147.5°.
The donor angle Asp272 OD1–HG1···N is95.6°, N···Asp272 OD1 is3.064 Å.
Thus both amine and protonated Asp-OH orientation need attention before
assuming the second transfer will proceed. The terminal N remains NH2.
The negative sum of N–CA/N–H unit vectors was also recorded as a geometrical
lone-pair proxy, explicitly not an electronic-orbital calculation.
Script: analyze_amine_orientation.py; output: data/amine_orientation.json.
Hydrogen bonds can favor the current local minimum; local minimization alone
does not establish whether an alternative rotated/inverted geometry is accessible.

### Why the energies rise

An exact printed-energy comparison of window00 to window13 yields(kcal/mol):
- Total:+49.670900.
- Artificial driving bias:+4.455684.
- MM positional restraints:+2.030016.
- Printed DFTBESCF term:+46.532200.
- Remaining printed terms:−3.347000.

This partition rules out the applied springs as the main accounting cause of
the rise. DFTBESCF is a printed QM term, not an Asp-specific or hydrogen-bond
energy, and cannot be assigned to a particular interaction or pKa shift.
N–H donation can stabilize carboxylate relative to neutral carboxylic acid,
making proton acceptance less favorable; the user's proposal is plausible.
A primary example of an enzyme hydrogen-bond network lowering an Asp pKa is
*Direct evidence that an extended hydrogen-bonding network influences activation
of pyridoxal 5′-phosphate in aspartate aminotransferase*(2017):
https://pmc.ncbi.nlm.nih.gov/articles/PMC5392587/.
That is supporting chemical precedent, not a measured Nyl12 pKa or proof of
this pathway. The opposing effects of other charges, solvation, geometry and
method error have not been separated.

At window13, Asp272–HG1=1.094 Å and OG1–HG1=1.504 Å, while OG1–C=2.053 Å
and substrate C–O=1.279 Å. Proton transfer has progressed strongly before
formation of the tetrahedral linkage. A costly alkoxide-like configuration
along an unfavorable imposed route is therefore another supported hypothesis.
The current linear target schedule need not follow the lowest-energy valley.
This is an inference from the geometry, not proof of the correct mechanism.

Next scientific decision: review the bounded phase1 endpoint, then compare
attack/proton-transfer coupling with a small local two-coordinate map and
prepare amine/Asp-OH orientations before phase2. Do not assign a real barrier,
remove an Asp hydrogen bond artificially, change protonation silently or claim
that increasing a tracking allowance fixes the energetic problem. Relative
DFTB proton affinities may also require matched higher-level spot checks before
relying on a barrier; none are launched by this diagnostic.
The analysis is scripted in audit_relay_energy_components.py; numeric results
and native source hashes are in data/relay_energy_component_audit.json.


### Asynchrony and DFTB3 proton-transfer validation — 2026-10-25

The user's correction is accepted: asynchronous proton transfer and attack are
not intrinsically problematic. Our two harmonic targets allow different actual
coordinates to advance at different rates; their linear schedule still samples
one imposed route, rather than finding the lowest-energy route automatically.
OG1–C has shortened from about 3.14 to 2.05 Å at window13, so approach has begun,
but no tetrahedral minimum has been demonstrated. Deprotonation increases
nucleophilicity without making attack necessarily downhill or barrierless.
A high intermediate energy could reflect relative proton affinities, attack
geometry, insufficient oxyanion stabilization, restrained environmental response,
or an unsuitable route. The scan alone does not identify the biochemical cause.

Primary literature reviewed:
- Arantes and Rezac (2025), DOI 10.1021/acs.jctc.5c00690, uses our exact
  DFTB3/3ob-3-1 set. Mean unsigned relative-energy errors are 3.63 kcal/mol
  overall, 1.37 for carboxylic-acid/carboxylate and 3.44 for amine profiles.
  These isolated proton-transfer benchmarks do not directly test coupled
  aliphatic-alcohol deprotonation and amide attack in Nyl12.
- Kromann et al. (2016), DOI 10.7717/peerj.1994, reports DFTB3 barrier and
  reaction-energy mean absolute deviations of 5.8 and 6.3 kcal/mol across five
  enzyme model reactions against its DFT reference. Individual errors can be
  larger and are chemistry dependent; this is not a universal error bound.
- Berselli et al. (2025), DOI 10.1021/acscatal.5c01263, obtains a PET4 acylation
  barrier around 21 kcal/mol with DFTB3, versus about 24 after B3LYP/RI-MP2
  single-point free-energy reweighting. Ester hydrolysis is useful precedent
  for reaction-specific validation, not validation of nylon amide hydrolysis.

These results justify checking DFTB3, but do not establish that our roughly
50 kcal/mol rise is a method artifact. Gas-phase proton affinities, potential
energy profiles and enzymatic pKas are different observables; pKas require
free-energy and reference-state treatment.

Recommended next diagnostic, NOT executed: matched higher-level single points
on ES, proton-shared, Asp-protonated and advanced-attack geometries, using the
same QM chemistry and protein electrostatic embedding. Compare relative method
corrections, not absolute energies between Hamiltonians. Use the agreed
omegaB97X-D without extra dispersion and include a diffuse-basis sensitivity
check for the anionic region. A small unconstrained gas cluster would not be an
equivalent test. Continue the already bounded phase1 run, but defer further
path variants until this diagnostic and endpoint/orientation review are planned.

### Phase1 completed; numerical completion versus TI1 — 2026-10-25

#### Two additional phase-I windows — latest user instruction

The user requested extending the original phase-I path, then explicitly kept
the original intervals and reduced the request to a few steps. The authorized
extension is exactly stage1_17 and stage1_18 from the original stage1_16
restart. Neither the previously discussed phase-II continuation nor ten finer
windows was prepared or launched. There is no amine rotation or change to the
reaction-coordinate definition.

Continue the original linear targets with global fractions 17/16 and 18/16:
attack target = a0 + f*(1.55-a0), and proton target = q0 + f*(0.8-q0),
where a0 and q0 are the original initial target values. The attack targets are
approximately 1.450560 and 1.351120 Å; the OG1–HG1 minus Asp272–HG1 targets
are approximately 0.912429 and 1.024858 Å. Retain the two k100 harmonic
springs, original MM reference/k10 positional restraints, DFTB3/QM region,
native optimizer tolerances and warning-only tracking policy. These are soft
target distances, not fixed atomic distances. Stop after the two windows for
energy/geometry review; do not extrapolate another eight original increments
toward an unphysical sub-angstrom attack target.

Numerical convergence and energy turnover remain separate observations. Both
total and bias-subtracted energy curves, all reaction distances and the endpoint
VMD movie are republished after each window. A downturn in a driven profile is
not alone proof of an unrestrained intermediate minimum. The unrestrained
diagnostic is a separate branch and is not the continuation's starting model.

#### User-requested unrestrained endpoint diagnostic — launched

The original stage1_16 endpoint is retained without amine rotation. A separate
child calculation `phase1_unrestrained/` starts from a byte-identical copied
restart and frozen topology. Both harmonic driving terms and all MM positional
restraints are removed (`nmropt=0`, `ntr=0`, no DISANG or reference argument).
All 103251 atoms are active, with no SHAKE (`ntc=1`, `ntf=1`). The periodic
cell remains fixed; QM98/link8/charge−2, DFTB3, embedding and the inherited
ff19SB/GAFF1/AM1-BCC/TIP3P Hamiltonian are unchanged. This is a geometry
optimization, not MD or an amine-orientation intervention.

Runner: `provenance/phase1_unrestrained_20261025/run_optimization.py`; the
executed frozen copy and exact command/environment/input hashes accompany
the child run. Session: `nyl12_phase1_unrestrained_20261025`; job metadata:
`output_nyl12_pa6_dftb_pilot_20261024/jobs/phase1_unrestrained_20261025/`.
Native LBFGS Cartesian criteria remain tol1e−4 Hartree/Bohr, tole1e−7;
maxcycle is 2000. `ntwx=5` saves force-evaluation snapshots, including trial
coordinates, rather than physical MD time. Engine convergence and endpoint
chemistry are recorded separately. Seven focused input-transformation checks
passed in tmux; the launch log already confirms RESTRAINT=0.

Inspect this result before deciding whether to sample a different terminal
amine/Asp-OH orientation. No manual rotation or further driving was launched
with this diagnostic. The separate stripped trajectory will be published below
`qmmm/pa6_asp272/phase1_unrestrained/`, preserving the original 17-window path.

An intermediate restart inspected at 2026-10-26 00:42:50 UTC has OG1–C
3.1500 Å, OG1–HG1 0.9921 Å, Asp272–HG1 1.8245 Å, N–HG1 2.5896 Å,
substrate C–O 1.2650 Å and C–N 1.3561 Å. Thus this checkpoint has returned
toward ES, with HG1 back on OG1. This is not the final optimized endpoint or
a general conclusion about native enzyme chemistry. Scripted measurements and
restart SHA-256 are in the child run's `checkpoint_geometry.json`.

The VMD renderer uses pair-specific DynamicBonds with distance cutoffs 1.25 Å
for the candidate proton bonds and 1.9 Å for OG1–C54. These are display cutoffs,
not computed bond orders. Original candidate topology edges are removed only
inside the VMD display molecule to prevent permanent ghost sticks; source
topologies are unchanged. HG1 is shown as a small VDW sphere, with persistent
O49–Phe113 backbone N and O49–Asn186 ND2 distance labels. Non-QM substrate
carbons remain orange and the protein NewCartoon remains available but off.

All 17 stage1 windows converged with the unchanged native optimizer tolerance.
Window16: OG1–C 1.795256 Å; OG1–HG1 1.793376 Å; Asp272 OD1–HG1
1.021580 Å; terminal N–HG1 2.731658 Å; carbonyl C–O 1.302931 Å;
substrate C–N 1.427586 Å. This passes the controller's stage1 geometric gate
but is not the intended NH3+ TI1 or a demonstrated unbiased minimum.
The per-window `chemical_state_validated:false` is intentional: convergence
and the stage gate do not certify a stationary intermediate. Total relative
energy is +66.1096 kcal/mol, including 6.0946 kcal/mol of driving bias;
bias-removed rise is +60.0173 kcal/mol, with MM position restraints retained.
These endpoint energies are not a measured activation barrier.

The user emphasizes that a numerically converged driven path remains useful
even if it is biochemically incorrect. Agreed: method accuracy must not be used
to excuse numerical failure or arbitrary tracking cutoffs. A sequence of
converged biased optimizations is attainable without asserting that the
requested unbiased intermediate exists under every Hamiltonian. Smoothness,
window convergence, endpoint identity, stationary-point character and native
mechanistic relevance remain separate acceptance records.

Next path work: review donor/amine orientation at window16; continue the
Asp–HG1 to terminal-N transfer while retaining the forming OG1–C contact and
allowing physical local reorientation. Add a declared orientation coordinate
only if needed, preserving the unmodified reference. Complete the forward
construction, reduce/remove driving bias to test TI1 stability, and reverse
the route from an accepted endpoint. Do not wait for full method benchmarking
to finish numerical path construction. Use continuous distance/energy plots
and VMD updates. Refine candidate saddle regions with DL-FIND P-RFO, then
check the unstable mode and downhill connections; a driven maximum is a seed,
not a TS. Multiple genuine minima imply multiple connecting TSs.

### Tight-binding alternatives and focused benchmark — planned, not executed

R35 isolated proton-transfer mean unsigned errors (kcal/mol):

| Method | All reactions | Carboxylate | Amine |
|---|---:|---:|---:|
| DFTB3/3ob-3-1 | 3.63 | 1.37 | 3.44 |
| GFN2-xTB | 3.23 | 2.39 | 5.31 |

GFN2-xTB is therefore a useful comparator, not a guaranteed upgrade for our
alcohol/carboxylate/amine chemistry. R39 also finds it competitive with DFTB3
for a ribozyme example, with both differing quantitatively from PBE0.

R38 supporting Table3 reports full-GMTKN55 WTMAD-2 values of 24.612
(DFTB3-D3(BJ)/3ob-3-1), 25.031 (spin-polarized GFN2-xTB), 9.254 (g-xTB),
7.366 (r2SCAN-3c), 8.032 (B3LYP-D4/def2-TZVP), and 3.246
(omegaB97M-V/aTZ), in kcal/mol. These weighted broad-set statistics are not
mean enzyme-barrier errors, not all necessarily identical supported subsets,
and not results for our exact DFTB setup or Gaussian omegaB97X-D. They refer
to the 2025 preprint model, not the subsequently changed g-xTB development
binary. No reaction-specific g-xTB advantage is established here.

Compatibility observed 2026-10-25:
- sander has a published direct GFN2-xTB library interface with electrostatic
  embedding/QM-MM PME [R39]. Local AmberTools26 `xtb_module.F90` implements
  GFN1/GFN2, but retained sander compile flags lack `-DXTB`; its disabled
  branch explicitly raises `XTB is not enabled`. Do not claim the current
  executable is ready for xTB or rebuild it silently.
- Current Grimme g-xTB README documents analytic gradients, but explicitly
  lists point-charge embedding as unsupported. It is not a drop-in sander
  Hamiltonian. ASH can drive g-xTB cluster calculations and also documents
  the lack of electrostatically embedded g-xTB QM/MM. Switching drivers does
  not repair that backend limitation.
- If a new driver becomes necessary, ASH/OpenMM with an appropriate QM backend,
  or ChemShell with ORCA and DL-FIND, merits evaluation. ORCA itself supports
  GFN2-xTB and QM/MM. Preserve boundary, force-field, electrostatic and periodic
  conventions explicitly; these are alternatives, not installed/validated
  replacements. Retaining sander is the least disruptive immediate choice.
- No MLP fitting is planned. Pretrained corrections would still need charge,
  element, reactive-chemistry and embedding validation; availability alone is
  not a reason to substitute one.

Small-model benchmark specification: use 2-propanol/acetate (Thr-to-Asp,
charge -1), acetic acid/methylamine (Asp-to-terminal-N, charge 0), and
2-propanol/methylamine (direct alternative, charge 0), all closed-shell.
Construct matched constrained proton-transfer profiles at declared donor-
acceptor separations. Preserve unstable ion-pair branches rather than allowing
unconstrained optimization to erase them. Compare DFTB3 with omegaB97X-D on
identical coordinates; separate geometry sensitivity from energy differences.
Use the agreed def2-SVP for economical preparation and larger/diffuse-basis
single-point checks for anions. Neither this DFT reference nor a gas-phase
cluster alone is an experimental pKa standard. Include a small amide-attack
test and a few matched embedded pathway points before transferring conclusions
to the enzyme. GFN2 is optional after explicit software setup approval;
g-xTB can first be evaluated only for clusters if installation is approved.
Report signed profile errors, reaction energies, barriers when present,
geometry/proton-state retention, convergence and model/version hashes.
# Corrected unrestrained construction protocol (2026-10-25)

The initial relay pilot used a 10 kcal mol-1 A-2 positional restraint on every
MM atom. That suppresses protein, substrate and solvent reorganization and is
not suitable as the relaxed E-S reference or the production construction path.
Those calculations are retained only as restrained diagnostics; they cannot
establish a mechanistic barrier.

The corrected workflow first converges the complete solvated E-S complex with
all atoms active and no positional or reaction-coordinate restraints. It then
uses 21 bounded windows for OG1 attack/HG1 transfer to Asp272 OD1 and 21 windows
for Asp272-to-terminal-N proton transfer. During the second stage, OG1-C is
biased at 1.55 A solely to construct a TI1 seed; this is explicitly artificial
and is not interpreted as an unbiased path. All biases are subsequently
released and TI1 must remain tetrahedral in a converged optimization before it
is accepted as a minimum candidate. Highest-energy driven structures will only
serve as starting guesses for later first-order saddle searches.

