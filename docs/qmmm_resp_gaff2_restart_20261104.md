# Nyl12 PA6 RESP/GAFF2 QM/MM restart checkpoint

Status: active, 2026-11-04. This checkpoint replaces the old AM1-BCC/GAFF
E--S seed as the production starting point; it does not certify a QM/MM result.

## Starting ensemble and selection

The noncovalent substrate is the released modular zwitterionic PA6 L4 model
`P6N-P6I-P6I-P6C` with DFT/RESP charges and GAFF2 bonded/Lennard-Jones terms.
The exact solvated topology SHA-256 is
`01e06919f919b0905a6fbe794ad5cf24463ff0b80c30593d0dc8fca6ebfd104f`.
WT replicas 1--3 completed 100 ns and passed the declared persistence rule.
Replica 3's full trajectory currently resides on azalea, so the first violet
pilot uses independent replicas 1 and 2; replica 3 remains a follow-up rather
than being silently omitted from the ensemble claim.

Snapshot selection does not rank eligible frames by catalytic distances. Three
nonoverlapping windows span 50--100 ns. Within each window, Tier A is preferred
to Tier B, the held-out Thr233 torsions must satisfy absolute native differences
of at most 35 degrees for psi and 20 degrees for chi1, and the eligible frame
nearest the window midpoint is selected. The all-frame screens found 68 and 135
joint Tier A/B plus native-torsion frames in replicas 1 and 2. Exact inputs,
hashes, extracted restarts and native overlays are under
`qmmm/pa6_resp_gaff2_wt/00_snapshot_selection/` and
`qmmm/pa6_resp_gaff2_wt/01_native_thr_snapshot_selection/`.

The first promoted pilot is replica 1 frame 1553 (77.60 ns). Before
minimization it is Tier A: d1/d2/d3 = 3.061/3.099/3.293 A, Burgi--Dunitz angle
= 91.72 degrees, and projected Herschlag phi = -94.88 degrees. Thr233 psi and
chi1 differ from `cif/n12c.cif` by -8.35 and +10.41 degrees. Its aligned Thr N
displacement is 0.845 A; OG1 displacement is 0.850 A, 0.100 A outside the
conservative 0.75 A gate. It is a near-pass seed, not an accepted endpoint.

## Conditional MM to QM/MM queue

The MM stage has `ntr=0`, `nmropt=0`, all atoms active, `maxcyc=10000`, and
`drms=0.01`. It is queued on violet GPU 0 and will not preempt the active PA6
production run. DFTB3/MM launches only if the MM endpoint converges, remains
Tier A/B, and passes the native Thr233 displacement/torsion gate.

Two conditional DFTB3/MM branches compare 8 and 12 A movable shells centered
on Phe113, Asn186, Thr233/234, Asp272/274 and PA6 residues 1313--1316. The
catalytic groups, substrate and nearby protein/water environment are active;
distant bulk atoms are fixed by DL-FIND's active mask. There are no harmonic
positional or reaction-coordinate restraints. This tests boundary sensitivity
without repeating the rejected all-MM positional-restraint protocol or asking
the optimizer to remove gradients from every bulk-water coordinate.

Active sessions are `nyl12_wt_r1f1553_mm_queue_20261104b`,
`nyl12_wt_r1f1553_qmmm8_queue_20261104`, and
`nyl12_wt_r1f1553_qmmm12_queue_20261104`.

## QM subsystem

The first comparison retains the reviewed chemistry of the prior primary
region: 98 real QM atoms plus eight link hydrogens, net charge -2, singlet. It
includes N-terminal Thr233 and its following backbone fragment, the Phe113
backbone oxyanion donor, Asn186 side chain, both Asp272 and Asp274 carboxylates,
and a 16-atom fragment spanning the scissile modular PA6 amide (P6I 1314 C/O
and P6I 1315 N). The mapping was rebuilt by residue/atom identity and audited;
old whole-ligand indices were not reused. No water was forced into QM. A water
may be added only under a declared structural/persistence rule or when testing
a specific water-mediated mechanism.

Retaining both Asp groups is conservative: previous trajectories showed strong
coupling to the N-terminal amine, so excluding either group would prejudge
proton transfer. Region-size sensitivity remains required. Andjelkovic et al.
reported that their hASNaseIII DFTB3 screening region contained at least the
whole substrate and N-terminal Thr, while larger regions adding interacting
residues did not materially change their profile. Their accepted acylation path
was stepwise: proton transfer to the substrate amino group preceded Thr attack;
DFTB3 was used for exploration and B3LYP-D3/6-31+G*/MM adaptive-string work for
the reported free-energy profile (DOI: 10.1021/acs.jcim.3c00900). This is an
analogy for substrate hydrolysis, not evidence about Nyl12 autocatalytic
cleavage and not proof that Nyl12 must use the identical base.

## Prior TI1 result and Hamiltonian gate

The previous 42-window coordinate drive produced a TI1-like seed, but a fully
unrestrained DFTB3/MM optimization returned to E--S: OG1--C changed from
1.70283 to 3.05524 A, HG1 returned to OG1 (0.99217 A), and carbonyl C--O/C--N
became 1.26529/1.35661 A. Thus TI1 was not a stable minimum under that old
Hamiltonian/seed; this does not prove that physical TI1 is absent.

DFTB3 is not assumed accurate for the decisive proton transfer. Arantes and
Rezac found chemistry-dependent errors and warned that no tested approximate
method is generally transferable without validation; DFTB3's isolated-reaction
MUE was 15.2 kJ/mol overall and 14.4 kJ/mol for amine reactions, while
omegaB97X was more accurate overall (DOI: 10.1021/acs.jctc.5c00690). The
planned gate uses minimal gas-phase alcohol/amine/carboxylate proton-transfer
models with fixed donor--acceptor distances and matched coordinates. Compare
DFTB3/3ob-3-1 against the agreed omegaB97X-D/def2-SVP reference, followed by a
diffuse/larger-basis sensitivity check for anionic states. Report signed profile
errors, proton localization, state ordering and geometry retention. A gas-phase
benchmark diagnoses the Hamiltonian; it does not replace embedded enzyme checks
or yield a pKa.

## Go/no-go decisions

1. Promote an E--S endpoint only after engine convergence, Tier A/B geometry,
   native Thr233 gates and visual inspection all pass.
2. Require the 8 and 12 A movable-shell endpoints to agree qualitatively before
   using an endpoint for path construction.
3. For the direct E-S-to-TI1 path, use the asynchronous composite coordinate
   `d(OG1-HG1) - d(N-HG1) - d(OG1-Ccarbonyl)`. Bias only this sum so proton
   transfer and attack can advance independently; coordinate driving supplies
   seeds, not a barrier or transition state.
4. If DFTB3 fails the small-model or embedded checks qualitatively, pivot the
   production path Hamiltonian rather than forcing a DFTB3 mechanism.
5. A retained TI1 minimum and reproducible forward/reverse path are required
   before DL-FIND P-RFO saddle refinement and DFT/MM validation.

## MM continuation outcome — 2026-11-05

The original 10,000-cycle unrestrained MM attempt reached its cycle limit at
RMS 0.089048 kcal mol-1 A-1. A new immutable continuation started from that
exact restart with a 20,000-cycle budget and no restraints. It converged after
9,317 additional cycles at RMS 0.009969, below the unchanged 0.01 target.

Convergence did not restore the initial Tier-A pose. The endpoint remains Tier
B: d1/d2/d3 = 3.492/2.957/2.899 A, attack angle = 92.40 degrees, projected phi
= -121.22 degrees and CB--OG1--C = 120.26 degrees. Thr233 torsions and N
displacement pass, but aligned OG1 displacement is 0.871 A versus the 0.75 A
gate. CB--OG1--C changed from 141.28 degrees initially to 120.29 degrees after
the first segment and 120.26 degrees after convergence, showing that the
reorientation stabilized rather than reversing.

No QM/MM optimization was launched from this endpoint. The retained acceptance
record is under
`qmmm/pa6_resp_gaff2_wt/02_mm_minimization/replica_01_frame_1553_continuation_20000_attempt1/assessment_attempt1/`.
Before another QM/MM attempt, screen the full eligible WT-frame set for native
Thr233 displacement rather than spending additional minimization on this
converged basin. The other six midpoint-selected snapshots are not immediate
replacements: their initial OG1 displacements are 1.865--2.753 A except for the
already tested frame 1553.

## Full WT native-Thr screen and selected seeds — 2026-11-06

The all-frame screen is complete for all three WT replicas. Replicas 1 and 2
contained 203 joint Tier A/B plus native-torsion frames; four passed the held-out
native Thr233 displacement gate, all from replica 1 between 77.65 and 80.10 ns.
Replica 3 contained 157 joint-eligible frames but no native-displacement pass.
The exact replica-3 calculation and provenance were run on azalea and copied to
`qmmm/pa6_resp_gaff2_wt/06_replica03_native_thr_screen_20261106_attempt1/`.

Three complementary seeds are promoted for unrestrained MM minimization:
replica-1 frames 1554 (77.65 ns, the only Tier-A pass), 1562 (78.05 ns, the
lowest Thr233 RMSD and nearly native chi1), and 1603 (80.10 ns, the latest pass
and least temporally redundant, although its d1 is near the Tier-B limit).
Frame 1572 remains an alternate. All promoted frames are correlated within a
2.45-ns interval, so agreement among them is not equivalent to independent
replica support. The machine-readable decision and exact source hashes are in
`qmmm/pa6_resp_gaff2_wt/07_selected_es_candidates_20261106/selection_decision.json`.
QM/MM remains conditional on converged MM endpoint Tier A/B and native-Thr
acceptance.

## Three-seed unrestrained MM outcome — 2026-11-06

Replica-1 frames 1554, 1562, and 1603 were extracted as complete solvated
restart files and minimized independently with no positional or reaction-
coordinate restraints (`ntr=0`, `nmropt=0`, all atoms active). Each used
immutable 10,000-cycle segments, a total ceiling of 30,000 cycles, and the
unchanged RMS-gradient target of 0.01 kcal mol-1 A-1. All three converged in
their second segment. The azalea frame-1603 first attempt used the generic
PMEMD executable and failed before minimization because that older binary had
no Blackwell kernel image; immutable attempt 2 used the previously validated
CUDA-12.8 Blackwell `pmemd.cuda_SPFP` and completed normally.

No endpoint passed the full scientific gate. Frame 1554 remained Tier B but
its aligned Thr233 OG1 displacement was 0.830 A. Frame 1562 lost Tier-B status
(`d1 = 3.576 A`) and had OG1 displacement 0.846 A. Frame 1603 remained Tier B
and was the nearest native-Thr pass, with OG1 displacement 0.782 A versus the
predeclared 0.750-A limit. Its d1/d2/d3 values were 3.426/3.106/2.904 A,
attack angle 89.48 degrees, and projected phi -114.71 degrees. The criterion
was not changed post hoc. Consequently, no 8- or 12-A DFTB3/MM shell was
launched. Complete centralized artifacts are under
`qmmm/pa6_resp_gaff2_wt/09_mm_qmmm_selected_candidates_20261106/`.

## Active frame-1603 near-pass QM/MM sensitivity — 2026-11-06

The user explicitly authorized QM/MM exploration from the converged frame-1603
near-pass. The standard gate remains recorded as failed. A narrowly scoped
override permits only the OG1-displacement criterion to increase from 0.750 to
0.800 A; convergence, Tier A/B geometry, N displacement, psi, and chi1 gates
remain mandatory. The observed OG1 displacement is 0.782 A. The exact decision,
source hashes, and user authorization are in
`qmmm/pa6_resp_gaff2_wt/10_frame1603_nearpass_qmmm_20261106/mm_override_authorization_status.json`.

Independent 8-A and 12-A movable-shell DFTB3/MM optimizations are active in
tmux sessions `nyl12_f1603_nearpass_qmmm08_20261106` and
`nyl12_f1603_nearpass_qmmm12_20261106`. Both `mm_acceptance.json` files record
standard acceptance as false and the explicit override as accepted. Both
Sander engines passed startup. Final scientific acceptance remains the standard
unmodified QM/MM endpoint gate; this override authorizes the starting point,
not automatic promotion of the optimized endpoint.
