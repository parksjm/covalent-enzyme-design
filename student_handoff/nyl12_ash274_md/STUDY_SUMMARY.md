# Protonation-state sensitivity of near-reactive PA6 binding in Nyl12

Working manuscript summary for a student computational study; 2026-11-19.
All numerical results below are exploratory existing results, not new student
simulations. Full references and source provenance accompany the parent
`manuscript/` drafts and the evidence JSON/CSV in this directory.

## Introduction

Nyl12 is being investigated as a promiscuous polyamide hydrolase and a scaffold
for improving PA6 conversion. Prior WT simulations suggested stable PA6 L4
binding across three independent 100 ns trajectories, whereas several single
mutants and PA66 systems were less persistent. Those observations motivate
model- and force-field-controlled comparisons rather than assuming a sequence
mutation directly caused every observed binding difference.

Related NylC work proposed a catalytic Asp–Asp–Thr network and protonated
Asp306, corresponding to Nyl12 Asp272, on structural/mechanistic grounds [1].
The assignment is circumstantial; mutations that impair maturation confound
some catalytic interpretations. Our fixed-charge MD and subsequent QM/MM
results motivate considering a protonated Asp274 as well. They do not assign
an experimental protonation state. Substrate hydrolysis must be distinguished
from precursor autocatalytic cleavage throughout this study.

The proposed student arm asks whether one fixed Asp-pair microstate increases
persistent near-reactive substrate poses and whether such pose persistence can
support a useful sequence-screening metric. It is complementary to ongoing
reaction-path calculations, which have not yet established a verified,
kinetically credible transition state or a required stable TI1 minimum.

## Existing computational methods

Six matched WT simulations used the flipped Model30 E–S starting model,
neutral mature Thr233 amino termini and three fixed microstates, each with
two independent 100 ns replicas: ASP272/ASP274, ASH272/ASP274 and ASP272/ASH274.
All four corresponding Asp copies had the declared state. Temperature was
300 K and protein/water treatment ff19SB/TIP3P [2,3]. PA6 L4 used the modular
multiconformer RESP/GAFF2 release `pa6_l4_modular_resp_gaff2_v1_20261031`,
assembled P6N–P6I–P6I–P6C, with zwitterionic termini and net substrate charge 0.

The release is accepted for this individual L4 MD application, not bulk nylon,
arbitrary lengths or catalytic adducts. It has documented validation limitations,
including a compact-conformer electrostatic holdout; use its manifest and
current acceptance statement rather than calling it universally validated.
RESP/GAFF method background is given in [4,5]; GAFF2 implementation provenance
is the exact released Amber parameter payload, not the original GAFF paper.

Production was unrestrained and classified at 10 ps cadence, 10,000 frames per
replica. Tier A required d1≤3.2 Å, d2/d3≤3.8 Å, Bürgi–Dunitz 90–120° and
70≤|projected Herschlag phi|≤110°. Positive signed phi identifies the flipped
local pose under the frozen implementation/atom convention. A separate native
Thr screen compares crystal-aligned geometry to the supplied unpublished
`n12c.cif` reference; citation information is pending from the coauthor/PI.

## Existing results

| Fixed Asp microstate | Tier A replica 1 | Tier A replica 2 | Pooled Tier A fraction |
|---|---:|---:|---:|
| ASP272/ASP274 |37/10,000|238/10,000|1.375%|
| ASH272/ASP274 |3,886/10,000|1,724/10,000|28.050%|
| ASP272/ASH274 |1,677/10,000|0/10,000|8.385%|

Thus ASH274 has an approximately 6.1-fold larger pooled fraction than ASP/ASP,
but one replica has no Tier A frames. The more complete structural screen is
less favorable: Tier A plus positive phi, intact PA6 and native-Thr gates yields
0+21 frames for ASP/ASP, 7+0 for ASH272 and 5+0 for ASH274. The wider Tier A/B
joint counts are 0+31, 12+0 and 8+0 respectively. Do not interchange these counts
or interpret the pooled Tier A increase as robust catalytic competence.

The closest Asp oxygen pairs were more compact in protonated-pair MD than in
ASP/ASP: median site-E O–O distances were approximately 2.56 Å for ASH272,
2.65–2.73 Å for ASH274 and 4.89–4.90 Å for ASP/ASP. Measure all four O–O
combinations because the two carboxylate oxygen labels are chemically
interchangeable. A short contact supports a possible hydrogen-bond arrangement,
not a measured pKa or proof of the correct proton site.

Water proximity alone also proved insufficient. A subsequent directional
screen found only 32 Tier A frames with a complete OG1→water→Thr-N arrangement
in the selected ASH274 replica, with no same-water two-leg streak longer than
three 10 ps samples. None of the three selected representative trajectories
provided a Tier A complete bridge to leaving amide N. These results neither
establish an ordered relay nor rule out rare reorganization/multiple waters.

The verified older sequence panel comprises WT, T53Q, F55N, K59D, K59E, F101E,
F101H, P105K and V111C, with 27 completed legacy runs. Their substrate treatment
was AM1-BCC/GAFF. The present inventory did not identify a completed modular-FF
mutant panel. Consequently, a new ASH274-versus-legacy comparison would confound
protonation, FF and possibly starting-model quality. Matched new controls are
required to test causation.

In QM/MM, some initially ASH272 models transferred their proton to Asp274.
This is a potential-energy behavior of those selected models, not experimental
proof of an ASH274 population. Parameter benchmarks and repeated TI1-like
releases also did not establish a stable TI1. A balanced acyl-enzyme candidate
has survived unbiased DFTB3/MM optimization, but the connecting NEB profiles
remain unconverged with high provisional maxima. None is a verified viable
activation barrier. These limitations do not prevent an explicitly exploratory
E–S MD screen, but they prevent calling it TS-guided activity optimization.

## pKa evidence and limits

Static protein-only PROPKA screens covered mature NylC 5XYG, uncleaved
mutant precursors 7YU0/7YU1 (not back-mutated), and three mature Nyl12
Tier A snapshots. Corrected terminal recognition in mature 5XYG gave
Asp306 pKa 2.19–2.21, Asp308 5.49–6.49 (default) or 7.24–8.45
(coupled alternative), and Thr267 N-terminal pKa 6.35–6.55 or 4.59–4.60.
These values supersede the earlier malformed-TER mature input.

In Nyl12's substrate-facing site, Asp272 predictions were 0.88–2.37,
Asp274 3.44–6.95 and Thr233 N-terminal 8.93–9.70 across the three
snapshots. Neither protonated Asp272 nor neutral Nyl12 Thr233 is established
by these estimates. The lack of PA6, conformation/coupling sensitivity and
mutant precursor structures limit interpretation. A possible elevated
Asp308/Asp274 pKa motivates testing ASH274; it does not demonstrate its
experimental population.

MCCE did not provide an accepted independent pKa result. Mature NylC's
extreme bounds were rejected because artificial internal termini corrupted
the input; complete Nyl12 failed terminal-recognition preflight before
titration. Do not interpret these rejected bounds as a credible disagreement
with PROPKA. Per-chain tables and rejection receipts are in `evidence/pka/`;
full Methods/Results tables are in the current repository manuscript.
The first handoff archives predate this documentation addition; obtain the
current repository documents/evidence alongside those frozen data archives.

## Interpretation and proposed next study

The defensible hypothesis is that local protonation influences the accessible
binding ensemble. The student should first test matched WT ensembles with at
least three, preferably four, independent 100 ns replicas per condition. Next,
rebuild the prior mutants with the same current FF and paired protonation
protocol. The appropriate primary uncertainty is among independent replicas;
10,000 correlated frames are not 10,000 independent experimental observations.

Report native-Thr orientation, binding persistence,signed phi and access alongside
Tier A occupancy. Preserve site-AE reproduction and add any-of-four-site analysis
before ranking designs. If a reproducible structural benefit survives those
gates, focused fixed-length LigandMPNN sampling can explore coordinated
substitutions while preserving catalytic, maturation and assembly constraints [6].
Designed sequences still require independent structural and dynamic validation;
sequence likelihood, fold confidence, binding affinity and catalytic efficiency
are distinct quantities.

## References

1. Negoro et al. FEBS Journal 2023. DOI: 10.1111/febs.16755.
2. Tian et al. ff19SB. JCTC 2020. DOI: 10.1021/acs.jctc.9b00591.
3. Jorgensen et al. Water-model comparison. JCP 1983. DOI: 10.1063/1.445869.
4. Bayly et al. RESP. J Phys Chem 1993. DOI: 10.1021/j100142a004.
5. Wang et al. GAFF. J Comput Chem 2004. DOI: 10.1002/jcc.20035.
6. Dauparas et al. LigandMPNN. Nature Methods 2025. DOI: 10.1038/s41592-025-02626-1.

The complete working bibliography, caveats and source ledger are in the parent
manuscript directory. Parameter-specific validation belongs to the separate
polyamide force-field manuscript; this summary does not duplicate its claims.
