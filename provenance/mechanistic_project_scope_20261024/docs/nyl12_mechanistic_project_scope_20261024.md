# Nyl12 mechanistic study: proposed scope and decision gates — 2026-10-24

## User constraints and existing evidence

The user has more than 20 years of computational experience and emphasizes the
nontrivial, error-prone nature of QM/MM. The proposal is a potentially separate
mechanistic paper, not an incidental extension of a redesign screen. Target
journals include ACS Catalysis, JCIM or J Phys Chem B; publication suitability
will depend on mechanistic insight, validation and evidence, not the journal
name or simply obtaining a reaction path.

Focus PA6, optionally poly(MAXA); do not expand PA66 now. NylB is mentioned as
a potentially simpler alternative, NOT authorization to pivot or start a second
enzyme campaign. Use a validated DFTB-level QM/MM Hamiltonian for free-energy
sampling; umbrella sampling or a suitable related approach is a candidate.
Hybrid DFT is restricted to selected PES optimizations/checks, not free-energy
simulations. No new method, dependency or simulation was started by this note.

User reports existing antechamber/AM1-BCC noncovalent substrate parameters and
three independent 100 ns WT PA6 L4 trajectories with stable binding; PA66 performed
less well. The October 9 synthesis independently supports persistence of the
selected PA6 pose in 3/3 WT replicas. It also says corrected PA66 can remain
bound globally but loses near-productive geometry/register. Do not rewrite this
as an affinity measurement or universal PA66 dissociation. Source:
/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/NYL12_WT_MUTANT_PROTOCOL_NARRATIVE_AND_AF3_RESTRAINT_TEST_20261009.md
and the September 28 final-interpretation note in that directory. Old PA66
chemistry-invalidated trajectories remain excluded. Exact topology/charge/
protonation and valid trajectory paths still need pinning before reuse.

## Recommendation: two linked scopes

A full chemical-cycle QM/MM free-energy study is valuable as a standalone
mechanistic project but is not a universal prerequisite for exploratory native-
enzyme sequence design. The proposed reference study follows WT PA6 from ES
through acylation/TI1/acyl enzyme and deacylation/TI2/product formation. It does
not initially promise complete polymer association, diffusion or product-release
free-energy profiles. Long-polymer/material effects are separate from local
chemical activation. The catalytic bottleneck must be assessed across both
half-reactions, not assigned to TS1 beforehand.

For redesign, concentrate expensive validation on a reference mechanism and
small finalist set, not every sequence. Early screening uses chemically credible
state ensembles labelled with their actual provenance. Hybrid DFT PES checks
can discriminate mechanisms and test DFTB geometries/relative energy profiles
at representative states. They do not convert a DFTB PMF into a hybrid-DFT
free-energy profile. A PMF maximum is an ensemble bottleneck, not automatically
a first-order saddle on a single potential-energy surface. Preserve these
labels when comparing with the previous Hessian-based TS1 plan.

## Go/no-go stages

1. **Recover validated ES evidence.** Pin the prior PA6 graph, termini, register,
   force-field versions, AM1-BCC data, protein protonation, restraints and three
   valid WT trajectories. Select structurally distinct near-reactive snapshots
   from separate replicas. Existing MD provides starting ensembles, not automatic
   equilibration for a changed Hamiltonian or confidence in chemistry. Avoid
   repeating 3 × 100 ns wholesale unless an input change or deficient ensemble
   justifies it. Preserve the assembly and relevant intersubunit environment;
   one active QM site does not imply truncating away the other subunits.
2. **Validate a local reactive Hamiltonian before production sampling.** Define
   the QM region from plausible catalytic/proton networks, oxyanion donors,
   reacting amide and waters actually present in the chosen state. Map residues
   to native n12c and the NylC literature before assigning Asp roles. Keep cuts
   away from chemistry and audit link atoms, embedding, MM charges/exclusions
   and water LJ parameters. Test bond addition/elimination as well as proton
   transfer, state order and geometries against selected hybrid DFT PES models.
   Check QM-region and parameter-set sensitivity. DFTB3/3ob availability is not
   reaction-specific validation. No production umbrella campaign if it changes
   the mechanism qualitatively or fails the intended chemical states.
3. **Small sampling pilot.** Choose coupled CVs for bond formation/breaking and
   proton motion as required by competing mechanisms; evaluate whether water,
   side-chain or substrate reorganization supplies a hidden slow coordinate.
   Check window overlap, hysteresis, independent starts and estimator stability.
   A drawn one-dimensional reaction coordinate is not sufficient evidence.
   Follow the promising mechanism first; avoid exhaustive proton-path enumeration.
4. **WT chemical-cycle PMF.** Expand only after the pilot passes. Resolve both
   acylation and deacylation, report convergence/replica variation and limitations
   of the approximate Hamiltonian, and perform targeted hybrid DFT PES checks.
   Do not claim a hybrid DFT free-energy barrier or derive precise rates from
   an unbenchmarked DFTB surface.
5. **Mechanistic discrimination and redesign connection.** Choose a small number
   of substrate/register variants or residue perturbations that distinguish the
   mechanism and can be checked experimentally. Poly(MAXA) is a later extension,
   not a requirement for the first PA6 result. Feed validated local geometries
   and multistate incompatibilities into focused Nyl12 sequence design. Neither
   unchanged bound poses nor confidence scores replace chemical-barrier evidence.

## Parameterization: reuse without conflating two tasks

The existing noncovalent AM1-BCC workflow is reusable and is not invalidated by
our current Gaussian/RESP work. A RESP calculation is not automatically a better
model for every oligomer, nor is a charge-method change merely a software detail.
Use a declared consistent parameter family and matched checks if changing GAFF,
protein/water models or substrate charges. Preserve old ES simulations as their
own dataset instead of silently relabelling them RESP-compatible simulations.

State-specific TI1/acyl-enzyme MM models support separate classical state
sampling/refinement. A reactive QM/MM calculation instead needs the reacting
chemistry treated by QM, rather than held by an unchanging MM bond. It need not
wait for every bonded parameter of the reactive QM core to be fitted as a
standalone classical residue. This does not waive topology/nonbonded audits:
water and QM atoms still need the MM interactions required by the implementation,
and boundary bonds, link atoms and exclusion behavior must be verified as
chemistry changes. Do not use this distinction to run unvalidated TI1 classical
relaxation, which remains explicitly gated by the earlier user requirement.

## NylB: useful comparator, not a pivot

NylB has a published reactive QM/MM mechanism: Kamiya et al. (2014),
*Nylon-Oligomer Hydrolase Promoting Cleavage Reactions in Unnatural Amide
Compounds*, JPCL 5, 1210–1216, DOI 10.1021/jz500323y,
https://pubmed.ncbi.nlm.nih.gov/26274473/ . The authors report a Ser112/Lys115/
Tyr215 mechanism with a role for Tyr170 and an acylation barrier near 21 kcal/mol.
That is evidence for a possible methodological benchmark, not a transferable
Nyl12 barrier. NylB chemical catalysis and substrate recognition are different
from Ntn Thr Nyl12. Prior substrate-complex structures show a substantial
substrate-induced loop rearrangement (JMB 2007, publisher record
https://www.sciencedirect.com/science/article/pii/S0022283607005347 ).

Monomeric assembly simplifies bookkeeping, but does not remove conformational
sampling, protonation or substrate-specificity problems. Logging, workflow
execution, reaction-coordinate analysis, sampling/validation and design-ranking
infrastructure can transfer. Catalytic topology, parameters, QM boundaries,
residue mapping, chemical states and CVs must be rebuilt and tested. A new
NylB paper must distinguish itself from the existing QM/MM study; it is not
a blank-slate alternative that should displace the better-developed Nyl12 case.

DFTB validation reference: *Benchmark of Approximate Quantum Chemical and
Machine Learning Potentials for Biochemical Proton Transfer Reactions* (2025),
DOI 10.1021/acs.jctc.5c00690,
https://pubs.acs.org/doi/10.1021/acs.jctc.5c00690 . The primary indexed abstract
reports chemical-group-dependent performance. This motivates testing our actual
reaction and does not itself validate amide bond cleavage by a chosen DFTB set.
