# Introduction — working draft

Editorial status, 2026-10-20: conceptual framing and study objectives. No
accepted transition state, validated new force field, or improved enzyme is
claimed. Reference keys resolve in [REFERENCES.md](REFERENCES.md).

Enzymatic polyamide hydrolysis requires recognition of both polymer chemistry
and the local sequence surrounding a cleavable amide. This distinction becomes
particularly important for copolyamides, for which cleavage at different
positions can generate different soluble products. Dishner and colleagues
demonstrated enzyme-dependent product distributions during hydrolysis of the
sequence-controlled copolyamide poly(X,AMA), with Nyl12 favoring the MAXA tetrad
[R1]. These observations motivate a product-informed modeling strategy in
which the substrate graph and scissile-bond register are defined before
structural prediction or active-site redesign.

Protein–ligand structure prediction provides candidate geometries for such
studies. AlphaFold 3 extends complex prediction to diverse biomolecular
components, including covalent connections [R2]. However, predicted structural
confidence and compatibility with an imposed bond do not establish a catalytic
reaction pathway. Flexible oligomers can occupy several registers, and favorable
contacts involving remote substrate segments need not organize the reacting
amide. Our earlier Nyl12 modeling work further identified limits of interpreting
the persistence of a pre-bound pose as an explanation of relative activity:
bound-start simulations address local structural stability, whereas measured
product formation also depends on productive recognition and subsequent
chemical steps [local sources S1–S3].

For a putative covalent hydrolase, a useful design model must accommodate changes
in connectivity, protonation and geometry over the catalytic cycle. A stable
tetrahedral intermediate is not a prerequisite for nucleophilic acyl
substitution: addition, proton transfers and leaving-group departure can be
coupled without a distinct metastable tetrahedral basin, as calculated for
FAAH-catalyzed anandamide hydrolysis [R48]. The existence of such a basin is
therefore a substrate- and environment-dependent hypothesis, not a mandatory
design target. Product-connected pathways and the relevant transition states
must be sought whether or not a persistent TI1 is found. Path-generation
methods must be followed by stationary-point and endpoint validation.
Established Hessian-based saddle-search algorithms and newer optimization
methods provide complementary approaches to this task [R3,R12,R13].

The present workflow is being developed to connect product-defined substrate
chemistry, audited covalent-state models and QM/MM reaction geometries to
structure-conditioned sequence design. Catalytic motifs will guide focused
LigandMPNN redesign and, where required, RFdiffusion3 backbone diversification
[R10,R11]. Candidate sequences will be assessed in both covalent and independent
noncovalent complexes, accounting for the four equivalent active sites of the
mature Nyl12 assembly. Additional checks on acyl-enzyme hydrolysis and product
accessibility will test whether improvement around the first transition state
creates an unfavorable later step. The objective is a reproducible route to
experimentally testable designs, with explicit separation between structural
plausibility, electronic-energy barriers and catalytic performance.

Mechanistic models for N-terminal nucleophile hydrolases do not establish a
universal requirement for a water bridge: water-assisted and direct-transfer
models have been reported in other family members, and a carboxylate-containing
network has been proposed for NylC [R23–R26]. We therefore treat relay geometry
and protonation as hypotheses requiring structural and chemical validation.
An unsuccessful placement search in one fixed conformation does not exclude
transient solvation or establish an alternative mechanism.

Negoro et al. proposed a triad-mediated substrate-hydrolysis mechanism for
NylC in which Asp308 activates the Thr267 hydroxyl through Asp306 and the
neutral alpha-amino group of Thr267 [R25, Fig. 10]. This proposal motivates
testing the corresponding Asp274–Asp272–Thr233 network in Nyl12. Their
assignment of protonated Asp306 is a mechanistic inference from the structural
environment and proposed acid/base chemistry, rather than a direct measurement
of its protonation state. Moreover, loss of activity in Asp mutants is partly
confounded by failed autocatalytic maturation: an inactive uncleaved precursor
does not independently establish the residue's catalytic role in mature-enzyme
substrate hydrolysis. These limitations motivate explicit protonation-state
comparisons and separate evaluation of maturation and substrate chemistry.

Negoro et al. also found that K189A in the H-GYAQ NylC background strongly
reduced Ahx-cyclic-oligomer activity: its cleavage-normalized specific activity
was 0.032 U mg⁻¹, about 0.6% of parental H-GYAQ, while processing reached only
33% after further incubation [R25, Table 1, Fig. 8D]. This supports an
important role in NylC enzyme function and hydrolysis beyond the maturation
defect, although the exact chemical contribution remains proposed rather than
isolated from every structural or maturation effect.

Human asparaginase type III provides a particularly relevant multiscale
substrate-hydrolysis precedent for an N-terminal threonine nucleophile [R41].
In the selected acylation pathway, proton transfer from the threonine hydroxyl
to its own neutral N-terminal amine precedes nucleophilic attack; the reported
free-energy barriers were 14.0 and 18.6 kcal/mol, respectively. The authors
used DFTB3-D3/MM to screen mechanisms and B3LYP-D3/6-31+G*/MM adaptive-string
calculations for the selected free-energy path. This result motivates testing
direct Thr233 self-activation in Nyl12 before imposing an Asp272 shuttle, while
neither proving mechanistic equivalence nor conflating mature-substrate
hydrolysis with the separate autocatalytic maturation reaction.

This stepwise precedent does not imply that the same intermediates must exist
for a flexible secondary polyamide amide. In contrast, the FAAH calculations
found that removing a constraint from an early proton-transfer geometry returned
the system to reactants, although the complete calculated path accomplished
acylation without a metastable intermediate [R48]. Moreover, “concerted” can
refer solely to deprotonation accompanying attack, while a tetrahedral minimum
still forms, as in a Zika serine-protease study [R46]. Mechanistic comparisons
must specify which bond changes are coupled. Amide nitrogen pyramidalization,
inversion and access to a proton donor can be as important as carbonyl attack
and oxyanion stabilization. Modeling and enzyme-engineering studies implicate
a scissile-NH hydrogen-bond acceptor in facilitating this chemistry [R44,R45].

The broad nylon activity found among diverse Ntn-hydrolases without evidence
of prior selection for nylon degradation supports a promiscuous-activity
framework [R47]. For Nyl12, a poor substrate pose, a rarely populated catalytic
microstate or inefficient proton delivery could each limit activity; low
polymer conversion alone cannot identify the rate-limiting chemical barrier.
The correspondence of NylC Tyr146 to Nyl12 Phe113 is mechanistically relevant:
the backbone NH oxyanion donor is retained, whereas the phenolic OH assigned
a possible additional hydrolytic role in NylC is absent [R25; local mapping
in docs/negoro_2023_mechanism_audit.md]. This is a hypothesis about alternative
catalytic networks, not proof that introducing Tyr will improve Nyl12.

Recent generative enzyme-design studies motivate complementing local chemical
models with structural assessment across reaction states [R27–R31]. These
approaches do not convert a predicted intermediate or transition-state analogue
into a verified electronic saddle point. Native-fold and polymer-hydrolase
engineering provide additional precedents for focused redesign and conformer-
aware evaluation [R32,R33]. For Nyl12, this distinction is essential: local
amide chemistry must be organized while long substrate segments retain access,
conformational diversity and product-specific cleavage registers. We propose
separating economical catalytic-state screening from detailed reference and
shortlist reaction-path calculations, while retaining the maturation and
oligomeric architecture of the existing enzyme.

Efficient QM/MM models make it practical to explore coupled proton transfer and
amide-bond chemistry within a solvated polymer-binding pocket, but their
energetic accuracy requires reaction-specific assessment. DFTB3/3ob-3-1 provides
a useful starting Hamiltonian: biochemical proton-transfer benchmarks show
chemical-group-dependent errors [R35], and enzyme-model benchmarks identify
non-negligible errors in both barriers and reaction energies [R36]. Successful
DFTB3 sampling with higher-level energy corrections in PETase illustrates its
potential for polymer-hydrolase studies, although ester chemistry does not
validate nylon amide hydrolysis [R37]. Accordingly, we propose comparing DFTB3
with omegaB97X-D for representative proton-transfer models and selected
embedded pathway geometries. Numerical convergence, intermediate stability and
biochemical plausibility will be assessed separately; an asynchronous or
energetically unfavorable driven path is not by itself evidence of optimizer
failure or of the native catalytic mechanism.
