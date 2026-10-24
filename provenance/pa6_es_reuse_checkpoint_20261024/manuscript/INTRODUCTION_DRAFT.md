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

For a covalent hydrolase, a useful design model must accommodate changes in
connectivity, protonation and geometry over the catalytic cycle. A tetrahedral
adduct can provide a chemically explicit starting point from which the
preceding transition structure and noncovalent enzyme–substrate complex can be
investigated. This construction is a mechanistic hypothesis to test, rather
than proof that an imposed intermediate is stable. Path-generation methods
must therefore be followed by stationary-point and endpoint validation.
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
