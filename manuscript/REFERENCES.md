# Working references

Started 2026-10-20. Stable DOI keys are used in the drafts. This is a working
bibliography; journal formatting and complete author lists remain to be added.
R7–R9 identifiers are consistent with prior project citations; publisher
full-text retrieval was incomplete during this checkpoint, so detailed claims
must be checked against the articles before submission. Software documentation
is cited as such, not represented as a peer-reviewed benchmark.

- **R1.** Dishner et al. *Circularity in Sequence-Controlled Copolyamides Enabled
  by Regioselective Enzymatic Hydrolysis*. JACS (2026).
  https://doi.org/10.1021/jacs.6c03466
- **R2.** Abramson et al. *Accurate structure prediction of biomolecular
  interactions with AlphaFold 3*. Nature (2024).
  https://doi.org/10.1038/s41586-024-07487-w
- **R3.** Kästner, Carr, Keal, Thiel, Wander and Sherwood. *DL-FIND: An
  Open-Source Geometry Optimizer for Atomistic Simulations*. J. Phys. Chem. A
  (2009), 113, 11856–11865. https://doi.org/10.1021/jp9028968
- **R4.** Chai and Head-Gordon. *Long-range corrected hybrid density functionals
  with damped atom–atom dispersion corrections*. PCCP (2008), 10, 6615–6620.
  https://doi.org/10.1039/B810189B
- **R5.** Weigend and Ahlrichs. *Balanced basis sets of split valence, triple
  zeta valence and quadruple zeta valence quality for H to Rn: Design and
  assessment of accuracy*. PCCP (2005). https://doi.org/10.1039/B508541A
- **R6.** Robertson, Tirado-Rives and Jorgensen. *Improved Peptide and Protein
  Torsional Energetics with the OPLS-AA Force Field*. JCTC (2015), 11, 3499–3509.
  https://doi.org/10.1021/acs.jctc.5b00356
- **R7.** Bayly, Cieplak, Cornell and Kollman. *A Well-Behaved Electrostatic
  Potential Based Method Using Charge Restraints for Deriving Atomic Charges:
  The RESP Model*. J. Phys. Chem. (1993), 97, 10269–10280.
  https://doi.org/10.1021/j100142a004
- **R8.** Tian et al. *ff19SB: Amino-Acid-Specific Protein Backbone Parameters
  Trained against Quantum Mechanics Energy Surfaces in Solution*. JCTC (2020).
  https://doi.org/10.1021/acs.jctc.9b00591
- **R9.** Gaus, Goez and Elstner. *Parametrization and Benchmark of DFTB3 for
  Organic Molecules*. JCTC (2013), 9, 338–354.
  https://doi.org/10.1021/ct300849w
- **R10.** Dauparas et al. *Atomic context-conditioned protein sequence design
  using LigandMPNN*. Nature Methods (2025).
  https://doi.org/10.1038/s41592-025-02626-1
- **R11.** RosettaCommons. *Enzyme Design in RFdiffusion3*, official Foundry
  documentation, accessed 2026-10-20.
  https://rosettacommons.github.io/foundry/models/rfd3/tutorials/enzyme_design_tutorial.html
  Additional official input/symmetry documentation:
  https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/input.md
  and the sibling docs/examples/symmetry.md.
  See R31 for the verified final 2026 primary-paper record. Any production
  design must additionally pin its actual software release/commit.
- **R12.** Hermes, Sargsyan, Najm and Zádor. *Sella, an Open-Source
  Automation-Friendly Molecular Saddle Point Optimizer*. JCTC (2022), 18,
  6974–6988. https://doi.org/10.1021/acs.jctc.2c00395
- **R13.** Steinmetzer et al. *pysisyphus: Exploring potential energy surfaces
  in ground and excited states*. Int. J. Quantum Chem. (2021).
  https://doi.org/10.1002/qua.26390
- **R14.** Allen, Payne and Cole. *Harmonic Force Constants for Molecular
  Mechanics Force Fields via Hessian Matrix Projection*. JCTC (2018), 14,
  274–281. https://doi.org/10.1021/acs.jctc.7b00785
  Publisher abstract and author-repository bibliographic record verified
  2026-10-22; full-text retrieval timed out. This is the modified-Seminario
  candidate method, not evidence that the present TI1 parameters are validated.

Before submission finalize software versions for production design/prediction,
the original Bürgi–Dunitz/Herschlag metric references and original Nyl12
discovery/structure/maturation references. The present projected-phi implementation is
identified by its actual local source rather than an unverified literature
attribution. These open bibliographic tasks do not justify inventing references.


## Force-field parameterization references added 2026-10-22

- **R15.** Wang, J.; Wolf, R. M.; Caldwell, J. W.; Kollman, P. A.; Case, D. A.
  *Development and testing of a general amber force field*. J. Comput. Chem.
  **2004**, 25 (9), 1157–1174. https://doi.org/10.1002/jcc.20035
  Original GAFF reference; do not identify this as a GAFF2 development paper.
- **R16.** Wang, J.; Wang, W.; Kollman, P. A.; Case, D. A.
  *Automatic atom type and bond type perception in molecular mechanical
  calculations*. J. Mol. Graph. Model. **2006**, 25 (2), 247–260.
  https://doi.org/10.1016/j.jmgm.2005.12.005
  Antechamber atom/bond-perception reference. The journal component is `j.jmgm`,
  not `j.theochem`; the latter DOI was rejected during bibliography checking.
- **R17.** Case, D. A.; Aktulga, H. M.; Belfon, K.; et al. *AmberTools*.
  J. Chem. Inf. Model. **2023**, 63 (20), 6183–6191.
  https://doi.org/10.1021/acs.jcim.3c01153
  Software citation recommended by https://ambermd.org/CiteAmber.php.
  This paper is the software citation, not a claim that its publication year
  identifies the installed AmberTools26 parameter-file revision.
- **R18.** Amber development team. *Ligand Force Fields* and *The Amber Force
  Fields*, official Amber documentation, accessed 2026-10-22.
  https://ambermd.org/AmberModels_organic.php and
  https://ambermd.org/AmberModels.php; linked *Amber 2026 Reference Manual*.
  GAFF2 implementation identifier: installed `gaff2.dat`, version 2.2.30,
  October 2025, distributed in the present AmberTools26 environment. The exact
  file hash is retained with manuscript-source evidence. The official page
  describes GAFF2 as an Amber/AmberTools distribution; no standalone GAFF2
  development paper was verified in this check.
- **R19.** Jorgensen, W. L.; Chandrasekhar, J.; Madura, J. D.; Impey, R. W.;
  Klein, M. L. *Comparison of simple potential functions for simulating liquid
  water*. J. Chem. Phys. **1983**, 79 (2), 926–935.
  https://doi.org/10.1063/1.445869
  TIP3P reference. Publisher full-text retrieval was unavailable during this
  check; the implemented TIP3P values are independently traced to local Amber
  source units, rather than inferred from an inaccessible article.
- **R20.** Frisch, M. J.; Trucks, G. W.; Schlegel, H. B.; et al.
  *Gaussian 16, Revision C.02*; Gaussian, Inc.: Wallingford, CT, **2019**.
  Revision, year and author-citation block were read directly from the actual
  calculation's "Cite this work as" header; build identifier
  ES64L-G16RevC.02, 7-Dec-2021. Official citation page:
  https://gaussian.com/citation/ (retrieval unavailable at this checkpoint).
- **R21.** Singh, U. C.; Kollman, P. A. *An approach to computing electrostatic
  charges for molecules*. J. Comput. Chem. **1984**, 5 (2), 129–145.
  https://doi.org/10.1002/jcc.540050204
- **R22.** Besler, B. H.; Merz, K. M., Jr.; Kollman, P. A. *Atomic charges
  derived from semiempirical methods*. J. Comput. Chem. **1990**, 11 (4),
  431–439. https://doi.org/10.1002/jcc.540110404
  R21/R22 document the ESP charge-fitting lineage; citing R22 does not mean
  semiempirical ESP replaces the specified HF/6-31G(d) calculation here.

R15–R18 and R21/R22 bibliographic/method identifiers were checked against
publisher or official software sources. This is not a claim that all cited
full texts were obtained. Full author-list and journal-style formatting can be
completed when the final manuscript and accepted software manifests are frozen.

## Mechanistic alternatives and recent enzyme design — added 2026-10-24

- **R23.** Wei et al.. *Fundamental reaction pathway for peptide metabolism by proteasome: insights from first-principles quantum mechanical/molecular mechanical free energy calculations*. J. Phys. Chem. B (2013).
  https://doi.org/10.1021/jp405337v
  Primary access record: https://pubmed.ncbi.nlm.nih.gov/24111489/
- **R24.** Chilov, Sidorova and Svedas. *Quantum chemical studies of the catalytic mechanism of N-terminal nucleophile hydrolase*. Biochemistry (Moscow) (2007).
  https://doi.org/10.1134/S0006297907050057
  Primary access record: https://pubmed.ncbi.nlm.nih.gov/17573703/
- **R25.** Negoro et al.. *X-ray crystallographic and mutational analysis of the NylC precursor: catalytic mechanism of autocleavage and substrate hydrolysis of nylon hydrolase*. FEBS Journal (2023).
  https://doi.org/10.1111/febs.16755
  Primary access record: https://febs.onlinelibrary.wiley.com/doi/10.1111/febs.16755
- **R26.** *A unified mechanism for proteolysis and autocatalytic activation in the 20S proteasome*. Nature Communications (2016).
  https://doi.org/10.1038/ncomms10900
  Primary access record: https://pmc.ncbi.nlm.nih.gov/articles/PMC4792962/
- **R27.** Lauko et al.. *Computational design of serine hydrolases*. Science (2025).
  https://doi.org/10.1126/science.adu2454
  Primary access record: https://pmc.ncbi.nlm.nih.gov/articles/PMC12288761/
- **R28.** Kim et al.. *Computational design of metallohydrolases*. Nature (2026).
  https://doi.org/10.1038/s41586-025-09746-w
  Primary access record: https://www.nature.com/articles/s41586-025-09746-w
- **R29.** Choi et al.. *De novo design of cysteine proteases*. Nature (2026).
  https://doi.org/10.1038/s41586-026-11143-w
  Primary access record: https://www.nature.com/articles/s41586-026-11143-w
- **R30.** Chen et al.. *Computational design of metalloproteases*. Nature (2026).
  https://doi.org/10.1038/s41586-026-11167-2
  Primary access record: https://www.nature.com/articles/s41586-026-11167-2
- **R31.** Butcher et al.. *Atomic functional site scaffolding in biopolymers with RFdiffusion3*. Science (2026).
  https://doi.org/10.1126/science.aeh2982
  Primary access record: https://pubmed.ncbi.nlm.nih.gov/42879533/
- **R32.** Wang et al.. *Computational loop reconstruction based design of efficient PET hydrolases*. Communications Biology (2025).
  https://doi.org/10.1038/s42003-025-08364-6
  Primary access record: https://pmc.ncbi.nlm.nih.gov/articles/PMC12174350/
- **R33.** Listov et al.. *Complete computational design of high-efficiency Kemp elimination enzymes*. Nature (2025).
  https://doi.org/10.1038/s41586-025-09136-2
  Primary access record: https://www.nature.com/articles/s41586-025-09136-2

R23/R24 support alternative relay/direct-transfer models in other enzymes,
not an assigned Nyl12 mechanism. R25 separates precursor autocleavage from
mature substrate hydrolysis. R26 discusses competing proteasome general-base
models. Full author lists and publication formatting remain editorial tasks.
R27 is the final Science article, not the 2024 ChemNet preprint. R28 appeared
online 2025-12-03 and in the January 2026 issue. R29/R30 appeared online
2026-10-14; R31's primary PubMed record gives online 2026-10-15. R30 detailed
motif provenance was checked in the authors' preprint, not final SI. Final
R31 full Methods/SI were inaccessible; final abstract and official code/docs
are distinguished from the older preprint. Do not transplant preprint metrics.
Official implementation sources: https://github.com/hojaec/RFD2_MI_cysteine_protease_pipeline
and https://rosettacommons.github.io/foundry/models/rfd3/tutorials/advanced_enzyme_design_tutorial.html .
Additional preserved relay structural reference (full bibliographic details
pending): *Crystal structures of penicillin acylase enzyme-substrate complexes:
structural insights into the catalytic mechanism* (2001), publisher record
https://www.sciencedirect.com/science/article/pii/S0022283601950430 .

## Proposed mechanistic-project references — 2026-10-24

- **R34.** Kamiya, K.; Baba, T.; Boero, M.; Matsui, T.; Negoro, S.; Shigeta, Y.
  *Nylon-Oligomer Hydrolase Promoting Cleavage Reactions in Unnatural Amide
  Compounds*. J Phys Chem Lett (2014),5,1210–1216.
  https://doi.org/10.1021/jz500323y ; primary metadata/abstract:
  https://pubmed.ncbi.nlm.nih.gov/26274473/ . This is NylB, not Nyl12, and
  is not evidence for a DFTB Hamiltonian until its actual methods are checked.
- **R35.** *Benchmark of Approximate Quantum Chemical and Machine Learning
  Potentials for Biochemical Proton Transfer Reactions*. JCTC (2025).
  https://doi.org/10.1021/acs.jctc.5c00690 . Primary indexed abstract verified;
  full PMC text retrieval returned CAPTCHA. Author-list formatting pending.
  Chemical-group-dependent proton-transfer benchmarking does not by itself
  validate a polyamide-cleavage free-energy surface.
