# Working references

## Primary-amine release assay references (verified 2026-11-26)

- **R61.** Negoro, S.; Shibata, N.; Tanaka, Y.; et al.
  *Three-dimensional Structure of Nylon Hydrolase and Mechanism of Nylon-6
  Hydrolysis*. J. Biol. Chem. **2012**, 287 (7), 5079–5090.
  https://doi.org/10.1074/jbc.M111.321992 .
  Direct NylC precedent: enzyme activity was measured as an increase in amino
  groups using TNBS; polymer-product analysis also used TLC/ninhydrin. These
  are different readouts, not a proof that all soluble products are monomers.
  Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3281642/ .
- **R62.** Shiraishi, Y.; Kato, D.-i.; Miyazaki, K.; Yonemura, M.; Furuno, Y.;
  Yokoyama, R.; et al. *Quantitative nylon monomerization by the combination
  of chemical pretreatment and enzymatic hydrolysis using nylon hydrolases*.
  PLOS ONE **2025**, 20 (2), e0318641.
  https://doi.org/10.1371/journal.pone.0318641 .
  Section 1.4 gives a 96-well TNBS assay in tetraborate at pH 9.5, with 420 nm
  detection after 40 °C derivatization and calibration against 6-aminohexanoate.
  Hydrolysis conditions are separate from these detection conditions.
  This is a nylon-specific adaptation precedent, not a validated Nyl12 assay
  in 200 mM Tris/crude lysate. First author is Shiraishi, not Negoro.

## Core computational references

- **R43.** Kästner, J.; Sherwood, P. *Superlinearly converging dimer method
  for transition state search*. J. Chem. Phys. **2008**, *128*, 014106.
  https://doi.org/10.1063/1.2815812
  Bibliographic details and method abstract verified 2026-11-13.

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
  Arantes, G. M.; Rezac, J. JCTC 21(14), 7149–7159.
  https://doi.org/10.1021/acs.jctc.5c00690 . Full primary text read:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12288007/ . Uses DFTB3/3ob-3-1;
  chemistry-dependent errors do not certify coupled amide attack in Nyl12.
- **R36.** Kromann et al. *Towards a barrier height benchmark set for
  biologically relevant systems*. PeerJ (2016), 4, e1994.
  https://doi.org/10.7717/peerj.1994 . DFTB3 enzyme-model barrier and reaction
  energy benchmarking against DFT; not a direct Nyl12 benchmark.
- **R37.** Berselli, A.; Menziani, M. C.; Piccini, G. M.; Muniz-Miranda, F.
  *Molecular-Level Mechanistic Insights into PETase-Catalyzed Plastics Hydrolysis
  from Accurate QM/MM Free Energy Calculations*. ACS Catalysis (2025), 15,
  10702–10721. https://doi.org/10.1021/acscatal.5c01263 . Primary publisher
  indexed methods/results checked; DFTB3 sampling with higher-level single-point
  reweighting. Ester chemistry, not amide hydrolysis.

- **R38.** Froitzheim, T.; Müller, M.; Hansen, A.; Grimme, S.
  *g-xTB: A General-Purpose Extended Tight-Binding Electronic Structure Method
  For the Elements H to Lr (Z=1–103)*. ChemRxiv (2025), version 1, June 24.
  https://doi.org/10.26434/chemrxiv-2025-bjxvt . Preprint, not treated as
  peer-reviewed validation. Primary indexed supporting Table 3 compares full
  GMTKN55 WTMAD-2; current development software differs from that model.
- **R39.** Giese, T. J.; Zeng, J.; Lerew, L.; McCarthy, E.; Tao, Y.; Ekesan, Ş.;
  York, D. M. *Software Infrastructure for Next-Generation QM/MM–ΔMLP Force
  Fields*. J Phys Chem B (2024), 128, 6257–6271.
  https://doi.org/10.1021/acs.jpcb.4c01466 . Primary author-hosted full text
  reviewed. Direct sander/GFN2-xTB interface with electrostatic embedding and
  QM/MM PME; the uncorrected Hamiltonian does not require MLP training.
- **R40.** Bannwarth, C.; Ehlert, S.; Grimme, S. *GFN2-xTB—An Accurate and
  Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method
  with Multipole Electrostatics and Density-Dependent Dispersion Contributions*.
  JCTC (2019), 15, 1652–1671. https://doi.org/10.1021/acs.jctc.8b01176 .
- **R41.** Andjelkovic, M.; Zinovjev, K.; Ramos-Guzmán, C. A.; Ruiz-Pernía,
  J. J.; Tuñón, I. *Elucidation of the Active Form and Reaction Mechanism in
  Human Asparaginase Type III Using Multiscale Simulations*. J Chem Inf Model
  (2023), 63, 5676–5688. https://doi.org/10.1021/acs.jcim.3c00900 . Primary
  full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC10852353/ . This is a
  mature-substrate hydrolysis precedent; its discussion of self-cleavage as an
  activation prerequisite must not be substituted for the modeled catalytic
  substrate-hydrolysis steps.
- **R42.** Tao, Y.; Giese, T. J.; Ekesan, Ş.; Zeng, J.; Aradi, B.; Hourahine,
  B.; Aktulga, H. M.; Götz, A. W.; Merz, K. M., Jr.; York, D. M. *Amber Free
  Energy Tools: Interoperable Software for Free Energy Simulations Using
  Generalized Quantum Mechanical/Molecular Mechanical and Machine Learning
  Potentials*. J Chem Phys (2024), 160, 224104.
  https://doi.org/10.1063/5.0211276 . Primary full text reviewed. Establishes
  direct Sander/xTB and Sander/DFTB+ interfaces with QM/MM Ewald and a
  DeePMD-kit correction interface; its G-T proton-transfer example also shows
  that uncorrected DFTB3 and GFN2-xTB profiles can differ substantially from
  the ab initio reference.

Implementation documentation checked 2026-10-25: g-xTB development README,
https://github.com/grimme-lab/g-xtb ; ASH xTB interface,
https://ash.readthedocs.io/en/latest/xTB-interface.html ; ORCA semiempirical
methods, https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/semiempirical.html .
The current g-xTB README explicitly excludes point-charge embedding; ASH
likewise states that electrostatically embedded g-xTB QM/MM is unavailable.
Analytic gradients are available in the current development binary, unlike
the older implementation. Do not conflate the two releases.

- **R43.** Kästner, J.; Sherwood, P. *Superlinearly converging dimer method
  for transition state search*. J Chem Phys (2008), 128, 014106.
  https://doi.org/10.1063/1.2815812 . Existing BibTeX key retained; reference
  list synchronized 2026-11-18.

## Amide hydrolysis without a required TI1 minimum — reviewed 2026-11-18

- **R48.** Palermo, G.; Campomanes, P.; Cavalli, A.; Rothlisberger, U.;
  De Vivo, M. *Anandamide Hydrolysis in FAAH Reveals a Dual Strategy for
  Efficient Enzyme-Assisted Amide Bond Cleavage via Nitrogen Inversion*.
  J Phys Chem B (2015), 119, 789–801. https://doi.org/10.1021/jp5052276 .
  Primary abstract: https://pubmed.ncbi.nlm.nih.gov/25205244/ ; author manuscript:
  https://infoscience.epfl.ch/record/208620/files/Full%20Text.pdf . Indexed primary
  full-text passages on pp. 793 and 795 explicitly report reversal of the early
  PT1 state on constraint release and no metastable intermediates within the
  calculated acylation/deacylation steps. Their conclusion is model-dependent
  and does not establish the Nyl12 mechanism. Direct PDF access returned
  429/redirect errors; the primary indexed passages and PubMed abstract were
  available. Distinguish a concerted full acylation from concerted addition only.
- **R44.** Syrén, P.-O.; Hult, K. *Amidases Have a Hydrogen Bond that Facilitates
  Nitrogen Inversion, but Esterases Have Not*. ChemCatChem (2011), 3, 853–860.
  https://doi.org/10.1002/cctc.201000448 . Primary modeling/experimental study
  proposes a scissile-NH hydrogen bond to an enzyme/substrate acceptor that
  facilitates nitrogen inversion. A mechanistic hypothesis and design precedent,
  not proof of a universal requirement or an observed Nyl12 contact.
- **R45.** Syrén, P.-O.; Hendil-Forssell, P.; Aumailley, L.; Besenmatter, W.;
  Gounine, F.; Svendsen, A.; Martinelle, M.; Hult, K. *Esterases with an
  Introduced Amidase-Like Hydrogen Bond in the Transition State Have Increased
  Amidase Specificity*. ChemBioChem (2012), 13, 645–648.
  https://doi.org/10.1002/cbic.201100779 . Publisher graphical abstract reports
  up to 50-fold improved amide-over-ester hydrolysis specificity after rational
  introduction of the proposed interaction; this is a specificity ratio, not
  evidence for a 50-fold absolute Nyl12 turnover improvement.
- **R46.** Nutho, B.; Mulholland, A. J.; Rungrotmongkol, T. *Quantum
  Mechanics/Molecular Mechanics (QM/MM) Calculations Support a Concerted
  Reaction Mechanism for the Zika Virus NS2B/NS3 Serine Protease with Its
  Substrate*. J Phys Chem B (2019), 123, 2889–2903.
  https://doi.org/10.1021/acs.jpcb.9b02157 . Primary author record/abstract:
  https://research-information.bris.ac.uk/en/publications/quantum-mechanicsmolecular-mechanics-qmmm-calculations-support-a- .
  Here “concerted” denotes nucleophile deprotonation plus attack; a tetrahedral
  intermediate still forms. It is not evidence for intermediate-free acylation.
- **R47.** Drufva, E. E.; Cahill, J. F.; Saint-Vincent, P. M. B.; et al.
  *Identification and characterization of substrate- and product-selective
  nylon hydrolases*. Chem Catalysis (2025), 5, 101418.
  https://doi.org/10.1016/j.checat.2025.101418 . Published primary study of diverse
  Ntn-hydrolases supports widespread nylon promiscuity and substrate/product
  selectivity without evidence of prior selection for this activity. Use this
  published version rather than the earlier 2024 preprint. It does not establish
  the native function or microscopic hydrolysis mechanism of Nyl12.

- **R48.** Gumbart, J. C.; Roux, B.; Chipot, C. *Standard Binding Free Energies
  from Computer Simulations: What Is the Best Strategy?* J Chem Theory Comput
  (2013), 9, 794–802. https://doi.org/10.1021/ct3008099 . Primary methodological
  reference for explicit standard-state and restraint corrections in binding/
  dissociation free energies. Relevant when comparing product-bound and
  product-released acylenzyme states; deleting a fragment is not a computed
  product-release free energy.

Additional primary historical counterexample, retained with limitations:
Harrison, M. J.; Burton, N. A.; Hillier, I. H.; Gould, I. R. *Mechanism and
transition state structure for papain catalysed amide hydrolysis, using a
hybrid QM/MM potential*. Chem Commun (1996), 2769–2770.
https://doi.org/10.1039/CC9960002769 . Its AM1/MM concerted, no-TI prediction
is a cysteine/thiol example using an early approximate model, not direct
validation of hydroxyl-nucleophile Nyl12 chemistry.


- **R49.** Rosta, E.; Woodcock, H. L.; Brooks, B. R.; Hummer, G.
  *Artificial reaction coordinate “tunneling” in free-energy calculations:
  the catalytic reaction of RNase H*. J Comput Chem (2009), 30, 1634–1641.
  https://doi.org/10.1002/jcc.21312 . Primary methodological evidence that an
  omitted slow coordinate can distort a reaction PMF; not direct Nyl12
  mechanistic precedent.
- **R50.** Rivas-Fernández et al. *Efficient QM/MM Modeling of Enzymatic
  Reactions Combining PathCV with OPES*. J Chem Theory Comput (2026), 22,
  5867–5881. https://doi.org/10.1021/acs.jctc.6c00492 . Potential future
  pathway/free-energy methodology; not used in the present calculations.
- **R51.** IUPAC Gold Book. *Transition-state theory*.
  https://goldbook.iupac.org/terms/view/T06470 . Defines the relationship among
  elementary rate constants, transmission coefficient and activation free energy.
- **R52.** Sigma-Aldrich. *TRIS buffer technical data*.
  https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/product/documents/142/961/93286dat.pdf .
  Manufacturer reference for temperature-sensitive Tris pH; not a measurement
  of the collaborators' assay pH at 75 °C.

- **R53.** Zoi, I.; Antoniou, D.; Schwartz, S. D. *Incorporating Fast Protein
  Dynamics into Enzyme Design: A Proposed Mutant Aromatic Amine Dehydrogenase*.
  J Phys Chem B (2017), 121, 7290–7298.
  https://doi.org/10.1021/acs.jpcb.7b05319 . Primary full text:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5763476/ . Uses committor analysis and
  15 uncorrelated TS configurations per enzyme in a theoretical engineering
  study; does not demonstrate a LigandMPNN workflow or experimental validation
  of the proposed mutant. Distinguish the 2017 publication from its 2018 PMC
  availability date. Together with R27/R29, motivates a proposed Nyl12 hybrid
  PES-reference/transition-region-ensemble comparison, not a claim it has run.

TS-template source recheck (2026-11-25): R27 primary full text distinguishes
the QM-derived motif from stochastic PLACER ensembles across catalytic states;
R29's indexed primary Extended Data Fig. 9 reports 50 PLACER models per state.
These are learned structural ensembles, not measured thermodynamic populations
or automatically committor-qualified TSEs. The R29 publisher full page was
inaccessible during this recheck; no additional final-SI claims were made.
R10 primary model description and the official LigandMPNN README
(https://github.com/dauparas/LigandMPNN, section `--pdb_path_multi`) distinguish
coordinate/element conditioning and efficient independent input batching from
a joint multistate optimization objective.

- **R54.** Lodola, A.; Branduardi, D.; De Vivo, M.; Capoferri, L.; Mor, M.;
  Piomelli, D.; Cavalli, A. *A Catalytic Mechanism for Cysteine N-Terminal
  Nucleophile Hydrolases, as Revealed by Free Energy Simulations*.
  PLoS ONE (2012), 7, e32397. https://doi.org/10.1371/journal.pone.0032397 .
  Primary QM/MM free-energy study with an N-terminal proton donor and
  chair-like transition-region geometry; the tetrahedral adduct is metastable
  rather than a well-defined free-energy minimum. This cysteine/thioester
  system is mechanistic precedent, not quantitative validation of the
  threonine/ester Nyl12 pathway. Full text:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3289653/ .

## Equilibration, truncation and sampling uncertainty (verified 2026-11-26)

- **R55.** Chodera, J. D. *A Simple Method for Automated Equilibration Detection
  in Molecular Simulations*. J. Chem. Theory Comput. **2016**, 12 (4),
  1799–1805. https://doi.org/10.1021/acs.jctc.5b00784 .
  Primary method for choosing a truncation time by maximizing an estimated
  effectively uncorrelated sample count; implemented in PyMBAR's
  `timeseries.detect_equilibration`. This does not validate our additional
  latest-of-thirteen/common-maximum-g aggregation. Full text:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4945107/ .
- **R56.** Clark, F.; Cole, D. J.; Michel, J. *Robust Automated Truncation Point
  Selection for Molecular Simulations*. J. Chem. Theory Comput. **2025**,
  21 (1), 88–101. https://doi.org/10.1021/acs.jctc.4c01359 .
  Online 23 December 2024; issue date 14 January 2025. Primary comparison
  documenting late/variable truncation and spurious optima from noisy
  correlation estimates. Distinguishes truncation selection from detecting
  insufficient sampling. Synthetic series modeled on protein–ligand
  free-energy calculations; not direct validation of Nyl12 QM/MM sampling.
  Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC11736681/ .
- **R57.** Grossfield, A.; Patrone, P. N.; Roe, D. R.; Schultz, A. J.;
  Siderius, D. W.; Zuckerman, D. M. *Best Practices for Quantification of
  Uncertainty and Sampling Quality in Molecular Simulations*, Article v1.0.
  Living J. Comput. Mol. Sci. **2019**, 1 (1), 5067.
  https://doi.org/10.33011/livecoms.1.1.5067 .
  General guidance on statistical uncertainty and sampling-quality assessment;
  not the source of the Chodera truncation objective or our aggregation rule.
  The journal places this article in the 2019 issue; its landing-page metadata
  gives online publication/citation year 2018, while the downloadable v1.0
  guide gives 2019. The issue-year convention is used here, retaining this
  metadata distinction for final reference-manager reconciliation. Publisher:
  https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067 .

- **R58.** Lu, X.; Fang, D.; Ito, S.; Okamoto, Y.; Ovchinnikov, V.; Cui, Q.
  *QM/MM free energy simulations: recent progress and challenges*.
  Mol. Simul. **2016**, 42 (13), 1056–1078.
  https://doi.org/10.1080/08927022.2015.1132317 .
  QM/MM-specific review emphasizing environmental sampling, enhanced sampling,
  collective-variable choice and limitations of approximate Hamiltonians.
  Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC4993472/ .
- **R59.** Jäckering, A.; van der Kamp, M.; Strodel, B.; Zinovjev, K.
  *Influence of Wobbling Tryptophan and Mutations on PET Degradation Explored
  by QM/MM Free Energy Calculations*. J. Chem. Inf. Model. **2024**, 64 (19),
  7544–7554. https://doi.org/10.1021/acs.jcim.4c00776 .
  Adaptive-string/DFTB3/MM polymer-hydrolase study including an orientation
  coordinate and umbrella integration; used a TS-region 95% confidence interval
  criterion of less than 1 kcal/mol. Also illustrates QM-method sensitivity of
  proton-transfer and tetrahedral-state energetics. This criterion is theirs,
  not sufficient alone to qualify the present MBAR calculation. Primary PDF:
  https://juser.fz-juelich.de/record/1031817/files/j%C3%A4ckering-et-al-2024-influence-of-wobbling-tryptophan-and-mutations-on-pet-degradation-explored-by-qm-mm-free-energy.pdf .
- **R60.** Zheng, M.; Chen, J.; Dong, W.; Wei, R.; Chen, J.; Tang, X.;
  Zhang, Q.; Wang, Q.; Wang, W.; Wang, G.; Li, Y.
  *Dynamic shift of internal electric field accelerates enzymatic polyethylene
  terephthalate depolymerization*. Commun. Chem. **2026**, 9, 85.
  https://doi.org/10.1038/s42004-026-01888-w .
  QM/MM US example using 30 ps/window, the last 15 ps for WHAM, 100 bootstrap
  data sets, profile-convergence checks and committor analysis. Protocol example,
  not a universal equilibration or sampling-time recommendation for Nyl12.

- **R63.** Politis, D. N.; Romano, J. P. *A Circular Block-Resampling Procedure
  for Stationary Data*. In *Exploring the Limits of Bootstrap*; LePage, R.;
  Billard, L., Eds.; Wiley, 1992; pp 263–270. Primary preprint: Stanford
  Department of Statistics, report EFS NSF 370, March 1991:
  https://statistics.stanford.edu/technical-reports/circular-block-resampling-procedure-stationary-data .
  Source for fixed-length circular block resampling. This is distinct from
  the authors' geometric-length stationary bootstrap (JASA 1994). The
  stationary/weak-dependence assumptions do not hold automatically for short,
  drifting QM/MM trajectories; the resampling scheme cannot repair missing
  conformations or initialization dependence.


- **R64.** Dunbrack, R. L. *Rēs ipSAE adhūc loquuntur: What the directionality
  of predicted aligned error in cofolding models can tell us about scoring
  diverse molecular complexes*. bioRxiv **2026**, posted 3 November, version 1.
  https://doi.org/10.64898/2026.11.02.763449 . Preprint, not peer reviewed.
  Ligand-aware directional ipSAE: align on protein and score ligand heavy atoms
  with fixed d0 = 2 Å. The supplied full text was inaccessible during this
  review; public author-abstract metadata and the primary software were
  examined. Do not describe full benchmark details as independently reviewed.
  Code: https://github.com/DunbrackLab/IPSAE2 ; reviewed commit
  fe560f4c2461eeb91aebded61ac8c1e83b264f60 (alpha 0.1.11). Nyl12/PA6
  scoring and thresholds remain unqualified.
