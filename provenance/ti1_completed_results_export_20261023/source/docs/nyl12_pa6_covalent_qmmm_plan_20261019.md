# Nyl12–PA6 L4z → L2: covalent-state and QM/MM plan

**Status:** implementation plan, updated 2026-10-20. Initial deliverable: audited
wild-type first tetrahedral-intermediate (TI1) models at one active site of the
mature Nyl12 assembly, for PA6 L4 → L2 and poly(X,AMA) → MAXA. This is a
follow-up methods project, separate from the commercial-polyamide manuscript.
No new production dependencies or predictor weights are installed by this plan.

Current execution status is in `nyl12_covalent_resume_20261019.md`. Both
substrates now have covalent AF3 predictions and TI1 topologies. Water placement
and geometry refinement remain unresolved; this plan does not establish that
either model is optimized or suitable for classical MD.

## What the new papers change

- Wu *et al.* (2026, AlphaProtein Novo; [preprint](https://doi.org/10.64898/2026.10.01.756017)) report that mechanism-specific geometry, multiple AF3 seeds, multiple catalytic states, and the fraction of LigandMPNN resequences passing filters improve experimental hit enrichment. Their DEHP esterase campaign explicitly modeled a covalent tetrahedral intermediate and screened several states. The authors also state that their predictors do not yet model catalytic physics adequately. This supports a **multi-state, multi-seed evaluation**, not interpreting a predicted pose or confidence as a reaction barrier. The paper is a preprint as of this plan.
- Butcher *et al.* ([Science 2026](https://doi.org/10.1126/science.aeh2982), including the supplied supporting material) show atomic motif scaffolding with RFdiffusion3 (RFD3), followed by LigandMPNN/sequence design and independent structure prediction. They experimentally validate a carbamylated-lysine, binuclear-metal phosphotriesterase and a zinc protease. Their supporting methods derive a theozyme from a QM transition-state calculation and, in a second round, refine a transition-state geometry with a machine-learning-force-field NEB calculation. Their reported NEB barrier is an **electronic** barrier without zero-point, thermal, or entropic corrections. This is a useful precedent for constructing motifs, not a validated Nyl12 reaction profile.
- Existing Nyl12 work already has mature alpha A–D (1–232) and beta E–H (1–96) chains, tested PA6 registers, a focused 21-position LigandMPNN mask, and matched AF3/Boltz controls. The selected Boltz L4z → L2 scissile atoms are **C54/O49/N59**. The new covalent AF3 custom CCD preserves these names exactly; **C12/O4/N3** belongs to an earlier mapped AF3 representation. Never transfer names without a graph isomorphism or apply an affinity-mode naming offset to the new covalent models. Existing bound-start MD and whole-ligand scores did not robustly explain variant activity. See [focused workflow](../manuscript_nyl12_wt_vs_mutants/NYL12_FOCUSED_REDESIGN_WORKFLOW_20261016.md) and [current interpretation](../manuscript_nyl12_wt_vs_mutants/NYL12_PA6_PA66_FINAL_INTERPRETATION_AND_NEXT_PLAN_20260928.md).
- Dishner *et al.* ([JACS 2026](https://doi.org/10.1021/jacs.6c03466)) report that Nyl12 predominantly releases the **MAXA** tetrad from hydrolysis of **poly(X,AMA)** under their assay conditions, while Nyl50 favors XAMA. Here M = hexamethylenediamine, X = *p*-xylylenediamine, and A = adipic acid. The paper also synthesizes **poly(MAXA)** from purified MAXA and finds similar repeat ordering but a different end-group distribution. It does not establish a Nyl12 hydrolysis rate or product distribution for poly(MAXA) itself. The reported soluble product yield was below 5 mol%, so polymer accessibility and product recovery matter alongside catalytic geometry.

## Correct the chemical target first

The proposed Thr233 OG1 attack and relay of its proton through water to the
neutral beta-chain N-terminal amine describes **formation of the first
tetrahedral intermediate**, if the scissile C–N bond remains intact. It does
not by itself produce a stable acyl enzyme. For the acyl enzyme, carbonyl C=O
must reform and the scissile C–N bond must cleave: the retained PA6 acyl
fragment is an **ester of beta Thr1 OG1**, and the departing PA6 fragment is a
separate amine product. Both product fragments, including their terminal
hydrogens, must be accounted for. Model the post-cleavage pocket both with and
without the departing product if release is relevant.

Define these separate, versioned graphs and proton inventories:

| State | Required reactive connectivity | Proton question |
| --- | --- | --- |
| ES | Intact substrate amide; Thr–OH; free beta Thr1 alpha amine | User-specified neutral NH2 starting state; water occupancy |
| TI1 | Thr OG1–scissile C **and** scissile C–N; tetrahedral oxyanion | User-specified NH3+ beta amine and neutral relay water |
| Acyl | Thr OG1–scissile C; C=O restored; scissile C–N broken | Which proton supplies the leaving amine, and is beta Thr1 NH2 or NH3+ afterward? |
| TI2/products, later | Water-derived O attacks ester; eventual Thr–OH regeneration | Separate deacylation proton relay and oxygen provenance |

These protonation states define the requested models. Whether the one-water
path connects them, and whether TI1 is a local minimum, remain questions for
the calculations. Track the proton inventory without asserting that an
individual water-mediated proton retains its identity from OG1 to the amine.
Additional microstates can be compared later with integer charge and atom
accounting. Avoid importing the old ES Thr-alkoxide topology into the acyl
state. Do not graft a serine–histidine–aspartate mechanism onto this N-terminal
threonine hydrolase.

For example, a neutral Thr1 NH2/OG1H starting state could give a TI1 with
Thr1 NH3+ and a substrate oxyanion after water-mediated proton transfer. A
subsequent proton donation to the scissile amide N could yield a neutral
departing amine and restore Thr1 NH2. This is one charge-balanced candidate
path, not an assigned mechanism or proof that one water suffices.

## Product-constrained poly(X,AMA) graph

The user narrowed the second substrate target to **poly(X,AMA)** and its
assay-observed **MAXA** product. The implemented explicit oligomer is
`M-A-X-A|M-A-X-A`; the cut is at `A|M`, with the A carbonyl attacked by
Thr1 OG1 and the M amide nitrogen departing. Graph cleavage and hydrolytic
termination give two MAXA products, each neutral formula C26H42N4O5 and exact
mass 490.31552 Da. Keep explicit terminal chemistry, stereochemistry, and atom
mapping in the molecular graph. A tetrad name alone is insufficient.

The competing XAMA product and assays on separately synthesized poly(MAXA)
are outside the current requested design target. PA6 likewise has a hard
product gate: L4 cleavage must yield the observed L2 product, neutral
C12H24N2O3, exact mass 244.17869 Da. The TI1 itself retains the scissile C–N
bond; product validation is a graph-derived prediction of subsequent cleavage,
not a claim that the intermediate already contains the released products.

The X aromatic ring creates a distinct packing and solvation problem from
PA6's uniform aliphatic repeat. Parameterize M, X, and A with explicit amide
bond orders and consistent GAFF charges and torsions. Keep monad and repeat
indices stable across conformers and predictors. Compare the paper's
poly(X,AMA) end-group distribution with poly(MAXA)'s heterotelechelic ends;
use internal, capped fragments to study local catalysis, then vary length and
realistic termini to measure end effects. Do not transfer a PA6 atom map,
force-field library, or register score to this copolyamide. A direct purified
Nyl12/poly(MAXA) assay would be needed to validate predictions for that
synthesized polymer.

## Large, flexible oligomers: representation and sampling

1. **Separate access from chemistry.** A solvated oligomer measures a pocket
   hypothesis. A solid or semicrystalline polyamide also requires chain
   exposure, adsorption, and extraction from polymer–polymer contacts. Record
   substrate morphology and assay temperature; compare enzyme variants in
   matched experiments before attributing a rate change to transition-state
   stabilization.
2. **Use a length and end-group ladder.** Start with a short internal fragment
   that includes the scissile bond and its neighboring repeats, then extend
   both tails until local active-site geometry, water access, and QM/MM
   reaction-core results stop changing materially. Include at least one
   deliberately different terminal chemistry. PA6 L4z and capped controls
   already provide an initial ladder; `MAXA|MAXA` starts the copolyamide
   ladder. Make the chemical graph, repeat direction, and scissile register
   explicit for every length.
3. **Sample conformations before and after binding.** Generate diverse
   solution conformers, including intrachain hydrogen-bonded states, then
   sample multiple binding registers, chain directions, and tail poses with
   independent predictor seeds. Cluster by the reacting core **and** tail
   arrangement, retain failed/unbound predictions, and avoid choosing only
   the best whole-ligand score. Inspect intersubunit threading, chain
   crossings, steric clashes, and unsatisfied buried amides. Structure
   predictors trained mainly on compact ligands may fail on extended
   polyamides; use short-fragment predictions as pose hypotheses and graft
   longer tails from the conformer ensemble only with explicit graph and
   clash checks.
4. **Relax tails with validated full MM parameters.** Minimize and run
   independent explicit-solvent MD from several distinct pose clusters;
   monitor core geometry, water occupancy, tail contacts, and excursions of
   terminal groups. Use enhanced conformational sampling only if ordinary
   replicas remain trapped and a specific slow coordinate can be defined.
   Compare states at the same length, caps, salt, temperature, and occupancy.
   Whole-ligand contact counts and `interface_delta` grow with oligomer size
   and must not rank unequal-length substrates directly.
5. **Keep QM local, retain the environment.** The long tails normally remain
   in the MM region, but their electrostatics, sterics, and conformations
   enter QM/MM through the solvated assembly. Reoptimize the same reacting
   core over representative tail/conformer clusters, and report the spread
   rather than one minimized structure. If tail length changes the local
   result, extend the model or sample more of the relevant conformations.

Add a versioned substrate-builder interface and a conformer/register ledger
to the proposed scripts. The builder must emit molecular graph, atom names,
monad/repeat labels, terminal groups, product fragments, charges, and an
unambiguous scissile C/O/N map; the ledger must connect each model and MD
replica to its length, cap, register, direction, and starting conformer.

## Staged, scripted workflow

1. **Freeze reference chemistry and coordinates.** Select one audited WT
   PA6 L4z → L2 pose from the existing three-structure focused set, retaining
   the full mature assembly and one substrate. Record source hash, chain map,
   atom map, terminal/capping chemistry, bond orders, explicit hydrogens,
   protonation, net charges, scissile register, and observed oxyanion contacts.
   Inspect Phe113 backbone NH and Asn186 ND2. Make an atom ledger that maps
   every substrate atom into the acyl or departing fragment, including newly
   added protons and water-derived atoms. Reject any missing or duplicated atom.
2. **Build covalent candidates from the audited graph.** Chemically edit the
   ES graph to make TI1 and the cleaved acyl state, then generate several local
   conformers and side-chain/water arrangements around the same full protein
   pose. A deterministic graph edit is the reproducible baseline. Preserve the
   physical alpha/beta chain break and the free beta Thr1 alpha-amino group.
   Reject invalid valence, stereochemistry, charge, product identity, or
   protein–substrate clashes before prediction.
3. **Small predictor comparison, only after the graph gate.** AF3 with
   `userCCD` and `bondedAtomPairs` is the preferred **covalent pose sampler**:
   its official interface supports named custom ligand atoms and inter-entity
   bonds ([input specification](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)). AF3-LA weights are specifically tuned for leaving-atom handling ([official README](https://github.com/google-deepmind/alphafold3/blob/main/README.md)), but are absent from the local install and are optional for the already-cleaved acyl graph. Use installed standard AF3 first, with a few fixed seeds, and compare against graph-built poses. Boltz-2 also declares a protein–CCD ligand `bond` constraint ([official schema](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)); use it as a matched secondary arm only if the local Lightning Boltz parser round-trips the custom CCD and linkage. Neither predictor determines proton transfer, charge, or a transition state. Avoid changing either upstream predictor until a minimal failing input shows an actual feature gap; any patch belongs in a separate, pinned fork with a regression case.
4. **Parameterize TI1 first.** Use ff19SB for unchanged protein, the
   project's GAFF treatment for substrate atoms, and a modified **N-terminal** Thr-TI1
   residue or equivalent cross-residue topology. The beta residue must still
   connect to beta residue 2. Generate Amber-compatible partial charges for
   each fixed protonation state using a documented ESP/RESP-like calculation
   on a capped model; preserve consistent ff19SB boundary charges and audit
   total integer charge after cap removal. Derive or validate atom types,
   Lennard-Jones terms, and all cross-junction bonds, angles, torsions,
   impropers, exclusions, and 1–4 interactions. Check the resulting MM minimum
   and key reactive-junction torsions against quantum reference data. Do not infer bond
   orders from a ligand PDB or accept `parmchk` completeness as proof of
   parameter quality.
5. **QM/MM geometry optimization.** Use the installed AmberTools26 `sander`
   with DFTB3 and its 3ob-3-1 Slater–Koster set for a bounded pilot. Keep the
   whole solvated assembly as MM environment, with electrostatic embedding and
   a QM region containing beta Thr1 free amine and nucleophilic side chain,
   scissile C/O/N with adjacent substrate groups, and candidate relay water(s).
   Put QM/MM cuts on ordinary nonreactive sigma bonds away from the ester,
   amide, and proton-transfer path; audit link atoms and boundary charges.
   Start with a restrained environment, then release local restraints in
   stages. Include oxyanion donors in a QM-region sensitivity calculation if
   their polarization affects the result. In QM/MM, internal QM bond energy is
   computed quantum mechanically, but Amber still needs a valid complete
   topology, MM charges/LJ and exclusions, and correct terms at QM/MM
   boundaries. **Only charges and LJ parameters is therefore insufficient as
   a general input recipe.** Record convergence, bond/angle geometry, water
   network, local protein shifts, and whether TI1 collapses on relaxation.
6. **Validate the electronic model.** DFTB3 is a practical geometry search
   method, not a trusted amide-hydrolysis barrier by default. Benchmark small
   capped active-site clusters and representative stationary points against
   a higher-level dispersion-aware DFT calculation, including protonation and
   water alternatives, before making barrier claims. Compare geometries and
   relative energies; test QM-region and MM-embedding sensitivity. A saddle
   candidate needs a reaction-coordinate/frequency check appropriate to the
   constrained system. If a TI minimum disappears under unbiased relaxation,
   report that outcome rather than holding it with arbitrary restraints.
7. **MD only when it answers a structural question.** After full classical
   parameter validation, run staged minimization, heating, restraint release,
   then independent explicit-solvent replicas for a WT acyl shortlist. Monitor
   ester geometry, product escape, local waters, protein assembly, and PA6-tail
   conformations. Existing bound-start controls and matched settings provide
   comparison. Classical MD relaxes tails; it cannot sample bond making or
   report turnover. Tail mobility alone is not failure if the reacting core
   remains well placed.
8. **Reaction profile after state validation.** Optimize ES, TI1, acyl,
   deacylation TI2, and products with the **same atom inventory and compatible
   QM/MM setup**. For acylation, sample at least Thr O–C formation, scissile
   C–N cleavage, and relevant proton-transfer coordinates; for deacylation,
   include water O–C approach and proton transfers. Begin with restrained
   scans and candidate saddle refinement. NEB may be a useful path finder once
   chemically correct endpoints exist, but its minimum-energy path is not a
   solution free-energy profile. For rate-limiting-step or mutant claims, use
   replicated QM/MM umbrella sampling or another validated free-energy method
   after the DFTB3 benchmark, with explicit uncertainty. Never compare raw
   energies of states with different atom counts, protonation, or solvation.

## Reproducible implementation and stop gates

Keep adapters small and separate: `substrates.py` (repeat chemistry, length,
caps, and register), `chemistry.py` (state graph/atom ledger),
`build_states.py` (coordinate editing), `predict_covalent.py` (AF3/Boltz input
and output maps), `amber_topology.py` (library/frcmod/tleap audit),
`qmmm.py` (Amber input builder and parsed diagnostics), `ensemble.py`
(conformer/register provenance), `profile.py` (reaction coordinates), and
`report.py` (state comparison). Put chemistry and run settings
in versioned YAML, and record software/weight hashes, seeds, exact commands,
input/output hashes, atom maps, charges, QM selections, restraints, and failures
in a machine-readable manifest. These names are proposed interfaces, not files
already present. Use `uv` for any Python environment changes, with confirmation
before adding a production dependency.

**Gate A, chemistry:** exact graph and hydrogen ledger, product mapping,
integer charges, correct mature chain boundaries, and the intended PA6 L2 or
copolyamide MAXA/XAMA-producing bond.
**Gate B, classical topology:** zero missing/ambiguous parameters, correct
ff19SB/GAFF boundary, integer total charge, valid ester geometry and basic
torsional behavior. **Gate C, QM/MM:** reproducible convergence from more than
one start, no unintended chemistry, and robustness to a reasonable QM-region
and protonation/water alternative. **Gate D, reaction claim:** benchmarked
electronic method and, for kinetics, a sampled free-energy calculation with
uncertainty. Failure at a gate stops downstream interpretation.

## Design after the chemistry gate

Start from the already prepared, symmetry-tied **focused LigandMPNN** Nyl12
mask and the protected maturation/catalytic positions, because preserving the
known alpha/beta fold and four-site assembly is central to this redesign.
Evaluate apo, ES, and audited covalent states over seeds and sequences, using
state-specific reacting-core geometry and water access rather than one global
confidence score. Use RFD3 for a later **local backbone/loop diversification**
branch if focused sequence sampling lacks viable geometry; test whether its
partial-diffusion controls preserve the processed chains, oligomer, and all
protected catalytic atoms. A full de novo scaffold campaign would answer a
different question and consume far more experimental screening. Select a
diverse small panel for purified-enzyme PA6 assays; the existing lysate data
alone cannot establish an intrinsic catalytic-rate gain.

## GPT model for implementation

Use **GPT-6 Astra at high reasoning effort** for the initial state-graph,
charge/proton inventory, covalent topology, QM/MM boundary, and reaction-path
implementation and review. Those choices couple chemistry and code, and an
internally consistent but wrong graph could invalidate every downstream run.
Once those interfaces and acceptance gates are frozen, **GPT-6.1 Sol at medium
or high effort** is suitable for routine adapters, manifests, batch scripts,
and analysis. If one model must own the whole implementation, choose Astra
high; increase to extra-high for disputed mechanism or barrier interpretation.
This division follows [official OpenAI model-selection guidance](https://developers.openai.com/api/docs/guides/model-selection), not a benchmark of either model's enzyme-chemistry accuracy. Chemistry review and experimental validation remain required.
