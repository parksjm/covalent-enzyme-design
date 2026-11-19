# Modular polyamide residue parameterization plan

Status: design proposal, 2026-10-28. No modular production parameters are yet
certified.

## Decision and scope

Build PA6 as a sequence of Amber residues rather than a separate whole-ligand
template for every oligomer length. Start with PA6 only. Extend the same schema
to poly(MAXA) only after its monomer structures, residue orientation, cleavage
bond and intended terminal chemistry are declared explicitly.

The natural PA6 residue is the repeat unit
`NH-(CH2)5-C(=O)`. Cut each inter-repeat amide C--N bond, as for peptide
residues:

- the head connection atom is the repeat nitrogen;
- the tail connection atom is the repeat carbonyl carbon;
- LEaP connects the preceding residue's tail C to the following residue's
  head N;
- the internal repeat has formal charge zero.

This boundary preserves the chemically important amide within the same bonded
parameter family while allowing any oligomer length to be assembled.

## Proposed residue library

Names are provisional until the implementation preflight confirms Amber/PDB
compatibility.

| Template | Chemistry | Charge | Use |
|---|---|---:|---|
| `P6I` | internal `NH-(CH2)5-C(=O)` | 0 | polymer repeat |
| `P6N` | N-terminal `NH3+-(CH2)5-C(=O)` | +1 | charged amino terminus |
| `P6C` | C-terminal `NH-(CH2)5-COO-` | -1 | charged carboxyl terminus |
| `P6U` | N-terminal `NH2-(CH2)5-C(=O)` | 0 | optional neutral amino terminus |
| `P6H` | C-terminal `NH-(CH2)5-COOH` | 0 | optional neutral acid terminus |
| `P6A` | acetyl-like N cap | 0 | mimic continuation before a modeled segment |
| `P6M` | N-methylamide-like C cap | 0 | mimic continuation after a modeled segment |

The neutral caps must be derived with the nylon library. Do not silently reuse
protein ACE/NME charges or atom types. Terminal protonation is a pH-dependent
model choice; the library should support charged and neutral alternatives
rather than asserting that all nylon termini are always charged.

For a soluble PA6 L4 zwitterion, the intended sequence is
`P6N-P6I-P6I-P6C`, total charge zero. A capped internal four-repeat segment is
`P6A-P6I-P6I-P6I-P6I-P6M`, also total charge zero. These represent different
physical systems and must not be interchanged.

## Charges and bonded parameters

Retain GAFF2-compatible atom types initially. The active QM torsion and
improper scans remain applicable because they test local repeat and terminal
chemistry. Do not interrupt or replace those calculations.

Do not derive repeat charges from an isolated monomer. Fixed charges depend on
the oligomer environment, and independently fitted monomers can produce poor
amide electrostatics or a length-dependent accumulated charge error. Instead:

1. Optimize multiple conformers of capped PA6 trimers or pentamers and charged
   terminal oligomers at the declared geometry level.
2. Compute gas-phase HF/6-31G(d) Merz--Kollman ESPs, consistent with the current
   RESP procedure.
3. Perform a simultaneous multiconformer RESP fit with graph-equivalent atoms
   constrained across conformers and chemically equivalent internal repeats.
4. Constrain every `P6I`, `P6A` and `P6M` template to net charge zero and the
   charged terminal templates to +1 or -1. Fit neutral terminal variants
   separately at net zero.
5. Inspect ESP residuals before accepting the integer-residue constraints. Do
   not hide a poor fit by distributing a residual charge arbitrarily.

The current 79-atom full-L4 RESP/GAFF2 model is the reference for electrostatic
and energetic comparison. Its per-atom charges should not simply be copied
into templates because the four repeats are context-dependent and were not
fit under residue-integer constraints.

## Boltz and AlphaFold 3 mapping

Predicted ligand atom order and atom names are not authoritative. Mapping must
be graph-based:

1. Build the canonical heavy-atom graph for the declared PA6 length and
   terminal/cap state from the residue templates.
2. Read the predicted ligand graph from mmCIF chemical-component atom and bond
   records.
3. Match graphs by element, formal charge, bond order and degree. Anchor chain
   orientation with the chemically distinct N and C termini.
4. Walk monotonically from the N terminus to the C terminus, assign repeat
   indices, then rename and split atoms into Amber residues.
5. Regenerate hydrogens from the Amber templates; do not depend on predicted
   hydrogen names or ordering.
6. Write an immutable JSON atom-map manifest containing source atom IDs,
   canonical residue/atom names, coordinates, connectivity and hashes.
7. Rebuild the Amber topology and require a graph round trip before relaxation.

For neutral, symmetrically capped chains, graph orientation may be ambiguous.
In that case orientation must come from the input/query atom map; it must not
be guessed from coordinates.

## Acceptance gates

The modular library is not production-ready until all of the following pass:

1. **Assembly:** lengths 1--8 build with expected atom counts, bonds, residue
   order and exact total charges for every declared terminal state.
2. **Mapping:** representative Boltz and AF3 CIF files round-trip to the same
   canonical heavy-atom graph without manual atom-name edits.
3. **Electrostatics:** constrained multiconformer RESP residuals are acceptable,
   and modular L2/L4/L6 electrostatic potentials and dipoles are compared with
   whole-oligomer QM references.
4. **Bonded terms:** the active proper/improper QM scans pass or flagged terms
   are refit and revalidated.
5. **Whole-L4 comparison:** the modular `P6N-P6I-P6I-P6C` topology is compared
   with the current whole-L4 candidate at identical coordinates for energy,
   force, geometry and short explicit-solvent MD behavior.
6. **Products:** separate L2 product topologies have the intended newly formed
   N and C termini and declared protonation states.
7. **Enzyme states:** modular noncovalent substrate templates remain distinct
   from TI1, acyl-enzyme and other covalently linked state topologies.

## Implementation sequence

Finish the current QM bonded-term benchmark first. Then implement, in order:

1. a chemistry-only residue schema and canonical PA6 graph builder;
2. graph-to-residue mapping with round-trip tests;
3. RESP constraint generation for capped and charged oligomer training sets;
4. Amber OFF/LIB and FRCMOD generation;
5. L2/L4/L6 assembly and comparison gates;
6. integration with Boltz/AF3 postprocessing and enzyme relaxation.

This work should replace per-length parameter duplication only after it matches
or improves upon the current validated whole-L4 reference.
