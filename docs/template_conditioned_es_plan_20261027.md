# Template-conditioned Nyl12 E-S pilot

## Decision

Do not interpret failure of one fully unrestrained, whole-system Model 23
optimization as failure of Nyl12. It shows that this seed relaxes into a
non-native Thr233 basin under the present global optimization protocol. Test a
small crystal-template-conditioned ensemble before considering an enzyme pivot.

## Pilot scope

1. Use `cif/n12c.cif` as the protein template and preserve the mature symmetric
   tetramer, catalytic N-terminal Thr233 residues, and canonical mapping.
2. Model one PA6 L4-to-L2 substrate at a declared active site using the minimum
   contact restraints required to set register and site. Generate independent
   seeds with both Boltz and the validated contact-restraint AF3 implementation
   when that implementation is ready. Contact restraints are generation aids,
   not acceptance evidence.
3. Before relaxation, reject molecular-graph errors, clashes, noncanonical
   numbering, wrong termini/protonation, and models that fail the established
   PA6 Tier A geometry. Rank retained seeds by crystal-aligned Thr233 heavy-atom
   RMSD and N/OG1 displacement, Thr chi1 and psi differences, and the
   CB–OG1–C54 angle. Report that angle without a hard cutoff initially;
   105–120 degrees is an exploratory range supported by the 116.95-degree
   hASNaseIII TS2 value, not a validated Nyl12 gate.
4. Perform staged explicit-solvent MM preparation: relax solvent and hydrogens,
   then the local pocket, while initially restraining the template backbone and
   ligand register. Release all ligand/contact restraints before acceptance.
   Use short independent relaxation/MD replicas and rescore distributions;
   do not infer stability from a single minimized endpoint.
5. For DFTB3/MM refinement, optimize the QM region and a declared inner MM
   shell while holding only the distant environment fixed or weakly restrained.
   Test at least two shell radii. Catalytic atoms, substrate, oxyanion donors,
   Asp272/Asp274, and nearby waters must be active. This is distinct from the
   rejected approach that restrained every non-QM atom and from the inefficient
   current attempt to optimize all bulk-water coordinates to a strict maximum
   gradient.
6. Release any local catalytic restraints and require retained Tier A/B
   geometry, acceptable Thr233 orientation, finite coordinates, and visual
   trajectory inspection. Promote only geometries reproduced across independent
   seeds; preserve failures in the denominator.

## Stop/go criterion

Begin with a small pilot rather than a production ensemble. If neither Boltz
nor AF3 produces multiple template-conditioned models that retain native-like
Thr233 and near-attack geometry after contact-restraint release and local
QM/MM refinement, pause Nyl12 redesign and run the same minimal pilot on NylB.
NylB is a useful lower-complexity control, not yet a replacement target.
