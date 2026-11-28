# RFD3 backbone symmetry qualification

Updated 27 November 2026. RFD3 remodeling must preserve the assembly symmetry
of the four coupled mature alpha/beta repeats: A/E, B/F, C/G and D/H.
LigandMPNN sequence ties alone do not meet this requirement.

## Audit and current hold

The [retained audit](../design/nyl12_relay_ts_ensemble_20261126/backbone_qualification/symmetry_qualification_20261127/symmetry_contract_audit.json)
found no symmetry specification, the default sampler and chain-A-only mobility
in the existing local RFD3 parents. Their earlier chemistry/graph acceptance
does not certify backbone symmetry. Preserve these outputs and let existing
frozen sequence-generation/scoring jobs finish as diagnostics. Do not advance
their RFD3-generated candidate lineages. Original MD-backbone lineages remain
eligible for review; cross-scores on unqualified remodels remain diagnostic.

Crystal CA correspondences give three nonidentity rotations near 180 degrees,
supporting approximate D2 noncrystallographic symmetry rather than C4. These
fits do not establish exact operator closure or qualify a native RFD3 run.

## Required qualification before another remodeled parent advances

1. Qualify D2 operators and the coupled alpha/beta correspondence in the actual
   input coordinate frame; do not treat eight mature chains as eight identical
   subunits. Use equivalent mobile-residue masks in all four repeats.
2. Use the installed native symmetry specification and symmetry sampler.
   Native Gaussian perturbations are independent; denoised coordinates are
   projected. Do not describe the implementation as tied raw noise.
3. Resolve compatibility between protein symmetry and the fixed single-site
   TS-like context before execution. Do not silently duplicate substrate atoms,
   alter the proton inventory or move the fixed reaction geometry.
4. Adapt and qualify the export coordinate transform for this sampler: the
   current mobile-center-of-mass exporter explicitly rejects symmetry inputs.
   Preserve canonical maps, identities, typed ligand graph and clash checks.
5. Declare backbone residual tolerances before the run; independently measure
   final symmetry residuals and fixed-context retention. Freeze the actual
   specification, sampler arguments, masks, operators, source and output hashes.

The admission gate in `rfd3_symmetry_acceptance.py`, called by
`relay_broad_runtime.verify_parent`, requires a new hash-bound
`backbone_symmetry_acceptance` JSON receipt. It checks schema version 1,
`accepted: true`, `point_group: D2`, `sampler_kind: symmetry`, the four
`coupled_repeats`, `equivalent_mobile_masks`, `operator_closure_accepted` and
`fixed_ts_context_retained`. The receipt must bind the exact parent PDB and
canonical mapping digests. It must report finite nonnegative
`backbone_rmsd_A` and `backbone_max_deviation_A`, each within its positive
predeclared limit (`<metric>_limit`). Synthetic test tolerances are not a
scientifically approved tolerance for a real run.

This gate consumes independent evidence; it does not produce operators or
measure coordinates. No existing parent has been issued a symmetry acceptance
receipt. Missing evidence fails explicitly. Include the new imported module
in every future frozen source release; never modify a running release.

## Experimental testing remains a separate decision

QM/MM PMFs and verified PES transition states are not prerequisites for PI
experimental nomination. Sequence, maturation/assembly, structural and available
matched MD evidence can support a bounded WT/design experiment. Optional
matched QM/MM calculations may be run before testing if resources permit, or
afterward to explain an observed activity difference. Neither geometry nor a
sequence score establishes a lower activation free energy.
