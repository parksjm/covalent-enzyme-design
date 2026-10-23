# Redesign backbone, symmetry and substrate access — 2026-10-22

## User requirements and current decision

The user explicitly requested consideration of whether RFD3 is necessary,
which backbone segments would change, whether WT lengths should be retained,
and how insertion/deletion variants would be postprocessed. They emphasized
that access loops must be able to adopt more open conformations, and that
both RFD3 and LigandMPNN must respect the symmetric tetramer. The user has
working LigandMPNN symmetry examples and will point to them when needed.
Do not ask for those examples as a prerequisite to independent preparation.

Recommended first campaign: focused LigandMPNN sequence redesign on retained
WT backbones, followed by independent AF3 predictions and physical relaxation,
with WT chain lengths retained. RFD3 is optional, not a prerequisite for a
working sequence-redesign pipeline. New backbone generation is a later,
controlled comparison if sequence redesign cannot provide the required
geometry/accessibility. Indels are not enabled in the first campaign.
This is a staged project recommendation; no new design campaign has launched.

| Route | Backbone and length policy | Role and limitation |
| --- | --- | --- |
| LigandMPNN -> AF3 -> relaxation | WT backbone at sequence-design input; WT alpha/beta lengths | First baseline; simplest correspondence and mature-fold preservation. AF3 and relaxation can move the backbone. Fixed input does not guarantee catalytic geometry or flexible access. |
| Local RFD3 -> LigandMPNN -> AF3/relaxation | Selected short loops move; exact WT segment lengths initially | Later comparison when backbone placement limits sequence-only designs. Retain the core, maturation architecture and essential interfaces. |
| RFD3 with loop indels | Explicit loop length ranges and matched changes in each equivalent copy | Separate exploratory branch only after full variant-mapping/postprocessing validation. Longer/shorter loops are not intrinsically more flexible or accessible. |

Compare routes with common WT controls, ligand states, screening criteria and
reported compute budgets. Preserve the parent design/backbone/sequence lineage.
AF3 confidence and a favorable closed pose are not activity measurements.

## Where RFD3 would be used

Do not regenerate the entire enzyme fold or arbitrarily select loops from a
sequence-number interval. Identify short pocket/entrance loops and second-shell
segments from the accepted WT assembly and state-specific geometry, including
substrate/product entry/exit routes. Their actual residue ranges have not yet
been selected; the validated structure and access analysis must supply them.

Preserve catalytic Thr identity, its maturation context and the processed beta
N terminus. Protect the core and required oligomer contacts. Use a compact,
validated catalytic atomic motif/ensemble rather than fixing an entire access
loop to one closed structure. If a selected loop carries a catalytic anchor,
restrict changes to compatible flanks and check all relevant states. Oxyanion
and proton-relay requirements remain explicit; a TS is not represented by
inventing an ordinary stable covalent adduct/bond order.

RFD3 documents atom-specific fixed-coordinate masks, independent sequence
masks, exact or ranged generated contigs and ligand exposure/burial hints [D1].
An exact-length generated segment between retained context regions is the
initial local-backbone concept. Exact loop/anchor/chain-break specifications
must be checked on this assembly before production. The installed RFD3
repository and relevant input examples are retained by source identity; this
is not a claim that this combined Nyl12 configuration is already validated.

Selected distal ligand atoms can potentially receive exposure conditioning
where consistent with the desired substrate register and terminal recognition.
This is a static geometric hint, not a continuous substrate channel or an
opening-rate constraint. Do not force all charged termini to remain solvent
exposed when terminal recognition is itself a design goal. No identified
RFD3 control guarantees a flexible gate or a low opening barrier.

## Length and residue-identity contract

Initial lengths are WT mature alpha 232 and beta 96 residues for each repeat.
The precursor-to-mature correspondence includes WT precursor Thr233 = beta1.
Retain chain breaks and maturation boundaries. A fixed-length design can still
change residue identities and backbone coordinates; that must not erase its
reference correspondence.

Current WT-specific code is not indel-ready. Examples include:

- `af3_postprocess.py` checks alpha/beta counts 232/96 explicitly.
- `inspection.py` hard-codes WT Amber residue blocks and Thr929/ligand1313.
- `local_repair_model.py` fixes dimer selections and a QM mask from WT indices.
- `launch_nyl12_af3_water_test.py` is a specific WT diagnostic launcher, not a
  generic variant predictor.

Keep these checks as guards for their historical workflows; do not silently
remove them or accept a length mismatch. Future variant-capable entry points
must consume an explicit per-design map and regenerate topology/QM selections.

Before enabling any indel branch, each residue must have an unambiguous mapping
between: design identity, chain family, biological copy, WT precursor/mature
position (nullable for insertions), designed sequence position, structure
label/auth residue identifiers including insertion codes, and Amber residue/
atom indices assigned after topology construction. Give inserted residues
explicit parent anchors and stable insertion identities. Record deletions.
Use generator correspondence and declared anchors where available; ambiguous
sequence alignment is an error rather than a silent positional fallback.

Generate catalytic selections, symmetry ties, contact restraints, maturation
checks, geometry analysis, QM masks, visualizations and mutation reports from
that map. Do not retain numeric selections such as Thr233 or Amber929 without
resolving their mapped identities. Prefer mmCIF as the identity-preserving
structure record, with PDB exports accompanied by correspondence tables.

Required future tests include identity maps, insertions/deletions before a
catalytic residue, insertion codes and reordered chains, missing anchors,
ambiguous repeats, inconsistent copy lengths/indel patterns and round trips
through AF3 and Amber. No current exporter is certified for these variants.
A same-length RFD3 output also requires checking generator residue order.

## Symmetric assembly and sequence ties

The mature enzyme is four alpha/beta repeats, (alpha-beta)4, represented by
alpha chains A/B/C/D and beta chains E/F/G/H, with sites AE/BF/CG/DH. It is
not eight identical protein chains. Tie homologous alpha positions across
A/B/C/D and homologous beta positions across E/F/G/H independently. Never
use an all-chain homooligomer shortcut that ties alpha to beta.
LigandMPNN explicitly supports arbitrary tied-position groups and weights [D2].
Use mapped identities rather than chain-letter assumptions for generated or
renumbered structures, and validate any user-provided examples against the
chosen source sequence/assembly.

Do not equate tetramer with C4. A scripted Kabsch audit of the retained
AF3-derived WT scaffold measured reference-to-copy rotation angles
179.9874, 179.9805 and 179.9291 degrees, CA RMSDs 0.1114, 0.0511 and 0.1157 A,
and nearly orthogonal rotation axes (largest off-diagonal absolute dot product
0.00157). These are D2-like geometric diagnostics, not a biological-assembly
assignment from a crystal. Verify the experimental assembly, copy pairing and
symmetry operators before choosing a production RFD3 group. D2 has four group
operations; it should not be confused with a two-copy oligomer.

RFD3 documents C/D symmetry and pre-symmetrized catalytic motifs, including
examples with ligands [D3]. The exact combination of two-chain alpha/beta
repeats, retained scaffold regions and the Nyl12 catalytic motifs still needs
an input preflight/pilot. Orient/center the assembly and motif consistently;
retain any canonicalization transform and copy map. Do not independently
generate four unrelated loops and average their coordinates afterward.

Sequence symmetry does not require identical instantaneous conformations.
Physical relaxation and accessibility tests should permit equivalent copies
to open asynchronously. Specify ligand occupancy separately from protein
symmetry: one occupied site and four occupied sites are different tests.
Retain all four site outcomes; symmetry-equivalent site recovery is one model
outcome, and sites within a model are not independent experimental replicates.

## Avoiding inaccessible or over-stabilized closed loops

The target is compatibility with both productive catalytic states and an
accessible apo/ES ensemble, not maximal flexibility everywhere. Keeping the
catalytic core organized while allowing a gate to open is a multistate problem.
Fixed length, low backbone RMSD, high confidence and low minimized energy do
not measure opening kinetics. Nor does high RMSF by itself imply useful access.

Use a staged screen:

1. Establish matched WT apo/ES and catalytic-state reference ensembles with
   explicit protonation, assembly and occupancy. Preserve available open
   structural conformations without calling prediction diversity a measured
   thermodynamic ensemble.
2. Design on representative compatible backbones and evaluate candidate
   sequences across them. LigandMPNN's multiple-PDB option is batch processing,
   not a documented joint multistate objective. Cross-backbone scoring and
   selection therefore require an explicit project wrapper; a good score on
   one closed scaffold is insufficient.
3. Independently predict apo and noncovalent ES models with multiple seeds,
   matched WT controls and no imposed active-site contact restraints or
   designed-pocket template. Preserve the separate covalent-state arm.
4. Compare entrance geometry, bottlenecks and steric accessibility for the
   actual substrate/product register. A small-probe tunnel does not prove
   passage of a PA6 oligomer or copolyamide. Check neighbouring subunits too.
5. After force-field validation, use matched independent relaxation/sampling
   from both open-compatible and catalytic conformations. Monitor gate
   distances, contacts, secondary-structure/hinge changes and ligand access,
   alongside catalytic geometry. Flag new closed-state locks for review;
   do not apply blanket Gly/Pro rules as substitutes for this analysis.
6. If short trajectories do not show opening, report it unresolved rather than
   conclude a loop cannot open. More detailed conformational sampling or
   opening free-energy calculations are a shortlist task. AF3 and ordinary
   minimization alone cannot establish an opening barrier or rate.

Retain state-specific reaction checks (ES/TI1, acyl enzyme and deacylation
compatibility), so improving access does not remove the catalytic arrangement
or over-stabilize only one intermediate. Detailed mechanism and enhanced
sampling remain gated on accepted structures/force fields.

## Immediate parallel work and provenance

Safe independent work includes sequence/symmetry manifest preparation, auditing
WT-specific index assumptions, planning matched-state/access screens and
preparing standard Amber typing/RESP/validation adapters without fitting an
unaccepted reference. No RFD3/LigandMPNN campaign, new AF3 run or physical loop
sampling has launched during this decision record. The user is independently
improving their AF3 contact-restraint fork; do not modify it here.

The repeat-geometry calculation is scripted in tmux with source hashes and
retained input under
`output_nyl12_parameter_validation_20261022/repeat_symmetry_audit_01/`.
Active journal: `output_agent_actions/20261021.jsonl`. Source review and
recommended procedures are distinct from implemented design capabilities.

## Primary software sources

- [D1] RFD3 input specification:
  https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/input.md
- [D2] LigandMPNN README and sequence-tying implementation:
  https://github.com/dauparas/LigandMPNN/blob/main/README.md and
  https://github.com/dauparas/LigandMPNN/blob/main/run.py.
  Scientific reference: https://doi.org/10.1038/s41592-025-02626-1.
- [D3] RFD3 symmetry and enzyme examples:
  https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/examples/symmetry.md
  and the sibling enzyme-design example. Installed small-molecule examples
  document exposure/burial controls. Installed commits inspected:
  Foundry `c97ad246877e2c7cb3b50e60aab7d76d02e8d393`, LigandMPNN
  `26ec57ac976ade5379920dbd43c7f97a91cf82de`.


## Recovered recent workflow after the user's source correction

The user identified `/data/jpc/enzyme-substrate-modeling/` as the more recent
source, superseding `/data/jpc/code/redesigner/` for this recovery. Found:

- `configs/nyl12_focused_redesign_pa6_20261016.json`.
- `scripts/lib/nyl12_focused_redesign.py`: `symmetry_tokens` and `mpnn_command`.
- `tests/test_nyl12_focused_redesign.py`: explicit separate alpha/beta tie tests.
- Current `output_nyl12_focused_redesign_pa6_20261016/protected_inputs/`, plus
  the dated workflow under `manuscript_nyl12_wt_vs_mutants/`.

The actual command uses explicit tied groups, weights0.25 for each of four
copies, `--homo_oligomer 0`, ligand-atom context and fixed-sidechain context.
The current preparation protects alpha113/186/222–232 and beta1–3/40/42.
These conservative maturation/anchor protections must be carried by mapped
identity, not treated as proof every protected residue is experimentally
essential. Original root and pilot snapshots are superseded and explicitly
marked DO NOT EXECUTE. The corrected snapshot prepared commands but states
no sequence-generation campaign was executed. Do not overstate that record.

Reuse the protected preparation and tested argument construction as the
starting pattern. Its mask is a candidate from an earlier noncovalent PA6
workflow, not an already-reviewed TI1/TS1 or loop-access mask. Its one-ligand
AE occupancy and equal symmetry weights are explicit assumptions; preserve
or compare them deliberately. The older redesigner used precursor-style
chains/automatic homooligomer mode, so copying it unchanged to eight mature
alpha/beta chains would be inappropriate.

The prior prose point-group declaration has not yet been located in the
bounded searches. The independent scaffold measurement supports D2-like
geometry; crystallographic space group P 2 21 21 in an old template is not
itself a biological point-group assignment. Keep these evidence types distinct.
Exact recovered source snapshots/hashes are under
`provenance/redesign_symmetry_and_gaff2_20261022/`.

## Parallel Amber preparation outcome

Standard antechamber GAFF2 typing and parmchk2 inventory were run on the
accepted S/-60 capped reference only. Charges were explicitly deleted with
`-c dc`; no RESP fit, native topology or protein refinement was performed.
Atom/bond order and the intended graph were preserved. The first audit
rejected serialization differences; its fresh successor checks the exact
three-decimal AC intermediate round trip. No chemical acceptance gate changed.

No ATTN placeholders were produced, but the inventory includes an empirical
angle and analog-derived torsions. The generic c3-os equilibrium distance is
1.427 A, versus approximately1.509 A for the QM reactive linkage. That generic
type pair also applies to Thr CB-OG1, so a global override could damage an
ordinary side-chain bond. These are targeted-review flags, not proof a fully
optimized MM model fails. Native RESP, complete MM validation and any justified
reactive-specific type/term fitting remain pending. Zero-charge typing files
are not eligible for energy calculations or MD.
