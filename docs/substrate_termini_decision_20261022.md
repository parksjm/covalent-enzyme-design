# Substrate termini and TI1 reference scope — 2026-10-22

## User clarification

The PA6 design objective is a specified oligomer product, particularly L2.
Retain a real N- or C-terminal group so designs can recognize its charge and
its distance from the scissile bond. The user is concerned that capping both
substrate ends would instead favor internal cleavage. This is a design
hypothesis to test, not proof that a capped model necessarily gives endo activity.

For the copolyamide, the user accepts considering a smaller capped model,
written `Ncap-XA/MA-Ccap`, with `/` indicating cleavage, and recalls using this
model in the Dishner JACS study. This notation is not yet an atom-level input:
record exact cap structures, protonation, monad direction and attacked C/O/N
before generating a graph. Do not infer cap identities from their labels.

## Separate two kinds of cap

The earlier proposed **Thr-NME cap** terminates the protein peptide fragment
after Thr. It does not cap either substrate end. The coupled TI1 reference
must retain Thr NH3+, unprotonated covalent OG1, intact scissile C–N, oxyanion
O− and neutral relay water. A larger protein fragment remains an alternative;
the user's substrate clarification does not resolve that fragment-size choice.

The existing PA6 loader selects C54/O49/N59 from the original full L4 MOL2.
Source inspection confirms NH3+ at N56 and COO− at C53/O48/O51, net formal
charge zero. Its strict hydrolysis check requires two L2 fragments. Preserve
that source and terminal chemistry as the **soluble-oligomer reference**.

## Oligomer turnover versus polymer-end recognition

The user further clarified that both charged ends describe a small soluble
oligomer, whereas one exposed end of a long polymer presents only one nearby
terminal charge. The user has detected both PA6 L4 and L2 products by MS and
therefore hypothesizes that L4 can itself be hydrolyzed to two L2 molecules.
This is a reasonable substrate hypothesis, not an established precursor–product
relationship from product detection alone.

Keep L4 → two L2 as the first PA6 reference for parameter work and soluble
oligomer turnover. It cannot by itself establish polymer-end recognition.
The separate polymer-end comparison would retain either a real N-terminus or
a real C-terminus two units from the cut, continue beyond the cut, and use a
neutral cap at the distant computational boundary. Compare the two orientations
and extend the distal tail to check boundary effects; a neutral cap alone does
not reproduce a long polymer. These are proposed later comparisons, not new
graphs, parameter transfers or jobs. Do not replace the existing L4 reference
silently or assume that the first detected product fixes the cleavage mechanism.

The MS-based activity assay detects only certain soluble, low-molecular-weight
oligomers, according to the user. Report the measured product panel and its
quantitation limits; do not equate this readout with all bond-cleavage events
or total polymer conversion. End-directed cleavage, internal cleavage followed
by further hydrolysis, and oligomer turnover can contribute to the observed
products. Direct incubation of L4 with monitoring of L4 depletion and L2
formation, with an enzyme-free control, would test the L4 → two L2 hypothesis.
This assay is proposed; no new experimental result is claimed.

## Copolyamide model roles

The existing product-defined reference is `M-A-X-A|M-A-X-A`, with the central
A carbonyl and M nitrogen selected. It is not the proposed short capped graph.
Keep it as a provenance-preserved product/register reference. A short capped
XA/MA model may support local parameter or catalytic-geometry calculations,
but it cannot by itself demonstrate production or selective release of the
complete MAXA tetrad. Its cut must be mapped by chemical connectivity rather
than relying on the written left-to-right sequence.

Primary reference: Dishner et al., *Circularity in Sequence-Controlled
Copolyamides Enabled by Regioselective Enzymatic Hydrolysis*, JACS (2026),
[DOI 10.1021/jacs.6c03466](https://doi.org/10.1021/jacs.6c03466).
The primary citation was located again on 2026-10-22; full-text/SI retrieval
was blocked by publisher access and a browser challenge. The exact capped
input has therefore **not** been independently verified in this resume.
The original modeling input is preferable to reconstructing caps from memory.

The user subsequently supplied the original capped-substrate SMILES verbatim:

```text
CC(NCC(C=C1)=CC=C1CNC(CCCCC(NCCCCCCNC(CCCCC(NC)=O)=O)=O)=O)=O
```

The scripted RDKit audit completed in tmux under
`output_nyl12_parameter_validation_20261022/substrate_termini_01/` and retained
the input, exact script, report, command, hashes and process status. It confirms
neutral **acetyl-X-A|M-A-N-methylamide**, formula C29H47N5O5, 39 heavy atoms and
five amides. The user-supplied provenance is distinct from independent SI
verification. Atom indices below refer to the **original supplied SMILES order,
zero-based**; they are not PA6 names or canonical-SMILES indices.

| Amide | Carbonyl C | Carbonyl O | Amide N |
| --- | --- | --- | --- |
| N-acetyl cap | 1 | 38 | 2 |
| X/A | 12 | 37 | 11 |
| A/M, selected slash reference | 17 | 36 | 18 |
| M/A, separate candidate bond | 26 | 35 | 25 |
| Terminal N-methylamide cap | 31 | 34 | 32 |

The slash in `Ncap-XA/MA-Ccap` selects C17–N18. Keep the other M/A bond as an
explicit alternative; do not silently substitute it. The old long-oligomer
atom names must not be copied onto this graph without a new atom crosswalk.

Frozen graph inputs, PA6 source MOL2, exact audit script, job records and
source digests are in `provenance/xa_ma_capped_reference_20261022/`.
The audit script's actual-input execution covers 30 of 31 executable lines
(96.77%); this is run coverage, not a new unit-test suite. The curated laptop
copy is under `references/xa_ma_capped_reference/`. The MOL2 is original
substrate reference chemistry, not a validated TI1 parameter set or minimum.

## Execution boundary

No substrate graph, coordinates, charges or QM input were changed by this
clarification. Next prepare an explicit coupled-state reference specification,
including how the charged intermediate is stabilized, before generating fitting
data. Keep the existing force-field validation gate: validated TI1 classical
parameters precede AF3-model MM relaxation, then DFTB3/MM refinement. No new AF3
prediction will run while the user is independently repairing the restraint
fork. Scripts, source hashes, tmux job records and the append-only action
journal remain required for subsequent calculations.
