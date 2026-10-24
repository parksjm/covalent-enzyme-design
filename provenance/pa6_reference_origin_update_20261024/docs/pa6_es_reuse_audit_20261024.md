# WT PA6 ES recovery and first validation checks — 2026-10-24

## Completed work

The accepted first stage of the proposed mechanistic study recovers the existing
WT PA6 L4 noncovalent ensemble. Original simulations and parameters remain
unchanged. No new production MD, umbrella sampling, TI1 fit or TS calculation
was started by this recovery step.

The source system is
`/data/jpc/enzyme-substrate-modeling/output_nyl12_pa6_md_pilot_100ns_20260818/systems/pa6_l4z_l2_nyl12_wt_model_41/`.
All three production logs reached 50,000,000 steps at 2 fs with normal completion
markers (100 ns each). Inputs declare distinct thermostat seeds 820117, 820134
and 820151. Different seeds do not establish statistically independent ensembles.
Each trajectory has 2,000 stored frames at 50 ps spacing.

Original LEaP input explicitly sources **ff19SB, GAFF and TIP3P**. GAFF here is
the first version; do not relabel these simulations GAFF2 or RESP. Original
parameter-cache library/frcmod hashes match the MD copies, and retained
AM1BCC intermediates support the reported AM1-BCC workflow. The exact historical
invocation/version was not recovered; current scripts alone cannot prove it.

Topology checks establish an intact noncovalent ligand and the declared
C54/O49/N59 scissile atom connections. Catalytic Thr929 has one OG1 hydrogen
and two N-terminal hydrogens: this dataset models **Thr–OH/NH2 ES**, not
Thr–NH3+/TI1. Phe113 backbone N and Asn186 ND2 retain their explicit donor
hydrogens. Manual N-terminal neutralization left a small Thr residue charge
residual of approximately −0.001 e; ligand charge is approximately −0.001001 e.
These values are preserved and recorded, not silently corrected or fitted.

## Native mapping and local environment

Direct parsing of the user-supplied `cif/n12c.cif` confirms these native author-A
identities. The declared full-topology residues have the same identities.

| Native author numbering | Full eight-chain Amber residue | Identity |
| --- | --- | --- |
| A113 | 113 | Phe |
| A186 | 186 | Asn |
| A233 | 929, mature beta E1 | Thr |
| A272 | 968, mature beta E40 | Asp |
| A274 | 970, mature beta E42 | Asp |

In the supplied reference, Asp272 OD1 is 2.83048 Å from Thr233 N and 4.92035 Å from OG1.
This supports examining the carboxylate/amine environment rather than treating
Asp272 solely as a steric obstacle. It does not establish protonation or a
catalytic mechanism. Correspondence to the Asp306/Asp308 network proposed for
another NylC remains a homology hypothesis pending sequence/structure mapping;
a constant residue-number offset is insufficient evidence. The user clarified
that this model comes from a coauthored manuscript currently in revision;
the formal citation will be supplied after publication. The CIF header says
`NYL12_CRYO-READY-FOR-MD` and lacks experiment/citation tags. It is retained as
the supplied structural reference with citation pending publication.

## Seed extraction and provenance

The versioned manifest is `configs/pa6_wt_es_reuse_20261024.json`. The reproducible
entry point is `scripts/audit_nyl12_pa6_md_reuse.py`, launched with
`scripts/run_nyl12_job.py` in tmux. Source hashes cover the original trajectories,
topologies, production inputs/logs, parameter files, reference CIF, CSV and code.
Original large trajectories remain in place; they are not laptop-export payloads.

Nine conditional ES seeds were selected: one per replica in each of the
predeclared 20–45, 45–70 and 70–100 ns windows. Candidates must satisfy historical
Tier B geometry; the eligible frame nearest the window midpoint is selected.
Geometry does not rank eligible frames. These seeds are starting conformations,
not a population estimate, an unbiased ensemble or optimized TI1 structures.

The legacy geometry CSV's time estimate starts at zero, while saved full frame 1
is at 0.05 ns relative to the production start. Ordinal CSV frame `k` maps to
full frame `10*(k−1)+1`. Both time conventions and their 0.05 ns offset are
retained. NetCDF absolute time may also include preceding equilibration.

Results:
`output_nyl12_pa6_es_audit_20261024/reuse_audit_attempt1/reuse_audit.json`.
Each seed retains all 103,251 atoms, 27,690 explicit waters, hydrogens and the
periodic box. CPPTRAJ images the complete assembly/solvent around the substrate.
The matching solvated topology remains the original unchanged source.

## Engine checks and next scientific gate

The installed Amber DFTB3 code was checked separately using its shipped water
minimization and 10-step periodic N-methylacetamide/water regression examples.
Both engine runs exited successfully; the compared printed energy series match
the shipped references exactly at their printed precision. Inputs, original
reference outputs, CHNO Slater–Koster hashes, executable hash and run logs are
retained under `output_nyl12_dftb_engine_checks_20261024/`.
These checks validate those installation regressions, **not nylon chemistry**.
The supplied NMA test's solvent is its own regression model, not a substitution
for the TIP3P model in the Nyl12 ES dataset.

Next inspect the observed water/carboxylate networks in the nine seeds, verify
their catalytic metrics independently, and define explicit QM fragments, boundary
cuts and formal charge budgets. Then run bounded DFTB3/MM PES/force checks and
selected hybrid-DFT PES comparisons. Native TI1 classical relaxation remains
blocked until its force field is validated. A successful ES recovery or engine
regression does not waive that requirement, authorize production free energies,
or demonstrate a first-order saddle point.
