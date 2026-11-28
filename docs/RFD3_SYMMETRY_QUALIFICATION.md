# RFD3 optional-arm qualification and retained symmetry audit

Updated 27 November 2026. RFD3 is outside the current critical path. Review the existing 12 original
MD/TS-backbone LigandMPNN designs from the two chemically coherent MD-derived
parents, retaining the 101-position protection mask and mandatory four-copy
sequence ties. No RFD3 GPU qualification or 12-design replacement is required.
Identical sequences need not have identical thermal conformations; do not force
exact D2 coordinates onto an MD snapshot. Multiple chain conformations alone
also do not establish TS compatibility: retain reaction mapping, proton inventory,
ligand graph and catalytic geometry for each parent.

Fixed-backbone LigandMPNN design followed by independent template-free,
paired-MSA Boltz prediction and matched Amber relaxation can assess
sequence-driven repacking and assembly; it does not promise improved activity.
A future optional RFD3 arm may evaluate coupled perturbations of A/E, B/F,
C/G and D/H, with independent physical-movement, clash, assembly, fixed-context
and chemical checks. It is not implemented or qualified for production, and
exact coordinate symmetry is not a universal acceptance requirement.

## Retained 27 November audit history

The [original audit](../design/nyl12_relay_ts_ensemble_20261126/backbone_qualification/symmetry_qualification_20261127/symmetry_contract_audit.json)
found default-sampler, chain-A-only local remodeling without symmetry inputs.
Earlier graph/coordinate receipts retain their stated scope; these remodels
remain diagnostics rather than current baseline parents. Approximate crystal
D2 NCS does not establish exact operators or TS compatibility.

The [exact-projection rejection](../design/nyl12_relay_ts_ensemble_20261126/backbone_qualification/symmetry_qualification_20261127/exact_projection_rejection.json)
records 10 inter-repeat contacts below 1 Å and a 17.5 Å displacement of alpha
residue 232 in the exact A/E-copy audit. Exact projection was rejected.
Two actual CPU preflights failed on an E233 canonical/local numbering mismatch
before feature generation. No GPU qualification launched. Prototypes are
parked and unqualified; preserve their evidence and all prior outputs.

## Scope of the optional exact-symmetry gate

`rfd3_symmetry_acceptance.py` validates a hash-bound acceptance receipt only
for a future explicitly symmetry-constrained experiment. The current
`relay_broad_runtime.verify_parent` does not invoke it. It is an optional
helper, not a universal MD-parent requirement. Its D2/native-sampler, equivalent-mask, operator-closure and
residual checks describe that specific experiment. They do not measure
coordinates or independently establish physical plausibility, fixed-TS
compatibility or chemistry. No production-qualified optional RFD3 arm is
established. Any later coupled-perturbation study needs a prespecified physical
hypothesis, motion masks, canonical maps, ligand/proton checks and independent
assembly/clash acceptance before execution; no new run is authorized here.

## Experimental testing remains a separate decision

QM/MM PMFs and verified PES transition states are not prerequisites for PI
experimental nomination. Sequence, maturation/assembly, structural and available
matched MD evidence can support a bounded WT/design experiment. Optional
matched QM/MM calculations may be run before testing if resources permit, or
afterward to explain an observed activity difference. Neither geometry nor a
sequence score establishes a lower activation free energy.
