# TI1 parameterization status — 24 October 2026

Native TI1 parameterization is paused at scientific validation, not accepted
for classical refinement or MD. The enzyme-specific reactive QM/MM pilot is
the next discriminating calculation, rather than another isolated-cluster fit.

Completed work retained in the repository and original output directories:

- Four constrained, 37-atom capped TI1 optimizations, exact Hessians and
  conditional-curvature assessments passed their individual state checks.
- Separate capped-surrogate HF/MK ESP calculations passed their ESP checks.
  These are not native RESP charges or a completed force-field fit.
- The native attack configuration is R. Its cap-sensitivity test failed:
  OG1–C* changed 0.05327 Å and C*–O− changed 0.03166 Å, exceeding the
  predeclared 0.03 Å threshold. Selecting the passing S configuration would
  not repair this result.
- Subsequent native/local cluster pilots either were interrupted or lost the
  intended TI1 state. Adding donors/water did not establish an accepted native
  reference. The rejected water placement and interrupted SMD experiment do
  not establish whether continuum stabilization would work.

Still required before native TI1 classical refinement:

1. Establish an enzyme-compatible, explicitly protonated intermediate and
   reference geometry, with selected hybrid-DFT checks of the reactive model.
2. Obtain appropriate native/capped parameter references, assess boundary/cap
   sensitivity and calculate the derivatives needed for bonded fitting.
3. Fit and validate native ESP/RESP charges separately from geometry and bonded
   parameters; test transfer, held-out geometries and important torsions.
4. Validate the assembled ff19SB/nonstandard-residue force field before using
   it to relax TI1 protein models or run classical TI1 MD.

The reactive DFTB3/MM pilot starts from the audited **noncovalent ES** topology.
The reacting atoms are QM, so this calculation does not require a fitted
classical TI1 bond. It retains the original ff19SB/GAFF1/TIP3P MM environment,
the full PA6 L4 substrate and its charged ends. This is distinct from the
prospective ff19SB/GAFF2/RESP TI1 parameterization workflow; no charge model
or parameter family is silently upgraded.

Two independently audited, charge −2 singlet QM regions contain 98 or 104
original atoms plus eight link hydrogens (106/112 QM nuclei). They include
neutral terminal Thr NH2/OH, the local reacting amide, the Phe backbone donor,
Asn side-chain donor and two deprotonated Asp side chains. The larger variant
moves the PA6 cuts outward. This starting protonation is an explicit assumption,
not an experimentally established catalytic microstate. No forced water is added.

The first gate is reproducible finite energies and finite analytic forces,
followed by a separately declared bounded ES relaxation and selected chemical/PES
tests. The optional finite-difference check has a recorded, unresolved small
step-dependent discrepancy; after the user's questions about analytic gradients,
root stopped further numerical-force work and removed it as a prerequisite
for this bounded exploratory relaxation. Its original failed record is retained.
Passing implementation checks alone does not validate DFTB3 for nylon
hydrolysis or establish a TI1 minimum, TS1 or reaction free-energy barrier.
If the relevant enzyme mechanism has no metastable TI1 basin, do not force a
classical minimum merely to complete the parameterization checklist. Revisit
which state-specific classical models are scientifically useful in that case.

See `docs/ti1_qm_results_20261023.md`,
`docs/pa6_dftb_validation_pilot_20261024.md` and the current handoff for actual
job status. The separate local-agent idea is deferred for later consideration;
Chinese-origin LLMs and their derivatives are excluded by the user's constraint.
