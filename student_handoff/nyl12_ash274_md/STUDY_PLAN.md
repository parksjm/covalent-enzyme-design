# Study plan and measurable gates

## Hypothesis and scope

Test a mutation-by-protonation interaction in the **noncovalent** mature
Nyl12–PA6 L4 complex. Primary comparison: ASP272/ASP274 versus
ASP272/ASH274, with neutral catalytic N termini. This is a fixed-microstate MD
experiment, not a calculation of equilibrium protonation populations.
A favorable result supports pose persistence under that model; it does not
establish a faster chemical step or experimentally improved conversion.

The first substrate is the released zwitterionic PA6 L4, representing the
soluble oligomer hypothesis L4→two L2 products. Both charged substrate termini
are therefore intentional. Do not treat it as a model of every polymer end or
bulk nylon. PA66, arbitrary PA6 lengths, other terminal caps, temperature changes
and catalytic adducts are outside the initial comparison.

## Stage 0 — reproduce and audit

Reproduce the cached six WT counts and SI panels before running new MD. Check
all parameter/system SHA values, protein sequence, four processed alpha/beta
pairs, scissile bond and per-site residue/atom mapping. Record the PA6 release
identifier and manifest SHA. Verify substrate bonds and charge, neutral beta
N termini, all four Asp microstates and whole-system charge/counterions.

The transferred WT starting coordinates are MD-prepared MM structures, not
certified stationary points. Use their exact acceptance receipts. Hardware
changes need a short stability check, not bitwise trajectory reproduction.
Never infer acceptance merely from the existence of a restart or exit code 0.

## Stage 1 — matched WT protonation test

Use the same starting heavy-atom E–S pose/model, protein/substrate force fields,
water model, box/salt recipe, equilibration schedule, production duration,
analysis cadence and scoring definitions in both microstates. Rebuild topology
and counterions when charge changes; do not rename a residue in a PDB and reuse
its old charged topology. Changes affect all four symmetry copies.

Start with at least **three independent 100 ns replicas per condition**;
prefer four given the observed 0-versus-16.77% disparity. Existing replicas can
count only if provenance confirms the same protocol. Different velocities/seeds
are required; continuing one trajectory does not create an independent replica.
Report all replicas, including dissociation and poor geometries. Do not select
only the successful seed. ASH272/ASP274 is an optional explicitly labeled
sensitivity control; it is already available for historical analysis.

The supplied primary/baseline starter configs contain replicas 01/02. Their
`tools/config_wt_asp272_ash274_wave2.json` and
`tools/config_wt_asp272_asp274_wave2.json` add replicas 03/04 with distinct seeds,
outputs and session names. Launch a second wave only after all jobs using its
GPUs 0/1 have finished; see `tools/README.md`. No wave is queued automatically.
Only WT systems are delivered ready to run; mutant rebuilding remains Stage 2.

The existing validated schedule is: accepted unrestrained MM preparation;
100 ps heating with 5 kcal/mol/Å² solute positional restraint; 100 ps density
relaxation at 2; 100 ps weak equilibration at 0.5; 200 ps unrestrained equilibration;
then unrestrained 10 ns+90 ns production. Preserve the actual generated mdins,
thermostat/barostat and integration settings rather than reconstructing them
from this prose. Production has no catalytic-contact or positional restraints.
Use 300 K initially. The assay used 75 °C, 200 mM Tris, pH 8,6 h; temperature of pH
adjustment is unknown. Assay-temperature work at 348.15 K is a separate arm.

## Stage 2 — prior single mutants

The verified historical panel is WT plus T53Q, F55N, K59D, K59E, F101E, F101H,
P105K, V111C, each with three completed legacy trajectories. Their PA6 treatment
was AM1-BCC/GAFF. Rebuild each chosen genotype with the **same current** modular
RESP/GAFF2 release and compare both Asp274 microstates where feasible.
Do not use old-FF-versus-new-FF differences as evidence for an ASH274 effect.
No completed new-modular-FF mutant panel was identified in this inventory;
add verified manifests if another project supplies it.

Within each genotype, use the identical starting structural model for the
paired protonation comparison. Across genotypes, record model confidence,
initial geometry and strain: previous single-mutant models were often poorer
than WT. A subsequent matched-WT-backbone mutation arm can separate sequence
effects from prediction/model-quality effects. Do not silently merge the arms.

Prepare new systems from genotype-specific coordinates and the released
parameter/topology-building recipe. The portable runner validates prebuilt
systems; it does not perform mutation modeling. Preserve all four copies,
processed termini, chain order and PA6 graph. Require explicit preparation
acceptance and matching restart/topology hashes before launch.

## Analysis and decision criteria

Historical Tier A: d1≤3.2 Å; d2,d3≤3.8 Å; OG1–C–O Bürgi–Dunitz 90–120°;
70≤|projected phi|≤110°. Tier B and the precise signed-phi implementation are
frozen in the supplied repository source. d1=ThrOG1–scissileC,
d2=Phe113N–carbonylO,d3=Asn186ND2–carbonylO. Report the signed value, not just |phi|.

Report per-replica full 100 ns and final 50 ns occupancies, 10 ns blocks,
geometry distributions/histograms and time series, binding/contact persistence,
Asp-pair minimum over all four OD–OD combinations, ASH hydrogen bonds,
protein/substrate RMSD, and Thr233 orientation/native displacement. The stricter
reference gates are native-aligned OG1≤0.75 Å, Nalpha≤1.25 Å, |Δpsi|≤35° and
|Δchi1|≤20°, plus intact substrate. Keep Tier A-only and Tier A/B joint counts
separate. Preserve failed/missing frames in explicit denominators.

Use independent replicas as the primary statistical unit. Frames are correlated;
four symmetry sites within one tetramer are also not independent replicas.
Report effect sizes and replica/block uncertainty, and assess early/late
agreement before extending sampling. A single percentage or a frame-wise
binomial significance test is insufficient. Declare a practical effect-size
threshold before ranking; calibrate it against WT replica variability.

The existing matched comparison script reproduces **site AE** counts. Before
using MD to rank designs, extend/report four-site-aware binding and catalytic
occupancy: any equivalent active site may be occupied. Preserve the historical
AE metric alongside the new any-site metric, record site identity/hops, and
count each frame once for any-site occupancy. Imaging and atom mapping must be
validated; apparent unbinding from AE alone is not proof of dissociation.

Proceed to broad mutation/design screening only if the microstate advantage
persists across independent replicas and late-time blocks without worsening
native-Thr geometry, substrate integrity or binding/access. Otherwise revise
the hypothesis or sampling; do not optimize sequences to a misleading proxy.

## Stage 3 — conditional LigandMPNN exploration

Start with fixed-length, focused sequence design on accepted E–S ensembles.
This is exploratory binding-pose design, not validated TS-guided catalysis.
Protect F113, K156, N186, T233, D272, D274 in all four copies; initially also freeze
cleavage-junction 232–233, Thr234, Gly273 and critical core/interface contacts.
A prior mutation at a protected position would be a labeled mechanistic
control, not an automatically eligible candidate.

Use explicit equal-weight ties among alpha A/B/C/D and separately beta E/F/G/H.
Do not blindly treat eight unequal alpha/beta chains as one homo-oligomer.
Validate every generated sequence against the fixed-residue and symmetry maps.
No indels or backbone redesign in this first arm. Use several accepted starting
poses rather than one unusually favorable snapshot.

Generate individual and coordinated substitutions, retain WT controls, model
all candidates with one declared Boltz/AF3 protocol, relax with the same released
FF and apply geometric/preparation gates before MD. Compare new designs with
matched WT at the same microstate and sampling budget. Maintain sequence
and structural diversity; do not select only increasingly closed or rigid
pockets. Assess substrate access/product release and, later, compatibility with
refined catalytic states. Final activity claims require experiments and/or
credible reaction-state energetics; Tier A occupancy alone is not sufficient.
