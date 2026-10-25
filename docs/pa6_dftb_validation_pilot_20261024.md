# First enzyme-specific DFTB validation pilot — prospective, 2026-10-24

The WT ES reuse audit and installed-engine regressions have passed. This note
specifies the next scientific gate; it is not a completed QM/MM calculation.

Update after user review: Sander supplies analytic QM/MM gradients. The extra
finite-difference sanity check found a small, step-dependent discrepancy and
is retained as failed evidence; its cause has not been established. The user
questioned making numerical derivatives a prerequisite. Root therefore permits
only a separately declared, bounded exploratory ES relaxation using analytic
gradients after source, inventory and realized-physics checks. This does not
certify the DFTB reaction surface or waive subsequent chemical validation.
No further numerical-force calculations are planned for this pilot.

1. Use the recovered single-ligand PA6 L4/tetramer seeds, with their original
   full MM environment. Preserve full substrate graph, terminal charges and
   the central L2-producing register. Do not transplant AF3 atom names.
2. Propose a QM core containing the catalytic N-terminal Thr, local reacting amide,
   backbone oxyanion donor, Asn side-chain donor, and the two nearby Asp groups.
   Resolve literature correspondence separately from native residue identities.
   Preserve donor amide chemistry: a link cut directly through the donor's
   C(O)–N bond would change the model being tested. Keep the complete PA6 L4
   and its charged termini in the enzyme/MM model, while truncating the QM
   substrate at declared saturated C–C bonds away from the reacting amide.
   The exact capped atom map must be verified before preparing QM inputs.
3. A useful *capped reference accounting model* is Thr–NME + a neutral capped
   local amide fragment +
   N-methylacetamide + acetamide + two acetate groups. With neutral Thr–OH/NH2,
   neutral local substrate and two deprotonated carboxylates, this is formally −2
   before adding neutral water. This is a candidate reference, not an accepted
   QM mask. Explicit caps, hydrogen placement, atom identities and MM boundary
   treatment must be audited before execution.
4. Include water only from declared observed configurations, or label any newly
   generated configuration as a separate preparation experiment. Retain explicit
   hydrogen identities and TIP3P MM Lennard-Jones bookkeeping. Select the water
   set once for a path; do not add/remove QM atoms dynamically by distance.
   A water contacting both N and OG1 does not by itself establish proton relay.
5. Check finite energies/SCF, actual link atoms and charge redistribution, and
   finite analytic forces before the bounded exploratory relaxation. The
   finite-difference diagnostic and its limitation remain separately recorded.
   Test at least one boundary/embedding variant before accepting the reactive
   Hamiltonian. The engine's
   passed water/NMA examples do not replace these enzyme-specific checks.
6. Start with bounded ES PES relaxations and proton-conserving local reaction
   models. Keep the total atom/proton inventory fixed when moving a proton
   between Thr, water and a carboxylate. Different protonation states require
   separately declared charge/atom inventories; changing only qmcharge is not
   a valid way to generate them.
7. Compare selected local states and candidate paths with the agreed
   omegaB97X-D/def2-SVP PES protocol, followed by larger/diffuse-basis checks
   where anion energetics require them. No additional dispersion correction
   is added to omegaB97X-D. Any continuum calculation is a labelled diagnostic
   with a different environment, not a substitute for embedded QM/MM or for
   the declared HF/6-31G* ESP/RESP charge protocol.

Only after this pilot passes should a short sampling/CV test be expanded into
production umbrella sampling. Both acylation and deacylation remain in scope.
No first-order TS, converged PMF, native TI1 classical force field or improved
enzyme has been established by the ES recovery work.

The prior MD's neutral terminal amine is a declared starting-state assumption,
not a measured dominant protonation state at assay pH. A conditional PMF must
state its protonation ensemble; accessing a reactive microstate may carry a
population cost in addition to its conditional chemical barrier.
