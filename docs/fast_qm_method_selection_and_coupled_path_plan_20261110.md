# Fast-QM selection and coupled Nyl12 reaction-path plan

## Decision prompted by the enzyme-path benchmark

The small gas-phase proton-transfer models favored GFN2-xTB over DFTB3, but
the five-point embedded enzyme test did not. For the exact 98-atom QM region,
eight link atoms and matched per-geometry electrostatic embedding, DFTB3 had a
1.49 kcal/mol double-difference RMSE against omegaB97X-D/def2-SVP, whereas
GFN2-xTB had a 3.91 kcal/mol RMSE and a +7.93 kcal/mol endpoint error. Method
quality is therefore reaction- and environment-dependent. Neither result is a
sufficient production-method decision because the five geometries contain
proton transfer but not the carbonyl polarization, nucleophilic attack and
oxyanion formation required for TI1.

This conclusion is consistent with recent biochemical proton-transfer
benchmarks: average semiempirical performance conceals substantial
functional-group dependence, and N--H/charge-separated reactions remain
challenging. Coupled enzyme coordinates can also give qualitatively similar
profiles but quantitatively different intermediate and barrier energies among
DFTB3, GFN2-xTB and hybrid DFT. The production method will consequently be
selected on the actual coupled Nyl12 chemistry, not on gas-phase proton
affinities alone.

## Immediate corrected attack branch

The chemically valid `asynchronous_24` endpoint is the immutable seed. It has
HG1, H1 and H2 bound to the terminal N and `d(OG1,C) = 2.104 A`. The invalid
continuation is not reused: it allowed H1 to transfer back to OG1 because the
old scalar coordinate tracked only HG1.

The corrected branch drives only `d(OG1,C)` from 2.104 to 1.50 A in increments
no larger than 0.075 A. Three flat-bottom N--H gates retain HG1/H1/H2 on N
without fixing their bond lengths: the unpenalized range is 0.90--1.20 A and
the outer walls are 0.80 and 1.35 A. The accepted 12 A DL-FIND movable shell
remains active; there are no global positional restraints. Each window must:

- converge natively;
- retain every N--H distance at or below 1.25 A and every corresponding
  OG1--H distance at or above 1.30 A; and
- preserve finite coordinates and the existing topology/cell.

Attack-target lag is recorded but is not a hard acceptance limit. Attempt 2
showed why: `attack_04` converged, retained every intended proton identity and
advanced OG1--C to 1.950 A, but lagged its nominal target by 0.115 A. The
predeclared 0.10 A cutoff stopped a chemically valid path for an arbitrary
numerical reason. Attempt 3 therefore records an immutable acceptance revision
for that window and resumes from its converged endpoint in the same calculation
directory. Convergence and chemical identity remain hard gates.

These biased endpoints are method-benchmark seeds, not stationary points or a
claimed TI1. Attempt 1 stopped before Sander because its minimal source bundle
omitted a transitive visualization module. Attempt 2 passed local and remote
source-hash/import checks and stopped only on the superseded tracking limit.
Attempt 3 is running on azalea in tmux session
`nyl12_corrected_attack_window24_attempt3_20261110`. The canonical directory is
`qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths/02_asynchronous_es_to_ti1/03_corrected_attack_from_window24/`.

## Reaction-specific Hamiltonian panel

Sparse geometries will be selected only after the corrected branch passes its
chemistry gates. The minimum set is E--S, proton-transferred, early attack,
late attack and TI1-like. Each exact geometry receives matched energy and
gradient evaluations with:

1. Amber DFTB3/3OB/MM (current baseline);
2. Amber PM6-DH+/MM, plus PM6/MM as a control;
3. standalone GFN2-xTB with the exact Amber QM/link geometry and adjusted
   point charges; and
4. Amber external omegaB97X-D/def2-SVP/MM through Gaussian 16 as the reference.

PM6 and PM6-DH+ are already supported by the installed AmberTools build, so
this panel introduces no production dependency. DFTB3-D3H5 is a later candidate
only if the installed DFTB3/PM6-DH+/xTB panel fails: it requires a DFTB+ build
and explicit user approval. The development `g-xTB` implementation is not a
production candidate because its public implementation does not support the
point-charge embedding required here.

The comparison uses method-minus-reference double differences, because the
finite external-charge membership can otherwise change the arbitrary
electrostatic zero for the net-charged QM region. Where possible, the same
charge indices and periodic images will be frozen across the sparse set. Raw
finite-cutoff profiles with changing charge membership are rejected.

## Promotion gates for a fast method

A method is promoted only if all of the following hold on the coupled set:

- the same proton identities and covalent connectivity are retained;
- E--S, transferred-proton, attack and TI1-like states have the same
  qualitative ordering as the DFT reference;
- relative-energy RMSE is no more than 3 kcal/mol and no endpoint error exceeds
  5 kcal/mol (reported as targets, not guaranteed outcomes);
- projected reactive-coordinate gradients have the same signs as DFT and a
  force-direction cosine of at least 0.8 on the declared reactive atoms; and
- short state-release tests do not create a spurious proton swap or destroy a
  DFT-retained basin.

No method is selected solely by gas-phase RMSE or the isolated direct-proton
coordinate. If no fast method passes, the fallback is DFTB3/MM sampling with
sparse DFT/MM correction/reweighting only after overlap is demonstrated, or a
GPU-capable DFT/MM implementation. It is not acceptable to force an incorrect
Hamiltonian through the desired mechanism.

## Reaction-path calculation after method selection

The productive path will be sampled in two coordinates:

```text
xi_PT     = d(Thr233 OG1, HG1) - d(terminal N, HG1)
xi_attack = d(Thr233 OG1, substrate carbonyl C)
```

A two-dimensional restrained grid or string seeded from the corrected path
allows proton transfer and attack to be stepwise, concerted or asynchronous
without collapsing both events into the chemically degenerate scalar used in
the failed continuation. H1 and H2 identity walls remain explicit safety gates;
HG1 is governed by `xi_PT`. Window overlap, proton identities, all forming and
breaking distances, carbonyl C--O elongation, oxyanion-hole contacts and the
Burgi--Dunitz angle are updated during sampling. A two-dimensional MBAR PMF is
accepted only after blockwise stability, adequate neighbor overlap and stable
minimum-free-energy-path location.

The E--S, any real intermediate, TI1 and PMF saddle seeds are then optimized
without reaction-coordinate restraints. Candidate transition states receive
DFT/MM DL-FIND/eigenvector-following refinement, a Hessian with exactly one
relevant imaginary mode, and forward/reverse connectivity validation. Those
validated stationary points—not biased coordinate-drive endpoints—become the
multi-state structures used for LigandMPNN redesign.

## References

- Tao, P.; Giese, T. J.; York, D. M.; Woodcock, H. L. *J. Chem. Phys.* **2024**,
  *160*, 224104. DOI: 10.1063/5.0211276.
- *Benchmark of Approximate Quantum Chemical and Machine Learning Potentials
  for Biochemical Proton Transfer Reactions*. *J. Chem. Theory Comput.*
  **2025**. DOI: 10.1021/acs.jctc.5c00690.
- Gaus, M.; Jin, H.; Demapan, D.; Christensen, A. S.; Goyal, P.; Elstner, M.;
  Cui, Q. *J. Chem. Theory Comput.* **2018**, *14*, 69--79. DOI:
  10.1021/acs.jctc.7b00600 (DFTB3-D3H5).
