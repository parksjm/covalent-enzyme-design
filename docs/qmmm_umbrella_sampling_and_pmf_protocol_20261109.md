# Nyl12 QM/MM umbrella-sampling and PMF protocol

## Status and interpretation of the initial campaigns

The 20-window direct-proton-transfer and 25-window asynchronous E--S-to-TI1
campaigns use 1 ps heating, 1 ps biased equilibration and 2 ps biased
production per window. These short campaigns are engine, stability and overlap
pilots. They are not production PMFs and must not be used to report a barrier.

## Estimator decision

Use MBAR as the primary estimator and binned WHAM as an independent
cross-check. MBAR uses all cross-window reduced bias energies without choosing
histogram bins and provides an overlap matrix, effective sample information and
uncertainty estimates. WHAM remains useful because it is conventional for
enzyme QM/MM umbrella sampling and offers an implementation-independent check,
but its result must be tested against bin width and convergence tolerance.

Do not interpret agreement between estimators as proof of convergence. Sampling
duration is adaptive: 30 ps/window is the present literature-aligned allocation,
not the definition of convergence. Analyze cumulative retained blocks after
automated equilibration detection and require the following preregistered project
gates in *two successive* 5 ps additions before PMF promotion:

1. continuous neighboring-window overlap with no unsampled gap and adjacent
   histogram-intersection overlap of at least 0.20;
2. integrated-autocorrelation/statistical-inefficiency analysis for every
   window and decorrelation or an equivalent correlated-data treatment, with
   at least 100 effective uncorrelated samples per window for the biased
   coordinate after discarded equilibration;
3. a connected MBAR overlap matrix with no isolated state or bottleneck that
   disappears under modest resampling;
4. changes of no more than 0.5 kcal/mol in the barrier and basin free-energy
   differences and no more than 0.10 A in their coordinate locations between
   successive cumulative analyses;
5. a contiguous-block-bootstrap 95% confidence-interval half-width of no more
   than 1.0 kcal/mol for the barrier and basin free-energy differences;
6. MBAR and bin-width-sensitive WHAM profiles that agree within 1.0 kcal/mol
   over the thermodynamically relevant coordinate range;
7. independent-start or forward/reverse estimates whose barrier and basin
   differences agree within 1.0 kcal/mol or their joint 95% uncertainty; and
8. stationary component-distance and important orthogonal-coordinate
   distributions without unresolved drift, hysteresis or hidden multimodality.

These numerical values are declared Nyl12 project acceptance thresholds, not
universal constants. Failure of an overlap gate calls for new/intermediate
windows; failure of effective-sample-size or uncertainty/stability gates calls
for longer sampling; failure caused by hidden multimodality calls for a better
reaction coordinate rather than merely longer one-dimensional trajectories.

The asynchronous scalar coordinate

`xi = d(OG1-HG1) - d(N-HG1) - d(OG1-Ccarbonyl)`

can map distinct combinations of proton transfer and nucleophilic attack onto
the same value. It is also blind to proton identity: the 2026-11-09 extension
satisfied the HG1-based coordinate by transferring a different ammonium proton
(Thr929 H1) back to OG1. If conditional component-distance distributions are
multimodal, show hysteresis or permit proton exchange, do not solve that problem
merely by extending the 1D run. Use independent proton-transfer and attack
coordinates with proton-ownership diagnostics, followed by a targeted 2D
umbrella patch or a path/string coordinate in the affected region.

## Sampling plan

Retain the current pilots. After their overlap and structural diagnostics pass:

1. continue each accepted window from its restart for 5 ps of additional
   biased equilibration;
2. accumulate production in separately retained 5 ps blocks;
3. analyze every cumulative 5 ps increment rather than only selected nominal
   durations;
4. treat 30 ps/window as the initial publication-quality allocation, excluding
   heating and discarded equilibration, while applying the gates above; and
5. extend selected or all windows beyond 30 ps when overlap, effective sample
   size, uncertainty, PMF stability or independent-start tests fail.

With one campaign using violet, the initial measured rate was approximately
64 minutes per simulated picosecond per window. When all 45 direct and
asynchronous windows shared violet, the observed rate slowed to approximately
1.4--1.9 wall-hours per simulated picosecond. Consequently, one 5 ps block is
expected to require approximately 5.3 hours on a dedicated host or 7--10 hours
under the present shared load. Production remains segmented into retained
blocks for restart safety and convergence analysis, but sampling continues to
at least 30 ps/window.

The direct-PT extension is assigned to violet and the asynchronous extension
to azalea. Each host uses its native AmberTools 26 installation; architecture-
dependent executable hashes are recorded rather than treated as an error. The
topology and scientific input hashes must match their declared sources.

## Literature basis

- Pan et al., *J. Chem. Theory Comput.* **2023**, 19, 8234--8246,
  DOI: [10.1021/acs.jctc.3c00973](https://doi.org/10.1021/acs.jctc.3c00973).
  Their QM/MM protocol used 80 windows, Hamiltonian replica exchange, 50 ps per
  window and the final 30 ps for MBAR/TI analysis.
- Paz et al., *J. Comput. Chem.* **2025**, DOI:
  [10.1002/jcc.70109](https://doi.org/10.1002/jcc.70109). Their DFTB3/MM tests
  used replicated profiles and 10--40 ps production per window depending on
  model size, analyzed by WHAM.
- Maschio et al., *Chem. Sci.* **2024**, DOI:
  [10.1039/D4SC02908A](https://doi.org/10.1039/D4SC02908A). Their SCC-DFTB
  enzyme protocol used 2 ps per window only with ten independent repeats, or
  20 ps aggregate sampling per window; it does not support a single 2 ps
  production trajectory as a converged PMF.
- Zheng et al., *Commun. Chem.* **2026**, DOI:
  [10.1038/s42004-026-01888-w](https://doi.org/10.1038/s42004-026-01888-w).
  This polymer-hydrolase QM/MM study used 30 ps per window, analyzed the final
  15 ps with WHAM, bootstrapped the PMF and performed convergence and committor
  checks.
- Andjelkovic et al., *J. Chem. Inf. Model.* **2023**, 63, 5676--5688,
  DOI: [10.1021/acs.jcim.3c00900](https://doi.org/10.1021/acs.jcim.3c00900).
  Their adaptive-string protocol used neighbor Hamiltonian exchange and then
  path-coordinate umbrella sampling with WHAM. Their reported 2 ps criterion
  concerns string convergence, not evidence that 2 ps/window is sufficient
  production sampling for this Nyl12 PMF.
- Chodera, *J. Chem. Theory Comput.* **2016**, 12, 1799--1805,
  DOI: [10.1021/acs.jctc.5b00784](https://doi.org/10.1021/acs.jctc.5b00784).
  The automated equilibration criterion selects the retained region that
  maximizes the number of effectively uncorrelated samples.
- Feddersen and Biggin, *J. Chem. Theory Comput.* **2026**, 22, 5299--5310,
  DOI: [10.1021/acs.jctc.6c00246](https://doi.org/10.1021/acs.jctc.6c00246).
  Their workflow benchmarks emphasize direct convergence/error analysis and
  warn that apparently smooth PMFs can remain undersampled.

## Analysis and visualization

Use `pymbar==4.0.3`, pinned through `uv`, for MBAR analysis rather than a
project-local reimplementation. Treat PMF stability, effective sample sizes,
autocorrelation and adjacent-window overlap as scientific acceptance criteria;
the presence of a reconstructed PMF alone is insufficient.

`scripts/analyze_umbrella_overlap.py` accepts either one production-stage
directory or a stage glob such as `*production*`. It concatenates only stages
whose `process.json` records a successful exit. This permits cumulative overlap
analysis across retained 5 ps blocks without merging or modifying source
trajectories.

Use `scripts/export_umbrella_vmd.py` to create an immutable, stripped VMD
bundle. It samples frames evenly from each completed window in campaign order,
then concatenates the windows in reaction-coordinate order. The output includes
source hashes, a frame manifest, stripped topology, DCD and VMD script. The VMD
display follows the coordinate-driving convention: the QM region is licorice,
only MM substrate carbons are orange, reactive bonds are dynamic, and the
oxyanion-hole distances and Burgi--Dunitz angle are displayed. A sequential
umbrella trajectory illustrates each window's sampled ensemble; it is not a
continuous dynamical reaction movie.
