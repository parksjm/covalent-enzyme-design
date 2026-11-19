# Parallel follow-up during TI1 DFT/MM relaxation

The first three DFT/MM evaluations took approximately 7.5 minutes with 28
Gaussian threads. Budget approximately 2--3 minutes per evaluation, or about
4 hours for 100 evaluations and 18--25 hours for 500; line searches and SCF
changes can increase this. Convergence timing is not yet predictable.

The scripted second-production-block audit and cumulative 12-ps/window MBAR
analysis completed in separate tmux sessions. Results are under
`07_direct_pt_umbrella_30ps_violet/analysis/mbar_12ps_complete/result/` and
`analysis/proton_state_block02_5ps/result/`. The cumulative analysis has 48,000
raw samples and 28,268 decorrelated samples. Minimum adjacent overlap is 0.183;
global overlap scalar is 0.0128. Compared with 7 ps/window, the common-domain
PMF changes by 0.414 kcal/mol RMS and 0.562 kcal/mol maximum. This is an interim
sampling-stability assessment, not proof of convergence or a barrier estimate.

Topology inspection identifies the scissile carbon as atom 19916, its oxygen
as 19917, and its bonded leaving-group nitrogen as 19929 (residue 1315).
The next acylation-coordinate set must therefore monitor OG1(13916)--C(19916),
C(19916)--N(19929), C(19916)--O(19917), all terminal N(13905)--H distances,
and proton delivery to N(19929). The proton donor/identity must be selected
from inspected geometry; do not assume that HG1 necessarily supplies the
leaving-group proton. No new cleavage-driving calculation has been launched.

Gradient audit: Gaussian benchmark logs contain analytic forces, and xTB
retained `gfn2.gradient`, `gfn2.pcgrad` and `gfn2.qm_link.engrad`. A matched
projection is still pending. It must preserve atom order, distinguish forces
from gradients and Hartree/bohr from kcal/mol/angstrom, and treat link atoms
consistently. Do not compare total hybrid forces directly against isolated
QM/link forces. DFTB3 reactive-gradient extraction must be established before
reporting a three-method direction comparison.
