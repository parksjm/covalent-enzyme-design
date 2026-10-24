"""Polish and retain the proposed mechanistic-study scope without simulation."""
from pathlib import Path

ROOT = Path('/data/jpc/covalent-enzyme-design')
path = ROOT / 'docs/nyl12_mechanistic_project_scope_20261024.md'
text = path.read_text()
replacements = {
    'more than20years': 'more than 20 years',
    'independent100ns WT PA6L4': 'independent 100 ns WT PA6 L4',
    'October9': 'October 9', 'September28': 'September 28',
    'in3/3WT': 'in 3/3 WT', 'ExistingMD': 'Existing MD',
    'repeating3x100ns': 'repeating 3 × 100 ns',
    'hybridDFT': 'hybrid DFT', 'etal(2014)': 'et al. (2014)',
    'JPCL5,1210': 'JPCL 5, 1210', 'DOI10.': 'DOI 10.',
    'near21kcal/mol': 'near 21 kcal/mol', 'NtnThr': 'Ntn Thr',
    'JMB2007': 'JMB 2007',
}
for old, new in replacements.items():
    text = text.replace(old, new)
path.write_text(text)
path = ROOT / 'docs/standalone_handoff_20261020.md'
text = path.read_text()
marker = '## Literature/manuscript checkpoint — 2026-10-24'
_, remainder = text.split(marker, 1)
path.write_text('''## Proposed separate WT PA6 mechanistic project — 2026-10-24

User supports considering a full chemical-cycle DFTB QM/MM free-energy study
as a standalone paper. Focus PA6, optionally poly(MAXA), without expanding PA66.
Use hybrid DFT only for selected PES optimizations/checks, not free-energy
sampling. QM/MM validation and convergence are substantive research tasks.
Prior antechamber/AM1-BCC WT noncovalent PA6 L4 MD (3 × 100 ns) already exists;
corrected records show PA66 mainly loses local geometry/register, without
necessarily dissociating. Pin valid topologies/trajectories before reuse;
do not repeat simulations or change charges without an explicit rationale.
NylB is a comparator, not an authorized pivot; its published 2014 reactive
QM/MM study is DOI 10.1021/jz500323y. Scope and validation gates, including the
distinction between a reactive QM core and state-specific MM parameters, are
in docs/nyl12_mechanistic_project_scope_20261024.md. No new simulation started.
A full-cycle mechanism would support a paper; it is not required for every
exploratory sequence. A PMF maximum is not automatically a PES saddle point.
Updated scope and references are included in the existing laptop pull package.

''' + marker + remainder)
path = ROOT / 'manuscript/SOURCE_LEDGER.md'
text = path.read_text()
marker = '## Proposed PA6 mechanistic project and prior MD — 2026-10-24'
prefix, _ = text.split(marker, 1)
path.write_text(prefix + marker + '''

The user's three 100 ns WT PA6 L4 persistence report and AM1-BCC provenance are
retained in ../docs/nyl12_mechanistic_project_scope_20261024.md with the
corroborating October 9 original-project synthesis. That source distinguishes
local PA66 geometry loss from global dissociation and preserves corrected versus
invalidated histories. Exact MD topology/protonation hashes await a reuse audit.
DFTB sampling with sparse hybrid DFT PES checks is a proposal only; no new MD
or QM/MM simulation was started. No NylB pivot is authorized. Primary references
R34/R35 document existing NylB QM/MM work and approximate-QM proton-transfer
benchmarking. Source access limits are recorded in REFERENCES.md.
''')
path = ROOT / 'manuscript/recent_design_and_relay_references.bib'
with path.open('a') as stream:
    stream.write('''
@article{R34,
  author = {Kamiya, Katsumasa and Baba, Takeshi and Boero, Mauro and Matsui, Toru and Negoro, Seiji and Shigeta, Yasuteru},
  title = {Nylon-Oligomer Hydrolase Promoting Cleavage Reactions in Unnatural Amide Compounds},
  journal = {The Journal of Physical Chemistry Letters},
  year = {2014},
  volume = {5},
  number = {7},
  pages = {1210--1216},
  doi = {10.1021/jz500323y},
  url = {https://pubmed.ncbi.nlm.nih.gov/26274473/},
  note = {Primary abstract verified; full computational methods not checked}
}

@article{R35,
  title = {Benchmark of Approximate Quantum Chemical and Machine Learning Potentials for Biochemical Proton Transfer Reactions},
  journal = {Journal of Chemical Theory and Computation},
  year = {2025},
  doi = {10.1021/acs.jctc.5c00690},
  url = {https://pubs.acs.org/doi/10.1021/acs.jctc.5c00690},
  note = {Working record; author list and full bibliographic metadata pending; primary indexed abstract verified}
}
''')
