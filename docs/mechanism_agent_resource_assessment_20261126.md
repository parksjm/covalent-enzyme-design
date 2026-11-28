# Mechanism-agent resources: source assessment, 26 November 2026

## Recommendation

Improve evidence selection and scientific verification before changing model
weights or orchestration frameworks. The completed project review exposed
source-version, benchmark-identity and numerical-interpretation errors even
when citation formatting passed. General chemistry benchmarks help diagnose
ability; they do not qualify enzyme-mechanism planning or scientific dispatch.
The current PI/GPT-led scientific workflow and November 30 priorities remain.

This assessment authorizes no installation, model startup, post-training,
framework migration or new scientific campaign. Existing QM/MM and umbrella
sampling jobs were not changed. Model servers remain off following the earlier
supervised test; this research did not restart or independently recheck them.

## LitmusBench v0.2

The [NVIDIA dataset card](https://huggingface.co/datasets/nvidia/Nemotron-RL-litmus-bench-v0.2/raw/main/README.md)
reports **47,118 training and 2,056 test questions: 49,174 total**, 64 MB,
126 properties and 12 assay endpoints. It emphasizes short-answer molecular
property/ADME tasks from ChEMBL, ExpansionRx and public Biogen data, including
RDKit descriptors. ChEMBL-derived questions use CC BY-SA 4.0; others use
OpenMDW 1.1. Preserve source-specific terms if reusing data.

Use a bounded, stratified diagnostic subset with a frozen question list,
separate tool/no-tool results and documented numeric tolerances. Neither its
test label nor public availability establishes independence from Nemotron's
training. Do not call it a clean held-out evaluation of our deployed checkpoint
without a contamination audit. Training and evaluation records must stay
separate. We have not evaluated or trained Nemotron on this dataset.

The card was read successfully through its raw README; direct dataset/API/viewer
access failed, so no individual records or full row schema were verified.
Do not infer reaction-mechanism or QM/MM coverage from the total question count.
[NeMo Gym's Litmus environment](https://github.com/NVIDIA-NeMo/Gym) documents
expected-answer checks for float/bool/string and optional sandboxed Python tools;
these check answer agreement, not the correctness of every reasoning step.

## Resources ranked by usefulness for our observed failures

| Resource | Proposed use | Qualification limit |
| --- | --- | --- |
| Existing frozen project index plus explicit state/version metadata | Retrieve exact current facts; distinguish original direct-path US from different protonation/NEB trials; retain original primary evidence | The 20 reviewed responses and answer-cued lookups are development data, not fresh held-out cases |
| [Llama Nemotron Embed 1B v2](https://huggingface.co/nvidia/llama-nemotron-embed-1b-v2) and [Rerank 1B v2](https://huggingface.co/nvidia/llama-nemotron-rerank-1b-v2) | Candidate local dense retrieval plus cross-encoder reranking alongside lexical search | Cards identify Meta Llama 3.2 1B ancestry; check model-specific organizational approval/terms. Compare retrieval recall and wrong-version selection against our existing lexical baseline before adoption |
| [PaperQA2](https://github.com/Future-House/paper-qa) | Optional literature component for local papers, evidence summaries, citations and contradictory claims; complementary to Paperclip search/download | Local models/embeddings can be configured; hosted defaults must be explicitly replaced. Retain source passages and inspect claims. It is not an enzyme-mechanism validator |
| [ChemBench](https://github.com/lamalab-org/chembench) | Broader chemistry diagnostic beyond property/ADME questions | Public benchmark scores cannot substitute for project-specific QM/MM and enzymology checks; training overlap may be unknown |
| [NeMo Gym](https://github.com/NVIDIA-NeMo/Gym) | Reusable evaluation/verifier patterns; later, small tool-grounded scientific tasks | Start with evaluation only. Installing an RL framework or post-training the 120B model is not the immediate priority |
| [NVIDIA OO Agents / NOOA](https://github.com/NVIDIA-NeMo/labs-OO-Agents) | Typed specialist interfaces, deterministic Python methods, agent composition and tracing | Research framework; CodeAct can execute generated Python. Evaluate only through bounded adapters with explicit PI checkpoints; do not replace the existing controller now |
| [NeMo Guardrails](https://github.com/NVIDIA-NeMo/Guardrails) | Optional checks around tool arguments/results and retrieved content, using local rules where possible | Useful policy enforcement, but grounding/safety checks do not prove proton balance, TS order or PMF convergence. Existing deterministic scientific checks and PI authority remain |

Both NOOA and Guardrails advertise Apache-2.0 licensing. Optional models,
dependencies and services have their own terms and approval requirements.
No Chinese-origin model or derivative is recommended here. Framework choice,
retrieval improvements and general benchmark success do not by themselves
improve intrinsic mechanistic reasoning.

## RFD4-Proteina: promising optional design arm, pending Nyl12 qualification

The [model card](https://huggingface.co/nvidia/RFD4-Proteina) describes joint
sequence/all-atom structure generation with ligand and catalytic-motif context.
Its current release pairs the 2026-10-19 flow model and autoencoder; pin both
together. The code/weights are BSD-3-Clause; optional refolding assets have
separate terms. Its published computational-pass statistics are not catalytic
activity measurements, and the card says they describe an earlier checkpoint.

[Partial diffusion](https://rfd4-proteina.org/latest/how-to/partial-diffusion.html)
can remodel an existing structure, which addresses a real limitation of
fixed-backbone LigandMPNN. Inverse folding is also supported. However, the
[current annotations reference](https://rfd4-proteina.org/latest/reference/annotations.html)
still lists **symmetry as planned**. A symmetry implementation was proposed in
[PR 98](https://github.com/RosettaCommons/RFD4-Proteina/pull/98), which is closed;
that does not establish released support. Do not assume fourfold sequence tying
or geometry symmetry works because multichain inputs can be supplied.

[Conditioning documentation](https://rfd4-proteina.org/latest/concepts/conditioning.html)
says ligands/modified residues are supplied atomized context; only the 20
standard amino acids are generated. Formal charges and discrete bond orders
are explicit inputs. This does not establish support for transferring protons
or fractional transition-state bonds. A barrier-region heavy-atom motif would
be a design constraint, not an electronic model or verified TS. Custom PA6,
covalent Thr adducts and mature-chain termini require a lossless mapping audit.

A future bounded comparison should:

1. Retain the matched symmetry-tied LigandMPNN baseline and qualified reaction
   ensemble; do not postpone PMF convergence for RFD4 setup.
2. Test restrained local backbone remodeling around a small selection of
   reaction-region geometries, preserving catalytic identities, chain breaks,
   interfaces and the maturation region. Input/output checks must demonstrate
   these properties; model conditioning alone is not acceptance.
3. Establish a reproducible route to identical sequences across all four
   protomers. If symmetry cannot be enforced and independently checked, RFD4
   is unsuitable as the primary generator for this campaign. Symmetry-aware
   LigandMPNN on an accepted remodeled backbone remains a possible alternative.
4. Validate the full mature assembly and precursor plausibility; refold with
   Boltz/AF3, relax, compare all four sites and assess multiple catalytic states.
   Preserve WT and known-active controls. Structural pass rates do not prove
   an activity improvement; native-Thr similarity remains diagnostic.

Generation and scoring environments must be isolated from working Amber/QM/MM
environments. The [installation guide](https://rfd4-proteina.org/latest/installation.html)
offers CUDA 12/13 profiles and documents qualified Linux x86_64 operation.
Neither a small weight bundle nor listed GPU support establishes memory/runtime
feasibility for the full Nyl12 assembly. No install or resource pilot occurred.

## Next decision

First compare source retrieval on fresh, PI-reviewed project questions with
deliberate historical/current distractors. Measure correct primary-source
retrieval, wrong-state/version errors, numerical accuracy, abstention and
runtime separately. Use existing tools to verify charge/proton inventory,
atom mapping, units and stationary-point/sampling claims. Consider small
domain post-training only if these interventions leave persistent failures
and a reliable disjoint training/evaluation set exists. RFD4 remains an
independent redesign-method pilot rather than a prerequisite for mechanism-agent.
