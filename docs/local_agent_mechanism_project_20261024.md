# Proposed local agents for the Nyl12 PA6 mechanism project

2026-10-24. **Standing constraint: no Chinese-origin LLMs or derivatives.**
Model provenance must satisfy that constraint before selection. Concept only:
no model downloads, installations, dependencies or scientific jobs authorized.

A local team can organize evidence and run approved reproducible workflows
alongside the staged DFTB QM/MM plan. It cannot guarantee a complete mechanism.
Personas are semantic roles, not independently qualified scientists; agreement
is not scientific evidence. Protonation, connectivity, QM-region/embedding,
method validation and sampling convergence remain the substantive work.

## Hardware and candidate models

Two RTX A6000 cards provide 48 GB ECC GDDR6 each, with optional two-way NVLink.
Bridge presence is unverified. Installed 96 GB is not automatically one pooled
allocation; software must shard weights/cache across devices.
[Retained official NVIDIA specifications](https://www.nvidia.com/en-us/products/workstations/rtx-a6000/)

| Candidate | Proposed use and conditions |
|---|---|
| **gpt-oss-20b** | Routine evidence extraction, bounded planning and scripting baseline: 21B total/3.6B active parameters; Apache-2.0; 131,072-token advertised context. First qualify actual Ampere runtime support. |
| **gpt-oss-120b** | Conditional higher-capacity scientific reviewer: 117B total/5.1B active parameters, Apache-2.0 and 131,072 context. A two-GPU quantized trial may be feasible; no measured local fit/speed claim. Same family as 20b, so not independent corroboration. |

[Official 20b documentation](https://developers.openai.com/api/docs/models/gpt-oss-20b),
[official 120b documentation](https://developers.openai.com/api/docs/models/gpt-oss-120b).
Neither establishes Nyl12 chemistry competence.

The official OpenAI local guide quotes about 16 GB for 20b with MXFP4 and
at least 60 GB or multiple GPUs for 120b, but its stated MXFP4 fast path requires
Hopper or newer; A6000 is Ampere. Do not transfer those memory/speed claims
unconditionally. The guide estimates approximately 48 GB for 20b BF16, leaving
inadequate single-card headroom. Before deployment, prove a pinned compatible
quantized backend or explicitly choose a two-card BF16 configuration; never
silently substitute precision/offload. Use the required Harmony chat template.
[Official OpenAI local execution guide](https://developers.openai.com/cookbook/articles/gpt-oss/run-transformers)

**Recommended initial European baseline: Mistral-Small-3.2-24B-Instruct-2506.**
Its official Apache-2.0 card specifies approximately 55 GB BF16/FP16 GPU memory
and a two-way tensor-parallel serving example. It inherits the 3.1 family's
128k-context features. This offers a plausible two-A6000 trial without relying
on a third-party quantization; real context/cache/latency still require testing.
Preserve official weights and Mistral tokenizer/tool templates. Do not assume
an unofficial derivative satisfies the origin constraint.
[Official 3.2 card](https://huggingface.co/mistralai/Mistral-Small-3.2-24B-Instruct-2506),
[inherited 3.1 features](https://huggingface.co/mistralai/Mistral-Small-3.1-24B-Instruct-2503).

Start with one shared server, sequential role requests and an 8–16k total
input/output budget; measure quality as well as memory. Multiple processes
can duplicate weights; concurrent contexts consume cache. Advertised context
is not a fit guarantee. Switch larger reviewers sequentially and preserve GPU
capacity for science. Record per-device peak VRAM, latency and tool-call errors.

## Roles and gates

| Role | Required artifact |
|---|---|
| Enzymology | Product constraints, atom/proton accounting and competing relay hypotheses. |
| Literature | DOI/page/figure evidence ledger; distinguish mature catalysis, precursor maturation, experiment and calculation. |
| QM/MM | Explicit method, parameters, QM boundary and controlled validation proposal. |
| PI/planner | Next discriminating calculation, dependencies and resource/stop limits; human scientific signoff. |
| Critic | Original-artifact checks, counterexamples and convergence challenges before reading the planner conclusion. |

1. Qualify read-only against known traps: invalid PA66, AF3/Boltz atom naming,
   frame/time offsets and proton/charge inconsistencies. Reject invented
   citations or silent substitutions; compare against a one-role baseline.
2. Human review fixes products, mapping, competing mechanisms and initial QM
   region. Include acylation, leaving-group protonation, deacylation and enzyme
   regeneration without assuming stable intermediates or one universal relay.
3. Planner proposes declarative jobs; deterministic preflight checks units,
   charge/topology and budget. Executor uses approved scripts and fresh tmux
   jobs. Reviewer separates process completion from scientific acceptance.
4. Retain engine regressions, representative DFTB-versus-hybrid-DFT chemical
   checks and region/boundary sensitivity before enzyme campaigns. Retain
   overlap, bidirectional sampling, replicate uncertainty and hidden-water/proton
   checks before larger PMFs. PMF maxima do not automatically establish PES
   saddles. Human review approves scientific claims and campaign expansion.

Prefer existing scripts plus a deterministic planner–executor–reviewer state
machine. Candidate skills: evidence retrieval, atom mapping, proton bookkeeping,
QM/MM preflight, scientific-tmux provenance and PMF diagnostics. Preserve model,
quantization/kernel hashes, prompts, evidence, environments, inputs and failures.
Start with one reviewed evidence ledger and correct next-job proposal. Continue
only if measured correctness or researcher time improves.

Apache-2.0 model licensing does not grant rights to restricted papers, DFTB
parameters or simulation engines; preserve applicable redistribution notices.
Ask before new production dependencies and use an approved isolated `uv`
environment if local deployment is later authorized. No external model API is
needed for this proposal.
