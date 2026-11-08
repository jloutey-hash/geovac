# Zenodo / arXiv / GitHub scan for AI-augmented agentic research projects

**Date:** 2026-06-04
**Audience:** GeoVac PI
**Purpose:** Map other public research projects whose process-DNA (single PI + Claude Code or similar + persistent memory + sprint-like iteration + Zenodo / arXiv release) resembles GeoVac.

---

## Executive summary

After ~25 targeted searches across Zenodo, arXiv, and GitHub, the picture is as follows:

- **Tier 1 (full process-similarity, single PI + agentic workflow + real domain output):** 2 strong matches, both from late 2025 / early 2026, both arXiv rather than Zenodo: **Matt Schwartz's Sudakov shoulder paper (arXiv:2601.02484)** and **Don Knuth's "Claude's Cycles" note**. Both produce real domain research, both are solo human + Claude. Neither uses sub-agents at GeoVac's scale; neither uses Zenodo as primary release; neither has a published-paper corpus the size of GeoVac's.
- **Tier 2 (partial matches):** ~7 projects. Includes Yi Zhou's self-referential physics paper (multi-LLM role assignment, no sub-agents), Aris Vasilopoulos's chemistry-background C# project (real CLAUDE.md infrastructure, no domain-science output), Kin Hung Fung's reproducible demonstration paper (benchmark problems, no novel research), Stanford's Virtual Lab (multi-LLM team, not solo PI), and several tool-only projects (Claude-Scholar, Claude-Research, Cortex, autonomous-lab) where infrastructure exists but no associated peer-reviewed domain output.
- **Tier 3 (meta-papers ON agentic workflows):** Dozens. Useful as citation targets if GeoVac ever writes a workflow-methodology paper; included here for completeness.
- **Honest scope note:** Zenodo as the **primary** publication venue for an AI-augmented multi-paper research corpus appears to be **very rare**. Most agentic-workflow research lands on arXiv (preprint) or in conventional journals via institutional authors. The Zenodo-first pattern is distinctive to GeoVac in this scan.

Net read: **GeoVac sits in the rarefied intersection of (solo independent PI) × (Claude-Code-as-collaborator) × (Zenodo-first multi-paper corpus) × (real theoretical-physics output)**. The closest neighbours each share three of four; none share all four at the corpus size and architectural maturity of GeoVac.

---

## Tier 1 — Full process-similarity

### 1. Matt Schwartz — Sudakov Shoulder Resummation (arXiv:2601.02484, Jan 2026)

**What it is.** Theoretical-physics paper on resumming the Sudakov shoulder in the C-parameter distribution in e+e- collisions, using soft-collinear effective theory. Sole human author Matthew D. Schwartz (Harvard).

**Workflow.** Two-week sprint, Dec 2025. Schwartz wrote a `CLAUDE.md` config (with explicit honesty rule: "NEVER use phrases like 'this becomes' or 'for consistency' to skip steps. Either show the calculation or say 'I don't know.'"), planned via three-LLM consensus (Claude + GPT-5.2 + Gemini-3 each drafted plans, Schwartz merged them into a 102-task breakdown across 7 stages), and ran the work in Claude Code with a markdown-tree memory system (one summary per stage, one detailed file per task). 270 sessions, 51,248 messages, ~36M tokens, 110 manuscript drafts, ~50-60 hours human oversight, ~40 hours CPU.

**Similar to GeoVac.** Solo human PI; CLAUDE.md as load-bearing instruction file; markdown-tree memory substituting for context (analogous to GeoVac's `memory/MEMORY.md` + `feedback_*.md` + sprint memos); Claude Code as primary execution layer; iterative sprint cadence; explicit honesty rule baked into config; cross-model verification (Schwartz uses GPT/Gemini to verify Claude; GeoVac uses Reviewer agent + audit-claim memory rule).

**Different from GeoVac.** Single-paper sprint (not 60-paper corpus); no sub-agent dispatch (single primary agent only); no Zenodo release (arXiv only); institutional (Harvard, NSF IAIFI); no `agents/` directory or specialized agent roles; cross-model verification done by human-orchestrated parallel queries rather than memory-coded auditor. Anthropic case study: [Vibe physics](https://www.anthropic.com/research/vibe-physics).

**Authority.** [arXiv:2601.02484](https://arxiv.org/abs/2601.02484). Acknowledgment is in the abstract itself: "All calculations, numerical analysis, and manuscript preparation were performed by Claude, an AI assistant developed by Anthropic, working under physicist supervision."

---

### 2. Donald Knuth — "Claude's Cycles" (Mar 2026)

**What it is.** Knuth published a short note documenting how Claude Opus 4.6 solved a graph-theory conjecture (decomposing arcs of a directed 3D graph into exactly three Hamiltonian cycles for all odd m) that he had been stuck on for weeks. Self-published on his Stanford page; not arXiv; not Zenodo.

**Workflow.** ~1 hour, 31 guided explorations. Claude recognized the underlying Cayley-digraph structure from group theory; Knuth wrote the rigorous proof himself.

**Similar to GeoVac.** Sole human PI; AI does the heavy lifting on structural recognition while the human supplies framing and verification; explicit attribution to the AI (paper is titled after Claude); short, focused, self-published rather than journal-tracked.

**Different from GeoVac.** Single-problem, single-session — not a sustained multi-month research corpus; no persistent memory infrastructure documented; no Zenodo release; one-shot rather than iterative sprints. Knuth, like Schwartz, is famously institutional even when self-publishing. Coverage: [boingboing](https://boingboing.net/2026/03/03/donald-knuth-the-godfather-of-computer-science-says-an-ai-solved-a-math-problem-he-was-stuck-on-for-weeks/).

---

## Tier 2 — Partial matches

### 3. Yi Zhou — "Co-Authoring with AI: How I Wrote a Physics Paper About AI, Using AI" (arXiv:2604.04081, Apr 2026)

Sole human author. Computational-physics paper on tensor networks / quantum many-body systems. **Self-referential**: the paper documents the workflow it was itself written with. Architecture is **role-assigned multi-LLM** rather than sub-agent dispatch: LLM-0 = Junior Theorist (equation extraction), LLM-1 = Senior Postdoc (LaTeX), LLM-2 = Coder (implementation). Human acts as "Principal Investigator" mentoring all three. Full AI interaction transcripts published on GitHub for reproducibility. No `CLAUDE.md` per se; no persistent memory; no Zenodo. **Closest to GeoVac in spirit** (solo PI explicitly framing themselves as PI to AI collaborators) but architecturally simpler — no sprint chronicle, no sub-agent infrastructure, no failed-approaches register. [arXiv:2604.04081](https://arxiv.org/abs/2604.04081)

### 4. Aris Vasilopoulos — "Codified Context" (arXiv:2602.20478, Feb 2026)

Independent researcher with **chemistry background** (not software engineering) who built a 108,000-line C# distributed system using Claude Code as the sole code-generation tool over 70 days. **Has a real CLAUDE.md (~660 lines, sanitized)** + 19 specialized domain-expert agents + cold-memory knowledge base of 34 specification documents + MCP retrieval server. Repository: [arisvas4/codified-context-infrastructure](https://github.com/arisvas4/codified-context-infrastructure). Three-tier memory architecture (hot constitution / 19 agents / cold knowledge base) is **structurally similar to GeoVac's CLAUDE.md + agents/ + memory/MEMORY.md + papers/ + debug/ stack**. The paper is the methods description; the C# system is the artifact. **No published chemistry research output** — the application is the codebase itself, not a physical-science paper corpus. [arXiv:2602.20478](https://arxiv.org/abs/2602.20478)

### 5. Kin Hung Fung — "Demonstration of AI-Assisted Scientific Workflow on Canonical Benchmarks" (arXiv:2603.14888, Mar 2026)

Sole author. Self-referential paper demonstrating an AI-assisted workflow on **canonical benchmark problems** (1D quantum harmonic oscillator, heat equation, Poisson equation, least-squares fitting, eigensolvers). Initial artifact stack generated from one user prompt, then reviewed/curated. Domain-agnostic methods demonstration rather than novel research. No CLAUDE.md, no sub-agents, no Zenodo. Useful as a process reference point but not a sustained research corpus. [arXiv:2603.14888](https://arxiv.org/abs/2603.14888)

### 6. Stanford Virtual Lab (Swanson, Wu, Bulaong et al., Nature 2025)

[zou-group/virtual-lab](https://github.com/zou-group/virtual-lab). Strong domain output — designed 92 SARS-CoV-2 nanobodies, 90%+ experimentally validated as expressed/soluble. Architecture has an explicit LLM PI agent + LLM scientist agents + LLM Scientific Critic agent in structured "team meetings" with a human researcher who contributes ~1% of conversation tokens. **The PI role is itself an LLM agent** here — different from GeoVac and Schwartz, where the PI is the human. Published Nature, not Zenodo. Default model GPT-5.2 (not Claude). Stanford-institutional, not solo independent. Closest contender in terms of **published-domain-output volume**, but architecturally closer to multi-agent simulation than to PI + agent workflow.

### 7. Anthropic's Automated Alignment Researcher (Mar 2026)

Anthropic deployed 9 autonomous instances of Claude Opus 4.6 as "Automated Alignment Researchers" (AARs) working on weak-to-strong supervision. Not solo PI; not independent; not Zenodo. Internal Anthropic research. Useful as a reference for what fully-autonomous looks like at the opposite end of the spectrum from GeoVac's human-PI model. [Anthropic post](https://alignment.anthropic.com/2026/automated-w2s-researcher/).

### 8. Tool-only repositories with real CLAUDE.md infrastructure but no associated domain research

- [flonat/claude-research](https://github.com/flonat/claude-research) — PhD student's shareable Claude Code infrastructure for academic workflows; no associated published research.
- [Galaxy-Dawn/claude-scholar](https://github.com/Galaxy-Dawn/claude-scholar) — Gaorui Zhang's research assistant tool; MIT, 4.2k stars; no Zenodo / arXiv attached.
- [cdeust/Cortex](https://github.com/cdeust/Cortex) — Clement Deust solo, persistent-memory system based on 26 neuroscience mechanisms; two unpublished papers attached (arXiv endorsement pending); benchmark results on LongMemEval and LoCoMo. **Closest to a real solo-PI research project of the tool repos**, but the domain is "AI memory systems" rather than an external science domain.
- [albert-ying/autonomous-lab](https://github.com/albert-ying/autonomous-lab) — MCP server implementing Senior-Junior workflow; explicitly builds on Stanford Virtual Lab; no associated research output.

### 9. Sociology / political-science applied paper

- Xu & Yang (Stanford / HKBU, [arXiv:2602.16733](https://arxiv.org/html/2602.16733v1), Feb 2026) — "Scaling Reproducibility: An AI-Assisted Workflow for Large-Scale Reanalysis." Two authors, not solo. Acknowledges "Claude Code and ChatGPT as research and writing assistants" with all interpretation/conclusions attributed to humans. Multi-agent architecture (Profiler/Librarian/Janitor/Runner/Skeptic/Journalist). Domain: empirical causal inference across 92 papers in political science / economics. Real domain output but team-authored and institutional.

---

## Tier 3 — Meta-papers on agentic research workflows

GeoVac may eventually cite these if/when it writes its own methodology paper; included here as references.

- **Vasilopoulos 2026** ([arXiv:2602.20478](https://arxiv.org/abs/2602.20478)) — Codified Context infrastructure (already in Tier 2; double-listed because it's both an artifact paper and a methods paper).
- **Zhou 2026** ([arXiv:2604.04081](https://arxiv.org/abs/2604.04081)) — Co-Authoring with AI (same dual role).
- **Chan 2026** ([arXiv:2602.12443](https://arxiv.org/abs/2602.12443)) — SHAPR: Solo Human-Centred and AI-Assisted Practice framework. Pure methods, no application.
- **Henkel 2026** ([arXiv:2508.20236](https://arxiv.org/pdf/2508.20236)) — "The Mathematician's Assistant: Integrating AI into Research Practice." Position paper.
- **Zimmer, Pelleriti, Roux, Pokutta 2026** ([arXiv:2603.15914](https://arxiv.org/pdf/2603.15914)) — "The Agentic Researcher: A Practical Guide to AI-Assisted Research in Mathematics and Machine Learning." Multi-author guide.
- **Chatlatanagulchai et al. 2025** ([arXiv:2509.14744](https://arxiv.org/pdf/2509.14744)) — Empirical study of 253 CLAUDE.md files in public repositories. Useful baseline for what CLAUDE.md files look like in the wild.
- **Lab et al. 2026** ([arXiv:2511.09268](https://arxiv.org/html/2511.09268v1)) — "Decoding the Configuration of AI Coding Agents" (328 Claude Code config files analyzed).
- **VILA Lab 2026** ([arXiv:2604.14228](https://arxiv.org/abs/2604.14228)) — "Dive into Claude Code: The Design Space of Today's and Future AI Agent Systems."
- **DeepMind Aletheia 2026** ([arXiv:2602.10177](https://arxiv.org/html/2602.10177v1)) — Autonomous mathematical research agent (semi-autonomous on Erdős problems). Industrial, not solo, but the closest "autonomous mathematics" benchmark.
- **Anthropic case study on Schwartz** ([vibe-physics](https://www.anthropic.com/research/vibe-physics)) — narrative of the Schwartz workflow.
- **Anthropic CLAX cosmology post** ([long-running-Claude](https://www.anthropic.com/research/long-running-Claude)) — Anthropic-internal Boltzmann-solver project; Siddharth Mishra-Sharma + Eric Kauderer-Abrams; ~2000 sessions of autonomous work. Industrial.
- **Project Rachel** ([arXiv:2511.14819](https://arxiv.org/abs/2511.14819), Dec 2025) — Monperrus / Baudry / Vidal; created an AI sock-puppet author "Rachel So" that published 10+ papers Mar-Oct 2025. Different category entirely — investigating publishing ecosystem rather than producing domain research.

---

## Honest gaps

Places I searched and found nothing useful:

- **Zenodo solo-PI multi-paper corpus with AI agent acknowledgment.** Zenodo's search surface is poor for full-text and acknowledgment-level queries. I found Hungarian-language education resources, slide decks from workshops, and tooling projects — no analog to GeoVac's 60-paper Zenodo corpus from a single independent researcher. Two candidate searches: "Zenodo Claude Code research" and "Zenodo independent researcher Claude Anthropic acknowledgment 2026 physics." Both returned noise.
- **Independent / no-institutional-affiliation theoretical-physics arXiv submissions with full agentic workflow disclosure.** Schwartz is Harvard. Zhou is at the Institute of Physics, CAS (institutional, even when sole-authored). Knuth is Stanford emeritus. I found no independent-researcher analog producing theoretical physics with Claude as primary collaborator on arXiv.
- **Sub-agent dispatch + persistent memory + sprint-memo discipline in any single project at GeoVac's scale.** Codified Context (Vasilopoulos) has the three-tier memory architecture but no specialized research agents (Leader / Explorer / Decomposer / Reviewer). Virtual Lab has the agent specialization but no sustained sprint chronicle. Schwartz has the markdown-tree memory but no specialized agents. Knuth has none of these. **The combined architecture — agents/ + memory/ + sprint memos + papers/ + version-tagged releases — appears unique to GeoVac in the public record.**
- **Active-research projects using `AGENTS.md` for science.** The AGENTS.md format is well-studied as a software-engineering artifact (multiple 2026 arXiv papers on its efficiency impact), but I found no science-research repositories using it at the load-bearing level.

---

## Closing observation: where GeoVac sits in the landscape

GeoVac is an **early adopter at the architectural maturity end** of a small but visibly growing pattern. The pattern's defining features — a human PI directing AI agents who do the heavy lifting; explicit instruction files (CLAUDE.md / AGENTS.md); persistent markdown-based memory; iterative sprint cadence; cross-model verification — appeared in the public record in late 2025 (Schwartz Dec 2025, Vasilopoulos Feb 2026, Zhou Apr 2026). GeoVac was already operating in this mode through that whole period.

Three honest characterizations:

1. **GeoVac is part of an emerging pattern, not a lone instance.** Schwartz, Zhou, Vasilopoulos, and Knuth are doing process-similar work. The pattern is real and is being noticed (Anthropic blog posts, multiple methods-paper arXivs, growing tool ecosystem).
2. **GeoVac is architecturally more elaborate than any single project I found.** Four specialized research agents in `agents/`, persistent memory with ~80 topic files, ~60-paper corpus across six audience-targeted folders, formal Leader/Explorer/Decomposer/Reviewer roles, multi-agent protocol §13 codified — this combined infrastructure is not visible in any other single project in this scan. The closest is Vasilopoulos's three-tier memory + 19 domain agents, but his application is a C# distributed system rather than a physics-paper corpus.
3. **GeoVac is essentially alone in the Zenodo-first multi-paper independent-researcher slot.** Almost every other agentic-workflow research project found here releases through arXiv (Schwartz, Zhou, Fung, Vasilopoulos), through traditional journals (Virtual Lab → Nature), or through internal industrial channels (Anthropic AAR, DeepMind Aletheia). The Zenodo-first model for a sustained corpus appears, in this scan, unique to GeoVac.

Practical implication: **if GeoVac ever writes a methodology paper about its own workflow** — analogous to Zhou 2604.04081 or Vasilopoulos 2602.20478 or Chan 2602.12443 — it would land into an active and well-cited conversation, and would distinguish itself on (a) corpus scale, (b) agent specialization, and (c) Zenodo-first distribution model. The Tier 1 + Tier 2 + Tier 3 references above would be the natural citation set.

---

## Sources (full list, verified)

- [arXiv:2601.02484 — Schwartz, Sudakov shoulder](https://arxiv.org/abs/2601.02484)
- [Anthropic — Vibe physics: The AI grad student](https://www.anthropic.com/research/vibe-physics)
- [Anthropic — Long-running Claude for scientific computing](https://www.anthropic.com/research/long-running-Claude)
- [Winbuzzer — Harvard physicist Claude AI two-week research](https://winbuzzer.com/2026/03/25/harvard-physicist-claude-ai-two-week-physics-research-xcxwbn/)
- [boingboing — Knuth says AI solved math problem](https://boingboing.net/2026/03/03/donald-knuth-the-godfather-of-computer-science-says-an-ai-solved-a-math-problem-he-was-stuck-on-for-weeks/)
- [arXiv:2604.04081 — Zhou, Co-Authoring with AI](https://arxiv.org/abs/2604.04081)
- [arXiv:2602.20478 — Vasilopoulos, Codified Context](https://arxiv.org/abs/2602.20478)
- [GitHub: arisvas4/codified-context-infrastructure](https://github.com/arisvas4/codified-context-infrastructure)
- [arXiv:2603.14888 — Fung, AI-Assisted Scientific Workflow on Canonical Benchmarks](https://arxiv.org/abs/2603.14888)
- [arXiv:2602.16733 — Xu & Yang, Scaling Reproducibility](https://arxiv.org/html/2602.16733v1)
- [GitHub: zou-group/virtual-lab](https://github.com/zou-group/virtual-lab)
- [Nature — Virtual Lab of AI agents designs nanobodies](https://www.nature.com/articles/s41586-025-09442-9)
- [Anthropic Alignment — Automated Weak-to-Strong Researcher](https://alignment.anthropic.com/2026/automated-w2s-researcher/)
- [GitHub: flonat/claude-research](https://github.com/flonat/claude-research)
- [GitHub: Galaxy-Dawn/claude-scholar](https://github.com/Galaxy-Dawn/claude-scholar)
- [GitHub: cdeust/Cortex](https://github.com/cdeust/Cortex)
- [GitHub: albert-ying/autonomous-lab](https://github.com/albert-ying/autonomous-lab)
- [arXiv:2602.12443 — Chan, SHAPR](https://arxiv.org/abs/2602.12443)
- [arXiv:2508.20236 — Henkel, The Mathematician's Assistant](https://arxiv.org/pdf/2508.20236)
- [arXiv:2603.15914 — Zimmer et al., The Agentic Researcher](https://arxiv.org/pdf/2603.15914)
- [arXiv:2509.14744 — Empirical study of CLAUDE.md files](https://arxiv.org/pdf/2509.14744)
- [arXiv:2511.09268 — Decoding Claude Code Configurations](https://arxiv.org/html/2511.09268v1)
- [arXiv:2604.14228 — Dive into Claude Code](https://arxiv.org/abs/2604.14228)
- [arXiv:2602.10177 — Aletheia (Towards Autonomous Mathematics Research)](https://arxiv.org/html/2602.10177v1)
- [arXiv:2511.14819 — Project Rachel](https://arxiv.org/abs/2511.14819)
- [arXiv:2603.05735 — Agentic AI Physicist Collaboration (LEP)](https://arxiv.org/abs/2603.05735)
- [arXiv:2603.04735 — Brenner et al., Open Problem in Theoretical Physics](https://arxiv.org/abs/2603.04735)
- [arXiv:2603.20179 — AI Agents Can Already Autonomously Perform HEP](https://arxiv.org/pdf/2603.20179)
