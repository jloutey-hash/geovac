# L3b-2 Preflight Concurrent Work Memo

**Date:** 2026-05-17
**Scope:** Fresh literature check before committing to original math.OA paper claiming Lorentzian propinquity convergence on truncated SU(2) × U(1)_T Krein spectral triples (weak-form, K⁺ state space, Wasserstein–Kantorovich distance, explicit rate).
**Window:** 2025-06-01 → 2026-05-17 (two months past prior May 2026 audit at `debug/l3_literature_audit_memo.md`).

## Summary verdict

**CLEAR — no direct overlap; one moderate adjacency to watch (Mondino–Sämann line on the Riemannian-geometry side, still NOT spectral triples).**

## Per-search results

**(1) "Lorentzian propinquity"** — Zero hits. Term not yet in the literature. Closest hits are Nieuviarts arXiv:2512.15450 (Dec 2025) and arXiv:2502.18105 (May 2025), both Riemannian → pseudo-Riemannian *twist morphisms* on twisted spectral triples — algebraic, NOT distance/convergence. **NO_OVERLAP.**

**(2) Latrémolière post-2025-06** — `arXiv:2504.11715` (Apr 2025, pre-window-edge), "Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics" — strictly Riemannian, analytic-path continuity result. No Krein/Lorentzian extension found. **NO_OVERLAP.**

**(3) Mondino–Sämann** — `arXiv:2504.10380` (Apr 2025), "Lorentzian Gromov–Hausdorff convergence and pre-compactness" — extends their `2209.14384` program (LMP 2024). Lives in Lorentzian pre-length spaces / synthetic metric-measure geometry via ε-nets of causal diamonds. **Does NOT bridge to spectral triples or NCG.** Adjacent at the conceptual level (both are "Lorentzian GH-type convergence"); structurally disjoint (synthetic geometry vs operator algebra). `arXiv:2510.13069` ("Gromov's Compactness Theorem for the Intrinsic Timed-Hausdorff Distance", Oct 2025) is a downstream pure-metric extension. **ADJACENT — flag in related-work; not a scoop.**

**(4) "Krein spectral triple convergence"** — No GH/propinquity-style convergence theorems for Krein spectral triples found. `arXiv:2402.05839` (Mar 2026 update, "Signature change by a morphism of spectral triples") is structural, not metric. **NO_OVERLAP.**

**(5) "noncommutative Lorentzian geometry" + convergence** — `arXiv:2404.00240` ("Collapse in Noncommutative Geometry and Spectral Continuity", Oct 2025) discusses convergence/approximation language but in standard Riemannian NCG; `arXiv:2601.07350` (Jan 2026) is algebra-construction of a novel NC spacetime, not convergence. **NO_OVERLAP.**

**(6) Lorentzian-NCG lineage (vdD, BBB, Franco–Eckstein, Devastato–Lizzi–Martinetti, Strohmaier)** — `arXiv:2601.22171` (Jan 2026, "Pseudo-Riemannian Spectral Triples for SU(1,1)") verifies a specific triple satisfies vdD–Paschke–Rennie + vdD–Rennie axioms; structural, no convergence. van den Dungen 2025 papers (Dirac–Schrödinger / Callias) are Riemannian index theory. `arXiv:2503.24192` (κ-Minkowski quantum causality) twisted-Lorentzian-triple, no propinquity. **NO_OVERLAP.**

**(7) Hekkelman / McDonald / Leimbach / vS** — Hekkelman–McDonald `arXiv:2412.00628` and Leimbach–vS (tori) remain Riemannian; no Lorentzian extensions found. **NO_OVERLAP.**

**(8) Bisognano–Wichmann + spectral triple + convergence** — `arXiv:2505.10530` (May 2025, lattice BW convergence) is condensed-matter / quant-ph, not spectral triple. No BW + spectral triple + convergence work found in window. **NO_OVERLAP.**

## Recommendations

**PROCEED with Paper 45 (L3b-2 propinquity convergence theorem).** No direct overlap; the Lorentzian-propinquity term is empty in the literature.

Required adjustments to the draft:
1. **Cite Mondino–Sämann `2209.14384` (LMP 2024) + `2504.10380` (Apr 2025)** in §Related Work, with one paragraph distinguishing: synthetic Lorentzian GH on pre-length spaces vs operator-algebra propinquity on Krein spectral triples — different objects, different distances, no implication either way.
2. **Cite Nieuviarts `2512.15450` + `2502.18105`** as the algebraic Riemannian→pseudo-Riemannian alternative (twist morphism, not metric); reiterate the odd-dim obstruction Paper 42 §3.2 already documents.
3. **Cite vdD–Rennie / vdD–Paschke–Rennie via `2601.22171`** to anchor the Krein-spectral-triple definition used.
4. **Cite Latrémolière `2504.11715`** as the most recent Riemannian propinquity continuity result and the natural Riemannian companion the L3b-2 theorem extends.

No scope modifications needed. No abort. **Window is open.**

## Sources

- [arXiv:2512.15450 — Nieuviarts (Dec 2025/May 2026)](https://arxiv.org/abs/2512.15450)
- [arXiv:2502.18105 — Nieuviarts (May 2025)](https://arxiv.org/abs/2502.18105)
- [arXiv:2504.11715 — Latrémolière (Apr 2025)](https://arxiv.org/abs/2504.11715)
- [arXiv:2504.10380 — Mondino–Sämann (Apr 2025)](https://arxiv.org/abs/2504.10380)
- [arXiv:2510.13069 — Intrinsic Timed-Hausdorff (Oct 2025)](https://arxiv.org/abs/2510.13069)
- [arXiv:2402.05839 — Signature change morphism (Mar 2026)](https://arxiv.org/abs/2402.05839)
- [arXiv:2404.00240 — Collapse in NCG (Oct 2025)](https://arxiv.org/abs/2404.00240)
- [arXiv:2601.07350 — Novel NC spacetime (Jan 2026)](https://arxiv.org/abs/2601.07350)
- [arXiv:2601.22171 — Pseudo-Riemannian ST for SU(1,1) (Jan 2026)](https://arxiv.org/abs/2601.22171)
- [arXiv:2412.00628 — Hekkelman–McDonald NC integral (Dec 2024)](https://arxiv.org/abs/2412.00628)
- [arXiv:2503.24192 — κ-Minkowski causality (Mar 2025)](https://arxiv.org/abs/2503.24192)
- [arXiv:2505.10530 — Lattice BW (May 2025)](https://arxiv.org/abs/2505.10530)
