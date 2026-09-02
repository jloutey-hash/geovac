# Multi-Observable Focal-Length Decomposition Program (retired directive)

> **Status: DORMANT, retired from CLAUDE.md §1.8 on 2026-09-01.**
>
> Queued 2026-05-09 as a PI directive with five named deliverables (Paper 34
> §V.C.2 through §V.C.6). As of 2026-09-01, **none of the five exists in
> Paper 34** — measured by grep, not assumed. The thread has a documented
> successor: the v4.71.0 chemistry-error projection sprint describes itself as
> "the chemistry-solver analog of the §1.8 focal-length program".
>
> It was moved here during the CLAUDE.md compaction because a dormant
> directive should not occupy 6.5 KB of a file loaded into every session and
> every sub-agent dispatch. **The text below is the directive verbatim.**
> Reviving it is a PI call: if revived, move it back to §1.8 or state here
> which of the five targets are actually in scope.

---

## 1.8. Multi-Observable Focal-Length Decomposition Program (Directive)

The May 9 r_Z extraction thread (four iterations, fully documented in §2 below) demonstrated that the framework's projection-chain machinery functions as a unique diagnostic instrument for precision physics: **a multi-observable global fit using a single structural framework exposes convention-mismatches in literature compilations and kernel approximation gaps that are invisible from any single-observable analysis.** This subsection codifies that finding as an ongoing research program with explicit targets.

**The directive.** When a sprint touches precision atomic, molecular, or nuclear observables, the PM should ask three questions in order:

1. **Is there a literature convention mismatch the framework can surface?** Different compilations (Eides, Karshenboim, Krauth, Pachucki–Yerokhin, Drake, Antognini) itemize Layer-2 corrections at sub-percent levels with different conventions for which sub-leading terms count as "Zemach" vs "polarizability" vs "recoil NLO" vs "multi-loop QED." Multi-observable global fits using GeoVac's single projection-chain dictionary expose these mismatches at numerical level (the W1b finding, 2026-05-09: ~0.01% of $\nu_F$ propagates to ~25 mfm in extracted $r_Z$, below the 30 mfm Eides-vs-lattice gap). Frameworks that consume single-observable Layer-2 inputs in their native conventions cannot see this.

2. **Is there a GeoVac kernel approximation gap the framework can identify?** Tightening precision via more observables tends to relocate the load-bearing systematic (the W1a-D and W1b iterations: cross-register V_eN was correct; kernel-leading-order at coarse precision; sub-leading recoil-mixing at 51 mfm; convention mismatch at 16 mfm). Each precision tightening uncovers a new systematic. PMs run per-observable consistency checks before accepting any global-fit result.

3. **Does the observable admit a focal-length decomposition that sharpens the §III dictionary's coverage?** Multi-component observables (Lamb shift, hyperfine, fine structure) deserve Roothaan autopsies in Paper 34 §V.C. Single-component observables get single rows in §V/§V.B. The autopsy format makes the projection-chain decomposition visible at the observable level.

**Three classes of problem the program targets.**

- **Literature convention mismatches**: differences in Layer-2 itemization between Eides 2024 / Krauth 2017 / Karshenboim 2005 / Pachucki–Yerokhin 2010 / Drake 1990 / Antognini 2013 compilations. Framework's multi-observable fits expose these structurally.
- **GeoVac kernel approximation gaps**: places where the framework's leading-order operator (Eides Zemach, Foldy–Friar contact, Friar moment) is too coarse for the precision target. Each gap is a named follow-up (W1a, W1b, etc., per CLAUDE.md §1.7 multi-focal-composition wall taxonomy).
- **General focal-length decomposition cataloguing**: §V.C Roothaan autopsies for multi-component observables. Discipline: every Layer-2 input gets a focal-length tag, every framework-native contribution gets a §III chain, the decomposition closes at sub-MHz (or sub-kHz) precision.

**Active Roothaan-decomposition targets** (queued 2026-05-09):

1. **§V.C.2 Hydrogen 21 cm hyperfine four-component autopsy** (placeholder fill). Bohr-Fermi + Schwinger $a_e$ + reduced-mass + Zemach decomposition at +18 ppm framework-residual. Tests §III.18 magnetization-density at the operator level.

2. **§V.C.3 Muonic hydrogen 2S–2P Lamb shift autopsy** (placeholder fill). Decomposition into full Uehling kernel (Antognini-style) + SE Bethe-log + Friar moment via §III.17 + Källén–Sabry two-loop VP + deuteron polarizability. Tests §III.17 + §III.18 + §III.16 in the muonic regime ($\beta = 1.475$).

3. **§V.C.4 Helium $2{}^3P$ fine structure full autopsy** (NEW). Decomposition into spin-orbit + spin-spin + spin-other-orbit (Drake combining coefficients) using bipolar harmonic $(k_1, k_2)$ decomposition. First operator-level test of Observation 3's angular compositional projection rule (internal multi-focal at $\alpha^2$). Reference: NIST + Pachucki–Yerokhin theory at sub-ppm.

4. **§V.C.5 Helium $2{}^1P \to 1{}^1S$ oscillator strength** (NEW). Multi-electron extension of Sprint Calc-L (Lyman α, +0.055% match). Tests vector-photon promotion + Wigner 3j + multi-electron correlation. Reference: Drake handbook $f \approx 0.276$. Verifies whether Sturmian's continuum-closing property (validated for hydrogen polarizability at exact at $N_\text{basis} = 2$, Sprint Calc-P) extends to multi-electron transitions.

5. **§V.C.6 Cesium $6S_{1/2}$ hyperfine (atomic clock)** (NEW, prospective). $Z = 55$ heavy-atom regime where relativistic spinor lift dominates and §III.17/§III.18 are tested at very different focal lengths than hydrogen. Atomic-clock-grade reference precision ($10^{-15}$ relative; the SI second is defined by this transition). Tests Z-scaling of the framework's projection-chain machinery in a regime where simple hydrogenic formulas break down. **Prospective angle**: heavy-atom Cs PNC (parity non-conservation) extractions have known atomic-structure uncertainties at the percent level; the framework's projection-chain decomposition might illuminate atomic-structure systematics that single-observable analyses don't expose.

These five targets together exercise §III.17, §III.18, §III.19, spinor lift (§III.7), Wigner 3j (§III.8), Wigner $D$ (§III.9), Hopf bundle (§III.2), vector-photon promotion (§III.11), spectral action (§III.6), Sturmian (§III.5), and rest-mass (§III.14) — most of the projection dictionary. Successful decomposition of all five at sub-percent framework-native precision would validate the directive's broader applicability across atomic precision physics. Where each lands also tests the directive's three problem-classes empirically: how often does multi-observable consistency reveal a literature convention issue vs a framework kernel gap vs neither?

**Sprint cadence.** 2 tracks per sprint (1 verification + 1 prospective), established in the May 9 thread. Continue this cadence for the new targets.

**Cross-references.** May 9 r_Z thread (§2 below): demonstrated all three problem classes in a single observable across four iterations. Paper 34 §V.C.1 (Lamb shift autopsy): founding instance of the cataloguing discipline. Audit memo `debug/paper34_v_base_unit_audit_memo.md`: methodological precedent for systematic analysis. Diagnostic-before-engineering memory `feedback_diagnostic_before_engineering.md`: the load-bearing discipline that made the trajectory cleanly informative.

---
