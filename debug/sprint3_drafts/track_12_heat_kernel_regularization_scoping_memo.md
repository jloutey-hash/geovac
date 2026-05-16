# Track 12 Scoping Memo: Heat-Kernel Regularization / Schwinger Proper-Time

**Sprint 3, Paper 34 dictionary-completion effort**
**Status:** scoping memo only — no §III subsection draft
**Date:** 2026-05-15
**Verdict (one-line):** **Absorbed into §III.6 with a cross-reference note to the master Mellin engine framing of Paper 32 §VIII / Paper 18 §III.7.** A standalone "Schwinger proper-time" projection would not add structural content beyond §III.6 plus the case-exhaustion theorem already cited from Paper 32, and would risk double-counting the same Mellin transform under two names.

---

## 1. Question

Schwinger proper-time integration is the operation
\[
   \mathrm{Tr}\, f(D/\Lambda) \;=\; \int_0^\infty dt\, f(t) \cdot \mathrm{Tr}\, e^{-t D^2},
\]
which, via the Mellin bridge
\(\zeta_A(s) = \frac{1}{\Gamma(s)} \int_0^\infty t^{s-1} \mathrm{Tr}(e^{-tA}) dt\),
converts the discrete spectrum of $D$ on a compact manifold into physical observables (heat-kernel coefficients, spectral zeta values, effective actions). This is the single mathematical operation that produces every M1 / M2 / M3 π-source in the framework, per the master Mellin engine sharpening of Paper 32 §VIII and Paper 18 §III.7. The case-exhaustion theorem (Paper 32) names the three sub-cases:
\[
  \pi\text{-source} = \mathcal{M}\!\left[\mathrm{Tr}(D^k \cdot e^{-tD^2})\right](s),
  \qquad k \in \{0, 1, 2\}.
\]

The scoping question: should this Mellin-transform engine sit as its own §III entry under a name like "Schwinger proper-time integration" or "Mellin / heat-kernel regularization," or is it correctly absorbed into the existing §III.6 (Connes–Chamseddine spectral action / heat kernel)?

The structural test: does the proper-time integration $\int_0^\infty dt\, f(t) \cdot \mathrm{Tr}\, e^{-tD^2}$ inject new physical variables, dimensions, or transcendental content beyond what §III.6 already encodes? If yes, it deserves its own slot. If no, it is a re-naming of §III.6's operational content.

---

## 2. Comparison to §III.6 spectral action

§III.6 (`sec:proj_spectral_action`) currently encodes:

- **Source:** Dirac operator $D$ on bare $S^3$ graph (or extension).
- **Target:** heat-kernel expansion $\mathrm{Tr}\, f(D^2/\Lambda^2)$ with cutoff function $f$ and scale $\Lambda^2$.
- **Variables introduced:** $\Lambda$ (UV cutoff), $\alpha$ (when gauge sector included).
- **Dimension introduced:** energy via $\Lambda$.
- **Transcendental signature:** $\sqrt{\pi} \times \mathbb{Q}$ Seeley–DeWitt coefficients on unit $S^3$ ($a_0 = a_1 = \sqrt{\pi}$, $a_2 = \sqrt{\pi}/8$); resulting observable coefficients are $\pi^{2k} \times \mathbb{Q}$ at one loop (T9 theorem, Paper 28).

Operationally, the Connes–Chamseddine spectral action $S_\mathrm{CC} = \mathrm{Tr}\,f(D/\Lambda)$ IS a Schwinger proper-time integral. The standard derivation reads $f(D/\Lambda)$ as an inverse Laplace transform $f(x) = \int_0^\infty d\tau\, \tilde{f}(\tau) e^{-\tau x^2}$ applied to the heat kernel $\mathrm{Tr}\, e^{-\tau D^2}$, and the Seeley–DeWitt expansion is the small-$\tau$ asymptotic of that heat trace. The target signature in §III.6 ($\sqrt{\pi}\cdot\mathbb{Q}$ from $a_k$ coefficients, $\pi^{2k}\cdot\mathbb{Q}$ from the Mellin output) is exactly what the Schwinger proper-time picture delivers when applied to $D^2$. There is no operational difference between "do the Connes–Chamseddine spectral action evaluation" and "do the Schwinger proper-time integration of $e^{-\tau D^2}$ against a test function $\tilde{f}$" — they are the same computation under two names from two different literatures.

The M2 sub-mechanism of the master Mellin engine is exactly §III.6's content: $k=2$, $\pi$-source from the Mellin transform of $\mathrm{Tr}(D^2 e^{-tD^2})$ evaluated on Seeley–DeWitt coefficients, producing $\sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2\cdot\mathbb{Q}$ ring elements per Sprint MR-B's closed-form modular residual (Paper 18 §III.7, eq. `dirac_modular_residual`). §III.6 already cites this lineage indirectly through "Paper 28 (QED on $S^3$), Paper 32 (spectral triple synthesis)," and the existence proof for the spectral-action / proper-time identification is implicit in those cross-references.

What §III.6 does *not* currently spell out is that M1 ($k=0$, Hopf-base measure, identified with §III.2/Paper 25) and M3 ($k=1$, vertex-parity Hurwitz, identified with the Camporesi–Higuchi spinor lift in §III.7 / Paper 28 §IV) are *also* outputs of the same Mellin-transform technology at different operator orders. A reader of §III.6 alone could be left thinking the heat-kernel expansion is specific to the spectral action; a reader of Paper 32 §VIII (case-exhaustion theorem) knows it is the master mechanism. A cross-reference inside §III.6 to the Paper 32 §VIII / Paper 18 §III.7 framing would close this gap without requiring a new §III entry.

---

## 3. Comparison to the master Mellin engine partition

Paper 32 §VIII's case-exhaustion theorem and Paper 18 §III.7's master-mechanism reading together establish that M1, M2, M3 are *three sub-cases of one mechanism* parameterized by operator order $k \in \{0, 1, 2\}$. The Mellin transform itself — the operation $\mathcal{M}[\cdot](s) = \frac{1}{\Gamma(s)} \int_0^\infty t^{s-1} \cdot (\cdot) \, dt$ applied to $\mathrm{Tr}(D^k e^{-tD^2})$ — is the *meta-mechanism* that produces each named §III projection's transcendental signature, not itself a §III projection.

The structural reading: §III projections in Paper 34 are the *named operations* that act on the bare graph and inject physical content. The Mellin transform is not such a named operation in the same sense; it is the *evaluator* that, given a §III projection (Hopf bundle for M1, spectral action / spinor lift for M2/M3), extracts the transcendental content the projection introduces. The bare graph plus Hopf bundle injects $\mathrm{Vol}(S^2)/4 = \pi$ via the M1 sub-mechanism; that "via M1" is the Mellin operation at $k=0$. But the projection is the Hopf bundle (§III.2), not the Mellin transform.

This reading is reinforced by §VIII apparatus identity (`sec:proj_apparatus_identity`), which Sprint TD Track 5 (May 2026) established as the *spectral-vs-state-side* divide: spectral-side projections all live inside the master Mellin engine $M_1 \cup M_2 \cup M_3$, while state-side reductions ($\mathrm{Tr}\,\rho\log\rho$) sit outside it. Paper 34's §III.28 (apparatus identity) treats the spectral-side ↔ state-side bridge as the projection. The Mellin engine itself is the *organizing principle* of all spectral-side §III entries, not a peer among them. Adding "Schwinger proper-time" as a §III peer would conflate the evaluator with the projections it evaluates.

The same logic applies to the analogy with Paper 18 §IV taxonomic tiers (intrinsic / calibration / embedding / etc.): these tiers are the *output classes* of the Mellin engine applied to different operator content, not themselves §III projections in Paper 34.

---

## 4. §X open-question (a) status

The §X open-question entry (a) — "the Mellin / heat-kernel evaluation at fractional $s$ (currently treated as part of spectral-action)" — is named in Paper 34 §X (the open-projections-on-review list, in the synthesis section reading "Are there projections in the framework that are still not yet in this list?"). The entry has been pending since pre-Sprint 1.

Does a standalone "Schwinger proper-time" §III entry resolve open-question (a)? **No, only partially.** Open-question (a) is specifically about *fractional* $s$ (i.e., $\zeta_A(s)$ at $s \notin \mathbb{Z}$, where the master Mellin engine output is not constrained by the T9 theorem to $\pi^{2k}\cdot\mathbb{Q}$ and can in principle access richer transcendental classes). Sprint Q (Paper 28) explored fractional-$s$ behaviour of the Dirac spectral zeta and found the transcendental class is discontinuous at integer $s$ vs fractional $s$ (memory file `fractional_order_spectral_zeta.md`). This is a genuine open question about the *output ring* of the Mellin engine evaluated off the integer-$s$ slice, not about whether the engine itself deserves its own §III slot.

Resolving open-question (a) requires either: (i) a clean statement of which transcendental ring $\mathrm{Mel}(s)$ lands in at non-integer $s$, with empirical anchors; or (ii) a structural argument that fractional-$s$ evaluation produces a *new sub-mechanism* M4 distinct from M1/M2/M3 (which would expand the master Mellin engine partition). Neither task is accomplished by promoting Schwinger proper-time to a §III entry — that promotion would simply rename existing content.

The right disposition for open-question (a) is to leave it on the §X list as a follow-up sprint candidate (rate this as a $\sim$2-week diagnostic to test whether Sprint Q's discontinuity finding generalizes), separate from the present scoping question. The fractional-$s$ question is genuinely open; the heat-kernel-as-its-own-projection question is structurally resolved by §III.6 + the Paper 32 §VIII cross-reference.

---

## 5. Verdict

**Absorbed into §III.6, with a cross-reference paragraph to the Paper 32 §VIII / Paper 18 §III.7 master Mellin engine framing.**

The structural test in §1 returns negative: Schwinger proper-time integration introduces no physical variables, dimensions, or transcendental content beyond what §III.6 already encodes. The proper-time integral $\int_0^\infty dt\, f(t)\, \mathrm{Tr}\, e^{-tD^2}$ and the Connes–Chamseddine spectral action $\mathrm{Tr}\, f(D/\Lambda)$ are the same operation under two names; both inject $\Lambda$ (UV cutoff) as the energy variable, both produce the $\sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2\cdot\mathbb{Q}$ ring on Seeley–DeWitt coefficients, both live in the M2 sub-mechanism of the master Mellin engine.

Promoting Schwinger proper-time to a standalone §III entry would:
- Double-count the M2 mechanism under two §III labels.
- Confuse the dictionary by treating the Mellin *evaluator* (a meta-mechanism organizing all spectral-side projections) as a *peer* of the projections it evaluates.
- Not resolve the genuinely-open §X question (a) about fractional-$s$ behaviour, which is about the *output ring* of the Mellin engine off-integer-$s$, not about whether the engine deserves its own slot.

**Recommended action (deferred to PI decision; this memo does not modify Paper 34):**

Add a short cross-reference paragraph at the end of §III.6 (~3–5 sentences) noting that:

1. The heat-kernel evaluation of §III.6 is the M2 sub-mechanism ($k=2$) of the master Mellin engine of Paper 32 §VIII / Paper 18 §III.7.
2. The same Mellin technology, applied at $k=0$ to $\mathrm{Tr}\, e^{-tD^2}$, recovers the Hopf-base measure (M1, §III.2 / Paper 25); at $k=1$ to $\mathrm{Tr}(D\, e^{-tD^2})$, recovers the vertex-parity Hurwitz content (M3, §III.7 spinor lift composed with vertex-parity restriction, Paper 28 §IV).
3. The Mellin transform itself is the evaluator/meta-mechanism, not a §III peer projection; consequently §III.6 covers the operational content of "Schwinger proper-time integration" without a separate entry.
4. The open question (Paper 34 §X (a)) about fractional-$s$ behaviour remains structurally separate and is the natural next sprint target on this mechanism.

This is a $\sim$10-line edit to §III.6, not a new §III entry, and keeps the dictionary's nineteen-projection structure intact while making the master-Mellin-engine reading explicit at the §III level rather than only at the Paper 32 / Paper 18 level.

**Counterargument considered and rejected:** one could argue that "Schwinger proper-time" deserves its own slot because it is a *continuous-parameter integration over $t \in (0, \infty)$*, while spectral action sums over modes with a test function $f$. This is a textbook distinction (Schwinger's proper-time formalism vs Wilsonian mode-cutoff), but on the spectral side they reduce to the same computation: the test function $f$ in $\mathrm{Tr}\, f(D/\Lambda)$ can always be rewritten as an inverse Laplace transform giving the proper-time integral, and conversely the proper-time integral against a cutoff function recovers the spectral action. The distinction is conventional, not structural, and Paper 34's dictionary tracks structural projections. (Paper 18's Mellin section already presents the two forms interchangeably: eq. `mellin_bridge` is the proper-time form; the spectral-action discussion uses the test-function form; they are the same.)

If, in a future sprint, a *non-Connes–Chamseddine* application of Schwinger proper-time emerges (e.g., a real-time Schwinger–Keldysh closed-time-path application on a Lorentzian wedge algebra, which Paper 34 §III.16-class extensions hint at), that would be the place to revisit the question. Until then, §III.6 + cross-reference is the cleanest disposition.

---

## Appendix: anchor citations consulted

- Paper 34, §III.6 (`sec:proj_spectral_action`, lines 319–335): existing spectral-action entry, target/variables/transcendental signature.
- Paper 34, §X open-question list (lines ~5995–6013): pending question (a) on Mellin / heat-kernel at fractional $s$.
- Paper 34, §III.15 / §III.28 (apparatus identity, Sprint TD Track 5 context, lines ~2470–2540): spectral-vs-state-side divide and the M1∪M2∪M3 disjointness from state-side observables.
- Paper 32 §VIII case-exhaustion theorem (cross-referenced throughout Paper 34 as `loutey_paper32` §VIII): every π in any GeoVac observable engages M1, M2, or M3.
- Paper 18 §III.7 master-mechanism reading (`sec:mellin`, lines 693–940): Mellin transform as taxonomic engine, master-mechanism reading, mechanism-as-domain sharpening.
- Paper 18, eq. `mellin_bridge` (line ~707): Mellin bridge $\zeta_A(s) = \Gamma(s)^{-1} \int_0^\infty t^{s-1} \mathrm{Tr}(e^{-tA}) dt$.
- Paper 18, eq. `master_mellin` (line ~853): π-source = $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})](s)$ with $k \in \{0,1,2\}$.
- Paper 18, eq. `dirac_modular_residual` (line ~885): MR-B closed-form modular residual confirming M2 ring.
- Memory `fractional_order_spectral_zeta.md`: Sprint Q finding that fractional-$s$ behaviour is discontinuous at integer $s$ — relevant for §X open question (a), distinct from this scoping question.
