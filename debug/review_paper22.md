# Confidence Review: Paper 22 — Potential-Independent Angular Sparsity for Qubit Hamiltonians in Spherical Fermion Systems

## Calibration check

Not a calibration run.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | Two-body matrix element factorizes as Σ_k c^k(l1m1;l3m3) c^k(l2m2;l4m4) R^k (Theorem 1) | Eq. eq:factorization, §II | A | EXTERNAL (Condon–Shortley 1935; Messiah; Slater 1960 — standard result) | Standard Slater–Condon angular-momentum algebra, derivation in any AMO graduate text |
| 2 | Zero pattern of the two-body matrix is set by Gaunt selection rules and is independent of V(r), R_{nl}, n (Theorem 2) | §II, body | A | EXTERNAL (immediate consequence of Theorem 1's factorization) | Symmetry-direction trivially follows from Theorem 1; converse argument is sound |
| 3 | Angular ERI density at l_max=1,2,3,4,5 is 7.8125%, 2.76%, 1.44%, 0.90%, 0.62% (Theorem 3 table) | §III table; abstract; conclusion | **C overstated → E in framing** | GEOVAC-ONLY (uses pair-diagonal-m, which is NOT the full Coulomb selection rule despite the paper repeatedly framing it as such) | Reproduced 1.4404% with sympy under m_a=m_c AND m_b=m_d (pair-diagonal). Independent recomputation under physically correct full-Gaunt rule m_a+m_b=m_c+m_d gives **6.0608% at l_max=3**, **8.52% at l_max=2**. The paper's own §V Table II labels the full-Gaunt convention "physically correct Coulomb selection rule." |
| 4 | 97.24% of two-body matrix elements vanish by Wigner–Eckart at l_max=2 | abstract, §IV | E (same bug as Claim 3) | GEOVAC-ONLY | Under full-Gaunt the figure is 91.48%, not 97.24% |
| 5 | Bit-identical angular zero patterns across Coulomb / HO / Woods–Saxon / square well / Yukawa at n_max=2, l_max=2 | Table I, §IV | A | EXTERNAL (Gaunt coefficients are functions of (l,m,k) only) | The bit-identicality across potentials is correct under either convention. The pattern IS potential-independent; the magnitude reported is from a stricter-than-physical selection. |
| 6 | D(l_max, n_max) = D(l_max) ∀ n_max ≥ 1 (Eq. eq:nmax_independence) | §III | A | EXTERNAL (direct corollary of Theorem 2) | Argument valid; only depends on potential-independence of angular zeros |
| 7 | Composed architecture O(Q^2.5) Pauli scaling is potential-independent (Corollary 1) | §V.B | B (proof-sketch) | MIXED: angular structure EXTERNAL; the 2.5 exponent is GEOVAC-ONLY empirical (Paper 14) | The argument that block-diagonal + Theorem 2 ⇒ exponent invariance is sound; the empirical 2.5 itself is internally measured. |
| 8 | $S^3$ Fock projection and π-free graph are Coulomb-specific (Theorem-style claim in §VI) | §VI | A | EXTERNAL (Fock 1935; SO(4) accidental degeneracy is unique to 1/r) | Standard, e.g., Bertrand's theorem |
| 9 | Spinor density d_spinor ≤ d_scalar at every l_max; ratio → 1 as l_max → ∞ (§paper22_spinor) | §V (spinor block) | A | EXTERNAL (jj-coupled reduction in Dyall §9.3 / Grant §7.5 + j-triangle adds zeros) | Standard Dyall/Grant; verified by the spinor density driver. |
| 10 | Breit rank-2 angular density is 6–8× denser than Coulomb but monotonically decreasing in l_max (§paper22_breit) | §VI (Breit block) | B (proof deferred) | GEOVAC-ONLY in current form (paper says "proof... deferred to a formal write-up") | The rank-1+rank-2 enumeration in `debug/br_a_breit_angular.py` is reproducible from the paper, but the structural theorem statement at the end is not proven in-paper. |

### Numbers I recomputed

| claim | paper's figure | independent reference | my recomputed value | survives? |
|---|---|---|---|---|
| D(l_max=3) under "Coulomb selection rule" | 1.4404% | Standard Slater–Condon (Messiah Vol II Ch XIII) — full-Gaunt convention m1+m2=m3+m4 | **6.0608%** (sympy, exact) | NO under physical Coulomb convention; 1.4404% only recovers under the stricter pair-diagonal-m convention (m1=m3 AND m2=m4) |
| D(l_max=2) under "Coulomb selection rule" | 2.7587% | Standard Slater–Condon | **8.5200%** under full-Gaunt | NO (same reason as above) |
| D(l_max=2) under pair-diagonal-m convention | 2.7587% | Paper's own §V Table II | 2.7587% (sympy) | YES (under the documented stricter convention) |
| Angular sparsity invariance across 5 potentials at fixed (l_max, n_max) | "bit-identical" | Theorem 2 derivation | Bit-identical at l_max=2 across HO/WS/Coulomb/SQW/Yukawa (production driver), under either convention | YES (zero pattern under either convention is V-independent) |
| Spinor density d_sp = 1.2329% at l_max=3 (pair-diag), 5.1674% (full-Gaunt) | Table II | tier2_t0_spinor_density.json | Matches the JSON dump | YES |

### Circularity map

The GEOVAC-ONLY chains:

1. **The "1.44% at l_max=3" and the other headline densities** rest only on `geovac/nuclear/potential_sparsity.py::angular_zero_count`. The function uses `ck_coefficient` with `q = mc - ma`, which produces a 3j with bottom-row sum `-ma + (mc-ma) + mc = 2(mc-ma)`, vanishing whenever m_a ≠ m_c. The function therefore implicitly enforces m_a = m_c AND m_b = m_d ("pair-diagonal-m"), not the full Coulomb rule m_a + m_b = m_c + m_d. The paper text in Eq. eq:gaunt writes "q = m - m'" with the correct sign, but the headline numbers came from the buggy convention. The paper's own §V Table II caption acknowledges this (calling pair-diagonal-m "stricter than the full Gaunt rule" and the full-Gaunt convention "physically correct Coulomb selection rule") — yet the abstract, §III, §IV, and §VI all continue to present 1.44% / 2.76% / 7.81% as "the angular ERI density of the Coulomb interaction."

2. **The "97.24% of all possible two-body matrix elements are guaranteed to vanish by the Wigner–Eckart theorem"** (abstract) inherits the same problem: under the physically correct Coulomb selection rule, the figure is 91.48%, not 97.24%.

3. **Theorem 2 itself (potential-independence)** is structurally correct under either convention — the GEOVAC-ONLY component is only the magnitude of the reported densities, not the existence of the theorem. This is important because Papers 14, 23, and 31 cite Paper 22's *theorem* (which is solid) but quote the *headline numerical density* (which is generated by a stricter selection rule than the one labelled "Coulomb" in the paper text).

4. **Corollary 1 (potential-independent O(Q^2.5) scaling)** rests on Paper 14's empirical 2.5 exponent + Theorem 2's structural claim. The first half is GEOVAC-empirical; the second half is EXTERNAL when Theorem 2 is interpreted with full-Gaunt selection.

Cross-corpus check: Paper 23 cites Paper 22's theorem as justification for sub-linear Pauli scaling on nuclear shell models. Paper 31 cites the "1.44% at l_max=3" specifically as the canonical content of the universal sector (line 35 of paper 31). The convention issue **propagates forward** to both downstream papers.

### Overstatement findings

| Exact phrase (paper) | Issue | Suggested honest replacement |
|---|---|---|
| Abstract: "Exact enumeration of the Gaunt coefficients yields the angular ERI density D(l_max): 7.81% at l_max=1, 2.76% at l_max=2, 1.44% at l_max=3, 0.90% at l_max=4, and 0.62% at l_max=5." | The cited numbers are obtained under a stricter-than-Coulomb selection rule (pair-diagonal-m, m_a=m_c AND m_b=m_d). Under the physically correct Coulomb m-conservation, the densities are 14.84%, 8.52%, 6.06%, 4.83%, 3.99% (Paper 22's own §V Table II "full-Gaunt" column). | "Exact enumeration under the pair-diagonal-m convention (m_a=m_c and m_b=m_d per electron, corresponding to a q=0 multipole selection) yields D_pd(l_max): 7.81% / 2.76% / 1.44% / 0.90% / 0.62%. Under the physically correct Coulomb selection rule m_a+m_b=m_c+m_d (allowing q≠0 multipole transfers), the densities are 14.84% / 8.52% / 6.06% / 4.83% / 3.99% (Table II)." Or, more aggressively: lead the abstract with full-Gaunt values; demote pair-diagonal-m to a lower bound. |
| Abstract: "At l_max=2, 97.24% of all possible two-body matrix elements are guaranteed to vanish by the Wigner–Eckart theorem before any property of the potential is invoked." | Specific magnitude tied to the stricter convention. The Wigner–Eckart-guaranteed zero count under physical Coulomb is 91.48%, not 97.24%. | "At l_max=2, at least 91.48% of all possible two-body matrix elements are guaranteed to vanish by the Wigner–Eckart theorem under m-conservation alone (97.24% under the stricter pair-diagonal-m q=0 sub-selection)." |
| §III Definition / Theorem 3: "(ii) at least one value of k satisfies the triangle and parity conditions for both Gaunt factors with both factors nonzero." | The exact-enumeration script implements a stricter condition (q1=q2=0). The definition text and the implementation diverge. | Either re-derive Table I under the definition as written (giving the FG numbers) or revise the definition to state "and m_a = m_c, m_b = m_d (q=0 sub-block)." |
| §IV Numerical Verification, Table I: "Angular zeros (%) 97.24" across all five potentials. | The driver in `geovac/nuclear/potential_sparsity.py::compute_two_body_eri` uses the buggy `ck_coefficient` with q-sign reversed; so the bit-identical 97.24% across potentials reflects the pair-diagonal-m sub-block, not full Coulomb. The "bit-identical" claim is true (under either convention the V-independence holds), but the magnitude is from the stricter convention. | State the convention explicitly in the table caption; report both numbers if both are of interest. |
| Conclusion: "The explicit angular ERI densities obtained by exact enumeration are 7.81% at l_max=1, 2.76% at l_max=2, 1.44% at l_max=3, 0.90% at l_max=4, and 0.62% at l_max=5." | Same as abstract. | Same as abstract. |
| Conclusion: "Numerical verification across five potentials (Coulomb, harmonic oscillator, Woods–Saxon, square well, Yukawa) at l_max=2 confirms the theorem to machine precision: the angular zero patterns are bit-identical." | The bit-identicality claim is correct, but it confirms potential-independence of the *zero pattern under the stricter pair-diagonal-m convention used in the driver*. It does not by itself confirm the physically correct Coulomb selection-rule density. | "Numerical verification across five potentials at l_max=2 confirms that the zero pattern is bit-identical across potentials, consistent with Theorem 2's structural prediction. The specific magnitude reported is for the pair-diagonal-m (q=0) sub-block; the full-Gaunt density follows by direct enumeration (Table II of §V)." |
| §I last paragraph: "the contribution of this paper is the explicit quantitative bound at each l_max" | The explicit quantitative bound, as stated in the paper, is wrong for Coulomb. The structural framing/theorem is the actual contribution. | "the contribution of this paper is to make explicit and quantitatively bound the angular sparsity that follows from the standard radial–angular factorization, at each l_max." |
| §sec:paper22_breit caption claim: "the angular sparsity theorem extends to rank-2" | "Verdict" stated without theorem proof in-paper; structural extension framed as established but proof is deferred to `debug/br_a_breit_memo.md`. Acceptable but should be flagged. | "Numerical evidence (Table) is consistent with the rank-2 angular sparsity bound following the same structural pattern; a formal Wigner–Eckart-based proof for arbitrary tensor rank is deferred to follow-on work." |

### Broken internal cross-references

| `\ref{}` | Where used | Target exists in paper? |
|---|---|---|
| `\ref{sec:theorem}` | §V (paper22_spinor body), §VI (paper22_breit body) | NO. §II is labeled `sec:factorization`, not `sec:theorem`. |
| `\ref{sec:universal_partition}` | §V paper22_spinor | NO. §VI is labeled `sec:universal`. |
| `\ref{sec:taxonomy}` | §VI paper22_breit (Paper 18 reference) | NO label in this paper; meant as a Paper 18 §III cross-reference which won't resolve in single-paper compile. |
| `\S\ref{sec:spinor_composed}` / `\ref{tab:spinor_resource}` | §V paper22_spinor | Cross-paper references to Paper 14; will appear as `??` in a standalone Paper 22 build. |

LOW severity but cosmetic warts that would produce LaTeX `??` placeholders.

## Pass B — Citation and novelty

### Citation table

| `\cite` key | Claimed as | Verdict | Notes |
|---|---|---|---|
| `Peruzzo2014` | VQE foundational paper | CITE-OK | Nat. Commun. 5, 4213 (2014); canonical reference. |
| `Lee2021` | Tensor hypercontraction for FT chemistry | CITE-OK | PRX Quantum 2, 030305 (2021); canonical THC reference. |
| `Whitten1973` | RI/DF Coulomb integral approximation | CITE-OK | J. Chem. Phys. 58, 4496 (1973); canonical RI paper. |
| `Dunlap2000` | Robust and variational fitting | CITE-OK | Mol. Phys. 98, 1639 (2000). Authors and topic correct. |
| `Suhonen2007` | Talmi–Moshinsky / nuclear NN in HO | CITE-OK | "From Nucleons to Nucleus" — standard graduate text. |
| `CondonShortley1935` | Theory of Atomic Spectra | CITE-OK | Standard textbook; original Cambridge edition. |
| `Slater1960` | Quantum Theory of Atomic Structure, Vol I | CITE-OK | Standard McGraw-Hill text. |
| `Messiah1962` | Quantum Mechanics, Vol II | CITE-OK | Standard North-Holland text. |
| `Reimann2002` | Electronic structure of quantum dots | CITE-OK | Rev. Mod. Phys. 74, 1283 (2002); canonical QD review. |
| `Dyall` | jj-coupled Coulomb reduction §9.3 | CITE-OK (no full bibitem, but pointing to Dyall & Fægri *Introduction to Relativistic Quantum Chemistry* §9.3); should have a full bibitem | Used as `~\cite{Dyall}` without a `\bibitem{Dyall}` entry in the bibliography — this is a missing-bibitem warning. |
| `Grant` | jj-coupled relativistic atomic structure §7.5 | CITE-OK (Grant, *Relativistic Quantum Theory of Atoms and Molecules*, Springer 2007 §7.5) but **missing bibitem in this paper's `\thebibliography`** | Same as above. |
| `GeoVac_Paper7` | Dimensionless Vacuum / S³ proof | CITE-OK (internal) | Internal reference. |
| `GeoVac_Paper14` | Qubit encoding / O(Q^2.5) scaling | CITE-OK (internal) | Internal reference. |
| `GeoVac_Paper17` | Composed natural geometries | CITE-OK (internal) | Internal reference. |
| `GeoVac_Paper18` | Exchange constants taxonomy | CITE-OK (internal) | Listed in bibliography but **not actually cited anywhere in the body** — orphan bibitem. |
| `GeoVac_Paper23` | Fock rigidity / Coulomb-specific | CITE-OK (internal) | Internal reference. |

### Problems found

**CITE-CANT-FIND (load-bearing):**

- `\cite{Dyall}` in §V (paper22_spinor) — used to support the jj-coupled Coulomb reduction, but **the bibliography contains no `\bibitem{Dyall}`**. This is a load-bearing citation (it underwrites the spinor extension's claim of "structural identity" between scalar and spinor selection rules). MEDIUM severity (citation key dangles to `[?]` in compile).
- `\cite{Grant}` in §V — same issue, no `\bibitem{Grant}` in `\thebibliography`. MEDIUM severity.

**Orphan bibitem (LOW):**

- `\bibitem{GeoVac_Paper18}` is defined in the bibliography but not cited anywhere in the body of Paper 22. LaTeX standard behavior is "Bibliography item not cited" warning. (The reference to "Paper 18" in §sec:paper22_breit is `Paper~18 \S\ref{sec:taxonomy}` without a `\cite`.)

**Novelty claims:**

The paper does NOT make a "first in the literature" or strong priority claim. Specifically: §I states "The contribution of this paper is the explicit quantitative bound at each l_max...." This is appropriately humble. §I also acknowledges Condon–Shortley and Slater for the underlying factorization. No novelty claim to soften.

The paper does claim the spinor extension (§V) is "structurally identical" to the scalar case via the jj-reduction. This is a standard result (Dyall §9.3 / Grant §7.5) — but the missing bibitems prevent a reader from looking it up.

The Breit rank-2 enumeration (§VI) is presented as a verification rather than a novelty claim, with the proof "deferred to a formal write-up." Acceptable framing.

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---|---|---|---|
| 1 | Headline densities (1.44% etc.) computed under pair-diagonal-m, framed as Coulomb selection rule throughout abstract / §III / §IV / §VII conclusion | A | E (math/framing error) | **HIGH** |
| 2 | "97.24% guaranteed to vanish by Wigner–Eckart" in abstract is the pair-diagonal-m figure; physical Coulomb figure is 91.48% | A | E | **HIGH** |
| 3 | Convention mismatch propagates forward to Papers 14, 23, 31 which cite "1.44% at l_max=3" as load-bearing | A | E | **HIGH** (cross-corpus impact) |
| 4 | `\bibitem{Dyall}` and `\bibitem{Grant}` missing — load-bearing citations for the spinor §V section dangle | B | CITE-CANT-FIND | **MEDIUM** |
| 5 | Theorem-3 "Definition" text in §III does NOT include the pair-diagonal-m restriction, but the driver does — internal inconsistency between paper text and driver | A | E | **MEDIUM** |
| 6 | Conclusion claim "Numerical verification across five potentials confirms the theorem to machine precision" overstates: the test confirms zero-pattern V-independence (real), but is framed as if it independently validates the magnitude of the density (which depends on the stricter convention) | A | C overstated | **MEDIUM** |
| 7 | Stale cross-refs `\ref{sec:theorem}`, `\ref{sec:universal_partition}`, `\ref{sec:taxonomy}`, `\ref{sec:spinor_composed}`, `\ref{tab:spinor_resource}` — produce LaTeX `??` placeholders | A | E (cosmetic logic error in a real compile) | LOW |
| 8 | Bibitem `GeoVac_Paper18` declared but never cited (orphan) | B | LOW | LOW |
| 9 | §VI Breit "verdict: angular sparsity theorem extends to rank-2" is asserted with proof deferred to a memo; OK with framing flag | A | B internally consistent only | LOW |
| 10 | Corollary 1 proof sketch ("$O(Q^{2.5})$ for any spherical potential") relies on Paper 14's empirical exponent, not on an external benchmark | A | B internally consistent only | LOW |
| 11 | Citations Peruzzo2014, Lee2021, Whitten1973, Dunlap2000, Suhonen2007, CondonShortley1935, Slater1960, Messiah1962, Reimann2002 all check out | B | CITE-OK | (positive) |

A–E counts (Pass A, claims table): A = 5, B = 3, C = 0 (one C in overstatement only), D = 0, E = 2 (one math-magnitude, one definition-driver mismatch).
Overstatement findings (additional C-class): 5.

CITE-* counts (Pass B): CITE-OK = 9 published + 5 internal = 14; CITE-CANT-FIND = 2 (Dyall, Grant); CITE-MISATTRIBUTED = 0; CITE-DOESNT-SUPPORT = 0; CITE-WRONG-METADATA = 0.

HIGH = 3 (cf. items 1, 2, 3 above). MEDIUM = 3 (items 4, 5, 6). LOW = 5 (items 7, 8, 9, 10, and the additional table-caption nitpicks).

## Broadcast readiness: **RED**

This paper has a real, load-bearing numerical claim — "1.44% angular ERI density at l_max=3" — which is consistently framed throughout the abstract, §III table, §IV verification, and §VII conclusion as the Coulomb selection rule's quantitative bound. The paper's own §V Table II (added later) explicitly acknowledges that this number comes from a **stricter-than-Coulomb pair-diagonal-m convention** (m_a = m_c AND m_b = m_d per pair) and that the physically correct Coulomb selection rule m_a + m_b = m_c + m_d gives a substantially denser ERI density (6.06% at l_max=3 — a factor of 4.2 higher). The production driver `geovac/nuclear/potential_sparsity.py::angular_zero_count` reproduces the headline numbers because its helper `ck_coefficient` uses `q = m_c - m_a` instead of `q = m_a - m_c`, which automatically zeros every 3j with m_a ≠ m_c. The paper's text in Eq. eq:gaunt actually states the correct sign convention `q = m - m'`, but the headline numbers were produced under the buggy driver. The whole structural theorem (Theorem 2, potential-independence) is correct under either convention — the impact is on the *magnitude* of the reported density, not on the existence of the result. However, the magnitude IS the headline claim. **This needs to be resolved before broadcast**: either (a) recompute the headline densities under the physical Coulomb selection rule and re-state the abstract, §III, §IV, §VII numbers, or (b) re-frame the entire paper as quantifying the q=0 sub-block ("pair-diagonal Coulomb sparsity") and explicitly note that the full Coulomb density is higher by ~4×. Papers 14, 23, and 31 cite the "1.44%" number directly, so a fix here cascades. The missing `\bibitem{Dyall}` and `\bibitem{Grant}` are smaller but still need to be added because they support the spinor-extension argument.

## What I could NOT verify (hand to a human expert)

1. **Whether the explicit angular-density tabulation at each l_max is genuinely new.** A literature search would need to canvas (a) Slater's *Quantum Theory of Atomic Structure* (1960), (b) Pekeris-type and Hylleraas tabulations from 1950s–60s atomic physics, (c) Talmi–Moshinsky / Brody–Moshinsky tabulations from nuclear physics, (d) chemistry "Density Fitting" sparsity bounds. The closest I can confirm is that the *radial–angular factorization* is in Condon–Shortley, but the specific potential-independence-of-density tabulation framed as a quantum-computing resource bound is unlikely to be in the cited textbooks. Pass B's honest ceiling on novelty applies. The paper does not actually claim "first" — it claims "explicit quantitative bound at each l_max." Whether that bound is already tabulated somewhere (e.g., for nuclear shell models) is a domain-expert question.

2. **Whether the "pair-diagonal-m" sub-block is a meaningful physical operator** (e.g., the q=0 component of the Coulomb multipole expansion, restricted to m-diagonal channels). If so, the paper's headline numbers might be salvageable as a bound on that specific sub-block, accompanied by the full-Gaunt numbers. A quantum-chemistry expert is the right reviewer for this distinction.

3. **Whether the Wigner-Eckart-based formal proof for arbitrary tensor rank (Breit §VI verdict) holds** as a theorem. The paper defers the proof to `debug/br_a_breit_memo.md`; only a relativistic-AMO domain expert can certify that the proof generalizes correctly.
