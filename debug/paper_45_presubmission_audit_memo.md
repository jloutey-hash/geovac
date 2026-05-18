# Paper 45 Pre-Submission Hardening Audit Memo

**Date:** 2026-05-18
**Paper:** `papers/standalone/paper_45_lorentzian_propinquity.tex`
**Title:** Lorentzian propinquity convergence on truncated SU(2) × U(1)_T Krein spectral triples
**Status entering audit:** drafted, three-pass clean compile, zero substantive LaTeX warnings
**Verdict:** **REQUIRES_FIXES** — five citation errors (two factually wrong attribution, one wrong title, one date inconsistency, one bibliography inconsistency with Paper 38). Once fixed, **ARXIV_READY**. No scope/framing changes required. No proof gaps. Tests pass.

---

## Track 1 — Reviewer-style critique

### 1.1 Scope claims (load-bearing K⁺-weak-form framing)

**PASS.** The K⁺-weak-form framing is consistently load-bearing across the paper:

- **Abstract** (lines 129–185): "restricted to the Krein-positive state space" stated up front; the dictionary "K⁺-restricted weak-form Lorentzian propinquity" appears in line 152; honest scope paragraph (lines 178–184) names strong-form as open and identifies Mondino–Sämann distinction.
- **§1.2** (line 232): "What this paper does" frames the construction as five-lemma argument transferred from Paper 38.
- **§1.3** (lines 269–295): Honest scope and named gaps — four explicit gaps G1–G4.
- **§7.3 Open questions** (lines 1296–1347): five named open questions Q1–Q5, with Q1 explicitly the strong-form gap.
- **§7.4 Structural-skeleton scope** (lines 1349–1358): defensible scope statement.

No softening of scope detected. The "first Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature" claim (line 174) is qualified by "weak-form" and "K⁺-restricted" everywhere it appears — defensible.

### 1.2 Proof rigor (per-lemma audit)

**§4.1 Lemma L2 (joint cb-norm):** PASS. Proof at lines 681–718 cites Bożejko–Fendler 1991 + Pisier 2001 Ch. 8 Thm 8.10, applies factor-wise via Paper 38 L2(c) for SU(2) and direct geometric-series for U(1), and combines via Takesaki IV.4.14 (commutative-factor tensor product). Memo §3 verifies this verbatim.

**§4.2 Lemma L3 (joint Lichnerowicz):** PASS. The headline structural identity (`[D_L, a_s ⊗ a_t] = i[D_GV, a_s] ⊗ a_t` bit-exact) is proved in lines 770–803 via two-term decomposition: cross-term vanishes because (i) [γ⁰, a_s] = 0 (chirality-block-diagonal), (ii) [∂_t, a_t] = 0 (both momentum-diagonal). Tensor norm factorization gives the bound. Sub-sprint A memo §3 verifies. The bound C_3^joint ≤ 1 asymptotically tight is inherited verbatim from Paper 38 Lemma L3.

**§4.3 Lemma L4 (joint Berezin):** PASS. Lines 909–984 prove five properties (a) positivity via UCP compression of non-negative convolution, (b) contractivity via Lemma L2's cb-norm 2/(n_max+1), (c) approximate identity with γ^joint = max(γ_SU(2), γ_U(1)), (d) L3 compatibility, (e) Krein-positivity preservation. The tensor-product factorization Remark 4.6 (lines 971–984) makes the SU(2) inheritance from Paper 38 L4 explicit. Memo C verifies.

**§5 L5 (propinquity assembly):** PASS. Lines 1018–1147 follow Latrémolière 2017 §4 weak-form construction verbatim on the K⁺ substrate. Proposition 5.2 (reach + height bounds) proves the four constituents; main theorem 5.3 combines via max(γ, γ, C_3·γ, 0) = C_3·γ. Memo D verifies.

**Recommended polish (NOT a fix):** Proposition 5.2 cites "Latrémolière 2017 Def. 3.5" for the Lipschitz-distortion height definition — fine, but a parenthetical "(Definition 3.5 of arXiv:1811.10843)" would be more reader-friendly. Acceptable to defer.

### 1.3 Definitions

**PASS.** Definition 2.3 (lines 553–568) defines the K⁺-restricted weak-form propinquity as Λ^L(T₁^L, T₂^L) := Λ(T₁⁺, T₂⁺) with Λ the standard Latrémolière metric-spectral-triple propinquity (citing latremoliere_metric_st_2017). The construction inherits all standard apparatus (Hilbert space, ‖[D, ·]‖_op Lipschitz seminorm). Definitions of K⁺ subspace (2.2), tunneling pair (5.1), and joint Berezin (4.1) are all precise enough for a careful reader to reconstruct.

### 1.4 Citation spot-checks (FIVE ERRORS FOUND)

| Cite key | Claim in Paper 45 | Verified | Verdict |
|----------|-------------------|----------|---------|
| `bizi_brouder_besnard2018` | J.Math.Phys. 59 (2018), 062303; arXiv:1611.07062 | YES via arXiv abstract page | OK |
| `vandungen2016` | Math.Phys.Anal.Geom. 19 (2016) Art.4 22pp; arXiv:1505.01939; "Krein spectral triples and the fermionic action" by Koen van den Dungen | YES verified verbatim | OK |
| `mondino_samann2025` (arXiv:2504.10380) | "Synthetic Lorentzian Gromov-Hausdorff convergence and pre-compactness" by Mondino, Sämann | YES verified | OK |
| **`mondino_samann2025b` (arXiv:2510.13069)** | **Attributed to Mondino, Sämann** | **WRONG — actual authors: Mauricio Che, Raquel Perales, Christina Sormani** | **REQUIRES_FIX** |
| **`vandungen2026_su11` (arXiv:2601.22171)** | **"Pseudo-Riemannian Spectral Triples for SU(1,1)" by K. van den Dungen** | **WRONG — actual author: Jort de Groot** | **REQUIRES_FIX** |
| `nieuviarts2025a` (arXiv:2502.18105) | (v3, May 2025) | YES — v3 dated May 6, 2025 | OK |
| **`nieuviarts2025b_proceedings` (arXiv:2512.15450)** | "Emergence of pseudo-Riemannian structures from twisted spectral triples within the almost-commutative framework (proceedings synthesis)" (v2, May 2026) | **TITLE WRONG — actual title: "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry"** by G. Nieuviarts, December 17, 2025 (revised May 11, 2026) | **REQUIRES_FIX** |
| `farsi_latremoliere2025` (arXiv:2504.11715) | "Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics" by Farsi & Latrémolière, 2025 | YES verified | OK |
| `bozejko_fendler1991` | "Boll. Un. Mat. Ital. A (6) 3 (1984), 297–302; extended version 1991" | **INCONSISTENT — Paper 38 (and the standard literature) cites this as "Arch. Math. 57 (1991), 290–298" with title "Herz–Schur multipliers and uniformly bounded representations of discrete groups"** | **REQUIRES_FIX** (align with Paper 38 to avoid confusion) |
| `latremoliere_metric_st_2017` | "preprint arXiv:1811.10843, 2017" | **MINOR DATE INCONSISTENCY — arXiv preprint was Nov 2018, not 2017. Published Adv. Math. 2023 is correct.** | **REQUIRES_FIX** (one-character year change) |

The five-error pattern is the standard "author/title drift" that happens when bibliography entries are populated from secondary sources. None of the errors invalidate the proof structure or the paper's claims. All five are mechanical fixes.

### 1.5 Related-work paragraph (§1.4)

**PASS with minor recommendation.** The Mondino–Sämann distinction is crisply stated at lines 338–346: "synthetic Lorentzian Gromov–Hausdorff program... lives on a categorically different mathematical object: pre-length spaces with causal-diamond geometry, not operator-algebraic spectral triples." Repeated more rigorously in §7.2 Remark 7.1 (lines 1284–1294). The Nieuviarts odd-dim obstruction is correctly stated at lines 332–336 referencing Paper 42 §3.2.

The Strohmaier 2006 / Franco–Eckstein 2014 / van den Dungen 2016 / BBB 2018 lineage is correctly enumerated at lines 321–336. The recent Riemannian-side continuity reference (Farsi–Latrémolière 2025) is cited at lines 302–307.

**Recommended polish:** the related-work paragraph could be tightened by adding one sentence noting that no published Lorentzian propinquity exists in the May 2026 literature (this is implicit in the negation but would be explicit and useful for a reviewer's quick scan). Acceptable to defer.

### 1.6 Numerical panel §6

**PASS.** Table 1 (lines 1174–1198) matches `debug/data/l3b_2_sub_sprint_D.json` verbatim:

| Cell | Paper | JSON | Match |
|------|-------|------|-------|
| (2,1) γ_SU(2) | 2.0746 | 2.0745510936998897 | ✓ |
| (3,1) γ_SU(2) | 1.6101 | 1.6100599680657361 | ✓ |
| (4,1) γ_SU(2) | 1.3223 | 1.3223327942828407 | ✓ |
| (2,3) γ_U(1) | 0.7220 | 0.7219831099305748 | ✓ |
| (3,5) γ_U(1) | 0.4956 | 0.49562415926901987 | ✓ |
| (4,7) γ_U(1) | 0.3841 | 0.38406187302086475 | ✓ |
| All N_t=1 Riem residuals | 0.0 bit-exact | match_bit_exact=true, gamma_residual=0.0 | ✓ |
| Ratio Λ(4,7)/Λ(2,3) | 0.6374 | 0.6374067133359951 | ✓ |

The Riemannian-limit bit-exact recovery at N_t = 1 (load-bearing falsifier) is verified in the JSON `riemannian_limit_check` block — all three cells have `match_bit_exact: true` and `gamma_residual: 0.0`.

**Minor observation (NOT a fix):** the JSON additionally reports `reach_B_panel` values (0.151, 0.212, 0.315) that are tighter than the bound `reach_B_bound` (1.20, 1.14, 1.02), and these are tighter than `gamma_joint_su2`. Paper 45 reports the theoretical bound `Λ^L ≤ C_3·γ`, which is dominated by γ_SU(2) at every cell — correct presentation for the theorem statement. The tighter panel values support the bound but aren't claimed as tight.

---

## Track 2 — Concurrent-work re-check

Fresh search performed 2026-05-18. Five web searches across distinct angles:

1. **`"Lorentzian propinquity" Krein spectral triple arXiv 2026`** — No new hits. Top results are Nieuviarts 2502.18105, 2512.15450, 2402.05839 (signature change morphism), 1710.04965 (Lorentz signature & twisted spectral triples), all pre-existing and cited in Paper 45.

2. **`Mondino Sämann spectral triple Lorentzian operator algebra 2026`** — No new spectral-triple results. New finding: Mondino–Sämann–Ryborz (May 2026) on stability of timelike Ricci curvature bounds under low-regularity limits — pure synthetic Lorentzian geometry, NOT spectral triples, NOT propinquity. Same category as prior Mondino–Sämann work; remains on synthetic side; no bridge to operator algebra.

3. **`Latrémolière propinquity Krein 2026 arXiv`** — No new Krein/Lorentzian extensions. arXiv:2602.23080 ("Noncommutative coarse metric geometry", Feb 2026) treats proper quantum metric spaces as noncommutative coarse spaces — orthogonal direction, no overlap.

4. **`"Krein spectral triple" convergence Gromov-Hausdorff 2026`** — Top results are pre-existing (Connes–vS 2021, Leimbach–vS 2024). No new Krein-spectral-triple GH-convergence papers found.

5. **`Hekkelman McDonald Leimbach van Suijlekom spectral truncation 2026`** — Hekkelman–McDonald arXiv:2412.00628 accepted to J. Funct. Anal. (July 2025); Leimbach defended PhD December 2025; no new papers from this group in May 2026 window changing Paper 45's scope.

**Verdict: CLEAR — no new direct overlap.** The Lorentzian-propinquity term remains empty in the published literature as of 2026-05-18. The May 2026 preflight memo's `NO_OVERLAP` verdict stands. **Window is still open.**

---

## Track 3 — Test verification

Command:
```
python -m pytest tests/test_lorentzian_propinquity.py tests/test_joint_berezin_compact_temporal.py tests/test_lorentzian_cb_norm.py tests/test_lorentzian_lichnerowicz.py tests/test_lorentzian_propinquity_foundation.py tests/test_krein_positive_state_space.py tests/test_operator_system_lorentzian.py -x --tb=short -q
```

Result: **316 passed, 9 skipped in 317.28s (5min 17s).**

Zero failures. Zero errors. The 9 skipped tests are `@pytest.mark.slow` tests (memory-budget panel cells beyond (4,7)) — expected behavior per memo §7 "memory note." All load-bearing claims in the paper are tested by this suite, including:

- Tunneling pair construction (`test_lorentzian_propinquity.py`)
- Joint Berezin map properties a–e (`test_joint_berezin_compact_temporal.py`)
- cb-norm = 2/(n_max+1) (`test_lorentzian_cb_norm.py`)
- Joint Lichnerowicz [D_L, a_s ⊗ a_t] = i[D_GV, a_s] ⊗ a_t bit-exact (`test_lorentzian_lichnerowicz.py`)
- L1' substrate from Paper 44 (`test_lorentzian_propinquity_foundation.py`, `test_operator_system_lorentzian.py`)
- K⁺ state-space construction (`test_krein_positive_state_space.py`)
- Riemannian-limit recovery at N_t = 1 (in `test_lorentzian_propinquity.py`)

**Test gate: PASS.**

---

## Required fixes (must address before arXiv push)

All five fixes are bibliography-only — no proof, framing, or content changes required.

### FIX 1 — `vandungen2026_su11` author and title

**Current** (lines 1656–1659):
```latex
\bibitem[van~den~Dungen(2026)]{vandungen2026_su11}
K.~van~den~Dungen,
\newblock \emph{Pseudo-Riemannian spectral triples for $\mathrm{SU}(1,1)$},
\newblock arXiv:2601.22171 (January 2026).
```

**Correct attribution:** Author is **Jort de Groot**, not van den Dungen. Title is "Pseudo-Riemannian Spectral Triples for SU(1,1)" (capitalize SU(1,1) correctly per paper).

**Proposed replacement:**
```latex
\bibitem[de~Groot(2026)]{de_groot2026_su11}
J.~de~Groot,
\newblock \emph{Pseudo-Riemannian Spectral Triples for $\mathrm{SU}(1,1)$},
\newblock arXiv:2601.22171, January 2026.
```

The cite key `vandungen2026_su11` is used at line 333 (`\cite{vandungen2026_su11}`); rename to `de_groot2026_su11` at both citation site and bibitem.

### FIX 2 — `mondino_samann2025b` author attribution

**Current** (lines 1588–1592):
```latex
\bibitem[Mondino \& Sämann(2025b)]{mondino_samann2025b}
A.~Mondino and C.~Sämann,
\newblock \emph{Gromov's compactness theorem for the intrinsic
timed-Hausdorff distance},
\newblock arXiv:2510.13069, 2025.
```

**Correct attribution:** Authors are **Mauricio Che, Raquel Perales, and Christina Sormani**, not Mondino & Sämann.

**Proposed replacement:**
```latex
\bibitem[Che, Perales \& Sormani(2025)]{che_perales_sormani2025}
M.~Che, R.~Perales, and C.~Sormani,
\newblock \emph{Gromov's compactness theorem for the intrinsic
timed-Hausdorff distance},
\newblock arXiv:2510.13069, October 2025.
```

The cite key `mondino_samann2025b` is used at line 339 (`\cite{mondino_samann2025,mondino_samann2025b}`); rename to `che_perales_sormani2025` at both citation site and bibitem.

### FIX 3 — `nieuviarts2025b_proceedings` title

**Current** (lines 1600–1605):
```latex
\bibitem[Nieuviarts(2025b)]{nieuviarts2025b_proceedings}
A.~Nieuviarts,
\newblock \emph{Emergence of pseudo-Riemannian structures from twisted
spectral triples within the almost-commutative framework (proceedings
synthesis)},
\newblock arXiv:2512.15450 (v2, May 2026).
```

**Correct title:** "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry" by G. Nieuviarts (initial A. is wrong too — actual first initial is G.).

**Proposed replacement:**
```latex
\bibitem[Nieuviarts(2025b)]{nieuviarts2025b_proceedings}
G.~Nieuviarts,
\newblock \emph{Emergence of time from a twisted spectral triple in
almost-commutative geometry},
\newblock arXiv:2512.15450, December 2025 (revised May 2026).
```

Also fix `nieuviarts2025a` author initial from `A.~Nieuviarts` to `G.~Nieuviarts` for consistency (line 1595):
```latex
\bibitem[Nieuviarts(2025a)]{nieuviarts2025a}
G.~Nieuviarts,
\newblock \emph{Emergence of pseudo-Riemannian spectral triples within
the almost-commutative framework},
\newblock arXiv:2502.18105 (v3, May 2025).
```

### FIX 4 — `bozejko_fendler1991` align with Paper 38

**Current** (lines 1483–1488):
```latex
\bibitem[Bo\.zejko \& Fendler(1984/1991)]{bozejko_fendler1991}
M.~Bo\.zejko and G.~Fendler,
\newblock \emph{Herz--Schur multipliers and completely bounded multipliers
of the Fourier algebra of a locally compact group},
\newblock Boll. Un. Mat. Ital. A (6) \textbf{3} (1984), 297--302;
extended version 1991.
```

**Issue:** Paper 38 cites this differently — see `papers/standalone/paper_38_su2_propinquity_convergence.tex` line 1338–1342:
```
\bibitem[Bo\.zejko \& Fendler(1991)]{bozejko_fendler1991}
M.~Bo\.zejko and G.~Fendler,
\newblock \emph{Herz--Schur multipliers and uniformly bounded
representations of discrete groups},
\newblock Arch.~Math. \textbf{57} (1991), 290--298.
```

For consistency across the math.OA series, **align with Paper 38**:
```latex
\bibitem[Bo\.zejko \& Fendler(1991)]{bozejko_fendler1991}
M.~Bo\.zejko and G.~Fendler,
\newblock \emph{Herz--Schur multipliers and uniformly bounded
representations of discrete groups},
\newblock Arch.~Math. \textbf{57} (1991), 290--298.
```

(The 1984 Boll. Un. Mat. Ital. paper is a different earlier work by the same authors; the 1991 Arch. Math. paper is the canonical cb-norm equality reference.)

### FIX 5 — `latremoliere_metric_st_2017` preprint year

**Current** (lines 1558–1563):
```latex
\bibitem[Latr\'emoli\`ere(2017)]{latremoliere_metric_st_2017}
F.~Latr\'emoli\`ere,
\newblock \emph{The Gromov--Hausdorff propinquity for metric spectral
triples},
\newblock Adv.~Math. \textbf{415} (2023), Paper No.~108876, 88pp;
\newblock preprint arXiv:1811.10843, 2017.
```

**Issue:** arXiv:1811.10843 was first submitted November 27, 2018, not 2017. The published Adv. Math. 2023 reference is correct.

**Proposed replacement:**
```latex
\bibitem[Latr\'emoli\`ere(2018/2023)]{latremoliere_metric_st_2017}
F.~Latr\'emoli\`ere,
\newblock \emph{The Gromov--Hausdorff propinquity for metric spectral
triples},
\newblock Adv.~Math. \textbf{415} (2023), Paper No.~108876, 88pp;
\newblock preprint arXiv:1811.10843, November 2018.
```

(The cite key `latremoliere_metric_st_2017` may be retained for cross-paper consistency or renamed to `latremoliere_metric_st_2018`; recommend retaining the key to avoid breaking other citation files.)

---

## Recommended polish (acceptable to defer)

1. **Latrémolière 2017 Def. 3.5 inline citation** in §5 Proposition 5.2 — add parenthetical `(Def. 3.5 of arXiv:1811.10843)` for reader convenience.
2. **Explicit "no published Lorentzian propinquity" sentence** in §1.4 (currently implicit). One sentence: "To our knowledge, no Lorentzian propinquity has been defined in the published math.OA literature as of May 2026."
3. **Wikipedia-style author redirect note** — none needed; arXiv handles redirects.
4. **Paper 38 cross-reference at line 174** — "first Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature" — defensible; the qualifier "K⁺-restricted weak-form" downstream makes the scope clear.

---

## Verdict

**REQUIRES_FIXES** at the bibliography level. Five mechanical fixes (FIX 1–5) listed above. Once applied:

- **All proofs** verified against source memos (Sub-sprints A, B, C, D).
- **All numerical claims** verified against `debug/data/l3b_2_sub_sprint_D.json`.
- **All 316 tests pass** (plus 9 skipped slow-mark tests, expected).
- **No concurrent-work scoop** as of 2026-05-18.
- **K⁺-weak-form framing** load-bearing throughout (abstract, §1, §7); not softened by any track.
- **All internal cross-references** (Papers 32, 38, 39, 40, 42, 43, 44) verified consistent.

After FIX 1–5 are applied, **Paper 45 is ARXIV_READY** pending PI metadata sign-off (math.OA primary, math-ph + gr-qc secondary).

---

**End of memo.**
