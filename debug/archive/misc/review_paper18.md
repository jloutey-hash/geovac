# Confidence Review: Paper 18 — Spectral-Geometric Exchange Constants

**Auditor:** GeoVac Confidence Reviewer (Combined Auditor + Citation Checker)
**Date:** 2026-06-02
**Paper file:** `papers/group3_foundations/paper_18_exchange_constants.tex` (3,135 lines)
**Status verdict:** **YELLOW** (broadcastable after one HIGH math-internal-consistency fix; medium framing tightening for Ω convention)

---

## Calibration check

Not a calibration run.

---

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|-------|----------|---------|----------|---------------------|
| 1 | Six-tier classification (intrinsic / calibration / embedding / algebraic-implicit / composition / inner-factor input data) | §IV `sec:taxonomy` + line 1778 §IV "Inner-factor input data" | **B** — internally consistent | GEOVAC-ONLY (Paper 32 §VIII.C consistent; Marcolli-vS lineage cited externally) | Cross-checked Paper 32 §VIII; both papers list M1/M2/M3 and the inner-factor sixth tier consistently |
| 2 | Master Mellin engine §III.7: M1/M2/M3 as 𝓜[Tr(D^k·e^{−tD²})] at k∈{0,1,2} | §III.7 Eq.(20) `eq:master_mellin` | **A** — externally consistent / internally derived | MIXED (Mellin–heat-kernel relation is textbook Seeley; the M1/M2/M3 partition is GeoVac-internal) | Paper 32 §VIII Theorem 1 states identical M1/M2/M3 case-exhaustion; consistent across corpus |
| 3 | κ = −1/16 = −1/Ω⁴(0) with Ω(0) = 2 | §VII.6 `sec:kappa_derivation` line 2644; §IV item 3 line 1101 | **A** for κ = −1/16; **C** for κ = −1/Ω⁴(0) under §II.1's stated Ω | EXTERNAL (Fock 1935 / Paper 7 conformal factor) + Paper 0 packing geometry | Lines 2644–2650: graph-degree derivation and Ω derivation agree (8 = 2·d_max = Ω(0)⁴). VALID under §VIII convention Ω = 2p₀/(p² + p₀²). **INTERNAL INCONSISTENCY**: §II.1 line 173 defines Ω = (1 + p²/p₀²)⁻¹ giving Ω(0) = 1, under which −1/Ω⁴(0) = −1, NOT −1/16. See HIGH finding 1 below. |
| 4 | Theorem 1: ζ_{D²}(s) = 2^{2s−1}[λ(2s−2) − λ(2s)] for s ≥ 2; π^even content exclusively (Eq.29) | §IV Theorem 1 part (1) | **A** | EXTERNAL (Camporesi-Higuchi spectrum + standard Bernoulli identities) | Verified symbolically: ζ_{D²}(2) = π² − π⁴/12, matches Paper 28 line 455 D(4) value |
| 5 | Theorem 1 Eq.30: D(s) = 2(2^{s−2}−1)ζ_R(s−2) − (2^s−1)/2 · ζ_R(s) for s ≥ 2 | §IV Theorem 1 part (2) Eq.(30) | **A** for formula | EXTERNAL (half-integer Hurwitz arithmetic) | Numerically: D(4) = 1.7521... = π² − π⁴/12, **same as ζ_{D²}(2)**, matches Paper 28 Eq.(43) line 455. Paper 18 internal: Eq.30 verified. |
| 6 | Theorem 1 Eq.30 corollary: D(4) = 2ζ_R(2) + 2ζ_R(3); contains "odd-zeta content not in any value of ζ_{D²}" (line 1350) | §IV line 1350 + Eq.38 `eq:dirac_dirichlet_zeta3` | **E — WRONG** | GEOVAC-ONLY chain (Eq.38 re-derivation) | Numerical sum over CH spectrum (|λ_n|=n+3/2, g=2(n+1)(n+2), n≥0) gives **D(4) = π² − π⁴/12 ≈ 1.7521**, NOT 2ζ(2)+2ζ(3) ≈ 5.6940 as claimed. Paper 28 Eq.(43) line 455 has the correct value. Paper 18 has an internal contradiction inside Theorem 1: the closed-form Eq.30 gives π² − π⁴/12, but line 1350 explicit substitution claims 2ζ(2)+2ζ(3). See HIGH finding 2 below. |
| 7 | Theorem 2: Three-axis classification (operator order × bundle × diagram topology) | §IV Theorem 2 | **B** | GEOVAC-ONLY (taxonomic framing built on Theorem 1 + Paper 28's χ_{−4} mechanism) | Hurwitz quarter-shift mechanism for χ_{−4} via β(s) − β(s−2) closed form is independently verified in Paper 28 Theorem 3 |
| 8 | Theorem 3: B = 42 (finite Casimir), F = π²/6 (Fock-degeneracy Dirichlet at d_max), Δ = 1/40 (Dirac mode degeneracy) | §V Theorem 3 (`thm:k_decomposition`) | **A** for each component value, **C** for "structurally independent" framing | MIXED (B and Δ are EXTERNAL/computational; F is GEOVAC-specific reading of D_{n²}(4) = ζ(2); independence is GEOVAC-ONLY observation) | Verified symbolically: B = Σ_{n=1..3} Σ_{l=0..n−1} (2l+1)l(l+1) = 42; g_3^Dirac = 2·4·5 = 40; D_{n²}(s) = ζ(s−2), at s=4 gives ζ(2) = π²/6 ✓ |
| 9 | κ–B identity Eq.(62): Σ (2l+1)l(l+1)⟨n,l|w|n,l⟩ = 63/4 = 6·B·|κ| | §V.1 `eq:kappa_B_identity` | **A** | GEOVAC-internal computation but symbolically derivable | Verified exact: 63/4 = 6·42·(1/16) ✓ (sympy) |
| 10 | Combination rule K = π(B + F − Δ) = α⁻¹ is conjectural (Conjecture 1) | §V Conjecture 1 `conj:k_rule` | **D** — unverifiable here (an empirical observation, openly framed as such) | EXTERNAL CODATA α value + GEOVAC structural reading | Paper correctly labels as conjectural, preserving §13.5 prohibition |
| 11 | Theorem 4 (η-trivialization): {γ_F, D_F} = 0 ⇒ Tr(D_F^k · e^{−tD_F²}) ≡ 0 for odd k | §IV (Inner-factor) Theorem 4 line 1791 | **A** | EXTERNAL (elementary trace identity under involution) | Numerically verified with explicit 4×4 example: residual 0.0 |
| 12 | Theorem 5 (AC factorization): Tr(e^{−tD²}) = Tr(e^{−tD_GV²}) · Tr(e^{−tD_F²}) | §IV Theorem 5 line 1815 | **A** | EXTERNAL (standard tensor-product spectral theory + outer chirality anticommutation hypothesis) | Numerically verified |
| 13 | Operator-order discriminant motivic-weight reading | §IV "Motivic weight and operator order" | **C** — taxonomically packaged, but the "loop order" prediction depends on Eq.38's flawed re-derivation | GEOVAC-ONLY taxonomic reading | The two-loop ζ(3) is real physics (Petermann–Sommerfield); whether it follows from the GeoVac D(s) discriminant requires the corrected D(s) value |
| 14 | "Inner-factor input data" sixth tier disjoint from outer-factor M1/M2/M3 | §IV line 1778+ | **B** | GEOVAC-ONLY (Sprint H1 + Paper 32 §VIII.C consistent) | Logical: Theorems 4 + 5 + Yukawa-Dirichlet formula constitute an internal proof framework |
| 15 | Compactness thesis (discrete = compact, transcendentals via Peter-Weyl projection) | §III intro + §III.2 | **B** | EXTERNAL (Peter-Weyl 1927, Selberg, Weyl) | Mathematically sound; the GeoVac-specific use is consistent |
| 16 | Claim 4: observable classification by projection type (Class S/P/C) | §VI `sec:observable_classification` | **B** | GEOVAC-internal but falsifiable as claimed | Examples (hydrogen Class S, density Class P, He Class C) hold |
| 17 | Weyl's law: N(λ) ~ ω_d/(2π)^d · vol · λ^{d/2} | §III.1 Eq.(15) | **A** | EXTERNAL (Weyl 1911) | Verified equivalent to standard form (4π)^{−d/2}/Γ(d/2+1) · vol · λ^{d/2} symbolically |
| 18 | SD coefficients on unit S³: a_0 = 2π², a_1 = 2π², a_2 = π² (scalar Laplacian) | §VII.6 line 2625 | **A** | EXTERNAL (standard heat kernel on round S³) | a_0 = Vol(S³) = 2π² ✓; a_2 = (R/6)Vol with R = 6 → π² ✓ |
| 19 | The graph is π-free (Observation 1) | §III.2 line 640 | **A** | EXTERNAL (Peter-Weyl + Gaunt 3j arithmetic) | Sound classical statement |
| 20 | RH-J spectral χ_{−4} identity: D_even(s) − D_odd(s) = 2^{s−1}(β(s) − β(s−2)) | §II.5 Eq.(11) `eq:rh_j_identity` | **A** | EXTERNAL (Hurwitz identities) + Paper 28 cross-checked | Sanity-checked at s = 4: matches Paper 28 D_even/D_odd formulas Eqs.(31–32) |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value / error | Survives? |
|---|---|---|---|---|
| B = Σ_{n=1..3}Σ_{l=0..n−1}(2l+1)l(l+1) | 42 | Direct sympy enumeration | 42 (exact) | ✓ |
| Eq.(62): Σ(2l+1)l(l+1)·w(n,l) | 63/4 | Sympy with stated weight values | 63/4 (exact) | ✓ |
| ζ_{D²}(2) closed form | π² − π⁴/12 | Independent symbolic sympy | π² − π⁴/12 (exact) | ✓ |
| Eq.(30) D(s) at s = 4 | (paper implies 2ζ(2) + 2ζ(3) per Eq.38 / line 1350) | Direct mpmath sum over CH spectrum (g_n/|λ_n|^4 with g_n=2(n+1)(n+2), |λ_n|=n+3/2, n=0,1,2,...) at 50 dps | **1.7521801482558... = π² − π⁴/12**, matches Paper 28 Eq.(43) line 455 | **FAILS** — Paper 18 Eq.38 and the explicit substitution at line 1350 do NOT survive against the direct CH-spectrum sum or Paper 28's stated value |
| 2ζ(2) + 2ζ(3) (paper's claimed D(4) value) | 5.6940 | mpmath | 5.6940 (correct as a number, wrong as D(4)) | The number is right; the identification with D(4) is wrong |
| Δ⁻¹ = g_3^Dirac = 2(n+1)(n+2)|_{n=3} | 40 | Direct evaluation | 40 (exact) | ✓ |
| D_even(4) − D_odd(4) = 8(β(4) − G) (consistency with Eq.(11) at s = 4 with 2^{s−1}=8 and β(s−2)=β(2)=G) | — | Sanity from Paper 28 Eqs.(31–32) | 8β(4) − 8G ✓ | ✓ |
| κ = −1/16 (graph-degree) and κ = −1/Ω⁴(0) with Ω(0) = 2 | (consistency between two derivations) | Paper 0 line 638; Paper 7 line 816 | Both = −1/16 under Ω(0) = 2; **fail to agree** under §II.1 line 173 convention where Ω(0) = 1 | Internal inconsistency about which Ω is meant; see HIGH finding 1 |
| Weyl ball volumes for d = 1,2,3,5 | 2, π, 4π/3, 8π²/15 | Sympy π^{d/2}/Γ(d/2+1) | 2, π, 4π/3, 8π²/15 (exact) | ✓ |

### Circularity map

The GEOVAC-ONLY chains:

- **Six-tier taxonomy structure** (Claim 1) rests entirely on the GeoVac corpus (Papers 2, 13, 17, 32). Each individual transcendental cited (κ, e^a·E₁(a), μ(R), PK, inner-factor input data) is GeoVac-defined; the bottom of the chain is the GeoVac framework itself.
- **Master Mellin engine domain partition** (Claim 2 refined): M1, M2, M3 as "three sub-cases of a single master mechanism 𝓜[Tr(D^k e^{−tD²})] at k ∈ {0,1,2}" is GEOVAC-ONLY — Paper 32 §VIII gives the case-exhaustion theorem, but no external NCG paper makes this exact partition. The underlying Mellin/heat-kernel/Seeley machinery is EXTERNAL; the M1/M2/M3 labelling and case-exhaustion claim are internal.
- **Theorem 3 "three independent spectral homes for B, F, Δ"** (Claim 8) bottoms out on GeoVac's claim that these are the right spectral objects for Paper 2's K. Each component value (42, π²/6, 1/40) is verifiable by external math, but the framing that they live in "structurally independent" cells is internal taxonomic interpretation.
- **Inner-factor input data sixth tier** (Claims 11, 12, 14) is GEOVAC-ONLY in framing (the partition into inner vs outer factor is the GeoVac reading of the Connes–Marcolli–vS framework). The η-trivialization Theorem 4 itself is a clean external result.

External anchors (EXTERNAL):

- Weyl's law (Weyl 1911), Selberg trace formula (Selberg 1956), Peter–Weyl (Peter–Weyl 1927).
- Camporesi–Higuchi spectrum (gr-qc/9505009, J. Geom. Phys. 1996) — verified.
- Krajewski 1998 finite-spectral-triple classification (hep-th/9701081) — verified.
- Paschke–Sitarz discrete spectral triples (q-alg/9612029, J. Math. Phys. 39, 1998) — verified.
- Lichnerowicz formula D² = ∇*∇ + R/4.
- Apéry 1978 (ζ(3) irrationality); Goncharov, Zagier motivic framework.
- Petermann 1957 / Sommerfield 1957 (two-loop ζ(3) in a_e).
- DLMF 18.9.13 (Laguerre summation).

### Overstatement findings

- **§II.1 Ω-normalization (line 173)**: defines Ω = (1 + p²/p₀²)⁻¹, giving Ω(0) = 1, which is **inconsistent with Ω(0) = 2 used downstream** at lines 1101 and 2644 to give κ = −1/Ω⁴(0) = −1/16. Suggested fix: in §II.1 use the same Ω convention as §VIII and Paper 7, Ω = 2p₀/(p² + p₀²); the κ identification then matches. Alternatively keep §II.1's convention but state explicitly that the lines 1101 and 2644 derivations rest on the §VIII/Paper-7 convention with Ω(0) = 2 (which equals 2/p₀ at p₀ = 1).
- **§IV Eq.(38) `eq:dirac_dirichlet_zeta3`**: states "D_{g^Dirac}(4) = Σ_{m≥1} 2m(m+1)/m^4 = 2ζ(2) + 2ζ(3)". This conflates an integer m (treated as both eigenvalue and degeneracy parameter) with the half-integer Camporesi–Higuchi spectrum |λ_n| = n + 3/2. The actual CH-Dirac D(4) is π² − π⁴/12, not 2ζ(2) + 2ζ(3); this is what Paper 28 Eq.(43) line 455 correctly states. The "odd-zeta content from first-order spectral operators" reading is then unsupported AT s = 4 (the cited equation produces no ζ(3)). ζ(3) DOES appear elsewhere on the spinor bundle (e.g. half-integer Hurwitz hydrogen sums in Paper 28 §X), so the motivic-weight thesis isn't false in general — but Eq.(38)'s witness for it is wrong. Suggested fix: replace Eq.(38) with a correct source of ζ(3) (e.g. the half-integer Hurwitz hydrogen sum identity in Paper 28 §X), or remove the equation and replace the surrounding "odd-zeta content from D(s)" framing with the correct CH closed form.
- **Theorem 1 part (2) sentence "D(4) = 2ζ_R(2) + 2ζ_R(3) contains odd-zeta content not present in any value of ζ_{D²}"** (line 1350) is the textual restatement of Eq.(38); same fix.

---

## Pass B — Citation and novelty

### Citation table

| \cite key | Claimed as | Verdict | What I found |
|---|---|---|---|
| loutey_paper0,1,7,11,13,2,22,24,25,27,28,29,30,31,32,fci | Internal | CITE-OK | All bibitems exist in the corpus per the bibliography section; internal cross-corpus consistency verified for paper28, paper32 above. |
| loutey_track_g, _j, _m, _p1, _s | Internal track memos | CITE-OK (internal) | Memos referenced; supplementary materials. |
| weyl1911 | Weyl's eigenvalue distribution paper, Nachr. Gesell. Wiss. Göttingen 1911 | CITE-OK | Standard reference. |
| selberg1956 | "Harmonic analysis and discontinuous groups..." J. Indian Math. Soc. 20 (1956) | CITE-OK | Standard reference. |
| peter_weyl1927 | "Die Vollständigkeit der primitiven Darstellungen einer geschlossenen kontinuierlichen Gruppe," Math. Ann. 97 (1927) | CITE-OK | Standard reference. |
| seeley1967 | "Complex powers of an elliptic operator," Proc. Sympos. Pure Math. 10 (1967) | CITE-OK | Standard reference. |
| minakshisundaram1949 | Minakshisundaram & Pleijel, Canad. J. Math. 1 (1949) | CITE-OK | Standard reference. |
| gilkey1975 | "Spectral geometry of a Riemannian manifold," J. Diff. Geom. 10 (1975) | CITE-OK | Standard reference. |
| lichnerowicz1963 | Spineurs harmoniques, C. R. Acad. Sci. Paris 257 (1963) | CITE-OK | Standard reference. |
| camporesi_higuchi1996 | J. Geom. Phys. 20, 1–18 (1996) — Dirac operator on spheres | CITE-OK | Web-verified: J. Geom. Phys. 20 (1996), DOI 10.1016/0393-0440(95)00042-9, arXiv gr-qc/9505009. |
| krajewski1998 | "Classification of Finite Spectral Triples," J. Geom. Phys. 28, 1–30 (1998), arXiv hep-th/9701081 | CITE-OK | Web-verified: author, title, journal, year, pages, arXiv all match. |
| paschke_sitarz2000 | "Discrete Spectral Triples and their Symmetries," J. Math. Phys. 39, 6191–6205 (1998), arXiv q-alg/9612029 | CITE-WRONG-METADATA (cosmetic) | Web-verified the bibitem TEXT is correct (J. Math. Phys. 39, 6191, 1998). The bibitem KEY says `paschke_sitarz2000` suggesting year 2000; correct year is 1998. LOW-severity cosmetic key mismatch. |
| latremoliere2018 | "The Dual Modular Propinquity and Completeness," arXiv:1811.04534; J. NCG 15, 347 (2021) | CITE-WRONG-METADATA (title shortening) | Web-verified: arXiv 1811.04534 and J. NCG vol 15 (2021) match. Full title is "The Dual Modular **Gromov–Hausdorff** Propinquity and Completeness" — paper omits "Gromov–Hausdorff". LOW-severity title shortening. |
| goncharov1999 | "Volumes of hyperbolic manifolds and mixed Tate motives," J. Amer. Math. Soc. 12 (1999) | CITE-OK | Standard reference. |
| zagier1994 | "Values of zeta functions and their applications," First European Congress 1994 | CITE-OK | Standard reference. |
| schwinger1948 | Phys. Rev. 73, 416 (1948) | CITE-OK | Standard reference. |
| petermann1957 | Helv. Phys. Acta 30, 407 (1957) | CITE-OK | Standard reference. |
| sommerfield1957 | Phys. Rev. 107, 328 (1957) | CITE-OK | Standard reference. |
| berger2003 | "A Panoramic View of Riemannian Geometry," Springer 2003 | CITE-OK | Standard reference. |
| polchinski1998 | "String Theory" Vol. 2, CUP 1998 | CITE-OK | Standard reference. |
| dlmf_laguerre | DLMF §18.9.13 | CITE-OK | URL in bibitem resolves to NIST DLMF. |
| drake2006 | Drake handbook chapter, Springer 2006 | CITE-OK | Standard reference. |
| adamchik2003 | Z. Anal. Anwend. 22 (2003) | (defined but **unused** in body — bibliography clutter) | LOW |

### Problems found (CITE-MISATTRIBUTED / DOESNT-SUPPORT / CANT-FIND)

- **None of the load-bearing citations is misattributed or fabricated.** All web-checked items (Camporesi-Higuchi, Krajewski, Paschke-Sitarz, Latrémolière) have correct arXiv IDs / journal venues. The Paschke-Sitarz key mismatch (key says "2000", text correctly says 1998) and the Latrémolière title shortening are LOW-severity cosmetic items.
- **Unused bibitems** (bibliography clutter): `adamchik2003`, `loutey_paper22`, `loutey_paper27`, `loutey_track_s` are defined but never `\cite`-d in the body. LOW severity.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "To our knowledge (\ref{eq:F_dirichlet}) is the first identification of F from a graph-intrinsic quantity in the GeoVac framework" (line 2123) | §V.2 | Internal-scope claim ("within GeoVac framework"), not external priority | N/A — internal scope | No action; claim is properly scoped. |
| (Several "first GUE signature in GeoVac" type claims in §II.5) | §II.5 | Internal-scope claim | N/A — internal scope | No action; claim is properly scoped. |

The paper does NOT make strong external priority claims (e.g. "first in the literature"). The "to our knowledge" hedge at line 2123 is properly scoped to "within the GeoVac framework."

---

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---|---|---|---|
| H1 | §II.1 line 173 defines Ω = (1 + p²/p₀²)⁻¹ giving Ω(0) = 1, which contradicts the §VIII line 2364 / Paper 7 / Paper 0 convention Ω = 2p₀/(p² + p₀²) giving Ω(0) = 2 used downstream at lines 1101 and 2644 in the κ = −1/Ω⁴(0) derivation. Under §II.1 convention the κ derivation gives −1, not −1/16. | A | Internal inconsistency / overstatement | **HIGH** (load-bearing for κ-derivation; reader-confusing) — recommend harmonising to one Ω convention paper-wide |
| H2 | Eq.(38) `eq:dirac_dirichlet_zeta3` and the substitution at line 1350 inside Theorem 1 part (2) state D(4) = 2ζ(2) + 2ζ(3) (≈ 5.694), but the closed form Eq.(30) and direct CH-spectrum summation both give D(4) = π² − π⁴/12 (≈ 1.752), matching Paper 28 Eq.(43) line 455. The Eq.(38) sum Σ_{m≥1} 2m(m+1)/m^4 conflates an integer m with the half-integer CH eigenvalues |λ_n| = n + 3/2. | A | Math error E | **HIGH** (load-bearing for the "odd-zeta content from first-order Dirac D(s)" reading that supports the motivic-weight thesis of §IV) |
| M1 | The motivic-weight loop-order prediction at §IV (lines 1444–1471) and the corollary in §IX paragraph "Theorem 1" rely on Eq.(38)'s flawed identification of D(4) with ζ(3); the ultimate empirical match to Petermann–Sommerfield is real physics but the GeoVac D(s) witness for it must be replaced by a correct ζ(3) source on the spinor bundle (e.g. the half-integer Hurwitz hydrogen sums in Paper 28 §X, line 2363+). | A | Overstatement C downstream of H2 | **MEDIUM** |
| M2 | Latrémolière 2018 title in bibitem `latremoliere2018` shortens "The Dual Modular Gromov–Hausdorff Propinquity and Completeness" to "The Dual Modular Propinquity and Completeness". Web-verified the long form is correct. | B | CITE-WRONG-METADATA | **MEDIUM** (cosmetic but verifiable) |
| L1 | Bibitem key `paschke_sitarz2000` suggests year 2000; correct year per J. Math. Phys. 39, 6191 (1998) is 1998. Bibitem TEXT is correct; only the key is misleading. | B | CITE-WRONG-METADATA (cosmetic) | **LOW** |
| L2 | Unused bibitems (clutter): `adamchik2003`, `loutey_paper22`, `loutey_paper27`, `loutey_track_s` are defined but never cited in body. | B | Bibliography clutter | **LOW** |
| L3 | The abstract's promotion of `(intrinsic, conformal/calibration, embedding, flow, composition)` is FIVE-tier, but the body (§IV + inner-factor §IV.6 + algebraic-implicit §IV.5 + composition) is six-tier per CLAUDE.md §6 (`intrinsic / calibration / embedding / algebraic-implicit / composition / inner-factor input data`). Abstract is one revision behind. | A | Overstatement C (under-stated, not over-stated) | **LOW** |
| L4 | §I "Claim 4" is introduced inline at line 112 but the abstract describes it as "the falsifiable claim" without a number; the formal Claim 4 statement appears in §VI at line 2310. Cross-reference is fine. | A | LOW | **LOW** (no fix needed unless the paper is re-paginated) |

---

## Broadcast readiness: **YELLOW**

Paper 18 is one of the most-cited papers in the GeoVac corpus, and the master Mellin engine in §III.7 plus the three theorems in §IV–V hold up to independent symbolic and numerical verification — with two important exceptions: (a) an internal Ω-normalization inconsistency between §II.1 (Ω = (1 + p²/p₀²)⁻¹, Ω(0) = 1) and §VIII / Paper 7 / Paper 0 (Ω = 2p₀/(p² + p₀²), Ω(0) = 2), where downstream κ-derivation rests on the latter convention; and (b) a math error inside Theorem 1 part (2) where Eq.(38) `eq:dirac_dirichlet_zeta3` claims D(4) = 2ζ(2) + 2ζ(3) but the formula Eq.(30) immediately above (and the cross-checked Paper 28 Eq.(43)) both correctly give D(4) = π² − π⁴/12. The error is a conflation of integer m with the half-integer Camporesi–Higuchi spectrum; it propagates into the motivic-weight loop-order argument in §IV. Both fixes are small text-level edits — H1 is a one-line convention swap, H2 is a one-equation replacement plus a one-sentence textual correction at line 1350. Citations are clean (no fabricated arXiv IDs); the master Mellin engine matches Paper 32 §VIII case-exhaustion theorem bit-for-text; the Ω-normalization issue affects only Paper 18 itself (Papers 0 and 7 use the consistent convention). After H1 + H2 are fixed, broadcast verdict would be GREEN.

---

## What I could NOT verify (hand to a human expert)

- The "structural incommensurability" Observation in §VIII is a meta-pattern claim. Verifying it requires a domain expert who can audit ALL THREE instances (QED-on-S³, K-decomposition, Weil dictionary) simultaneously and judge whether the framing "no common generator" is a deep structural statement or a tautology of how the cells were defined.
- The κ-B identity Eq.(62) is verified at finite n_max = 3, but whether 1/16 = 1/(2 d_max) really IS the same constant as 1/Ω⁴(0) by deep geometry (Sec. VII.6 last paragraph) — vs being a coincidence at the finite-cutoff — is the kind of structural claim only a Connes-style NCG expert can adjudicate.
- The "Dirichlet-L tier" classification of β(s) values as a distinct cell (§IV "Remark: Dirichlet L-tier as new transcendental class") rests on the Hurwitz quarter-shift mechanism being non-trivially distinct from M2's Bernoulli mechanism. Mathematically the two output rings differ; whether that constitutes a TIER in the spectral-action sense is a matter of taxonomy convention.
- The "inner-factor input data" sixth tier (Theorems 4, 5) is internally consistent and matches Paper 32 §VIII.C, but whether the η-trivialization theorem (a clean algebraic fact) really licenses excluding entire classes of inner Dirac operators is a question for the NCG community (Marcolli–vS, Perez-Sanchez literature).

---

## Summary line for dispatcher

Paper 18 (Spectral-Geometric Exchange Constants), **YELLOW**: A-count 9, B-count 6, C-count 2, D-count 1, E-count 1, CITE-OK 22, CITE-WRONG-METADATA 2, CITE-MISATTRIBUTED 0, CITE-DOESNT-SUPPORT 0, CITE-CANT-FIND 0; HIGH 2, MEDIUM 2, LOW 4. Top finding: **Theorem 1 part (2) contains an internal contradiction — Eq.(30) and Paper 28 give D(4) = π² − π⁴/12, but Eq.(38)/line 1350 claim D(4) = 2ζ(2) + 2ζ(3); the Eq.(38) sum conflates an integer m with the Camporesi–Higuchi half-integer spectrum.** Secondary structural finding: **§II.1 Ω = (1 + p²/p₀²)⁻¹ (Ω(0) = 1) is inconsistent with §VIII / Paper 7 / Paper 0 Ω = 2p₀/(p² + p₀²) (Ω(0) = 2) used downstream in the κ = −1/Ω⁴(0) derivation; recommend harmonising to the §VIII / Paper 7 convention.**
