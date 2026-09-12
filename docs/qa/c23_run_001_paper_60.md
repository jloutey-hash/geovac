# C23 run #001 — Paper 60 (2026-09-12)

**The first exercise of the inverse-citation criterion** (`docs/qa/criteria.md`, C23),
run against the paper whose `eq:sigma_law` created it. Not a certification; a
single-criterion pass. Scans: `debug/lit_scan/c23_paper60_{linalg,analysis}_memo.md`.

**Headline: 6 of 8 audited claims returned PRIOR ART, including one written the
day before the run.** No claim was found false — C23 re-tiers attribution, not truth.
Every load-bearing identity was re-verified numerically here before any edit, per the
criterion's own first hard rule.

## Verdicts

| # | Claim as the paper had it | Verdict | Disposition |
|:--|:--|:--|:--|
| A1 | `spec[[I,C],[C^T,I]] = {1±σ}`, `cond = (1+σ)/(1−σ)`, σ = cosines of principal angles | **PRIOR ART** | Jordan's principal angles; Jordan–Wielandt spectrum; two-block CBS constant. Named in prose as classical, used not derived. Primaries **owed**. |
| A2 | `‖[P_A,P_B]‖ = max σ√(1−σ²)`, cited to Halmos | **PRIOR ART** | Halmos primary re-read: Thm 2 is the canonical form *only*, no norm. Loring 2014 added (verified); the inline derivation stands. |
| A3 | "Proposition D" — no block-diagonal congruence orthogonalizes a non-block-diagonal metric | **PRIOR ART (content)** / ABSENT (exact contrapositive) | The catch. Its positive direction is the Löwdin symmetry-preservation property, known since Slater–Koster (1954); in operator terms, block-diagonals are a commutant and inverse-closed. **Demoted** from Proposition to a specialization; the ℓ-vs-m application is what the paper now claims. Primary **owed**. |
| A4 | `2/(1+min j₀) = 2.555041…` | **ABSENT** (constant) | No name, no tabulation. But the *mechanism* is classical (self-adjoint Toeplitz spectrum = convex hull of the symbol's essential range) and is now named. Primary **owed**. |
| B1 | `\|c_j\| ~ j^{-5/4}` by stationary phase | **PRIOR ART (chain)** / ABSENT (statement) | The asymptotic needed no deriving: the model integral is DLMF 10.32.10 at ν=2, and DLMF 10.40.2 delivers constant *and* phase. Claim **upgraded** from an exponent to `eq:chirp_decay` in full. Two corrections fell out (below). DLMF cited (verified verbatim with phase conditions). |
| B2 | `j₀` composed with a cotangent, coefficient asymptotics | **ABSENT** | Nobody has done it. Recorded trap: "Toeplitz operators on the *Fock space*" in the OA literature is Bargmann–Fock, not Fock's projection — it makes the search look populated when it is not. |
| B3 | the `S^{-1/2}` non-locality mechanism | **PRIOR ART, and our reading is correct** | Jaffard's class: polynomial decay of order >1 plus bounded invertibility ⟹ decay-preserving inversion; inverse-closedness makes holomorphic and Riesz calculi agree. So the escape is the spectrum touching zero, **not** the decay class. Stated in the paper. Primaries **owed**. |
| B4 | rank-`M−1` all-ones degeneracy at `χ=π` (written 2026-09-12) | **PRIOR ART** | This is the **flat limit** of the RBF/kernel literature (Barthelmé–Usevich 2021, verified; lineage Driscoll–Fornberg 2002, Schaback 2005). What survives as ours: that the block symbols realize it at a *symbol point* rather than a shape-parameter limit — which is why a fixed rank-`M−1` rotation removes it. |

## Two corrections to our own reasoning, from B1

1. **The `π/4` is a branch phase**, from the `(π/2z)^{1/2}` prefactor of DLMF 10.40.2 — *not* the stationary-phase `sign(φ'')·π/4`. Same number, wrong mechanism. The backing guard now pins the phase by sign agreement, so the right reason is tested rather than the right value.
2. **`Σ|c_j|` converges.** `5/4 > 1`, so the symbol is in the Wiener algebra; what diverges is `Σ j|c_j|`, which is Böttcher–Widom's hypothesis and a different condition. The paper said the right thing but the distinction was not drawn, and it matters for B3.

## Citations added (verified at primary source in this session, and only these)

- **DLMF 10.32.10 / 10.40.2** — quoted verbatim with their phase conditions (`|ph z| < π/4`; `|ph z| ≤ 3π/2 − δ`).
- **Loring**, *Principal angles and approximation for quaternionic projections*, Ann. Funct. Anal. **5**(2), 176 (2014), arXiv:1306.1923 — title, venue and topic verified. The `‖PQ−QP‖ = ½sin 2θ` line was read by the scan as a displayed proof line, not a numbered theorem; cited accordingly.
- **Barthelmé & Usevich**, *Spectral properties of kernel matrices in the flat limit*, SIAM J. Matrix Anal. Appl. **42**(1), 17 (2021), arXiv:1910.14067 — title, venue and topic verified.

## OWED — named in prose, deliberately given no bibitem

The session's WebSearch budget (200) was exhausted partway through; `WebFetch` remained, so anything with a known URL was reachable and anything requiring a *search* was not. Per C23's second hard rule these were **not** invented:

| Owed primary | For |
|:--|:--|
| Jordan (1875), *Bull. Soc. Math. France* **3**, 103 | principal angles (A1) |
| Jordan–Wielandt form; Stewart–Sun / Horn–Johnson | block spectrum (A1) |
| Eijkhout–Vassilevski, *SIAM Review* **33**, 405 (1991) | CBS constant (A1) |
| Slater & Koster, *Phys. Rev.* **94**, 1498 (1954), Appendix | symmetry preservation (A3) |
| Rokob–Szabados–Surján, note on Löwdin orthogonalization symmetry | the `[T,S]=0` hypothesis (A3) |
| Hartman–Wintner | Toeplitz spectrum = convex hull of essential range (A4) |
| Jaffard (1990), *Ann. IHP Anal. Non Linéaire* **7**, 461 | the decay class (B3) |
| Gröchenig–Leinert, *TAMS* **358**, 2695 (2006) | matrix inverse-closedness (B3) — **note:** the JAMS 2004 paper is the *Gabor* one and is the wrong cite for matrices |
| Driscoll–Fornberg (2002) | coined "flat limit" (B4) |
| Björck–Golub (1973) | principal-angle computation (A1) — AMS PDF returned 403 |
| Böttcher–Spitkovsky (2010) | **already cited by Paper 60 for content nobody here has read** — paywalled, no preprint. Its title matches its use (two-projections theory), so it is retained rather than dropped; flagged as unverified. |

**CLOSED 2026-09-12 (v5.11.2).** The WebSearch budget was raised (project `env`, 200 -> 500) and took effect immediately. Ten of the eleven owed primaries were verified at source and are now cited: Jordan 1875, Björck–Golub 1973, Eijkhout–Vassilevski 1991, Hartman–Wintner 1954, Löwdin 1950, Slater–Koster 1954 (**primary read** — Sec. II, p. 1500), Rokob–Szabados–Surján (existence + abstract; PDF cert-blocked), Jaffard 1990, Gröchenig–Leinert 2006, Driscoll–Fornberg 2002. **Two items were deliberately not closed:** the "Jordan–Wielandt" label stays prose-only (Stewart–Sun and Horn–Johnson are books, unopened), and the scan's claim that Slater–Koster's *Appendix* states the symmetry theorem is **dropped, not cited** — what was read is their Sec. II use of Löwdin, which is what the citation carries. Böttcher–Spitkovsky remains cited-but-unread and flagged. A1/A3/A4/B3 verdicts stand as PRIOR ART, now with primaries.

## Process findings

- **C23 works, and it is cheap relative to its yield.** Two scans caught six prior-art items in one pass on a paper that had been through three DELTA runs and a FULL run. The two most valuable catches were the *newest* claims (A3 written 2026-09-11, B4 written 2026-09-12) — which argues for running C23 close to authorship, not only at certification.
- **The criterion's priority ordering held.** Both catches came from category (3), "derivations under a page", and A1/A4 from category (1), "a clean closed-form constant". Category (2) — a well-developed external field entered sideways — described every one of them.
- **`UNVERIFIABLE` earned its place in the vocabulary.** Eleven items landed there. Without it the honest options would have been to invent citations or to silently drop real prior art.
