# C23 inverse-citation scan — the contraction seam: is `j_0` as an E(3) zonal function prior art?

**Date:** 2026-09-12. **Target:** Paper 60 `sec:molecular` (`papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex`, near L832 and L925–942, L1093–1101).
**Trigger:** C23 trigger 2 (authorship) — three group-theoretic readings about to be added to an already-certified identification.
**Read-only.** Nothing under `papers/`, `geovac/`, `tests/`, `CLAUDE.md`, `CHANGELOG.md` was modified.

**Prior scans consumed first, not repeated:** `frames_riesz_overcompleteness_memo.md` (Q4.3 already states the translation-phase mechanism), `toeplitz_finite_section_memo.md` (Q4A already settles "the operator is prior art, the operator-theoretic reading is not"), `c23_paper60_analysis_memo.md` (B4 already returned PRIOR ART on the M-centre all-ones degeneracy), `c23_paper60_linalg_memo.md`, `sturmian_conditioning_prior_art_memo.md`. Kac–Murdock–Szego, Bottcher–Widom, Serra, Jaffard, flat-limit/RBF, Halmos, Lowdin/Slater–Koster, DLMF Bessel are settled and untouched here.

---

## Consolidated verdict

| Claim | Verdict | One-line justification |
|:--|:--|:--|
| **C1** zonal/contraction: `j_0` = E(3) zonal function, Gegenbauer `C^{(1)}` = SO(4) zonal function, related by the Inonu–Wigner contraction realized as Mehler–Heine | **PRIOR ART** | A named theorem with a dedicated lineage: Mehler–Heine (1861/1868) -> Clerc (1976) for symmetric spaces -> Dooley–Rice (1983) for rotation groups -> Diaz Martin–Pacharoni (2018) states *exactly* the `(SO(n+1),SO(n)) -> (M(n),SO(n))` contraction of zonal spherical functions. Our `n=3` case is the textbook special case. |
| **C2** trivial-character reading of the degeneracy + M-centre rank-one corollary | **PRIOR ART** (twice over), and **already written and already attributed in the paper** | The degeneracy content was ruled PRIOR ART on 2026-09-11 (flat limit, `c23_paper60_analysis_memo.md` B4) and the paper already cites Driscoll–Fornberg + Barthelme–Usevich at L1093–1101. The *quantitative* M-centre version is a sharper named theorem we do not yet cite: Batenkov–Demanet–Goldman–Yomdin's clustered-node bound, `lambda_min ~ SRF^{-2(l-1)}`. |
| **C3** V_0-independence | **ABSENT as a literature statement — and the claim as posed in the brief is WEAKER AND PARTLY WRONG relative to what the PM has already written** | No source says "the ill-conditioning of multicentre bases is a property of the translation, not of the metric". The live (uncommitted) Paper 60 text at L925–942 already states it *correctly*: measured across five weights, scoped to momentum-space multipliers, and concluding that the constant **is** preserved. The brief's wording ("the constant is not") contradicts the paper's own measurement. |

**State of play at scan time (checked, not assumed).** The working tree is *not* clean: `paper_60...tex` (+51 lines), `geovac/sturmian_sigma_law.py`, `debug/qa/numeric_registry.py`, `debug/p60_contraction_seam_probe.py` and `tests/test_paper60_contraction_window.py` are all modified/untracked. So **C3 has already landed** (uncommitted, and landed correctly); **C2 landed earlier** and is committed at L1093–1101 with its flat-limit attribution; **C1 has not landed** — no zonal/contraction/Mehler–Heine/E(3) string appears anywhere in the `.tex`. The C1 verdict below is therefore the only one that is still ahead of the edit, and it is the one that matters most.

**Headline:** C1 is genuine, verifiable, citable prior art with a named lineage. C2 is a relabel of something already found and already cited. C3 as posed would *undo* a correct, measured, carefully-fenced passage that is already in the paper.

---

## C1 — the zonal/contraction statement. **PRIOR ART.**

### What the literature actually says

The single most on-point source, **read**:

> R. Diaz Martin and I. Pacharoni, *Mehler–Heine formula: a generalization in the context of spherical functions*, arXiv:1807.03904 (submitted 10 July 2018).
> Abstract, verbatim: *"In this article, using the notion of group contraction, we obtain the spherical functions of the strong Gelfand pair (M(n),SO(n)) as an appropriate limit of spherical functions of the strong Gelfand pair (SO(n+1),SO(n)) and also of the strong Gelfand pair (SO_0(n,1),SO(n))."*

From their introduction (quoted via the ar5iv rendering of the arXiv source):

- the classical formula `lim_{N->inf} P_N(cos(z/N)) = J_0(z)`, attributed to *"Heine in 1861 and ... Mehler in 1868"*;
- at `alpha = beta = (n-2)/2`, *"on the left side we have the Gegenbauer polynomials that are orthogonal polynomials that correspond to the spherical functions associated with the Gelfand pair (SO(n+1),SO(n))"*;
- the right side is `J_{(n-2)/2}(z)/(z/2)^{(n-2)/2}`, which *"is a spherical function associated with the Gelfand pair (SO(n) x R^n, SO(n))"*.

For `n = 3` that reads: Gelfand pair `(SO(4),SO(3))` = the Fock sphere `S^3`, zonal functions `C_k^{(1)}`; Gelfand pair `(E(3),SO(3))`; limit `J_{1/2}(z)/(z/2)^{1/2} ∝ sin z / z = j_0(z)`. **That is C1, verbatim, as the `n=3` instance of a published general theorem.**

Their own reference list (read, as printed) gives the lineage:

| Source | Verified to | What it is |
|:--|:--|:--|
| E. Inonu, E. P. Wigner, *On the contraction of groups and their representations*, **Proc. Nat. Acad. Sci. USA 39** (1953) 510–524 | BIB-ONLY (as printed in [IW] of 1807.03904) | the contraction itself |
| J. L. Clerc, *Une formule asymptotique du type Mehler–Heine pour les zonales d'un espace riemannien symetrique*, **Studia Math. 57** (1976) 27–32 | BIB-ONLY (as printed in [Cl]) | **the named theorem**: a Mehler–Heine asymptotic for zonal spherical functions on a Riemannian symmetric space |
| A. H. Dooley, J. W. Rice, *Contractions of rotation groups and their representations*, **Math. Proc. Camb. Phil. Soc. 94** (1983) 509–517 | BIB-ONLY (as printed in [DR1]) | contraction at the level of representations of rotation groups |
| A. H. Dooley, J. W. Rice, *On contractions of semisimple Lie groups*, **Trans. Amer. Math. Soc. 289** (1985) 185–202 | BIB-ONLY (as printed in [DR2]); DOI not retrieved | the general contraction theorem |

Neighbouring lineage for the *physics* framing (sphere -> Euclidean contraction and separation of variables), all **BIB-ONLY**, retrieved from search indices and publisher landing pages, primaries not opened:

- A. A. Izmest'ev, G. S. Pogosyan, A. N. Sissakian, P. Winternitz, *Contractions of Lie algebras and separation of variables. The n-dimensional sphere*, **J. Math. Phys. 40** (1999) 1549.
- E. G. Kalnins, W. Miller Jr., G. S. Pogosyan, *Contractions of Lie algebras: applications to special functions and separation of variables*, **J. Phys. A 32** (1999) 4709–4732.
- G. S. Pogosyan, A. Yakhno, *Lie-algebra contractions and separation of variables. Three-dimensional sphere*, **Phys. At. Nucl. 72** (2009) 836–844, DOI 10.1134/S1063778809050123. Landing page 303-redirected to an SSO endpoint; abstract read only through a search index, which reports the contraction as *SO(4) -> E(3)* and *"six systems of coordinates on the three-dimensional sphere contract to nine systems"* on Euclidean space.

Also textbook and worth one clause: `j_0(pR)` as the angular average of `exp(i p.R)` is the `l=0` term of the **Rayleigh plane-wave expansion**, not something the paper needs to derive.

One source checked and **cleared as not relevant**: Subag, Baruch, Birman & Mann, *On the contraction of so(4) to iso(3)*, arXiv:1210.6023 — read; it does the contraction abstractly through matrix coefficients and contains no zonal-function / Gegenbauer / Bessel content.

### Independent re-verification (C23 hard rule)

Both legs re-derived here, not taken on report (`scratchpad/mh.py`, `anti.py`, `par.py`):

- `C_k^{(1)}(cos chi) = sin((k+1)chi)/sin chi` to `<= 5e-15` over `k in {1,3,7,12}`, `chi in {0.3,1.1,2.7}`.
- `n^{-1} U_{n-1}(cos(z/n)) -> j_0(z)`: errors `8.2e-3, 8.0e-5, 8.0e-7, 8.1e-9` at `n = 10,100,1000,10000` (`z=5`), successive ratios `0.0097, 0.0100, 0.0101` — the limit holds and its own remainder is `O(n^-2)`.
- Index bookkeeping: `S^3` gives Jacobi `alpha=beta=1/2`, Gegenbauer `lambda=1`, Bessel order `1/2`. Consistent.

### Two technical cautions the paper must carry (both measured here)

1. **The contraction applies at the ANTIPODAL end, and the paper's chart puts the degeneracy there.** In Paper 60's convention `p = k cot(chi/2)`, so `chi -> 0` is `p -> infinity` (the chirp end) and `chi -> pi` is `p = 0` (the degenerate end). Mehler–Heine is classically stated at `cos chi -> +1`, i.e. `chi -> 0`. The limit *does* also hold at `chi = pi - z/n` — verified — but it is a different endpoint of the Jacobi asymptotic and must be written as such, or a referee reads the citation as pointing at the wrong corner.
2. **There is a parity factor.** `U_{n-1}(cos(pi-theta)) = (-1)^{n+1} sin(n theta)/sin theta`. Measured at `z=5`: `n=10 -> +0.19178565`, `n=11 -> -0.19178565`, `n=100 -> +`, `n=101 -> -`. So `n^{-1}U_{n-1}(cos(pi - z/n))` converges to `(-1)^{n+1} j_0(z)` — it **does not converge**, it converges along parities. This is the Hankel half of the Toeplitz-minus-Hankel structure showing up, and a bare statement "the S^3 zonal functions contract to `j_0` at the degeneracy" is false as written without it.

### The sentence to use (attribution wording, not a rewrite)

> The identification of the two objects is classical and we use rather than derive it: `sin(n chi)/(n sin chi)` is the zonal spherical function of the Gelfand pair `(SO(4),SO(3))` and `j_0(pR)` that of `(E(3),SO(3))`, and the first degenerates into the second under the Inonu–Wigner contraction `SO(4) -> E(3)` [Inonu–Wigner 1953], realized analytically as the Mehler–Heine limit — a theorem stated for a general Riemannian symmetric space by Clerc [Clerc 1976], at the level of representations of the rotation groups by Dooley and Rice [Dooley–Rice 1983], and in exactly the Gelfand-pair form used here, `(SO(n+1),SO(n)) -> (M(n),SO(n))`, by Diaz Martin and Pacharoni [DM–P 2018]. What is ours is only the placement: that this contraction is where the finite-section asymptotics of the molecular metric localize, at the antipodal corner `chi -> pi` of the Fock chart and along the scaling `pi - chi ~ 1/n` (with the parity factor `(-1)^{n+1}`).

---

## C2 — the trivial-character reading. **PRIOR ART, and largely already in the paper.**

### It is already written, and already attributed

Paper 60 L1093–1101 already says: *"a kernel matrix whose entries all tend to a common value degenerates to a rank-one all-ones matrix with an (M-1)-fold null space, the flat limit of the radial-basis-function literature [driscoll_fornberg2002, barthelme_usevich2021]. What is ours is that the block symbols realize it at a symbol point --- chi=pi, i.e. p=0 --- rather than in a shape-parameter limit."* That is the correct C23 disposition and it was applied yesterday (B4).

The addition on offer — "`p=0` is the **trivial character** of the translation group" — is a **relabel of an already-attributed statement**, not a new claim. It is worth making, because it names the mechanism in one word, but it must not re-open the novelty ledger. `frames_riesz_overcompleteness_memo.md` Q4.3 already stated the same mechanism in coordinate-free form on 2026-09-11.

### The corollary is a sharper named theorem that we do NOT yet cite

The claimed corollary — *"for M centres the symbol degenerates to the rank-one all-ones matrix, so the null space has fixed dimension M-1 and a geometry-independent direction"* — is the qualitative shadow of a quantitative theorem in the super-resolution literature:

> D. Batenkov, L. Demanet, G. Goldman, Y. Yomdin, *Conditioning of partial nonuniform Fourier matrices with clustered nodes*, arXiv:1809.00658; **SIAM J. Matrix Anal. Appl. 41** (2020) 199–220 `[theorem READ via ar5iv; volume/pages from secondary sources plus the SIAM manuscript number 18M1212197 appearing in the MIT DSpace filename — NOT confirmed on the publisher page]`.
> Theorem 3.2, as read: `sigma_min(V_N(x,Omega)) >= C_main(s) (Delta Omega)^{l-1}`, equivalently `lambda_min(G) >= C_main^2 (Delta Omega)^{2(l-1)}`, with `l` the maximal number of nodes clustered within one inverse-bandwidth. In SRF terms `lambda_min ~ SRF^{-2(l-1)}`.

And the bridge between that and the flat-limit literature we already cite is itself published:

> N. Diab, D. Batenkov, *Spectral properties of infinitely smooth kernel matrices in the single cluster limit, with applications to multivariate super-resolution*, arXiv:2407.10600 (July 2024, rev. Jan 2026), no journal ref listed `[abstract READ]`: *"We study the spectral properties of infinitely smooth multivariate kernel matrices when the nodes form a single cluster. We show that the **geometry of the nodes plays an important role** in the scaling of the eigenvalues of these kernel matrices."*

### Independent re-verification — and it contradicts the naive corollary

I built the M-centre symbol matrix `[j_0(p d_ij)]` directly and measured the eigenvalue orders in `p` (`scratchpad/mcent.py`):

| configuration | orders of the `M-1` small eigenvalues in `p` | `lambda_min` scaling |
|:--|:--|:--|
| `M=2`, `d=1` | 2.00 | `SRF^-2` |
| `M=3` collinear equispaced | 2.00, 4.00 | `SRF^-4` |
| `M=3` collinear uneven (0, 1, 3.7) | 2.00, 4.00 | `SRF^-4` |
| `M=3` equilateral | 2.00, **2.00** | `SRF^-2` |
| `M=4` collinear | 2.00, 4.00, 6.07 | `SRF^-6` |
| `M=4` tetrahedron | 2.00, **2.00, 2.00** | `SRF^-2` |

The `(M-1)`-dimensional null space at the exact point is confirmed in every case (each null vector's entries sum to `<= 3.4e-16`, i.e. the null space is exactly `1^perp`, and that part *is* geometry-independent). **But the orders are not.** A linear molecule pays `SRF^{-2(M-1)}`; a non-degenerate `d`-dimensional arrangement pays only `SRF^{-2}` with multiplicity. This is exactly Batenkov's `l-1` exponent and exactly Diab–Batenkov's "geometry plays an important role", and it is the flat-limit order structure (Schaback's eigenvalue orders, set by polynomial degrees available in `d` dimensions) seen in the Fourier variable.

**This is the thing the PM is about to over-claim.** "Geometry-independent direction" is true and trivial (`ker J_M = 1^perp`). "Fixed rank-`(M-1)` rotation removes it" is true for `M=2`, is *measured* for water's `A_1` block (a 2-dimensional block, so effectively the same case), and is **not established for a linear polyatomic**, where the residual after removing the `(M-1)` space has eigenvalue orders `2,4,...,2(M-1)` rather than a single order. Paper 60's L1093 sentence "spanned by a fixed vector set that does not move with geometry or with basis size" is defensible; the sentence at L856 "the fixed, geometry-independent rank-`(M-1)` rotation ... removes it" is the one that needs a scope clause, because *what is left after the rotation* does move with geometry.

### The sentence to use

> The `M`-centre degeneracy is the trivial character `p = 0` of the translation group, at which every block symbol equals `j_0(0) = 1` and the matrix symbol is the rank-one all-ones matrix `J_M` with null space `1^perp`. Both the degeneracy and its `(M-1)`-fold null space are known — as the flat limit of kernel matrices [Driscoll–Fornberg; Barthelme–Usevich], and quantitatively as the clustered-node conditioning bound `lambda_min ~ SRF^{-2(l-1)}` of Batenkov, Demanet, Goldman and Yomdin. We claim only the placement (a symbol point rather than a shape-parameter limit) and we adopt their warning: the *null space* is geometry-independent, the *rates along which it opens* are not — we measure orders `(2,4)` for three collinear centres against `(2,2)` for three centres in general position.

---

## C3 — the V_0-independence corollary. **ABSENT in the literature; and the brief's wording is weaker than what the PM already wrote.**

### The literature search

No source states "the ill-conditioning of multicentre bases is a property of the translation, not of the metric". Searched: shift-invariant/systems-of-translates (Ron–Shen lineage — already settled by the frames scan, not re-litigated), Schoenberg matrices, kernel conditioning, and the QC basis-set linear-dependence canon.

The nearest thing found, and it is genuinely adjacent:

> L. Golinskii, M. Malamud, L. Oridoroga, *Schoenberg matrices of radial positive definite functions and Riesz sequences in L^2(R^n)*, arXiv:1403.2234 `[READ via ar5iv; journal reference NOT verified]`. Their Theorem 1.9, as read: for a radial `f in L^2(R^n)` with non-vanishing Fourier transform, `{f(. - x_j)}` is a Riesz sequence **iff the node set X is separated** — i.e. the Riesz property is decided by the geometry of the translations and not by which radial `f` generates them. That is the closest published statement in the direction of C3, and it is a *different* theorem (infinite node set, one generator per node, separation as the criterion; ours is two nodes with an infinite per-node basis and truncation order as the asymptotic parameter). Cite as an analogue, never as the source.

QC side: nothing. Aissing–Monkhorst's "linear dependence is intrinsic to three-dimensionally extended systems" is a dimension-counting argument about density of functions, not a translation-vs-metric statement.

The mathematical content of C3 — that the pencil `T_n(W sigma) v = lambda T_n(W) v` has its spectrum governed by the *ratio* symbol `sigma = (W sigma)/W` — is standard preconditioned-Toeplitz theory (the eigenvalues of `T_n(f)^{-1}T_n(g)` lie in the essential range of `g/f`, and the one-line proof is monotonicity of `T_n` in the symbol: `m f <= g <= M f` pointwise implies `m T_n(f) <= T_n(g) <= M T_n(f)`). I could **not** open a primary stating it, so it is flagged in the weakest-link section rather than cited.

### The wording gap — this is the most important finding in the scan

**Paper 60 already contains C3, at L925–942 of the live working tree (uncommitted), stated better than the version posed to me in the brief.** Verbatim from the paper: the collapse `(1-sigma_max)(n/kR)^2` reads `0.4072, 0.4123, 0.4107, 0.4164` for `W = 1, 1 + (4/5)cos chi, 2 + sin chi, e^{-chi}` against `pi^2/24 = 0.4112`; the control `W = 1 + cos chi`, which *vanishes* at the degeneracy, moves the constant to `0.828`; and the paper then explicitly fences the overstatement: *"This is not independence of the weighting potential V_0: a position-space-local V_0 acts on momentum space by convolution rather than multiplication and therefore leaves this class entirely. What is shown is narrower and still worth having --- the law is carried by the translation phase, not by any momentum-space weight multiplying it."*

The C3 *as posed in the brief* says two things the paper text has already, correctly, refused:

1. **"for any radial local weighting potential `W(p)`"** — a *local* `V_0` is a multiplication operator in **position** space, hence a **convolution** in momentum space, hence not of the form `W(p) x`. The paper flags exactly this conflation as "the natural overstatement [that] is close by". C3 as posed commits it.
2. **"the CONSTANT is not [independent]"** — this is backwards for the class in which the argument works. I re-measured it independently (`scratchpad/c3.py`, `c3b.py`), five weights, `kR=2`, `n = 20..320`, generalized symmetric eigenproblem in the sine basis:

| weight `W(chi)` | exponent (last window) | `(1-lambda_max) n^2` at `n=320` |
|:--|:--:|:--:|
| `1` (control) | 1.993 | 1.637 |
| `1 + 0.5 cos chi` | 1.992 | 1.636 |
| `2 - cos chi` | 1.982 | 1.624 |
| `(1 + p^2/k^2)^{-1}` | 1.977 | 1.619 |
| `sin^2(chi/2) + 0.3` | 1.981 | 1.623 |
| **`cos^2(chi/2)` — vanishes at `p=0`** | 1.989 | **3.338** |
| **`cos^4(chi/2)` — double zero at `p=0`** | 1.982 | **5.483** |
| `1/(cos^2(chi/2)+0.01)` — peaked at `p=inf` | 1.953 | 1.592 |

(reference `pi^2 (kR)^2/24 = 1.64493`.)

So: the exponent is `2` for **every** weight tested, including ones degenerate at the trivial character; and the constant is invariant too — to within the pre-asymptotic `O(1/n)` residue — for every weight that is continuous and **non-zero at `chi = pi`**. It moves only when `V_0` degenerates *exactly at the trivial character*, and then by a factor `2.04` for a simple zero and `3.35` for a double zero. My `cos^2` factor `2.04` independently reproduces the paper's own `0.828/0.4112 = 2.01` control.

**The honest statement is stronger than the one on offer, and localizes the freedom precisely:** both the exponent and the constant are carried by the translation character alone; the only escape inside the class is a weight that vanishes at `p = 0`, which is the same escape route the frames scan already identified (bound the Fock-momentum content away from the trivial character).

### The sentence to use

Do not add a second C3 and do not restate it in the brief's wording. Keep the L925–942 paragraph as written (its guards in `tests/test_paper60_contraction_window.py` — `test_smooth_weights_preserve_the_collapse_constant`, `test_vanishing_weight_control_moves_the_constant`, `test_weight_independence_holds_at_the_exponent_too` — already fence it correctly, including the vanishing-weight control). If anything, tighten it to:

> Within the class of momentum-space multipliers the law is carried entirely by the translation character: exponent *and* constant are the same for every weight continuous and non-vanishing at `p = 0` (measured across eight weights), and the constant moves only for a weight that degenerates at the trivial character itself (factor `2.04` for a simple zero, `3.35` for a double zero). A position-space-local `V_0` is a convolution in momentum space and lies outside this class; whether any such `V_0` preserves the constant remains open.

---

## Weakest link

What I could not verify, and what would close it:

1. **Shibuya & Wulfman (1965), the primary.** `Proc. R. Soc. A 286, 376–389` — royalsocietypublishing.org returned **HTTP 403**; no preprint exists for a 1965 paper; Semantic Scholar carries metadata only. **UNVERIFIABLE.** So the question the task posed — do SW themselves state a translation operator on the Fock sphere, or any group-theoretic reading of their integrals? — is *unanswered*, exactly as the 2026-09-11 scan left it. Secondary summaries consistently say only that they "derived an expansion of a plane wave involving the four-dimensional spherical harmonics". **Closing move:** an institutional PDF of Proc. R. Soc. A 286, or Avery's book chapter *Electronic Structure Theory in Momentum Space* (Springer, `10.1007/978-94-009-2329-4_10`, also paywalled). Until then the paper must not assert what SW do or do not say.
2. **Avery's canon, priority 4.** Every route tried was paywalled or 403 (ResearchGate "Request PDF" for *Many-center Coulomb Sturmians and Shibuya–Wulfman integrals* and for Red & Weatherford's *Derivation of a general formula for the Shibuya–Wulfman matrix*). **UNVERIFIABLE**, unchanged from the prior scan's verdict.
3. **The Toeplitz-pencil / ratio-symbol theorem** that is C3's mathematical backbone. I verified it numerically and it has a one-line proof, but I could not open a primary (Serra Capizzano's BIT/LAA papers, Chan–Ng survey — all paywalled; search returned only GLT-era secondary material). Do **not** add a Serra citation on this pass; either derive the monotonicity line in-paper or open the primary first.
4. **Clerc (1976) and Dooley–Rice (1983/1985)** are **BIB-ONLY** — read as printed in Diaz Martin–Pacharoni's reference list, primaries not opened. Diaz Martin–Pacharoni itself is `[READ]` (ar5iv rendering of the arXiv source, direct quotes above) but its **published** venue is unknown; cite it as arXiv:1807.03904 unless someone confirms a journal version.
5. **Batenkov et al. journal reference.** Theorem 3.2 is `[READ]`; the `SIAM J. Matrix Anal. Appl. 41 (2020) 199–220` coordinates are secondary-source plus the SIAM manuscript id in a DSpace filename. Confirm before printing volume/pages.
6. **Not attempted, and it is the real remaining gap:** whether the *momentum-space quantum-chemistry* literature (Aquilanti–Cavalli–Coletti; Avery) has ever used the SO(4) -> E(3) contraction to discuss the `p -> 0` end of a Sturmian basis. The generic physics statement ("the hydrogen atom's SO(4) contracts to E(3) at the ionization threshold") is well known; whether anyone joined it to multicentre overlap conditioning is unverified. This is where a surviving novelty claim for the *seam* would live or die.

---

## What survives as ours

The corpus's standing position — that the **identification** of the two-centre Shibuya–Wulfman metric with a Toeplitz-minus-Hankel finite section of symbol `j_0(kR cot(chi/2))` is ours, even though the Toeplitz asymptotics are Kac–Murdock–Szego — **is not changed by C1–C3.** Nothing found here touches the identification. What C1–C3 change is the *size* of what may be added on top of it, and the answer is: less than proposed, in all three cases.

- **C1 adds nothing to the novelty ledger and should be written purely as attribution.** The contraction, the Mehler–Heine realization, and even the exact Gelfand-pair form are a named, published theorem chain. The only sentence that is ours is the *placement* — that this contraction is where the finite-section asymptotics localize, at the antipodal corner and along `pi - chi ~ 1/n` — and even that is a reading of our own already-published law rather than a new result. It is worth one clause, not a subsection, and it must carry the endpoint and parity cautions above or it is wrong as written.
- **C2 is a relabel of a claim already found, already re-tiered, and already cited in the paper.** Adding "trivial character" is good expository physics; presenting it as a new group-theoretic result would be the second rediscovery of the same degeneracy in two days.
- **C3 should not be added again** — the live tree already has it, measured, better scoped, guarded, and with the constant's behaviour the right way round. The only action is to not let the brief's weaker wording overwrite it.

**One thing NOT re-litigated, and it should be:** the landed paragraph before C3 tags the `pi^2` as truncation-side M2 and the `sqrt(pi)` / `pi/4` as continuum-side M2. That is a Paper-18 / Paper-34 classification claim, it is new, it is `[SYMBOLIC + MEASURED]`, and it was outside this scan's three-claim remit. It carries its own C23 trigger-2 obligation (and the `feedback_tag_transcendentals` memory rule). Flagging, not scanning.

The net effect of this scan is therefore subtractive on claims and additive on citations, which is the ordinary C23 outcome. The one genuinely new *finding* is the M-centre order table: the `(M-1)` null space is geometry-independent but the orders along which it opens are not (`2,4` collinear vs `2,2` in general position), which is both a correction to the proposed corollary and a concrete, cheap prediction the composed-polyatomic arc can test.

---

## Files

- This memo.
- Scratchpad drivers, not in the repo: `mh.py` (Gegenbauer/Chebyshev identity + Mehler–Heine convergence), `anti.py` (which sphere end the contraction sits at), `par.py` (the `(-1)^{n+1}` parity factor), `mcent.py` (M-centre eigenvalue orders), `c3.py` / `c3b.py` (pencil ratio-symbol, eight weights).
- Nothing under `papers/`, `geovac/`, `tests/`, `CLAUDE.md` or `CHANGELOG.md` was modified.
