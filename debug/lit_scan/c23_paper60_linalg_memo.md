# C23 inverse-citation scan — Paper 60, linear-algebra / operator-theory claims

**Date:** 2026-09-12 · **Type:** LITERATURE SCAN (read-only; no paper/code/doc edits) ·
**Criterion:** `docs/qa/criteria.md` C23 (inverse citation — is an UNcited result already named elsewhere?)
**Scope:** claims A1–A4 of the Paper 60 conditioning / composition-wall section
(`papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex`, L786–L892).

**Verification discipline.** Every citation below is either (a) read in primary full text in this
session — marked **[READ]** — or (b) bibliographic data confirmed from a primary source that cites it —
marked **[CITED-BY-PRIMARY]** — or (c) explicitly marked UNVERIFIED. Nothing is asserted from a
search-engine summary. The mathematics of A1/A2/A4 was independently re-derived numerically here
(`scratchpad/c23_verify.py`, see Sec. 5) per C23's "re-verify the identification yourself" rule.

**Search budget note.** The session's shared WebSearch quota (200) was exhausted partway through;
the second half of the scan ran on direct WebFetch/curl of primary sources and on citation-chaining
out of papers already obtained. Where that left a leg unreachable it is recorded as UNVERIFIABLE
rather than as a clearance (C23 hard rule).

---

## Verdict table

| # | Claim | Verdict | Load-bearing citation |
|:--|:------|:--------|:----------------------|
| A1 | `spec[[I,C],[C^T,I]] = {1±sigma_k}`, `cond = (1+s_max)/(1-s_max)` | **PRIOR ART** (standard, assembled from three named results; no single source states the composite) | Jordan–Wielandt; Jordan 1875; two-projections canonical form |
| A2 | `norm[P,Q] = max_k sigma_k sqrt(1-sigma_k^2) <= 1/2` | **PRIOR ART** — and the *existing* Halmos attribution is **DEFECTIVE** | Loring, Ann. Funct. Anal. **5** (2014) 176–187 |
| A3 | block-diagonal congruence cannot create block structure ("Proposition D") | **ABSENT** as the stated contrapositive / **PRIOR ART** for its content | Slater–Koster, Phys. Rev. **94**, 1498 (1954), Appendix |
| A4 | `2/(1+min_x j_0(x)) = 2.555041...` | **ABSENT** (constant); mechanism is standard Toeplitz theory | — (mechanism UNVERIFIED, see Sec. 4) |

---

## 1. A1 — Gram spectrum `{1 ± sigma_k}` and `cond = (1+s_max)/(1-s_max)`

**Paper 60 text (L786–L794), tiered `[SYMBOLIC + MEASURED]`:** "The conditioning law is in fact
*derived* ... Because the intra-center block is exactly the identity, `S = [[I,C],[C^T,I]]` has spectrum
`{1 +- sigma_k}` ... so `cond(S) = (1+s_max)/(1-s_max)` *exactly*." Paper 60 *does* cite `amos_hall1961`
and `king1967` — but only for the sigma_k / principal-angle identification. The spectral identity and
the condition-number formula are presented as the paper's own derivation.

### Verdict: **PRIOR ART.** Not new; standard in three named pieces.

**(i) The spectrum.** `[[0,C],[C*,0]]` is the **Jordan–Wielandt matrix** (a.k.a. the Hermitian
dilation), whose eigenvalues are `+-sigma_k(C)` together with zeros; `S = I +` that matrix, so
`spec(S) = {1 +- sigma_k}` union `{1}`. The naming and its standard attributions (Wielandt 1955;
Stewart & Sun, *Matrix Perturbation Theory*, Academic Press 1990; Horn & Johnson, *Matrix Analysis*,
CUP 2nd ed. 2012) are stated in **[READ]** "Higher-Order Singular-Value Derivatives of Real Rectangular
Matrices", arXiv:2506.03764, Sec. 1.2 and "Theorem 1.2 (Spectrum of Jordan–Wielandt Embedding)".
**I did not open Stewart–Sun or Horn & Johnson myself** — a search summary offered "Theorem 7.3.3 /
I.4.2" and that number is **UNVERIFIED**; cite the books without a theorem number, or open them first.

**(ii) sigma_k = cosines of principal angles.** Origin: **C. Jordan**, "Essai sur la geometrie a n
dimensions", *Bull. Soc. Math. France* **3**, 103–174 (1875) — **[CITED-BY-PRIMARY]** as the origin in
both Kim & Kim arXiv:1811.10518 Sec. 1 (read) and Loring arXiv:1306.1923 ref. [8] (read). Computational
form: A. Bjorck and G. H. Golub, "Numerical methods for computing angles between linear subspaces",
*Math. Comp.* **27**, 579–594 (1973) — the AMS PDF returned **HTTP 403**, so its *content* is
**UNVERIFIABLE** here; do not add it as backing for the Gram-spectrum statement specifically.

**(iii) The composite, in two-projections language.** The Gram operator of the union of two orthonormal
sets is unitarily `P + Q`; decomposing into Jordan planes gives, per plane, `P|J = [[1,0],[0,0]]` and
`Q|J = [[c^2, cs],[cs, s^2]]`, whence `(P+Q)|J` has eigenvalues `1 +- c`. Read directly in **[READ]**
J. Kim and Y. Kim, "Jordan plane and numerical range of operators involving two projections",
arXiv:1811.10518 [math.FA] (2018) — the displayed 2x2 representations in Sec. 3, and their explicit
`norm(P_M + P_N) = 1 + cos theta_1` obtained by combining their Prop. 2.3 (`norm(P+Q) = 1 + norm(PQ)`)
and Prop. 2.4 (`norm(P_M P_N) = cos theta_1`).

**(iv) The condition-number form is a *named constant* in numerical analysis.** `D^(-1/2) A D^(-1/2) =
[[I,C],[C*,I]]` is exactly the block-Jacobi-preconditioned form of a 2x2 block SPD matrix, and
`s_max = norm(C)` is exactly the **CBS constant** (strengthened Cauchy–Bunyakowski–Schwarz), which is
*defined* as the cosine of the abstract angle between the two subspaces:
`gamma := cos(V1,V2) := sup_{u in V1, v in V2} A(u,v)/sqrt(A(u,u)A(v,v))` — **[READ]** verbatim in
J. Kraus, M. Lymbery, S. Margenov, "Robust algebraic multilevel preconditioning in H(curl) and H(div)",
arXiv:1301.3269, Sec. 3 (the display introducing their Eq. (3.19)). That paper attributes the surrounding
theory to V. Eijkhout and P. S. Vassilevski, "The role of the strengthened Cauchy–Buniakowskii–Schwarz
inequality in multilevel methods", *SIAM Review* **33**, 405–419 (1991) (their ref. [15]).
**The specific equality `kappa = (1+gamma)/(1-gamma)` I could not open a source for** —
Eijkhout–Vassilevski, Axelsson's *Iterative Solution Methods* (CUP 1994), and Vassilevski's *Multilevel
Block Factorization Preconditioners* (Springer 2008) are all paywalled/offline here. That leg is
**UNVERIFIABLE**; do not cite it until someone opens one of them.

**Recommended remediation.** Re-tier A1 from "derived" to "standard, recalled here", cite
Jordan 1875 plus the Jordan–Wielandt fact (books, no theorem number until checked), and keep
Amos–Hall/King where they are. GeoVac's actual contribution in this passage survives intact and is
*sharper* once A1 is recalled rather than derived: it is the identification of `C` as the finite
section of a multiplication operator with symbol `j_0(kR cot(chi/2))`, and the resulting
`1 - s_max ~ (kR)^2 pi^2 / 24 n^2` law. That is the part nobody else has.

---

## 2. A2 — commutator norm `norm[P,Q] = max_k sigma_k sqrt(1-sigma_k^2) <= 1/2`

**Paper 60 text (L847–L856):** the norm identity is "a one-line consequence of the two-subspaces
canonical form~\cite{halmos1969,bottcher_spitkovsky2010}, which resolves the pair into 2x2 blocks at the
principal angles theta_k, in which the commutator has norm sin theta_k cos theta_k."

**Credit where due:** the paper already attributes only the *canonical form* to Halmos and derives the
norm itself — i.e. the 2026-09-06 scan's recommendation was applied. What remains is (a) the norm
identity is still presented as GeoVac's own one-line derivation when a citable primary source exists,
and (b) Bottcher–Spitkovsky is cited for content nobody in this corpus has read.

### Verdict: **PRIOR ART for the identity; Halmos alone cannot back it (confirmed independently).**

**Halmos — re-verified against the primary PDF this session.** P. R. Halmos, "Two subspaces",
*Trans. Amer. Math. Soc.* **144**, 381–389 (1969). **[READ]** Theorem 2, verbatim:

> "**Theorem 2.** If M and N are subspaces in generic position in a Hilbert space H, with respective
> projections P and Q, then there exists a Hilbert space K, and there exist positive contractions S and
> C on K, with S^2 + C^2 = 1 and ker S = ker C = 0, such that P and Q are unitarily equivalent to
> `[[1,0],[0,0]]` and `[[C^2, CS],[CS, S^2]]` respectively."

That is the **canonical form only**. Halmos states **no** commutator norm anywhere in the paper. The
earlier scan (`debug/lit_scan/toeplitz_finite_section_memo.md` Sec. D) was right, and this scan
independently confirms it.

**Loring is the correct primary citation for the norm identity.** T. A. Loring, "Principal angles and
approximation for quaternionic projections", *Ann. Funct. Anal.* **5** (2) 176–187 (2014),
DOI 10.15352/afa/1396833512; preprint arXiv:1306.1923. **[READ]** in full text. His Theorem 1.1 gives
the block reduction with `P = [[1,0],[0,0]]`, `Q = [[cos^2 t, cos t sin t],[cos t sin t, sin^2 t]]`;
in Sec. 3, in the argument establishing Theorem 3.2, the text states verbatim:

> "For all theta we find `norm(PQ - QP) = (1/2) sin(2 theta)`."

and his **Corollary 3.1** imposes the relation `norm(pq - qp) <= epsilon` only for
`0 <= epsilon < 1/2`, with `C_epsilon = (1/2) arcsin(2 epsilon)` — i.e. the `1/2` ceiling is built into
the statement. Since `(1/2) sin 2t = sin t cos t = sigma sqrt(1 - sigma^2)`, this is exactly GeoVac's
per-block identity; the `max_k` over blocks then follows from the direct-sum structure of Halmos Thm 2 /
Loring Thm 1.1.

**One honest caveat to carry into the paper:** Loring's sentence is a displayed line inside a proof, not
a numbered theorem. Cite it as "Loring (2014), Sec. 3 (proof of Thm. 3.2)"; do not invent a theorem number.

**Bottcher–Spitkovsky: bibliographic data confirmed, content still UNVERIFIABLE.** A. Bottcher and
I. M. Spitkovsky, "A gentle guide to the basics of two projections theory", *Linear Algebra Appl.*
**432** (6), 1412–1459 (2010), DOI 10.1016/j.laa.2009.11.002 — the volume/issue/pages are confirmed
**[CITED-BY-PRIMARY]** by two independent reference lists read in full this session (Kim & Kim
arXiv:1811.10518 ref. [2]; arXiv:1212.1996 ref. [BS10]). Whether it *states*
`norm[P,Q] = max sigma sqrt(1-sigma^2)` remains unverified — paywalled, no preprint found.
**Do not attribute the identity to it.** Either open the PDF or drop it from that particular cite pair
and let Halmos (canonical form) + Loring (norm) carry it.

**Recommended remediation.** Keep `\cite{halmos1969}` where it is (canonical form — correct as written);
add Loring 2014 for the norm identity and change "a one-line consequence of" to "recorded for the 2x2
blocks by Loring, and summed over blocks by the canonical form of Halmos"; Bottcher–Spitkovsky either
verified by opening the PDF or demoted to a general "see also" for two-projections theory.

---

## 3. A3 — "Proposition D": block-diagonal congruence cannot create block structure

**Paper 60 text (L858–L874), tiered `[SYMBOLIC]`:** "If `X` is invertible and block diagonal with
respect to `H = (+)_l H_l`, and `X^dagger S X` is block diagonal, then
`S = X^(-dagger) (X^dagger S X) X^(-1)` is block diagonal too, being a product of block-diagonal factors.
Contrapositively: if `S` is not l-block diagonal, *no* block-diagonal congruence — Lowdin, canonical, or
Cholesky — orthogonalizes it."

### Verdict: **ABSENT as stated / PRIOR ART for the content it packages.** Split, deliberately.

**ABSENT — the exact contrapositive.** I found no source stating this proposition. Searched: the
quantum-chemistry orthogonalization/symmetry literature (Lowdin symmetric vs. canonical
orthogonalization; basis-set linear dependence); numerical-linear-algebra structured-congruence and
structure-preserving-canonical-form material; two-projections and frame theory. Nothing states
"block-diagonal congruence cannot create block-diagonality the metric lacks" as a proposition.

**But the mathematical content is not new, and GeoVac cites nobody for it.** Two observations, and the
second is the load-bearing one for remediation:

**(a) It is one line from a textbook fact about algebras.** Block-diagonal matrices w.r.t. a fixed
grading are precisely the commutant `{P_l}'` of the grading projections. A commutant is a von Neumann
algebra, hence a unital *-algebra closed under adjoints *and* under inverses of its invertible elements.
`S = X^(-dagger) (X^dagger S X) X^(-1)` is then a product of three elements of that algebra. GeoVac's
proof *is* this argument; it just does not name the structure. This is exactly C23's "anything whose
derivation took under a page" bucket.

**(b) The chemistry-side statement is the SLATER–KOSTER THEOREM, and Paper 60 does not cite it.**
J. C. Slater and G. F. Koster, "Simplified LCAO method for the periodic potential problem",
*Phys. Rev.* **94**, 1498 (1954) — Appendix. **[CITED-BY-PRIMARY, quoted in full text]**:
T. A. Rokob, A. Szabados, P. R. Surjan, "A note on the symmetry properties of Lowdin's orthogonalization
schemes" (Zahradnik-80 dedication; read in full at
`https://coulson.chem.elte.hu/surjan/Zahradnik80.pdf`) states:

> "In the appendix of their seminal paper, Slater and Koster[5] proved the following theorem: Let T be a
> symmetry operator of the system. Then, the transformation properties of symmetrically orthogonalized
> vectors (1) are the same as those of the original nonorthogonal set, i.e. the matrices representing T
> in both sets are identical. This orthogonalization thus preserves the symmetry of the basis."

their ref. [5] being `J. C. Slater and G. F. Koster. Phys. Rev., 94:1498, 1954`.

That note is **more useful to GeoVac than Slater–Koster itself**, because it isolates the hypothesis
that is exactly GeoVac's mechanism. Rokob–Szabados–Surjan show, reading the operator relation
`T^dagger T = I` in an overlapping basis, that `t = S^(-1) T` is unitary — and hence that
Lowdin/canonical orthogonalization preserves symmetry — **if and only if `[T, S] = 0`**; they exhibit
the failure for redundant Cartesian 6d/10f sets, where `[T,S] != 0`, the S-eigenvalues lose their
degeneracies, and canonically orthogonalized AOs no longer span representations. This is the same
"the orthogonalization can only inherit the structure the metric already has" content, stated in the
quantum-chemistry register, in 1954 and 2010. They also give the `s`-only case explicitly (t a
permutation matrix, hence `[T,S] = 0`), which is the homonuclear-diatomic case Paper 60 works in.

**Also-ran, flagged UNVERIFIED.** The tight-binding/LMTO literature asserts the l-specific form —
"Lowdin orthogonalization is a symmetry transformation: angular symmetry is preserved, orbitals with
angular momentum l are not mixed with l' != l" — I saw this in search summaries only and did **not**
read a primary source for that wording. Do not cite it on this memo's authority.

**Recommended remediation.** Keep the proposition — it is correctly proved and it does real work in the
argument — but stop presenting it as unattributed new mathematics. Demote "Proposition D" to a Remark,
state it as the contrapositive of the Slater–Koster theorem, cite Slater & Koster 1954 and
Rokob–Szabados–Surjan for the `[T,S] = 0` hypothesis, and note in one clause that the linear algebra is
the inverse-closedness of the commutant. **The contribution that survives — and it is the real one — is
the application:** that the two-center Shibuya–Wulfman metric couples l while preserving `m`, so
`m`-selection survives orthogonalization and within-`m` l-selection provably cannot, *at every*
`cond(S) > 1`. That specialization is not in the prior art and should carry the emphasis.

**This is the C23 catch the criterion was built for** — a brand-new, one-line, uncited proposition whose
positive direction has had a name in this paper's own field for seventy years.

---

## 4. A4 — the constant `2/(1 + min_x j_0(x)) = 2.555041...`

**Paper 60 text (L890):** "`cond(I+C) -> 2/(1+min_x j_0(x)) = 2.555041...`, an exact constant of the
sinc symbol independent of R and of basis size".

### Verdict: **ABSENT.** The ingredients are tabulated; the combination is not.

**Ingredients, verified numerically here (Sec. 5):** `x_0 = 4.493409457909064` (first positive root of
`tan x = x`) and `min_x sinc(x) = -0.21723362821122166`. These are standard tabulated values of the sinc
function; they are *not* named constants in any source reached, and the search-summary claim that
MathWorld tabulates them did **not** survive an attempt to read the page (the fetched MathWorld
`SincFunction` page contains no such tabulation), so treat the "well-known constant" framing as
UNVERIFIED provenance even though the numbers themselves are elementary and were recomputed here to
16 digits.

**The combination `2/(1 + min j_0) = 2.5550407785526077...`:** no name, no source, no tabulation found.
OEIS could not be queried (HTTP 403 on the search endpoint), so "not in OEIS" is **UNVERIFIABLE**, not
established. Everything else reached was negative.

**But the *mechanism* is textbook and should be named.** `cond(I + T(a))` for a self-adjoint Toeplitz
operator is `(1 + max a)/(1 + min a)` because the spectrum of a self-adjoint Toeplitz operator is the
closed convex hull of the essential range of its symbol. That is a classical theorem (Hartman–Wintner);
**I did not open a source for it in this session and the corpus's own Toeplitz memo does not cite it
either** — flagged UNVERIFIED and owed. With `max j_0 = 1` and `min j_0` the sinc minimum, the constant
follows immediately. So A4 is best presented as "the Hartman–Wintner value for this symbol", i.e. a
worked instance rather than a discovery — which also explains its R- and n-independence in one clause
instead of a paragraph.

**Recommended remediation.** Keep the number (it is right, and it is the paper's own measurement), but
present it as the symbol's convex-hull value rather than as a standalone constant. Verify
Hartman–Wintner before citing it.

---

## 5. Independent re-verification (C23 hard rule)

`scratchpad/c23_verify.py`, run this session:

```
A1  max abs err (spectrum + cond), 300 random subspace pairs, dims 4-20 : 2.02e-10
A2  max abs err norm[P,Q] vs max_k sigma_k sqrt(1-sigma_k^2), 300 pairs : 1.31e-14
    ceiling  max_{s in [0,1]} s sqrt(1-s^2)                             : 0.5000000
A4  root of tan x = x : 4.493409457909064
    min j_0            : -0.21723362821122166
    2/(1 + min j_0)    :  2.5550407785526077
```

All four identities are mathematically correct. **Nothing here is a retraction** — per C23, finding
prior art re-tiers attribution, not truth.

---

## 6. Dependents to sweep if A1–A3 are re-attributed

Not swept here (read-only scan); listed so the remediation pass does not have to rediscover them.
`eq:sigma_law` and the composition-wall commutator are cited outside Paper 60 — check at minimum
`tests/test_paper60_sigma_law.py` (the backing test, named at both loci),
`memory/composition_wall_non_commuting_projections.md`,
`memory/avery_method_and_prior_art_gaps.md`, and Paper 58's decompactification-front passage, which
uses the same principal-angle object. `debug/lit_scan/projector_angle_bonding_memo.md` (2026-09-06) and
`debug/lit_scan/toeplitz_finite_section_memo.md` are the sibling scans; Sec. D of the latter reached the
same Halmos verdict independently and its recommendation is hereby confirmed rather than superseded.
