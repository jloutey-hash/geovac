# Lit scan: two primaries the corpus had twice failed to reach

Date: 2026-09-12. Read-only scan (nothing under `papers/`, `geovac/`, `tests/`,
`CLAUDE.md`, `CHANGELOG.md` was touched). Targets:

1. **Shibuya & Wulfman (1965)** — do *they* give a group-theoretic (translation-operator)
   reading of their own integrals?
2. **The Toeplitz pencil / ratio-symbol theorem** behind Paper 60's new
   weight-independence result (CHANGELOG v5.11.4) — is it a named theorem, and does the
   standard statement reach the *constant*?

| Target | Verdict |
|:--|:--|
| 1. Shibuya–Wulfman 1965 | **SECONDARY-QUOTED** (body UNREACHED; the paper's *own published abstract* reached verbatim, plus its complete reference list; the attribution question is settled by adjacent primaries) |
| 2. Toeplitz pencil / ratio symbol | **REACHED** — it is a named, published, **open-access** theorem, opened in full, with theorem number, hypotheses, and an *exact* τ-algebra identity that fits GeoVac's matrices better than the Toeplitz statement does |

---

## Method note: what changed this pass

The two previous scans were run through a summarising fetch tool, which 403s on
Royal Society / Wiley / Elsevier and cannot extract PDF text. This pass used
**direct network access from the shell** (`curl`) plus local `pdftotext`. That
opened three things the earlier passes could not reach:

- the **publisher's own abstract** of SW 1965 and its **complete reference list**,
  via the OpenAlex inverted index (raw JSON, reconstructed locally);
- **full text of two Toeplitz primaries** (one open-access journal article, one
  arXiv preprint), downloaded and read as text rather than summarised;
- **Semantic Scholar citation contexts** — verbatim sentences from citing papers.

Where a PDF was fetched and grepped locally, it is marked **[PRIMARY, full text
read]**. Where only the publisher-supplied abstract was obtained, it is marked
**[ABSTRACT, verbatim]**. Nothing in this memo is a search-engine paraphrase
unless explicitly labelled as such.

---
---

# TARGET 1 — Shibuya & Wulfman (1965)

**T. Shibuya and C. E. Wulfman, "Molecular orbitals in momentum space",
*Proc. R. Soc. Lond. A* **286** (1406), 376–389 (1965), DOI 10.1098/rspa.1965.0151.**
Authors' affiliation at the time: University of the Pacific. 79 citations (S2).

## Verdict: SECONDARY-QUOTED

The **body is still unreached** — the 1965 text is closed everywhere, and this is now
a well-characterised negative rather than a guess (routes and failure modes below).
But two things were obtained that the previous two scans did not have, and together
they answer the question the corpus was asking.

### 1a. What SW 1965 says about it — their own abstract, verbatim

Reconstructed from the publisher-deposited abstract (OpenAlex inverted index,
`abstract_inverted_index` for the DOI above; reconstruction done locally by
sorting the word positions). **[ABSTRACT, verbatim]**

> "The troublesome problem of developing cusps in ordinary molecular wave functions
> can be avoided by working with momentum-space wavefunctions for these have no
> cusps. The need for continuum wavefunctions can be eliminated if one works with a
> hydrogenic basis set in Fock's projective momentmn [*sic*, OCR] space. This basis
> set is the set of R4 spherical harmonics and as a consequence one may obtain,
> **solely by the ordinary angular momentum calculus**, algebraic expressions for all
> the integrals required in the solution of the momentum space Schrödinger equation.
> A number of these integrals and a number of **R4 transformation coefficients** are
> tabulated. The method is then applied to several simple united-atom and l.c.a.o.
> wavefunctions for H+2 and ground state energies and corrected wavefunctions are
> obtained. It is found in this numerical work that the method is most appropriate at
> internuclear distances somewhat less than the equilibrium distance. In Fock's
> representation both l.c.a.o. and unitedatom approximations become exact as the
> internuclear distance approaches zero. The united-atom expansion can be viewed as
> an eigenvalue equation for the root-mean-square momentum, p0 = √( — 2E).
> **In the molecule, the matrix operator corresponding to p0 is related to the
> operator for the united-atom by a sum of unitary transformations, one for each
> nucleus in the molecule.**"

Read this carefully, because it cuts both ways.

- **It IS an operator statement, one unitary per centre.** The last sentence says the
  molecular `p0` matrix is built from the united-atom operator by *a sum of unitary
  transformations, one for each nucleus*. That is, in substance, the translation
  structure: one unitary per nuclear displacement. The corpus cannot claim to have
  been the first to see that the many-centre object is the united-atom object acted on
  by a per-centre unitary.
- **It is NOT stated group-theoretically.** No group is named. The stated method is
  "solely by the ordinary angular momentum calculus"; the deliverable is "algebraic
  expressions for all the integrals" and tabulated "R4 transformation coefficients".
  The framing is *integral evaluation*, exactly as the earlier corpus scan concluded.
- **The word "translation" does not appear in the abstract**, and neither does
  "representation", "group", or `exp(ip·R)`.

### 1b. The corroborating evidence: SW 1965's complete reference list

The paper cites **exactly 8 works** (OpenAlex `referenced_works`, each resolved):

| Year | Authors | Title | Venue |
|:--|:--|:--|:--|
| 1935 | V. Fock | Zur Theorie des Wasserstoffatoms | Z. Phys. |
| 1929 | Podolsky & Pauling | The Momentum Distribution in Hydrogen-Like Atoms | Phys. Rev. |
| 1937 | B. L. Hicks | The Shape of the Compton Line for Helium and Molecular Hydrogen | Phys. Rev. |
| 1941 | Coulson | Momentum Distribution in Molecular Systems | Proc. Camb. Phil. Soc. |
| 1941 | Coulson & Duncanson | Momentum Distribution in Molecular Systems | Proc. Camb. Phil. Soc. |
| 1941 | Duncanson & Coulson | Momentum distribution in molecular systems | Proc. Camb. Phil. Soc. |
| 1949 | R. McWeeny | The Computation of Wave Functions in Momentum Space – II: H2+ | Proc. Phys. Soc. A |
| 1953 | Bates, Ledsham & Stewart | Wave functions of the hydrogen molecular ion | Phil. Trans. R. Soc. A |

**Not one group-theory reference.** No Wigner, no Bargmann, no Biedenharn, no
Talman, no Bander–Itzykson (which postdates it anyway). The bibliography is
momentum-space quantum chemistry plus Fock 1935. This is consistent with the
abstract's own "ordinary angular momentum calculus" and strongly against SW 1965
being a group-representation paper.

### 1c. The decisive adjacent primary: Wulfman does it himself — two years later

**C. E. Wulfman and Y. Takahata, "Noninvariance Groups in Molecular Quantum
Mechanics. I", *J. Chem. Phys.* **47** (2), 488–498 (1967), DOI 10.1063/1.1711921.**
**[ABSTRACT, verbatim]**

> "The problem of one electron moving in the field of many stationary nuclei is
> formulated in terms of the operations of the Lie algebras of **E4**, R5, and O4,1,
> all noninvariance groups of quantal electrostatics. **This makes it possible to
> carry through and completely analyze one-electron molecular calculations using only
> concepts and operations of the theory of continuous groups.** The approach is
> designed for the study of general laws of behavior of simple molecular systems..."

`E4` is the Euclidean group in four dimensions — rotations **and translations** of the
Fock hypersphere's ambient space. So the group-theoretic reading of the many-centre
one-electron momentum-space problem *is* Wulfman's, but it is **1967, not 1965**, and
it is **Wulfman & Takahata, not Shibuya & Wulfman**. The phrasing "*This makes it
possible to* carry through ... using only concepts and operations of the theory of
continuous groups" reads as announcing a new capability, which is further evidence
that the 1965 paper had not done it that way.

### 1d. The translation reading in the later literature — verified verbatim

This is the part that matters most for the corpus's attribution, and it is
unambiguous. The "SW matrix = the translation operator in a Coulomb-Sturmian basis"
reading is **explicit, titled prior art by 2002**:

- **C. A. Weatherford**, "*Diagonalization of the **Translation Operator** in a Coulomb
  Sturmian Basis: Application to Multicenter Integrals Over Exponential-Type
  Orbitals*", APS DAMOP Meeting Abstracts (2002). **[TITLE, verbatim]**
- **C. A. Weatherford and E. Red**, "*Representation of the **Translation Operator** in a
  Soft-Slater Basis: Application to the **Translation of Coulomb Sturmians***", APS
  DAMOP Meeting Abstracts **34** (2003). **[TITLE, verbatim]**
- **E. Red and C. A. Weatherford**, "Derivation of a general formula for the
  Shibuya–Wulfman matrix", *Int. J. Quantum Chem.* **100** (2), 208–213 (2004),
  DOI 10.1002/qua.20122. **[ABSTRACT, verbatim]**: "An explicit and general formula
  for the Shibuya–Wulfman matrix, which may easily be programmed, is derived. The
  derivation closely follows the Fourier transform methods of Weniger. A graph of a
  number of typical matrix elements is given to show the dependence of the
  Shibuya–Wulfman matrix on **the translation distance multiplied by the screening
  parameter**."
  → Note that last clause: Red & Weatherford already plot SW matrix elements against
  **(translation distance) × (screening parameter)** — i.e. against GeoVac's own `kR`,
  as the single combined variable. That specific parameterisation is prior art.
- **P. Hoggan**, "How Exponential Type Orbitals Recently Became a Viable Basis Set
  Choice in Molecular Electronic Structure Work and When to Use Them",
  arXiv:1010.5425, Conclusions. **[PRIMARY, full text read — PDF downloaded and
  grepped locally]**: "...some mathematical structure has been emerging in the BCLFs
  used to translate Slater type orbitals and even more in the **Shibuya-Wulfman matrix
  used to translate Coulomb Sturmians**."
- **J. Avery**, "Many-center Coulomb Sturmians and Shibuya–Wulfman integrals",
  *Int. J. Quantum Chem.* **100** (2), 121–130 (2004; online 2 Dec 2003),
  DOI 10.1002/qua.10820. **[ABSTRACT, verbatim]**: "When momentum space is projected
  onto the surface of a unit 4-D hypersphere by means of Fock's mapping, Coulomb
  Sturmian basis functions can be simply represented in terms of hyperspherical
  harmonics. The properties of these harmonics can be used to evaluate
  Shibuya–Wulfman integrals and other integrals that arise when the Sturmian basis
  functions are used in molecular calculations."
  *Correction owed to the corpus:* the earlier scan recorded this as "IJQC **98**
  (2004)". Crossref and OpenAlex both give **vol. 100, issue 2, pp. 121–130**. The
  DOI is right; the volume/pages are wrong wherever the corpus carries them.
- The sentence "**the Shibuya-Wulfman integrals ... form a representation of the group
  of translations**" is attested repeatedly in search indexing of Avery's book/chapter
  material (with an "Eq. (1.122)" numbering pointing at a monograph chapter, most
  likely *Hyperspherical Harmonics and Generalized Sturmians* or *Hyperspherical
  Harmonics and Their Physical Applications*). **This one I could NOT verify in
  primary** — every candidate host (Wiley, Elsevier, World Scientific, ResearchGate,
  Academia.edu, Copenhagen CURIS/researchprofiles, Aalborg VBN) refused. Flagged as
  search-index-level, do not quote it as verified. It is not load-bearing: the
  Weatherford/Red titles and the Hoggan sentence already carry the point.

### 1e. Answer to the question as asked

> *Do Shibuya and Wulfman themselves give a group-theoretic reading of their
> integrals — specifically, do they describe the two-centre SW matrix element as the
> action of a translation operator on the Fock hypersphere?*

**On the evidence reachable: no, not as group theory — but yes, as an operator
statement, and the identification is in any case not GeoVac's.** Precisely:

1. SW 1965's own abstract frames the method as integral evaluation "solely by the
   ordinary angular momentum calculus", and its bibliography contains no group theory.
   No group is named, "translation" is not used.
2. SW 1965's own abstract *does* state that the molecular `p0` matrix operator is the
   united-atom operator acted on by "**a sum of unitary transformations, one for each
   nucleus in the molecule**" — which is the translation structure in substance, one
   unitary per centre, stated by them.
3. The explicit continuous-group formulation is **Wulfman & Takahata 1967** (E4, R5,
   O(4,1)), i.e. the same author, two years later.
4. The explicit "SW matrix = the translation operator represented in a
   Coulomb-Sturmian basis" reading is **Weatherford & Red, 2002–2004**, in paper
   titles, and is routine in the review literature by 2009 (Hoggan).

**Therefore the corpus's attribution must move.** Paper 60 currently writes
(`paper_60_sturmian_secular_quantum.tex`, ~L830) that what is ours is "the
*identification* — that the two-center Shibuya–Wulfman metric in the sine basis is
such a finite section ... **equivalently that the Shibuya–Wulfman operator is
multiplication by the translation phase `e^{ip·R}` on the Fock sphere**". The clause
after "equivalently" is prior art. What survives as GeoVac's is narrower and should
be stated as such: **not** that SW is a translation operator, but that *in the sine
basis on the Fock polar angle it is the finite section of a multiplication operator
with symbol `j0(kR cot(χ/2))`*, i.e. the **Toeplitz-minus-Hankel/symbol** reading and
everything downstream of it (`eq:sigma_law`, the ratio-symbol argument, the
conditioning law). The symbol itself is the part the earlier Q4 scan already found
had "never appeared", and that finding stands.

### 1f. Routes tried for the 1965 body, and how each failed

| Route | Result |
|:--|:--|
| `royalsocietypublishing.org/doi/10.1098/rspa.1965.0151` (WebFetch) | HTTP 403 |
| same URL via `curl` + browser UA | HTTP 403, Cloudflare "Just a moment..." interstitial |
| Unpaywall API (`10.1098/rspa.1965.0151`) | `is_oa: false`, `oa_locations: []` — no OA copy anywhere |
| Semantic Scholar Graph API | `openAccessPdf.status: "CLOSED"`, abstract "elided by the publisher" |
| OpenAlex | **abstract recovered verbatim + 8-item reference list recovered** (this is the win) |
| S2 citation-contexts API (79 citing papers) | only 6 contexts stored; none quotes SW's own framing beyond "found momentum-space solutions to the one-electron many-center wave equation" (Hoggan) |
| arXiv full-text/metadata search `all:"Shibuya-Wulfman"` | zero hits |
| ar5iv full text of math-ph/0606062, arXiv:2007.00698, arXiv:1811.05777 | no occurrence of Shibuya/Wulfman |
| HathiTrust full-text search, phrase "Shibuya-Wulfman", full-view limit | **"No results"** |
| archive.org full-text search API | connection failed (HTTP 000) |
| Google Books API (snippet search) | HTTP 429, anonymous daily quota exhausted |
| Wiley `pdfdirect` for Avery 2004 (bronze OA per OpenAlex) | HTTP 403 (Cloudflare) — via curl *and* WebFetch |
| ScienceDirect PDF for Avery 2010 JCAM (bronze OA per Unpaywall) | HTTP 403 — via curl *and* WebFetch |
| World Scientific sample chapter (10.1142/10690) | HTTP 403 |
| Copenhagen CURIS landing page (submitted version of Avery 2004) | page served, **no file attached** |
| `researchprofiles.ku.dk` search | HTTP 403 |
| Aalborg VBN (Avery 1996 ×2) | pages served, **no PDFs attached** |
| ResearchGate / Academia.edu | 403 |
| J. E. Avery's DIKU page | only links out to `sturmian.kvante.org`, which returns empty |

**What would close it:** a library copy of *Proc. R. Soc. A* **286** (1965), or of
Wulfman & Takahata, *J. Chem. Phys.* **47**, 488 (1967) — the latter is the more
valuable read, since it is where the group structure is actually developed and it
would fix, in the authors' own words, how much of it was already in the 1965 paper.

### 1g. Two by-catch items the corpus should know about

1. **H. J. Monkhorst and B. Jeziorski, "No linear dependence or many-center integral
   problems in momentum space quantum chemistry", *J. Chem. Phys.* **71**, 5268
   (1979).** Found while chasing SW. The *title alone* is a direct claim on the
   territory of the repo memory note `avery_method_and_prior_art_gaps.md`
   ("conditioning is a blind spot in the whole Avery canon" / "no prior art for
   secular-matrix norm growth"). I did **not** open it, and the title is about
   *absence* of linear dependence rather than a scaling law, so it may well not
   contradict the GeoVac claim — but a claim of that shape, published in JCP in 1979,
   must be read before the corpus repeats "no prior art". **Owed.**
2. The earlier conditioning scan's Q1a flag (arXiv:1208.4156, nuclear Coulomb-Sturmian
   "no linear dependency problem") is the same claim family. Both should be closed in
   one pass.

---
---

# TARGET 2 — the Toeplitz pencil / ratio-symbol theorem

## Verdict: REACHED

It is a named result with a clean lineage, and the single best citation is **open
access and was read in full**. Better: the *exact* form of the theorem lives in the
**τ (DST-I) algebra**, whose matrices are **Toeplitz-minus-Hankel** — which is
GeoVac's structure, not the pure-Toeplitz one. The corpus's object is closer to the
exactly-solvable case than to the asymptotic one.

## 2a. The primary, opened

**F. Ahmad, E. S. Al-Aidarous, D. A. Alrehaili, S.-E. Ekström, I. Furci and
S. Serra-Capizzano, "Are the eigenvalues of preconditioned banded symmetric Toeplitz
matrices known in almost closed form?", *Numerical Algorithms* **78** (3), 867–893
(2018), DOI 10.1007/s11075-017-0404-z. Open access (CC-BY).**
**[PRIMARY, full text read — PDF downloaded and text-extracted locally.]**

Setting (their §1, verbatim structure): `f, g` real-valued **cosine trigonometric
polynomials** (RCTPs) on `[0,π]`; `M_g = max g > 0`, `m_g = min g ≥ 0` so `T_n(g) > 0`;
`P_n(f,g) = T_n^{-1}(g) T_n(f)`; and — their words — "**we define the new symbol
r = f/g**".

### The exact identity — their Eq. (20)–(22)

Inside the proof they decompose

> `T_n(f) = τ_n(f) + H_n(f)`,  `T_n(g) = τ_n(g) + H_n(g)`   — Eq. (20)

where `τ_n(φ) = Q diag_{1≤j≤n}(φ(jπ/(n+1))) Q` with `Q` the **DST-I** matrix
`Q_{ij} = sqrt(2/(n+1)) sin(ijπ/(n+1))`, `Q = Qᵀ = Q⁻¹`, citing **Bini & Capovani,
*Linear Algebra Appl.* **52–53** (1983) 99–126** for the `τ` class; and `H_n(φ)` is a
corner Hankel matrix built from `φ̂_2 … φ̂_m` with `rank H_n(φ) ≤ 2(deg φ − 1)`.
Then, verbatim:

> `P̃_n = τ_n⁻¹(g) τ_n(f) = Q diag(f(jπ/(n+1))/g(jπ/(n+1))) Q = Q diag(r(jπ/(n+1))) Q`
>
> "Hence, for `j = 1, …, n`:  **λ_j(P̃_n) = r(jπ/(n+1))**."   — Eq. (22)

**That is the corpus's theorem, exactly, with no error term at all**, for matrices in
the τ algebra: the generalized eigenvalues of the pencil are the **ratio symbol
sampled on the grid `θ_{j,n} = jπ/(n+1)`**, and the weight `g` has cancelled
identically. Rate and constant both follow, because they are just the behaviour of
`r` near its extremum evaluated at the last grid point.

### Theorem 1 — the Toeplitz version, with its error bound

**Theorem 1** (Appendix, verbatim statement):

> "Let `f`, `g` be real-valued cosine trigonometric polynomials (RCTP) on `[0,π]` with
> `M_g = max g > 0` and `m_g = min g ≥ 0`. If `r = f/g` is **monotone** on `[0,π]` then
> `∃ C > 0` such that
>
>   `|λ_j(P_n(f,g)) − r(jπ/(n+1))| ≤ C h`  ∀ j, n,
>
> where ... `h = 1/(n+1)` and `θ_{j,n} = jπ/(n+1) = jh`."

Proof route: Eq. (22) is exact for the τ pencil; the Hankel parts are a perturbation
of rank `≤ R_{f,g} = 2(max{deg f, deg g} − 1)`, handled by min–max plus the
interlacing theorem for Hermitian matrices.

**Remark (immediately after the proof, verbatim):** "With regard to Theorem 1, the
case where `r` is bounded and **nonmonotone** is even easier. If we consider `r̂`, the
monotone nondecreasing rearrangement of `r` on `[0,π]`, taking into account that the
derivative of `r` has **at most a finite number `S` of sign changes** ... the
eigenvalues of `τ_n(r)` are exactly given `r(jπ/(n+1))` so that, by ordering these
values nondecreasingly, ... the proof follows exactly the same steps as in Theorem 1."

### The higher-order expansion — and its status

**Conjecture (their Eq. (1)), numerically supported, not proved in that paper:**
for every integer `α ≥ 0`,

> `λ_j(P_n(f,g)) = r(θ_{j,n}) + Σ_{k=1}^{α} c_k(θ_{j,n}) h^k + E_{j,n,α}`,
> `E_{j,n,α} = O(h^{α+1})`, with `{c_k}` "a sequence of functions ... **which depends
> only on r**".

Status, from the second primary I opened —
**M. Bogoya, S. Serra-Capizzano and P. Vassalos, "Fast Toeplitz eigenvalue
computations, joining interpolation-extrapolation matrix-less algorithms and
simple-loop conjectures: the preconditioned setting", arXiv:2203.11338 (2022).
[PRIMARY, full text read.]** There it is stated as **Conjecture 1.1**, and:

> "When `g ≡ 1` and `l` satisfies further technical additional assumptions, those of
> the **simple-loop** method, Conjecture 1.1 was formally proved by Bogoya, Böttcher,
> Grudsky, and Maximenko in a series of papers [6, 8, 10, 12]. **For a positive
> function `g`, relation (1.1) was proven, using only purely matrix-theoretic tools,
> and only for `K = 1` in [1].**" ([1] = Ahmad et al. 2018, i.e. Theorem 1 above.)

The same paper also records the two *proved* ingredients, with their homes:

- **Theorem 2.3 (localization)**, from **S. Serra-Capizzano, "The extension of the
  concept of the generating function to a class of preconditioned Toeplitz matrices",
  *Linear Algebra Appl.* **267** (1997) 139–161**: for `g, l ∈ L¹(−π,π)`, `g ≥ 0` a.e.
  and not a.e. zero, `f = l/g`, `m = ess inf f`, `M = ess sup f` with `m < M`, then
  `λ_j(T_n⁻¹(g)T_n(l)) ∈ (m, M)` **for all j and all n ≥ 1**.
- **Theorem 2.2 (distribution)**, from **S. Serra-Capizzano, "An ergodic theorem for
  classes of preconditioned matrices", *Linear Algebra Appl.* **282** (1998) 161–183**:
  `{T_n⁻¹(g)T_n(l)} ∼_λ (f, Q̃)` — the whole sequence is Weyl-distributed as the ratio
  symbol. (Reduces to Szegő / Tyrtyshnikov–Zamarashkin when `g ≡ 1`.)

## 2b. Answers to the questions as asked

**Is it a named theorem?** Yes, but it is three statements, not one, and they are not
equally strong:

| Statement | Status | Home |
|:--|:--|:--|
| Extreme generalized eigenvalues **contained in** `[inf r, sup r]`, `r = f/g`, for every `n`, `g` arbitrary ≥ 0 | **PROVED**, classical | Di Benedetto–Fiorentino–Serra 1993 / Serra-Capizzano LAA **267** (1997) Thm as restated Thm 2.3 of arXiv:2203.11338; also Serra, *Math. Comp.* **66** (1997) Thm 2.2 (the corpus read this one in primary in the earlier scan) |
| Extreme generalized eigenvalues **converge to** `inf r` / `sup r` | **PROVED** | Serra, *Math. Comp.* **66** (218) (1997) 651–665, Thm 2.2 ("lim λⁿ₁ = r, lim λⁿ_n = R") |
| Eigenvalues **equal the ratio symbol on the grid** `jπ/(n+1)`, to `O(h)` — this is what carries the **rate** | **PROVED** for RCTP `f,g` (banded), `g ≥ 0`, `r` monotone or boundedly nonmonotone | **Ahmad et al. 2018, Theorem 1** |
| Same, **exactly, with no error**, in the τ algebra | **PROVED, exact** | **Ahmad et al. 2018, Eq. (22)**; τ class from Bini & Capovani, LAA **52–53** (1983) 99–126 |
| Higher-order expansion `r(θ_{j,n}) + Σ c_k(θ)h^k` — what would pin the **constant** for a Toeplitz (non-τ) pencil | **CONJECTURE** for `g ≢ 1` beyond `K=1`; proved for `g ≡ 1` under simple-loop hypotheses | Ekström–Garoni–Serra-Capizzano, *Exper. Math.* **27**(4) (2018) 478–487 (origin); Bogoya–Böttcher–Grudsky–Maximenko (`g≡1` proofs); arXiv:2203.11338 Conj. 1.1 |

**Under exactly what hypotheses on `W` and `g`?** For the containment/limit
statements: only `g, l ∈ L¹`, `g ≥ 0` a.e. and not a.e. zero. **No smoothness, no
monotonicity, no bandedness.** For the grid statement (Theorem 1): `f` and `g` must be
**cosine trigonometric polynomials** (hence banded), `max g > 0`, `min g ≥ 0`, and
`r = f/g` monotone on `[0,π]` — relaxed by the Remark to bounded `r` whose derivative
has **finitely many** sign changes.

**Does the standard statement cover the CONSTANT or only the exponent/distribution?**
For a *pure Toeplitz* pencil: **only down to `O(h)`**, which is one order too coarse to
pin a constant sitting at `O(h²)`. The constant is covered only by the conjectural
`K ≥ 2` expansion. For a *τ* pencil: **exactly, constant and all**, by Eq. (22) — there
is no error term to beat.

**Does it cover Toeplitz-minus-Hankel (the τ/DST algebra) as well as pure Toeplitz?**
**Yes — and that case is the *stronger* one.** τ matrices *are* Toeplitz-minus-Hankel;
Eq. (22) is exact there, and it is the pure-Toeplitz case that needs the perturbation
argument. Neither arXiv:2203.11338 nor the GLT line addresses Toeplitz-minus-Hankel as
a target class in its own right (I grepped the extracted text of 2203.11338: **zero**
occurrences of "Hankel" — the `τ` glyph itself is not greppable after PDF extraction,
so that half is unchecked), but Ahmad et al. use it as the exactly-solvable core.

## 2c. Where GeoVac's matrices actually sit — measured this pass

`geovac/sturmian_sigma_law.py::weighted_blocks` builds
`A_ab = (2/π)∫₀^π sin(aχ)sin(bχ) W(χ) dχ` and
`B_ab = (2/π)∫₀^π sin(aχ)sin(bχ) W(χ) j₀(kR cot(χ/2)) dχ`,
i.e. `c_{a−b} − c_{a+b}` with `c_k` the cosine coefficients of the symbol — Toeplitz
minus **one** Hankel. The exact τ matrix is Toeplitz minus **two** Hankels:

    tau_n(phi) = T_n(phi) - [Hankel phi_hat_{i+j}] - [Hankel phi_hat_{2(n+1)-i-j}]

(verified numerically this pass by forming `Q diag(φ(θ_j)) Q` and differencing against
`T_n(φ)`: the residual is exactly `+φ̂_2` at `(1,1)` and at `(n,n)` for a degree-2
symbol). So GeoVac's blocks are **exactly τ when the symbol has cosine degree ≤ 1**,
and approximately τ otherwise, with the discrepancy governed by the coefficient tail.
Measured relative off-diagonality of `QAQ`: `2e-15` for `W = 1` and `W = 1+0.8cosχ`;
`8e-3` (`W = 2+sinχ`) and `2e-2` (`W = e^{−χ}`) at `n = 12`.

**The measurement that matters.** Generalized eigenvalues of the pencil `(B, A)`
against `max_j r(θ_{j,n})`, `θ_{j,n} = jπ/(n+1)`, at `kR = 2`:

| `n` | `W=1` | `W=1+0.8cosχ` | `W=2+sinχ` | `W=e^{−χ}` | control `W=1+cosχ` (vanishes at `χ=π`) |
|--:|--:|--:|--:|--:|--:|
| 10 | 5.3e-4 | 3.3e-3 | 1.2e-3 | 3.3e-3 | 1.33e-2 |
| 40 | 9.7e-6 | 6.0e-5 | 3.3e-5 | 7.8e-5 | 9.98e-4 |
| 160 | **1.6e-7** | 9.6e-7 | 7.0e-7 | 1.6e-6 | **6.6e-5** |

(relative difference `|σ_max − max_j r(θ_j)| / max_j r(θ_j)`.)

Three readings, all load-bearing:

1. **For smooth positive `W` the deficit from the grid value falls like `n^{-3}`**
   (ratio ≈ 7.9 per doubling), i.e. **one order below** the `O(n^{-2})` main term. That
   is why the constant is `W`-independent — and it is two orders better than Theorem 1's
   proved `O(h)` bound, consistent with the conjectural `K = 2` expansion.
2. **The vanishing-weight control behaves differently in exactly the right way:** its
   deficit falls like `n^{-2}` (ratio ≈ 3.9), *the same order as the main term*, so it
   shifts the constant by a factor ≈ 2 — reproducing Paper 60's measured `0.828` against
   `π²/24 = 0.4112`. The control is doing real work; it is not an artifact.
3. **`n → n+1` is the corpus's own unexplained `O(1/n)` residue.** Paper 60 notes that
   the `~1%` gap at `n = 160` "is the asymptotic's own `O(1/n)` term rather than
   scatter". It is: `1 − max_j r(θ_{j,n}) = (kR)²π²/(24(n+1)²)` — at `n = 160`,
   `kR = 2`, that predicts `6.3459e-5` against the measured `6.3462e-5`. Collapsing with
   `(n+1)²` instead of `n²` moves the residue from **−0.98% to +0.26%**, a 3.7×
   reduction. `eq:sigma_law` would be sharper written with `(n+1)`.

## 2d. Hypothesis gaps — stated flatly, because they are real

GeoVac's ratio symbol `r(χ) = j₀(kR cot(χ/2))` **fails the stated hypotheses of
Theorem 1 on three counts**:

- it is **not** a cosine trigonometric polynomial (coefficients decay as `j^{-5/4}`,
  per `eq:chirp_decay`), so the Hankel corrections are not low rank;
- `r′` has **infinitely many** sign changes as `χ → 0` (the chirp), so neither the
  monotonicity hypothesis nor the "finite number `S` of sign changes" of the Remark
  holds;
- `min r < 0` (`j₀` dips to `−0.2172`), though this bites only the `m_g ≥ 0` condition
  on the *preconditioner* symbol, not on `r`, so it is the least serious of the three.

What survives regardless: the **containment and limit** results (Serra-Capizzano LAA
267/282; Serra *Math. Comp.* 66 Thm 2.2) need only `L¹` and non-negativity, which
GeoVac satisfies; and the corpus's own basis-free positivity argument (recorded in
`toeplitz_finite_section_memo.md` Q3.1) gives the containment for any Gram/moment
matrices in any orthonormal basis, with no Toeplitz theory at all.

## 2e. Bibliographic corrections owed

- The candidate the corpus had declined to cite, **"LAA 270 (1998)", is the wrong
  paper.** *Linear Algebra Appl.* **270** (1998) 15–27 is **Tyrtyshnikov & Zamarashkin,
  "Spectra of multilevel Toeplitz matrices: advanced theory via simple matrix
  relationships"** — a Szegő-type distribution theorem, not the pencil/ratio result.
  The Serra-Capizzano papers wanted are **LAA 267 (1997) 139–161** (localization) and
  **LAA 282 (1998) 161–183** (distribution). Good that it was not added.
- **Page-number conflict in the literature itself** for Di Benedetto–Fiorentino–Serra,
  "C.G. preconditioning for Toeplitz matrices", *Comput. Math. Appl.* **25**(6) (1993):
  Ahmad et al. ref. [11] gives **33–45**; Bogoya–Serra-Capizzano–Vassalos ref. [15]
  gives **35–45**. Check against the journal before citing; do not copy either blindly.

## 2f. The one citation to add

If Paper 60 adds exactly one reference for the ratio-symbol mechanism, it should be
this one — it is open access, it was read in full, and its Eq. (22) is the mechanism
in exact form for GeoVac's own Toeplitz-minus-Hankel structure:

> F. Ahmad, E. S. Al-Aidarous, D. A. Alrehaili, S.-E. Ekström, I. Furci and
> S. Serra-Capizzano, "Are the eigenvalues of preconditioned banded symmetric Toeplitz
> matrices known in almost closed form?", *Numerical Algorithms* **78** (3) (2018)
> 867–893, DOI 10.1007/s11075-017-0404-z. (Open access.) — Theorem 1 and Eq. (22).

with **Bini & Capovani, *Linear Algebra Appl.* **52–53** (1983) 99–126** for the τ
(DST-I) algebra if the exact identity is used, and
**Serra-Capizzano, *Linear Algebra Appl.* **267** (1997) 139–161** for the
hypothesis-free localization. *Provenance discipline:* I opened Ahmad et al. and
arXiv:2203.11338 in full. I did **not** open Bini–Capovani, Serra-Capizzano LAA
267/282, or Di Benedetto–Fiorentino–Serra; their statements here are as restated in
the two primaries I did open, and are marked as such.

## 2g. Verification artifacts

Scratch scripts (session scratchpad, not committed): `tau_check.py` (τ-membership of
`weighted_blocks` output under DST-I), `tau2.py` (pencil `σ_max` vs `max_j r(θ_{j,n})`
across five weights and `n = 10…160`), `tau5.py` (recovery of the exact τ correction
convention `T_n − H_{i+j} − H_{2(n+1)−i−j}`). Downloaded primaries:
`webfetch-…-i2l6bq.pdf` → `ahmad.txt` (Ahmad et al. 2018),
`webfetch-…-qfttgb.pdf` → `out2203.txt` (arXiv:2203.11338), `hoggan.pdf`/`hoggan.txt`
(arXiv:1010.5425).

---

## Two sentences the corpus should now write

**Target 1.** *The reading of the Shibuya–Wulfman matrix as the translation operator
in a Coulomb–Sturmian basis is prior art — Shibuya and Wulfman's own 1965 abstract
already builds the molecular `p₀` operator from "a sum of unitary transformations, one
for each nucleus in the molecule", Wulfman and Takahata gave the explicit
continuous-group formulation in terms of the Lie algebras of `E₄`, `R₅` and `O(4,1)`
two years later [J. Chem. Phys. **47**, 488 (1967)], and Weatherford and Red titled
papers on the diagonalization and representation of that translation operator in
2002–2004 — so what is ours is not the translation identification but the
**symbol**: that in the sine basis on the Fock polar angle the operator is the finite
section of multiplication by `j₀(kR cot(χ/2))`.*

**Target 2.** *The ratio-symbol law is a published theorem, not an observation: for
matrices in the τ (discrete-sine) algebra the generalized eigenvalues of the pencil are
**exactly** `r(jπ/(n+1))` with `r = f/g`, so the weight cancels identically [Ahmad
et al., Numer. Algorithms **78** (2018) 867–893, Eq. (22); τ algebra: Bini & Capovani,
LAA **52–53** (1983) 99–126], and for the Toeplitz pencil the same grid formula holds
to `O(1/n)` (their Theorem 1, for cosine-polynomial symbols with `r` monotone or
boundedly nonmonotone) — the localization and limit statements requiring nothing but
`L¹` positivity [Serra-Capizzano, LAA **267** (1997) 139–161]; what remains
conjectural, and is what our `n^{-3}` residue actually exhibits, is the higher-order
expansion that would pin the **constant** for a non-τ pencil.*
