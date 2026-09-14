# C23 run #3 — the M2 tagging paragraph of Paper 60 (at-authorship trigger)

**Date:** 2026-09-12 · **Criterion:** `docs/qa/criteria.md` C23, trigger 2 (at
authorship) · **Target:** the `[SYMBOLIC + MEASURED]` paragraph beginning "The two
transcendentals of this section sit on opposite sides of the compactness boundary",
`papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex` ~L892;
`debug/sprint_contraction_seam_memo.md` §3; CHANGELOG v5.11.4.

**Read-only.** Nothing under `papers/`, `geovac/`, `tests/`, `CLAUDE.md`,
`CHANGELOG.md` was modified. Local verification scripts were run from the session
scratchpad only.

**Not re-covered** (settled by earlier scans): KMS itself, Böttcher–Widom
arXiv:math/0412269, DLMF 10.32.10 / 10.40.2 and the π/4-as-branch-phase correction,
Serra *Math. Comp.* **66** 651 (1997), Jaffard, Gröchenig–Leinert, the RBF flat limit,
Ron–Shen, Balian–Low, Halmos, Löwdin/Slater–Koster. Prior scans read first:
`toeplitz_finite_section_memo.md`, `frames_riesz_overcompleteness_memo.md`,
`c23_paper60_analysis_memo.md`, `contraction_seam_e3_memo.md`.

---

## Verdict table

| # | Claim | Verdict | One line |
|:--|:--|:--|:--|
| **T1** | `c_1 = pi^2` "is the first Dirichlet eigenvalue of `-d^2/dx^2` on the unit interval" | **PRIOR ART** | Explicit in substance in Böttcher–Widom — the source Paper 60 already cites — as the `alpha = 1` case of a boundary-value problem they write out and whose `c_1 = pi^2` they state. Not folklore, not ours. The words "Dirichlet"/"Laplacian" do not appear there; the BVP does. |
| **T2** | `n^2 min<theta^2> -> pi^2` over `span{sin(a theta)}_{a<=n}` is a band-limited concentration result "containing no symbol at all" | **PRIOR ART** (twice) | (a) It is the `b == 1` case of the *same* KMS/Böttcher–Widom statement, i.e. the variational characterisation of `c_1` — so it is **not an independent route**. (b) Externally it is the classical minimum **second-moment (Gabor) bandwidth/duration** problem, extremal = the **half sinusoid**, sharp constant = the Dirichlet/Wirtinger eigenvalue of the band; `M`-signal version Nuttall–Amoroso 1965. **NOT Slepian/Landau/Pollak** — their functional is energy concentration, eigenvalues `lambda_n(c)`, not `pi^2`. |
| **T3** | The two-sided truncation-side / continuum-side classification, read as a predictor of which wall is removable | **ABSENT** as a taxonomy — but see the weakest-link section: two of its three legs are prior art, and its operational inference misnames its own mechanism | No "where the constant comes from tells you whether the wall is removable" accounting found in the numerical-analysis or asymptotic-analysis literature. The NA literature *does* classify removability — keyed on **symbol zeros vs winding number vs smoothness**, not on constant provenance. |

---

## T1 — the KMS constant as a Dirichlet eigenvalue

### What the primary says

**Source (verified, full-text render):** A. Böttcher and H. Widom, *From Toeplitz
eigenvalues through Green's kernels to higher-order Wirtinger–Sobolev inequalities*.
Read at https://ar5iv.labs.arxiv.org/html/math/0412269 (full text), abstract
re-verified at https://arxiv.org/abs/math/0412269, bibliographic record verified at
Crossref (DOI `10.1007/978-3-7643-7980-3_4`).

Quoted from the rendered full text:

> `lambda_min(T_n(a)) ~ c_alpha / n^{2 alpha} * b(1)` as `n -> infinity`, with `c_1 = pi^2`
> [attributed to] Kac, Murdock, and Szegő

> For a natural number `alpha`, consider the boundary value problem
> `(-1)^alpha u^{(2 alpha)}(x) = v(x)`, `u(0) = u'(0) = ... = u^{(alpha-1)}(0) = 0`,
> `u(1) = u'(1) = ... = u^{(alpha-1)}(1) = 0`

> The minimal eigenvalue of this boundary value problem can be shown to be just `c_alpha`.

> [Wirtinger–Sobolev] the best constant `C` for which (6) is true for all
> `u in C^alpha[0,1]` satisfying `u^{(j)}(0) = u^{(j)}(1)` for `0 <= j <= alpha-1`
> is equal to `1/c_alpha`.

At `alpha = 1` the displayed BVP **is** `-u'' = lambda u` on `[0,1]` with
`u(0) = u(1) = 0` — the Dirichlet problem for `-d^2/dx^2` on the unit interval, first
eigenvalue `pi^2`. So our sentence is a one-line specialisation of an explicitly
stated formula in a source the paper already cites.

### "Written in these words" vs "an expert would recognise it"

Text-searched the rendered full text for `Dirichlet`, `Laplacian`, `Poincaré`:
**none occurs.** `Wirtinger` and `boundary value problem` do. So:

- **written in these words: NO** — Böttcher–Widom never say "Dirichlet";
- **stated in substance, explicitly, in the cited source: YES** — they write the BVP
  and give `c_1 = pi^2`;
- **the general framing "the finite section behaves as a Dirichlet problem on an
  interval": standard**, not folklore. For the model symbol it is exact and textbook:
  `T_n(2 - 2 cos theta)` *is* the `(-1, 2, -1)` discrete Dirichlet Laplacian with
  eigenvalues `4 sin^2(j pi / (2(n+1)))`. The modern Toeplitz/GLT literature routinely
  states that finite sections of symbols with a zero of order `2k` are finite-difference
  discretisations of order-`2k` differential operators with the matching boundary
  conditions.

### Re-verified here (C23 hard rule: reproduce the load-bearing leg)

Exact tridiagonal section, no fitting:

| `n` | `lambda_min * (n+1)^2` | closed form `4 sin^2(pi/(2(n+1))) (n+1)^2` |
|---:|---:|---:|
| 40 | 9.864776 | 9.864776 |
| 160 | 9.869291 | 9.869291 |
| 640 | 9.869585 | 9.869585 |

against `pi^2 = 9.869604`. The identification is exact at every `n`, not asymptotic.

### Exact attribution sentence to use

> The constant is the Kac–Murdock–Szegő constant $c_1=\pi^2$. Böttcher and
> Widom~\cite{bottcher_widom2005} characterise $c_\alpha$ as the minimal eigenvalue of
> the boundary value problem $(-1)^\alpha u^{(2\alpha)}=\lambda u$ on $[0,1]$ with
> $u^{(j)}(0)=u^{(j)}(1)=0$ for $0\le j\le\alpha-1$; at $\alpha=1$ that is the
> Dirichlet problem for $-d^2/dx^2$ on the unit interval, and it is that
> identification --- not ours --- which places the constant on the truncation side.

### Two citation defects found while verifying (both live)

1. **Missing KMS bibitem (C20 class).** Paper 60 L823 attributes `c_1 = pi^2` inline
   to "Kac, Murdock and Szegő" with **no resolvable `\bibitem`**; `grep -i` for
   `kac|murdock|szeg|parter|wirtinger` over the `.tex` returns **nothing** in the
   bibliography. Full data, verified (Böttcher–Widom's own reference [11]→[8] list at
   ar5iv, plus the IUMJ record and Semantic Scholar):

   > M. Kac, W. L. Murdock and G. Szegő, "On the eigen-values of certain Hermitian
   > forms," *J. Rational Mech. Anal.* **2** (1953), 767–800.
   > (Reprinted in *Indiana Univ. Math. J.* 2, No. 4.)

2. **Böttcher–Widom cited as arXiv-only.** The published coordinates exist and are
   verified at Crossref:

   > A. Böttcher and H. Widom, "From Toeplitz eigenvalues through Green's kernels to
   > higher-order Wirtinger–Sobolev inequalities," in *The Extended Field of Operator
   > Theory* (M. A. Dritschel, ed.), Operator Theory: Advances and Applications **171**,
   > Birkhäuser, Basel, 2007, pp. 73–87. DOI 10.1007/978-3-7643-7980-3_4;
   > arXiv:math/0412269.

   *Verified:* title/authors/pages/publisher from Crossref DOI record; volume title and
   editor from the book's own Crossref record (DOI `10.1007/978-3-7643-7980-3`, year
   2007). The volume number 171 comes from a search summary, **not** from a primary
   record — state it or drop it, but do not present it as verified.

---

## T2 — the minimal mean-square spread of a band-limited function

### Re-verified here, independently of the corpus probes

Exact matrix elements, no quadrature. With `f(m) = pi^3/3` at `m = 0` and
`2 pi (-1)^m / m^2` otherwise, `M_ab = [f(a-b) - f(a+b)]/2`, Gram `= (pi/2) I`:

| `n` | `n^2 * min<theta^2>` |
|---:|---:|
| 160 | 9.7619256 |
| 320 | 9.8155435 |
| 640 | 9.8425183 |
| 1280 | 9.8560474 |

Richardson in `1/n` on the last pair: **9.869576** against `pi^2 = 9.869604`.
Consistent with the corpus's 9.86949. **The claim is true.**

### But it is the same theorem, not a second one

The minimising coefficient vector at `n = 320` is, to `1.3e-3` max deviation and
correlation `0.9999991`, `c_a = sin(pi a / (n+1))` — **the Dirichlet ground state in
the band index**. That is the whole content: rescaling `theta = x/n` sends the problem
to "minimise `int x^2 |g|^2 / int |g|^2` over functions whose sine transform is
supported in `[0,1]`", and by Plancherel that is `int_0^1 |G'|^2 / int_0^1 |G|^2` with
`G(0) = G(1) = 0`, i.e. the first Dirichlet eigenvalue `pi^2` — the **same** variational
problem Böttcher–Widom write down, at `b == 1`.

So the paper's parenthetical "a band-limited concentration problem containing no symbol
at all" is accurate as a description and **misleading as evidence**: `b == 1` *is* the
canonical KMS case, and the prior scan's `verify2.py` already ran exactly it
(9.864776 / 9.869291 / 9.869585). The measurement re-derives the constant the same
sentence attributes to KMS. It is a representation change, not an independent route.

### The external name we are not using

**It is not Slepian territory.** Slepian–Landau–Pollak maximise *energy concentration*
in an interval; the extremals are prolate spheroidal wave functions and the eigenvalues
`lambda_n(c)` are transcendental functions of the time–bandwidth product — never `pi^2`.
The second-moment branch of the same lineage is a different, elementary problem with an
exact Dirichlet answer, and the literature states the contrast explicitly:

> "In [6], Gabor presented an alternative set of optimal band-limited functions with a
> **second moment weighting** in the time domain, **the half sinusoids**, whose tails
> decay asymptotically like `1/t^2`." — L. Wei, R. A. Kennedy, T. A. Lamahewa,
> "Further results on signal concentration in time-frequency," *Proc. IEEE ICASSP 2010*,
> §1. (Primary read: PDF text extracted from
> https://users.cecs.anu.edu.au/~rod/papers/2010/05495736.pdf.)

Their [6] is **D. Gabor, "Theory of communication," *J. Inst. Electr. Eng.* 93, Part III,
No. 26 (1946), 429–457** — the origin of the second-moment ("Gabor") bandwidth measure.
*Caveat: the half-sinusoid attribution to Gabor is at secondary-source level; I did not
open Gabor 1946.*

The extremal problem itself, in the form closest to ours (the `M` lowest modes rather
than just the ground state), is:

> **A. Nuttall and F. Amoroso, "Minimum Gabor bandwidth of M orthogonal signals,"
> *IEEE Trans. Inf. Theory* **11** (3), 440–444 (1965).** DOI 10.1109/TIT.1965.1053803.
> Abstract: "The minimum Gabor bandwidth of M orthogonal real signals which together
> occupy a given time interval is derived, and for special values of M, one form of the
> optimum set of M signals is specified. Interrelations of the Gabor definition of
> bandwidth to other definitions are also discussed."

*Verified:* bibliographic record at Crossref (DOI above); abstract via the ACM DL/IEEE
listing. **Abstract level — I did not reach the article text**, so do not attribute a
specific formula to it.

On the pure-mathematics side the constant's name is **Wirtinger's (Poincaré's)
inequality** — which is exactly the reading Böttcher–Widom already supply ("the best
constant in a Wirtinger–Sobolev inequality is `1/c_alpha`"), specialised to `alpha = 1`.

### Exact attribution sentence to use

> The same constant is the sharp Wirtinger--Poincar\'e constant on the band, and the
> variational problem it solves is classical in signal analysis as the minimum
> second-moment (Gabor) bandwidth problem, whose extremals are the half
> sinusoids~\cite{gabor1946,nuttall_amoroso1965}. It is \emph{not} the
> Slepian--Landau--Pollak concentration problem, whose eigenvalues depend
> transcendentally on the time--bandwidth product. What the measurement shows is
> therefore that the near-null direction realises the classical extremal, not that
> $\pi^2$ arises here by a route independent of~\cite{bottcher_widom2005}.

---

## T3 — the two-sided transcendental classification

### Verdict: ABSENT as an accounting

Searched for a stated "where does the constant come from, and does that tell you
whether the wall is removable" principle in: Toeplitz finite-section and
extreme-eigenvalue literature; band-Toeplitz and circulant preconditioning surveys
(Serra-Capizzano lineage, R. Chan–Ng); GLT theory; boundary-layer/matched-asymptotics
literature. **Nothing of that shape was found, and I am not going to manufacture one.**
The classification as a *classification* — provenance of a constant used as a predictor
of removability — is corpus vocabulary (Paper 18), with no external analogue I could
reach.

### But two of its three legs are already someone else's

- **"The `pi^2` is truncation-side, not a property of the rest of the symbol."**
  This is Böttcher–Widom's own statement: formula (1) carries "a certain constant
  `c_alpha` ... **independent of `b`**" — i.e. independent of everything about the
  symbol except the order of the zero. (Quoted from the arXiv full text in
  `toeplitz_finite_section_memo.md`, consistent with the ar5iv render read here.) The
  paper should credit that sentence rather than present the placement as a reading
  obtained here.
- **"The `(2pi)^{-1/2}` and `pi/4` are symbol asymptotics."** True by construction and
  DLMF-sourced (already priced by run #1; not re-covered).

### What the numerical-analysis literature *does* classify

It has its own removability taxonomy, and it is keyed on different features:

| source of ill-conditioning | removable? | mechanism |
|:--|:--|:--|
| a zero of finite (even) order in the symbol | **yes** | divide by the trigonometric polynomial `g` of least degree carrying the same zero — band-Toeplitz preconditioning; `f = g h`, `h` bounded away from `0` and `infinity`, spectrum trapped in `[min h, max h]` |
| nonzero winding number with no zeros on the circle | **no** | topological; condition numbers grow exponentially |
| non-locality of inverses / functional calculus | governed by the symbol's smoothness and Wiener-algebra membership | Jaffard, Gröchenig–Leinert (already covered) |

(Rows 1 and 2 from search-level reading of the band-Toeplitz preconditioning survey
literature; **row 1 is independently confirmed inside this corpus** by the
`verify9/verify10` measurements recorded in `toeplitz_finite_section_memo.md` Q3.0b —
`cond` flat at `7.67` with the limit equal to `max h / min h`. Row 2 is search-summary
level only and should not be cited.)

---

## Weakest link

**W1 (the important one). T3's operational sentence names the wrong mechanism, and its
second half is false as written.** The paper says:

> a truncation-side price is a property of the matrix, which a preconditioner reaches,
> while a continuum-side price is a property of the symbol, which no congruence of the
> finite section can touch.

Both halves fail on inspection against the corpus's own measurement:

- The preconditioner that breached the conditioning wall is `g(theta) = 2 - 2 cos theta`,
  **chosen because it carries the symbol's zero**. It is built from the symbol, not from
  the truncation, and the theory that licenses it is `f = g h` with `h` bounded — a
  statement entirely about the symbol. So truncation-sidedness is not what makes the
  price reachable; carrying the symbol's zero is.
- "No congruence of the finite section can touch a property of the symbol" is **false in
  general**: preconditioning *is* a congruence, and it changes the symbol (to `f/g`).
  The defensible narrow statement is: a **banded** congruence multiplies the symbol by a
  trigonometric polynomial, which can cancel a zero but cannot improve the symbol's
  smoothness class — so it cannot restore locality at the chirp end. That version is
  arguable from the symbol algebra and does not need the constant-provenance framing at
  all.

**W2. Internal tension the new paragraph silently resolves.**
`toeplitz_finite_section_memo.md` Q3.0 measured the non-locality of `S^{-1/2}`
(decay length `L ~ 0.06 sqrt(cond) ~ 0.134 n`) and concluded "the conditioning and the
non-locality are **one fact, not two** ... set by the symbol's quadratic touch. A
structure-preserving fix and a conditioning fix are the same request." The same memo's
item 5b then warns the opposite: "the chirp causes non-locality and no ill-conditioning;
the quadratic touch causes all the ill-conditioning and no non-locality. Conflating them
is easy and the paper may do it." The new paragraph's clean two-sided split adopts the
second reading without saying that the first exists. **This is a corpus-internal
consistency item, not a C23 finding** — but T3's entire operational payoff rests on it,
so it is the load-bearing uncertainty, not the attribution.

**W3. The "independent route" framing.** Memo §3 offers "Second route for `c_1`: the
exact tridiagonal spectrum" and the paper offers the band-limited minimisation. All
three objects — the tridiagonal spectrum, the band-limited second-moment minimum, and
`c_1` itself — are the **same Dirichlet eigenvalue problem** in matrix, Rayleigh-quotient
and BVP form. The memo's own phrasing ("two representations, one constant") is honest;
the paper's ("containing no symbol at all") reads as independent corroboration and is
not. Under `memory/feedback_independent_route_crosscheck.md` this does not satisfy the
second-algorithmic-route requirement.

**W4 / W5.** The two citation defects in the T1 section above: no KMS bibitem (C20
class), and a load-bearing arXiv-only citation whose published coordinates exist.

**Directions that failed.** I could not open (a) the Nuttall–Amoroso article text (IEEE
paywall; abstract only), (b) Gabor 1946 (secondary attribution only), (c) the
Böttcher–Widom published chapter behind Springer's IDP redirect (used the ar5iv render of
the identical preprint instead, plus Crossref for metadata), (d) Serra-Capizzano's
"Practical band Toeplitz preconditioning and boundary layer effects" (*Numerical
Algorithms*, Springer IDP redirect) — which is the one source most likely to state a
truncation-vs-symbol separation in the literature's own words, and whose title alone
suggests it. **That is the single open direction for T3: if anyone later reaches that
paper, re-run T3 against it before the ABSENT verdict is relied on.**

---

## What survives as ours

1. **The identification of the symbol.** `a(chi) = j0(kR cot(chi/2))` as the generating
   function of the Shibuya–Wulfman cross block, and `b(1) = (kR)^2/24` as its curvature
   at the extremum. Unchanged from run #1; still the real contribution.
2. **The sharp constant surviving a failed hypothesis.** Böttcher–Widom require
   `sum |k| |b_k| < infinity`; our symbol measurably fails it (`|b_k| ~ k^{-1.25}`), and
   the constant holds anyway. No theorem found covers this symbol class. Correctly
   already stated in-paper as open.
3. **The weight-independence of §4** (run #2 C3), unchanged: **ABSENT**.
4. **The chirp asymptotic applied to this symbol**, constant and branch phase, verified
   by three quadrature routes — DLMF-sourced law, our application.
5. **The M2 placement itself.** It is bookkeeping in Paper 18's vocabulary, and there is
   nothing external to collide with. It should be presented as tagging, not as a
   discovery — and, per W1, its operational corollary should be re-derived from the
   symbol algebra rather than from the side of the seam a constant sits on.

---

## Sources

- Böttcher & Widom, *From Toeplitz eigenvalues through Green's kernels to higher-order
  Wirtinger–Sobolev inequalities*: https://arxiv.org/abs/math/0412269 ·
  https://ar5iv.labs.arxiv.org/html/math/0412269 ·
  https://link.springer.com/chapter/10.1007/978-3-7643-7980-3_4
- Kac, Murdock & Szegő (1953), record:
  http://www.iumj.indiana.edu/IUMJ/fulltext.php?artid=52034&year=1953&volume=2
- Nuttall & Amoroso, *Minimum Gabor bandwidth of M orthogonal signals*:
  https://ieeexplore.ieee.org/document/1053803/ · https://dl.acm.org/doi/10.1109/TIT.1965.1053803
- Wei, Kennedy & Lamahewa, *Further results on signal concentration in time-frequency*
  (ICASSP 2010): https://users.cecs.anu.edu.au/~rod/papers/2010/05495736.pdf
- Slepian & Pollak / Landau & Pollak, PSWF I–III (for the disambiguation):
  https://ieeexplore.ieee.org/document/6773660 ·
  https://onlinelibrary.wiley.com/doi/abs/10.1002/j.1538-7305.1962.tb03279.x
- Serra-Capizzano, *Practical band Toeplitz preconditioning and boundary layer effects*
  (NOT reached): https://link.springer.com/article/10.1023/B:NUMA.0000005355.94096.bc
