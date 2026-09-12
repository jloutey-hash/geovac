# Lit scan: is Paper 60's `eq:sigma_law` a known Toeplitz finite-section theorem, and is there a local escape?

Date: 2026-09-11. Read-only scan. Target: Paper 60 `eq:sigma_law`
(`papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex`, ~L780-830).

**Headline: GO on Q1. `eq:sigma_law` is Kac-Murdock-Szego (1953), constant and all.**
The constant pi^2/24 factors as `c_1 * b(1)` where `c_1 = pi^2` is the classical
Kac-Murdock-Szego constant and `b(1) = (kR)^2/24` is our symbol's curvature at its
extremum. Verified numerically here on both legs. Paper 60 currently cites **zero**
Toeplitz-family references (`grep -i -e toeplitz -e szeg -e parter -e widom -e kac
-e circulant` over the .tex returns nothing, against 26 bibitems).

---

## Consolidated verdict

| Q | Verdict | One line |
|---|---|---|
| Q1 | **GO** | `eq:sigma_law` is Kac-Murdock-Szego (1953), `c_1 = pi^2`, applied to `b(1) = (kR)^2/24`. Rediscovery. Must cite. |
| Q2 | **GO (rates fully known)** | `n^{-2k}` for order-`2k` zeros with sharp constants `c_1 = pi^2`, `c_2 = 500.5467`, `c_3 = 61529`; non-integer order via Rambour-Seghier; the rate survives our chirp under Serra-Capizzano's `L^1` hypotheses. `n^{-2}` is NOT an artifact. |
| Q3 | **GO on conditioning (proven theorem, we meet its hypotheses); FORBIDDEN on locality of `S^{-1/2}`** | Serra 1997 Thm 2.2/4.1: a **bandwidth-1** preconditioner gives `O(1)` conditioning — measured, `1.2e6 -> 7.67`. But `S^{-1/2}` decay is *excluded by hypothesis* from every decay theorem (holomorphic calculus needs `0 notin spec`), and `(1-a)^{-1/2} notin L^1`. **Constructive escape found and independently reproduced:** all the non-locality is `P^{-1/2}` for `P` tridiagonal, which is an exact DST. |
| Q4 | **STOP on the symbol, GO on citation defects** | The symbol has never appeared; the *operator* is Shibuya-Wulfman. Found 3 citation defects in Paper 60, one substantive (Halmos). |

**Overall: GO.** Q1 alone forces it — `eq:sigma_law` must be recast from a derivation to
a cited classical theorem plus an identification. Q3 additionally returns a real,
previously-unconsidered lever (banded preconditioning) *and* a proof-shaped obstruction
for the quantum side.

### What Paper 60 should do (highest value first)

1. **Recast `eq:sigma_law` as an application of KMS, not a derivation.** The sentence
   "the band-limited concentration rate at a quadratic symbol maximum gives" should
   become "identifying the symbol and applying the classical extreme-eigenvalue
   asymptotics of Kac-Murdock-Szego [ref] with `c_1 = pi^2` gives". **Paper 60 currently
   cites zero Toeplitz-family references** (`grep -i -e toeplitz -e szeg -e parter
   -e widom -e kac -e circulant` over the .tex returns nothing, against 26 bibitems).
   The surviving novelty is *identifying the symbol* `j0(kR cot(chi/2))` and reading
   its curvature `(kR)^2/24` — which is genuine and should be stated as the claim.
2. **Fix "Toeplitz" -> "Toeplitz minus Hankel"** wherever the section is described
   (independently flagged by both scan legs). Costs one sentence.
3. **Fix the Halmos attribution** — Halmos 1969 Thm 2 gives the canonical form, NOT
   `‖[P,Q]‖`; derive the commutator norm in one line or cite Loring 2014.
4. **Report the banded-preconditioner result and cite Serra 1997** (Q3.0b, Q3.1). It
   converts "the metric retains a polynomial conditioning multiplier ... it cannot
   simply be transformed away" into "the metric is spectrally equivalent to the 1-D
   discrete Dirichlet Laplacian, and a tridiagonal preconditioner gives `O(1)`
   conditioning" — truer, more useful, and it retires the classical half of the cost
   argument outright. **The current sentence "it cannot simply be transformed away" is
   too strong** and should be narrowed to the `S^{-1/2}`/quantum case where it is right.
5. **Do NOT claim it fixes the quantum resource count** — but DO report the constructive
   whitening `W = G^{-1/2}P^{-1/2}` (Q3.3), where all the non-locality is concentrated in
   `P^{-1/2}` (an exact DST) and the residual `G^{-1/2}` is uniformly local and uniformly
   conditioned. That is the most valuable new result in this scan and it is
   independently reproduced. Flag the uniform locality of `G` as measured, not proved.
5b. **Re-examine the §sec:obstruction framing** against Q3.7: the chirp causes
   non-locality and no ill-conditioning; the quadratic touch causes all the
   ill-conditioning and no non-locality. Conflating them is easy and the paper may do it.
5c. **Fix the numerics**: never compute `1 - sigma_max` by subtracting `sigma_max` from 1
   (Q3.5) — catastrophic cancellation, ~6 digits lost at `n = 512`.
6. **State the Fock chart convention** (`p = k cot(chi/2)` vs the literature's
   `p = p_0 tan(chi/2)`) — otherwise the symbol's maximum at `chi = pi` reads as an error.
7. Minor: Shibuya-Wulfman is pp. **376-389**, not just p. 376.

---

## 0. The object, restated in the literature's normalisation

Paper 60's cross block `C` is the n x n section of multiplication by

    a(chi) = j0( kR cot(chi/2) ),   j0 = sinc,

in the sine basis `sin(n chi)/sin(chi)` on `chi in (0,pi)` with weight `sin^2 chi`.
Put `theta = pi - chi`, so the extremum `chi = pi` sits at `theta = 0`, and set

    f(theta) = 1 - a  =  1 - j0( kR tan(theta/2) ).

Then `f >= 0`, `f(0) = 0`, and `1 - sigma_max = lambda_min` of the section of `f`.
In the Toeplitz normalisation `f(t) = |1-t|^{2 alpha} b(t)` with `|1-e^{i theta}|^2
= 4 sin^2(theta/2) ~ theta^2`, our case is **alpha = 1** with

    b(1) = lim_{theta->0} f(theta)/theta^2 = (kR)^2 / 24.

Verified to 8 digits (`verify2.py`): at kR = 0.25/0.5/1/2/4 and theta = 1e-3, 1e-4,
`(1-a)/theta^2` reproduces `(kR)^2/24` (e.g. kR=2: 0.1666666861 vs 0.1666666667).

One structural note the paper does not make: because the basis is the **sine** basis,
`C` is not Toeplitz but **Toeplitz minus Hankel**. Writing `f = sum_k c_k e^{i k theta}`
(`c_k` real, `c_k = c_{-k}`, since `f` is real and even), the section is

    M_{nm} = c_{n-m} - c_{n+m},   n, m = 1..n,     i.e.  M = T_n(f) - H_n(f).

*I initially assumed this was the Bini-Capovani tau algebra (Toeplitz-minus-Hankel
diagonalised by DST-I, eigenvalues exactly `f(j pi/(n+1))`) and **that is wrong** —
I tested it (`verify5.py`) and it fails from degree 2 onward. Truncation and
multiplication do not commute, so the sine-basis section of multiplication by `f` is
not `p(tau(2 cos theta))`; the two differ by a corner term (for `f = cos 2 theta` the
tau matrix carries `-1/2` at both corners, the section only at `(1,1)`). Recording the
error because the false version is tempting and would have made the law look trivial.*

What is true, and measured (`verify6.py`, `kR = 1`, coefficients from a 2^22-point FFT):
the Hankel part moves the `O(1/n)` correction but **not the leading constant**.

| n | `lambda_min(T - H) * n^2` (our matrix) | `lambda_min(T) * n^2` (pure Toeplitz) |
|---:|---:|---:|
| 40 | 0.398577 | 0.405558 |
| 80 | 0.404736 | 0.408288 |
| 160 | 0.407942 | 0.409730 |

both against `pi^2/24 = 0.411234`. So the KMS theorem is stated for `T_n(f)`, our
object is `T_n(f) - H_n(f)`, and the `pi^2` constant is common to both. Note the left
column reproduces the independent momentum-space build of `verify3.py` digit for digit
(0.398577 / 0.404736 / 0.407942) — two unrelated routes to the same matrix.

The relevant literature family for the Hankel correction is "Toeplitz plus Hankel" /
"flipped Toeplitz" spectral theory (e.g. Ekstrom & Serra-Capizzano, arXiv:2203.06992).
**I did not verify that any of it covers the extreme-eigenvalue-at-a-symbol-zero
question**, so for Paper 60 the honest line is: cite KMS for the constant, and report
the Hankel-part invariance as measured.

---

## Q1. Is this a named theorem, with a sharp constant? YES.

### The theorem

**Bottcher & Widom (2006), Eq. (1)** — read from the arXiv full text
(arXiv:math/0412269, `bw_raw.txt` lines 40-50):

> Suppose `a` is of the form `a(t) = |1-t|^{2 alpha} b(t)` (`t` in T) where `alpha`
> is a natural number and `b` is a positive function on T whose Fourier coefficients
> are subject to the condition `sum_k |k| |b_k| < infinity`. Then the matrix
> `T_n(a)` is positive definite and its smallest eigenvalue satisfies
>
>     lambda_min(T_n(a))  ~  (c_alpha / n^{2 alpha}) * b(1)   as n -> infinity
>
> with a certain constant `c_alpha` in (0, infinity) independent of `b`.
> **Kac, Murdock, and Szego [8] proved that `c_1 = pi^2`**, and Parter [11] showed
> that `c_2 = 500.5467`.

Later in the same paper (line 366): `c_1 = pi^2 = 9.8696, c_2 = 500.5467, c_3 = 61529`.

The constant is **sharp and explicitly characterised**: `c_alpha` is the minimal
eigenvalue of the boundary value problem `(-1)^alpha u^{(2 alpha)} = lambda u` on
[0,1] with `u = u' = ... = u^{(alpha-1)} = 0` at both ends. For `alpha = 1` that is
the Dirichlet Laplacian on the unit interval, whose minimal eigenvalue is `pi^2`.
Bottcher and Widom also show `c_alpha` is simultaneously the norm of a Green kernel,
the best constant in a higher-order Wirtinger-Sobolev inequality, and the conditioning
constant of a least-squares problem — five equivalent problems.

**Attribution** (Bottcher-Widom Section 2, line 212): *"For alpha = 1, formula (1)
goes back to Kac, Murdock, Szego [8]. In the late 1950's, Seymour Parter and the
second of the authors started tackling the general case, with Parter embarking on the
Toeplitz case and the second of us on the Wiener-Hopf case. In [11] (alpha = 2) and
then in [10], [12] (general alpha), Parter established (1)."*

### So `eq:sigma_law` is exactly `c_1 * b(1)`

    lambda_min ~ c_1 b(1) / n^2 = pi^2 * ((kR)^2/24) / n^2 = (kR)^2/24 * pi^2/n^2

which is `eq:sigma_law` **character for character**. It is a rediscovery.

### Numerical confirmation of both legs (done here)

`verify2.py`, canonical symbol `a = |1-t|^2 = 2 - 2 cos theta`, `b == 1`:

| n | lambda_min * n^2 | lambda_min * (n+1)^2 |
|---:|---:|---:|
| 40 | 9.389436 | 9.864776 |
| 160 | 9.747072 | 9.869291 |
| 640 | 9.838814 | 9.869585 |

against `pi^2 = 9.869604`. KMS `c_1 = pi^2` confirmed.

`verify3.py`, **our actual matrix** `C` (built in momentum space via `p = k cot(chi/2)`,
where the chirp is regular and `j0(kR p)` is slowly oscillating; 90 geometric panels,
120 Gauss points each, `p` from 1e-7 to 1e7):

| kR | n | `1 - sigma_max` | `x n^2/(kR)^2` | `x (n+1)^2/(kR)^2` |
|---:|---:|---:|---:|---:|
| 1.0 | 40 | 2.491107e-04 | 0.398577 | 0.418755 |
| 1.0 | 80 | 6.323993e-05 | 0.404736 | 0.414917 |
| 1.0 | 160 | 1.593523e-05 | 0.407942 | 0.413057 |
| 2.0 | 160 | 6.362191e-05 | 0.407180 | 0.412286 |

against `pi^2/24 = 0.411234`. Converges, from below in `n^2` and from above in
`(n+1)^2`.

### A free refinement Paper 60 can take

The paper reports the collapse "confirmed to ~1% at n=160". That residual is **not
noise, it is the known O(1/n) term**: the deficit is -0.80% at n = 160 and -1.58% at
n = 80, i.e. exactly halving. Fitting `(n+c)^{-2}` gives a stable `c`: at kR = 1,
`c = 0.63, 0.64, 0.65` for n = 40, 80, 160; at kR = 0.5, `c ~ 0.57`; at kR = 2,
`c ~ 0.80`. So `c` is symbol-dependent (not simply `n -> n+1`), which is precisely
what the modern simple-loop expansions supply in closed form:

- Bogoya, Bottcher, Grudsky & Maximenko (2015) give higher-order **uniform individual**
  asymptotics `lambda_j(T_n(a))` under the simple-loop (SL) condition — symbol has
  exactly two intervals of monotonicity, first derivative nonvanishing on them, second
  derivative nonzero at the min and max. The leading grid is `j pi/(n+1)` with computable
  corrections in powers of `1/(n+1)`.
- Ekstrom, Garoni & Serra-Capizzano (2018) turn the same expansion into the
  computational "almost closed form" `lambda_j = f(j pi/(n+1)) + sum_k c_k(theta_j) h^k`.

**Caveat on applying SL to our symbol:** our `a` has *infinitely many* intervals of
monotonicity (the `chi -> 0` chirp), so the SL condition **fails globally**. The
expansions are not literally licensed here. That is why the `c ~ 0.6-0.8` above is
reported as measured, not derived.

### The one hypothesis we genuinely fail (report this honestly)

Bottcher-Widom require `sum_k |k| |b_k| < infinity`. **Our symbol fails it.**
Measured (`verify4.py`, 2^21-point FFT): `|b_k| ~ k^(-1.246)` at kR = 1 and
`k^(-1.254)` at kR = 2 (and the symbol itself `|f_k| ~ k^(-1.24)`, consistent with
Paper 60's reported `j^(-1.3)`). Hence `sum k|b_k| ~ sum k^(-0.25)` **diverges**, and
visibly so — partial sums 9.07, 50.3, 283, 803 at K = 1e3, 1e4, 1e5, 4e5, growing
like `K^0.75`.

So: the *law and its constant* are KMS's; but no theorem I could find states the
**sharp constant** under hypotheses our symbol verifiably meets. What survives under
weak hypotheses is the **rate only** (next section). This is a real, small, honest
residue of originality — and it should be stated as "the KMS constant is observed to
survive at a merely continuous symbol with power-law Fourier decay", not as a derivation.

---

## Q2. How does the rate change?

### (b) Higher-order / flat maximum — fully known, with constants

Rate becomes `n^{-2 alpha}` for a zero of order `2 alpha`. Constants
(Bottcher-Widom, above): `c_1 = pi^2 = 9.8696`, `c_2 = 500.5467` (Parter),
`c_3 = 61529` (Bottcher-Widom, numerically). Their **main theorem** gives the
large-alpha asymptotics and two-sided bounds (arXiv:math/0412269, Eq. (9)-(10)):

    c_alpha = alpha sqrt(8) (4 alpha / e)^{2 alpha} (1 + O(1/alpha))   as alpha -> infinity

    (4 alpha - 2)/(4 alpha^2 - alpha) * (4 alpha)! (alpha!)^2 / [(2 alpha)!]^2
        <= c_alpha <=
    (4 alpha + 1)/(2 alpha + 1) * (4 alpha)! (alpha!)^2 / [(2 alpha)!]^2

The paper also explicitly kills the tempting conjecture `c_alpha ~ ((alpha+1)/2)^{2 alpha}`
(it matches to 3 figures at alpha = 3 by coincidence and is wrong).

**Non-integer order:** Rambour & Seghier (2012) extend Bottcher-Widom to
`phi_alpha(theta) = |1 - e^{i theta}|^{2 alpha} f_1(e^{i theta})` with
`alpha in (1/2, infinity) \ N`, obtaining `lambda_min ~ c_alpha N^{-2 alpha} f_1(1)`
with asymptotic expansions for the non-integer `c_alpha`. (Read from the numdam
landing page for the article; **abstract/summary level, I did not read the proof**.)

**General rate of attaining the extremum:** Novosel'tsev & Simonenko,
"Dependence of the asymptotics of extreme eigenvalues of truncated Toeplitz matrices
on the rate of attaining the extremum by the symbol", Algebra i Analiz 16 (2004)
146-152 (St. Petersburg Math. J. translation). **I have this by title and
bibliographic record only — I did not obtain the text and cannot state its theorem.**

**Flat in the strongest sense (extremum on a set of positive measure):** the decay
stops being polynomial and the relevant statement is a *lower* bound —
`lambda_min(T_n) >= exp(-c n)` with no smoothness assumption on `f`; and if
`a >= epsilon > 0` on some arc of T then `1/lambda_min = O(e^{alpha n})`. Reported by
Serra-Capizzano, "How bad can positive definite Toeplitz matrices be?" (2000).
**Obtained via search summaries of that paper and of citing works, not from the
primary text** — treat as directional, not quotable.

### (c) Symbol continuous but non-smooth elsewhere (our chi -> 0 chirp) — rate survives

This is the decisive point for us, and it is settled: **Serra-Capizzano (1998)**
proves the rate under `L^1` alone. From the paper's own abstract (retrieved via the
Semantic Scholar API for DOI 10.1016/S0024-3795(97)00231-0 — **abstract only, I did
not read the paper**): the rate of convergence of `lambda_0(n)` to `inf f` *"depends
only on the order rho"* of the zero, proved under weaker conditions than the prior
Kac-Murdock-Szego / Widom / Parter / R. H. Chan results, which had assumed `f` in
`C^k` at least locally; and the paper establishes lower bounds for minimal eigenvalues
of Toeplitz matrices generated by nonnegative `L^1` functions plus upper bounds on the
Euclidean condition number.

So our chirp — which lives at `chi -> 0`, i.e. as far from the extremum as the domain
allows — **does not change the exponent**. `n^{-2}` is not an artifact of our
particular symbol; it is forced by "nonnegative, single zero, order 2". It costs us
only the sharp-constant hypothesis, and the measurement above shows the sharp constant
survives anyway.

### (a) Extremum at an endpoint rather than an interior point — a non-question, twice over

- **In the Toeplitz setting there are no endpoints.** The symbol lives on the circle
  `T`; every point is interior, and `T_n(a(-t)) = D T_n(a) D` with `D = diag((-1)^j)`
  is unitarily equivalent, so `theta = 0` and `theta = pi` are interchangeable. The
  classical statements are normalised to a zero at `t = 1` with no loss.
- **In our sine-basis setting `chi = 0` and `chi = pi` are symmetric.** Under the even
  periodic extension both map to points of the circle (`theta = 0` and `theta = pi`),
  and the argument above applies to each. Our symbol does not attain its supremum at
  `chi -> 0` anyway: there `a` oscillates to zero with amplitude `~ chi/(2 kR)`, so
  `f = 1 - a -> 1`, bounded away from the zero. *(My own elementary reading, not a
  quoted theorem.)*

So `n^{-2}` here is not an endpoint artifact. The exponent is fixed by "nonnegative,
isolated zero, order 2" and nothing else.

A caution about the tempting inverse idea: one might hope a different transform whose
sampling grid **includes** the extremum point would make `lambda_min = min f` exactly
and kill the divergence. That is not available to us — the sine basis is forced by the
Fock geometry — and in any case the section of a multiplication operator is not
diagonalised by the grid (see the corrected note in Section 0, where exactly this
assumption was tested and failed).

---

## Q3. The prize: a structure-preserving escape

### Q3.0 — Measured first, on the real object: `S^{-1/2}` has NO effective locality

Before any literature: I built `S = [[I,C],[C^T,I]]` for our actual symbol, formed
`S^{-1/2}` by exact diagonalisation, and fitted the within-block off-diagonal decay
`|(S^{-1/2})_{r, r+d}| ~ exp(-d/L)` along a central row (`verify8.py`, `kR = 1`):

| n | cond(S) | decay length L | L/n | L/sqrt(cond) |
|---:|---:|---:|---:|---:|
| 20 | 2 065.27 | 2.83 | 0.142 | 0.0623 |
| 40 | 8 027.56 | 5.51 | 0.138 | 0.0615 |
| 80 | 31 624.59 | 10.88 | 0.136 | 0.0612 |
| 160 | 125 507.05 | 21.52 | 0.135 | 0.0607 |
| 320 | 500 027.54 | 42.78 | 0.134 | 0.0605 |

Two readings, both clean:

1. **`L` grows linearly in `n`** (`L/n -> 0.134`). The bandwidth needed to represent
   `S^{-1/2}` to any fixed accuracy is `Theta(n)`, so a banded truncation has
   `Theta(n^2)` nonzeros — **dense up to a constant. There is no banded approximation
   of `S^{-1/2}` at any fixed accuracy.**
2. **`L` is a constant times `sqrt(cond(S))`** (`L/sqrt(cond) = 0.0605-0.0623`, flat to
   3%). This is exactly the Demko-Moss-Smith-type scaling: the decay ratio
   `q = (sqrt(kappa)-1)/(sqrt(kappa)+1)` gives `q^d ~ exp(-2d/sqrt(kappa))`, i.e.
   decay length proportional to `sqrt(kappa)`.

As a side check, `cond(S) = 500 027` at `n = 320` against `2/(1-sigma_max) =
2/(0.411234/320^2) = 498 000` from `eq:sigma_law` — the conditioning law and this
table are the same object, confirmed independently.

**The structural conclusion, which is stronger than a literature answer:** the
conditioning and the non-locality are *one fact*, not two. The decay length of
`S^{-1/2}` is set by `sqrt(cond(S))`, which is set by `1 - sigma_max`, which is set by
the symbol's quadratic touch of its supremum. Any transformation that restores locality
must raise `1 - sigma_max`, i.e. must change the symbol — and the symbol is fixed by
the Fock geometry and `kR`. **A structure-preserving fix and a conditioning fix are the
same request, and it cannot be granted by a similarity transform.** This also
re-derives, from the metric side, why Lowdin orthogonalisation densifies: it is
`S^{-1/2}`, and `S^{-1/2}` is provably non-local here.

### Q3.0b — Measured second: a BANDWIDTH-1 preconditioner flattens the conditioning completely

This is the positive half, and it is large. Take the ill-conditioned sector
`M = ` the sine-basis section of `f = 1 - a` (the ungerade block `I - C`; Paper 60
already localises all the divergence there). Precondition with the section of the
**trigonometric polynomial** `g(theta) = 2 - 2 cos theta = 4 sin^2(theta/2)` — i.e. the
**tridiagonal** matrix, bandwidth 1 — which has *the same quadratic zero at
`theta = 0`* as `f`. Then (`verify9.py`, `kR = 1`):

| n | cond(M) | cond(G^{-1} M), G Toeplitz-Hankel | cond(G^{-1} M), G pure tridiagonal |
|---:|---:|---:|---:|
| 20 | 1 213.1 | 7.3425 | 7.3425 |
| 40 | 4 828.1 | 7.5760 | 7.5760 |
| 80 | 19 180.2 | 7.6470 | 7.6470 |
| 160 | 76 315.8 | 7.6677 | 7.6677 |
| 320 | 304 254.7 | 7.6733 | 7.6733 |
| 640 | 1 214 781.5 | 7.6747 | 7.6747 |

`cond(M)` grows as `n^2` (x1001 over a x32 range in `n`); **the preconditioned
condition number is FLAT — it converges to a constant.** And the constant is exactly
the one the theory predicts, `max h / min h` with `h = f/g` (`verify10.py`):

| kR | min h | max h | max/min | measured limit |
|---:|---:|---:|---:|---:|
| 0.5 | 0.010289 | 0.308102 | 29.945 | — |
| 1.0 | 0.041551 | 0.319801 | **7.6965** | **7.6747 at n=640, rising** |
| 2.0 | 0.166601 | 0.371632 | 2.231 | — |

with `h(0+) = (kR)^2/24` (our `b(1)`) and `h(pi) = 1/4` recovered analytically.
This is **Serra-Capizzano's band-Toeplitz preconditioning working exactly as
advertised**: `f = g h` with `g` a trigonometric polynomial carrying the zero and `h`
continuous and bounded away from 0 and infinity, and the preconditioned spectrum is
trapped in `[min h, max h]`. Note our chirp does not obstruct it — `h` stays positive
and bounded because `f -> 1` and `g -> 4` at the chirp end.

### Q3.0c — But read what the preconditioner actually does: it MOVES the problem

The preconditioner `G` is the discrete Dirichlet Laplacian, and its *own* condition
number is `4n^2/pi^2` (`verify10.py`: n=320 gives 41 760 against `4n^2/pi^2 = 41 501`).
So `G^{-1}M` being `O(1)` does not mean the ill-conditioning was removed — it means

> **our metric's ill-conditioning is spectrally equivalent to the 1-D discrete
> Dirichlet Laplacian's**, with a bounded, `kR`-dependent distortion.

That is the sharpest structural statement available here, and it is worth more to
Paper 60 than a "fix". Consequences, stated honestly:

- **Classically the problem is solved.** Solves with `S` need `O(1)` PCG iterations
  with a *tridiagonal* preconditioner (an `O(n)` direct solve), or standard multigrid.
  The metric is no longer an obstacle to anything classical.
- **For the quantum block-encoding it is NOT a fix, and Paper 60 should not claim it
  is.** The deliverable there is `S^{-1/2}`, and Q3.0 shows `S^{-1/2}` is provably
  non-local (decay length `~ 0.06 sqrt(cond) ~ 0.134 n`). A preconditioner accelerates
  *iterative solves*; it does not hand you a sparse `S^{-1/2}`. Whether the
  spectral-equivalence-to-the-Laplacian can be exploited in a QSVT circuit (e.g.
  block-encode `G` and `H = G^{-1/2} S G^{-1/2}` separately) is an **open question I
  did not resolve** — `G^{-1/2}` is itself non-local, so the naive route does not close.

### Q3.0c2 — Circulant preconditioners FAIL here, and for the textbook reason

The complementary control (`verify12.py`, pure Toeplitz section `T_n(f)`, `kR = 1`):

| n | cond(T_n(f)) | cond(Strang^{-1}T) | cond(TChan^{-1}T) | min eig of Strang circulant |
|---:|---:|---:|---:|---:|
| 16 | 682.6 | 92.27 | 51.27 | +1.574e-03 |
| 32 | 2 925.6 | **indefinite** | 138.92 | -1.370e-03 |
| 64 | 12 016.1 | 152.09 | 386.26 | +1.012e-04 |
| 128 | 48 466.5 | **indefinite** | 1 087.18 | -1.985e-06 |
| 256 | 194 160.9 | **indefinite** | 3 057.50 | -2.933e-06 |

- **Strang's circulant goes indefinite** (negative eigenvalues) from `n = 32` — unusable.
- **T. Chan's optimal circulant stays positive definite but does not cluster**: the
  preconditioned condition number *grows* like `n^{1.5}` (ratio ~2.8 per doubling).
  It buys a constant (194 161 -> 3 058 at n=256) and then loses ground.

This is exactly the known dichotomy: circulant preconditioners deliver spectral
clustering only for symbols **bounded away from zero**; a symbol with a zero defeats
them, while a **banded** preconditioner carrying the same zero succeeds (Q3.0b). So the
answer to "is there a *structured* fix" is **banded yes, circulant no** — and that is
a literature-predicted split, reproduced here on the actual object.

### Q3.0d — The Cholesky lever is closed too

The generalized problem `W c = k S c` does not actually require `S^{-1/2}` — any
factorization `S = L L^T` works (`L^{-1} W L^{-T} y = k y`). So: is the Cholesky factor
local even though `S^{-1/2}` is not? Measured (`verify11.py`, decay of `L^{-1}` along a
row inside the second block):

| n | decay length of `L^{-1}` | /n |
|---:|---:|---:|
| 20 | 5.28 | 0.264 |
| 40 | 8.94 | 0.224 |
| 80 | 16.02 | 0.200 |
| 160 | 29.83 | 0.186 |
| 320 | 57.03 | 0.178 |

`L^{-1}`'s decay length grows like `n^{0.86}` over this window and the ratio `/n` is
still drifting slowly downward — but it is nowhere near bounded. **No fixed bandwidth
suffices for the Cholesky route either.** Both whitening factorizations are non-local;
the non-locality is a property of `S`, not of the choice of square root.

### Q3.1 — The literature side: the positive half is a PROVEN CLASSICAL THEOREM

**S. Serra, "Optimal, quasi-optimal and superlinear band-Toeplitz preconditioners for
asymptotically ill-conditioned positive definite Toeplitz systems", *Math. Comp.* **66**
(218) (1997) 651-665**, AMS id S 0025-5718(97)00833-8. **[PRIMARY, full text read]**

Verbatim from the introduction:

> "When f has zeros, i.e., ess inf f = 0 we know [20] that the Euclidean condition
> number mu_2(A_n(f)) ... grows to infinity ... **In the case where r = 2k is an even
> number only tau preconditioners [15] and band-Toeplitz preconditioners [7, 16] are
> shown to be able to reduce the condition number from O(n^{2k}) to O(1).**"
>
> "The main idea (see [16]) is to find a trigonometric polynomial g for which
> r < f/g < R where r, R are positive constants. The associated band-Toeplitz matrix
> A_n(g) results to be the desired preconditioner in the sense that **the spectrum of
> A_n^{-1}(g) A_n(f) lies in (r, R) for any dimension n.**"

**Theorem 2.2** (Serra's statement, attributed to Di Benedetto-Fiorentino-Serra 1993
Thms 3.1/3.2): for `f, g` essentially nonnegative in `L^1`, the eigenvalues of
`A_n(g)^{-1} A_n(f)` lie in `(r, R)` with `r = ess inf f/g`, `R = ess sup f/g`; their
union over `n` is **dense** in the essential range of `f/g`; and `lim lambda^n_1 = r`,
`lim lambda^n_n = R`.

**Theorem 4.1** gives the explicit minimal preconditioner for `f` nonnegative continuous
**with zeros of even order**: `z_k = prod_i (2 - 2 cos(x - x_i))^{l_i}` over the zeros
`x_i` of order `2 l_i`.

**Our symbol is squarely inside the hypothesis class** — this is worth stating flatly,
because I had expected otherwise from the Q1 analysis:
- `f = 1 - a` nonnegative: yes (`|j0| <= 1`, equality only at argument 0);
- continuous on the closed interval: yes (the chirp has amplitude `-> 0`, so `a(0+) = 0`);
- exactly one zero, of order exactly **2**: yes;
- **the sign change of `j0` is irrelevant** — it lives in `I + C`, where it is bounded
  and harmless (see Q3.3);
- **the power-law `j^{-1.3}` coefficient decay violates nothing.** Serra needs `f`
  *continuous*; smoothness enters only in Thm 4.2 (a Jackson estimate) to bound how fast
  the preconditioner improves with *bandwidth*. At minimal bandwidth no smoothness is used.

So the hypothesis gap that bites the *constant* in Q1 does **not** bite the
*preconditioning* result. Different theorem, weaker hypotheses, and we satisfy them.

**Even better, we do not need the Toeplitz machinery at all.** Because our matrix is the
moment matrix of a multiplication operator, positivity is immediate and elementary:

> For any `f >= 0`, `v^T M_n(f) v = (2/pi) integral_0^pi f(chi) |sum_j v_j sin(j chi)|^2 dchi >= 0`,
> so `M_n` is an order-preserving positive map. Hence `r <= f/g <= R` pointwise with
> `g >= 0` gives `r M_n(g) <= M_n(f) <= R M_n(g)`, so
> **`spec(M_n(g)^{-1} M_n(f)) subset [r, R]` for every `n`** — no Szego theory, and it
> applies verbatim to *any* Gram/moment matrix in *any* orthonormal basis.

*(This two-line argument is the Q3 scan's, not a quotation; I checked it and it is
correct. It is the cleanest way to state the result in Paper 60, and it sidesteps the
Toeplitz-vs-Toeplitz-minus-Hankel issue entirely.)*

Supporting primary sources:
- **R. H. Chan, "Toeplitz preconditioners for Toeplitz systems with nonnegative
  generating functions", *IMA J. Numer. Anal.* **11**(3) (1991) 333-345,
  DOI 10.1093/imanum/11.3.333** — the original banded-preconditioner-for-zeros theorem;
  band-widths independent of `n`, spectra of `C_n^{-1}T_n` uniformly bounded.
  *[PRIMARY, abstract only.]*
- **S. Hon, S. Serra-Capizzano, A. Wathen, "Band-Toeplitz preconditioners for
  ill-conditioned Toeplitz systems", *BIT* **62**(2) (2022) 465-491,
  DOI 10.1007/s10543-021-00889-6** *[PRIMARY, full text read]* — restates Serra Thm 4.1
  verbatim and confirms the history: *"While most circulant preconditioners do not work
  well, ... R. Chan [7] and Di Benedetto et al. [13] proposed band-Toeplitz
  preconditioners ... constructed by using certain polynomial that matches the zeros of f."*
- **F. Di Benedetto, G. Fiorentino, S. Serra, "C.G. preconditioning for Toeplitz
  matrices", *Comput. Math. Appl.* **25** (1993) 35-45** — origin of Thm 2.2.
  *[SECONDARY — statement obtained only through Serra 1997 and Hon et al.]*
- **F. Di Benedetto, "Analysis of preconditioning techniques for ill-conditioned
  Toeplitz matrices", *SIAM J. Sci. Comput.* **16** (1995) 682-697** — the **tau/DST**
  version: preconditioners `S_n Lambda S_n` with `S_n` the discrete sine transform,
  bounded condition number given only the position and order of the zeros.
  ***[SECONDARY, via search index — NOT read. This is the single most important item
  for anyone to verify, since it is the DST result matching our basis.]***
- Multigrid alternative, also `O(n)` and also handling zeros: **S. Serra Capizzano,
  "Convergence analysis of two-grid methods for elliptic Toeplitz and PDEs
  matrix-sequences", *Numer. Math.* **92** (2002) 433-465, DOI 10.1007/s002110100331**
  — convergence "under the sole assumption that the generating function is nonnegative
  and has a zero at x0 = 0 of finite order". *[abstract only.]* Our order-2 zero is
  inside even the earlier, weaker versions.

### Q3.2 — The negative half: decay for `S^{-1/2}` is FORBIDDEN, not merely unproven

Every decay theorem excludes our case **by hypothesis**, and the exclusion is the same
one each time: the operator must be **boundedly invertible** / `0 notin spectrum`.

- **Demko, Moss & Smith, "Decay rates for inverses of band matrices", *Math. Comp.*
  **43**(168) (1984) 491-499** *[PRIMARY, full text read]*. Thm 2.4:
  `|A^{-1}(i,j)| <= C lambda^{|i-j|}` with
  `lambda = [(sqrt(cond A) - 1)/(sqrt(cond A) + 1)]^{2/m}`,
  `C = ‖A^{-1}‖ max{1, (1 + sqrt(cond A))^2/(2 cond A)}`.
  Fails on **two** counts for us: `A` must be **banded** (ours is power-law, `r ~ 1.3`),
  and boundedly invertible (ours is not, in the limit).
- **Benzi & Golub, "Bounds for the entries of matrix functions with applications to
  preconditioning", *BIT* **39**(3) (1999) 417-438, DOI 10.1023/A:1022362401426**
  *[PRIMARY, full text read]* — and it treats `A^{-1/2}` **explicitly** (p. 421), with
  `chi_bar = (sqrt(kappa)+1)^2/(kappa-1)`. Their own verbatim caveats:
  > "**This constant can be very large for a near zero, as is the case if A is nearly
  > singular.**" ... "**Conversely, the same formulas show that decay can be arbitrarily
  > slow as the condition number of A increases to infinity, as in this case q -> 1^-
  > and K_0 grows without bound.**"
  and p. 423:
  > "**because we are dealing here with matrices of finite order n, the decay rate may be
  > so slow that no entry in f(A) is actually small, particularly in the ill-conditioned
  > case.**"
  Note `1/chi_bar` inverts to `(sqrt(kappa)-1)/(sqrt(kappa)+1)` — **Benzi-Golub's
  limiting `q` for `x^{-1/2}` is identical to Demko's `q` for `x^{-1}`.** Both degrade at
  the same rate.
- **Groechenig's functional calculus is the decisive exclusion.** From Groechenig,
  "Wiener's Lemma: Theme and Variations" (lecture notes; published in *Four Short
  Courses on Harmonic Analysis*, Birkhauser 2010, DOI 10.1007/978-0-8176-4891-6_5)
  *[PRIMARY, read]*: the decay-preserving calculus is the **holomorphic** (Riesz) one,
  requiring **`f` analytic on an open neighbourhood of the spectrum**, and its corollary
  "⟹ square roots, powers, pseudoinverse". `x^{-1/2}` is singular at 0; our limit
  operator's spectrum **contains** 0. **Excluded hypothesis, not an open case.**
  Same exclusion in **Groechenig & Klotz, *Constr. Approx.* **32** (2010) 429-466,
  arXiv:0904.0386** *[PRIMARY, read]* — every theorem carries "If A is invertible on
  ell^2(Z)"; and in **Q. Sun, *Constr. Approx.* **33** (2011) 317-342, arXiv:1001.1457**
  (ell^q_w-stability = uniform invertibility).
- **Jaffard (1990)**, *Ann. IHP Anal. Non Lineaire* **7**(5) 461-476: polynomial
  off-diagonal decay is inverse-closed — **but only for invertible `A`**. *[Statement via
  Groechenig-Klotz section 1, NOT from Jaffard itself.]* And even if it applied,
  `r = 1.3` is barely above the 1-D threshold `r > 1`: `sum_{d>b} d^{-1.3} ~ b^{-0.3}/0.3`,
  so 99% of the `ell^1` mass needs `b ~ 2e8`. **Jaffard-class locality at exponent 1.3 is
  theoretical, not practical.**
- **Groechenig & Leinert, *Trans. AMS* **358**(6) (2006) 2695-2711** — *[BIBLIOGRAPHIC
  RECORD ONLY, not read; no theorem numbers claimed.]*
- **No published theorem of the form "f(A) banded for all banded A ⟹ f polynomial" was
  found. Do not attribute one.** But an elementary Fejer-Riesz argument settles it for
  multiplication operators: `A^{-1/2}` is banded iff `f^{-1/2}` is a trigonometric
  polynomial, so for nonconstant nonnegative `f` with a zero, **never**. For us,
  `(1-a)^{-1/2} ~ c/|chi - pi|`, which is **not in `L^1`** — the limiting operator's
  entries do not exist as an absolutely summable family. *(Argument is the Q3 scan's,
  not a citation.)*

**Consistency with my own measurement (Q3.0):** the scan measured `(I-C)^{-1/2}` entries
at **fixed** offsets and found they *grow like log n* (increments ~+0.539 per doubling),
with a 99%-mass bandwidth that is a fixed **fraction** (0.44) of `n`. That is the same
phenomenon as my "decay length `= 0.134 n`", from the other side, and it is the stronger
statement. **Lowdin whitening of this metric is not expensive — it is non-local by
construction, at every n.**

**The arithmetic requested (how the bounds degrade as `kappa ~ n^2`):** with
`kappa ~ 48 n^2/((kR)^2 pi^2)`, `sqrt(kappa) ~ 6.93 n/(kR)`:
`q = 1 - 2/sqrt(kappa)`, `lambda = q^{2/m} ~ exp(-(4/m)/sqrt(kappa))`, so the
**decay length is `m sqrt(kappa)/4 ∝ n`** — i.e. `exp(-c|i-j|/n)`, no effective locality
— while the prefactor `C = ‖A^{-1}‖ ~ n^2` diverges. At `kR = 2, n = 512, m = 2`:
`kappa = 3.20e5`, `q = 0.996464`, `C = 1.60e5`, so the bound drops below 1 only at
`|i-j| >= 3386` — **6.6x the matrix dimension. Vacuous for every entry.**

### Q3.3 — The constructive result (verified independently by me)

The scan proposed, and I reproduced from scratch (`verify13.py`), a **structure-preserving
whitening**. Lowdin's symmetric `S^{-1/2}` is only one of many whitenings, and it is the
non-local one. Instead put

    P := tridiagonal preconditioner (symbol 2 - 2 cos theta, the zero-matching polynomial)
    G := P^{-1/2} (I - C) P^{-1/2},        W := G^{-1/2} P^{-1/2},   so  W (I-C) W^T = I.

Measured at `kR = 2` (my run; the scan's independent numbers in parentheses):

| n | cond(G) | G_rr | d=1 | d=2 | d=4 | d=8 | d=16 | d=32 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 128 | 2.2286 (2.2286) | 0.2503 | 3.411e-2 | 2.041e-2 | 5.398e-3 | 1.906e-3 | 1.325e-3 | 8.136e-4 |
| 256 | 2.2295 (2.2295) | 0.2500 | 3.382e-2 | 2.067e-2 | 5.573e-3 | 1.990e-3 | 1.579e-3 | 6.356e-4 |
| 512 | 2.2297 (2.2297) | 0.2500 | 3.382e-2 | 2.069e-2 | 5.608e-3 | 1.919e-3 | 1.499e-3 | 6.924e-4 |

and `G^{-1/2}` inherits it — diagonal `2.0488, 2.0500, 2.0500`; `d=1: 0.1495, 0.1483,
0.1482`; `d=4: 1.191e-2, 1.265e-2, 1.280e-2`. **Stable in `n` to 3 digits.**

So: **the entire non-locality of Lowdin whitening for this metric is the inverse square
root of a tridiagonal discrete Laplacian**, and `P` — unlike `I - C` — really *is* in the
Bini-Capovani tau algebra, so `P^{-1/2}` is **exactly diagonalized by the DST**:
`P^{-1/2} = S diag((2 - 2cos(j pi/(n+1)))^{-1/2}) S`. Classically `O(n log n)`; as a
circuit it is a Fourier-type primitive. Everything else (`G`, `G^{-1/2}`) is uniformly
local with `n`-independent decay — and now **inside** the Benzi-Golub/Jaffard hypotheses,
because `cond(G) = 2.23` is bounded.

**Caveat, stated plainly:** the uniform locality of `G` and `G^{-1/2}` is **measured, not
proved**. The literature supplies the `O(1)`-conditioning theorem (Serra 1997) and the
decay theorems that *would* apply to `G` given bounded conditioning; nobody has proved
`G` lies in a decay class with `n`-independent constants. Do not state it as a theorem.

### Q3.4 — Obstruction list

1. **Serra 1997 §3 (primary, verbatim):** band preconditioners reduce the condition
   number to `O(1)` *"only when the zeros of f have even order"*; for odd/fractional
   order the preconditioned spectra **cannot** be contained in a positive interval.
   **Does not bite us (order 2) — but it is the exact boundary.**
2. **Serra 1997 Thm 2.3(2) (primary):** the band-preconditioned spectrum is **dense** in
   `[1/(1+h), 1/(1-h)]` — *"we have no clusters ... but practically a uniform
   distribution of the spectrum."* **Band preconditioners are optimal (`O(1)` iterations)
   but never superlinear.** Do not expect clustering.
3. **Serra Capizzano & Tyrtyshnikov, "Any circulant-like preconditioner for multilevel
   matrices is not superlinear", *SIAM J. Matrix Anal. Appl.* **21**(2) (1999) 431-439;
   and "How to prove that a preconditioner cannot be superlinear", *Math. Comp.* **72**
   (243) (2003) 1305-1316, DOI 10.1090/S0025-5718-03-01506-0** *[abstracts/records only]*.
   **These are MULTILEVEL results and do NOT bar the present 1-level problem — do not
   cite them as such.** But they bite the moment the construction goes multi-index
   (three-centre / polyatomic), where circulant-, tau-, Hartley- and DCT-algebra
   preconditioners are all provably non-superlinear. **That is the real warning for the
   next step of the program.**
4. **Groechenig/Jaffard/Sun:** decay-preserving functional calculus is holomorphic and
   requires `0 notin spectrum`. Excluded hypothesis.
5. **Benzi-Golub p. 423 (primary):** finite-`n` decay bounds can be so slow that no
   entry is small — with a worked tridiagonal example whose inverse has *no* decay.
6. **Fejer-Riesz:** `A^{-1/2}` banded ⟹ `A = cI`, for multiplication operators.

### Q3.5 — A numerical-hygiene point that bites Paper 60 directly

`cond(S) = (1 + sigma_max)/(1 - sigma_max)` is a **catastrophic-cancellation formula**.
At `n = 512`, `1 - sigma_max ~ 6e-6`, so forming `sigma_max` in double precision and
subtracting loses ~6 digits; by `n ~ 1e4` it loses all of them. **Never compute
`1 - sigma_max` from `sigma_max`.** Get it as the smallest eigenvalue of `I - C`
(backward stable — this is what all the probes above do), or use the sine-based
principal-angle formulation: **Bjorck & Golub, *Math. Comp.* **27**(123) (1973) 579-594**
*[record only, not read]*; **Knyazev & Argentati, *SIAM J. Sci. Comput.* **23**(6) (2002)
2008-2040** *[abstract only]* — *"All popular software codes for canonical correlations
compute only cosine of principal angles, thus making impossible, because of round-off
errors, finding small angles accurately."*; **Drmac, *SIAM J. Matrix Anal. Appl.* **22**(1)
(2000) 173-194** *[abstract only]* (mixed stability of Bjorck-Golub).

### Q3.6 — Chan & Ng survey: NOT verified

**R. H. Chan, M. K. Ng, "Conjugate gradient methods for Toeplitz systems", *SIAM Review*
**38**(3) (1996) 427-482, DOI 10.1137/S0036144594276474** — **[BIBLIOGRAPHIC RECORD
ONLY.]** Full text inaccessible (SIAM 403; the author's own copy also 403). **Nothing is
claimed about its contents.** The "circulant clustering requires a positive symbol"
statement is attested instead by Serra 1997 §1 and Hon et al. §1, both read in primary.

### Q3.7 — The two pathologies are independent (important for the write-up)

- The **chirp** at `chi -> 0` (`j^{-1.3}` coefficients) causes **non-locality** but *no*
  ill-conditioning (`a -> 0` there, so `1 -/+ a -> 1`).
- The **quadratic touch** at `chi = pi` causes **all** the ill-conditioning but *no*
  non-locality (it is exactly matched by a tridiagonal).

Fixing one does not fix the other, and only the second has a structure-preserving fix.
**Any claim of the form "the metric is badly conditioned *because* the basis is
overcomplete / non-smooth" conflates them** — Paper 60's §sec:obstruction framing should
be checked against this.

Corollary, measured: the `I + C` block needs nothing at all. `1 + a` ranges in
`[1 - 0.21723, 2]`, so `cond(I+C) <= 2/(1 - 0.21723) = 2.5549` — exactly Paper 60's own
gerade constant `2.555041...`, here re-derived as a *sign-change* effect of `j0` rather
than a symmetry effect. (And applying `P` to the `+` block is a blunder: it manufactures
a spurious zero; the scan measured `cond(P^{-1}(I+C))` growing 938 -> 229 459.)

---

## Q4. Prior appearance of this symbol / AO overlap as a Toeplitz section

### (A) The symbol `j0(kR cot(chi/2))` itself — NO, but the *operator* is prior art

Nothing returns on "Fock projection" + Toeplitz, or Sturmian + Toeplitz. No source
writes the two-centre Sturmian overlap as a multiplication operator with a named symbol.

But the operator is old and named. **The two-centre Coulomb-Sturmian overlap in the
`1/r`-weighted metric IS the Shibuya-Wulfman integral**, and Avery's own framing is
that the Shibuya-Wulfman expansion expands a plane wave in four-dimensional spherical
harmonics — i.e. the kernel is `e^{i p . R}` on S^3, which makes the SW matrix literally
the matrix of multiplication by `e^{i p . R}` in the hyperspherical-harmonic basis.
The `l = 0` reduction to `j0(pR)` is the elementary angular average. *(Avery's IJQC 2004
and the books were inaccessible — 403 — so this is abstract + concordant secondary
summaries, not a verbatim equation. Flag as secondary-source.)*

So: **the operator is prior art; the operator-theoretic reading of it — symbol, finite
section, extreme-eigenvalue asymptotics — is not.** That is the correct seam for the
novelty claim.

Two adjacent items:
- Red & Weatherford's general SW formula states the matrix depends on *the translation
  distance times the screening parameter*, i.e. on `kR` alone — matching our symbol's
  single parameter, though not the symbol.
- Meremianin (arXiv:math-ph/0510080, Eq. 41) has the 4-D plane-wave to
  Gegenbauer-times-Bessel expansion, explicitly motivated by Sturmian bases in
  many-centre Coulomb problems, but does **not** connect it to Shibuya-Wulfman or to
  two-centre overlaps *(verified from the paper)*.

**Convention warning (a referee will catch this).** Every source checked writes Fock's
projection as `p = p_0 tan(chi/2)` (Meremianin & Rost Eq. 3; Elesin, Podlivaev & Openov,
arXiv:quant-ph/0403195, `alpha = 2 arctan(|p|/|p_0|)`; Aquilanti). Paper 60 uses
`p = k cot(chi/2)` — the antipodal chart, which is what moves the symbol's maximum to
`chi = pi`. **State the convention explicitly** or the formula will read as wrong.

### (B) AO overlap as a Toeplitz finite section — PARTIAL, and thin

**No one has done this for a two-centre overlap.** Two precedents:

1. A **single-centre** even-tempered AO overlap matrix has been called Toeplitz once:
   Wang, Dowdle & Whitfield, *APL Comput. Phys.* **2** (2026) 016106, state the
   normalised even-tempered overlap "forms a Toeplitz matrix only dependent on beta
   and N". **Caveat: that sentence is from a search snippet of the published AIP PDF
   (403); the arXiv v1 (2511.03579), which was read in full, contains no occurrence of
   "Toeplitz"** and reports condition-number growth purely empirically — no symbol, no
   Szego, no finite-section theory. The word appears; the analysis does not.
2. The `n^2` law is Parter-Widom, as established in Q1/Q2 above. This was found
   **independently** by the Q4 scan, which is a useful cross-check on the Q1 verdict.

Everything else in the QC literature is non-symbolic — verified: Lehtola 2019
(arXiv:1911.10372, read in full: purely algorithmic pivoted Cholesky, no scaling law,
no Toeplitz/Szego/principal angles); Aissing & Monkhorst 1992 (linear dependence as
intrinsic to 3-D extended systems — a dimensional argument, abstract-level);
Hoyvik 2020 (the one QC paper on the overlap *spectrum* as such — abstract only,
Taylor & Francis 403); Klahn & Bingel 1977 / Klahn 1981 (Gram determinants, no
asymptotics); Lowdin's orthogonalisations (remedy, not analysis).

**On principal angles for AO overlaps:** the corresponding-orbital machinery
(Amos-Hall, King et al.) *is* the SVD/principal-angle construction, but it is applied
to the orbital sets of **two determinants**, not to the **two nuclear centres' basis
subspaces**. No paper was found that writes the diatomic AO overlap as `[[I,C],[C^T,I]]`,
identifies `sigma_k` as canonical correlations between per-centre subspaces, and reads
off `cond = (1+sigma_max)/(1-sigma_max)`. The identity is elementary and standard in
two-subspace theory; **this application appears to be ours.**

### (C) A "Szego theorem for principal angles" — PARTIAL; a mature neighbouring literature

The angles between the **past and future** of a stationary process are governed by the
spectral density — a symbol. This is the Helson-Szego / "angle between past and future"
literature (verified from N. H. Bingham's open survey *"Szego's theorem and its
probabilistic descendants"*, section 6.2): Helson & Szego 1960 (positive angle iff
`log w = u + v~` with `‖v‖_inf < pi/2`, equivalently Muckenhoupt A_2); Helson & Sarason
1967; Jewell & Bloomfield 1983 (canonical correlations of past and future, explicitly);
Peller & Khrushchev 1982 and Pavon (canonical correlations = **Hankel** singular values).

**But it is not our theorem.** That literature is Hankel, and its asymptotic parameter
is the **time lag**, not the truncation order `n`. Nothing was found of the form
"asymptotic distribution of the principal angles between two per-centre finite sections,
governed by a symbol, as `n -> infinity`". The nearest generic result is random rather
than symbol-driven (Aubrun 2022: principal angles between random half-dimensional
subspaces are asymptotically uniform on `[0, pi/2]`).

**Position this as a transfer from prediction theory, not an invention.**

### (D) Citation defects found in Paper 60's existing five

| key | real? | says what Paper 60 claims? |
|---|---|---|
| `amos_hall1961` | YES | YES — corresponding orbitals; proves diagonalisation of a rectangular matrix by two unitaries (= the SVD). *Abstract-level; full text 403.* |
| `king1967` | YES | YES — "for any two sets of spin orbitals there exist equivalent sets with a diagonal overlap matrix". *Abstract-level.* |
| `west_ruedenberg2013` | YES | **PARTIALLY — over-characterised.** Uses SVD to maximise overlap between MO and free-atom subspaces, which *is* a principal-angle computation in substance, but no evidence it uses "corresponding orbitals", "principal angles" or "canonical correlations". Soften to "cf. the SVD-based subspace-alignment construction of…". *Full text inaccessible (AIP + ISU repo 403) — UNRESOLVED.* |
| `halmos1969` | YES — **primary PDF read** | **PARTIALLY. Halmos does NOT state `‖[P,Q]‖`.** Theorem 2 gives the canonical form: for `M, N` in generic position with projections `P, Q`, there are positive contractions `S, C` with `S^2 + C^2 = 1`, `ker S = ker C = 0`, such that `P ~ [[1,0],[0,0]]` and `Q ~ [[C^2, CS],[CS, S^2]]`. The commutator norm is a one-line corollary: `PQ - QP = [[0, CS],[-CS, 0]]`, so `‖[P,Q]‖ = ‖CS‖ = max sigma sqrt(1-sigma^2) <= 1/2`. **Cite Halmos for the canonical form and derive the commutator norm in one line — do not attribute it to him.** |
| `bottcher_spitkovsky2010` | YES (LAA **432**(6) (2010) 1412-1459, DOI 10.1016/j.laa.2009.11.002) | **UNVERIFIED** that it states `‖[P,Q]‖ = max sigma sqrt(1-sigma^2)` or the `1/2` ceiling — paywalled everywhere, no preprint. **Do not attribute the identity to it without opening the PDF.** |

**The identity itself is correct** — I checked it directly (`verify7.py`): over 200
random subspace pairs in dimensions 4-12, `max | ‖[P_A,P_B]‖ - max_k sigma_k
sqrt(1-sigma_k^2) | = 1.3e-14`, and `max_{s in [0,1]} s sqrt(1-s^2) = 1/2` at
`s = 1/sqrt(2)`. So only the *attribution* is defective, not the mathematics.

**Verifiable substitute:** T. A. Loring, "Principal angles and approximation for
quaternionic projections", *Ann. Funct. Anal.* **5** (2014) 176-187,
DOI 10.15352/afa/1396833512 — the proof of Theorem 3.3 states
`‖P_theta Q_theta - Q_theta P_theta‖ = (1/2) sin(2 theta)` for the 2x2 blocks.
Halmos Thm 2 (direct sum over blocks) + Loring gives the full statement with two
citations that can both be stood behind.

### (E) One more bibliographic fix

Shibuya & Wulfman, *Proc. R. Soc. Lond. A* **286** (1965) is pages **376-389**;
Paper 60 currently gives only p. 376.

---

## Citations (Q1/Q2), with what each actually proves

1. **M. Kac, W. L. Murdock, G. Szego**, "On the eigenvalues of certain Hermitian
   forms", *J. Rational Mech. Anal.* **2** (1953) 767-800.
   JSTOR: http://www.jstor.org/stable/24900353
   — **PROVES `c_1 = pi^2`**, i.e. `lambda_min(T_n(a)) ~ pi^2 b(1)/n^2` for a
   second-order zero. *Attribution taken from Bottcher-Widom Section 2, which I read
   in full text; I did not read KMS 1953 itself.* This is the paper `eq:sigma_law`
   must cite.

2. **A. Bottcher and H. Widom**, "From Toeplitz eigenvalues through Green's kernels to
   higher-order Wirtinger-Sobolev inequalities", *Operator Theory: Advances and
   Applications* **171** (2006) 73-87. Preprint arXiv:math/0412269.
   DOI: 10.1007/978-3-7643-7980-3_4
   — **READ IN FULL TEXT.** States Eq. (1) with hypotheses; records `c_1 = pi^2`,
   `c_2 = 500.5467`, `c_3 = 61529`; **proves** the large-alpha asymptotics
   `c_alpha ~ alpha sqrt(8) (4 alpha/e)^{2 alpha}` and the two-sided bounds (10);
   proves the equivalence of the five problems. Best single citation: it carries both
   the statement and the attribution.

3. **S. V. Parter**, "On the extreme eigenvalues of truncated Toeplitz matrices",
   *Bull. Amer. Math. Soc.* **67** (1961) 191-196.
   **S. V. Parter**, "Extreme eigenvalues of Toeplitz forms and applications to
   elliptic difference equations", *Trans. Amer. Math. Soc.* **99** (1961) 153-192.
   **S. V. Parter**, "On the extreme eigenvalues of Toeplitz matrices",
   *Trans. Amer. Math. Soc.* **100** (1961) 263-276. DOI 10.1090/S0002-9947-1961-0138981-6
   — general `alpha`; `c_2 = 500.5467`. *Bibliographic data and role taken from the
   Bottcher-Widom reference list, which I read; I did not read Parter.*

4. **H. Widom**, "On the eigenvalues of certain Hermitian operators",
   *Trans. Amer. Math. Soc.* **88** (1958) 491-522.
   **H. Widom**, "Extreme eigenvalues of translation kernels", *Trans. Amer. Math. Soc.*
   **100** (1961) 252-262.
   **H. Widom**, "Extreme eigenvalues of N-dimensional convolution operators",
   *Trans. Amer. Math. Soc.* **106** (1963) 391-414.
   — the Wiener-Hopf/continuous side of the same programme. *Bibliographic data from
   the Bottcher-Widom reference list.*

5. **S. Serra Capizzano**, "On the extreme eigenvalues of Hermitian (block) Toeplitz
   matrices", *Linear Algebra Appl.* **270** (1998) 109-129.
   DOI: 10.1016/S0024-3795(97)00231-0
   — the rate under `L^1` only, no smoothness. *Abstract only (Semantic Scholar API).*
   This is the citation that licenses the **exponent** for our chirpy symbol.

6. **S. Serra Capizzano**, "On the extreme spectral properties of Toeplitz matrices
   generated by L^1 functions with several minima/maxima", *BIT* **36** (1996) 135-142.
   DOI: 10.1007/BF01740550 — *title/venue verified via the Springer landing page; not read.*

7. **J. M. Bogoya, A. Bottcher, S. M. Grudsky, E. A. Maximenko**, "Eigenvalues of
   Hermitian Toeplitz matrices with smooth simple-loop symbols", *J. Math. Anal. Appl.*
   **422** (2015) 1308-1334. DOI: 10.1016/j.jmaa.2014.09.057
   — higher-order uniform individual asymptotics under the SL condition. *Bibliographic
   data verified from two independent reference lists; not read. Our symbol violates SL.*

8. **J. M. Bogoya, S. M. Grudsky, E. A. Maximenko**, "Eigenvalues of Hermitian Toeplitz
   matrices generated by simple-loop symbols with relaxed smoothness",
   *Oper. Theory Adv. Appl.* **259** (2017) 179-212. DOI: 10.1007/978-3-319-49182-0_11
   — *bibliographic data verified from a reference list; not read.*

9. **S.-E. Ekstrom, C. Garoni, S. Serra-Capizzano**, "Are the eigenvalues of banded
   symmetric Toeplitz matrices known in almost closed form?", *Experimental Mathematics*
   **27**:4 (2018) 478-487. DOI: 10.1080/10586458.2017.1320241
   — the computational form of the `j pi/(n+1)` expansion. *Bibliographic data via
   search; full text not obtained (diva-portal fetch timed out).*

10. **P. Rambour, A. Seghier**, "Inversion des matrices de Toeplitz dont le symbole
    admet un zero d'ordre rationnel positif, valeur propre minimale",
    *Ann. Fac. Sci. Toulouse Math.* (6) **21** no. 1 (2012) 173-211.
    DOI: 10.5802/afst.1332. Preprint arXiv:1005.4073.
    — non-integer `alpha`. *Landing page read; proof not read.*

11. **A. A. Novosel'tsev, I. B. Simonenko**, "Dependence of the asymptotics of extreme
    eigenvalues of truncated Toeplitz matrices on the rate of attaining the extremum by
    the symbol", *Algebra i Analiz* **16** (2004) 146-152.
    — **title and bibliographic record only; text not obtained, theorem not seen.**

12. **D. Bini, M. Capovani**, "Spectral and computational properties of band symmetric
    Toeplitz matrices", *Linear Algebra Appl.* **52/53** (1983) 99-126.
    DOI: 10.1016/0024-3795(83)80009-3
    — origin of the tau algebra (Toeplitz-minus-Hankel matrices diagonalised by DST-I,
    eigenvalues = symbol sampled at `j pi/(n+1)`). *Bibliographic data verified; paper
    not read; the DST-I eigenvalue formula comes from secondary descriptions of the tau
    class.* **Listed for completeness only — I checked and our sine-basis section is
    NOT a tau matrix** (Section 0), so this is background, not a citation Paper 60 needs.

16. **S.-E. Ekstrom, S. Serra-Capizzano**, "Theoretical results for eigenvalues,
    singular values, and eigenvectors of (flipped) Toeplitz matrices and related
    computational proposals", arXiv:2203.06992.
    — the Toeplitz-plus-Hankel / flipped family our matrix actually belongs to.
    *Title and arXiv id verified via search; **not read**, and I did not verify it
    addresses the extreme eigenvalue at a symbol zero.*

13. **A. Bottcher, S. M. Grudsky**, *Spectral Properties of Banded Toeplitz Matrices*,
    SIAM, Philadelphia, 2005. — *not consulted directly in this scan.*

14. **A. Bottcher, B. Silbermann**, *Introduction to Large Truncated Toeplitz Matrices*,
    Springer, New York, 1999. — *not consulted directly in this scan.*

15. **S. Serra Capizzano**, "How bad can positive definite Toeplitz matrices be?" (2000)
    — the `exp(-cn)` universal lower bound. *Search-summary level only; do not quote.*

## Files written during this scan (scratchpad, not the repo)

`verify2.py` (KMS constant + symbol curvature), `verify3.py` (direct `1 - sigma_max`
in momentum space), `verify4.py` (Fourier decay of `b`, Bottcher-Widom hypothesis test)
in the session scratchpad. Downloaded full texts: arXiv:math/0412269 (Bottcher-Widom),
arXiv:1903.10551 (Batalshchikov et al., SL survey with the reference tree),
arXiv:2104.12394 (Rambour).
