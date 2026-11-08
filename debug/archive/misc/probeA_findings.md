# Probe A — does knowing the SW metric's spectrum EXACTLY reduce the QSVT degree for S^(−1/2)?

Pre-registered diagnostic probe. Target: Paper 60 `sec:resource`, which prices the metric
penalty as a QSVT degree `d_inv ~ kappa*ln(kappa/eps)` for a polynomial approximating
`x^(-1/2)` on `[1/kappa, 1]`, `kappa = cond(S)`. Question: the two-center Shibuya–Wulfman
metric's spectrum is `{1 +/- sigma_k}` and is *known* — classically by an `n x n` SVD, and
asymptotically in closed form by the sigma law `1 - sigma_k ~ c_sym*(k*pi)^2/n^2`
(Paper 60 `eq:sigma_law`, `geovac/sturmian_sigma_law.py`). Does knowing it buy a real degree
reduction versus the generic (continuum-accuracy) requirement?

Drivers: `debug/probeA_specaware_lp.py` (LP core) · `debug/probeA_specaware_scan.py` (main scan)
· `debug/probeA_validate_lp.py` (V1–V6) · `debug/probeA_robustness.py` (leg a: law-vs-true) ·
`debug/probeA_legs_bc.py` (legs b, c) · `debug/probeA_band.py` + `debug/probeA_band2.py`
(leg a2; `band2` is the ladder version actually used) · `debug/probeA_water.py` ·
`debug/probeA_analyze.py`. Data: `debug/data/probeA_*.json`; logs `debug/probeA_*.log`; full
tables `debug/probeA_tables.txt`.
**Diagnostic only — no paper, test, or CLAUDE.md edits.**

---

## VERDICT — CONSTANT-FACTOR (about 8–10x), and the exploit needs the exact eigenvalues

1. **Spectrum-awareness is a constant factor, not a scaling change — proved, not just measured.**
   Restricting the accuracy requirement from the continuum `[lam_min, lam_max]` to the `2n`
   actual eigenvalues reduces the minimal degree by a factor rising from ~3x at small `kappa` to
   **8.5–10x**, but *both* arms stay `Theta(kappa)`. In the one window that reaches the
   asymptotic regime (`s = 1.4`, `kappa` = 51–864) the fitted laws are
   `d_aware = 0.192*kappa^1.040` (R^2 = 0.9992) against `d_generic = 2.33*kappa^0.996`
   (R^2 = 0.9999) — **the same exponent**. A Bernstein-inequality argument (Table 5) makes this
   a theorem rather than a fit: any admissible polynomial obeys `d >= h*kappa/6` on the aware
   problem and `d >= h*kappa/2` on the generic one — both linear in `kappa`, fixed ratio 3.
   No knowledge of the spectrum can change the exponent.

2. **A *different* exploit does change the exponent — the positive-definite shift — and the
   sigma law is exactly what enables it.** Paper 60's convention leaves `lam_min` in the
   *interior* of the QSP domain `[-1,1]`, where the resolution scale is `1/d` and Bernstein
   costs `Theta(kappa)`. Because `S` is PSD with *known* `lam_min`, `lam_max` (both delivered
   in closed form by `eq:sigma_law`), one may instead block-encode the affinely shifted
   `(2S - (lam_max+lam_min)I)/(lam_max-lam_min)`, putting `lam_min` at the *endpoint*
   `-1 + 2/kappa`, where the scale is `1/d^2`. Measured: `d = 2.78*kappa^0.505`
   (R^2 = 0.9996, `kappa` = 12–864) — a clean `sqrt(kappa)`. **But** the naive LCU that
   implements the shift has 1-norm ~3, which squeezes the spectrum back into `[-1/3, 1/3]`
   (interior) and destroys the gain outright — measured, leg (b) below. The exponent change is
   real only if the shifted encoding is re-amplified to unit normalization, which costs ~3x in
   queries. **This is a resource-model claim, partially measured; it should not be asserted
   without pricing the amplification.**

3. **The aware polynomial needs the eigenvalues essentially exactly; the closed-form law alone
   is useless for it.** Designing on law-predicted eigenvalues and evaluating on the true
   spectrum gives relative errors of `1.4e+1` to `5.8e+1` — 1000%–5800%, versus a `1e-3` target.
   The mechanism is structural: the aware polynomial is a *bounded interpolant* that oscillates
   at full amplitude between nodes, so an eigenvalue displaced by the law's 11–21% error lands
   on an oscillation. The obvious repair — design on the law with a safety band wide enough to
   contain the truth — was measured and **costs the full generic degree** (leg a2): a band wide
   enough to survive the law's error is wide enough to be the local continuum. Exact eigenvalues
   are cheap classically (an `n x n` SVD, and only eigen*values* — no basis rotation, so unlike
   Löwdin it never touches the integrals or their sparsity), but the law is not a substitute for
   them, which means the exploit's price is an `O(N^3)` classical pre-step.

**Gate reading.** The pre-registered gate was "SCALING-CHANGE if `d_aware` exponent < 1.7 while
generic ~2". The honest answer is **CONSTANT-FACTOR**. The fitted `N`-exponents *do* come out
below 1.7 at `s = 2.0` (1.565) and `s = 3.0` (1.294) — but that is a pre-asymptotic artifact,
and the artifact is identified rather than guessed: for small `kappa` the aware degree is
limited by the *node count* (`d ~ 2*dim`) rather than by `kappa`, with crossover at
`dim ~ h*kappa/12`. Only the `s = 1.4` window crosses over, and there the exponent is **1.960**
against a generic **~1.8** — no exponent change, exactly as the Bernstein floor demands in
advance.

**Recommendation.** Finding 1 does **not** warrant a Paper 60 follow-on: it is a constant factor
plus a proof that it can only ever be a constant factor. Finding 2 *might* — it is an exponent
change on the paper's own headline resource number and needs only `lam_min`/`lam_max` — but it
is gated on the amplification cost that leg (b) shows is load-bearing, so it is a scoped
follow-on question ("price the amplified PSD shift"), not a result.

---

## The object actually solved

For a Hermitian eigenvalue transformation only the polynomial's values *on the spectrum* matter
for correctness, but for QSP/QSVT realizability the polynomial must additionally satisfy
`|p(x)| <= 1` on the whole rescaled interval. The honest constrained problem is therefore

> minimal `d` such that there is a `p`, `deg p = d`, with `|p(x)| <= 1` on `[-1,1]` and
> `|p(x_i) - c*lam_i^(-1/2)| <= eps*c*lam_i^(-1/2)` on the accuracy set.

Solved as a linear program per degree — `min r` s.t. `|p(x_i) - y_i| <= r*y_i` and `|p| <= 1` on
a fine grid — whose optimum `r*(d)` is the best achievable relative error at degree `d`. Then
`d(eps) = min{d : r*(d) <= eps}` by bisection (`r*` is monotone in `d`). Chebyshev basis
(`T_j(x) = cos(j*arccos x)`), for which `|p| <= 1` on `[-1,1]` implies `|a_j| <= 2`, so the LP
is well scaled.

**Four arms**, one LP machine:

| arm | accuracy set | rescaling |
|---|---|---|
| `aware / p60` | the `2n` eigenvalues | `x = lam/lam_max` in `[1/kappa, 1]` (Paper 60's convention) |
| `generic / p60` | dense grid on `[lam_min, lam_max]` | same |
| `aware / shift` | the `2n` eigenvalues | `x = affine([lam_min, lam_max] -> [-1,1])` |
| `generic / shift` | dense grid | same |

**Subnormalization convention (fixed, and it is not a gauge).** `c = h*sqrt(lam_min)` with
`h = 1/2`, identical across all four arms, so the degrees are directly comparable and `h`
cancels in every ratio. `h` is a real parameter: at `h = 1` the target *touches* the QSP ceiling
at `lam_min`, and that tangency degrades the convergence from geometric to algebraic — measured
`r* ~ 1/d`, so `d(1e-3)` would be ~2000 even at `kappa = 54`. `h = 1/2` is the standard factor-2
head-room of QSVT inversion constructions. Downstream a smaller `c` costs `1/c`
amplitude-amplification rounds; that cost is common to all arms here.

**Parity caveat.** QSP proper requires definite parity. The `p60` arm's polynomial can be taken
even (the spectrum is PSD and `|x|^(-1/2)` is even); the `shift` arm's cannot, and would be
realized as an even+odd LCU (~2x, one ancilla) or by generalized QSP. The `shift` degrees below
should therefore carry an extra factor ~2 when compared to `p60` degrees. This does not change
any exponent.

## Solver validation (`debug/probeA_validate2.log`)

* **V1** — the LP reproduces the exact Chebyshev minimax `min_{deg<=k-1}||x^k - q||_inf = 2^(1-k)`
  to 1.5e-15 (k=5), 1.1e-14 (k=8), 2.3e-5 (k=12), 6.7e-3 (k=20, at an answer of size 1.9e-6 —
  the LP tolerance floor).
* **V2** — exact interpolation recovered: `r*(d = #nodes-1) = 0` with boundedness relaxed.
* **V3** — grid independence: `r*` moves `<0.3%` between 12 and 24 boundedness points per degree.
  Production setting is 6/degree, whose Ehlich–Zeller guarantee is `sup|p| <= 1/cos(pi/12) =
  1.035`; every reported solution's true sup was re-checked on a 10x finer grid (all `<= 1.034`).
* **V5 — independent, non-LP confirmation of the minimal degree.** Parametrize the affine family
  of exact interpolants through the `2n` targets (particular solution + SVD null space) and
  minimize `sup|p|` over the family by Nelder–Mead on a smooth-max — no linear programming
  anywhere. 3/3 agree with the LP:

  | case | `d_min-1`: min sup\|p\| | `d_min`: min sup\|p\| | LP `d_min` |
  |---|---|---|---|
  | s=2.0, n=4 | 2.1017 (infeasible) | 0.9882 (feasible) | 11 |
  | s=3.0, n=4 | 1.3977 (infeasible) | 0.7578 (feasible) | 11 |
  | s=2.0, n=6 | 1.1612 (infeasible) | 0.8223 (feasible) | 21 |

* **V6** — generic-arm accuracy-grid refinement: `r*` stable to 0.2% between `n_grid` 130 and 400.
* **Leg (c), QSP head-room** — the aware degrees are not an artifact of letting `|p|` graze the
  ceiling. Re-solving with `|p| <= 0.9` (real head-room for phase-factor synthesis, where
  `1-|p|^2` must not be tiny) changes almost nothing: `s = 2.0`, `n = 6/10/14/18` give
  `19->21`, `40->41`, `71->73`, `116->118` (ratios 1.11, 1.02, 1.03, 1.02).

---

## Table 1 — minimal degree `d`, `eps = 1e-3`

`kappa = (1+sigma_max)/(1-sigma_max)`; `N = 2n` is the metric dimension;
`d_model = kappa*ln(kappa/eps)` is Paper 60 `sec:resource`. `--` = outside the arm's compute
window (`generic/p60` runs at `d ~ 2.3*kappa` and is the cost driver).

| s | n | N | kappa | aware/p60 | generic/p60 | aware/shift | generic/shift | d_model |
|---|---|---|---|---|---|---|---|---|
| 1.4 | 4 | 8 | 50.85 | 12 | 116 | 7 | 20 | 551 |
| 1.4 | 6 | 12 | 107.32 | 24 | 243 | 11 | 30 | 1243 |
| 1.4 | 8 | 16 | 183.78 | 42 | -- | 14 | 39 | 2228 |
| 1.4 | 10 | 20 | 280.21 | 66 | -- | 18 | 48 | 3515 |
| 1.4 | 12 | 24 | 396.41 | 96 | -- | 21 | 57 | 5110 |
| 1.4 | 14 | 28 | 532.39 | 132 | -- | 25 | 66 | 7020 |
| 1.4 | 16 | 32 | 688.24 | 174 | -- | 28 | 75 | 9251 |
| 1.4 | 18 | 36 | 864.04 | 222 | -- | 32 | 84 | 11811 |
| 2.0 | 4 | 8 | 25.77 | 11 | 59 | 6 | 14 | 262 |
| 2.0 | 6 | 12 | 53.91 | 19 | 124 | 9 | 21 | 587 |
| 2.0 | 8 | 16 | 91.84 | 27 | 209 | 12 | 27 | 1050 |
| 2.0 | 10 | 20 | 139.42 | 40 | -- | 15 | 34 | 1652 |
| 2.0 | 12 | 24 | 196.78 | 54 | -- | 18 | 40 | 2399 |
| 2.0 | 14 | 28 | 263.92 | 71 | -- | 20 | 47 | 3295 |
| 2.0 | 16 | 32 | 340.75 | 91 | -- | 23 | 53 | 4341 |
| 2.0 | 18 | 36 | 427.27 | 116 | -- | 26 | 59 | 5540 |
| 3.0 | 4 | 8 | 12.10 | 9 | 29 | 5 | 10 | 114 |
| 3.0 | 6 | 12 | 24.91 | 13 | 57 | 7 | 14 | 252 |
| 3.0 | 8 | 16 | 42.04 | 20 | 96 | 9 | 18 | 448 |
| 3.0 | 10 | 20 | 63.55 | 26 | 146 | 11 | 23 | 703 |
| 3.0 | 12 | 24 | 89.34 | 33 | 203 | 13 | 27 | 1019 |
| 3.0 | 14 | 28 | 119.48 | 42 | -- | 15 | 31 | 1397 |
| 3.0 | 16 | 32 | 153.96 | 51 | -- | 17 | 36 | 1839 |
| 3.0 | 18 | 36 | 192.74 | 61 | -- | 19 | 40 | 2346 |

**`eps = 1e-5`**: the aware/p60 degree barely moves — 12/24/42/66/96/132/**175**/**223** at
`s = 1.4` versus 12/24/42/66/96/132/174/222 at `eps = 1e-3`. That is structural, not luck: once
the degree admits a *bounded exact interpolant* of the `2n` eigenvalues, `r* = 0` identically and
`eps` drops out of the aware problem. The generic arms pick up the usual `ln(1/eps)`
(`generic/shift` roughly doubles, 20->36 and 84->147; `generic/p60`'s prefactor goes
2.33 -> 5.13). Full table in `debug/probeA_tables.txt`.

## Table 2 — pooled log-log fits, `d = A*kappa^p` (`eps = 1e-3`)

| arm | #pts | kappa range | exponent `p` | `A` | R^2 | max\|resid\| (dex) |
|---|---|---|---|---|---|---|
| aware/p60 | 24 | 12–864 | 0.786 | 0.893 | 0.9564 | 0.213 |
| generic/p60 | 10 | 12–107 | **0.982** | 2.461 | 0.9998 | 0.008 |
| aware/shift | 24 | 12–864 | 0.450 | 1.560 | 0.9587 | 0.116 |
| generic/shift | 24 | 12–864 | **0.505** | 2.781 | 0.9996 | 0.011 |

The two *generic* arms are near-perfect power laws landing on the two theoretical values:
`kappa^1` (interior Bernstein, Paper 60's convention) and `kappa^(1/2)` (endpoint, shifted
convention). The two *aware* arms fit worse (R^2 ~ 0.956, residuals to 0.21 dex) **because they
are not single power laws** — they cross over from node-count-limited to `kappa`-limited inside
the window. Pooling them is the wrong fit; Table 6 does it per `s`.

## Table 3 — the pre-registered gate variable: exponent in `N = 2n` at fixed `s` (`eps = 1e-3`)

| arm | s | #pts | N range | exponent | R^2 | max\|resid\| (dex) |
|---|---|---|---|---|---|---|
| aware/p60 | 1.4 | 8 | 8–36 | **1.960** | 0.9986 | 0.027 |
| aware/p60 | 2.0 | 8 | 8–36 | 1.565 | 0.9932 | 0.043 |
| aware/p60 | 3.0 | 8 | 8–36 | 1.294 | 0.9939 | 0.037 |
| generic/p60 | 2.0 | 3 | 8–16 | 1.825 | 1.0000 | 0.001 |
| generic/p60 | 3.0 | 5 | 8–24 | 1.778 | 0.9994 | 0.011 |
| aware/shift | 1.4 | 8 | 8–36 | 0.996 | 0.9989 | 0.013 |
| aware/shift | 2.0 | 8 | 8–36 | 0.968 | 0.9991 | 0.010 |
| aware/shift | 3.0 | 8 | 8–36 | 0.891 | 0.9994 | 0.007 |
| generic/shift | 1.4 | 8 | 8–36 | 0.948 | 0.9998 | 0.005 |
| generic/shift | 2.0 | 8 | 8–36 | 0.956 | 0.9997 | 0.005 |
| generic/shift | 3.0 | 8 | 8–36 | 0.932 | 0.9985 | 0.014 |

The `generic/p60` `N`-exponent is 1.78–1.83, not 2 — and that is *not* a discrepancy with
`d ~ kappa`: over this window `cond(S) ~ N^1.85` rather than its asymptotic `N^2`, which is
exactly what Paper 60 says about its own `N^1.85`/`N^1.97` fits (`0.982 x 1.85 = 1.82`).
For the gate, the like-for-like comparison at the *same* `s` in the one window that reaches the
asymptotic regime is 1.960 (aware) against ~1.8 (generic) — **no exponent change**.

## Table 6 — per-`s` fits in `kappa`, and the ratio that is the deliverable (`eps = 1e-3`)

| arm | s | fit | R^2 | max\|resid\| (dex) | `d/kappa` at largest `kappa` |
|---|---|---|---|---|---|
| aware/p60 | 1.4 | `0.192*kappa^1.040` | 0.9992 | 0.022 | 0.257 (kappa=864) |
| aware/p60 | 2.0 | `0.671*kappa^0.838` | 0.9946 | 0.040 | 0.271 (kappa=427) |
| aware/p60 | 3.0 | `1.453*kappa^0.702` | 0.9955 | 0.031 | 0.316 (kappa=193) |
| generic/p60 | 2.0 | `2.326*kappa^0.996` | 0.9999 | 0.002 | 2.276 |
| generic/p60 | 3.0 | `2.503*kappa^0.977` | 0.9998 | 0.008 | 2.272 |
| aware/shift | 1.4 | `0.899*kappa^0.529` | 0.9986 | 0.014 | 0.037 |
| generic/shift | 1.4 | `2.817*kappa^0.503` | 0.9997 | 0.007 | 0.097 |
| generic/shift | 2.0 | `2.691*kappa^0.512` | 0.9996 | 0.007 | 0.138 |
| generic/shift | 3.0 | `2.783*kappa^0.506` | 0.9991 | 0.011 | 0.208 |

The `s = 1.4` aware/p60 exponent is **1.040** — the Bernstein-limited regime; the lower `s = 2, 3`
exponents (0.838, 0.702) are the node-count-limited regime that has not yet crossed over.
`d_aware/kappa` is climbing monotonically toward the asymptote (0.224 -> 0.257 across
`kappa` = 107 -> 864), i.e. the constant, not the exponent, is what is still moving.

**The ratio table (`s = 1.4` branch; generic from its own fit where not directly computed):**

| kappa | d_aware | d_generic | ratio gen/aware | d_model | d_model/d_generic | d_model/d_aware |
|---|---|---|---|---|---|---|
| 50.9 | 12 (0.236 k) | 117 (2.29 k) | 9.71 | 551 | 4.73 | 45.9 |
| 107.3 | 24 (0.224 k) | 243 (2.26 k) | 10.11 | 1243 | 5.13 | 51.8 |
| 183.8 | 42 (0.229 k) | 411 (2.24 k) | 9.79 | 2228 | 5.42 | 53.0 |
| 280.2 | 66 (0.236 k) | 622 (2.22 k) | 9.43 | 3515 | 5.65 | 53.3 |
| 396.4 | 96 (0.242 k) | 875 (2.21 k) | 9.11 | 5110 | 5.84 | 53.2 |
| 532.4 | 132 (0.248 k) | 1169 (2.20 k) | 8.85 | 7020 | 6.01 | 53.2 |
| 688.2 | 174 (0.253 k) | 1504 (2.18 k) | 8.64 | 9251 | 6.15 | 53.2 |
| 864.0 | 222 (0.257 k) | 1880 (2.18 k) | 8.47 | 11811 | 6.28 | 53.2 |

**Decomposition of the apparent ~53x.** `d_model / d_aware ~ 53` splits into
`(d_model/d_generic ~ 5-6) x (d_generic/d_aware ~ 8.5-10)`. The first factor is *not* an
exploit — it is the `O(1)` conservatism of `kappa*ln(kappa/eps)` versus the honest minimax at the
*same* convention, which Paper 60 already flags ("the `O(1)` constant in `d_inv` cancels in the
cross-metric ratios but not in the absolute counts"). Only the second factor is spectrum
awareness. At the paper's tabulated point (`N = 16`, `kappa = 91.8`, `d_inv = 1050`) the honest
same-convention generic minimax is `d = 209` and the spectrum-aware degree is `d = 27`.

## Table 5 — rigorous Bernstein lower bounds (why it can only ever be a constant)

For `|p| <= 1` on `[-1,1]`, Bernstein gives `|p'(x)| <= d/sqrt(1-x^2)`, hence between any two
accuracy points `|p(x_b) - p(x_a)| <= d*|arcsin(x_b) - arcsin(x_a)|`. Applied to the two
smallest eigenvalues (`lam_2 ~ 4*lam_1` by the *quadratic* sigma law, targets `h` and `h/2`)
this is a proof, not a fit:

| s | n | kappa | aware/p60 floor | generic/p60 floor | aware/shift | generic/shift | floor ratio |
|---|---|---|---|---|---|---|---|
| 1.4 | 4 | 50.85 | 3.6 | 12.7 | 0.5 | 0.7 | 3.51 |
| 1.4 | 8 | 183.78 | 14.6 | 45.9 | 1.0 | 1.3 | 3.14 |
| 1.4 | 12 | 396.41 | 32.3 | 99.1 | 1.4 | 1.9 | 3.06 |
| 1.4 | 18 | 864.04 | 71.3 | 216.0 | 2.1 | 2.8 | **3.03** |
| 2.0 | 18 | 427.27 | 35.3 | 106.8 | 1.5 | 2.0 | 3.02 |
| 3.0 | 12 | 89.34 | 7.4 | 22.3 | 0.7 | 0.9 | 3.02 |

`aware/p60 floor -> h*kappa/6`, `generic/p60 floor -> h*kappa/2`: **both `Theta(kappa)`**, ratio
exactly 3. The shift-convention floors are `O(sqrt(kappa))` and two orders of magnitude weaker —
consistent with that arm's measured `kappa^0.5`. Note the *quadratic* clustering of the sigma law
is what makes the aware floor only 3x below the generic one: the bottom eigenvalues are spaced
`lam_1, 4*lam_1, 9*lam_1, ...`, i.e. by `O(lam_1)` — the same scale on which `c*lam^(-1/2)`
varies. Clustering that is *quadratic* is exactly clustering that does **not** help.

---

## Robustness leg (the pre-registered leg 4): does the closed-form law suffice?

**(a) Law-designed polynomial, evaluated on the true spectrum — catastrophic failure.**
Same degree, accuracy imposed at the law-predicted eigenvalues instead of the true ones
(`debug/probeA_robust.log`):

| s | n | d | lam_min true / law | rel. err (true-designed) | rel. err (law-designed, on true) |
|---|---|---|---|---|---|
| 1.4 | 6 | 24 | 1.846e-2 / 2.239e-2 | 1.1e-13 | **1.4e+1** |
| 1.4 | 10 | 66 | 7.112e-3 / 8.060e-3 | 9.5e-13 | **3.3e+1** |
| 1.4 | 14 | 132 | 3.750e-3 / 4.112e-3 | 1.9e-12 | **4.7e+1** |
| 1.4 | 18 | 222 | 2.312e-3 / 2.488e-3 | 3.1e-4 | **5.8e+1** |
| 3.0 | 14 | 42 | 1.660e-2 / 1.888e-2 | 6.0e-4 | **2.1e+1** |
| 3.0 | 18 | 61 | 1.032e-2 / 1.142e-2 | 1.7e-4 | **2.7e+1** |

Target is `1e-3`; the law-designed polynomial is off by 1000%–5800% on the true spectrum. The
law's own worst per-eigenvalue error is 11% (`s=3, n=18`) to 21% (`s=1.4, n=6`), always at
`lam_min` — small in absolute terms, fatal here, because the aware polynomial is a bounded
interpolant that swings to `+/-1` between its nodes.

**(a2) Band design (law nodes + safety band) — the deployable version, and it buys nothing.**
Demand accuracy on `[lam_k^law/(1+w), lam_k^law*(1+w)]` for every `k`, never touching the true
spectrum, then measure the resulting polynomial's relative error *at the true eigenvalues*
(`debug/probeA_band2.log`; `eps` target `1e-3`):

| case (kappa) | d_aware | d_generic | w | err on TRUE at `d_aware` | at `2 d_aware` | at `4 d_aware` | at `d_generic` |
|---|---|---|---|---|---|---|---|
| s=1.4, n=6 (107) | 24 | 243 | 0.05 | 3.1e-1 | 1.5e-1 | 2.3e+0 | 1.5e+1 |
| s=1.4, n=6 (107) | 24 | 243 | 0.30 | 1.3e-1 | 5.7e-2 | 1.9e-2 | **6.8e-4** |
| s=2.0, n=8 (92) | 27 | 209 | 0.05 | 9.0e-2 | 6.2e-2 | 4.1e-2 | 6.9e-2 |
| s=2.0, n=8 (92) | 27 | 209 | 0.30 | 1.0e-1 | 4.3e-2 | 8.5e-3 | **7.0e-4** |
| s=3.0, n=12 (89) | 33 | 203 | 0.05 | 4.7e-2 | 1.5e-1 | 9.7e-2 | 3.3e-1 |
| s=3.0, n=12 (89) | 33 | 203 | 0.30 | 7.0e-2 | 3.3e-2 | 4.9e-3 | **8.9e-4** |

Two clean readings. (i) `w = 0.05` **never** reaches tolerance at any degree, because the law's
worst per-eigenvalue error is 16–21% (always at `lam_min`) and a 5% band simply does not contain
the true eigenvalue — the design is accurate on a band the truth is not in. (ii) `w = 0.30`
covers the truth and does reach `~7e-4` — but only at `d = d_generic` (243/209/203), i.e. **the
band design costs the full generic degree; the entire 8–10x spectrum-awareness gain is gone.**
That is the expected mechanism: a band wide enough to survive the law's error is wide enough to
be the local continuum.

*(The `r*_band` column in the raw log — the LP's own error on the band grid — is not used for
this verdict and does not behave monotonically in `d`, because the band sampling density is tied
to `d`; the `err_on_TRUE` column above is evaluated exactly at the eigenvalues and is free of
that subtlety.)*

## Leg (b) — the shift's `sqrt(kappa)` advantage vs the block-encoding subnormalization

Shifting to `[-1,1]` needs a block-encoding of `B = (2S - (lam_max+lam_min)I)/(lam_max-lam_min)`.
The natural LCU of `S` and `I` has 1-norm `(2*lam_max + lam_max + lam_min)/(lam_max - lam_min) ~ 3`,
so without re-amplification the spectrum lands in `[-1/3, 1/3]` — interior again. Measured
(`s = 2.0`, `eps = 1e-3`, `debug/probeA_legs_bc.log`):

| kappa | aware, squeeze 1 | aware, squeeze 3 | generic, squeeze 1 | generic, squeeze 3 | generic/p60 |
|---|---|---|---|---|---|
| 25.77 | 6 | 16 | 14 | **77** (2.99 k) | 59 (2.29 k) |
| 53.91 | 9 | 29 | 21 | **168** (3.12 k) | 124 (2.30 k) |
| 91.84 | 12 | 42 | 27 | **289** (3.15 k) | 209 (2.28 k) |
| 139.42 | 15 | 60 | 34 | — | 318 (fit) |

Squeezing restores `Theta(kappa)` exactly — `d = 3.1*kappa`, i.e. the squeezed shift is **1.37x
worse than Paper 60's own convention**, not better. The entire `sqrt(kappa)` gain lives in the
spectrum reaching the endpoints `+/-1`. So finding 2 stands or falls on uniform singular-value
amplification of `B/3` back to `B`, at ~3x queries (times ~2 for the even+odd parity LCU). As a
resource model: at `kappa = 864` that is `6 x 84 = 504` effective queries against
`d_generic/p60 ~ 1880` — still ~3.7x, with the margin growing like `sqrt(kappa)`.
**The amplification cost is not measured here; it is the gating question for any follow-on.**

## Optional leg — water `A_1` block (a spectrum that is *not* of `1 +/- sigma` form)

The `C_2v`-adapted `A_1` block of water's three-center SW metric (the block holding the ground
state and the symmetry-irremovable O-H coupling, Paper 60 `sec:molecular`), `eps = 1e-3`:

| nmax | dim | kappa | aware | generic | ratio | d_model | d_model/d_aware |
|---|---|---|---|---|---|---|---|
| 3 | 6 | 46.62 | 10 (0.215 k) | 107 (2.29 k) | 10.7 | 501 | 50 |
| 4 | 8 | 82.85 | 18 (0.217 k) | 189 (2.28 k) | 10.5 | 938 | 52 |
| 6 | 12 | 182.77 | 36 (0.197 k) | 410 (2.24 k) | 11.4 | 2214 | 62 |
| 8 | 16 | 318.92 | 66 (0.207 k) | (not run: `d ~ 725`) | — | 4042 | 61 |

The ratio ~10.5–11.4 matches the two-center asymptote and the `d/kappa` constants match too
(0.20–0.22 aware, 2.24–2.29 generic), so the verdict is **not** an artifact of the `1 +/- sigma`
spectral shape — it is a property of finite dimension against a `Theta(kappa)` Bernstein floor,
and it transfers to the polyatomic block that Paper 60 identifies as the hard case.

---

## Honest caveats

* This is a **resource-model probe**. The deliverable is the ratio and the exponent, never the
  absolute degree: the `O(1)` constant in any `d_inv` model (Paper 60's caveat) applies here too,
  and the LP's answer is a *lower* bound on what a practical phase-factor synthesis achieves.
* Boundedness is enforced on a finite Chebyshev grid; the Ehlich–Zeller guarantee at the
  production density is `sup|p| <= 1.035`, and every solution was re-checked on a 10x finer grid.
  Tightening to `|p| <= 0.9` changes the degrees by 2–11% (leg c).
* Parity: `p60`-arm polynomials can be taken even; `shift`-arm ones cannot and carry a ~2x
  even+odd LCU factor not included in the tabulated degrees.
* The `generic/p60` arm is measured only to `kappa ~ 107` (it costs `d ~ 2.3*kappa` per LP);
  larger-`kappa` generic values in the ratio table come from its own `R^2 = 0.9999` fit, an
  extrapolation of at most 8x in `kappa` on a law whose residuals are `<= 0.008` dex.
* Machine was heavily loaded by other jobs throughout; this affects wall time only, not results.
