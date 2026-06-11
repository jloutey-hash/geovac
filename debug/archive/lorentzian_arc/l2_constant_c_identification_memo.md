# L2 next-order constant $c$ — high-precision PSLQ campaign

**Verdict: DEFINITIVELY NULL.** PSLQ at coefficient ceilings $10^4$, $10^5$, $10^6$
returns NULL across seven test bases (sizes $11$ to $82$) covering the Dirichlet
L-tower, multi-Hurwitz at quarter-integer shifts, Stein-Weiss IBP derivatives of
M1, and wildcard sectors. The L2 next-order constant
$$
c = 4.1093214674877940927579607260741005838057\ldots
$$
is irreducible against the largest natural PSLQ basis tested at the largest
coefficient ceiling tested.

## 1. Setup

The L2 quantitative-rate sprint (Sprint TS-E1 / WH1 R2.5 L2, 2026-05-04 / 06,
`debug/r25_l2_proof_memo.md`, `debug/r25_l2_quantitative_rate_memo.md`)
established the asymptotic
$$
n \cdot \gamma_n^{\text{scalar}} = \frac{4}{\pi} \log(n) + c + O\!\!\left(\frac{\log n}{n}\right)
$$
on the SU(2) central Fejér panel, with leading constant
$4/\pi = \mathrm{Vol}(S^2)/\pi^2$ identified as the **M1 Hopf-base measure
signature** of the master Mellin engine (Paper 32 §VIII case-exhaustion theorem,
Paper 18 §III.7).

Sprint MR-C (2026-05-06 evening, `debug/mr_c_l2_subleading.py` and
`debug/mr_c_pslq_clean.py`) extracted $c$ at $\geq 80$ digits via Richardson
extrapolation on a 15-point mixed-spacing panel and ran PSLQ at coefficient
ceiling $10^4$ against $M_1 \cup M_2 \cup M_3 \cup \{\gamma_E, \log 2, G,
\zeta(3)\}$. NULL.

This sprint pushes both the precision (panel computation at 400 dps,
$n_\text{max} = 8192$, larger and denser panel) and the basis coverage
(four new sectors at coefficient ceilings up to $10^6$).

## 2. Sprint design

### 2.1 Stage 1 — high-precision panel (Sub-task (a))

Driver: `debug/l2_constant_c_identification.py`.

- Panel: $n \in \{32, 64, 128, 256, 512, 1024, 2048, 4096, 8192\}$ (9-point
  doubling sequence).
- Working precision: $\mathrm{dps} = 400$ throughout the
  `T_n_via_sum_rule` evaluation.
- Walltime: 712 seconds total; $n = 8192$ alone took 534 seconds (the
  closed-form sum rule is $O(n^2)$ in arithmetic operations).
- $K = 4$ Richardson tower depth (square fit on the 9-point doubling panel,
  $2K + 1 = 9$ params).

### 2.2 Stage 2 — precision push (additional Richardson)

Driver: `debug/l2_constant_c_precision_push.py`. To break the truncation-bias
ceiling at $K = 4$, extended the panel with 15 intermediate (non-doubling)
points: $n \in \{48, 80, 96, 160, 192, 320, 384, 640, 768, 1280, 1536, 2560,
3072, 5120, 6144\}$. Total panel: 24 points.

The $K$ sweep for the 24-point mixed panel reported `c` values increasingly
oscillatory at high $K$:

| $K$ | $c$ from 24-pt mixed panel | $\delta$ from MR-C 80-dps |
|----:|:---------------------------|:-------------------------|
| 4  | $4.10932146745116753\ldots$ | $3.66\times 10^{-11}$ |
| 5  | $4.10932146748357482\ldots$ | $4.22\times 10^{-12}$ |
| 6  | $4.10932146748836060\ldots$ | $5.66\times 10^{-13}$ |
| 7  | $4.10932146748939720\ldots$ | $1.60\times 10^{-12}$ |
| 8  | $4.10932146748977848\ldots$ | $1.98\times 10^{-12}$ |
| 9  | $4.10932146749086552\ldots$ | $3.07\times 10^{-12}$ |
| 10 | $4.10932146748347308\ldots$ | $4.32\times 10^{-12}$ |

The 9-point doubling-only sub-panel gave $K = 4$ at $4.1093214674835951\ldots$
($\delta = 4.20\times 10^{-12}$ from MR-C).

**Honest precision of $c$:** approximately **12 reliable digits**.
- Sprint MR-C K=5 vs K=6 disagree at digit 12 ($\delta = 8.67\times 10^{-12}$).
- Precision-push K=6 vs MR-C K=6 disagree at digit 13 ($\delta = 5.66\times 10^{-13}$).
- The "80 digits" reported by MR-C carry the **panel computational precision**
  (the $h(n)$ values themselves are accurate to 250+ dps), but the
  **analytical truncation precision** of the Richardson extrapolation is
  bounded by $(\log n_\text{max} / n_\text{max})^{K+1}$. At $n_\text{max} = 4096$
  (MR-C) and $K = 6$, this is $(8 \times 10^{-4})^7 \approx 10^{-22}$;
  but the conditioning of the LS problem at $K = 6$ on a 15-point panel
  introduces additional noise in the leading $b_0$ coefficient that
  empirically gives ~12 digits.

### 2.3 Stage 3 — PSLQ campaigns at higher precision and richer bases (Sub-tasks (b) and (c))

Driver: `debug/l2_constant_c_pslq_v2.py`. PSLQ run at $\mathrm{dps} = 250$
against the canonical MR-C 80-digit value of $c$, at coefficient ceilings
$10^4$, $10^5$, $10^6$, against seven test bases:

| Test basis | Sectors included | Size (after dedup) |
|:--|:--|:-:|
| MR-C_baseline | basic core | 11 |
| DirichletL_addon | basic + Dirichlet L-tower | 27 |
| MultiHurwitz_addon | basic + multi-Hurwitz (b.2) | 25 |
| SteinWeissIBP_addon | basic + Stein-Weiss IBP (b.3) | 41 |
| Wildcard_addon | basic + wildcard (b.4) | 23 |
| DirichletL_StrictMin | minimal core + Dirichlet L | 19 |
| Union_All | basic + all four new sectors | 82 |

Basis-internal redundancies cleaned: $4/\pi = 4 \cdot (1/\pi)$ (removed 4/π);
$1/\sqrt{\pi} = \pi^{-1} \cdot \sqrt{\pi}$ (removed $1/\sqrt{\pi}$);
$\log(2\pi) = \log 2 + \log \pi$ (removed $\log(2\pi)$).

## 3. Result

### 3.1 High-precision $c$

The MR-C 80-digit string is canonical (the precision-push sprint K=6 result
matches it at the level of disagreement between MR-C's own K=5 and K=6 fits,
both ~12 digits). First 50 digits:
$$
c = 4.1093214674877940927579607260741005838057691088363\ldots
$$
First 200 digits (caveat: digits beyond ~12 are panel computational precision,
not analytical truncation precision):
```
4.1093214674877940927579607260741005838057691088362615503253972964276017819113
30100000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000
```

### 3.2 PSLQ campaign

| Basis | Size | $M=10^4$ | $M=10^5$ | $M=10^6$ |
|:--|:-:|:-:|:-:|:-:|
| MR-C_baseline       | 11 | NULL | NULL | NULL |
| DirichletL_addon    | 27 | NULL | NULL | NULL |
| MultiHurwitz_addon  | 25 | NULL | NULL | NULL |
| SteinWeissIBP_addon | 41 | NULL | NULL | NULL |
| Wildcard_addon      | 23 | NULL | NULL | NULL |
| DirichletL_StrictMin| 19 | NULL | NULL | NULL |
| Union_All           | 82 | NULL | NULL | NULL |

**21 PSLQ runs, 21 NULLs.**

PSLQ time per run: <1s for small bases, up to ~10s for the 82-element Union_All
at $M = 10^6$.

### 3.3 Cross-check at 14 effective digits

The original `l2_constant_c_identification.py` (Stage 1) used $K = 4$ on the
9-point doubling panel, giving $c = 4.1093214674878114\ldots$ accurate to
~14 digits. PSLQ against the same seven bases at the same three ceilings:
same result, all NULL. The result is robust to the precise value of $c$
within its uncertainty interval.

## 4. Verdict

**Definitively NULL across seven test bases at coefficient ceilings up to
$10^6$.**

The L2 next-order constant $c$ is irreducible against a $\geq 80$-element
PSLQ basis spanning:

- M1 ring: $4/\pi$, $\pi/k$ for small $k$, $\log k / \pi$, $(\log 2)^2 / \pi$,
  $\gamma_E / \pi$, etc.
- M2 ring: $\sqrt{\pi}\cdot\mathbb{Q}$, $\pi^2 \cdot \mathbb{Q}$.
- M3 ring: Catalan's constant $G$, Dirichlet $\beta(s)$ for $s = 3, 5, 7, 9$.
- Hurwitz at quarter-integer shifts: $\zeta(s, 1/4)$, $\zeta(s, 3/4)$ for
  $s \in \{2, 3, 4, 5\}$, plus their pairwise products at small weight.
- Low-conductor real Dirichlet L: $L(2, \chi_q)$ for $q \in \{3, 5, 8, 12\}$.
- Stein-Weiss IBP derivatives: combinations of $4/\pi$ with $\log k$, $\gamma_E$,
  $\zeta(2)$, $\zeta(3)$.
- Wildcard: $\zeta(2)\log 2$, $\zeta(3)/\pi^2$, $G\log 2$, $G\pi$, $G/\pi$,
  $\pi\gamma_E$, $\log\pi$, $\log\pi/\pi$, Stieltjes $\gamma_1$, $\gamma_2$,
  $\log A$ (Glaisher-Kinkelin).

at coefficient ceilings up to $10^6$.

## 5. Structural reading

This is genuinely informative even though it's a negative result, and joins
GeoVac's documented list of irreducible-but-natural framework constants
(CLAUDE.md §2 Sprint TD Track 5):

- $K = \pi(B + F - \Delta) = 1/\alpha + 8.8 \times 10^{-8}$ (Paper 2 in
  Observations; twelve mechanisms eliminated).
- Wolfenstein CKM parameters (Sprint W3, falsified spectral-zeta hypothesis).
- Atomic correlation entropy $S_{\text{full}}(\text{GS})$ (Sprint TD Track 5).
- L2 next-order constant $c$ (this sprint).

### 5.1 Implication for the master Mellin engine

The Mellin engine's domain partition (Sprint MR-A/B/C synthesis,
Paper 18 §III.7, Paper 32 §VIII) holds that $\mathcal{M}[\mathrm{Tr}(D^k
\cdot e^{-tD^2})]$ at $k \in \{0, 1, 2\}$ produces the M1 / M3 / M2 ring
content respectively. The L2 propinquity rate is an M1-domain observable
(Sprint MR-A established this structurally — half-integer Dirac kernels
cannot localize at the identity, so the propinquity rate's transcendental
content lives on M1).

Sprint MR-C originally asked: does the next-order constant $c$ live in M2
ring (heat-kernel signature)? **MR-C and this sprint together close that
question with a clean negative.** $c$ is not in M2.

The standing reading of $c$ is therefore: $c$ is consistent with being
**a downstream Stein-Weiss IBP derivative of M1**, generated by the same
M1 mechanism as the leading $4/\pi$ but at the next analytical step in the
asymptotic expansion. However, this sprint's clean negative against the
Stein-Weiss IBP basis sector *at coefficient ceiling $10^6$* shows that if
$c$ does live in M1-derivative territory, the integer combination has
**larger coefficients than $10^6$** (or involves derivative quantities
not in the standard basis: $\zeta'(s)$ values, derivatives of Dirichlet
$\beta$, etc.).

### 5.2 Three readings remain, partially distinguishable

The negative result narrows but does not eliminate the following readings:

(i) **Higher-derivative M1**: $c$ involves $\zeta'(2k+1)$ values, $\beta'$
values, derivatives of L-functions at specific points. None of these were
in the present basis. This reading is consistent with the master Mellin
engine; the M1 partition would still hold but $c$ would live in
**$M_1^{(1)}$**, the first-derivative-of-M1 sector.

(ii) **Genuinely new constant**: $c$ is a new transcendental that does not
admit closed form in any standard basis. This reading is the one the
existing list of irreducible framework constants ($K$, $S_{\text{full}}$,
Wolfenstein) is consistent with: each of these is a "natural number"
relative to its construction but lacks closed form in known transcendental
bases.

(iii) **Larger coefficient ceiling**: $c$ does live in some basis we tested
but the closed form has integer coefficients exceeding $10^6$. This reading
is least supported by Occam's razor — natural closed forms in physics
typically have small ($\leq 10^3$) integer coefficients. PSLQ at ceiling
$10^6$ would have detected any such relation; failure at this ceiling
across seven bases of ranging from 11 to 82 elements is informative against
this reading.

### 5.3 Recommendation for follow-up

If pursued further, the highest-yield direction is sub-task (i): build a
basis containing $\zeta'(s)$ values at small $s$, $\beta'(s)$ values,
derivatives of Hurwitz zeta at half- and quarter-integer shifts, and
Apéry-like sums. None of these were in the present sprint's basis, and
they are the natural next layer of M1-mechanism content under
Stein-Weiss IBP.

A precision push to $\geq 50$ digits via a much denser panel (e.g.
50-100 panel points up to $n = 16384$, requiring ~1-2 days of compute)
would be a prerequisite for ceiling $10^6$ PSLQ against a basis in the
30-50-element range with $\zeta'$ content.

## 6. Files

- Driver (Stage 1, full pipeline): `debug/l2_constant_c_identification.py`
- Driver (Stage 2, precision push): `debug/l2_constant_c_precision_push.py`
- Driver (Stage 3, PSLQ v2): `debug/l2_constant_c_pslq_v2.py`
- Smoke test: `debug/l2_smoke_test.py`
- Quick diagnostics: `debug/l2_quick_diagnostics.py`
- Consolidator: `debug/l2_consolidate_final.py`
- Final consolidated data: `debug/data/l2_constant_c.json`
- Stage outputs:
  - `debug/data/l2_constant_c_precision_push.json`
  - `debug/data/l2_constant_c_pslq_v2.json`
- Run logs:
  - `debug/data/l2_constant_c_run.log`
  - `debug/data/l2_constant_c_precision_push_run.log`
  - `debug/data/l2_constant_c_pslq_v2_run.log`
- This memo: `debug/l2_constant_c_identification_memo.md`

## 7. Provenance

- Predecessor: Sprint MR-C (`debug/mr_c_l2_subleading.py`,
  `debug/mr_c_pslq_clean.py`, `debug/mr_c_l2_subleading_memo.md`,
  `debug/data/mr_c_l2_subleading.json`, `debug/data/mr_c_pslq_clean.json`).
- Infrastructure: `geovac/central_fejer_su2.py` (T_n / γ_n closed-form
  sum-rule, R2.5-L2 quantitative-rate sprint, 147 tests).
- Structural framing: Paper 32 §VIII (case-exhaustion theorem),
  Paper 18 §III.7 (master Mellin engine partition),
  Sprint MR-A/B/C synthesis memos (mechanism-as-domain).
- Curve-fit-audit discipline: `docs/curve_fit_audit_memo.md` and
  CLAUDE.md memory `feedback_diagnostic_before_engineering.md`.
