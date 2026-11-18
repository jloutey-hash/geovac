# Zhou and Matuzas against Paper 63: two light checks (2026-11-18)

Scripts: `debug/lit_scan_zhou_check.py` (7 min, one process), `debug/lit_scan_matuzas_check.py`
(under 1 min). Data: `debug/data/lit_scan_zhou_check.json`, `debug/data/lit_scan_matuzas_check.json`.
No paper, record or git state was touched. Both sources were read from the arXiv PDFs
(Zhou pp. 1-6; Matuzas pp. 1-8), not from abstracts alone.

## Part A. Zhou, arXiv:1911.04141 (Ramanujan J. 58, 723)

### What Zhou proves (read from the PDF)

- Object: IKM(a,b;2k+1) = int_0^oo I0(t)^a K0(t)^b t^(2k+1) dt. **Odd moments only.**
- Thm 1.2: for 1 <= a < b, the span over Q of {IKM(a,b;2k+1)} has dimension <= floor((a+b-1)/2).
  For a = 0: {1} with {IKM(0,n;2k+1)} has dimension <= floor((n+1)/2), the same number as the
  Borwein-Salvy recurrence gives; the new part there is only which moments form the basis.
- Thm 1.3 (exceptional): IKM(2,6;1) = 72 IKM(2,6;3), with companions for (4,4), (3,5), (1,7) whose
  right sides are 7 log 2 / 2, -5 pi^2/12, 7 pi^4/48, and one quadratic rule, Eq. (1.16).
- Tool: Vanhove's operator L~_n, of order n, acting on the off-shell diagram
  int I0(sqrt(u) t) I0^(a-1) K0^b t dt with a + b = n + 2 Bessel factors, evaluated at u = 1.

### Which objects coincide with ours

Paper 63: F_n(t) = int_0^oo K0(tx) K0(x)^(n-1) dx (prop:tree_bound, line 347), weight x^0.

1. **Zhou's moments are not ours.** F_n(1) = c_{n,0} = IKM(0,n;0), an **even** moment with a = 0.
   Zhou's theorem covers odd moments. The Borwein-Salvy recurrence links moments of one parity
   only, so no relation of Zhou's touches F_n(1) or its derivatives.
2. What does line up is the operator count. Zhou's L~_n has order n for n + 2 Bessel factors, that
   is the two-dimensional banana with n + 1 propagators. Paper 63 line 1083 says F_n "sits on the
   rung of the banana with n+1 propagators" with an operator of order n. Same order, same rung.
   Zhou's moments are that banana at the on-shell point (I0(t) = J0(it), u = 1); F_n is the
   banana "integrated over a line of external momentum" (line 1082). Different slices of one family.
3. The Taylor data of F_n at t = 1 is spanned by even moments. With theta = t d/dt,
   theta^j K0(tx) is a polynomial combination of x^(2k) K0 and x^(2k) (x K1), and
   int x^(2k+1) K1 K0^(n-1) dx = ((2k+1)/n) c_{n,2k}. So theta^j F_n(1), j < n, lies in the span
   of c_{n,0}, ..., c_{n,2 floor((n-1)/2)}: at most floor((n+1)/2) numbers, against an operator
   of order n.

### Numerical verification (mpmath, 50 digits, plain quadrature)

Our evaluators compute F_n and do not apply to I0-weighted odd moments, so these used direct
quadrature. Relative residual = |lhs - rhs| / (sum of |terms|).

| Relation | Residual |
|:--|:--|
| (1.4) int I0 K0^4 t(16 - 228t^2 + 45t^4) = 0 | 4e-52 |
| (1.5) int K0^5 t(16 - 228t^2 + 45t^4) = 24 | < 1e-50 |
| (1.6) int I0^2 K0^4 t(2 - 85t^2 + 72t^4) = 0 | 6e-52 |
| (1.7) int I0 K0^5 t(2 - 85t^2 + 72t^4) = 0 | 1e-51 |
| (1.8) int K0^6 t(2 - 85t^2 + 72t^4) = 15/2 | < 1e-50 |
| BBBG: 45 IKM(2,3;5) = 228 IKM(2,3;3) - 16 IKM(2,3;1) | < 1e-50 |
| **(1.14) IKM(2,6;1) = 72 IKM(2,6;3)** (exceptional) | < 1e-50 |
| (1.13) IKM(3,5;1) - 72 IKM(3,5;3) = -5 pi^2/12 | 1e-51 |
| (1.15) IKM(1,7;1) - 72 IKM(1,7;3) = 7 pi^4/48 | 6e-52 |
| (1.16) quadratic rule, right side 45 pi^6/16 | 6e-52 |

Controls: coefficient 71 in place of 72 on (2,6) leaves 7e-3; (1.15) with its right side dropped
leaves 0.39. (1.12), the (4,4) rule with the "honorary" moment, was not run (algebraic tails).

Borwein-Salvy recurrence (Zhou Eq. 1.2-1.3) on pure K0 moments, including the even class ours
lives in: n = 5, 6 at k = 0, 1, 2 hold to 2e-47 or better (even, k = 0: c_{5,0} - 731 c_{5,2} +
2435 c_{5,4} - 225 c_{5,6} = 0 and c_{6,0} - 2736 c_{6,2} + 24000 c_{6,4} - 9216 c_{6,6} = 0).

### Dimension bound against our operator orders, n = 3..8

Using Paper 63's own operator (`_t_ode(n)` from `tests/test_paper63_four_density_rung.py`),
evaluated at t = 1 with the jet written exactly in the symbols c_{n,2k}:

| n | operator order | rank of the jet {theta^j F_n(1)}, j < n | floor((n+1)/2) | operator's linear form at t = 1 |
|:--|:--|:--|:--|:--|
| 3 | 3 | 2 | 2 | identically 0 |
| 4 | 4 | 2 | 2 | identically 0 |
| 5 | 5 | 3 | 3 | identically 0 |
| 6 | 6 | 3 | 3 | identically 0 |
| 7 | 7 | 4 | 4 | identically 0 |
| 8 | 8 | 4 | 4 | identically 0 |

Reading:

- The jet rank equals floor((n+1)/2) exactly, the Borwein-Salvy number. The drop from n to
  floor((n+1)/2) comes from the elementary integration by parts above (all n masses equal at
  t = 1), not from Zhou.
- The operator at t = 1 imposes **no** rational relation among the even moments, for odd n
  (regular point) and even n (singular point) alike. This is a real check on the operators: a
  wrong coefficient would generically leave a nonzero form, which would assert a false Q-relation
  among c_{n,0}, c_{n,2}, .... Because the form is identically zero, the numerical evaluation of
  it adds nothing.
- **Zhou's bound floor((a+b-1)/2) is about a different thing** from "order <= n". It counts
  rational relations among odd on-shell moments with at least one I0. It neither reproduces nor
  sharpens the tree rung; it is consistent with it only through the shared banana order in
  item 2. Paper 63's credit line (line 1994: "rational linear relations among Bessel moments,
  from Vanhove's operators, are Zhou's") is accurate as written.

## Part B. Matuzas, arXiv:2610.00036

### The object (read from the PDF)

I(alpha) = int over R^12 of exp(-sum alpha_i r_i) / (r1 r2 r3 r4 r12 r23 r34 r41). One centre.
Four Coulomb links in a cycle; each electron carries exp(-alpha r)/r (a "Slater spoke" with a
1/r), not a 1s density. Result (Thm 3.1): (4 pi)^4/P times (7/8 zeta(3) + S/2), S a sum of
eighteen blocks of log, Li2, Li3 of rational functions of the exponents, weight at most three.
Equal exponents (Cor. 4.2): (4 pi)^4/lambda^4 times
(Li3(1/4)/2 - 7 zeta(3)/24 + pi^2 log 2/12 - 2 log^3 2/3). The two chorded topologies are excluded.

### (1) Same graph?

Same topology, different weights and different question.

- Our ring of four (record `docs/predictions/2026-10-17_ring_extensions.md`, QR): S(t) =
  <rho_a rho_b rho_c rho_d / (r_ab r_bc r_cd r_da)>, a square of four Coulomb kernels, with
  fingerprint-fibre densities (transform (c k^2 + 1)^(-3/2), i.e. rho proportional to K0(a r)) at
  scales (1,1,1,t). QR1 asks for the rung of S(t) as a function of t.
- His is the same square at one centre with weights exp(-alpha r)/r. A 1s Slater density
  exp(-alpha r) is minus the alpha-derivative of his weight, so the one-centre Slater square of
  our QR2 is the fourth mixed alpha-derivative of his I(alpha). Differentiating logs and
  polylogs of rational arguments stays in the same class.
- So his result is the one-centre limit, QR2, not the fibre rung, QR1.

### (2) His result class against what Paper 63 says

Paper 63 says nothing about the class of a ring of four:

- line 426: "how a ring of four behaves;\ the second test stopped for feasibility."
- line 2047: "Fourth, among rings only the triangle was computed"
- line 2053 lists "the ring of four" as open.
- lines 1136-1141 already cite the paper: "evaluates the cyclic integral without the diagonal in
  closed form, in logarithms, dilogarithms and trilogarithms". That matches the PDF.

The record's QR2 ("At one centre with 1s Slater densities the square is not rational, by analogy
with the triangle's dilogarithms"), UNSCORED, is **supported** by Matuzas: weight at most three,
with a trilogarithm at equal exponents. It fits the triangle pattern in sec:density_ring
(chain rational, triangle weight two at one centre, line 1190: "The ring raises the weight from
zero to two"): the square reaches weight three. QR1 is untouched: for the triangle the
one-centre value is a dilogarithm while the fibre sits on the elliptic rung, so the one-centre
class does not fix the rung. Verdict: consistent with QR2, orthogonal to QR1, contradicts nothing.

### (3) Numerical comparison

The formula is reachable. Checked at 40 digits:

- Cor. 4.2 at (1,1,1,1) against his Table 1 (3159.60899251365163791156042906108): all 32 printed digits.
- His one-dimensional form, Prop. 2.1, at (1,1,1,1), (1,2,3,5), (1,1,2,2) against Table 1: all 32 printed digits each.
- Independent route (angular reduction and the 16/8 sector count re-derived here; Gauss-Legendre, float64, equal exponents): 2e-15 relative against Cor. 4.2.

The eighteen-block Theorem 3.1 itself was not coded; Prop. 2.1 and Cor. 4.2 were.

Comparison with our own numbers was **not possible**. Our stored ring-of-four values
(S(0.7) = 0.0220238, S(1.6) = 0.0491398, `debug/data/ring_ext_QR_mom.json`) are for K0-type
densities at scales (1,1,1,t). Mapping his formula to them needs a four-fold integral over a
fourth derivative of his closed form; that is not a light computation. QR2 was never computed on
our side, so there is no one-centre number of ours to compare.

Follow-up within reach, not done: the one-centre Slater square from the same sector form, then
an integer-relation search in his letters. That would score QR2.

## Sentences in Paper 63

No change is required by either check.

One existing sentence is worth a second look, though this check does not decide it. Lines
1131-1134: "one-center three- and four-electron correlated integrals have closed
forms~\cite{fromm_hill1987,padhy2017}, consistent with the rung needing two-center densities."
Three lines later the paper reports Matuzas's other result, a one-centre ring with a diagonal tied
to a Calabi-Yau operator (AESZ-330). The two sit together only if "the rung" is read as the tree
ladder's. The Zenodo source for that result was not read here.

Optional, PI's call: at line 2053 the open item "the ring of four" could note that its one-centre
limit is now known in closed form at weight three (Matuzas), with the fibre rung still open.
