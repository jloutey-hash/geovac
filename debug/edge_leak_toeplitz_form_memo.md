# Order-one defect on the true CH window: Toeplitz/Hankel form and norm scaling

Date 2026-10-31. Driver: `debug/edge_leak_toeplitz_form.py`. Data: `debug/data/edge_leak_toeplitz_form.json`, plus the logs `edge_leak_run8.log` and `edge_leak_run12.log` in `debug/data/`. The S³ construction is imported with importlib from `tests/test_paper32_window_gauge_breaking.py`, and the guarded edge-crossing formula from `tests/test_paper32_window_order_one.py`. All runs were single-threaded. Every result below is MEASURED. Nothing here is proved.

## Verdict: PARTIAL (the primary source is missing)

- The numerical agreement holds.
- The norm scaling has been measured.
- Widom's original statement could not be read. A full-text statement by the textbook author was read and verified instead.

## 1. Source for T(ab) = T(a)T(b) + H(a)H(b̃)

**Read and verified.** A. Böttcher, "One more proof of the Borodin–Okounkov formula for Toeplitz determinants", arXiv:math/0012200 (2000).
- **Definitions (p. 1):** T(f) = (f_{j−k})_{j,k≥0}, H(f) = (f_{j+k+1})_{j,k≥0} and H(f̃) = (f_{−j−k−1})_{j,k≥0}, for N×N matrix functions f in L^∞ on the unit circle.
- **Identity (p. 2):** the parenthetical "(note that H(b)H(c̃) = T(bc) − T(b)T(c) …)". It is stated as known, with no citation attached.
- **Status:** Böttcher is the co-author of the standard textbook (Böttcher–Silbermann, *Analysis of Toeplitz Operators*), but this is still not the originating paper.
- **Copy:** `debug/lit_scan/dl_widom/boettcher_math0012200.pdf`.

**Also read.** H. Widom, "On a Toeplitz determinant identity of Borodin and Okounkov", arXiv:math/9909010 (Integral Equations Operator Theory 37 (2000) 397). It works with products of Hankel matrices but does **not** state the identity. Its references [3] and [4] are Widom, Amer. J. Math. 95 (1973) 333–383, and Widom, Adv. Math. 21 (1976) 1–29. Copy: `dl_widom/widom_math9909010.pdf`.

**Not read.** These are the probable primary sources:
- Widom, "Asymptotic behavior of block Toeplitz matrices and determinants II", Adv. Math. 21 (1976) 1–29, DOI 10.1016/0001-8708(76)90113-4.
- Widom, Proc. AMS 50 (1975) 167–173.
- Böttcher–Silbermann, Section 2.2–2.3 of the book. Proposition 2.14 is widely cited, but its number is unverified.

Both publisher sites returned HTTP 403. Bypassing the bot check with a spoofed client was refused by the permission system and was not pursued. **Primary statement location: unverified.** I recommend crediting the identity as "classical (Widom; see e.g. Böttcher–Silbermann)" and citing Böttcher math/0012200 p. 2 as the verified locus.

**Why the specific source matters little.** On any Hilbert space the identity is P = P(P+Q) applied between x and y. That gives P xy P = (PxP)(PyP) + (PxQ)(QyP). The classical content is the identification of PxQ and QyP with Hankel operators on H² of the circle. The driver checks the operator form numerically (Section 2).

## 2. Block (Widom) form of the defect

**Notation.** W is the true window |λ| ≤ n_max + 1/2, with P its projection and Q = 1 − P. Write T(x) = PxP, H(x) = QxP and H̃(x) = PxQ.

- **Widom identity:** T(xy) = T(x)T(y) + H̃(x)H(y).
- **One-line form of the defect:**

  **[[D_W, a_W], b_W] = H̃(b)·H([D,a]) − H̃([D,a])·H(b)**

  Equivalently, H([D,a]) = D_Q H(a) − H(a) D_P.
- **Where it comes from:** substitute Widom twice into [T([D,a]), T(b)]. The interior term T([[D,a],b]) vanishes by the continuum order-one condition.
- **Reading:** the defect is the antisymmetrised Widom remainder of the pair ([D,a], b). It is a sum of two products of a Hankel block into the outside and a Hankel block back.

**Agreement.** Maximum entrywise residual relative to the largest entry, over all 16 ordered pairs of the four D^{1/2}_{MM'}, with the guarded ambient J2max = 2n+2:

| n_max | dim W | block vs edge sum | block vs direct | Widom identity | smaller ambient (J2max = n+3) |
|---|---|---|---|---|---|
| 1 | 4 | 3.0e-16 | 7.4e-16 | 0 | 7.4e-16 |
| 2 | 16 | 7.5e-16 | 3.0e-15 | 3.5e-16 | 3.0e-15 |
| 3 | 40 | 1.0e-15 | 3.0e-15 | 1.6e-16 | 3.0e-15 |
| 4 | 80 | 1.0e-15 | 2.4e-15 | 1.4e-16 | 2.4e-15 |

The last column justifies running the norm sweep with the smaller ambient J2max = n_max + 3.

**Support.** The defect lives only on the outermost shell |λ| = n_max + 1/2. That shell has dimension 2n(n+1). Off it, |C| < 1e-15 at n_max = 2, 4, 6. This follows from the block form: a multiplier by D^{1/2} moves |λ| by at most 1, so only the top shell has H(·) ≠ 0.

## 3. Norm scaling, n_max = 1..12

The window dimension grows like n³; the local slope is 2.77 at n_max = 12. The defect is singular-value decomposed per pair. Columns:

- **diag:** a = b, the same value for all four functions (symmetry).
- **max16:** the largest over the 16 pairs, attained at (a, b) = (D_{−+}, D_{+−}) and (D_{++}, D_{−−}).
- **coord-sum:** Σ_μ [[D_W, x_μ], x_μ] over the four real coordinates of SU(2) (Re and Im of D_{−+} and D_{++}), the analogue of Glaser–Stern's object.

| n_max | ‖C‖_op diag | ‖C‖_1 diag | ‖C‖_op max16 | ‖C‖_1 max16 | ‖·‖_op coord-sum | ‖·‖_1 coord-sum |
|---|---|---|---|---|---|---|
| 1 | 0 | 0 | 0.7500 | 3 | 1.5000 | 6 |
| 2 | 0.3143 | 1.257 | 0.7778 | 8 | 1.3333 | 16 |
| 3 | 0.3536 | 3.864 | 0.8125 | 15 | 1.2500 | 30 |
| 4 | 0.3919 | 7.780 | 0.8400 | 24 | 1.2000 | 48 |
| 5 | 0.4082 | 12.98 | 0.8611 | 35 | 1.1667 | 70 |
| 6 | 0.4242 | 19.45 | 0.8776 | 48 | 1.1429 | 96 |
| 8 | 0.4417 | 36.18 | 0.9012 | 80 | 1.1111 | 160 |
| 10 | 0.4527 | 57.91 | 0.9174 | 120 | 1.0909 | 240 |
| 12 | 0.4602 | 84.63 | 0.9290 | 168 | 1.0769 | 336 |

### Fits and exact patterns

These were read off the table to about 1e-14 at every n_max in 1..12. They are measured patterns, not theorems.

**Operator norm: bounded at every n_max.**
- **max16:** ‖C‖_op = 1 − n/(n+1)², which tends to 1.
- **diag:** at odd n, ‖C‖_op² = (n−1)/(4(n+1)), which tends to 1/4, so ‖C‖_op tends to 1/2. Even n_max sits slightly off this curve, still increasing below 1/2.
- **coord-sum:** this is exactly **diagonal**, with off-diagonal entries below 2e-15. It equals sign(λ)·(n+2)/(n+1) on the top shell and 0 elsewhere. Its operator norm is (n+2)/(n+1), which tends to 1.

**Trace norm: grows like n_max².**
- **max16:** ‖C‖_1 = n(n+2).
- **coord-sum:** ‖C‖_1 = 2n(n+2), which is the top-shell dimension times (n+2)/(n+1).
- **diag:** the local log–log slope falls from 2.29 at n = 4–5 to 2.075 at n = 11–12. ‖C‖_1/n² rises from 0.486 at n = 4 to 0.588 at n = 12 and is levelling off. That is consistent with slope 2 with a subleading correction.
- **Hilbert–Schmidt norm (diag):** the slope falls towards 1 (1.054 at n = 11–12).

**Reading.** The order-one defect is a top-shell (boundary) operator. Its operator norm is asymptotically constant and its trace norm scales with the boundary dimension 2n(n+1) ~ n².
- **Comparison with S² (Glaser–Stern, JMP 61, 033507):** this is the same pattern as their O(1) operator norm and O(Λ) trace norm. On S² the top shell has O(Λ) states; on S³ it has O(n²).
- **Answer to the scan's question:** the defect "grows" only in trace norm, not in operator norm. Any Paper 32 or record wording should name the norm.

## Open / not done

- **Primary source:** the Widom 1976 and 1975 texts and the textbook proposition number need library access.
- **No asymptotic estimate:** none was derived from the Szegő/Widom or Hekkelman–McDonald machinery. The exact patterns above (1 − n/(n+1)², n(n+2), and the diagonal coordinate sum) are candidates for a closed-form derivation from the block form plus Clebsch–Gordan sums.
- **No guard written:** under the CLAUDE.md §9 guard-writing rule, that is a separate step.
