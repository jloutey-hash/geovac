# HG occupation weights: second route from the published parabolic-spherical coefficients

**Date:** 2026-10-31. **Driver:** `debug/hg_mardoyan_second_route.py` (exact `Fraction` arithmetic; no CG routine
in the book route). **Target:** record HG (`docs/predictions/2026-10-31_hopf_grouping.md`), occupation
p(n,l,m; mu) = |<j, m-mu; j, mu | l, m>|^2, j = (n-1)/2, computed by sympy CG as in
`tests/test_paper25_hopf_grouping.py::hopf_M`.

**Verdict: AGREE.** The book route reproduces every occupation exactly at n_max = 1, 2, 3, 4. Two planted errors are
caught. Two other planted errors are not caught, and the reason is a symmetry, explained below.

## Source and transcription

Mardoyan, Pogosyan, Sissakian, Ter-Antonyan, arXiv:2310.17336, Sec. 2.2 ("Solution to the Park-Tarter problem"),
text re-extracted with `pdftotext -layout` and `-raw`.

- Expansion (2.2.17): Psi_{n1 n2 m} = sum_{l=|m|}^{n-1} W^{nl}_{n1 n2 m} Psi_{nlm}, with n = n1 + n2 + |m| + 1.
  Spherical functions as in (2.2.18)-(2.2.21) (Condon-Shortley-type phase (-1)^{(m-|m|)/2}); parabolic functions
  (2.2.19), (2.2.22), with xi = r(1 + cos theta), eta = r(1 - cos theta) (2.2.23).
- CG form (2.2.25): W = (-1)^{l + n2 + (m-|m|)/2} C^{l,|m|}_{(n-1)/2, (n1-n2+|m|)/2 ; (n-1)/2, (n2-n1+|m|)/2}.
- **Route used: the 3F2 form (2.2.29)** (Tarter's matrix up to the phase (-1)^{(m+|m|)/2}, as the book notes):

  W = (-1)^{l+(m-|m|)/2} sqrt[(2l+1)(n1+|m|)!(n2+|m|)!(l+|m|)! / (n1! n2! (l-|m|)! (n-l-1)! (n+l)!)]
      x (n-|m|-1)!/|m|! x 3F2(-n2, -l+|m|, l+|m|+1; |m|+1, -n+|m|+1; 1).

  The layout extraction does not show where the square root ends. Read with (n-|m|-1)!/|m|! outside the root, the
  matrix W^2 has every row and every column summing to exactly 1 for every (n, m) with n <= 4, which is checked in
  the driver. That fixes the reading. The phase drops out of the square.

## Label map

In (2.2.25) the first spin-j factor has projection (n1-n2+|m|)/2 and the second (n2-n1+|m|)/2. For m >= 0 this is
<j, m-mu; j, mu | l, m> with **mu = (m - n1 + n2)/2**. For m < 0, (2.2.25) uses |m|. Two symmetries hold when
j1 = j2 = j: CG^2 is unchanged by flipping all projections, and unchanged by exchanging m1 and m2. Either one turns
the m < 0 case into the same formula, mu = (m - n1 + n2)/2. So the map holds for all m: in the CG ordering the
parabolic state (n1, n2, m) has A_z eigenvalue n1 - n2, and it is the J-_z = mu eigenstate. Within a shell, mu runs
over -j..j, and the m+ = m - mu label is collapsed. The group (n, mu) is therefore the set of Stark states with
n2 - n1 = 2mu - m, which is the regrouping identified by the prior-art scan.

## Result

| n_max | (node, group) entries | nonzero | max abs(book - record) |
|:-----:|:---------------------:|:-------:|:----------------------:|
| 1 | 1  | 1  | 0 |
| 2 | 7  | 7  | 0 |
| 3 | 26 | 25 | 0 |
| 4 | 70 | 67 | 0 |

The agreement is exact rational equality, including the exact zeros. One example of a zero is l = 1 at n = 3,
m = 0, mu = 0, where <1,0;1,0|1,0> = 0.

## Controls (n_max = 4)

| Planted error | Max diff | Discriminates? |
|:--|:--:|:--:|
| C-b: mu = (n2-n1)/2 (the m term dropped) | 1, at (3,2,2; mu=1) | yes |
| C-e: mu -> -mu at fixed m (L_z sign error) | 1, at (4,3,3; mu=-3/2) | yes |
| C-d: 3F2 lower parameter -n+abs(m)+1 -> -n+abs(m)+2 | 1332/5 | yes (row sums also break) |
| C-a: n1 <-> n2 swapped at fixed m | 0 | **no** |
| C-c: the swapped map for m < 0 only | 0 | **no** |
| C-f: -n2 -> -n1 in the 3F2, with the swapped map | 0 | **no** (same as C-a) |

**Why the n1 <-> n2 swap suggested in the brief cannot discriminate.** At fixed m, swapping n1 and n2 sends mu to
m - mu. Because j1 = j2, CG^2(j, a; j, m-a | l, m) = CG^2(j, m-a; j, a | l, m). The driver verifies this identity
for every case with n <= 4. So the occupation of every node in group mu equals its occupation in group m - mu, and
the swap reproduces the record's occupation matrix exactly. This is not an accident of the test. It is the same
fact as the record's F3: grouping by m+ instead of m- gives identical results, because the occupations cannot tell
J+ from J-, or the sign of A, apart. For the same reason no convention on the sign of A_z can be tested through
squared weights.

The convention errors that do change the squared weights are a lost m term in the map (C-b) and a sign error on L_z
(C-e). Both are caught at the first shell where they can act, which is n = 3 for C-b. A corrupted coefficient
formula (C-d) is caught both by the comparison and by the unitarity self-check. Agreement is therefore not
automatic: the route checks the map's dependence on m and the coefficient values. By construction it cannot check
the A-sign or J+/J- labelling, and HG does not depend on either.

## Scope

This check covers n <= 4, which includes both of HG's gated cutoffs (n_max = 3 and 4). It confirms the occupations
only. It does not check the base-graph scoring, which uses the occupations together with the adjacency.
