# BU: the one-electron bond-angle operator derived bottom-up (case N1)

Pre-registration: `docs/predictions/2026-10-28_bottom_up_angle_operator.md` (gate GO / BORDERLINE / STOP, P(GO) = 0.6).

## Outcome: GO

The derived operator D, of order 9 and degree 30 in the ansatz sum a_ij x^j T^i (T = (1 - x^2) d/dx), is a left
multiple of NA's guessed operator G, of order 8 and degree 31. In standard form the exact Ore right division
D = M o G over F_p(x) leaves remainder 0, with M a first-order operator. This holds mod 2^61 - 1 and mod
4611686018427387847. Under the gate's second GO clause, the guessed operator right-divides the derived one.

## Route (no fitting to coefficient data)

1. **Radial families and IBP** (`bu_ibp.py`, sympy, L symbolic). The three radial pieces of F_L are integrals of
   rank-2 families u_L = e^{-zX r} r^{L+2} i_L(zH r) on [0, D1], e^{-zX r} r^{L+2} k_L(zH r) on [D1, D2], and
   e^{-zX r} r^{1-L} k_L(zH r) on [D2, oo). Each family's d/dr system and its L -> L+2 shift come from the
   three-term and derivative rules of i_L and k_L. Integration-by-parts identities int d/dr(r^k u_j) = [r^k u_j]
   reduce every integral needed to masters plus boundary values. These identities hold for na_exact's formal
   integrals as well, because its monomial rule b J_n + n J_{n-1} = [r^n e^{br}] is the fundamental theorem of
   calculus.
   - With k >= 0 only, two masters per piece are needed.
   - With negative k allowed, every family has a first-order relation of the form int u_L = rho(L) int u_{L+1}
     plus boundary terms, so a single master per piece is enough (the "rank-one" reduction). On [0, D1], the r = 0
     end of the boundary terms vanishes for the k that occur, since u_j ~ r^{2L+2+j}.
   - Every relation was checked by mpmath quadrature, including on [0, 1.8] and [2.6, oo), at L = 0, 1, 2, 5.
     Worst relative error: 3e-38.
2. **First-order system in L** (`bu_system.py`). The coordinates U_L are 14 formal sequences: multiplier Bessel
   pairs tensored with masters, plus the boundary products at r = D1 and r = D2 (12 + 6 such terms in the two-master
   version). Transcendental constants sit inside the coordinates, so all coefficients lie in Q(L). The result is
   U_{L+1} = A(L) U_L and F_L / (2 pi N) = c(L) . U_L. Checked numerically against quadrature of the defining radial
   integral for L = 0..5: the system to 6e-48 and the output to 4e-49.
3. **Recurrence in L** (`bu_derive.py rec`). The cyclic vector rows r_k = c(L+k) A(L+k-1)...A(L) become dependent
   at order m. The coefficients come from kernels at random points of F_p, then rational reconstruction; they were
   stable at 300 and 500 points.
   - Two-master system (18-dim): order 11, degree 62 in L.
   - Rank-one system (14-dim): order 9, degree 40 in L (coefficient degrees 36..40).
4. **x-operator** (`bu_derive.py angle`). Unknown constants a_ij were solved for so that the coefficient-space image
   of sum a_ij x^j T^i annihilates c . U for **every** solution U of the system. Two kinds of rows impose this:
   - at random symbolic points lam, 250 points, 9 independent rows per point;
   - at the integer centres l < R + Q, propagating from L = 0 only. Coefficient-space X and T never couple index
     n < 0 to n >= 0. Without these rows the search returned spurious order-7 operators: they held for generic L but
     failed at the lowest ~30 L on the r = 0 component. The cause is that A(L) has poles at L = -1, -2.

   Order/degree trade, with nullity in parentheses: order 8 has nothing up to degree 256; 9 at degree 30 (1);
   10 at 22 (2); 11 at 18 (1); 12 at 17 (4). The (9, 30) operator is unique (degree 29 gives nullity 0) and the same
   on both primes. In standard form its coefficient degrees are [30, 31, 32, 33, 34, 35, 35, 36, 36, 37].

## Checks (data used only here)

- `bu_data.py` rebuilds the generic N1 sequence exactly as na_run.py does, for L = 0..1000 on both primes.
- `bu_check.py rec`: both recurrences (11/62 and 9/40) give zero on all rows (992 for the 9/40 one), for the
  generic sequence and for each of the six constant components, on both primes.
- `bu_check.py angle 9 30`: D gives zero on all 962 rows, generic and each component, on both primes.
- `bu_compare.py` is the only place G is computed, with na_run's call G.nullity(B, 8, 31, 30), nullity 1. The
  exact right division of D by G has remainder 0 on both primes.
- `bu_compare_ctrl.py`: the same division FAILS (non-zero remainder) for G with one coefficient changed, for a random
  operator of G's shape, and for D's own order-8 truncation, on both primes. The test can therefore reject.

## Reading

The derived family is gen's solution space plus one extra solution. The derived trade is NA's trade shifted up by
exactly one order (NA: (8, 31), (9, 22), (10, 18), (11, 17), (12, 16)). The explanation is that the system's
solution family is slightly larger than the module of the specific generic sequence. Removing that excess would
need initial conditions, which is data, so it was not attempted. The rank-one IBP relations are what brought the
family down to within one solution of gen. With two masters per piece, no operator of order 8 exists up to degree
256, and the scan was not continued.

## Errors found

None in NA's code or result: the derivation confirms the guessed operator. Pitfall recorded above: imposing the
identity only at generic L gives operators that miss the L = 0 boundary.

## Files

`debug/bu_ibp.py`, `bu_system.py`, `bu_derive.py`, `bu_data.py`, `bu_check.py`, `bu_compare.py`, `bu_compare_ctrl.py`;
`debug/data/bu_system_v{1,2}.json`, `bu_recurrence_v{1,2}_<p>.json`, `bu_angle_v2_<p>_9_30.json`,
`bu_angle_scan_v2_2305843009213693951.json`, `bu_check_seq_<p>.json`, `bu_compare.json`, `bu_data.log`.
