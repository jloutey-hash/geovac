# B2: the two-electron bond-angle structure, built bottom-up (case W1)

Pre-registration: `docs/predictions/2026-10-28_two_electron_bottom_up.md` (commit ab16914a).
Case W1: D1 = 17/10, D2 = 5/2, zX = 7/5, zH = 11/10.

## Outcome: BORDERLINE

- **B2-1 holds.** A recurrence in L was derived and passes B2-0.
- **B2-2 holds.**
- **B2-3 fails as worded.** At lam = 17/25 the exponents include beta = 3 and beta = 1 alongside -4, -5 and -6. The
  wording requires every beta to be an integer <= -4. The beta = -4 that B2-3 requires is present, and so are -5 and
  -6.

The gate reads GO only if B2-1, B2-2 and B2-3 all hold, and BORDERLINE otherwise.

## Construction (no fitting to coefficient data)

`debug/b2_engine.py` writes F = I_C + I_T. Here I_C = int_0^oo a1(t) Phi2(t) dt and I_T = int_0^oo a2(t) Phi1(t) dt,
with a_j = g_j t^{-L-1} and Phi_j(t) = int_0^t g_j s^L. For t > D_j, Phi_j(t) = Phi_j(oo) - int_t^oo (outer form), so
every inner integral starts at 0 or ends at infinity. That means no constant terms enter the inner functions.

Modules, each with a d/dr matrix and an L-shift matrix whose entries are Laurent polynomials in r:
- U[Y] = (u_L, u_{L+1}), with u_L = e^{beta r} r^{L+2} Y_L(zr).
- V[Y] = (v_L, v_{L+1}), with v_L = e^{beta r} r^{1-L} Y_L(zr).
- The nested modules N[VY, UY, side] = span{v_c psi, v_c u_d}, of rank 6, where psi = int_0^t u_{L+1} or
  int_t^oo u_{L+1}. They are built from d psi/dt = +-u_{L+1} and from U's rank-one shift psi(L+1) = rho psi + boundary
  terms. The variants used are N[i,i,0], N[k,i,0], N[k,k,oo] and N[i,k,oo].

IBP reduction:
- The identities int d/dr(r^k e) = [r^k e] are put in row-echelon form mod p at a point L = lam of F_p.
- Windows run over k in [-6, 6]; for U, k >= -1.
- Integrals reduce to masters plus boundary values at D1 and D2. The ends 0 and infinity vanish.
- On side 0, k = -6 is valid only for L >= 2. The checks below show the result holds from L = 0 anyway.
- U and V reduce to one master each (J(0,1)), as in BU. The nested modules have about 7 non-edge masters each.

Coordinates:
- Each coordinate is a monomial e^{tag} D1^{aL} D2^{bL} times a product of atoms.
- The atoms are the Bessel pairs (Y_L, Y_{L+1})(z t0) at t0 in {D1, D2}, and the integral masters (module, interval,
  master).
- The shift of a coordinate is the product of the shifts of its atoms. An integral atom's shift is its reduced
  master shift plus boundary monomials.
- Closing under the shift, starting from the 140 monomials of the output, gives **a system of dimension 180**:
  U_{L+1} = A(L) U_L and F_L = c(L) . U_L.

`debug/b2_derive.py` finds the recurrence with the cyclic-vector method:
- The rows r_k = c(L+k) A(L+k-1) ... A(L) first become dependent at **order 41**.
- Normalised kernel vectors were computed at 4200 consecutive points of F_p, which share matrices. Each p_k/p_41 was
  then reconstructed as a rational function.
- Result: **order 41, degree 926 in L.** The coefficient degrees are 906, 908, ..., 926 (plateau), ..., 908, 906,
  symmetric, and identical on both primes.
- 500 and 1600 points were too few: the reconstruction returned degree n - 1, meaning it failed. Those outputs were
  deleted.

## B2-0 (checks; data used only here)

`debug/b2_data.py` rebuilds the exact sequence with `ao_exact.coefficients` for L <= 300 on both primes. It finds 31
constants, sorts the keys by repr, and assigns random.Random(2026) weights in that order, as `ao3_run.py` does.
`debug/b2_check.py` then runs two checks.

| Check | Prime 2^61 - 1 | Prime 4611686018427387847 |
|:--|:--|:--|
| Pointwise: kernel of the cyclic rows at an integer L0 (no reconstruction) applied to the data | L0 = 2, 3, 5, 10, 20, 100, 200, 259: one kernel vector each; the generic sequence and all 31 components give 0 | L0 = 2, 3, 50, 150, 259: same |
| Reconstructed order-41 recurrence applied to the data | 0 failures on all 260 rows (L = 0..259), generic and all 31 components | 0 failures, same |

At L0 = 1 the pointwise kernel has dimension 9 (a degenerate integer point). At L0 = 0 the IBP pivot structure
changes, so that point was skipped. The reconstructed recurrence nonetheless holds from L = 0.

The radial pieces were not checked by quadrature. The exact identity on the data, for every component, is the
decisive check.

## b2_asym.py output (identical on both primes)

```
{'root': '-1875/2057', 'x': '-3873437/3856875', 'mult': 1, 'betas': [('-1/2', 1)]}
{'root': '-2057/1875', 'x': '-3873437/3856875', 'mult': 1, 'betas': [('-1/2', 1)]}
{'root': '-121/51', 'x': '-8621/6171', 'mult': 1, 'betas': [('-1/2', 1)]}
{'root': '-51/121', 'x': '-8621/6171', 'mult': 1, 'betas': [('-1/2', 1)]}
{'root': '25/17', 'x': '457/425', 'mult': 4, 'betas': [('-4', 1), ('-5', 1), ('3', 1), ('1', 1)]}
{'root': '17/25', 'x': '457/425', 'mult': 5, 'betas': [('-4', 1), ('-5', 1), ('-6', 1), ('3', 1), ('1', 1)]}
{'root': '0', 'x': None, 'mult': 10, 'betas': []}
{'root': 'irreducible deg 2', 'mult': 2}
{'root': 'irreducible deg 2', 'mult': 2}
```

The order of the lines differs between the two primes; the content is the same. The quadratic factors, from
`debug/b2_charpoly.py`, are the same on both primes:

- lam^2 + (13550/2057) lam + 625/289, mult 2 (product of roots (25/17)^2)
- lam^2 + (9214/3025) lam + 289/625, mult 2 (product of roots (17/25)^2)

Their roots are real and negative, and map to x < -1.

## Reading

- **B2-2.** Integer beta occurs only at 17/25 and 25/17, both mapping to x* = 457/425. The other rational roots map to
  two points with x < -1, and each has only beta = -1/2, which corresponds to an integer alpha = 0, a log/pole-type
  term. The root 0 (multiplicity 10) carries no exponents in this analysis.
- **B2-3.** The contact exponent beta = -4, i.e. alpha = 7/2, is present at 17/25, together with -5 and -6. The
  exponents 3 and 1, i.e. alpha = -7/2 and -3/2, are also present at both roots.
- **Open: the two extra exponents.** The recurrence annihilates every solution of the 41-dimensional cyclic family. In
  BU the analogous family was larger than the target sequence's own module by one solution. So beta = 3 and beta = 1
  may belong to solutions that the actual coefficients do not contain. Settling that would need either the minimal
  recurrence of the data itself or a large-L numerical test of F_L. Both use data, and neither was done.

## Files

- Code in `debug/`: `b2_engine.py` (modules, IBP, atoms, closure, output), `b2_derive.py` (system and recurrence),
  `b2_data.py` (check data), `b2_check.py`, `b2_charpoly.py`, and this memo.
- Data in `debug/data/`:
  - `b2_recurrence_<p>_4200.json`, in the format b2_asym.py reads
  - `b2_asym_<p>.txt`
  - `b2_check_seq_<p>.json`
  - `b2_point_check.log`, `b2_rec_4200.log`, `b2_rec_4200_p2.log`, `b2_data.log`
