# Paper 60 σ-basis SW conditioning: independent extended-precision cross-check (2026-10-09)

**Verdict: GO.** Every number Paper 60 quotes for the full σ basis (sec:molecular plateau paragraph
and the "[MEASURED] The full σ basis" paragraph) is reproduced to its printed digits by an
independent exact-arithmetic route at 60 digits. The large condition numbers are real. They are
not artefacts of float64 or the quadrature. No paper number has to change.

Driver `debug/p60_sigma_precision_crosscheck.py`; data `debug/data/p60_sigma_precision_crosscheck.json`
(+ `.log`); wall time 31 s; exit 0.

## Independent route

Every σ Coulomb–Sturmian on centre A or B has the form poly(ξ,η)·e^{−k r_c} in prolate coordinates,
because r_c = h(ξ±η), z_c = h(ξη±1), and r^l P_l(z/r) is a polynomial in (z, r²). The factor
dV/r_c = 2πh²(ξ∓η) cancels the singularity exactly. Each S and W entry is therefore a finite sum
of products X_p(α)·Y_q(b), where X_p = ∫₁^∞ ξ^p e^{−αξ} and Y_q = ∫₋₁¹ η^q e^{−bη}, with b = α, 0 or −α
for AA, AB and BB. Both are computed in closed form or by convergent series in mpmath at dps 60.
The polynomial coefficients are exact (built from the Legendre and Laguerre coefficients), there is
no quadrature, and eigenvalues come from `mp.eigsy` at the same precision.

Three checks on this route:
- Intra-centre S block equals the identity to 2.3e-56 (Sturmian orthonormality with 1/r weight).
- Third route: one hard cross entry, S[A(4,3), B(6,5)] at kR=3, by mpmath 2-D adaptive quad =
  0.36990981745016… This matches the closed form to every digit; float64 gives 0.369909817450969.
- For N=6, kR=3 (42×42), the largest entry discrepancy against the float64 production builder is
  |ΔS| = 2.1e-12 and |ΔW| = 1.3e-12. The float route's noise floor is therefore about 1e-12
  absolute. This is consistent with its own intra-block deviation of 1.7e-12, and plausibly comes
  from the w·e^{t} Gauss–Laguerre weight product.

## Quoted vs independent

| Quantity (paper) | Quoted | Independent (dps 60) | float64 (production) |
|---|---|---|---|
| cond(S), kR=3, 3 fns/centre (N=2) | 18 | 18.262 | 18.262 |
| ″ 6 (N=3) | 5.1e3 | 5.1071e3 | 5.1071e3 |
| ″ 10 (N=4) | 1.2e7 | 1.2305e7 | 1.2305e7 |
| ″ 15 (N=5) | 1.1e11 | 1.1057e11 | 1.1057e11 |
| gerade cond, kR=3, N=2 | 2.3 | 2.3167 | 2.3167 |
| ″ N=3 | 21 | 21.121 | 21.121 |
| ″ N=4 | 7.4e3 | 7.3966e3 | 7.3966e3 |
| ″ N=5 | 2.1e7 | 2.1009e7 | 2.1009e7 |
| Bessel deficit kR=1, N=2..6 | 0.158, 0.090, 0.064, 0.048, 0.037 | 0.15791, 0.09025, 0.06378, 0.04791, 0.03686 | identical to 5 dp |
| soft eigs (<1e-2), kR=3, N=2..6 | 0,1,3,5,8 | 0,1,3,5,8 | 0,1,3,5,8 |
| H2+ R=2 Galerkin err, 3 fns | 8.5e-4 | 8.48403e-4 | 8.48403e-4 |
| ″ 6 fns | 9.1e-5 | 9.12013e-5 | 9.12013e-5 |
| ″ 15 fns | 5.4e-6 | 5.43575e-6 | 5.43575e-6 |

True H2+ energies at R=2 (Galerkin, fixed-R isoenergetic posing, no orthogonalisation tolerance):
−1.10178581074458 (N=2), −1.10254301269743 (N=3), −1.10262877824858 (N=5).

**Smallest eigenvalues, kR=3, N=5** (mp vs float):
1.80877e-11 / 1.8088e-11, 9.51994e-8 / 9.5199e-8, 3.08115e-6 / 3.0812e-6, 1.77710e-4 / 1.7771e-4.
They agree to at least 4 significant figures, including λ_min. The N=5 cond quoted in the paper is
inside float64's range: entry noise of 1e-12 against λ_min of 1.8e-11 shifts it by less than 1e-4 relative.

## Where float64 is at its floor (none of these are quoted in the paper)

| Case | float64 cond | true cond |
|---|---|---|
| kR=3, N=6 (21/centre) | 2.458e15 (λ_min 8.14e-16) | **2.645e15** (λ_min 7.561e-16) |
| kR=1, N=4 | 3.633e13 | 3.645e13 |
| kR=1, N=5 | −4.1e15 (λ_min negative) | **2.83e19** |
| kR=1, N=6 | −5.3e15 | **5.72e25** |
| gerade, kR=1, N=6 | 5.0e16 | 8.33e19 |
| gerade, kR=3, N=6 | 2.1052e11 | 2.1053e11 |

- The paper's kR=3 sentence stops at 15 functions/centre, and that endpoint is float-safe.
- `debug/data/p60_sigma_basis_sw.json` records the kR=1, N≥5 and kR=3, N=6 conditions. Those
  entries are float-floor values; the kR=1 ones are even negative. Anyone who later quotes them
  should use the true values above.
- The soft-direction counts and the Bessel deficits are unaffected at every N, including N=6,
  because neither depends on the smallest eigenvalue's digits.
- The qualitative "smallest eigenvalue falling faster than geometrically" holds on exact values.
  At kR=3, λ_min = 1.04e-1, 3.92e-4, 1.63e-7, 1.81e-11, 7.56e-16, so the step ratios are 265, 2400,
  9000 and 24000.

## Canonical-orthogonalisation tolerance

Float Galerkin error at R=2 against CANON_TOL ∈ {0, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6}:
- **N=2 and N=3:** bit-identical at every tolerance.
- **N=5:** identical for tolerance ≤ 1e-12 (5.435751421e-6). The error moves by 3.6e-15 Ha at 1e-10
  and by 6.2e-9 Ha at 1e-6.

The mp solve uses no tolerance at all and gives 5.435751416e-6. The energies are therefore
insensitive to CANON_TOL = 1e-12 at every printed digit.

## Optional hygiene (not required)

`tests/test_paper60_sigma_basis.py` and the paper use only float-safe cases (kR=3, N≤5, and N=6
only for soft counts and deficits), so nothing is owed. If the paper ever extends the kR=3 ladder to
21 functions per centre, it should quote 2.6e15, not the float value 2.5e15.
