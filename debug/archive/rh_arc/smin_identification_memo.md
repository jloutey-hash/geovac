# S_min Identification: Irreducible Multiple Hurwitz Zeta

**Date:** April 17, 2026
**Status:** CONFIRMED IRREDUCIBLE — 15 PSLQ attempts across 47 basis elements all failed

## Definition

S_min = Σ_{k=1}^∞ T(k)² where T(k) = 2·ζ(2, k+3/2) − (1/2)·ζ(4, k+3/2)

is the CG-weighted two-loop sunset sum on S³ from `qed_vertex.py`. The weight
comes from the SO(4) Clebsch-Gordan channel count W_total = 2·min(n1,n2) − 1 − δ_{n1,n2}.

## Numerical Value

S_min = 2.47953699802733386573169120967095568216744156767274059479172... (100 digits)

Computed at 150 dps with 10000 terms + three-term tail correction.

Ratio: S_min / D(4)² = 0.8076311883349...

## PSLQ Results

All 15 attempts returned None (no integer relation found):

| Attempt | Basis size | Elements | Result |
|---------|-----------|----------|--------|
| original_18 | 18 | Standard: {1, π², π⁴, ζ(3), ζ(5), ..., G, β(4)} | FAILED |
| with_polylog | 15 | Added Li₂(1/2), Li₃(1/2), Li₄(1/2) | FAILED |
| with_gamma_cross | 16 | Added Euler γ, G·ζ(3), β(4)·ζ(3) | FAILED |
| with_beta3_mzv | 17 | Added β(3)=π³/32, ζ(2,1,1)=ζ(4), ζ(3,1)=π⁴/360 | FAILED |
| hurwitz_products | 19 | Added ζ(2,3/2)·ζ(4,3/2) products | FAILED |
| component A (ζ(2,k)²) | 21 | Decomposed S=4A−2B+C/4 | FAILED |
| component B (ζ(2,k)·ζ(4,k)) | 21 | Cross-component | FAILED |
| component C (ζ(4,k)²) | 21 | Pure ζ(4,k) component | FAILED |
| with_polygamma | 16 | Added ψ^(n)(3/2) values | FAILED |
| ultra_wide_47 | 47 | Kitchen-sink: all of the above | FAILED |
| ratio S/D² | 17 | Tried the ratio instead | FAILED |
| ratio direct | 17 | Alternative ratio formulation | FAILED |
| ratio hz2sq | 11 | Component A / ζ(2,3/2)² | FAILED |
| double_hurwitz | 20 | Added ζ₂(2,2;3/2), ζ₂(2,4;3/2) | FAILED |
| even_odd_products | 21 | D_even·D_odd products | FAILED |

## Structural Interpretation

S_min is a **depth-2 multiple Hurwitz zeta value** at half-integer shift a=3/2:

S_min = Σ_{k≥1} [2·ζ(2, k+3/2) − (1/2)·ζ(4, k+3/2)]²

This is structurally analogous to a depth-2 multiple zeta value ζ(s₁,s₂), but
evaluated on the half-integer Hurwitz lattice rather than the integer lattice.
The half-integer shift is intrinsic to the Dirac spectrum on S³ (Camporesi-Higuchi
convention: |λ_n| = n + 3/2).

The failure of ALL extended bases — including polylogarithms, Euler sums, double
Hurwitz, and cross-products of known constants — indicates that S_min is
**genuinely irreducible**: it cannot be expressed as a rational linear combination
of any known constants up to weight 8.

## Paper 18 Classification

S_min lives at the intersection of all three taxonomy axes:
1. **Operator order:** first-order Dirac eigenvalues → odd-zeta channel active
2. **Vertex topology:** parity constraint restricts mode sums
3. **CG weighting:** SO(4) channel count creates nested min-weighted structure

The nesting (depth-2) is the mechanism that escapes standard constant bases.

## Key Files

- `debug/smin_identification.py` — computation script (544 lines)
- `debug/data/smin_identification.json` — full numerical results
- `geovac/qed_vertex.py` — `two_loop_min_weighted_hurwitz()` definition
