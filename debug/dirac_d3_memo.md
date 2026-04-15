# Track D3 Memo — Dirichlet Search for a Dirac Analog of F = π²/6

**Date:** 2026-04-15  
**Agent:** Worker (Track D3 of Dirac-on-S³ Tier 1)  
**Method:** Pure sympy symbolic (no numerical evaluation beyond float sanity checks).

## Reference identity (Phase 4F, α-J)

On the scalar Fock lattice of S³, with degeneracy g_m = m² at Fock level
m = n_CH + 1 ≥ 1:

```
D_{n²}(s) = Σ_{m≥1} m² · m^{-s} = ζ(s - 2)
D_{n²}(s = d_max = 4) = ζ(2) = π²/6 = F      (EXACT sympy equality)
```

The weight m² is the Fock scalar degeneracy; the exponent s = 4 is the
packing dimension d_max (Paper 0). Both are graph-intrinsic.

## Candidate series tested (all exact sympy closed forms)

With Dirac degeneracy `g_m^Dirac = 2m(m+1)` (Fock indexing, m = n_CH + 1):

| # | Series | Closed form (symbolic in s) |
|---|--------|-----------------------------|
| a | D_dirac_fock = Σ 2m(m+1)/m^s | **2·ζ(s−2) + 2·ζ(s−1)** |
| b | D_dirac_spec = Σ 2m(m+1)/(n+3/2)^s (CH-eigenvalue denom) | 2·ζ(s−2, 3/2) − (1/2)·ζ(s, 3/2) |
| c | Weyl (half) | ζ(s−2) + ζ(s−1) |
| d | Dirac − 2·scalar = Σ 2m/m^s | 2·ζ(s−1) |
| e | D_B_dirac = Σ m/m^s | ζ(s−1) |

## Evaluations at s ∈ {3, 4, 5, 6}

**Primary candidate (a), D_dirac_fock:**
- s = 3: diverges (ζ(1) pole)
- **s = 4: 2F + 2·ζ(3) = π²/3 + 2·ζ(3)** — mixed with an independent transcendental
- s = 5: 2·ζ(3) + π⁴/45
- s = 6: π⁴/45 + 2·ζ(5)

**(b) Spectral / Hurwitz form** — no reduction to π² · rational at any tested s.

**(d,e) Weight-m series** hit F = ζ(2) trivially at s = 3 (since Σ m·m^{−3} = Σ m^{−2} = ζ(2)). This is the weight-1 / dim(S³) = 3 coincidence, NOT a Dirac-specific lift: any Dirichlet series with weight m^1 gives ζ(2) at s = 3.

## F-extraction tests

Two "hits" on F emerge, both reducible:

1. **(½)·D_dirac_fock(4) − ζ(3) = F** — exact, but requires hand-subtraction of the independent transcendental ζ(3). Not a clean D_{Dirac}(d_max) = rational · F identity.

2. **D_dirac_fock(4) − 2·D_B_dirac(4) = 2F** — exact. But algebraically this is `2[2m(m+1)/m⁴] − 2[2m/m⁴] = 2·(2m²/m⁴) = 2·D_{n²}(4) = 2F`, i.e. subtracting the "m" part of 2m(m+1) recovers the scalar n²-weighted series. This is not a new Dirac identity; it merely strips the Dirac-specific 2m content and returns to the scalar case.

## Structural reason F does not lift

The Dirac degeneracy `g_m = 2m(m+1) = 2m² + 2m` is a **mixture** of two homogeneity classes:
- The **m² part** produces ζ(s−2) — the scalar F channel
- The **2m part** produces 2·ζ(s−1) — an independent channel with a different "effective dimension"

At s = d_max = 4, these evaluate to ζ(2) and 2·ζ(3). ζ(2) and ζ(3) are algebraically independent (Apéry; in any case not Q-linearly related), so no rational combination of the Dirac series collapses to a pure rational multiple of F. To isolate F from D_dirac_fock one must subtract the ζ(s−1) channel by hand — which is equivalent to projecting out the 2m part of g_m, recovering the scalar lattice.

## Verdict for D5

**F does NOT lift to the Dirac sector on S³.**

Phase 4F's identification F = D_{n²}(d_max) = ζ(s−2)|_{s=4} is *specific to the scalar Fock degeneracy m²*. The Dirac-degeneracy analog D_dirac_fock(s = d_max) = 2F + 2·ζ(3) mixes F with an independent transcendental, and no natural subtraction isolates F without reducing the Dirac series to its m²-subpart (i.e. back to the scalar case). The triviality of `D_B_dirac(3) = ζ(2)` arises from weight m at s = dim(S³), which is neither Dirac-specific nor tied to d_max and holds for any m-weighted series.

**Implication for Paper 2 / K-rule:** F entering K = π(B + F − Δ) is structurally tied to the scalar lattice's n² degeneracy at s = d_max; the Dirac sector introduces ζ(3) as an additional independent transcendental rather than reproducing F. If the Dirac sector enters the alpha construction (per Phase 4H's Δ⁻¹ = g_3^Dirac = 40 identification), it contributes ζ(3) at the same Dirichlet exponent, not F.

## Deliverables

- `debug/dirac_d3_dirichlet.py` — sympy driver (all five candidate series, four s-values, extraction tests)
- `debug/data/dirac_d3_dirichlet.json` — full symbolic table
- `debug/dirac_d3_memo.md` — this file
