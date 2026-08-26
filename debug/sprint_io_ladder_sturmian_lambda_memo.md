# I/O ladder — Sturmian-λ check: λ does NOT preserve on the genuine Coulomb-Sturmian basis

**Date:** 2026-08-17 · **Type:** diagnostic (no geovac/ / papers / CHANGELOG / version edits)
**Driver:** `debug/io_ladder_sturmian_lambda.py` · **Data:** `debug/data/io_ladder_sturmian_lambda.json`
**Follows:** `sprint_io_ladder_lambda_sweep_memo.md` (LG-3, hydrogenic — favorable).

## Why this run
The LG-3 λ result was measured on GeoVac's **hydrogenic** builder. To put the headline
"on James's basis" I rebuilt λ in the **genuine shared-k Coulomb-Sturmian** basis. Controlled
s-sector atomic build (He, Z=2): only the radial convention is flipped —
hydrogenic R_n0 (a=Z/n, L²-orthonormal) vs Coulomb-Sturmian S_n0 (shared scale k, L²-NON-
orthogonal). Anchored k=Z so n=1 is identical (fair). Same code both families: L² overlap S →
Löwdin S^{-1/2} → transform h1 + ERI tensor → JW → λ=Σ|c_i| (excl identity). Validated:
F⁰(1s,1s)=1.252 (exact 5Z/8=1.250, grid), ⟨1s|h|1s⟩=−2.000 (exact −Z²/2).

## Result — the genuine Sturmian basis INFLATES λ
| N | Q | λ hydrogenic | λ Sturmian | cond(S) Sturmian |
|--:|--:|--:|--:|--:|
| 1 | 2 | 1.687 | 1.687 | 1.0 |
| 2 | 4 | 2.248 | 5.669 | 3.0 |
| 3 | 6 | 3.571 | 24.27 | 5.8 |
| 4 | 8 | 5.101 | 60.73 | 9.5 |
| 5 | 10 | 6.679 | 119.7 | 13.9 |

**Scaling:** hydrogenic λ~**Q^1.19** (flat/healthy); genuine Sturmian λ~**Q^3.33** — inflates,
**17.9× worse at Q=10**.

## Mechanism — the cost-conservation wall, in the λ currency
The Coulomb-Sturmian basis is L²-non-orthogonal (orthonormal only in the 1/r weight); at fixed
k the higher-n functions grow linearly dependent, so cond(S) climbs (1→14 over N=1..5). A
second-quantized qubit encoding needs an L²-orthonormal single-particle basis, so it must
Löwdin-orthogonalize — and Löwdin of an ill-conditioned Gram spreads the coefficients, inflating
λ. This is the **same wall** the corpus documented as "Löwdin retrofit → 17.9× Pauli" and
"the overlap metric never disappears" (composition-wall / two-kinds-of-sparsity) — now measured
in λ. The 17.9× echo is a striking (if numerically coincidental) confirmation, not a new wall.

## What it means (the honest reframe)
- **GeoVac's ACTUAL encoding is hydrogenic** (per-center, a=Z/n), forced by the Papers 8–9
  dual-p₀ theorem (no shared-k Sturmian exists for a heteronuclear pair). Hydrogenic is
  L²-orthonormal ⇒ no Löwdin penalty ⇒ the earlier **LG-3-favorable λ result HOLDS for GeoVac's
  real qubit Hamiltonians** (Paper 14 / ecosystem). This run explains *why* hydrogenic is the
  right radial encoding: not only dual-p₀, but λ-health.
- **The genuine shared-k Coulomb-Sturmian basis — James's basis — does NOT preserve λ** under
  standard second-quantized Löwdin encoding. The cleanest angular-generation basis is the one
  whose non-orthogonality costs λ.

## Scope / what stays open (do not overclaim)
- s-only, He, k=Z anchor, 5 shells. l>0 / other k / molecular could shift magnitudes (not the
  qualitative direction — cond(S) growth is intrinsic to the shared-k basis).
- Measured for **second-quantized JW + Löwdin**. James's actual lever (Rung 3) is
  **first-quantized** simulation, which does not second-quantize into orthonormal spin-orbitals
  and so may not pay the Löwdin λ penalty at all. **The 2Q negative may not transfer to 1Q** —
  untested, and the precise open question for the call.

## Verdict — LG-3 is basis-specific
- Hydrogenic (GeoVac's real encoding, 2Q): **FAVORABLE** (λ~Q^1.19).
- Genuine Sturmian (2Q + Löwdin): **NEGATIVE** (λ~Q^3.33, non-orthogonality cost).
- Genuine Sturmian (first-quantized): **OPEN** — the measurement to settle it, and the natural
  James topic.

## Cross-refs
`sprint_io_ladder_{accounting,radial_seeds,costmodel,lambda_sweep}_memo.md`;
composition-wall / cost-conservation (`sprint_commutator_and_explorer_memo.md`).
