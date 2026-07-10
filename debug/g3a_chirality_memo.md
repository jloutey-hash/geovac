# G3-A: GeoVac-side Chirality Grading γ_GV

**Sprint:** G3 (Electroweak chirality co-location on S³), Track A (GeoVac side).
**Status:** CLOSED, POSITIVE.
**Date:** 2026-05-06.
**Module:** `geovac/chirality_grading.py` (~310 lines).
**Tests:** `tests/test_chirality_grading.py` (50/50 passing).
**Data:** `debug/data/g3a_chirality.json`.

## TL;DR

The natural NCG chirality grading on the GeoVac full-Dirac truncated
operator system is

> **γ_GV = σ_x ⊗ I_{spinor_dim(n_max)}**

in the basis ordering `[chi=+1 block, chi=-1 block]` from
`full_dirac_basis`. It satisfies all three Z_2-grading axioms
(γ² = I, {γ, D_truthful} = 0, [γ, M] = 0 for every scalar multiplier
M ∈ O) **bit-exactly at n_max ∈ {1, 2, 3}**, and γMγ = M holds as a
matrix equality (not merely as a span containment). No SDP is involved
on the GeoVac side; the construction is direct linear algebra.

The diagonal candidate γ_diag = σ_z ⊗ I (the "D-eigenvalue sign")
satisfies γ² = I and γMγ = M but **commutes** with D_truthful instead
of anticommuting; it is therefore not the chirality grading, only an
auxiliary Z_2 sign. We expose it as a diagnostic to make the negative
control explicit.

For G3-C the recommendation is unambiguous: **use γ_GV = σ_x ⊗ I with
the truthful CH Dirac.** The offdiag CH Dirac (R3.5 device for the
Connes-distance SDP) is *not* compatible with γ_GV in the
{γ, D} = 0 sense, and that is consistent with the existing reading
that offdiag is an SDP-bounding tool, not a spectral-action object.

## Convention chosen and why

### Setup

In `geovac/full_dirac_operator_system.py`, the full-Dirac sector is
the Weyl basis doubled by an explicit chirality label:

- `full_dirac_basis(n_max)` orders all chi=+1 labels first, then all
  chi=-1 labels (each block has length spinor_dim(n_max)).
- `camporesi_higuchi_full_dirac_matrix(basis)` is diagonal with
  eigenvalue chi · (n_fock + 1/2). In block form:

  D_truthful = σ_z ⊗ D_+,    where D_+ = diag(n_fock + 1/2 on Weyl basis).

- Scalar multipliers M_full = M_Weyl ⊕ M_Weyl (block-diagonal). In
  block form:

  M_full = I_2 ⊗ M_Weyl.

### Three Z_2-grading axioms

The chirality grading γ on a spectral triple (Connes 1995,
van Suijlekom 2015 §3) must satisfy:

| Axiom  | Statement | Geometric meaning |
|--------|-----------|------------------|
| (i)    | γ² = I    | Z_2 grading      |
| (ii)   | {γ, D} = 0 | γ swaps left/right components of D |
| (iii)  | [γ, a] = 0 ∀ a ∈ A | γ commutes with the algebra (not with D) |

### σ_x ⊗ I satisfies all three

| Quantity | Calculation | Result |
|----------|-------------|--------|
| γ²       | (σ_x ⊗ I)² = σ_x² ⊗ I = I_2 ⊗ I  | I  |
| {γ, D}   | (σ_x ⊗ I)(σ_z ⊗ D_+) + (σ_z ⊗ D_+)(σ_x ⊗ I) = (σ_x σ_z + σ_z σ_x) ⊗ D_+ = 0 ⊗ D_+ | 0 |
| [γ, M]   | (σ_x ⊗ I)(I_2 ⊗ M_Weyl) − (I_2 ⊗ M_Weyl)(σ_x ⊗ I) = σ_x ⊗ M_Weyl − σ_x ⊗ M_Weyl | 0 |
| γMγ      | (σ_x ⊗ I)(I_2 ⊗ M_Weyl)(σ_x ⊗ I) = σ_x² ⊗ M_Weyl = I_2 ⊗ M_Weyl | M  |

In particular γMγ = M holds as a *matrix equality*, the strongest
possible "γ preserves O" statement: γ acts trivially on the operator
system. There is no nontrivial action of γ on M, but there is a
nontrivial action on D, which is the structurally informative outcome.

### σ_z ⊗ I does NOT satisfy (ii)

Computing:

  {σ_z ⊗ I, σ_z ⊗ D_+} = σ_z² ⊗ D_+ + σ_z² ⊗ D_+ = 2 (I_2 ⊗ D_+)
                       = 2 |D_truthful|.

So σ_z ⊗ I anticommutes with the *off-diagonal* part of D (which is
zero in the truthful case) but commutes with the diagonal part. This
is the "D-eigenvalue sign" Z_2, not the chirality grading. We retain
σ_z as a diagnostic only.

The audit table below (n_max = 1, 2, 3) confirms numerically.

## Residual table

All values are absolute residuals; the spectrum lives in [0, dim_H]
so any value below 1e-14 is bit-exact zero in float128 arithmetic.

### Convention σ_x (NCG chirality)

| n_max | dim_H | dim_O | \|gen\| | γ²−I | {γ, D_truthful} | {γ, D_offdiag} | γMγ − M (worst) | [γ, M] (worst) |
|------:|------:|------:|--------:|-----:|---------------:|---------------:|----------------:|---------------:|
| 1     | 4     | 1     | 1       | 0    | 0              | 1.0e-2         | 2.5e-16         | 0              |
| 2     | 16    | 14    | 14      | 0    | 0              | 2.0            | 4.9e-16         | 0              |
| 3     | 40    | 55    | 55      | 0    | 0              | 2.0            | 1.3e-14         | 0              |

The "γMγ − M (worst)" column reports `op_sys.contains(γMγ)` residual,
which checks span-membership. The exact equality `γMγ = M` is verified
separately (test `test_sigma_x_gMg_equals_M_exactly`); the small
nonzero values shown are floating-point noise in lstsq, not a real
deviation.

### Convention σ_z (D-eigenvalue sign — diagnostic)

| n_max | dim_H | dim_O | \|gen\| | γ²−I | {γ, D_truthful} | {γ, D_offdiag} | γMγ − M (worst) | [γ, M] (worst) |
|------:|------:|------:|--------:|-----:|---------------:|---------------:|----------------:|---------------:|
| 1     | 4     | 1     | 1       | 0    | 3.0            | 3.01           | 2.5e-16         | 0              |
| 2     | 16    | 14    | 14      | 0    | 5.0            | 5.23           | 4.9e-16         | 0              |
| 3     | 40    | 55    | 55      | 0    | 7.0            | 7.45           | 1.3e-14         | 0              |

The σ_z anticommutator with D_truthful is exactly 2 · (n_max + 1/2)
= 2n_max + 1, the magnitude of the top Dirac eigenvalue (paired with
its antichiral partner): n_max=1 → 3, n_max=2 → 5, n_max=3 → 7.
This is the signature of "{γ, D} = 2D" rather than "{γ, D} = 0".

## {γ, D_offdiag} ≠ 0: honest disclosure

The σ_x convention does NOT anticommute with the offdiag CH Dirac
used in WH1-R3.5 / Sprint-TS R3.5 for the Connes-distance SDP. The
nonzero residuals at n_max = 2, 3 are 2.0 in absolute units (one
unit per E1 ladder coupling, of which there are O(n_max²) per
chirality). At n_max = 1 the residual is 1e-2 because the offdiag
generator only adds the small l_lift, m_lift diagonal perturbations
(no E1 partners exist at n_max = 1).

### Why this happens

The offdiag CH Dirac has the structure

  D_offdiag = D_diag + ΔD_within + ΔD_cross,

where:

- D_diag is diagonal with chi · (n_lift · (n + 1/2)) + l_lift · l + m_lift · 2m_j.
  The chirality-odd part anticommutes with σ_x, the chirality-even
  parts (l_lift, m_lift) commute with σ_x.
- ΔD_within has E1 selection-rule entries (Δn = ±1, Δl = ±1, |Δm_j| ≤ 1)
  with strength `offdiag_alpha` *equal on both chirality blocks*. So
  ΔD_within = I_2 ⊗ X_within in block form, which **commutes** with
  σ_x ⊗ I rather than anticommuting.
- ΔD_cross has chirality-flipping entries with strength
  `chirality_coupling`. Of cross-chirality contributions, those with
  `chirality_coupling = 1` (default) appear at the off-diagonal
  Pauli structure σ_x ⊗ Y, which commutes with σ_x ⊗ I.

So the "natural" offdiag construction is built to be approximately
chirality-symmetric, not chirality-antisymmetric. This is a feature
of the SDP-bounding device, not a bug: the offdiag perturbation was
designed (in R3.2 / R3.5) to lift the n-degeneracy of the Connes
distance kernel, NOT to be a J- or γ-symmetric Dirac.

The offdiag CH Dirac is consistent with the existing CLAUDE.md
reading (Sprint TS-C #14, real_structure_finite_nmax_memo.md):
"offdiag is an SDP-bounding device for Connes distance, NOT
spectral-action-foundational." γ_GV inherits the same boundary.

### Implication for G3-C

For the tensor-product test (γ_GV ⊗ I_F) ≟ (I_GV ⊗ γ_F), G3-C should
use γ_GV with the **truthful CH Dirac**. The offdiag CH Dirac is
spectrally informative for Connes-distance work but does not respect
the γ-grading by construction.

## What γ_GV does (and does not) entail for the AC extension

γ_GV is the chirality grading on the GeoVac side of the
almost-commutative extension A_GV ⊗ A_F (Sprint H1, Paper 32 §VIII.C).
It plays two structural roles in inner-fluctuation theory:

1. **It splits the spinor bundle into left/right.** The +1 eigenspace
   of γ_GV at finite n_max is exactly the chi=+1 block of
   `full_dirac_basis`, which is the Weyl spinor sector
   (`spinor_operator_system.py`). The −1 eigenspace is the anti-Weyl
   sector. Together they reconstruct the full Dirac sector with
   γ_GV-grading the Z_2 split.

2. **It is the source of the chirality factor in inner fluctuations.**
   In Connes-Marcolli, the inner fluctuation D → D + ω + γ ω γ
   contains both gauge-1-form and Higgs contributions, and γ is what
   distinguishes them. With γ_GV = σ_x ⊗ I and ω built from scalar
   multiplers, γ_GV ω γ_GV = ω exactly (since γMγ = M for our O), so
   the gauge-1-form sector is unmodified. The Higgs sector (Sprint H1
   verdict POSITIVE-THIN) emerges instead from the off-diagonal
   D_F block on the AC factor side; γ_GV's role there is to multiply
   the cross-block contribution by the appropriate sign in the
   tensor-product spectral triple — that is exactly the G3-C question.

γ_GV does NOT, on its own, force any particular Higgs structure or
Yukawa selection. The Sprint H1 verdict that GeoVac admits Higgs but
does not autonomously emit Y is unchanged by G3-A.

## Comparison to J on truthful CH (Step 2 of three-step PI commitment)

Sprint Connes-Step-2 (`debug/real_structure_finite_nmax_memo.md`,
2026-05-05/06) verified that the real structure J on truthful CH
satisfies the three load-bearing Connes axioms (J² = -I, JD = +DJ,
J·O·J⁻¹ = O) **exactly** at n_max ∈ {1, 2, 3}, and that JD = +DJ
fails on offdiag CH (residual 2.0). The same residual structure
appears here for the chirality grading γ:

| Axiom              | Truthful CH residual | Offdiag CH residual |
|--------------------|---------------------:|--------------------:|
| J² = −I            | 0 exact              | (J unchanged)       |
| JD = +DJ           | 0 exact              | 2.0 (Step 2)        |
| γ² = I             | 0 exact              | (γ unchanged)       |
| {γ, D} = 0         | 0 exact              | 2.0 (this sprint)   |
| γ·O·γ ⊆ O          | 0 exact              | 0 exact             |
| J·O·J⁻¹ ⊆ O        | 0 exact              | 0 exact             |

The pattern is consistent: **truthful CH respects J and γ; offdiag CH
respects neither.** Both J and γ are spectral-action-foundational
objects, and the offdiag construction sacrifices both axioms in
exchange for SDP-tractable Connes distances. This sharpens the WH1
reading: offdiag is a metric/SDP tool, truthful CH is the
spectral-triple object.

## Recommendation for G3-C

Use γ_GV with `convention="sigma_x"` and the **truthful CH Dirac**:

```python
from geovac.chirality_grading import build_gamma_GV
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
)
from geovac.almost_commutative import build_F_factor  # G3-B

n_max = 3
gamma_GV = build_gamma_GV(n_max, convention="sigma_x").matrix
op_GV = FullDiracTruncatedOperatorSystem(n_max)
D_GV = camporesi_higuchi_full_dirac_matrix(op_GV.basis)

# G3-C tensor product
F = build_F_factor()  # delivers gamma_F, dim_F, A_F, etc.
H_total_dim = op_GV.dim_H * F.dim_F
gamma_GV_x_I = np.kron(gamma_GV, np.eye(F.dim_F))
I_GV_x_gamma_F = np.kron(np.eye(op_GV.dim_H), F.gamma_F)
residual = gamma_GV_x_I - I_GV_x_gamma_F   # operator on H_GV ⊗ H_F
```

The G3-C diagnostic is then the operator residual `residual` (or
appropriately physical-sector-restricted version) at n_max = 2, 3.

## Files

- `geovac/chirality_grading.py` (310 lines).
- `tests/test_chirality_grading.py` (50/50 tests passing).
- `debug/data/g3a_chirality.json` (full audit at all n_max × both
  conventions).
- `debug/g3a_chirality_memo.md` (this file).
- Pre-existing tests still pass: `test_full_dirac_operator_system.py`
  (29/29 + 2 skip), `test_spinor_operator_system.py` (29/29 + 1 skip),
  `test_real_structure.py` (44/44 + 1 skip).

## Open items

- G3-A complete and self-contained.
- G3-B (γ_F audit on AC factor) and G3-C (tensor-product
  identification) are separate sprint tracks; this memo's
  recommendation feeds G3-C directly.
- Optional follow-up: build a γ_GV-symmetric **offdiag** CH Dirac
  variant (anticommutator-respecting offdiag couplings; opposite-sign
  E1 ladders on the two chirality blocks) for joint Connes-distance +
  spectral-action work. This would be a new device, not a fix to the
  current offdiag.

## Status verdict

**G3-A POSITIVE.** The natural NCG chirality grading γ_GV = σ_x ⊗ I
exists at finite n_max on the GeoVac full-Dirac sector and satisfies
all three Z_2 axioms bit-exactly with the truthful CH Dirac. Provides
the GeoVac-side input for G3-C (tensor-product chirality co-location
test).
