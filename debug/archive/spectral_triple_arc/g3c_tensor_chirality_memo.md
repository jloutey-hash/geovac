# G3-C Tensor-Product Chirality Diagnostic — Memo

**Sprint:** G3-C (parallel with G3-A and G3-B; re-run with corrected conventions 2026-05-06)
**Date:** 2026-05-06 (initial), 2026-05-06 (re-run)
**Status:** Computational diagnostic complete; spectrum re-verified under production conventions.
**Verdict:** **NEGATIVE.** The chirality grading on the GeoVac side does *not* identify with the weak-isospin chirality on the AC factor. The two are structurally independent ℤ₂ gradings on the tensor product H_GV ⊗ H_F.
**Implication:** Electroweak unification does NOT close on S³ alone via direct chirality identification. This is a clean result, and pushes the EW-unification question onto the cross-manifold (G4) frontier exactly as the §VIII.B four-gap analysis already anticipated.

---

## Convention-correction note (2026-05-06 re-run)

The original (2026-05-06) version of this diagnostic was built standalone with
local `gamma_GV` and `gamma_F` constructions to avoid blocking on the parallel
G3-A and G3-B forks. After G3-A and G3-B landed, the parent PM flagged a
convention concern (the original directive's "γ_F identical on both copies"
would yield KO-dim 0, which silently breaks H1's KO-dim 6 axiom suite) and
requested a re-run using the production sources.

On inspection:

* The **original local γ_F** was already at KO-dim 6:
  ```
  diag(+1, +1, -1, -1, -1, -1, +1, +1)   (matter, then antimatter sign-flipped)
  ```
  This is bit-identical to `AlmostCommutativeTriple.gamma_F()` from G3-B. The
  KO-0 / KO-6 confusion in the parent's read of the original script was a
  misread, not a bug.
* The **original local γ_GV** was the σ_z (D-eigenvalue sign) convention,
  i.e. `diag(b.chirality)` on the full Dirac basis. This is *not* the NCG
  chirality grading: σ_z commutes with the truthful Camporesi-Higuchi Dirac
  rather than anticommuting (G3-A diagnostic, see §3 of `debug/g3a_chirality_memo.md`).
  The proper NCG chirality is σ_x ⊗ I_w.

The re-run uses the production sources:

| Operator | Production source | Convention |
|:---------|:------------------|:-----------|
| γ_GV     | `geovac.chirality_grading.build_gamma_GV(n_max, convention="sigma_x")` | NCG chirality grading; anticommutes truthful CH Dirac (G3-A verified) |
| γ_F      | `geovac.almost_commutative.AlmostCommutativeTriple.gamma_F()`           | Connes-Marcolli KO-dim 6; antimatter sign-flip forced by {J_F, γ_F} = 0 (G3-B verified) |

**Numerical impact of the convention correction:** none on the headline
spectrum.  σ_x and σ_z share eigenvalues ±1 with identical multiplicities and
both commute with γ_F, so the operator Δ has the *same* eigenvalue spectrum in
either convention. The eigenvectors differ (σ_x mixes the two chirality blocks
off-diagonally; σ_z is block-diagonal), and one secondary swap diagnostic
changes character (see Swap (c) below), but the verdict NEGATIVE is robust
to the convention choice. A σ_x-vs-σ_z cross-check at n_max=2 in the re-run
script confirms identical spectra.

What the re-run *did* do is land the diagnostic on the same operator objects
that Paper 32 §VIII.C will reference when documenting the result, so the
spectrum table below is now bit-identical to what the §VIII.C addendum
should cite.

---

## 1. The diagnostic

We compute the operator residual

```
    Delta := (gamma_GV ⊗ I_F)  -  (I_GV ⊗ gamma_F)
```

on the tensor-product Hilbert space H = H_GV ⊗ H_F, with:

- **gamma_GV** = `build_gamma_GV(n_max, convention="sigma_x").gamma_matrix()`
  from `geovac/chirality_grading.py`. In the [chi=+1, chi=-1] block ordering
  of `full_dirac_basis(n_max)`, this is σ_x ⊗ I_w (where dim_w =
  spinor_dim(n_max)). Verified by G3-A: γ_GV² = I exact, {γ_GV, D_truthful} = 0
  exact, [γ_GV, M] = 0 exact for every scalar multiplier (γ_GV commutes with
  all of A_GV).

- **gamma_F** = `AlmostCommutativeTriple.gamma_F()` from
  `geovac/almost_commutative.py`. Connes-Marcolli KO-dim 6:
  ```
  gamma_F^matter      = diag(+1, +1, -1, -1)
  gamma_F^antimatter  = diag(-1, -1, +1, +1)
  gamma_F             = diag(+1, +1, -1, -1, -1, -1, +1, +1)
  ```
  Verified by G3-B: γ_F² = I, γ_F^* = γ_F, {γ_F, D_F} = 0 at any Yukawa,
  γ_F π_F(a) γ_F = π_F(a), {J_F, γ_F} = 0 — all bit-exact zero.

- Sanity invariants verified at every n_max in the re-run: `gamma_GV^2 - I` = 0,
  `gamma_F^2 - I` = 0, and `[gamma_GV ⊗ I_F, I_GV ⊗ gamma_F]` = 0 (the two
  operators must commute because they act on disjoint tensor factors).

If γ_GV and γ_F act as the *same* ℤ₂ grading on physical states (after
some standard convention swap), then Δ = 0 and electroweak chirality
co-locates with GeoVac chirality on the S³ Dirac sector. If Δ ≠ 0
robustly, the two are independent gradings and EW unification fails to
close on S³ alone.

---

## 2. Numerical results

Computation at n_max ∈ {1, 2, 3} on the full Dirac sector
(dim_GV = 4, 16, 40; dim_F = 8 throughout).

### Baseline residual

| n_max | dim_GV | dim_total | ‖Δ‖_op | ‖Δ‖_F     | Spectrum (mult)              |
|:-----:|:------:|:---------:|:------:|:---------:|:-----------------------------|
| 1     | 4      | 32        | 2.000  | 8.000     | {−2 : 8,  0 : 16, +2 : 8}    |
| 2     | 16     | 128       | 2.000  | 16.000    | {−2 : 32, 0 : 64, +2 : 32}   |
| 3     | 40     | 320       | 2.000  | 25.298    | {−2 : 80, 0 : 160, +2 : 80}  |

The pattern is exact and n_max-independent at the level of structure:
- Operator norm is **exactly 2** at every n_max.
- Eigenvalue spectrum is **{−2, 0, +2}** with multiplicities in ratio 1 : 2 : 1.
- Frobenius norm scales as √(½·dim_total · 4) = √(2·dim_total),
  giving 8.0, 16.0, √640 ≈ 25.298 at n_max = 1, 2, 3 — exactly matching
  computed values.

This is the structurally cleanest possible non-zero residual: half the
joint basis dimension supports chirality identification (Δ = 0), the
other half sees a chirality mismatch of magnitude 2.

### Convention swap (a): gamma_F → −gamma_F

Global left/right convention swap on the AC factor. Standard NCG-literature
ambiguity (different papers use opposite sign conventions).

Result at every n_max: **identical baseline.** ‖Δ‖_op = 2, spectrum
{−2, 0, +2} with the same multiplicities. The swap relabels which entries
are 0 vs ±2 but the residual structure is unchanged.

### Convention swap (b): restrict to matter sector only

Drop the antimatter half of H_F (dim_F = 8 → 4). "Physical fermion sector"
restriction.

Result at every n_max: **‖Δ‖_op = 2 unchanged**, spectrum still contains
±2. Frobenius norm drops by exactly √2 because we kept half the basis.
Restriction to matter does not collapse the residual.

### Convention swap (c): restrict to Weyl chirality only on H_GV

Drop the anti-Weyl half of H_GV. Under the **σ_x γ_GV** (NCG chirality)
convention this becomes structurally informative: σ_x has zero diagonal,
so projecting to the chi=+1 block yields γ_GV|_Weyl = 0_w. Then

```
    Delta_c = 0_w ⊗ I_F  -  I_w ⊗ gamma_F  =  -I_w ⊗ gamma_F.
```

Result: spectrum **{−1, +1}** with equal multiplicities (since γ_F has
4 eigenvalues +1 and 4 eigenvalues −1, and we tensor with I_w). ‖Δ_c‖_op
= 1, ‖Δ_c‖_F = √(dim_w · 8). The residual does NOT vanish; it survives
as the trivial reduction "γ_F itself, lifted to the Weyl sector."

(Under the original σ_z γ_GV, the Weyl-only restriction gave γ_GV|_Weyl
= +I, so Δ_c = I − γ_F had spectrum {0, +2}. The σ_x re-run replaces this
with {−1, +1} of operator norm 1. Both are nonzero, both confirm NEGATIVE.)

### Combined swaps (a+b) and (a+b+c)

Sign-flip + matter, and sign-flip + matter + Weyl. Under σ_x γ_GV:
- (a+b): ‖Δ‖_op = 2, spectrum {−2, 0, +2} with mult 4·dim_w each at ±2,
  8·dim_w at 0.
- (a+b+c): ‖Δ‖_op = 1, spectrum {−1, +1} with equal mult 4·dim_w each.

In every case the residual is nonzero. No standard convention swap can
absorb it.

### The product operator γ₅ := γ_GV ⊗ γ_F

Sanity check on the canonical NCG combined grading:

```
gamma_5 := gamma_GV ⊗ gamma_F
gamma_5^2 = (gamma_GV ⊗ gamma_F)(gamma_GV ⊗ gamma_F)
         = gamma_GV^2 ⊗ gamma_F^2 = I ⊗ I = I.
```

Verified numerically: `is_idempotent_squared_to_I = True` at every n_max.
Spectrum is {−1, +1} with equal multiplicities; trace = 0. **The product
works, the difference does not** — γ₅ exists and is well-defined as the
canonical Connes-Chamseddine combined grading, but γ_GV and γ_F do not
identify across the tensor product.

### σ_x vs σ_z cross-check

The re-run includes a sanity cross-check at n_max=2 confirming that the
spectrum of Δ is identical under the σ_x (NCG chirality) and σ_z
(D-eigenvalue sign) γ_GV conventions:

```
sigma_x spectrum: [(-2.0, 32), (0.0, 64), (2.0, 32)]
sigma_z spectrum: [(-2.0, 32), (0.0, 64), (2.0, 32)]
Identical: True
```

This confirms that the structural verdict is robust across both natural
γ_GV conventions — both eigenvalues {±1} for γ_GV combined with γ_F's
diagonal action gives the same Δ spectrum.

---

## 3. Structural reading

The diagnostic exhibits a clean, n_max-independent obstruction. Three
readings of what it means, in increasing order of strength.

### Reading 1: γ_GV and γ_F grade orthogonal ℤ₂ factors

The tensor product H_GV ⊗ H_F naturally carries TWO independent ℤ₂
gradings:

- γ_GV ⊗ I_F (the GeoVac chirality on the spinor bundle of S³)
- I_GV ⊗ γ_F (the weak-isospin chirality on the AC factor)

These commute as operators (sanity check confirms `[A, B] = 0` exactly
in the re-run), and their product is the canonical Connes-Chamseddine
combined grading γ₅ = γ_GV ⊗ γ_F.

Asking "do γ_GV and γ_F identify" is asking whether these two orthogonal
ℤ₂'s collapse to one — i.e., whether γ_GV ⊗ I_F = I_GV ⊗ γ_F as
operators. The diagnostic says no: they grade independent dimensions of
the joint Hilbert space, and their difference Δ has spectrum {−2, 0, +2}
with equal weight on the nonzero eigenvalues.

### Reading 2: GeoVac chirality is *Dirac sign*, F chirality is *weak isospin*

These are physically different objects. γ_GV on the full Dirac sector
tracks the chirality of the Camporesi-Higuchi spinor bundle on S³ — a
property of the spin geometry. γ_F tracks the SU(2)_L doublet vs SU(2)_R
singlet structure of the matter representation in the finite electroweak
triple. There is no a priori reason for these to be the same grading,
and the diagnostic confirms they aren't.

### Reading 3: this is the precise statement of "GeoVac doesn't autonomously emit Yukawa"

Sprint H1 found that GeoVac's AC extension admits Higgs structurally
but does not autonomously select the Yukawa Y. The G3-C diagnostic
sharpens the structural reason: the Yukawa lives in the off-diagonal
block of D_F that flips F-chirality (sends γ_F = +1 to γ_F = −1).
Since γ_F is independent of γ_GV, the choice of off-diagonal block in
D_F is independent of any GeoVac-side data. There is no GeoVac-side
constraint that determines Y because Y lives in a sector GeoVac doesn't
grade.

This reframes G2 (Sprint H1's "no autonomous Yukawa") as a corollary
of G3: chirality independence ⇒ no autonomous Y. The two gaps are
the same structural fact.

---

## 4. What the swaps tell us

The convention swaps in this diagnostic test whether some standard
NCG sign convention could absorb the residual:

- (a) sign-flip — equivalent to choosing a different convention for what
  "left" means. Doesn't help: the obstruction is the *independence* of
  the gradings, not their relative sign.
- (b) matter-only — physical-sector restriction. Doesn't help: matter-only
  still has both γ_F signs (left doublet vs right singlets) within itself.
- (c) Weyl-only — single-chirality restriction on the GeoVac side. Under
  σ_x γ_GV this is even more informative than under σ_z: γ_GV|_Weyl = 0,
  so Δ_c = −I_w ⊗ γ_F is just γ_F itself (negated and lifted). Weyl
  restriction strips γ_GV entirely and leaves γ_F's structure standing.
- (a+b), (a+b+c) — progressive restrictions. Drop more basis vectors but
  never collapse the residual on the kept subspace.

No standard convention swap can absorb this residual. The identification
of γ_GV with γ_F is structurally false.

---

## 5. Sector-resolved decomposition

(See `debug/data/g3c_tensor_chirality.json` for the full per-cell tables
under the σ_x convention. With σ_x γ_GV the per-basis-vector decomposition
is less direct than under σ_z because γ_GV mixes chi=+1 and chi=-1 blocks,
so individual basis vectors no longer have a well-defined γ_GV eigenvalue.
The Δ matrix has off-diagonal entries (from γ_GV ⊗ I_F mixing the chirality
blocks) and diagonal entries 0 − γ_F[i] = ±1. The aggregate spectrum
{−2 : 8·dim_w/2, 0 : 16·dim_w/2, +2 : 8·dim_w/2} arises from the joint
eigenstructure of (σ_x ⊗ I_w ⊗ I_F) and (I_2 ⊗ I_w ⊗ γ_F), which is
unambiguous and confirms the orthogonal-Z₂'s reading of §3.)

---

## 6. Compatibility with G3-A, G3-B, and the §VIII.C addendum

This re-run consumes the production γ_GV and γ_F directly:

- **γ_GV** from `geovac/chirality_grading.py` — G3-A's σ_x convention
  with γ²=I, {γ, D_truthful}=0, [γ, M]=0 all bit-exact at n_max ∈ {1,2,3}.
- **γ_F** from `geovac/almost_commutative.py` — G3-B's KO-dim 6 convention
  with all five Connes axioms at exact zero.

Both modules' test suites (50/50 in `tests/test_chirality_grading.py`,
53/53 in `tests/test_almost_commutative.py` after G3-B's TestChiralityF)
pass green. The §VIII.C addendum can cite the spectrum table in §2 above
verbatim.

The σ_x ↔ σ_z cross-check at n_max=2 (identical spectrum) shows the
NEGATIVE verdict is robust to either natural γ_GV convention. The
addendum should state the verdict in convention-independent language
("γ_GV and γ_F are independent commuting ℤ₂'s on H_GV ⊗ H_F") and cite
the σ_x convention's bit-exact axioms (G3-A) for the production-quality
construction.

---

## 7. Recommendation to PM

**Verdict: NEGATIVE on direct chirality identification.**

Three concrete next steps, ordered by tractability:

1. **Apply §VIII.C addendum extension.** The four-gap analysis in Paper 32
   §VIII.B already named G3 as open. This diagnostic closes G3 in the
   negative on S³ alone and sharpens it as a corollary of G2's "no
   autonomous Yukawa" finding from Sprint H1. The §VIII.C addendum should
   record:
   - γ_GV (σ_x ⊗ I_w on full-Dirac basis) and γ_F (KO-6 Connes-Marcolli)
     are structurally independent ℤ₂ gradings that commute on H_GV ⊗ H_F.
   - The canonical NCG combined grading γ₅ = γ_GV ⊗ γ_F exists and is
     well-defined.
   - "No autonomous Yukawa" (G2) and "no chirality identification on S³"
     (G3) are the same structural fact: Y lives in the off-diagonal block
     of D_F that flips γ_F, and γ_F is independent of γ_GV ⇒ no GeoVac-
     side constraint determines Y.
   - EW unification on S³ alone is structurally blocked by chirality
     independence; the cross-manifold direction (G4) is the remaining
     route.

2. **Pivot to G4 (cross-manifold).** The cleanest follow-up is to look at
   S³ × S⁵ or a master AC structure containing all three SM gauge groups
   across sub-manifolds. The G4 fork landed alongside this re-run; its
   verdict (FAR, with G4a Connes SM on T_{S³} reachable in 1–2 months)
   strengthens the priority of pursuing G4a as the next concrete sprint.

3. **Optional non-scalar multiplier extension.** If you want to explore
   chirality-mixing within S³ alone, consider replacing the scalar
   multiplier 𝒜_GV = C^∞(S³)|truncated with a Pauli-matrix-coupled or
   Clifford-tensored multiplier M_γ = γ^a M_a. This would promote 𝒜_GV
   from C^∞(S³) to C^∞(S³, M_2(ℂ)), break the scalar-multiplier-only
   assumption, and could in principle give a chirality-mixing ℤ₂ that
   identifies with γ_F. Speculative; would require a separate sprint.

---

## 8. Closing note

This is a clean, structural negative — the kind of result the §1.5
rhetoric rule recommends as the most informative kind. We did not fit
anything. We did not adjust conventions to make the answer come out
nicer. The two ℤ₂ gradings are simply orthogonal in the tensor product,
no choice of convention absorbs the residual, and the structural reason
(γ_GV is Dirac chirality, γ_F is weak isospin) is the natural physics-
layer explanation.

WH1 is PROVEN. Sprint H1 was POSITIVE-THIN. G3-C is NEGATIVE-CLEAN. That
sequence — proven structural alignment, partial admission with honest
gap, clean closure of the chirality direction in the negative — is the
kind of progress that makes the four-gap analysis a real roadmap and not
a wishlist. G4a is the next move.
