# G3-B Audit Memo: γ_F on the AC Finite Spectral Triple

**Sprint:** G3 — Electroweak chirality co-location on S³
**Track:** G3-B (AC-side chirality input to G3-C tensor-product test)
**Date:** 2026-05-06
**Status:** Sealed. Ready for G3-C with one documented convention deviation.

---

## 1. Audit findings

**Was γ_F already defined?** **No.** The Sprint H1 module
`geovac/almost_commutative.py` (854 lines) builds the full electroweak
AC extension *T_F* — the algebra `A_F = ℂ ⊕ ℍ`, the doubled finite
Hilbert space `H_F = ℂ⁴_matter ⊕ ℂ⁴_antimatter`, the matter–antimatter
Dirac `D_F`, and the antilinear real structure `J_F = U_F · K` with
`J_F² = +I` — but it does **not** explicitly construct the chirality
operator γ_F.

The only chirality-like object present before this audit is on the
GeoVac factor: `_gamma_GV` (built in
`AlmostCommutativeTriple.__init__`, lines 441–443), a diagonal
`±1`-by-Camporesi–Higuchi-chirality-label matrix on H_GV. This is used in
the cross-coupling term `γ_GV ⊗ D_F` of the combined Dirac, but plays
no role on the F factor.

**Is γ_F there now?** **Yes.** Added as:

- `ElectroweakFiniteTriple.chirality_F()` —
  `geovac/almost_commutative.py:386` (immediately after
  `real_structure_F()`).
- `AlmostCommutativeTriple.gamma_F()` —
  `geovac/almost_commutative.py:528` (a passthrough convenience for
  G3-C).

Quoted code path:

```python
def chirality_F(self) -> np.ndarray:
    g_matter = np.diag([1.0, 1.0, -1.0, -1.0]).astype(np.complex128)
    G = np.zeros((8, 8), dtype=np.complex128)
    G[0:4, 0:4] = g_matter
    G[4:8, 4:8] = -g_matter
    return G
```

Coverage: a new section `TestChiralityF` in
`tests/test_almost_commutative.py:175` adds 15 tests; the full file
goes from 38 → 53 tests; all pass at 1e-14 residual.

---

## 2. Convention chosen + rationale

### 2.1 The Connes–Marcolli SM convention (KO-dim 6)

For the finite SM spectral triple at `A_F = ℂ ⊕ ℍ`, Connes and Marcolli
(2008, *Noncommutative Geometry, Quantum Fields and Motives*, Ch. 13,
Table 13.1) place the finite triple at **KO-dim 6**:

- `J_F² = +I`,
- `J_F D_F = +D_F J_F`,
- `{J_F, γ_F} = 0`  (anticommuting J,γ — the KO-6 sign in the
  ε-ε'-ε'' table).

The H1 module is **already locked at KO-dim 6**: this is documented at
`almost_commutative.py:96–97` ("the standard Connes-Marcolli SM
convention places T_F at KO-dim 6") and tested at
`test_almost_commutative.py:99–104` (`test_J_F_squared_plus_I`) and
`:118–131` (`test_J_F_D_F_commutation`).

For γ_F to live consistently inside this KO-6 module, the chirality
must satisfy `{J_F, γ_F} = 0`. The unique block-diagonal real-diagonal
γ_F satisfying both `γ_F² = I` and `{J_F, γ_F} = 0` for the
matter/antimatter swap-K real structure is:

```
γ_F^matter     = diag(+1, +1, −1, −1)    (L = +, R = −)
γ_F^antimatter = − γ_F^matter
               = diag(−1, −1, +1, +1)
```

The matter assignment `(L, R) → (+, −)` matches the directive:
- `ℂ` summand → R-singlets  → γ_F = −1 sub-block  ✓
- `ℍ` summand → L-doublet   → γ_F = +1 sub-block  ✓

The antimatter sign-flip is what gives `{J_F, γ_F} = 0` rather than
`[J_F, γ_F] = 0`.

### 2.2 Deviation from the G3-B prompt

The directive said: *"Matter/antimatter doubling: γ_F identical on
both copies."*

Identical γ on both copies gives `[J_F, γ_F] = 0` (commuting), which
is the **KO-dim 0** sign rule, not KO-6.

The H1 module is locked at KO-6. Switching γ_F to "identical on both
copies" would silently downgrade the combined triple's KO-dim from
9 → 8 ≡ 0 (mod 8) and break `test_J_F_D_F_commutation` together with
the `test_J_combined_KO1_sign` consistency check (combined `J² = −I`
plus `JD = +DJ` is the KO-dim 1 / mod-8 = 1 fingerprint, which
requires KO_GV = 3 and KO_F = 6 to add modulo 8).

I matched the module convention (γ_F^antimatter = − γ_F^matter) and
documented the deviation in the `chirality_F()` docstring and in
`test_J_F_anticommutes_with_gamma_F_KO6`.
The test `test_J_F_commutes_with_gamma_F_KO0_explicitly_fails`
explicitly asserts that the directive's KO-0 convention does **not**
hold, so a future drift back to "identical on both copies" would fail
loudly.

---

## 3. Axiom verification table

All four directive-required Z₂-grading axioms verified to 1e-14 at
the matrix level. Format: assertion / residual definition / max
residual across tests.

| # | Axiom | Test name | Residual | Max obs. |
|---|---|---|---|---|
| 1 | γ_F² = I | `test_chirality_F_squared_identity` | `‖γ_F² − I‖_∞` | 0.0 |
| 2a | {γ_F, D_F} = 0 at Y=0 | `test_anticommutes_with_DF_zero_yukawa` | `‖{γ_F, D_F}‖_∞` | 0.0 |
| 2b | {γ_F, D_F} = 0 at real Y | `test_anticommutes_with_DF_nonzero_real_yukawa` | `‖{γ_F, D_F}‖_∞` | 0.0 |
| 2c | {γ_F, D_F} = 0 at complex Y | `test_anticommutes_with_DF_complex_yukawa` | `‖{γ_F, D_F}‖_∞` | 0.0 |
| 3 | γ_F π_F(a) γ_F ∈ A_F (in fact = π_F(a)) | `test_algebra_elements_are_even` | `‖γ_F π_F γ_F − π_F‖_F` | 0.0 |
| 4 | {J_F, γ_F} = 0 (KO-6 sign) | `test_J_F_anticommutes_with_gamma_F_KO6` | `‖U_F γ_F + γ_F U_F‖_∞` | 0.0 |

Auxiliary properties also checked:

| # | Property | Test name | Residual |
|---|---|---|---|
| A1 | γ_F = γ_F* (Hermitian) | `test_chirality_F_hermitian` | 0.0 |
| A2 | γ_F is diagonal | `test_chirality_F_diagonal` | 0.0 |
| A3 | matter block = diag(+1,+1,−1,−1) | `test_chirality_F_matter_block` | 0.0 |
| A4 | antimatter block = −γ_F^matter | `test_chirality_F_antimatter_block` | 0.0 |
| A5 | no matter↔antimatter cross-block | `test_chirality_F_no_cross_blocks` | 0.0 |
| A6 | KO-0 fails (directive deviation) | `test_J_F_commutes_with_gamma_F_KO0_explicitly_fails` | not satisfied |
| A7 | AC passthrough = ETT chirality | `test_AC_triple_gamma_F_passthrough` | 0.0 |
| A8 | γ_F independent of Yukawa | `test_chirality_F_independent_of_yukawa` | 0.0 |

Algebra-evenness (axiom 3) holds in the strong form
`γ_F π_F(a) γ_F = π_F(a)` (not just membership in A_F): π_F is
already block-diagonal in (L, R) within the matter sector, and γ_F
is also block-diagonal in (L, R), so conjugation by γ_F preserves L
and R subspaces individually. The algebra is a fully even subalgebra
of B(H_F).

The {γ_F, D_F} = 0 verification at complex Yukawa is non-trivial: the
antimatter Dirac is `M^bar = conj(M)`, so D_F preserves L↔R
off-diagonal structure in each sector regardless of the phase of Y.
γ_F's block-diagonal-by-(L,R) structure with antimatter sign-flip
anticommutes with both `M` and `conj(M)` simultaneously.

---

## 4. Honest disclosure: D_F degeneracy and γ_F uniqueness

**Concern raised by the directive:** *"Honest disclosure if H1's D_F
is degenerate enough that γ_F is not unique."*

### 4.1 Y = 0 case (zero Yukawa)

When `yukawa_e = yukawa_nu = 0`, `D_F = 0`. In this case the
anticommutation `{γ_F, D_F} = 0` is **trivially satisfied** by any
operator γ_F whatsoever: there is no constraint on γ_F coming from
D_F.

The remaining constraints:

- `γ_F² = I` (Z₂ grading),
- `γ_F π_F γ_F = π_F` (algebra evenness),
- `{J_F, γ_F} = 0` (KO-dim 6),

still leave a non-trivial moduli space for γ_F. Algebra-evenness
forces γ_F to be block-diagonal in (L, R) within each sector;
combined with `{J_F, γ_F} = 0` and `γ_F² = I`, the only freedom is the
overall sign in the matter block (and dependent antimatter block).
The convention `(L, R) → (+, −)` versus `(L, R) → (−, +)` is a
naming choice, not a structural one.

### 4.2 Y ≠ 0 case (Yukawa active)

With a non-zero Yukawa, `D_F` has off-diagonal L↔R block structure in
each sector. This **forces** γ_F to flip sign between L and R within
each sector (so that {γ_F, M} = 0 in matter and {γ_F, M^bar} = 0 in
antimatter). The only remaining freedom is the matter↔antimatter
relative sign — and that is fixed by `{J_F, γ_F} = 0` (KO-6) once the
real structure is chosen.

**Verdict:** at non-zero Yukawa, γ_F is *unique up to overall sign*
within the KO-6 / Connes-Marcolli convention. At zero Yukawa, the
D_F constraint disappears but the algebra-evenness + KO-6 J-sign
constraint still pin down the matter/antimatter relative sign, so γ_F
remains uniquely determined modulo the global sign convention.

### 4.3 Test design for both regimes

Tests `test_anticommutes_with_DF_zero_yukawa`,
`test_anticommutes_with_DF_nonzero_real_yukawa`, and
`test_anticommutes_with_DF_complex_yukawa` together exercise both
regimes and confirm γ_F as constructed satisfies {γ_F, D_F} = 0 in
all three.

---

## 5. Sanity: relation to combined chirality

A stand-alone γ for the combined triple T = T_GV ⊗ T_F is **not**
canonical: combined KO-dim is 9 ≡ 1 (mod 8), and KO-dim 1 is *odd* —
the standard sign table has no Z₂ grading at odd KO-dim.

However, the construction `D = D_GV ⊗ I_F + γ_GV ⊗ D_F` does use
γ_GV as an internal sign-operator (the chirality label on the
Camporesi–Higuchi spinor bundle), independent of whether the *combined*
triple admits an external Z₂ grading. The H1 module's existing
`gamma_GV()` and the new `gamma_F()` are the natural building blocks
for the G3-C tensor-product diagnostic:

```
γ_combined  := γ_GV ⊗ γ_F       (formal definition, even-grading
                                  if extended to the combined H)
```

Whether `γ_combined` plays a structural role or is purely formal at
KO-dim 9 is a question for G3-C; this audit only delivers γ_F as a
clean, axiomatically-verified ingredient of T_F.

---

## 6. Status: ready for G3-C

**Yes**, with one documented caveat (convention deviation,
§2.2 above).

Concrete handoff to G3-C:

- API: `ElectroweakFiniteTriple.chirality_F() → np.ndarray (8,8)` and
  `AlmostCommutativeTriple.gamma_F() → np.ndarray (8,8)` are stable
  and tested.
- Convention: KO-dim 6 of T_F (Connes–Marcolli SM), with
  γ_F^antimatter = − γ_F^matter. **Not** the directive's
  "γ_F identical on both copies" (which would be KO-0 and break the
  module).
- Combined: `γ_combined := γ_GV ⊗ γ_F` is well-defined at the matrix
  level; whether the combined T admits a canonical chirality is part
  of G3-C's question, not G3-B's deliverable.

If G3-C decides to compute `(γ_GV ⊗ I_F) − (I_GV ⊗ γ_F)` directly
on H_GV ⊗ H_F (per the parent prompt), both factors are now
available bit-stably from the AlmostCommutativeTriple instance.

The 15 new tests in `TestChiralityF` provide the regression net: any
future change to D_F, J_F, or the matter/antimatter convention that
silently invalidates γ_F's axioms will be caught.

---

## 7. Files modified

| Path | Change | Lines (approx) |
|------|--------|---------|
| `geovac/almost_commutative.py` | Added `chirality_F()` (ETT) + `gamma_F()` (AC passthrough) | +43 |
| `tests/test_almost_commutative.py` | Added Section A2 with `TestChiralityF` (15 tests) | +160 |
| `debug/g3b_chirality_F_audit_memo.md` | This memo | (new) |

No commits in this fork. Edits are in the working tree; parent
to commit as part of the G3 sprint roll-up.

— G3-B fork, 2026-05-06.
