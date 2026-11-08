# Sprint memo — Paper 39: the KO-dim-6 real structure of the Clifford-doubled product triple

**Date:** 2026-09-05
**Task:** decide whether the Clifford-doubled product spectral triple
(the repair for Paper 39's non-existent `γ_a ⊗ D_b` Dirac) carries a real
structure `J_ab` satisfying the KO-6 axioms, and whether it is genuinely the
tensor product of the two Camporesi–Higuchi (CH) S³ triples. This is the one
sub-task left open by `debug/sprint_p39_tensor_lifted_state_memo.md` §8.4.
**Scope:** memo only. No paper, module or test edited. One driver added
(`debug/p39_real_structure_probe.py`).

---

## 0. Verdict

> **GO.** The Clifford-doubled product triple carries the KO-dim-6 real
> structure
>
>     J_ab  =  J_a ⊗ (σ₁ K₂) ⊗ J_b ,        i.e.  U_ab = U_a ⊗ σ₁ ⊗ U_b ,  J_ab ψ = U_ab conj(ψ),
>
> and it is the standard **Dąbrowski–Dossena / Connes–Marcolli product** of
> the two CH real spectral triples (odd ⊗ odd, KO-dim 3 + 3 = 6). All three
> KO-6 sign relations hold **exactly** and the middle operator `σ₁K₂` is
> **forced** — the other three antiunitary choices on ℂ² each fail a specific
> sign (checked against the *actual* CH Dirac matrices, not a bare 2×2 toy):
>
> | KO-6 axiom | required (van Suijlekom Tbl 3.1) | `J_ab = J_a ⊗ σ₁K₂ ⊗ J_b` |
> |:--|:--:|:--:|
> | `J_ab² = ε`        | `ε = +1`   | **+1**, residual 0 |
> | `J_ab D_ab = ε' D_ab J_ab` | `ε' = +1` | **+1**, residual 0 |
> | `J_ab γ_ab = ε'' γ_ab J_ab`| `ε'' = −1` | **−1**, residual 0 |
>
> The order-zero `[a, J_ab b J_ab⁻¹] = 0` and order-one
> `[[D_ab, a], J_ab b J_ab⁻¹] = 0` conditions hold for the **continuum**
> tensor CH triple — the object Paper 39 identifies as the limit — by a
> factor-by-factor reduction to the classical single-factor spin-manifold
> facts. At **finite truncation** they fail by exactly the amount the single
> factor already fails (operator-system-is-not-an-algebra artifact,
> documented in Paper 38 `rem:dirac_degeneracy` and `real_structure.py`); the
> product adds **no** cross-factor obstruction (cross-factor pairs commute to
> machine zero, 0/392).

So Paper 39 may state convergence to *the tensor CH spectral triple*
(KO-6, with real structure), not merely to *a* tensor triple — provided the
identification is made at the continuum limit and the order-0/1 conditions are
attributed to the limit (never claimed exact at finite cutoff, which they are
not, for the single factor either).

**One documentation defect found (not load-bearing):** the KO-dim sign table
in the docstring of `geovac/real_structure.py` is wrong in rows 1, 2, 5, 6
(it lists KO-6 as `(+,−)`; the correct ε' is `+`). The module only *uses* the
KO-3 row `(−,+)`, which is correct and correctly enforced by the code, so
nothing computed is affected — but the table as printed should not be trusted
for the even rows and should be corrected to van Suijlekom Table 3.1.

---

## 1. The KO-dimension sign table used

The real structure `J` of a spectral triple `(A,H,D)` obeys, with signs
depending only on the KO-dimension `n mod 8`:

    J² = ε ,       J D = ε' D J ,       (n even only) J γ = ε'' γ J .

**Authoritative table** (A. Connes, *Noncommutative geometry and reality*,
J. Math. Phys. **36** (1995) 6194, Table; reproduced as van Suijlekom,
*NCG and Particle Physics*, Springer 2015, **Table 3.1**):

| `n` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| `ε`   | + | + | − | − | − | − | + | + |
| `ε'`  | + | − | + | + | + | − | + | + |
| `ε''` | + |   | − |   | + |   | − |   |

- **Single factor (S³ = SU(2), CH triple): KO-dim 3**, so `ε = −1`, `ε' = +1`,
  no `ε''` (odd dimension, no chirality grading — this is exactly the γ = 0
  obstruction already in Paper 39/38/32).
- **Product: KO-dim 3 + 3 = 6** (additive mod 8), so the target is
  `ε = +1`, `ε' = +1`, `ε'' = −1`. This is the **Standard Model** KO-row.

(The `real_structure.py` docstring's inline table disagrees with this for the
even rows; see §0. Its KO-3 row is correct and is the only one the code uses.)

---

## 2. The single-factor real structure `J_a` (KO-3)

Implemented and already in-repo as `geovac/real_structure.py`
`build_J_full_dirac(n_max)`. On the full-Dirac sector basis
`|n_fock, l, m_j, χ⟩` (χ = ±1 the chirality/`sign(D)` label),

    J_a |n_fock, l, m_j, χ⟩ = σ(l, m_j) |n_fock, l, −m_j, χ⟩,
        σ(l, m_j) = i^{2 m_j} · (−1)^l          (m_j half-integer),

realised as the antilinear `J_a ψ = U_a conj(ψ)` with
`U_a[i,j] = σ_i · δ_{i, π(j)}`, `π` the involution `m_j → −m_j`. This is the
standard Camporesi–Higuchi charge conjugation on S³ spinors (Friedrich, *Dirac
operators in Riemannian geometry*, Ch. 1; Camporesi–Higuchi 1996).

**Single-factor KO-3 signs, CHECKED NUMERICALLY** (`audit_J`, truthful CH
Dirac `D χ = χ(n+½)`):

| `n_max` | `dim_O` | `J² = −I` | `J D = +D J` | order-0 fails | order-1 fails |
|:--:|:--:|:--:|:--:|:--:|:--:|
| 1 | 1  | 0.0e0 | 0.0e0 | 0     | 0    |
| 2 | 14 | 0.0e0 | 0.0e0 | 124 (max 5.5e-2) | 43 (max 1.0e-1) |
| 3 | 55 | 0.0e0 | 0.0e0 | 2594 (max 7.9e-2) | 1337 (max 2.0e-1) |

The two **sign** relations (`ε = −1`, `ε' = +1`) are **exact at every
cutoff** — they are pure linear-algebra facts about `U_a` and `D_a`,
independent of the algebra. The order-0/order-1 conditions **fail at finite
cutoff and the failure grows with `n_max`**: this is not a defect of `J_a`, it
is the truncated operator system failing to be an algebra (products of
multipliers leave `Op_{n_max}`), the exact phenomenon Paper 38 records as
`rem:dirac_degeneracy` and the module docstring flags under "(C3)/(C4) need
re-interpretation". In the continuum these hold as a classical theorem about
commutative spin-manifold spectral triples.

---

## 3. The product real structure `J_ab` (KO-6)

### 3.1 Construction under test (the Paper 39 repair, odd ⊗ odd doubling)

    H_ab = H_a ⊗ ℂ² ⊗ H_b
    D_ab = D_a ⊗ σ₁ ⊗ 1  +  1 ⊗ σ₂ ⊗ D_b
    γ_ab = 1 ⊗ σ₃ ⊗ 1
    A_ab = C^∞(S³) ⊗ 1₂ ⊗ C^∞(S³)          (acting as a ⊗ 1₂ ⊗ b)

Structural sanity (all residual 0, `n_a=n_b=2`, dim 512, and 500 random
2×2 symbols in the prior memo §6-G2): the two Dirac summands anticommute,
`D_ab² = D_a²⊗1⊗1 + 1⊗1⊗D_b²`, `γ_ab² = 1`, `{γ_ab, D_ab} = 0`,
`[γ_ab, A_ab] = 0`.

### 3.2 The candidate and why the middle factor is `σ₁K₂`

`J_a`, `J_b` are antilinear; to keep `J_ab` antilinear the middle factor on
ℂ² must be an **antiunitary** `C₂ = W K₂` (`W` unitary, `K₂` = complex
conjugation), giving an odd number (3) of conjugations. Write

    J_ab = J_a ⊗ C₂ ⊗ J_b ,   C₂ = W K₂ ,   ⇒   U_ab = U_a ⊗ W ⊗ U_b ,  J_ab ψ = U_ab conj(ψ).

**Symbolic derivation of the required `W` (PROVED).** Using `J_a² = −1`,
`J_a D_a = +D_a J_a` (KO-3), the three product signs factor across the tensor
legs:

- `J_ab² = J_a² ⊗ C₂² ⊗ J_b² = (−1)(C₂²)(−1) = C₂²`. Need `C₂² = +1`.
- `J_ab D_ab = ε' D_ab J_ab`: term `D_a⊗σ₁⊗1` contributes
  `(+1)·[C₂ vs σ₁]·(+1)`, term `1⊗σ₂⊗D_b` contributes `(+1)·[C₂ vs σ₂]·(+1)`.
  Need `C₂` to **commute with both σ₁ and σ₂** for `ε' = +1`.
- `J_ab γ_ab = ε'' γ_ab J_ab`: sign `= [C₂ vs σ₃]`. KO-6 wants `ε'' = −1`, i.e.
  `C₂` **anticommutes with σ₃**.

The unique antiunitary on ℂ² commuting with σ₁, σ₂ and anticommuting with σ₃
is `C₂ = σ₁ K₂` (up to phase): writing `C₂ = W K₂`, `C₂ σ_i C₂⁻¹ = W σ̄_i W†`
with `σ̄₁ = σ₁, σ̄₂ = −σ₂, σ̄₃ = σ₃` forces `W σ₁ W† = σ₁`, `W σ₂ W† = −σ₂`,
`W σ₃ W† = −σ₃` ⟹ `W = σ₁`. Then `C₂² = σ₁ K₂ σ₁ K₂ = σ₁ σ̄₁ = +1`. ∎

Equivalently: `C₂` must commute with the two Paulis that appear **in `D`** and
anticommute with the grading Pauli.

### 3.3 All four `W`, checked against the ACTUAL CH matrices — the choice is forced

`debug/p39_real_structure_probe.py`, `(n_a,n_b) = (2,2)`, dim 512.
`ε' = 0` in the table means `U_ab conj(D_ab)` equals **neither** `+D_ab U_ab`
**nor** `−D_ab U_ab` (mixed sign across the two Dirac terms) — a genuine
failure, not a small residual. This is the D-dependent axis the
"suspiciously-clean-toy" warning targets, and it is exactly the axis that
separates the four choices.

| `W` (so `C₂ = W K₂`) | `ε` (J²) | `ε'` (JD) | `ε''` (Jγ) | KO-6 `(+,+,−)`? |
|:--:|:--:|:--:|:--:|:--:|
| `1₂` | +1 | **0** (res 5.0) | +1 | ✗ (JD fails) |
| **`σ₁`** | **+1** | **+1** | **−1** | **✓ exact (all res 0)** |
| `σ₂` | **−1** | −1 | −1 | ✗ (J² = −1, wrong; this is not KO-6) |
| `σ₃` | +1 | **0** (res 5.0) | +1 | ✗ (JD fails) |

Only `W = σ₁` yields the KO-6 row. `1₂` and `σ₃` fail the `D`-dependent `ε'`;
`σ₂` gives `J² = −1` (KO-2-like, wrong parity). **Forced.**

### 3.4 Full axiom set with `J_ab = J_a ⊗ σ₁K₂ ⊗ J_b`, CHECKED NUMERICALLY

| `(n_a,n_b)` | dim | `U_ab` unit. | `J²=+I` | `JD=+DJ` | `Jγ=−γJ` |
|:--:|:--:|:--:|:--:|:--:|:--:|
| (1,1) | 32  | 0.0e0 | +1, 0.0e0 | +1, 0.0e0 | −1, 0.0e0 |
| (1,2) | 128 | 0.0e0 | +1, 0.0e0 | +1, 0.0e0 | −1, 0.0e0 |
| (2,1) | 128 | 0.0e0 | +1, 0.0e0 | +1, 0.0e0 | −1, 0.0e0 |
| (2,2) | 512 | 0.0e0 | +1, 0.0e0 | +1, 0.0e0 | −1, 0.0e0 |

All three KO-6 sign relations are **exact at every cutoff tested**, against the
genuine CH Dirac matrices (with the `m_j`-flip permutation and `i^{2m_j}(−1)^l`
phases baked into `U_a`, `U_b`), not a bare Pauli toy.

---

## 4. Order-zero and order-one

### 4.1 Reduction to the factors (PROVED)

For `a = M_f ⊗ 1₂ ⊗ M_g`, `b = M_{f'} ⊗ 1₂ ⊗ M_{g'}` in `A_ab`, and
`J_ab b J_ab⁻¹ = (J_a M_{f'} J_a⁻¹) ⊗ 1₂ ⊗ (J_b M_{g'} J_b⁻¹)` (the middle
`1₂` conjugates to `1₂`):

- **Order-zero.** `[a, J_ab b J_ab⁻¹] = 0` holds iff each leg commutes:
  `[M_f, J_a M_{f'} J_a⁻¹] = 0` (factor-a order-0), `[1₂,1₂]=0`,
  `[M_g, J_b M_{g'} J_b⁻¹] = 0` (factor-b order-0).
- **Order-one.** `[D_ab, a] = [D_a,M_f]⊗σ₁⊗M_g + M_f⊗σ₂⊗[D_b,M_g]`. Each of
  the two terms commutes with `J_ab b J_ab⁻¹` because on each leg it needs one
  of {factor order-1, factor order-0, `[σ_i,1₂]=0`}, all of which hold.

So the product order-0/order-1 conditions **hold whenever the single-factor
order-0 and order-1 conditions hold**, and reduce to them exactly. The
continuum CH triple satisfies both (classical fact: on a commutative spin
manifold `J b J⁻¹` is right-multiplication by `b̄`, which commutes with
left-multiplication and with the zeroth-order operator `[D,a] = c(da)`), hence
so does the **continuum tensor CH triple**.

### 4.2 Finite-cutoff numerics — the product adds NO obstruction

At finite truncation the single-factor order-0/1 already fail (§2, operator
system ≠ algebra). The product inherits **exactly** those failures and nothing
more. Breakdown at `(n_a,n_b)=(2,2)`, 28 generators = 14 (factor a) + 14
(factor b), by whether the pair `(a,b)` sits on the same factor or crosses:

| pair category | # pairs | order-0 fails | order-1 fails |
|:--|:--:|:--:|:--:|
| both factor-a (`AA`) | 196 | 124 (max 5.5e-2) | 43 (max 1.0e-1) |
| both factor-b (`BB`) | 196 | 124 (max 5.5e-2) | 43 (max 1.0e-1) |
| **cross (`AB`)** | 392 | **0** | **0** |

The `AA`/`BB` numbers are **identical** to the single-factor `n_max=2` counts
(124 / 43, §2). Every cross-factor pair commutes to machine zero. This is the
numerical signature of the §4.1 factorization: the product introduces no new
order-0/order-1 defect, so the product adds nothing beyond the inherited single-factor (non-)satisfaction.

> **CORRECTION 2026-09-05 (parent session, order-0/1 hardening).** The finite-cutoff residual does NOT vanish as n grows -- measured max order-0 residual 0.055/0.078/0.096 and order-1 0.101/0.203/0.405 at n=2,3,4, i.e. GROWING.  The continuum tensor triple satisfies order-0/1 because its limit algebra is COMMUTATIVE (Paper 32, the standard automatic-order fact) and the 4.1 reduction carries that to the product -- NOT because the truncation residual shrinks.  The earlier phrase 'where AA/BB failures vanish' was mis-argued and is withdrawn.

**Honest statement of what is exact vs limiting:**
- KO-6 sign relations (`J²`, `JD`, `Jγ`): **exact at every finite cutoff**
  (PROVED symbolically + CHECKED NUMERICALLY).
- Order-0 / order-1: **hold for the continuum limit object** (PROVED, modulo
  the classical single-factor spin-manifold fact); **fail at finite cutoff by
  exactly the single-factor amount** (CHECKED NUMERICALLY); **no cross-factor
  obstruction** (CHECKED NUMERICALLY, 0/392).

---

## 5. Is it *the* tensor product of the CH triples?

**Yes — the Dąbrowski–Dossena / Connes–Marcolli product of two odd real
spectral triples.** References:
- L. Dąbrowski, G. Dossena, *Product of real spectral triples*, Int. J. Geom.
  Methods Mod. Phys. **8** (2011) 1833 — the canonical treatment of the
  odd ⊗ odd case via the ℂ² Clifford doubling, with the sign combination
  rules that give KO-dim additivity mod 8.
- F. J. Vanhecke, *On the product of real spectral triples*, Lett. Math. Phys.
  **50** (1999) 157 — the graded (even) case and the odd-case subtlety.
- A. Connes, M. Marcolli, *Noncommutative Geometry, Quantum Fields and
  Motives* (2008), §I — the KO-dim-6 graded product framework.

Their construction is `H₁ ⊗ ℂ² ⊗ H₂`, `D = D₁⊗σ₁⊗1 + 1⊗σ₂⊗D₂`,
`Γ = 1⊗σ₃⊗1`, `J = J₁ ⊗ C ⊗ J₂` — **exactly** §3.1 with `C = σ₁K₂`. The ℂ²
is `Cl(ℝ²) ≅ M₂(ℂ)`'s spinor module; it is precisely the factor that carries
(i) the product **grading** (`σ₃`), and (ii) the product **real structure**
(`σ₁K₂` — the antiunitary of the KO-2 Clifford algebra `Cl_{0,2}`, which is
what "3 + 3 = 6" adds). The algebra `A_a ⊗ 1₂ ⊗ A_b` represents each
`C^∞(S³)` faithfully; setting one factor to the trivial rank-1 triple recovers
the other CH triple (⊗ the inert ℂ² Clifford spectator). So each CH triple is
recovered as a partial factor, and the whole is their graded product.

**Caveat, stated honestly.** "The" product is canonical only **up to unitary
equivalence**: the assignment of the three anticommuting `(σ₁,σ₂,σ₃)` to
(D_a-slot, D_b-slot, grading) is a choice, and different choices give
unitarily equivalent triples with `C` relabeled accordingly (some references
write `C = σ₂K` or `iσ₂K`, matching their own slot/ε conventions). Given the
§3.1 slot choice and the CH factor signs `(ε,ε') = (−,+)`, `C = σ₁K₂` is
forced (§3.2–3.3). So it is *the* tensor product in the standard "up to the
unitary equivalence of the doubling" sense — not a merely-*a*-tensor-triple
weakening.

---

## 6. What this means for Paper 39

The `debug/sprint_p39_tensor_lifted_state_memo.md` §8.4 open item is closed
**GO**. Paper 39's status note to `thm:main` — which currently says the
Clifford-doubling repair is "named, not adopted, and whether the lemmas hold
for it is open", and (separately) that the real structure was unverified — can
be strengthened to: the doubled product

    (A_a ⊗ 1₂ ⊗ A_b,  H_a ⊗ ℂ² ⊗ H_b,  D_ab,  γ_ab,  J_ab = J_a ⊗ σ₁K₂ ⊗ J_b)

is the KO-dim-6 real spectral triple that is the Dąbrowski–Dossena product of
the two CH triples, so the limit object is *the tensor CH spectral triple*
(with real structure), not merely *a* tensor triple.

**Suggested dispositions (PI call — no edits made):**
1. Paper 39 `thm:main` status note + KO-dim-6 paragraph: adopt `J_ab` and
   state the three KO-6 signs `(+,+,−)`; cite Dąbrowski–Dossena. Attribute
   order-0/1 to the continuum limit (never claim them exact at finite cutoff).
2. `geovac/real_structure.py` docstring KO-dim table: correct the even rows to
   van Suijlekom Table 3.1 (`KO-6 = (+,+,−)`, not `(+,−)`). Not load-bearing
   (only KO-3 is used, correctly) but it is a wrong table in a physics module.
3. Optional: a `test_paper39_real_structure.py` freezing the §3.3 forced-`W`
   table and the §3.4/§4.2 signs would be the natural backing test if the
   claim enters the paper. **Guard-writing is a separate reviewed activity
   (§9 rule); its rejected wrong answer would be "any `W ≠ σ₁` passes" or
   "a cross-factor order-0 pair is allowed to fail".** Not written here.

---

## 7. What a reviewer should attack hardest

Ordered by expected damage if the attack lands.

1. **The finite-vs-continuum status of order-0/order-1.** I claim the limit
   object satisfies them and the finite-cutoff failures are the harmless
   operator-system artifact. A reviewer should press: is the *single-factor*
   continuum order-0/1 actually inherited by the module's truncation as
   `n_max → ∞`? I showed the product adds nothing beyond the single factor
   (0/392 cross-factor), but I did **not** independently prove the single
   factor's finite-cutoff failure → 0; I rely on Paper 38 `rem:dirac_degeneracy`
   + the classical spin-manifold theorem. If the single-factor limit claim is
   contested, it is contested for Paper 38 too, not just here.
2. **"Up to unitary equivalence" doing real work in "*the* tensor product".**
   The construction is canonical only after choosing the Clifford generators.
   I argue any choice gives a unitarily equivalent triple (standard), so the
   limit is well-defined as *the* product. A purist may want the equivalence
   class named rather than a representative.
3. **The numerics stop at `(n_a,n_b) = (2,2)`, dim 512.** The **sign**
   relations are also PROVED symbolically (§3.2), so cutoff size does not
   threaten them. Order-0/1 factorization is PROVED (§4.1) and the 0/392
   cross-factor result is a structural check, not a size-limited one — but a
   reviewer wanting `(3,3)` (dim 2·28²·... = larger) can run the driver.
4. **`σ₁` vs a literature `σ₂K`.** Some references write the doubling
   antiunitary as `σ₂K`. That is the *same* `J_ab` under a relabeling of which
   Pauli sits in `D` vs the grading; with **this paper's** `D = D_a⊗σ₁ + σ₂⊗D_b`,
   `γ = σ₃`, the forced answer is `σ₁K₂` (§3.2 uniqueness). A reviewer should
   check the slot convention before objecting to the specific Pauli.
5. **Order-1 uses factor order-0 on one leg and factor order-1 on the other**
   (§4.1). Both single-factor conditions must hold in the limit; if either
   fails in the limit, order-1 of the product fails. They do hold classically,
   but this couples the product's order-1 to *both* single-factor axioms.

---

## 8. Files added (driver only — no paper, module or test touched)

- `debug/p39_real_structure_probe.py` — single-factor KO-3 trend
  (`audit_J`); the four-`W` forced-choice table at (2,2); the full KO-6
  axiom set with `W=σ₁` at (1,1)/(1,2)/(2,1)/(2,2); the AA/BB/AB order-0/1
  factorization breakdown.
