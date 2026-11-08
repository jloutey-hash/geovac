# Door 4c — Does the combined real structure J = J_GV ⊗ J_F SELECT ℍ, or merely ADMIT it?

*Closes (negatively) the ℍ-vs-M₂(ℂ) fork left open by Door 4b.*
*Structure-only / sign-table probe. NO Yukawa value selected (guardrail).*
Started + closed 2026-06-01. Driver: `debug/door4c_j_signtable_audit.py`.
Data: `debug/data/door4c_j_signtable_audit.json`.

---

## TL;DR verdict — **ADMITS ℍ (does NOT force ℍ)**

The conjectured GeoVac-native handle — that the outer triple's
`J_GV² = −1` at KO-dimension 3 (the pseudoreal Spin(3)=SU(2) two-spinor
sign) selects ℍ over M₂(ℂ) at the n=2 inner factor — **does NOT close the
fork through the combined-J / KO-dim sign table.** The fork stays open at
the same place Door 4b left it. ℍ is a **literature import** (the
Chamseddine–Connes–Marcolli complex-fermion-rep / second-order-condition
argument on the finite factor), **not a GeoVac forcing** derivable from
`J_GV`.

This is the honest, likely outcome the decision gate anticipated. The
**inner algebra status remains PARTIAL** (Door 4b unchanged): factor count,
n=1 = ℂ, n=3 = M₃(ℂ) forced; n=2 = ℍ-vs-M₂(ℂ) admitted-not-forced; N_gen
and inner KO-dim free.

---

## The one structural fact that settles it

Over ℂ, **ℍ and M₂(ℂ) are the same complex algebra**: ℍ ⊗_ℝ ℂ ≅ M₂(ℂ).
"ℍ vs M₂(ℂ)" is therefore **not** a choice of complex algebra — it is a
choice of **real form**, fixed by an antilinear involution `j` on M₂(ℂ):

| Real form | Involution `j(m)` | Structure op `J_int` | `J_int²` | Fix(j) | Division ring? |
|:----------|:------------------|:---------------------|:--------:|:-------|:--------------:|
| **ℍ** | `(iσ₂) m̄ (iσ₂)⁻¹` | `J_q = (iσ₂)·K` | **−1** | span_ℝ{1, iσ₁, iσ₂, iσ₃} | **YES** (every nonzero element invertible; min \|det\|/‖·‖² = ½) |
| **M₂(ℂ)/M₂(ℝ)** | `m̄` | `J_r = K` | **+1** | M₂(ℝ) | **NO** (contains diag(1,0), singular) |

(Part 1 of the driver, all bit-exact; division-ring witness `diag(1,0)`
projects to a singular nonzero element of Fix(K) but never of Fix(j_q).)

So the n=2 division ring is decided by **one sign** — the square of the
*internal* C² antilinear involution: `J_int² = −1 ⇔ ℍ`, `J_int² = +1 ⇔ M₂(ℂ)`.
This internal involution lives on the **algebra side** (it constrains which
2×2 matrices are admissible). It is **a different object** from the
finite-triple real structure `J_F` (the matter↔antimatter swap-conjugation),
and a different object again from the outer `J_GV`.

---

## Why the combined-J / KO-dim sign table cannot see the fork

The combined real structure is `J = J_GV ⊗ J_F`. Its sign table is governed
by:

- `J_GV` (outer, GeoVac S³): `J_GV² = −1`, KO-dim 3 — **fixed by Paper 32
  Prop. `reality`** (genuine GeoVac datum, not free).
- `J_F` (finite, doubled): the **matter↔antimatter swap-conjugation**
  `U_F = [[0,I],[I,0]]`, `J_F² = +1`, KO-dim 6 — the **same for both
  candidate real forms** (Part 5).

The combined values (product rule, KO-dim adds mod 8):

| Candidate | internal `J_int²` | div-ring | `J_F²` | finite KO | combined KO | combined (ε, ε') |
|:----------|:-----------------:|:--------:|:------:|:---------:|:-----------:|:----------------:|
| **ℂ⊕ℍ⊕M₃** (ℍ) | **−1** | YES | +1 | 6 | 1 | (−1, +1) |
| **ℂ⊕M₂(ℂ)⊕M₃** (M₂) | **+1** | NO | +1 | 6 | 1 | (−1, +1) |

**The combined J², combined ε', and combined KO-dim are IDENTICAL for both
candidates.** They differ ONLY in the internal C² involution sign — which
the combined-J sign table never references, because the combined J is built
from `J_F` (the doubled swap), whose square `+1` is blind to the internal
quaternionic structure on the C² doublet.

Numerically confirmed on the actual operators (Part 3, ℂ⊕ℍ triple, real
`real_structure.build_J_full_dirac` ⊗ `ElectroweakFiniteTriple.J_F`):

```
n_max=1: J² → ε=−1 (|J²+I| = 0.00e0),  JD → ε'=+1 (|U·D̄ − D·U| = 0.00e0)
n_max=2: J² → ε=−1 (0.00e0),           JD → ε'=+1 (0.00e0)
n_max=3: J² → ε=−1 (0.00e0),           JD → ε'=+1 (0.00e0)
```

bit-exact, matching the table-predicted combined KO-dim 1. But this is the
**same** sign whether the n=2 factor is ℍ or M₂(ℂ).

---

## The KO-dim "forcing" is real but does NOT reach the division ring

Part 4: holding the outer KO-dim = 3 fixed, the SM-convention requirement
"combined KO-dim ≡ 1 (mod 8)" gives **`d_F = 6` uniquely**:

```
d_F giving combined KO-dim = 1 (mod 8): [6]
```

This looks like a forcing, and Door 4b hoped it would close the fork. It
does **not**, for two compounding reasons:

1. **The target "combined = 1 mod 8" is an IMPOSED SM input, not a GeoVac
   output.** Nothing in the Bertrand × Hopf-tower construction selects the
   combined KO-dim (Door 4b Q3: inner KO-dim is FREE). So even the `d_F = 6`
   conclusion is downstream of importing the SM total.

2. **Even granting `d_F = 6`, KO-dim 6 admits BOTH real forms.** KO-dim 6 is
   the property of the *doubled* `J_F` (`J_F² = +1`, swap-conjugation). It is
   satisfied identically whether the C² doublet carries `j_q` (ℍ) or `j_r`
   (M₂). The division-ring distinction is invisible to `d_F`.

The free link Door 4b flagged (inner KO-dim free) is therefore **not the
only gap**: even with the link supplied (`d_F = 6` imposed), the sign table
still does not separate ℍ from M₂(ℂ).

---

## The literature's actual selector (order-one with Yukawa) — also does not fire here

The CCM selection of ℍ over M₂(ℂ) is **not** a combined-J statement. It is
the **second-order / order-one condition** applied to a **nonzero
off-diagonal `D_F`** (Yukawa-carrying finite Dirac): with a nonzero `D_F`,
order-one `[[D, a], J b J⁻¹] = 0` is supposed to constrain the algebra
action.

Part 3c tested this directly (generic nonzero Yukawa y = 0.3, used ONLY to
probe the structural constraint — never to fix a value):

```
kind=H  n_max=1: order1_max = 0.00e0     kind=M2 n_max=1: order1_max = 0.00e0
kind=H  n_max=2: order1_max = 0.00e0     kind=M2 n_max=2: order1_max = 0.00e0
```

**Order-one is satisfied by both ℍ and M₂(ℂ)** in this construction.
Mechanism (verified directly): `J b J⁻¹` is supported entirely on the
**antimatter** sector (`‖(JbJ⁻¹)_matter‖ = 0`, `‖(JbJ⁻¹)_antimatter‖ ≠ 0`),
while `a` and `[D, a]` are supported on the **matter** sector
(`D_F` is block-diagonal in matter/antimatter, off-block norm 0). So
order-one vanishes by **matter/antimatter decoupling**, regardless of
whether the L-block is quaternionic or general. The general M₂ action and
the quaternionic ℍ action are genuinely different (`H action == M2 action?
False`), yet both pass.

This *reinforces* the ADMITS verdict: in the present (standard,
matter-doubled) construction, **neither** the combined-J sign table **nor**
the order-one condition selects ℍ. The selection would require the finer
CCM machinery (second-order condition + complex chiral fermion
representation classification, à la Chamseddine–Connes 2008 "Why the
Standard Model"), which:

- operates **entirely on finite-side data** (`D_F`, `J_F`, the algebra
  representation), and
- **never references `J_GV`**.

So even if that finer argument fires, it is the **literature's** ℍ-selection
transported into the finite factor — **not** a GeoVac-native forcing from
`J_GV² = −1`. The probe's specific target ("close the fork FROM GeoVac's
native `J_GV`") is unmet either way.

---

## The decisive distinction (audit-discipline: "consistent" ≠ "forces")

> **FORCES ℍ** would require: combined-J / KO-dim consistency, with outer
> `J_GV² = −1` fixed by GeoVac, admits ONLY the quaternionic completion and
> EXCLUDES M₂(ℂ).
>
> **ADMITS ℍ** is what we find: both real forms produce the identical
> combined-J sign table, identical combined KO-dim 1, and both pass
> order-zero AND order-one. The outer `J_GV² = −1` is **consistent with
> both**; it does not exclude M₂(ℂ).

The conflation the audit guards against would be: "the outer `J_GV² = −1` is
the same quaternionic sign that selects ℍ in CCM, therefore GeoVac forces
ℍ." This is exactly the move that fails. The signs *coincide* (both are
"the pseudoreal Spin(3) −1"), and that coincidence is genuinely suggestive —
but the combined-J bookkeeping built from `J_F` (the doubled swap) cannot
transmit the outer `J_GV` quaternionic sign to a constraint on the
*internal* C² involution that actually distinguishes the division ring. The
two −1's live on different objects (`J_GV` on the S³ spinor; `j_q` on the
finite C² doublet) and the construction provides no morphism welding them.

---

## J sign-table for both candidates (the deliverable)

```
                       internal   div-    doubled    finite   combined   combined
candidate    n=2 ring   J_int²    ring    J_F²       KO-dim   KO-dim     (eps, eps')
-------------------------------------------------------------------------------------
C(+)H(+)M3     H         -1       YES      +1          6         1        (-1, +1)
C(+)M2(+)M3    M2(C)     +1       NO       +1          6         1        (-1, +1)
-------------------------------------------------------------------------------------
DIFFER ONLY in the internal C^2 involution sign (and hence division-ring-ness).
IDENTICAL in everything the combined-J / KO-dim sign table can see.
```

All sign-table entries integer / ±1; **no transcendental appears**
(integer KO-dim and ±1 sign bookkeeping throughout — transcendental-tag
discipline: nothing to tag).

---

## Precise KO-dim consistency argument

1. GeoVac fixes outer KO-dim = 3, `J_GV² = −1` (Paper 32 Prop. `reality`,
   bit-exact at n_max ∈ {1,2,3} on Weyl and full-Dirac sectors). This is a
   genuine GeoVac datum.
2. The finite real structure `J_F` is the matter↔antimatter swap-conjugation,
   `J_F² = +1`, KO-dim 6. This is **identical** for the ℍ and M₂(ℂ) real
   forms — both double the *same* 4-dim matter Hilbert space the *same* way.
3. Combined KO-dim = 3 + 6 = 9 ≡ 1 (mod 8); combined (ε, ε') = (−1, +1).
   **Same for both candidates** (Parts 4, 5; bit-exact numerics Part 3).
4. The ℍ/M₂(ℂ) distinction is `J_int²` on the internal C² doublet
   (`−1` ⇒ ℍ, `+1` ⇒ M₂(ℂ)), which is **not** a factor of the combined
   `J = J_GV ⊗ J_F` and **not** constrained by the combined KO-dim.
5. ∴ outer `J_GV² = −1` **does not exclude** M₂(ℂ). The fork stays open.

---

## Inner-algebra status after Door 4c

```
FORCED (Bertrand × Hopf-tower, same principle as gauge content):
  • factor count = 3
  • Wedderburn / matrix-algebra type of each factor
  • n=1 factor = C
  • n=3 factor = M_3(C)

ADMITTED-NOT-FORCED (Door 4c verdict; the J_GV handle does NOT close it):
  • n=2 factor = H (vs M_2(C))
      handle tested:  combined J = J_GV ⊗ J_F sign table + KO-dim consistency
      result:         IDENTICAL sign table for both; order-one passes for both;
                      H is the CCM complex-fermion-rep import, NOT a J_GV forcing
      why it fails:   the division-ring distinction is an INTERNAL C^2 involution
                      sign (j_q^2 = -1 vs j_r^2 = +1) on the ALGEBRA side, which
                      the combined-J (built from the doubled J_F^2 = +1) cannot see

FREE (no construction reaches them):
  • generation count N_gen   (invisible to inner automorphisms)
  • inner KO-dimension       (imposed to hit SM total 1 mod 8)
  • Yukawa eigenvalues       (seam theorem; out of scope)
```

Net: **the inner algebra remains PARTIAL.** Door 4b's "sprint-reachable
door" (Falsifier A) is now **tested and closed negative** at the sign-table
level. The fork is not closable by `J_GV` through the combined-J / KO-dim
bookkeeping. Whether a *different* GeoVac-native object (not the combined-J
sign table) could weld the outer pseudoreal sign to the internal C²
involution is a deeper, un-handled question — but the specific handle Door 4b
proposed (combined-J / KO-dim consistency) does not do it.

---

## Honest scope / what was NOT shown

- This does **not** prove M₂(ℂ) is the *right* choice — it isn't; the
  physical SM uses ℍ, and the CCM finite-triple second-order-condition
  argument selects ℍ. What is shown is that **GeoVac's `J_GV` does not
  supply that selection**; GeoVac would import it.
- This does **not** rule out some *other* (non-sign-table) GeoVac-internal
  argument forcing ℍ. It rules out the specific Door 4b candidate (combined-J
  / KO-dim consistency). The deeper "second packing axiom for inner-factor
  data" (Paper 18 §IV open question) is untouched and remains the real wall.
- The order-one test used a generic nonzero Yukawa **only to probe the
  structural order-one constraint** — no Yukawa value was selected,
  proposed, or read off. Guardrail-clean (H1 / W3 / Koide negatives,
  CLAUDE.md §3, untouched).

---

## Files / sources

- Driver: `debug/door4c_j_signtable_audit.py`
- Data: `debug/data/door4c_j_signtable_audit.json`
- Prior memos: `debug/door4b_inner_algebra_forcing_memo.md`,
  `debug/door4_gauge_yukawa_boundary_memo.md`
- Code: `geovac/real_structure.py` (`build_J_full_dirac`, `J_GV`, KO-dim 3),
  `geovac/almost_commutative.py` (`ElectroweakFiniteTriple`, `J_F` KO-dim 6,
  `quaternion_to_matrix`), `geovac/standard_model_triple.py`
- Paper 32 §VIII: `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
  — Prop. `reality` (`J_GV² = −1`, KO-dim 3, pseudoreal Spin(3) spinor,
  lines 893–933); Q1 "order one beyond rank 1" (ℍ ≅ M₂(ℂ) over ℂ,
  lines 4605–4622); §VIII.C G3 + H1.
- Literature: Chamseddine–Connes 2008 "Why the Standard Model" J. Geom. Phys.
  58, 38 (finite-geometry classification, ℍ selection); Chamseddine–Connes–
  Marcolli 2007 Adv. Theor. Math. Phys. 11, 991; Connes 1995 / van Suijlekom
  2015 (KO-dim sign tables; ℍ selected by complex-rep + second-order
  condition).
