# Sprint memo — decided census, transcendental tagging, polyatomic scoping

**Date:** 2026-08-13
**Versions:** v4.78.0 (decided census + tagging), v4.79.0 (polyatomic scoping)
**Branch:** `work/sparsity-boundary` (local only; not pushed — repo frozen, PI controls release)
**Owning doc:** `docs/neumann_general_m_build_plan.md`; papers 58, 18, 34

One memo, three subsections. They are one arc: the two-center ERI engine closed at
transcendence weight one and π-free (v4.77.0), and everything here is either
cashing that in or asking what it does not reach.

---

## 1. Step 2 — the (AA|BB) block goes from *counted* to *decided*

**The gap.** Paper 58's `g` row is a **counted** census: it states which entries
the selection rules *permit*, and the paper is explicit that these are "claims
about how many entries the selection rules permit, not verifications that each
permitted entry is nonzero." Its Gaussian corroboration decides zeros by a
`1e-10` float threshold, which is a measurement, not a decision.

**Why it is now closable.** Increment 1c put (AA|BB) in closed form, every π
cancels, and the entry has the shape

```
(ab|cd) = A_0(R) + Σ_j A_j(R) e^{−λ_j R},    A_j rational, λ_j rational
```

For algebraic `R` the family `{1, e^{−λ_j R}}` is linearly independent over the
algebraics (Lindemann–Weierstrass), so the entry vanishes **iff** every `A_j`
does. Group the expanded form by exponential rate, ask whether each rational
coefficient vanishes. That is a decision procedure.

**Three outcomes were pre-separated deliberately**, because the third is
invisible to *both* prior methods:

| outcome | what it is | who can see it |
|:---|:---|:---|
| nonzero | decided nonzero | both |
| zero, Gaunt | permitted by the crude count, but the Gaunt coefficient vanishes | thresholding only |
| zero, **accidental** | Gaunt nonzero, radial combination cancels anyway | **neither** |

**Result** (Paper 58's own census configuration: Z_A=3, Z_B=1, n_max=2, R=3 bohr;
driver `debug/step2_decided_census.py`, ~110 s):

```
entries in the block                                   625
permitted by the counting rules                        195

  decided nonzero                                      195   (100.0%)
  zero, vanishing Gaunt coefficient                      0   (  0.0%)
  zero, ACCIDENTAL                                       0   (  0.0%)
```

All 195 permitted entries are genuinely nonzero. **On this block the counted
density is not an overcount — it is exactly the true density**, with no threshold
anywhere.

This also settles what the Gaussian corroboration could only suggest: it found
29.8% against 29.4% counted and attributed the gap to the real-Cartesian vs
complex-`m` basis convention rather than to over-permissive counting. There is no
over-permissiveness here to find.

**Scope.** The `g` row's overall tier is **unchanged**. (AA|BB) plus the
one-centre classes are the 20.5% carrying residual `l`-selection; the three
genuinely cross classes are 80% of the tensor and stay *counted* until an
independence argument over {E₁, ln, γ} exists. What changed is that the sector
Thm. `abelian_residue` makes a positive prediction about is now decided.

---

## 2. Transcendental tagging — discharged, and load-bearing

`feedback_tag_transcendentals` fires whenever a transcendental appears. `ln` and
γ_E appeared across increments 2, 3b, 3c, 3d; the rule was flagged four times and
discharged zero times. Closed here.

| object | Paper 18 tier | why it is not an independent import |
|:---|:---|:---|
| e^a E₁(a) | **embedding** | already catalogued, Level-2 Stieltjes seed |
| ln | **embedding** | *is* the Neumann kernel Q₀(ξ) = ½ln[(ξ+1)/(ξ−1)] — the same prolate-spheroidal curvature seen from the kernel side rather than the moment side |
| γ_E | **embedding** | not independent of E₁ at all; arrives with it via E₁(z) = −γ_E − ln z + Ein(z) |

**Mellin sub-mechanism: none applies.** M1/M2/M3 are the calibration tier's
machinery, and there is no calibration constant here to explain.

**Paper 34 chain:** Fock conformal → Sturmian → the multipole/Gaunt projection
evaluated *across two centres*. Three-axis tag (L, dimension-preserving,
{E₁, ln, γ_E} at weight one).

**The two absences are the result.**

- **No π survives.** It enters through the spherical normalisations and cancels
  identically against the 4π of the Coulomb kernel.
- **No weight-two object appears** — no Li₂, no ζ(2) — although the ordered
  ξ_</ξ_> integral is formally an iterated integral over a simplex, exactly the
  shape that generically produces dilogarithms.

Together those are *precisely* what makes §1 possible. With no π and no
weight-two content the closed form is rational-coefficient exponential-sum and
Lindemann applies. Had either survived, the decided census could not have been
run. **Cataloguing the seed set bought the zero test rather than merely
labelling it.**

**Scope correction shipped.** Paper 34's multipole/Gaunt entry asserted flatly
that "no transcendental is introduced at any step." That is true on one centre,
where the radial factor is a Slater integral; it is false across two centres,
where it is not. Sentence scoped; two-centre continuation stated.

---

## 3. Polyatomic scoping — Poly-0 and Poly-1

Diagnostic before engineering, on the PI's question "can we plan to solve water?"

### Poly-0: the three-centre burden

Water has three nuclei, so its ERI tensor splits by how many *distinct* centres
the four orbital indices touch. Driver `debug/poly0_three_center_burden.py`.

| system | 1-centre | 2-centre | 3-centre | **cost of dropping 3-centre** |
|:---|---:|---:|---:|---:|
| H₂ (2 nuclei, null control) | 18.0% | 82.0% | 0.0% | **+0.00000000 Ha** |
| BeH₂ (3 nuclei, linear control) | 25.1% | 64.9% | 10.0% | −0.54180387 Ha |
| **H₂O (3 nuclei, bent)** | 31.1% | 55.2% | 13.7% | **−2.34720266 Ha** |

The percentage columns are the misleading ones; the last is the gate. Water's
3-centre block is 13.7% of Σ|g| and 420 of 2401 entries — but dropping it costs
**1467× chemical accuracy**. The sign is a variational catastrophe: the energy
goes *down*, below the true value, because the 3-centre terms are largely
repulsive and removing repulsion over-binds. **There is no truncation story.**

The H₂ row returning exactly +0.0 is the control on the partition machinery.

Structural consolation: water has only three nuclei, so it has **no 4-centre
integrals at all**. Water needs exactly one new capability, not two.

Method note: both columns come from the *same* McMurchie–Davidson reference
tensor, partitioned. No fit error enters the comparison — the dropped-energy
figure is a statement about the partition alone, not about the basis.

### Poly-1: the rotation gate

The engine is inherently axial (prolate spheroidal puts both centres on z); water's
O–H bonds sit at ±52.25°. Testing the standard rotation route — evaluate in the
frame where the pair axis *is* z, rotate each index back with a Wigner-D, the same
trick `shibuya_wulfman` already uses for the **one**-body cross-centre integral.
Driver `debug/poly1_rotation_gate.py`.

```
max |g_lab − (D ⊗ D ⊗ D ⊗ D)·g_axis|  =  2.2e-16    over 5 orientations
```

with an l=1 shell in the basis, so the (x,y,z) index bookkeeping is genuinely
exercised — a wrong-basis bug could not hide.

**Net: 1981 of 2401 entries (86.3% of Σ|g|) are reachable exactly today, at any
geometry, with what this arc already built.** The gap is precisely the 420
three-centre entries, and it cannot be dropped.

### Why 3-centre is hard

Prolate spheroidal coordinates are built from exactly *two* foci. A third nucleus
has nowhere to sit. Gaussians escape this because two Gaussians on different
centres multiply into a single Gaussian on a third point; Slater functions have no
product theorem. This is the 70-year-old reason Gaussians won quantum chemistry,
not a GeoVac-specific wall.

### Routes, priced

| route | verdict |
|:---|:---|
| **A.** One-centre re-expansion (Löwdin / Barnett–Coulson) | ⚠️ **GUARDRAIL** (§3.5, Papers 8–9). Truncation is in `l`, destroying the exact angular sparsity the framework is built on — the polyatomic replay of the Löwdin-retrofit dead end already in §3. Also re-enters Cor. `dual_p0` (no shared p₀ for heteronuclear; water is O + H). Works numerically, costs the framework its identity. |
| **B.** Gaussian transform (Shavitt–Karplus) | Exact, no `l`-truncation, preserves the basis — but leaves a numerical integral per ERI, forfeiting closed form, Lindemann decidability, and the compiled-evaluation speed. Sound chemistry; converts the distinctive product into an ordinary one. |
| **C.** Momentum space / Fourier (Sturmian-native) | Coulomb kernel is 4π/k², translation is a phase e^{ik·R}; three centres = three phases. This *is* the Fock projection, and the one-body 3-centre analog is already solved in-repo. Whether the two-body case closes is **genuinely open**. Most GeoVac-native route; the right question for the Avery call. |
| **D.** Decide rather than solve | Paper 58's Prediction `angular` (C₂ᵥ abelian of order 4 ⇒ **2-bit spatial grading, two qubits, not a multiplicative factor**) is a *support* claim, not an energy claim. Cheapest real deliverable. **Recommended first.** |

Note on D: Paper 58 marks that prediction "not falsifiable on the present
builder" because the composed builder carries no bond angle (`MolecularSpec.nuclei
= None`, Obs. `no_angle`). That obstruction is about the *composed builder*, not
about the framework — a geometry-carrying reference tensor was built here.

---

## 4. Files

**Created:** `debug/step2_decided_census.py`, `debug/poly0_three_center_burden.py`,
`debug/poly1_rotation_gate.py`, this memo.

**Modified:** `papers/group2_quantum_chemistry/paper_58_abelian_residue.tex`
(decided-census passage), `papers/group3_foundations/paper_18_exchange_constants.tex`
(new §Level-2 seed set), `papers/group6_precision_observations/paper_34_projection_taxonomy.tex`
(multipole/Gaunt two-centre scope), `tests/test_two_center_eri_aabb.py` (+1 test),
`docs/claim_test_matrix.md` (`g` row rewritten, +2 rows), `docs/neumann_general_m_build_plan.md`
(§10 polyatomic scoping), CHANGELOG, CLAUDE.md.

---

## 5. Verification

- Regression: 110 passed / 4 skipped (`test_two_center_eri_aabb.py` +
  `test_fock_projection.py` + `test_fock_laplacian.py`), 78 s.
- New test: `test_decided_census_aabb_block_has_no_accidental_zeros` — 4-entry
  sample chosen where accidental cancellation is likeliest (l>0 both sides, m≠0);
  the full 195-entry sweep is ~110 s and lives in the driver.
- Papers 18, 34, 58 all compile three-pass EXIT=0. The undefined-reference
  warnings in the P18/P34 logs are **pre-existing** and untouched by these edits
  (surfaced and reported separately to the PI).
- No production `geovac/` code changed in §3; Poly-0/1 are pure diagnostics.
- Hard prohibitions (§13.5): nothing touched.

---

## 6. Honest scope

**Closed at theorem grade.**
- (AA|BB) zero-decidability. The Lindemann separation is a proof, not a
  measurement, and it rests on the π-cancellation and weight-one closure, both
  themselves proven earlier in the arc.
- The rotation transport of the axial engine (Poly-1) is a representation-theory
  identity; 2.2e-16 is verification of correct implementation, not of the claim.

**Structural sketch / classification, not theorem.**
- The transcendental tier assignments. "γ_E is not independent of E₁" is exact
  (the Ein identity); "`ln` is the same embedding object as the Stieltjes seed"
  is a *structural reading* of where the two come from, not a proof that they
  generate the same ring.
- The route pricing in §3. Routes A/B/C are literature-informed judgments about
  cost, not measured results.

**Numerical observation.**
- The 195/195 census is exhaustive **on one configuration** (Z_A=3, Z_B=1,
  n_max=2, R=3). It is a decided result *there*; that no accidental zeros exist
  at other Z, n_max, or R is not established and was not claimed.
- Poly-0's percentages and energies are basis-dependent (minimal Slater-shape,
  8-Gaussian fits). The *qualitative* verdict — 3-centre cannot be dropped — is
  robust across all three systems and the 1467× margin; the specific 2.35 Ha is
  not a converged number.

**Named open follow-ons.**
1. The three cross classes (80% of the tensor) remain **counted**. Deciding them
   needs an independence argument over {E₁, ln, γ} — the natural next
   decidability target, and strictly harder than (AA|BB).
2. Accidental zeros at other (Z, n_max, R). Cheap to sweep; not swept.
3. **Three-centre ERIs** — the polyatomic wall. Route C (momentum space) is open
   and is the Avery-call question.
4. Paper 58 Prediction `angular` (C₂ᵥ 2-bit grading) — testable now, not tested.
5. Step-1 driver still not migrated into `tests/` (one OWED row in the claim
   matrix).
6. Benchmarking Rule item 5 ("H2 Full CI < 1.0%") names no live artifact;
   `production_suite` H₂ sits at 2.80%. Flagged to PI, unresolved.

**What this sprint does NOT establish.** That native water is worth building. This
session's own results say the engine has no accuracy advantage over Gaussians, no
qubit/Pauli advantage (QC-1 negative), and the 420× speed figure (QC-2) is against
in-repo pure Python, not a production C code. The case for water is *structural*
— an exact, decidable sparsity statement for a bent molecule — not computational.
