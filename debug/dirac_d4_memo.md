# Track D4 Memo — Orduz S¹-equivariant Dirac decomposition on the Hopf fibration

**Sprint:** Dirac-on-S³ Tier 1 (Track D4).
**Date:** 2026-04-15.
**Dependencies:** D1 (`geovac/dirac_s3.py`).
**Algebraic-first:** all results are sympy-exact. No floating-point diagonalization. Charges and multiplicities are integers / half-integers.

---

## 1. Setup

The Hopf fibration S¹ → S³ → S² admits a U(1) action by rotation along the fiber. The Spin(4) = SU(2)_L × SU(2)_R decomposition of the Dirac spinor bundle on S³ at Camporesi–Higuchi level n_CH = n assigns the two chirality sectors to irreps

```
(+)-chirality :  (j_L, j_R) = ((n+1)/2,  n/2 )     dim (n+1)(n+2)
(-)-chirality :  (j_L, j_R) = ( n/2,  (n+1)/2)     dim (n+1)(n+2)
```

The Hopf circle sits as the maximal torus of SU(2)_R (the opposite choice sits in SU(2)_L; the resulting decompositions are mirror images and give the same answers to the structural questions below). The U(1)-charge of an irrep (j_L, j_R) runs

```
q = m_R  ∈  { -j_R, -j_R+1, ..., +j_R }    (step 1)
mult at each q  =  (2 j_L + 1)
```

Half-integer charges arise naturally: at odd n the (+)-chirality has j_R = n/2 (half-integer); at even n, (-)-chirality has j_R = (n+1)/2 (half-integer). The full Dirac support at level n_CH therefore has charges

```
q  ∈  ½ℤ  ∩  [-(n+1)/2,  (n+1)/2]
```

with multiplicities alternating between (n+1) and (n+2) depending on whether q is integer (from (-)-chirality, mult n+1) or half-integer (from (+)-chirality, mult n+2) — or vice versa at even n.

This is the Orduz / Bär–Dahl "tower of charge-q twisted Dirac operators on S²": each charge sector on S³ equals (up to the ℏω_Hopf shift) a twisted Dirac operator on S² with line-bundle twist q.

## 2. Explicit (q, multiplicity) partitions for n_CH = 0..5

Keying by 2q (integer) for type safety; q = (2q)/2 ∈ ½ℤ. Full Dirac sector only (sum of both chiralities).

| n_CH | g_n^Dirac | partition (2q → mult)                                                    | charge support |
|:----:|:---------:|:-------------------------------------------------------------------------|:---------------|
|  0   |  4        | {−1: 1, 0: 2, 1: 1}                                                      | q ∈ {−½, 0, ½} |
|  1   | 12        | {−2: 2, −1: 3, 0: 2, 1: 3, 2: 2}                                         | q ∈ {−1,…,1}   |
|  2   | 24        | {−3: 3, −2: 4, −1: 3, 0: 4, 1: 3, 2: 4, 3: 3}                            | q ∈ {−3/2,…,3/2} |
|  3   | 40        | {−4: 4, −3: 5, −2: 4, −1: 5, 0: 4, 1: 5, 2: 4, 3: 5, 4: 4}               | q ∈ {−2,…,2}   |
|  4   | 60        | {−5: 5, −4: 6, −3: 5, −2: 6, −1: 5, 0: 6, 1: 5, 2: 6, 3: 5, 4: 6, 5: 5}  | q ∈ {−5/2,…,5/2} |
|  5   | 84        | {alt 6,7,…,7,6 over 2q ∈ [−6,6]}                                         | q ∈ {−3,…,3}   |

Sanity check: Σ mult = g_n^Dirac = 2(n+1)(n+2) at every level (verified symbolically).

Hopf-charge Casimir moment at each level:

| n_CH | Σ q²·mult |
|:----:|:---------:|
|  0   |  1/2      |
|  1   |  11/2     |
|  2   |  23       |
|  3   |  65       |
|  4   |  295/2    |
|  5   |  581/2    |

---

## 3. Tests against Paper 2 targets

### 3.1 Δ = 1/40

**Reproduced, no new content.** Phase 4H SM-D already identified Δ⁻¹ = g_3^Dirac = 40. The Hopf decomposition refines the reading as follows:

- **By chirality:** 20 (+) ⊕ 20 (-) = 40.
- **By charge parity:** 20 half-integer-charge states ⊕ 20 integer-charge states = 40.
  - Half-integer: 4 charges × 5 mult = 20.
  - Integer: 5 charges × 4 mult = 20.

Paper 2's canonical factorization is Δ⁻¹ = |λ₃|·N(2) = 8·5 = 40 (scalar Laplacian eigenvalue |λ| = n²−1 at n_Fock=3 gives 8; N(2) = cumulative state count). The Dirac-Hopf factorizations are 4·5 + 5·4 or 2·20. Neither reproduces 8·5 as a simple product — the "8" in Paper 2 is a scalar eigenvalue, not a Dirac/spinor object. **No new structural handle on Δ from D4.**

### 3.2 B = 42 — NUMERICAL NEAR-MISS, NO CLEAN IDENTITY

The cleanest truncated Dirac Casimir trace at Paper 2's m=3 cutoff (n_CH ≤ 2):

```
Σ_{n=0}^{2} |λ_n|² · g_n^Dirac
 = (3/2)² · 4  +  (5/2)² · 12  +  (7/2)² · 24
 = 9 + 75 + 294
 = 378
 = 9 · 42.
```

**This is a numerical multiple of 42.** But the factor 9 has no clean spectral interpretation in the Hopf/Dirac setting — it is not any single |λ|², j_R, or 2j+1. Extending to n_CH ≤ 3 gives

```
Σ_{n=0}^{3} |λ_n|² · g_n^Dirac  =  378 + (9/2)² · 40  =  378 + 810  =  1188  =  198 · 6  =  not a multiple of 42.
```

The "9·42" is a coincidence at the specific cutoff n_CH ≤ 2. It does NOT promote to a B = 42 identity.

No Hopf-charge-truncated sum (Σ_{|q|≤Q, n≤N}) reproduces 42 as an exact hit — confirmed by the brute-force search in the script.

Conclusion: **B does not lift to the Hopf-Dirac sector as a structural identity.**

### 3.3 F = π²/6 — CLEAN NEGATIVE

The Dirac Dirichlet series at the packing exponent s = d_max = 4 admits a closed form via half-integer Hurwitz ζ. Letting k = n + 3/2 run over half-integers ≥ 3/2:

```
D_dirac(s)  =  Σ_{n≥0} 2(n+1)(n+2) · (n + 3/2)^{−s}
           =  2 · ζ_{3/2}(s−2)  −  (1/2) · ζ_{3/2}(s)
```

where ζ_{3/2}(s) = (2^s − 1)·ζ(s) − 2^s (half-integer Hurwitz identity).

At s = 4 this evaluates **symbolically** to

```
D_dirac(4) = π²  −  π⁴/12  ≈  1.7522
```

vs F = π²/6 ≈ 1.6449. The expression contains ζ(2) = π²/6 and ζ(4) = π⁴/90 in an inseparable rational mixture; there is **no sympy simplification isolating F as a multiple.**

At s = 3, the series diverges (ζ(1) pole). At s = 5, 6, closed forms exist (mixtures of ζ(3),ζ(5) and ζ(4),ζ(6) respectively); none are multiples of F.

Per-charge-sector Dirichlet series: each fixed q sector supports only n ≥ 2|q|, so the sum is a tail-shifted variant of the same half-integer Hurwitz. No charge-weighted sum produces F either.

Conclusion: **F does not lift to the Hopf-Dirac sector.**

---

## 4. Categorical distinction from Phase 4B α-D (scalar Hopf graph morphism)

This is not a re-derivation. The distinction is structural and must be preserved in the D5 writeup.

**Phase 4B α-D** tested the scalar Hopf quotient S³ → S² at the graph level: nodes (n, l, m) mapped to S² sectors labeled by scalar spherical-harmonic content. The Laplace–Beltrami is second-order and bosonic; the quotient is a deformation-retract-like averaging; no line bundle can appear. All charge content is trivial (only q = 0). Negative result: no spectral invariant of the quotient hits K targets.

**D4** tests the Dirac decomposition. The Dirac operator is first-order; the spinor bundle is nontrivial; the U(1) Hopf action gives a **genuine direct-sum decomposition over half-integer charges q**, with each charge sector a twisted Dirac operator on S² carrying a nontrivial Baum–Friedrich / Bär spectrum and a nontrivial line-bundle O(2q). This is genuinely new content unavailable in the scalar case.

The categorical difference manifests concretely as:

- The charge support at level n_CH is {−(n+1)/2, ..., (n+1)/2} with both integer AND half-integer values. The scalar case has only q = 0.
- The chirality/charge-parity split (+ ↔ half-int ; − ↔ int, or vice versa at even n) reproduces Δ⁻¹ = 40 as 20 + 20 — a genuinely new factorization mirror of the trivial 2·20 chirality split, impossible in the bosonic sector.

**However**, while D4 opens a new decomposition, it does NOT open new structural handles on B or F. The numerical "378 = 9·42" at n_CH ≤ 2 is a coincidence, not a theorem; D_dirac(4) = π² − π⁴/12 has no clean relationship to ζ(2).

## 5. Verdict for D5

| Target | D4 result | Status |
|:------:|:----------|:-------|
| Δ = 1/40 | Reproduced; new 20+20 charge-parity split is structurally cleaner than Paper 2's 8·5, but not a new identity | Known since Phase 4H |
| B = 42  | 9·42 near-miss at n_CH ≤ 2; no clean factor interpretation; drops at n_CH ≤ 3 | **Negative** |
| F = π²/6 | D_dirac(4) = π² − π⁴/12, inseparable mixture | **Negative** |

**D5 verdict recommendation:** NO Hopf-Dirac lift of the K combination rule.

Mapped onto the sprint-plan verdict table (dirac_s3_tier1_sprint_plan.md §D5):

- If D2 (Dirac B = 42) and D3 (Dirac F = π²/6) also come back negative, the verdict is **"All negative" → three-tier coincidence formally documented; Paper 2 §IV closure statement; α remains conjectural with structural-mystery framing.**
- If D2 and/or D3 are positive, D4 alone being negative maps to **"D2 + D3 positive, D4 negative: B, F both lift to Dirac spectrum; Hopf structure not the link."**

In all cases D4 rules out the Hopf-equivariant Dirac sector as the common generator.

**Recommended next steps if this verdict holds jointly with D2/D3:**

1. Close the alpha sprint series at Paper 2 §IV with the three-tier structural-mystery framing.
2. Move on to Tier 2 (spin-ful composed qubit encoding), which benefits from D1 regardless of whether the α program closes positively or negatively.

---

## 6. Deliverables produced

- `debug/dirac_d4_hopf_equivariant.py` (symbolic, sympy-exact).
- `debug/data/dirac_d4_hopf_equivariant.json` (all partitions, moments, tests, verdict).
- This memo.

## 7. What was NOT done (guardrails respected)

- No extension to S⁵ / S⁷ (Section 3 Phase 4E, α-I).
- No discrete / graph-level Dirac construction (Ginsparg–Wilson; Section 3 D1 guardrail).
- No re-run of scalar Hopf graph morphism (Section 3 Phase 4B, α-D).
- No edits to `geovac/dirac_s3.py` (D1 is frozen).
- No edits to Paper 2 (deferred until D5 per PI decision 3).
