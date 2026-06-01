# Door 1b — S¹¹ graduation of the F-theorem ladder recursion + analytic c_d

**Date:** 2026-06-01
**Driver:** `debug/door1b_s11_graduation.py`
**Data:** `debug/data/door1b_s11_graduation.json`
**Predecessor:** `debug/door1_ftheorem_odd_d.py` (Door 1, found the borderline recursion)
**Verdict:** **WALL** — the ladder recursion is a single-rung (S⁵) coincidence, not a
dimensional recursion. It breaks at S⁷, S⁹, S¹¹ under any *uniform* Dirac
normalisation. The original "c_d = 1/2 at d ≥ 7" was an artifact of an internally
inconsistent multiplicity convention. The breakdown is *derived analytically* (top-atom
mismatch 2^{−(d−5)/2}), so this is a clean, understood wall rather than an empirical
failure.

---

## 1. What Door 1 claimed

Door 1 reported a BORDERLINE structural recursion among the odd-d F-theorem coefficients
`F = −½ ζ'_Δ(0)`:

```
scalar_{S^d}(0)  −  2 · Dirac_{S^d}(0)  =  c_d · Dirac_{S^{d-2}}(0)
```

with c₅ = 1, c₇ = c₉ = 1/2 (exact rational, 3 rungs). The named graduation test was
S¹¹ + an analytic derivation of c_d.

(Here `scalar` = conformally coupled scalar, eigenvalue v²−1/4, v = n+(d−1)/2;
`Dirac` = single-chirality Camporesi–Higuchi Weyl Dirac, eigenvalue |λ| = n+d/2;
both quantities are ζ'_Δ(0) on the round **unit** S^d.)

---

## 2. The load-bearing convention problem (caught by the curve-fit audit)

Door 1's Dirac multiplicity normalisation was hardcoded
`DIRAC_CMAP = {3:2, 5:12, 7:360, 9:20160}`, written `g_n = prod_{j=1}^{d-1}(n+j)/c_d`.
**This is not a single formula.** It is

| d | Door 1 c_d | genuine Weyl (d−1)!/2^{⌊d/2⌋−1} | (d−1)!/2 |
|---|-----------|--------------------------------|----------|
| 3 | 2         | **2** ✓                        | 1        |
| 5 | 12        | **12** ✓                       | 12       |
| 7 | 360       | 180                            | **360** ✓ |
| 9 | 20160     | 5040                           | **20160** ✓ |

i.e. Door 1 used the genuine Weyl multiplicity at d = 3,5 but switched to *twice* the
Weyl (= (d−1)!/2) at d ≥ 7. That factor-of-2 switch is an **untracked free parameter**,
and it is exactly the parameter that produced c₇ = c₉ = 1/2. The recursion constant c_d
is *defined relative to* the Dirac normalisation; doubling the Dirac multiplicity halves
c_d. So "c_d = 1/2" carries no convention-independent meaning.

The Paper 50 anchor: the S⁵ Dirac F-coefficient matches `c₅ = 12 = (d−1)!/2` bit-exactly,
which equals the genuine Weyl value at d=5 (coincidence: (4)!/2 = 12 = 4!/2² ). So the
physically-anchored uniform choice is the **single-chirality Weyl multiplicity**
`c_d^{weyl} = (d−1)!/2^{⌊d/2⌋−1}` = {2, 12, 180, 5040, 226800} for d = 3,5,7,9,11.

This sprint runs the test under **two uniform conventions** applied to every rung:
- `weyl`: c_d = (d−1)!/2^{⌊d/2⌋−1} (single-chirality CH; matches Paper 50 at S⁵)
- `full`: c_d = (d−1)!/2^{⌊d/2⌋} (both chiralities)

---

## 3. S¹¹ closed forms (skeleton-exact / Layer-2 identification)

All quantities computed in **exact arithmetic** for the skeleton (eigenvalues,
multiplicities, multiplicity polynomials) and identified by **PSLQ at 320 dps** for the
Layer-2 spectral-zeta value. PSLQ is gated by a reconstruction-error check (<1e-140) to
reject the spurious small-denominator relation that fooled an earlier loose-tolerance run.

**Dirac S¹¹ (weyl convention), symbolic-exact, residual = 0:**
```
ζ'_{Dir,S^11}(0) = 63/32768 log2
                 + 117469/34406400 · ζ(3)/π²
                 + 17281/2064384 · ζ(5)/π⁴
                 + 4389/327680 · ζ(7)/π⁶
                 + 935/65536 · ζ(9)/π⁸
                 + 1023/131072 · ζ(11)/π¹⁰
```

**Scalar S¹¹ (conformally coupled), Hurwitz-series + PSLQ, recon err ≈ 1e-321:**
```
ζ'_{sc,S^11}(0) = −7/131072 log2
                − 3897/45875200 · ζ(3)/π²
                − 485/8257536 · ζ(5)/π⁴
                + 609/1310720 · ζ(7)/π⁶
                + 425/262144 · ζ(9)/π⁸
                + 1023/524288 · ζ(11)/π¹⁰
```
(The scalar lead coefficient is 825753600 ≈ 8.3×10⁸; PSLQ needs maxcoeff ≥ 10¹⁰ —
the reason the earlier coarse run returned a spurious /1314007 relation.)

The odd-d F-theorem closed forms therefore **continue up the ladder** (S⁷, S⁹, S¹¹ all
have clean rational-coefficient closed forms in {log2, ζ(2k+1)/π^{2k}}). That part of
Door 1 (closed forms exist for all odd d) graduates cleanly. It is only the *cross-species
ladder recursion* that fails.

---

## 4. The recursion BREAKS at S¹¹ (and at S⁷, S⁹ under uniform convention)

Solving for c₁₁ from the precise values, c₁₁ = (scalar − 2·Dirac)/Dirac_{S⁹} ≈ 0.91893 —
not 1/2, not any clean rational near it. Per-atom ratios of `scalar_{S^11} − 2·Dirac_{S^11}`
to `Dirac_{S^9}` (weyl):

| atom | ratio |
|------|-------|
| log2 | 29/32 |
| ζ3/π² | 378239/413312 |
| ζ5/π⁴ | 92327/94752 |
| ζ7/π⁶ | 221/192 |
| ζ9/π⁸ | 57/32 |
| ζ11/π¹⁰ | **no counterpart in Dirac_{S^9}** |

The ratios all disagree, and the top atom ζ11/π¹⁰ has no home on the right-hand side.
**Under a uniform convention the recursion holds only at S⁵ (weyl, c₅ = 1) and at no
other rung.** Full convention: fails at every rung including S⁵.

| d | c_d (weyl uniform) | c_d (full uniform) |
|---|--------------------|--------------------|
| 5 | **1**              | none               |
| 7 | none               | none               |
| 9 | none               | none               |
| 11| none               | none               |

---

## 5. Analytic derivation of the breakdown (the GO-criterion answer)

The recursion `scalar − 2·Dirac = c_d·Dirac_{S^{d-2}}` requires the residual's **top atom**
ζ_d/π^{d-1} to vanish, because Dirac_{S^{d-2}} has no such atom (its top atom is
ζ_{d-2}/π^{d-3}). Necessary condition:

```
scalar_top(d) = 2 · Dirac_top(d).
```

Dirac_top(d) is fully analytic (half-integer ladder; the half-integer Hurwitz factor
(2^{−(d−1)} − 1) ζ_R'(−(d−1))). scalar_top(d) is on the integer ladder. The diagnostic
ratio is, **exactly, across all four rungs:**

```
scalar_top(d) / (2 Dirac_top(d)) = 2^{−(d−5)/2}      (weyl convention)
                                 = 2^{−(d−3)/2}      (full convention)
```

verified bit-exact at d = 5, 7, 9, 11 (PSLQ-confirmed scalar, analytic Dirac):

| d | weyl ratio | full ratio |
|---|-----------|-----------|
| 5 | 1         | 1/2       |
| 7 | 1/2       | 1/4       |
| 9 | 1/4       | 1/8       |
| 11| 1/8       | 1/16      |

The top-atom cancellation (ratio = 1) happens at **d = 5 only** (weyl) and at **no d**
(full). The mismatch is a clean power of 2 in d, growing without bound — a genuine
structural obstruction, not a near-miss. The mechanism: the scalar lives on the *integer*
ladder (eigenvalue v²−1/4), the Dirac on the *half-integer* ladder (eigenvalue u); the
relative half-vs-integer Hurwitz factor matches the "2·" prefactor at exactly one
dimension.

This is the analytic content asked for: it explains *why* d=5 is special and *why* the
sequence cannot be a closed-form constant — the would-be recursion fails the necessary
top-atom condition for every d ≥ 7.

---

## 6. Why the original saw c₇ = c₉ = 1/2

At d = 7, 9 the original used c_d = (d−1)!/2 = **twice** the genuine Weyl multiplicity.
Doubling the Dirac multiplicity halves every Dirac value, which is exactly the factor
2^{−(d−5)/2}·2 = ... that brings the d=7 top-atom ratio from 1/2 to 1 — i.e. the doubled
normalisation *artificially supplies the missing factor of 2* that lets the top atom
cancel at d=7 (and a leftover at d=9 read as a still-consistent 1/2 across the surviving
atoms). It is a normalisation artifact, not a recursion.

---

## 7. Curve-fit audit summary (full audit in the sprint thread)

- **Claim:** c_d = 1, 1/2, 1/2 is a closed-form dimensional recursion constant.
- **Free-parameter count:** the hidden, uncontrolled binary free parameter is the Dirac
  normalisation (factor-of-2 per rung). It is precisely the parameter that set c₇=c₉.
- **Selection bias:** the convention was not frozen; it was an inconsistent hardcoded dict.
  Under a frozen uniform convention the recursion fails at the first independent rung (S⁷).
- **Robustness:** the relation evaporates under a factor-of-2 normalisation change (a
  perturbation far larger than the 320-dps test precision). Maximally fragile.
- **Independent test:** S¹¹ under both uniform conventions — killed.
- **Verdict: STOP / WALL.** Demote to "single-rung S⁵ top-atom coincidence." The
  convention-invariant structural finding is the *breakdown* and its analytic cause.

---

## 8. Transcendental tagging (discipline)

All transcendentals appearing: **log 2** and **ζ(2k+1)/π^{2k}** (odd-zeta tower,
k = 1..(d−1)/2). These are the Layer-2 content of the F-theorem coefficient `−½ζ'_Δ(0)`
on odd spheres:
- **Paper 18 tier:** spectral-action / spectral-zeta (intrinsic, observation-side).
- **Master Mellin engine:** M2/M3. The odd-zeta/π^{2k} tower on odd S^d is the
  **M3 (half-integer Hurwitz / vertex-parity)** signature — Dirac on the half-integer
  ladder produces ζ(odd) via ζ_R'(−2j). log 2 is the (2^x − 1) half-integer Hurwitz
  prefactor, also M3-side.
- **Paper 34 projection:** F-theorem / spectral-zeta-derivative (Paper 50 §3 spectral-zeta
  side). No anonymous transcendentals; every atom is pinned to the half-integer Hurwitz
  mechanism.

---

## 9. What graduates and what dies

**Graduates (clean, keep):**
- The odd-d F-theorem closed forms continue to S⁷, S⁹, S¹¹ — every dimension has a clean
  rational-coefficient closed form in {log2, ζ(2k+1)/π^{2k}}. This *is* a real, verified
  continuation of Paper 50's S³/S⁵ closed forms, the engine GENERATES here.
- The analytic top-atom mechanism (ratio 2^{−(d−5)/2}) is a clean structural result
  worth one line: it pins exactly why d=5 is the only rung where the cross-species
  subtraction drops a dimension.

**Dies:**
- The cross-species ladder recursion `scalar − 2 Dirac = c_d Dirac_{S^{d-2}}` as a
  dimensional recursion. It is a single-rung (S⁵) coincidence. c_d is NOT a closed-form
  constant; under uniform normalisation the sequence is {1, —, —, —}.

**No paper edits** (per task). If captured later, the right home is a one-line remark in
Paper 50 §7 (the F-theorem-on-odd-spheres extension): the closed forms continue, but the
scalar/Dirac dual-basis / ladder reading does NOT extend past S⁵ — consistent with Paper
50 Prop 7.4 (dual-basis projection over-determined for d ≥ 5).
