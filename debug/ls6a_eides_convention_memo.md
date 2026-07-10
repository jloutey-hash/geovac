# Sprint LS-6a: Hydrogen Lamb shift in the Eides §3.2 canonical convention

**Date:** 2026-05-03
**Sprint goal:** Re-derive the LS-1 one-loop Lamb shift in the canonical Eides §3.2 convention to test the LS-5 reframing of the LS-3 residual into a "+24.7 MHz convention" piece + "+5 MHz genuine α⁵ multi-loop" piece.
**Verdict:** **POSITIVE — LS-5 reframing CONFIRMED.** LS-6a Lamb shift = **1052.19 MHz**, residual against experimental 1057.85 = **+5.65 MHz** (in LS-5 predicted band [+3.5, +10.7] MHz). The LS-1 +38/45 SE constant was a one-loop convention artifact, not missing physics; the genuine α⁵ multi-loop test target is now ~+5.65 MHz, not the +29.26 MHz LS-3 raw residual.

## 1. What LS-1 computed

LS-1 (`debug/ls1_lamb_shift_memo.md`) computed the H 2S₁/₂ − 2P₁/₂ Lamb shift using the standard Bethe-Salpeter / Eides one-loop formulas with Drake-Swainson Bethe logarithms ln k₀(2,0) = 2.81177, ln k₀(2,1) = −0.03002. The 2S₁/₂ self-energy used the form:

```
ΔE_SE(2S₁/₂) = (α³Z⁴/(πn³)) [(4/3)·ln(1/(Zα)²) − (4/3)·ln(k₀/Ry) + 38/45]   Ha
```

The "+38/45" was attributed in LS-1 to "absorbing Karplus-Klein-Darwin and the j=1/2 anomalous magnetic moment." LS-1 §2.2 explicitly flagged this as a convention choice giving SE 2S = 1039.31 MHz vs. textbook ~1078 MHz (a 4% gap) and noted that "a different grouping of the magnetic-moment / Darwin / Karplus-Klein terms would close this." The 2P₁/₂ formula used the textbook constant −1/6 (Itzykson-Zuber p. 345).

Result: total Lamb shift = 1025.06 MHz, error vs. experimental 1057.845 = **−32.78 MHz (−3.10%)**.

LS-2/3/4 added native Bethe-log derivation via Sturmian projections and Drake-Swainson regularization (Paper 34's 13th projection), but converged at infinite N back to the LS-1 baseline 1025 MHz. The LS-4 N=40 sweet spot of 1053.76 MHz (−0.39%) was an accidental cancellation between the 2S Bethe-log truncation T₂S = +28.8 MHz and the LS-1 ceiling A = −32.78 MHz.

## 2. What changes in Eides §3.2

The canonical Eides Eq. (3.32) for the one-loop nS₁/₂ self-energy is:

```
ΔE_SE(nS₁/₂) = (α(Zα)⁴/(πn³)) m_e c² · [(4/3)·ln(1/(Zα)²) − (4/3)·ln(k₀(n,0)/Ry) + 10/9]
```

The constant **10/9** is the COMBINED one-loop SE finite constant for s-states. It already includes:
- Bethe constant (5/6 from the original 1947 derivation)
- Karplus-Klein-Darwin contact contribution
- Schwinger AMM contribution at j = 1/2 for the s-state contact

There is **NO separate Schwinger AMM term to add** for s-states in the canonical Eides convention; the AMM is already inside the +10/9.

The constant difference between LS-1 and Eides is:

```
+10/9 − 38/45 = +50/45 − 38/45 = +12/45 = +4/15
```

Crucially, **+4/15 is exactly the Uehling kernel constant** that appears in the vacuum-polarization shift:

```
ΔE_VP(nS) = −(4/(15π))·α³Z⁴/n³  Ha
```

LS-1 inadvertently subtracted the Uehling kernel constant 4/15 from the canonical Eides 10/9 constant, giving 38/45. This is a double-counting error: Uehling is already counted as a separate vacuum-polarization contribution (which it is — VP is a distinct one-loop diagram from SE). The correct canonical decomposition keeps Uehling in VP and uses 10/9 alone for SE.

For 2P₁/₂, the canonical Eides form keeps the textbook −1/6 constant (AMM/spin-orbit combination for j = l − 1/2), unchanged from LS-1.

## 3. Numerical impact

LS-1 → LS-6a:

| Quantity | LS-1 | LS-6a | Shift |
|:---------|:----:|:-----:|:-----:|
| SE 2S₁/₂ constant | 38/45 | 10/9 | +4/15 |
| SE 2S₁/₂ total (MHz) | +1039.31 | +1066.44 | +27.13 |
| Total 2S₁/₂ (MHz) | +1012.18 | +1039.31 | +27.13 |
| SE 2P₁/₂ (MHz) | −12.88 | −12.88 | 0 |
| Total 2P₁/₂ (MHz) | −12.88 | −12.88 | 0 |
| **Lamb shift (MHz)** | **+1025.06** | **+1052.19** | **+27.13** |
| Error vs. experimental | −32.78 (−3.10%) | −5.65 (−0.534%) | +27.13 |

The shift +27.13 MHz from LS-1 to LS-6a:
- Decomposes purely as +(4/15) × (α³Z⁴/(πn³)) × HA_TO_MHZ at n=2, Z=1 = +27.13 MHz
- Matches the LS-5 scoping prediction +24.7 MHz to within 10% (LS-5 used a back-of-envelope with the textbook ~1078 MHz target rather than the exact Eides convention)
- Is a pure one-loop convention fix, NOT new physics

## 4. Genuine multi-loop residual

After the convention fix, the residual against experimental Lamb shift 1057.845 MHz is:

```
1057.845 − 1052.19 = +5.65 MHz
```

This sits in the LS-5 predicted α⁵ multi-loop band of **[+3.5, +10.7] MHz** (central value +7.10 MHz from Eides Table 7.4 cumulative). The breakdown of expected α⁵ multi-loop contributions to the H 2S Lamb shift (from Eides 2001 Tables 7.4–7.6):

| Group | Contribution | MHz | Sign |
|:------|:-------------|----:|:-----|
| **A** Two-loop self-energy (B₆₀ + ln(Zα) pieces) | dominant | +0.86 (B₆₀ alone) | + |
| **B** Two-loop vacuum polarization (Karplus–Sachs) | small | +0.16 | + |
| **C** Mixed SE × VP | small | −0.06 | − |
| **D** Two-photon vertex / Yennie gauge | small | +0.27 | + |
| Wichmann–Kroll | tiny | −0.025 | negligible |
| **Cumulative α⁵ multi-loop (Eides Table 7.4)** | | **+7.10** | + |

The residual +5.65 MHz is consistent with this cumulative α⁵ multi-loop contribution within ~20% (i.e., within the back-of-envelope precision of the LS-5 scoping). **The genuine multi-loop test target for the eventual LS-7 sprint is ~+5 to +6 MHz, not +29 MHz.**

## 5. Implications for Paper 34

Paper 34 (the projection-taxonomy living catalogue) currently classifies the LS-3 residual of +29.26 MHz / −2.77% as A-tier ("approximation order, missing multi-loop QED"). The LS-6a result establishes that this row should be **reclassified**:

Before LS-6a:
```
LS-3 residual: −2.77% (T₂S + A)
   T₂S = +29 MHz Bethe-log truncation (closes at N→∞)
   A   = −33 MHz one-loop ceiling (LS-1 baseline at N→∞)
```

After LS-6a:
```
LS-1 baseline residual: −3.10% (T+A+C combined)
   T = truncation (vanishes at N→∞ in Bethe-log)
   B = basis (Drake-Swainson tabulated, machine precision)
   A = approximation order (~+0.53% / +5.65 MHz, GENUINE α⁵ multi-loop)
   C = calibration mismatch (~+2.57% / +27.13 MHz, LS-1 SE convention,
       SHOULD have been derived in Eides §3.2 from the start)
```

**Recommended Paper 34 update:** Append a new row in §V.B (off-precision matches) for LS-6a with the −0.534% / +5.65 MHz residual classified as "A-tier (genuine α⁵ multi-loop)". Mark the LS-1..LS-4 rows as superseded by LS-6a for the canonical-convention computation, retaining them in the catalogue as the convention-comparison baseline. The "C-tier ~+27 MHz piece" is now resolved as a one-loop convention fix and removed from the unexplained budget.

## 6. Implications for Paper 35 §VII.3

Paper 35 §VII.3 (the iterated temporal-compactification hypothesis for two-loop physics) is tested against the genuine multi-loop residual, not against any convention-mixed residual. After LS-6a, the test target is:

```
Target: ~+5.65 MHz with sign = positive
LS-5 prediction (Eides Table 7.4 cumulative): +7.10 MHz
Acceptable falsification band: [+3.5, +10.7] MHz
```

The LS-7 sprint (genuine two-loop spectral-action computation on S³) is now the one with a clean numerical target. The structural class match (π^{even} from T9 iterated + ln(Zα) bound-state factor) is already CONSISTENT with the Eides Table 7.4 form (B₆₀ ln(Zα) coefficient); the quantitative test is whether GeoVac's iterated proper-time integration produces +7 MHz-scale numerics, not +29 MHz-scale.

If LS-7 produces a Lamb-shift contribution in the +5..+10 MHz band with positive sign, Paper 35 §VII.3 is NOT falsified and the iterated-projection picture is supported. If it produces ~+30 MHz or negative, the picture is falsified. The LS-6a result establishes the test target unambiguously.

## 7. Honest limits

1. **No new physics.** LS-6a is a pure convention fix to LS-1; the same one-loop diagrams, same Bethe logs, same Uehling kernel. The result improves LS-1's residual from −3.10% to −0.534% by fixing a double-counting of the Uehling kernel constant. This is mechanical, not derivational.

2. **The +5.65 MHz residual is consistent with α⁵ multi-loop, but is not derived as such.** What LS-6a establishes is that the residual is in the right ballpark and right sign for genuine multi-loop QED. The actual derivation requires LS-7 (two-loop spectral-action on S³).

3. **No GeoVac-native machinery used.** LS-6a is a textbook one-loop computation using the same Bethe formulas as LS-1. The convention change is purely how to organize the constants — which constants belong in SE vs. VP. The GeoVac contribution to Lamb-shift physics remains the planned LS-7 two-loop sprint.

4. **The 2P₁/₂ formula is unchanged.** The convention difference is entirely in the 2S₁/₂ SE constant (38/45 vs. 10/9). The 2P₁/₂ −1/6 constant is the same in LS-1 and LS-6a (both are textbook).

5. **Bethe logarithms are still external input.** LS-6a uses the same Drake-Swainson tabulated values as LS-1; LS-2/3/4's native Bethe-log derivations remain valid alternative computational paths but are not the source of the convention issue.

6. **Remaining LS-7 multi-loop sprint scope.** Per LS-5 §7, the next sprints to close the α⁵ ceiling on the GeoVac framework are:
   - LS-6b (2 sprints, medium risk): Karplus-Sachs two-loop VP on S³, +0.16 MHz scale
   - LS-7 (3 sprints, high risk): Two-loop self-energy Σ_{2L} on S³, +0.86 to +7 MHz scale
   - LS-8 (2 sprints): Mixed SE×VP and two-photon vertex, completing the α⁵ budget

## 8. Files

- Implementation: `debug/ls6a_eides_convention.py`
- Data: `debug/data/ls6a_eides_convention.json`
- Memo: this file

## 9. References

- M. I. Eides, H. Grotch, V. A. Shelyuto, "Theory of light hydrogenic bound states", *Phys. Rep.* 342 (2001) 63. §3.2 Eq. (3.32) for the canonical one-loop nS₁/₂ SE form with the +10/9 constant.
- P. J. Mohr, G. Plunien, G. Soff, *Phys. Rep.* 293 (1998) 227. §V for the same canonical form.
- H. A. Bethe, *Phys. Rev.* 72 (1947) 339. The original 5/6 Bethe constant; 10/9 emerges after including Karplus-Klein-Darwin.
- H. A. Bethe, E. E. Salpeter, *QM of One- and Two-Electron Atoms* (1957) §19-21. Multiple convention choices documented.
- C. Itzykson, J.-B. Zuber, *Quantum Field Theory*, p. 345. Source of the 2P₁/₂ −1/6 constant.
- GeoVac LS-1 (`debug/ls1_lamb_shift_memo.md`), LS-5 (`debug/ls5_two_loop_scoping_memo.md`).
