# Precision Catalogue: Helium 2³P Fine Structure

**Sprint:** Precision catalogue extension to multi-electron fine structure (post-PK-cross-center round, May 2026).
**Status:** **POSITIVE** — bears partial fruit. Dominant intervals (P0-P1, P0-P2) sub-percent; small P1-P2 interval at -2.6% (partial-cancellation amplification of an absolute residual that sits inside the LS-8a one-loop QED budget).

---

## 1. Headline numbers

| Interval | GeoVac (MHz) | NIST (MHz) | Pachucki 2006 (MHz) | Abs err (MHz) | Rel err vs NIST |
|:---|---:|---:|---:|---:|---:|
| P₀ - P₁ | +29,612.91 | +29,616.951(6) | +29,616.952 | -4.04 | **-0.014%** |
| P₁ - P₂ | +2,231.16 | +2,291.176(15) | +2,291.180 | -60.02 | -2.62% |
| P₀ - P₂ | +31,844.07 | +31,908.131(6) | +31,908.132 | -64.06 | **-0.201%** |

The ABSOLUTE residuals on each interval are similar (~4–60 MHz). The fractional residuals differ because P₁-P₂ is the small partial-cancellation difference between two ~30 GHz numbers, amplifying the relative residual by ~10× over P₀-P₂.

---

## 2. Architecture

He (1s)(2p) ³P_J is the **first multi-electron precision test in the catalogue**. The two electrons sit on the same nucleus but at different effective Z's: 1s sees Z_eff ≈ 2 (full nuclear), 2p sees Z_eff ≈ 1 (full-shield). This is the *internal multi-focal* cell of the catalogue — distinct from cross-register multi-focal (Sprint MH, Sprint HF, Mu/D HFS) which couples electronic and nuclear/leptonic registers at vastly different focal lengths.

The Breit-Pauli operator decomposition follows Drake 1971:

$$E({}^3P_J) = \tfrac{\zeta_{2p}}{2} X(J) + A_{SS}\, f_{SS}(J) + A_{SOO}\, f_{SOO}(J)$$

with:
- **J-pattern** (Sprint 4 DD, sympy-exact from rank-k 6j algebra):
  - $f_{SS}(J) = (-2, +1, -\tfrac{1}{5})$ for $J = 0, 1, 2$ from $(-1)^{L+S+J} \cdot 6j\{L,S,J;S,L,2\}$ at $L=S=1$
  - $f_{SOO}(J) = (+2, +1, -1)$ from $(-1)^{L+S+J} \cdot 6j\{L,S,J;S,L,1\}$
- **Spin-tensor selection** (Sprint 4 DD, sympy-exact from 9j):
  - $\langle S{=}1\|[s_1\otimes s_2]^{(2)}\|S{=}1\rangle = \sqrt{5}/2$ (non-zero, SS uses rank-2 spin tensor)
  - $\langle S{=}1\|[s_1\otimes s_2]^{(1)}\|S{=}1\rangle = 0$ (SOO uses sum $\vec s_1 + 2\vec s_2$, Bethe-Salpeter §38.15)
- **Drake combining coefficients** (Sprint 3 BF-D, intrinsic-tier rationals):
  - $A_{SS} = \alpha^2(3/50\, M^2_{\text{dir}} - 2/5\, M^2_{\text{exch}})$
  - $A_{SOO} = \alpha^2(3/2\, M^1_{\text{dir}} - M^1_{\text{exch}})$
- **Bipolar harmonic structure** (Sprint 5 DV characterization):
  - Direct path: $(k_1=0, k_2=2)$ via $\langle 0\|C^{(0)}\|0\rangle\langle 1\|C^{(2)}\|1\rangle$
  - Exchange path: $(k_1=1, k_2=1)$ via $\langle 0\|C^{(1)}\|1\rangle$
  - Higher $(k_1, k_2)$ channels Gaunt-forbidden
- **Convention**: $\zeta_{2p} = \alpha^2 Z_{\text{nuc}} Z_{\text{eff}}^3 / (n^3 l(l+\tfrac{1}{2})(l+1))$ with $Z_{\text{nuc}}=2$, $Z_{\text{eff}}=1$. Sprint 5 CP §2.1 documents this is a He-specific coincidence (Z_nuc/Z_val = 2 = factor of 2 cancels with explicit /2 in $E_{SO}$). For Li/Be the convention requires Z_val=1 + Slater Z_eff (Sprint 5 CP).

---

## 3. Component decomposition

| Interval | SO (MHz) | SS (MHz) | SOO (MHz) | Total (MHz) | NIST (MHz) |
|:---|---:|---:|---:|---:|---:|
| P₀-P₁ | -29,198 | +23,748 | +35,063 | **+29,613** | +29,617 |
| P₁-P₂ | -58,396 | -9,499 | +70,126 | **+2,231** | +2,291 |
| P₀-P₂ | -87,594 | +14,249 | +105,189 | **+31,844** | +31,908 |

Note the structural cancellation in P₁-P₂: SO (-58 GHz) + SS (-9.5 GHz) + SOO (+70 GHz) = +2.2 GHz, a factor 30+ smaller than any single contribution. The framework's small absolute residual (60 MHz) on a small physical interval (2.3 GHz) yields the amplified -2.6% fractional residual.

---

## 4. Residual attribution (LS-8a wall)

The framework reproduces the leading-order $\alpha^2(Z\alpha)^2$ Breit-Pauli operator. The residuals attribute to LS-8a-class corrections that are inherited from underlying field theory and require renormalization counterterms the framework does not autonomously generate (Sprint H1, LS-8a, May 2026):

- **$\alpha^3/\pi$ one-loop QED on fine structure**: ~0.23% budget on dominant intervals (Drake 1971 §IV; Pachucki 2006 PRA 74, 022512)
- **$\alpha^3$ recoil corrections** ($m_e/m_{4\text{He}}$): ~0.04% budget
- **$\alpha^2(Z\alpha)^4$ second-order Breit-Pauli**: ~0.005% budget
- **$\alpha^4$ multi-loop QED**: LS-8a wall regime, sub-MHz on these intervals

P₀-P₁ residual -4 MHz / -0.014% sits below the one-loop QED budget — surprisingly clean for the dominant interval. P₀-P₂ residual -64 MHz / -0.20% sits within the one-loop QED budget. P₁-P₂ residual -60 MHz absolute is within the same budget but amplified to -2.6% by partial cancellation.

The residual pattern is consistent with Drake 1971's analysis: the dominant J-splittings are well-reproduced by leading Breit-Pauli, while the smaller P₁-P₂ interval is sensitive to higher-order corrections that only enter at α³ and beyond.

---

## 5. Multi-electron internal multi-focal architecture: VERIFIED

This catalogue entry is the first verification of the framework's architecture for **internal multi-focal** systems — two electrons on the same nucleus at distinct effective Z's, coupled via two-body operators with bipolar harmonic structure.

What's verified:
1. **Distinct effective Z's coexist**: 1s at full nuclear, 2p at full-shield, both sympy-exact with their respective ⟨1/rᵏ⟩ Pochhammer rationals.
2. **Bipolar harmonic decomposition is angular-content-only**: only two channels contribute (k₁=0,k₂=2 direct and k₁=1,k₂=1 exchange), independent of focal-length ratio.
3. **Drake J-pattern** is sympy-exact and intrinsic (Paper 18 intrinsic tier).
4. **Drake combining coefficients** $(\tfrac{3}{50}, -\tfrac{2}{5}, \tfrac{3}{2}, -1)$ are pure rationals with no convention-dependent transcendental content — Paper 18 intrinsic tier.

What's *not* introduced:
- No new transcendental class beyond α² (one-loop α³/π corrections sit in known LS-8a budget).
- No new projection mechanism beyond Wigner 3j/6j/9j (already cataloged in Paper 34 §III).

This consolidates the multi-focal architecture's coverage matrix:

| Multi-focal regime | Test cases | Residual range |
|:---|:---|:---|
| Cross-register (electron + leptonic nucleus) | Mu 1S-2S, Mu 2S-2P, Mu HFS, μH | ppm to <1% |
| Cross-register (electron + I=1/2 nucleus) | H 21cm | +18 ppm |
| Cross-register (electron + I=1 nucleus) | D 1S HFS | +40 ppm |
| Cross-register (e⁺e⁻, equal-mass) | Ps 1S HFS | +0.49% (with Layer-2) |
| **Internal (two electrons, same nucleus)** | **He 2³P** | **-0.014% / -0.20% / -2.6%** |

---

## 6. Cross-references

- **Sprint 3 BF-D** (April 2026): original Drake-pattern He 2³P at -0.20% on span. `tests/test_breit_integrals.py::test_drake_combining_coefficients_reproduce_nist`.
- **Sprint 4 DD** (April 2026): J-pattern derived sympy-exact from rank-k 6j algebra. `debug/dd_drake_derivation.py`, `debug/dd_drake_derivation.md`.
- **Sprint 5 CP** (April 2026): convention work for Li 2²P and Be 2s2p ³P; documents He as the special case where Z_nuc=2 coincides with Z_val=1 by structural numerology. `debug/cp_fine_structure_memo.md`.
- **Sprint 5 DV** (April 2026): bipolar harmonic structure for (1s)(2p) configuration. `debug/dv_drake_bipolar_memo.md`.
- **Track 1 precision catalogue** (Mu 1S-2S, Ps HFS) and **Sprint MH** (μH Lamb shift, μH HFS) — cross-register multi-focal context.
- **Sprint H1 + LS-8a** (May 2026): LS-8a wall framing for multi-loop QED corrections that the framework does not autonomously generate.

---

## 7. Paper updates applied

- **Paper 34 §V** (machine-precision): two new rows added under the multi-projection chain section for the dominant intervals P₀-P₁ (-0.014%) and P₀-P₂ (-0.20%).
- **Paper 34 §V.B** (off-precision): one new row for the partial-cancellation P₁-P₂ interval (-2.62%, error code A — approximation order, leading-order Breit-Pauli misses the higher-order corrections that survive in the small-difference regime).
- **Paper 14 §V**: cross-reference paragraph noting He 2³P as the **first internal multi-focal precision catalogue entry** and connecting Sprint 3 BF-D / Sprint 4 DD / Sprint 5 CP/DV to the post-MH catalogue arc. The existing Sprint 5 CP closure language for He/Li/Be is preserved.
- **CLAUDE.md §2**: one-paragraph entry summarizing the catalogue extension to multi-electron systems.

---

## 8. Verdict

The framework's **internal multi-focal architecture is verified at sub-percent on dominant fine-structure intervals**. Sprint 3 BF-D's Drake-pattern construction (which has been in the codebase since April 2026 with sub-percent on the span) is now formally entered into Paper 34's living matches catalogue and connected to the broader multi-focal architecture arc.

The catalogue now covers:
- **Mass-hierarchy axis**: H, μH, Mu, Ps (4 mass-ratio regimes)
- **Nuclear-spin axis**: I=1/2 (proton, antimuon), I=1 (deuteron)
- **Multi-focal kind axis**: cross-register electron-nucleus, cross-register equal-mass leptonic, **internal multi-electron** (this entry)
- Seven precision systems sub-percent on framework-native parts

---

## 9. Files

- `debug/precision_catalogue_he_2_3p.py` — main compute script
- `debug/precision_catalogue_he_2_3p_memo.md` — this memo
- `debug/data/precision_catalogue_he_2_3p.json` — full numerical results
- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` — 2 §V rows + 1 §V.B row
- `papers/group4_quantum_computing/paper_14_qubit_encoding.tex` — cross-reference paragraph
- `CLAUDE.md` — §2 sprint outcome paragraph

No production `geovac/` modifications. No git commits (worker fork).
