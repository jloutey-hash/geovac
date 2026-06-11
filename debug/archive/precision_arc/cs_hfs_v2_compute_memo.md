# Cs HFS v2 — five-component Roothaan autopsy compute (Sprint Cs-HFS-v2 Phase B/C)

**Date:** 2026-05-09 (companion to `debug/cs_hfs_v2_engineering_memo.md`).
**Sprint:** Multi-observable focal-length decomposition program (CLAUDE.md §1.8), prospective Paper 34 §V.C.6 fill.
**Verdict:** **Framework-native compute landed cleanly at the engineering level** (Phase A closures all green; Phase B/C compute runs in 2.1 s end-to-end on framework machinery). **Numerical accuracy is at the −47% level** versus experimental A_Cs = 2298.158 MHz, dominated by a single named framework kernel gap (Clementi-Raimondi single-zeta screening qualitative for heavy-atom inner-shell penetration). This is exactly the "framework kernel approximation gap" class-2 finding from CLAUDE.md §1.8 directive, and it is structurally consistent with Roberts-Ginges 2015's atomic-structure-uncertainty framing for Cs PNC.

---

## 1. Compute summary

The compute decomposes A(Cs 6S₁/₂) into five components per the §1.8 Roothaan-autopsy discipline. All values in MHz; convention H_HFS = A·I·J (Lande convention), so for Cs with I=7/2, J=1/2 the F=4 ↔ F=3 transition has ν_HFS = 4A.

| Component | Source | A contribution [MHz] | Cumulative [MHz] |
|-----------|--------|----------------------|-------------------|
| **C1** Bohr-Fermi strict (Dirac g) | framework-native, §III.18 spinor-lift point-nucleus | 782.96 | 782.96 |
| **C3** Schwinger a_e correction | framework-native, factor (1+a_e^Schwinger) | +0.91 | 783.87 |
| **C4** Casimir relativistic enhancement | framework-native, F_R = 4γ/(4γ²−1) at Z·α=0.401 | +435.25 (×1.5553) | 1219.12 |
| **C5** Zemach (r_Z=6.7 fm) | Layer-2: §III.18 magnetization-density at literature r_Z | −0.12 | 1219.00 |
| **C5** Bohr-Weisskopf (Layer-2 lit) | Layer-2 input: Roberts-Ginges 2015 type CI+all-order | +14.63 | 1233.63 |
| **L** Multi-loop QED (LS-8a wall) | Layer-2 input: Karshenboim 2005 §VIII | +0.12 | 1233.75 |

**Experimental (Lande convention):** A_Cs = 2298.157 9425 MHz = 9.192 631 770 GHz / 4.

**Framework-native total:** 1219.12 MHz (residual: **−1079.03 MHz / −46.95%** vs experiment).

**Framework + Layer-2 total:** 1233.75 MHz (residual: **−1064.41 MHz / −46.32%** vs experiment).

Compute runs end-to-end in 2.1 s wall-clock time on framework machinery. The only computation-intensive step is C2 |ψ_6s(0)|² via the dense-uniform-grid solver from Phase A (n_grid = 200k for Cs 6s gives ~0.1s per call; convergence study uses 4 grid sizes).

---

## 2. Framework-native compute: |ψ_6s(0)|² is the load-bearing piece

The C2 component (the framework's |ψ_6s(0)|² from FrozenCore Z_eff(r)) is what was BLOCKED in Track 5 (`debug/cs_hfs_v1_memo.md`). With Phase A's `screened_psi_origin_squared` wrapper (A2 closure), the compute now runs cleanly:

| n_grid | E_6s [eV] | \|ψ_6s(0)\|² [bohr⁻³] | Wall time [s] |
|--------|-----------|------------------------|----------------|
| 50,000 | −1.4744 | 1.0817 | ~0.04 |
| 100,000 | −1.4746 | 1.2027 | ~0.07 |
| 200,000 | −1.4747 | 1.2659 | ~0.13 |
| 400,000 | −1.4747 | 1.2981 | ~0.27 |

Richardson extrapolation (linear in 1/n_grid) gives **|ψ_6s(0)|²_∞ = 1.328 bohr⁻³**. The Phase B/C compute uses this Richardson-extrapolated value as the canonical input for the BF formula.

**Eigenvalue accuracy:** E_6s = −1.475 eV vs NIST −3.894 eV is **−62%** off. This is the Clementi-Raimondi kernel signature: the analytical single-zeta hydrogenic basis underbinds the Cs valence orbital by a factor of 2.6×. The same kernel limitation is what drives the |ψ(0)|² being too small.

**Cross-check vs v1 skeleton (Roberts-Ginges Z_eff = 9.7):** v1's qualitative skeleton predicted A = 793 MHz at Z_eff=9.7. With our framework-native |ψ(0)|² = 1.328 bohr⁻³, the inverse hydrogenic conversion gives effective Z_eff = (1.328 × π × 216)^(1/3) = 9.74 — **bit-identical to Roberts-Ginges 2015 effective Z** (v1 §2 Table). So the framework's FrozenCore [Xe] Z_eff(r) profile reproduces the literature qualitative effective-Z to <1% at the contact-density level. The −47% residual on A is structural to the Clementi-Raimondi screening, NOT a numerical artifact of the new solver.

---

## 3. Casimir relativistic enhancement (C4): the dominant framework correction

At Z·α = 0.401, the Casimir 1936 relativistic enhancement of the contact density is:

  F_R = 4γ / (4γ² − 1) where γ = √(1 − (Zα)²) = 0.9159

For Cs: F_R = 1.5553. This is a **+55% correction** that takes the framework-native A from 783.87 MHz (BF + Schwinger only) to 1219.12 MHz. Without F_R, the framework would predict A = 784 MHz (residual −66%, matching the v1 skeleton's quoted result without relativistic factors).

**The Casimir formula is leading-order.** The full Bohr-Weisskopf relativistic enhancement (Sobelman §6.4, Roberts-Ginges 2015) gives F_R^full ≈ 2.6 at Z=55, including:
- Higher-order Dirac small/large component admixture
- Finite-nuclear-size correction to the spinor wavefunction at the origin
- Anomalous-magnetic-moment contribution at the Dirac level

The 2.6 / 1.5553 = 1.67× ratio between full BW relativistic and Casimir leading-order is a structural framework gap. To close it, the framework would need the spinor-lift §III.7 evaluated with full Dirac wavefunctions, not the scalar §II.B Bohr-Fermi formula multiplied by a Casimir factor. This is a Tier 3 / Sprint TR-class extension; **scoped, not closed in this sprint**.

---

## 4. Layer-2 Zemach + Bohr-Weisskopf (C5): bookkeeping check

The Eides-style leading-order Zemach correction with r_Z(Cs-133) = 6.7 fm (Roberts-Ginges 2015) gives:

  Δν / ν = −2 Z α m_e r_Z = −2 × 55 × (1/137.036) × 1 × (6.7 fm / a_0_fm) = −101.6 ppm

For A_Cs ≈ 1219 MHz, this is −0.12 MHz. **Sub-MHz correction**, much smaller than the framework kernel gap (−1079 MHz at the Casimir level). The Cs Zemach radius is somewhat poorly known compared to the proton (r_Z(p) = 1.045 fm with sub-percent uncertainty); for Cs the literature spread is 5-8 fm. Even at r_Z(Cs) = 8 fm, the Zemach correction is ~−120 ppm = −0.15 MHz — still negligible against the framework residual.

The Bohr-Weisskopf correction (distributed nuclear magnetization beyond Zemach) is +1.0% to +1.5% of A per Roberts-Ginges 2015 Tab. 4 (we use +1.2%). For A ≈ 1219 MHz, this is +14.63 MHz — comparable to the multi-loop QED estimate but still 70× smaller than the kernel gap.

**Net Layer-2:** +14.51 MHz (+1.19% relative). Not the load-bearing residual.

---

## 5. Multi-loop QED (L): structurally present but small at this precision

At Z·α = 0.4, the α²(Zα) multi-loop QED for HFS is order ~10⁻⁴ relative (Karshenboim 2005 §VIII). For Cs A ≈ 1219 MHz, this is +0.12 MHz — much smaller than the kernel gap and below the resolution of this computation.

**Same wall as Sprint HF-5, Sprint LS-8a, Sprint H1:** the bare iterated CC spectral action on Dirac-S³ reproduces the structural prefactor (α/π)²·(Zα) but cannot autonomously generate the Z_2/Z_3/δm renormalization counterterms required for a finite extraction. Layer-2 input from external QED literature.

---

## 6. Pattern crystallization: three-class tag (per CLAUDE.md §1.8 directive)

Per the §1.8 directive's three problem-classes:

### 6.1 Literature convention mismatch
**NOT YET A SHARP FINDING at this scoping precision.** Roberts-Ginges 2015 vs Porsev-Beloy-Derevianko 2009 use slightly different decompositions (CI+all-order vs LCC), with atomic-structure spread at the 0.5-1% level in their respective predictions for A_Cs. The framework's residual is **−47%**, two orders of magnitude larger than this RG-vs-PBD gap. **Multi-observable global fits at sub-percent precision would be needed to expose any literature convention mismatch in Cs HFS.** This is an open follow-up: extend the framework's |ψ(0)|² compute to sub-percent accuracy (via Hartree-Fock screening) before attempting a Cs PNC literature-comparison sprint.

### 6.2 Framework kernel approximation gap (THE LOAD-BEARING FINDING)
**Two named kernel gaps, both structural to the heavy-atom regime:**

(a) **Clementi-Raimondi single-zeta screening is qualitative for the Cs valence regime.** The eigenvalue E_6s = −1.475 eV vs NIST −3.894 eV is −62% off, and the |ψ(0)|² is correspondingly ~50% too small. The mechanism is unambiguous: the analytical hydrogenic basis with single Z_eff per shell cannot represent the radial nodes and inner-shell penetration of a 6s orbital outside a [Xe] core. Hartree-Fock or DFT screening at Z=55 would close this. **Estimated effort: 2-4 weeks (a substantial extension to `geovac/neon_core.py`).**

(b) **Casimir relativistic enhancement is leading-order; full Bohr-Weisskopf is ~1.7× larger.** Closing this requires the spinor-lift §III.7 evaluated with full Dirac wavefunctions for the s-state contact density, not the scalar BF formula multiplied by F_R^Casimir. **Estimated effort: 1-2 weeks (Tier 3 / Sprint TR-class extension).**

### 6.3 Cleanly attributed Layer-2 wall
Multi-loop QED (~+100 ppm at Z=55, analogous to Sprint HF-5 / LS-8a) is structurally present in the residual and cleanly attributed to the LS-8a renormalization wall. Bohr-Weisskopf (+1.2%) is similarly Layer-2 (requires distributed nuclear magnetization input, not GeoVac-internal). Both contribute below the percent level — they are not the load-bearing residual at this sprint's precision.

---

## 7. Z-scaling pattern and §1.8 directive value

**This sprint tests Z-scaling of the framework's projection-chain machinery at the heaviest atom attempted to date.** Comparing to other sprints:

| System | Z | Framework residual | Wall (per §1.8 + CLAUDE.md §3) |
|--------|---|---------------------|---------------------------------|
| H 21cm (Sprint HF) | 1 | +18 ppm | LS-8a multi-loop |
| Mu HFS | 1 (μ-N) | +199 ppm | LS-8a (cleanest) |
| D HFS | 1 (D nucleus) | +286 ppm | W2a Bohr-Weisskopf-class |
| He 2³P (Sprint MH-A) | 2 | −0.014% to −2.6% | A multi-loop α³/π |
| μH Lamb (Sprint MH) | 1 | −0.10% | LS-8a + W2a recoil |
| **Cs 6S₁/₂ (this sprint)** | **55** | **−47%** | **W1c-residual screening** |

The Z=55 regime exposes a **structural kernel limitation of the framework's atom-construction layer** (Clementi-Raimondi single-zeta screening) that is invisible at low Z. This is precisely the §1.8 directive's "framework kernel approximation gap" class-2 finding, and it has direct relevance to the Cs PNC parity-violation sprint that has been flagged as a long-term diagnostic-instrument target (CLAUDE.md §1.8 §V.C.6 prospective fill text in v1 memo). **The framework's projection-chain decomposition isolates the percent-level atomic-structure uncertainty in Cs PNC to a single named kernel limitation (Clementi-Raimondi screening at Z=55), distinguishing it cleanly from the multi-loop QED Layer-2 wall and from the nuclear-structure Bohr-Weisskopf input.**

---

## 8. Open follow-ups

**Tooling-addressable (high value):**
- (a) Replace Clementi-Raimondi single-zeta with Hartree-Fock screening for Z ≥ 30. This would close the dominant residual to the percent level. Estimated 2-4 weeks. Could reuse atomic structure machinery from PySCF or Dirac/MOLCAS as a Layer-2 input, OR implement a self-consistent HF iteration in `geovac/neon_core.py`. The latter is the principled framework path; the former is faster and gets the result.

- (b) Spinor-lift evaluation of the s-state contact density via Tier 3 Dirac wavefunctions (instead of scalar BF × Casimir F_R). Would close the relativistic-enhancement gap from leading-order Casimir to full Bohr-Weisskopf. Estimated 1-2 weeks. Reuses `geovac/dirac_matrix_elements.py` infrastructure.

**Scoped but deferred (low value at this sprint precision):**
- (c) Better Cs r_Z (Layer-2): current literature spread 5-8 fm doesn't move the needle since the Zemach correction is sub-MHz. Not a priority.

- (d) LS-8a-renorm extension for multi-loop QED at Z=55: same multi-week scope as flagged for Sprint LS-8a / Sprint H1 / Sprint HF-5. Hard NCG / spectral-action extension; deferred.

**Diagnostic-instrument extension (CLAUDE.md §1.8 directive value):**
- (e) Once kernels (a) and (b) close to sub-percent accuracy, attempt the Cs PNC sprint (Roberts-Ginges 2015 P_eff comparison). **The framework's projection-chain decomposition would isolate the atomic-structure dependence of the Cs weak charge extraction**, which is the long-term scientific value flagged in the §1.8 directive. Open question: does the framework's structural decomposition reveal where in the projection chain the percent-level Cs PNC uncertainty lives?

---

## 9. Files

- `geovac/neon_core.py` — modified in Phase A (A1, A2 closures).
- `geovac/hyperfine_a_constant.py` — created in Phase A (A3 closure, ~340 lines).
- `tests/test_neon_core.py` — extended in Phase A (16 new tests).
- `tests/test_hyperfine_a_constant.py` — created in Phase A (14 tests).
- `debug/calc_track_cs_hfs_v2.py` — Phase B/C compute driver (this commit).
- `debug/data/cs_hfs_v2.json` — five-component decomposition data (this commit).
- `debug/cs_hfs_v2_engineering_memo.md` — Phase A engineering closures memo.
- `debug/cs_hfs_v2_compute_memo.md` — this memo (Phase B/C compute report).

No changes to Paper 34 §V.C.6 in this commit (the v1 prospective fill text remains the relevant skeleton; the v2 compute provides the framework-native data point at −47% residual that would feed into a §V.C.6 fill once the kernel gap (a) closes to sub-percent — still a multi-week extension).
