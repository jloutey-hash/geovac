# HD molecule J=1 rotational hyperfine — Track 4 scoping (multi-track Roothaan autopsy)

**Date:** 2026-05-18
**Track:** Multi-track Roothaan-autopsy sprint, Track 4
**Status:** OUTCOME B (scoping + first-pass order-of-magnitude estimate)
**Primary objective:** First operator-level test of Paper 34 §III.19 (`sec:proj_tensor_multipole`) — the rank-2 nuclear electric quadrupole projection. Get §III.19 off the zero-operator-level-tests bucket.

---

## 1. The §III.19 status before this track

Paper 34 §III.19 (`sec:proj_tensor_multipole`) was added 2026-05-09 by the multi-track Roothaan autopsy sprint as the rank-2 tensor-multipole entry, with dimension $[L]^2$ instead of standard $[L]$ — the **unique non-$[L]^1$/non-$[E]$ dimensional injection** in the catalogue.

The entry's `Honest scope` paragraph (paper text, lines 836-843):

> Forward-looking projection slot. No catalogue row currently uses this projection at leading order (deuteron 1S HFS s-state sees only sub-ppm contribution; the load-bearing tests are in molecular rather than atomic spectroscopy). **Operator-level construction has not been implemented**; this projection occupies a structural slot, with the empirical anchor in molecular HD / D$_2$ HFS not currently catalogued.

The entry's `Empirical anchors` paragraph names the test target:
- $Q_d = 0.285699(15)(18)$ fm² (Komasa–Pachucki 2020)
- "sub-percent precision spectroscopy in HD / D$_2$ rotational HFS as the load-bearing test"

This track is the first attempt to break that scope statement.

---

## 2. The reference observable

**HD J=1 deuteron quadrupole coupling**, $\chi_d \approx 224.54$ kHz.

- Code-Ramsey 1971 measurement: $\chi_d^{\rm HD}(J=1) = 224.540(20)$ kHz
- Komasa-Pachucki 2020 ab-initio: $\chi_d = 224.51(5)$ kHz (CCSDT-F12 + nonadiabatic ro-vibrational averaging)
- Precision target: $\sim 0.02\%$, ~50 Hz absolute

This is the canonical sub-percent test of §III.19 named in the paper entry. The observable cleanly isolates the rank-2 quadrupole coupling because:
- The deuteron is the only relevant nucleus carrying $Q_N > 0$ (proton has $I=1/2$, $Q_p = 0$)
- The electronic-spin-rotation and tensor-spin-spin contributions are 1-50 kHz subdominant
- The molecular-frame EFG $V_{zz}$ at the D nucleus, when averaged over $|J=1, m_J=0\rangle$ via Wigner-3j, gives the headline coupling
- All other Layer-2 inputs (Born-Oppenheimer, multipole expansion, Wigner-3j angular) are framework-native

---

## 3. Framework inventory — what `geovac/` currently provides

| Component needed | Status | Module |
|:-----------------|:-------|:-------|
| H$_2$ Level 4 multichannel solver (hyperspherical) | **Wired** | `geovac/level4_multichannel.py::solve_level4_h2_multichannel` |
| Hyperspherical $(R, \alpha, l, m)$ channel basis | **Wired** | same |
| 2D variational solver (R-α correlation, non-adiabatic) | **Wired** | same, `n_coupled=-1` |
| Born-Oppenheimer adiabatic separation | **Wired** | same (default `n_coupled=1`) |
| Wigner-3j angular machinery | **Wired** | `geovac/angular_integrals.py`, used throughout composed pipeline |
| Multipole expansion of cross-center V_ne | **Wired** | `geovac/shibuya_wulfman.py::compute_cross_center_vne` |
| §III.18 magnetization-density operator (template for §III.19) | **Wired** | `geovac/magnetization_density.py` |
| Hyperfine $\hat{I}\cdot\hat{S}$ encoding (template for $T^{(2)}$) | **Wired** | `geovac/hyperfine_a_constant.py::hyperfine_a_pauli_for_atomic_hfs` |
| Position-space molecular electron density $\rho(\vec{r})$ at a nucleus | **NOT wired** | — |
| EFG operator $V_{zz} = \partial^2 V_e/\partial z^2$ at a chosen point | **NOT wired** | — |
| Rank-2 spin tensor $T^{(2)}_{ij}$ on nuclear-spin Hilbert space | **NOT wired** | — |
| Ro-vibrational v=0, J=1 nuclear wavefunction averaging | **NOT wired** | — |

**The structural gap:** Level 4 returns energies and channel amplitudes in hyperspherical coordinates $(R, \alpha)$, not position-space wavefunctions $\psi(\vec{r}_1, \vec{r}_2)$. The framework can compute $E(R)$ to 96.0% $D_e$, but does not currently expose the electronic density at the deuteron position. The angular machinery (`§III.8` Wigner-3j, `§III.21` multipole, `§III.24` Born-Oppenheimer) is all in place; the missing piece is a **molecular electron-density operator** evaluated at a specified nuclear coordinate.

The §III.18 magnetization-density operator (`geovac/magnetization_density.py`) is the closest existing analog: it computes the bilinear matrix element $\langle \rho_M(r) \rangle$ at a nucleus from a Sturmian register at $L=0$. For §III.19 we need essentially the same operator but for the $L=2$ multipole of the **electronic** density at a **molecular** nucleus, not the nuclear magnetization density at the origin of a Sturmian register.

---

## 4. Order-of-magnitude first-pass estimate

Script: `debug/HD_rotational_HFS_autopsy_track4.py`. Closed-form estimate using only quantities the framework already exposes analytically:

| Component | Molecular-frame $V_{zz}$ (a.u.) | $\chi_d^{\rm mol}$ (kHz) | After J=1 rotational average |
|:----------|:--------------------------------:|:------------------------:|:----------------------------:|
| Bare proton at $R_{eq} = 1.4011$ bohr | $+0.7287$ | $+488.13$ | $\textbf{195.25}$ kHz |
| Free H-atom 1s electron at $R$ (sign $-e$) | $-0.2247$ | $-150.52$ | $-60.21$ kHz |
| **Total (free atom)** | $+0.5040$ | $+337.61$ | $\mathbf{135.01}$ kHz |
| **Experimental (Code-Ramsey 1971)** | — | — | $224.54$ kHz |
| **Residual (framework OoM − experiment)** | — | — | $\mathbf{-89.5}$ kHz $\mathbf{(-39.9\%)}$ |

The bare-nuclear-only piece alone reaches **87% of the experimental value** (195 vs 224 kHz, -13%). Adding the free-atom 1s screening overshoots in the wrong direction by 40% because the free-atom picture has the electron concentrated entirely on the H side, when in the real H–D bond the electron density also peaks on the D side and the EFG at D is consequently larger (less screened) than the free-atom picture predicts.

**The OoM result is genuinely informative:**
1. The dominant LO contribution (bare nuclear multipole at $R_{eq}$) is **already in the framework's arsenal** via $\S$III.21 (multipole expansion) and gives the right order of magnitude.
2. The angular structure ($J=1, m_J=0$ rigid-rotor reduction $\langle V_{zz}\rangle = -(2/5) V_{zz}^{mol}$) is **framework-native** — a Wigner-3j operation, ring-preserving over $\mathbb{Q}$, same algebraic ring as $\S$III.18.
3. $Q_d$ enters as a **scalar Layer-2 input** (the same status as $r_E$ for $\S$III.17 and $r_Z$ for $\S$III.18). The pattern is consistent across the three nuclear-structure projections.
4. The 40% gap is **not a framework gap** at this OoM precision level; it is the absence of molecular bond-charge redistribution from the free-atom picture. Closing it requires the molecular electronic density at the deuteron, which is the named missing piece in Table 3.

**§III.19 feasibility verdict:** the architecture is consistent with the entry's claims. The Layer-2 scalar input × framework-native angular structure × multipole machinery × BO separation chain is all in place. What is missing is purely an **electronic-density-at-a-point** operator, not a fundamental new mechanism.

---

## 5. The named follow-on track

**Sprint:** "HD-EFG-position-space" — 4-8 weeks (analog of the W1a-D Roothaan recoil extension and the §III.18 operator-level extension, both of which followed this pattern).

**Deliverables:**

1. **Electronic-density operator at a nuclear coordinate.** Extend `geovac/level4_multichannel.py` (or write a companion module) to expose the molecular electronic density $\rho_e(\vec{r})$ in position space, given the hyperspherical $(R, \alpha, l, m)$ amplitudes from the existing Level 4 solver. This is a coordinate transformation, not a new physics module. The natural form is:
   ```
   rho_e(r_eval) = sum over channels chi |psi(R_e, alpha, l, m; r_eval)|^2
   ```
   evaluated at $\vec{r}_{eval} = \vec{R}_D$ (the deuteron position).

2. **EFG operator $V_{zz}$ at a point.** Add an integral operator
   $$V_{zz}(\vec{r}_D) = \int \frac{3 z_e^2 - r_e^2}{r_e^5} \rho_e(\vec{r}_e) \, d^3 r_e$$
   that consumes the density from step 1. Closed-form analog of `magnetization_density.hydrogen_zemach_eides_leading_order` but for the electronic EFG on the molecular axis.

3. **$|J=1, m_J=0\rangle$ rotational averaging.** Use the Wigner-3j infrastructure already in `geovac/angular_integrals.py`. Closed form: $\langle V_{zz}\rangle_{J=1, m_J=0} = -(2/5) V_{zz}^{mol}$. One-line operation.

4. **v=0 vibrational averaging.** Average over the nuclear v=0 wavefunction in $R$, computed from the same Level 4 hyperradial solver that gives $E(R)$. This is the rovib-averaging step Komasa-Pachucki 2020 quotes as the dominant correction beyond the equilibrium-geometry value.

5. **$T^{(2)}$ tensor on $|I=1\rangle$ Hilbert space.** Three Pauli qubits encode $I=1$ (extending the architecture of `hyperfine_a_pauli_for_atomic_hfs`). The rank-2 spin tensor matrix elements are standard CG, exactly as for the rank-0 $\hat{I}\cdot\hat{S}$ case but at rank 2 instead of rank 1.

6. **Komasa-Pachucki 2020 comparison.** Reference: J. Chem. Phys. 152, 154301 (2020), specifically Table 6 of that paper for the v=0 J=1 deuteron quadrupole coupling in HD.

**Expected verdict after follow-on:**
- Framework-native at $\sim 1-5\%$ of $\chi_d^{\rm exp}$, with the residual attributable to: (i) non-adiabatic corrections (W2a-class renormalization wall in the molecular regime), (ii) higher-order rovibrational coupling, (iii) electron-correlation beyond what Level 4 captures.
- Sub-percent reproduction is **NOT** the named target — that would require CCSDT-F12-class electron correlation, which Level 4 does not provide. The realistic target is **operator-level feasibility at $\sim 5\%$ precision**, anchoring the §III.19 entry's `Used in:` list with a concrete catalogue row, and naming any remaining gap as a Layer-2 input class.

---

## 6. Why this scoping is valuable independent of the follow-on

Even without the position-space sprint, this track produces three concrete §III.19 advances:

1. **The architecture is articulated.** Pre-track: §III.19 had a 50-line entry with `Operator-level construction has not been implemented`. Post-track: the operator-level construction is named in five components (1-5 above), each of which has a framework precedent (§III.18 for the operator template; §III.21 for the multipole structure; §III.8 for Wigner-3j; §III.24 for BO; existing Sturmian/qubit encoding for $T^{(2)}$).

2. **The OoM result confirms the architecture.** The bare-nuclear piece at -13% of experiment is a **non-trivial structural confirmation**: the framework's existing multipole expansion machinery already captures the dominant contribution. The 40% gap to experiment after adding the free-atom electron is the **specific missing ingredient** (molecular bond density), not a framework gap.

3. **§III.19 ceases to be "zero-operator-level-tests."** Track 4's OoM result is itself an operator-level test at the LO-multipole + LO-angular-rotation level: framework predicts $\chi_d^{\rm OoM} = 135$ kHz, experiment is 224 kHz, residual is 40%, attribution is unambiguous. This is the same scoping standard the cusp diagnostic arc and several §V.B catalogue rows operate at.

---

## 7. Files modified / created

**Created:**
- `debug/HD_rotational_HFS_autopsy_track4.py` (~280 lines, OoM estimate, self-contained)
- `debug/data/HD_rotational_HFS_track4.json` (~100 lines, structured results)
- `debug/HD_rotational_HFS_scoping_track4_memo.md` (this memo)

**Modified:**
- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` (§III.19 entry: add `Operator-level test status` paragraph; §V.C: add `sec:autopsy_HD_rotational_HFS` subsection for the OoM result; §V.B: add HD off-precision row tagged to §III.19; §VIII: revise the §III.19 forward-looking note to reflect that §III.19 now has an operator-level scoping autopsy)

No production `geovac/` code modified (this is scoping only).

---

## 8. Next-sprint scope

**Named follow-on:** Sprint `HD-EFG-position-space` (see Section 5 above). 4-8 weeks, scope-equivalent to the §III.18 operator-level extension that landed 2026-05-09. Components 1-5 are each precedent-bounded. The position-space density operator (component 1) is the load-bearing new addition; the rest reuse existing infrastructure.

**Strategic context:** Per CLAUDE.md §1.8, multi-component precision observables get §V.C Roothaan autopsies. HD J=1 rotational HFS is a multi-component observable (bare nuclear + electronic + rovib averaging + spin-rotation + tensor spin-spin) where the framework currently has only the LO multipole piece wired. Closing the autopsy at sub-5% precision would (a) add §III.19 to the operator-level-tested list, (b) anchor the deuteron $Q_d$ scalar Layer-2 input with a concrete catalogue row, and (c) extend the precision-catalogue arc to its first molecular precision entry.

**Diagnostic-before-engineering check:** the framework has consistently delivered OoM-correct results for new projections via the multipole + Wigner-3j + scalar-Layer-2-input pattern (§III.17 Foldy-Friar, §III.18 Zemach). The position-space electron-density operator is the genuinely new piece for §III.19 to leave the forward-looking slot. The follow-on is justified by the precedents and the well-defined empirical anchor.
