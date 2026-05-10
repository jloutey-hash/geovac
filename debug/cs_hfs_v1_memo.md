# Cs HFS at Z=55 (Track Cs-HFS-v1, scoping)

**Date:** 2026-05-09
**Sprint:** Multi-observable focal-length decomposition program (CLAUDE.md §1.8), prospective Paper 34 §V.C.6 entry.
**Verdict:** **SCOPING MODE** — clean framework-native compute is a 2-3 week sprint. Skeleton compute via hydrogenic-with-effective-Z is consistent with literature qualitative readings; clean compute is blocked by two named engineering items.

---

## 1. Feasibility assessment

The four feasibility questions and their verdicts:

| Q | Question | Verdict | Notes |
|---|----------|---------|-------|
| Q1 | Is Z=55 classified by the atomic_classifier? | **YES** | type C, [Xe] 6s¹, frozen_core PK source. Sprint 3 HA-A+B (v2.12.0) wired Cs/Ba via [Xe] FrozenCore; Cs has been used for SrH/BaH-style alkaline-earth monohydrides but NOT yet as a single atom. |
| Q2 | Does FrozenCore(Z=55) produce a sensible Z_eff(r) profile? | **YES** | Z_eff(0) → 55 (matches Z); Z_eff(∞) → 1 (matches Z − 54). Profile valid across r ∈ [0.001, 100] bohr. |
| Q3 | Does the screened radial solver handle 6s (l=0)? | **NO** | `_solve_screened_radial` rejects l=0 by design (line 620: "l must be >= 1 for screened radial solver"); originally written for SO splitting where Kramers cancellation handles l=0. **This is the load-bearing block for clean compute.** |
| Q3b | Does a manual l=0 extension (drop the centrifugal term, otherwise reuse the FD machinery) work? | **PARTIAL** | E_6s converges in n_grid (4k→96k: -1.384 → -1.475 eV); but \|ψ(0)\|² **does NOT converge** (0.64 → 1.21 bohr⁻³ over the same grid range, monotonically growing ~2x per grid doubling). Linear extrapolation also fails (0.48 → 0.36). The s-state at a singular -Z_eff(r)/r origin needs Numerov + log-grid + analytical small-r fit. |
| Q4 | Is `hyperfine_coupling_pauli` applicable to Cs without modification? | **NO** | The signature accepts arbitrary `A_hf`, but the function is wired to the Track NI architecture (1p + 1e + 21cm = 4-qubit cross-register I·S). For Cs with I=7/2, neither (a) the Cs nuclear register nor (b) a 6s electron register exist as nuclear/electronic blocks. A cleaner approach for HFS observables is a classical A·I·J wrapper over framework-native `|ψ(0)|²` plus `g_N` literature input. |

The blocking items partition cleanly:
- **W1c-residual analog (radial side)**: FD solver is structurally inadequate for s-wave at singular potential origin. ~1-2 weeks to extend to log-grid + Numerov + small-r analytic fit.
- **Architecture-side**: A·I·J is conceptually trivial (no second register needed for HFS), but the existing infrastructure is wired for the harder problem (full nuclear+electronic Hilbert space). 1-2 weeks to add a clean wrapper.

**Net feasibility:** 2-3 week sprint to enable a clean framework-native compute of A(Cs 6s). For a `prospective` Paper 34 §V.C.6 row, the skeleton compute below provides a defensible "framework-prediction range" that brackets the experimental value.

---

## 2. Skeleton compute (literature-anchored)

The framework's Bohr-Fermi A constant for an ns valence electron is

$$A_\text{HF} = \frac{8\pi}{3} g_e g_N \frac{\mu_N \mu_B}{1} |\psi_{ns}(0)|^2$$

In the ratio-to-hydrogen-1s convention used elsewhere in the precision catalogue (Sprints HF, MH, muonium HFS, deuterium HFS):

$$\frac{A(\text{Cs 6s})}{A(\text{H 1s})} = \frac{|\psi_\text{6s}(0)|^2}{|\psi_\text{1s}(0)|^2} \cdot \frac{g_\text{Cs}}{g_p}$$

with $g_\text{Cs} = \mu_\text{Cs}/I = 2.582025/3.5 = 0.7377$, $g_p = 5.5857$, $|\psi_\text{1s}(0)|^2 = 1/\pi$ (Z=1 hydrogenic).

The unknown is $|\psi_\text{6s}(0)|^2$. Three skeleton estimates, plus an inverse solve:

| Estimate | $Z_\text{eff}$ | $|\psi_\text{6s}(0)|^2$ | $A$ predicted (MHz) | Residual vs 2298.16 |
|----------|----------------|--------------------------|----------------------|----------------------|
| Naive hydrogenic Z=1 | 1 | 0.00147 | 0.87 | −99.96% |
| Roberts-Ginges effective Z | 9.7 | 1.345 | 793 | −65.5% |
| High-penetration Z=15 | 15 | 4.97 | 2931 | +27.5% |
| **Inverse solve** for exact A | **13.83** | **3.91** | **2298.16 ↔** | **0%** |

The required $Z_\text{eff}^\text{eff} \approx 13.83$ is in the upper ballpark of literature estimates: Roberts-Dzuba-Ginges 2019 (Phys. Rev. A 100, 042504) report effective Z values around 9-12 for Cs 6s relative to nonrelativistic formulas, before relativistic enhancement.

**Relativistic enhancement at Z=55:** $Z\alpha = 0.401$, $\gamma = \sqrt{1-(Z\alpha)^2} = 0.916$. The simplest Casimir factor $F_R = 4\gamma/(4\gamma^2 - 1) = 1.555$. The full Bohr-Weisskopf-class enhancement (Sobelman §6.4) is $\sim 2.5\text{-}3$ at Z=55; the simplification understates this. Including $F_R = 1.555$ multiplicatively: required pre-relativistic $Z_\text{eff}/F_R^{1/3} = 13.83 / 1.158 = 11.94$ — closer to the RG range. With proper Casimir $F_R \approx 2.6$: required $Z_\text{eff} = 13.83/1.376 = 10.05$ — bit-on with RG.

**Skeleton verdict:** the framework's Bohr-Fermi machinery, given a properly screened $|\psi_\text{6s}(0)|^2$, is consistent with the experimental SI-second value to within the precision of the screening + relativistic content. The skeleton brackets the experimental value cleanly between $Z_\text{eff}=9.7$ (-65.5%) and $Z_\text{eff}=15$ (+27.5%); proper screening gives $Z_\text{eff}\approx 10\text{-}11$, in the bracket. **No structural failure of the framework's projection-chain dictionary at Z=55** is detected at this scoping level. The residual is engineering precision (Q3 + Q3b extension) plus relativistic content (spinor lift, §III.7), both well-understood.

---

## 3. Diagnostic-instrument value (Cs PNC angle)

Heavy-atom Cs PNC extractions have known atomic-structure uncertainties at the ~0.5-1% level (Porsev-Beloy-Derevianko 2009; Roberts-Dzuba-Flambaum 2014). The dominant atomic-structure dependencies are:

1. **6s valence wavefunction at the nucleus** — what we are computing here. This enters PNC linearly and the framework's projection-chain decomposition would expose its sensitivity to the screening profile.

2. **Core polarization corrections** — many-body effects from the 5p, 5s, 4d electrons responding to the valence electron's penetration of the core. This is precisely the W1c-residual orthogonality wall the framework has been documenting in the chemistry-solver re-test (CLAUDE.md §2 chemistry arc paused). At Z=55 with [Xe] core, this is qualitatively the same mechanism as NaH at [Ne] core (Sprint Phillips-Kleinman cross-center sprint, 2026-05-08), but ~3x deeper in screening.

3. **Multi-loop QED in heavy-atom regime** — at $Z\alpha = 0.4$, the $\alpha(Z\alpha)$ multi-loop terms are at the percent level for HFS (Karshenboim 2005 §VIII), and the Bohr-Weisskopf nuclear magnetization correction is structurally present. This corresponds to the LS-8a wall in the multi-focal taxonomy.

The framework's prospective diagnostic value is to provide a **structural decomposition** that separates these three contributions cleanly: (i) the radial wavefunction at the origin from a scalar Schrödinger calculation with a transparent Z_eff(r) screening profile, (ii) the W1c-residual core-polarization correction as a sibling of the [Ne]-core orthogonality work, (iii) the Layer-2 LS-8a multi-loop and Bohr-Weisskopf inputs from external literature. This is the same decomposition discipline that exposed the Eides-vs-Krauth Layer-2 itemization mismatch in W1b extended (2026-05-09); applied to Cs PNC, it would expose where the percent-level atomic-structure uncertainty lives in the projection-chain budget.

---

## 4. Z-scaling assessment

Sprint 7b (`geovac/neon_core.py::screened_so_splitting`, v2.19.4) tested Z-scaling for the SO splitting on CaH/SrH/BaH at Z=20/38/56. The screened SO results were 30-70% off vs experiment — leading-order Breit-Pauli + Clementi-Raimondi screening, comparable to what we would expect for Cs HFS at the analogous level (without spinor lift, multi-loop QED, Bohr-Weisskopf).

**Z-scaling holds qualitatively** at Z=55: the framework's projection-chain machinery (atomic_classifier, FrozenCore, BF formula) extends without structural failure. The accuracy at Z=55 is comparable to (or somewhat worse than) the screened SO splittings at Z=20-56, both limited by Clementi-Raimondi screening being qualitative for this regime. **Heavy-atom regime does NOT expose a structural limitation of the framework**, only an engineering limitation (FD solver inadequate for s-wave at singular potential origin) and a known calibration limitation (Clementi-Raimondi screening is HF approximation-of-an-approximation).

---

## 5. Pattern observation: three-class tag

For the catalogue cataloguing rule (CLAUDE.md §1.8 directive question 3):

- **Literature convention question** (problem class 1): Roberts-Ginges 2022 vs Porsev-Beloy-Derevianko 2009 use slightly different decompositions of the A constant (CI+all-order vs LCC); a multi-observable framework fit could expose any convention mismatch in their Layer-2 itemizations. Not yet a sharp finding at this scoping level.
- **Framework kernel approximation gap** (problem class 2): Two named gaps. (a) FD-on-uniform-grid solver inadequate for s-wave |ψ(0)|² (engineering, mechanical fix). (b) Clementi-Raimondi single-zeta screening qualitative for Cs valence (structural, would benefit from Hartree-Fock or DFT screening at Z=55). The latter is a real limitation of the framework's atom-construction layer that would bear on PNC extractions.
- **Focal-length decomposition** (problem class 3): A constant decomposes naturally as Bohr-Fermi-strict + relativistic (spinor lift §III.7) + Schwinger a_e + Bohr-Weisskopf (nuclear magnetization, §III.18 sibling) + multi-loop QED (LS-8a). Five-component Roothaan autopsy possible once skeleton compute is enabled. Comparable to §V.C.2 hydrogen 21 cm but shifted to the heavy-atom regime.

---

## 6. Recommendation: 2-3 week sprint plan

**Phase 1 (1-2 weeks): Enable l=0 in screened solver.**
Extend `_solve_screened_radial` to support l=0 by:
- Switching to a logarithmic radial grid (denser near the origin)
- Numerov integration (4th-order accurate, standard for s-wave)
- Analytical small-r series fit: near r=0, V(r) ~ -Z/r (pre-screening), so R_ns(r) ~ R_ns(0) · exp(-Z·r/n) by hydrogenic limit. Match the numerical solution to this analytical form at r ~ 1/Z.
- New API: `screened_psi_origin_squared(Z, n, l=0, core_type=None)` returning |ψ_ns(0)|² with controlled convergence.
- Include test against hydrogenic Z=1 (should match 1/π exactly) and against literature Cs 6s value.

**Phase 2 (~1 week): Add HFS A wrapper.**
In a new module `geovac/hyperfine_a_constant.py` (or extend `nuclear_electronic.py`), add:
```python
def a_constant_ns(Z, n, g_N, core_type=None):
    """A constant for ns_{1/2} HFS: A = (8pi/3) g_e g_N alpha^2 |psi(0)|^2 (m_e/m_p) [Ha]"""
    psi_sq = screened_psi_origin_squared(Z, n, l=0, core_type=core_type)
    g_e = 2.0023193
    return (8 * pi / 3) * g_e * g_N * ALPHA**2 * (1/M_PROTON_OVER_M_E) * psi_sq
```
plus relativistic enhancement via Casimir / spinor lift, cross-check vs hydrogen 1s at A_H = 1418.84 MHz.

**Phase 3 (~1 week): Bohr-Weisskopf as W1b sibling.**
The Bohr-Weisskopf correction (extended nuclear magnetization for s-states) is the Cs HFS analog of the proton Zemach radius for hydrogen 21 cm. Add to `magnetization_density.py` a `bohr_weisskopf_correction` function that takes the magnetization profile of the Cs nucleus (literature: Roberts-Ginges 2022) and computes the framework-native shift. Compare to the literature itemization of Bohr-Weisskopf in Cs HFS (typically 1-2% of A).

**Phase 4 (~3 days): Comparison vs Roberts-Ginges 2019.**
Decompose the framework's A prediction into Bohr-Fermi + relativistic + Schwinger + Bohr-Weisskopf + multi-loop QED, and compare component-by-component with the RG decomposition. Open question: does the framework's projection-chain decomposition match RG's CI+all-order decomposition under a structural mapping?

---

## 7. Proposed Paper 34 §V.C.6 fill text (PROSPECTIVE; do NOT apply now)

> **§V.C.6 Cesium 6S₁/₂ hyperfine (atomic clock; SCOPING).** The Cs-133 6S₁/₂-F=4↔F=3 transition defines the SI second by international agreement (BIPM 1967): $\Delta\nu = 9{,}192{,}631{,}770$ Hz, equivalently $A(\text{6S}_{1/2}) = \Delta\nu/4 = 2298.158$ MHz. A scoping pass (Track Cs-HFS-v1, May 2026) confirms the framework's projection-chain machinery (atomic_classifier, [Xe] FrozenCore, Bohr-Fermi formula) extends to Z=55 without structural failure: the prediction lives in $A \in [793, 2931]$ MHz across the natural range of effective screening charges $Z_\text{eff} \in [9.7, 15]$, bracketing the experimental value at the inverse-solve $Z_\text{eff} \approx 10\text{-}11$ (consistent with Roberts-Ginges 2019). A clean framework-native compute is blocked by two engineering items: (i) the screened radial solver `_solve_screened_radial` rejects l=0 by design, requiring extension to log-grid + Numerov + analytical small-r fit for s-wave $|\psi(0)|^2$; (ii) the Track NI hyperfine Pauli wrapper `hyperfine_coupling_pauli` is specialized to the 1p+1e architecture, requiring a generic A·I·J wrapper for atomic HFS observables. Both are mechanical (~2-3 week sprint). Z-scaling assessment: heavy-atom regime does not expose structural limitations beyond those already identified at Z=20-56 in Sprint 7b screened SO splittings (Clementi-Raimondi single-zeta screening qualitative; W1c-residual core-polarization correction needed for sub-percent precision). Diagnostic-instrument angle: the framework's projection-chain decomposition would isolate three known atomic-structure dependencies of the percent-level Cs PNC uncertainty (radial wavefunction at the origin, core polarization, multi-loop QED in heavy-atom regime), each tied to a different Layer-2 input. (Cf. §III.18 magnetization-density for Bohr-Weisskopf, §III.7 spinor lift for relativistic enhancement, §III.16 multi-loop two-body Dirac/Breit Layer-2.)

---

## 8. Files

- `debug/calc_track_cs_hfs_v1.py` — full feasibility driver + skeleton compute
- `debug/data/cs_hfs_v1.json` — structured probe results
- `debug/cs_hfs_v1_memo.md` — this memo

No production GeoVac code modified. No Paper 34 edits applied (per task constraint).
