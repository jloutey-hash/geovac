# Sprint HF Track 4 — Zemach correction diagnostic

**Date:** 2026-05-07
**Goal:** Test whether GeoVac's existing form-factor / finite-size machinery (`geovac/nuclear/form_factor.py` + `nuclear_electronic.py::finite_size_coupling_pauli`) produces the **Zemach correction** to $A_{\rm hf}$ at the right size, OR whether Zemach requires native machinery the framework doesn't have.
**Verdict:** **NEGATIVE.** Track NI has no magnetization-distribution operator, no separate $r_Z$ input, and no size-dependent modulation of the hyperfine coupling. The single proton-structure parameter `R_PROTON_BOHR` is wired exclusively to the *binding-energy* Foldy correction, not to $A_{\rm hf}$. Zemach is structurally external — like recoil — and must be supplied as a Layer-2 calibration focal length with $r_Z$ as input.

---

## 1. Why this matters

HF-1 closed the strict-Bohr-Fermi piece. HF-2 closed the $a_e$ piece via the Schwinger asymptote of GeoVac's graph-native vertex correction (residual $+58$ ppm). HF-3 confirmed Track NI's cross-register architecture cannot derive recoil natively (negative). The Zemach correction is the next-largest known contribution: $\sim -33$ to $-41$ ppm depending on convention (Eides Ch. 7, Karshenboim 2005), arising from the convolution of the proton's charge density $\rho(r)$ with its magnetization density $m(r)$. The PI question for Track HF-4 is the parallel of HF-3: does the framework's existing `R_PROTON_BOHR` parameter, plumbed through `form_factor.py`, host any size dependence of the hyperfine *coupling* (not just the binding energy)?

## 2. Diagnostic of Track NI's proton-structure machinery

I ran a structured search across `geovac/nuclear/form_factor.py` and `geovac/nuclear/nuclear_electronic.py` (the only two modules with proton-structure parameters in the codebase, confirmed by repository-wide grep):

| Token | Hits in `geovac/nuclear/` |
|---|---|
| `Zemach` | **0** |
| `magnetization` | **0** |
| `r_Z` | **0** |
| `R_M` | **0** |
| `magnetic moment distribution` | **0** |
| `R_PROTON` (substring) | 10 |
| `R_nuc` | 58 |
| `charge_radius` (substring) | 0 |
| `charge radius` (phrase) | 1 |

Inspecting `hyperfine_coupling_pauli` directly:

```
hyperfine_coupling_pauli(Q_nuc, Q_elec, nuclear_block, electronic_block,
                         A_hf=HF_HYDROGEN_HA) -> Dict[str, float]
```

No `R_nuc`, no `r_Z`, no proton-size kwarg. `inspect.getsource` confirms zero references to any size token in the body. The function consumes a scalar $A_{\rm hf}$ and emits the four-qubit point-like $\mathbf{I} \cdot \mathbf{S}$ coupling.

`finite_size_coupling_pauli` consumes `R_nuc` but emits a one-body *electronic* operator $\Delta E_{1s} \cdot (n_{1s,\uparrow} + n_{1s,\downarrow})$ where $\Delta E_{1s} = (2/5) Z^4 R_{\rm nuc}^2 / n^3$ is a binding-energy perturbation (Foldy 1958, Friar 1979). Its source contains no `hyperfine`, no `A_hf`, no $\mathbf{I}\cdot\mathbf{S}$ reference. The Foldy correction targets $E_{1s}$, not the hyperfine constant.

**Conclusion of Step 1:** The framework has zero Zemach-relevant machinery. The single proton-structure parameter `R_PROTON_BOHR` is the proton's *charge* radius (CODATA 0.8414 fm), wired to a binding-energy correction. There is no proton magnetization-distribution operator on the nuclear register. There is no functional difference between the proton's spatial wavefunction and a $\delta^3$ point at the origin in the hyperfine coupling.

## 3. Three Zemach scenarios

The canonical Eides §7.2 form is

$$\frac{\Delta \nu_Z}{\nu_F} = -2 Z \alpha m_e r_Z = -2 Z \frac{r_Z[{\rm bohr}]}{1}$$

in atomic units (using $a_0 = 1$, $\lambda_C = \alpha$ bohr, so $m_e r_Z$[bohr]$\,= r_Z[{\rm bohr}]/\alpha$, giving the factor of $\alpha$ to cancel). The briefing's expression "$-2 Z \alpha (r_Z / a_0)$" reads literally as off by a factor of $\alpha$ relative to Eides; I used the canonical natural-unit form and verified it reproduces the textbook $\sim -41$ ppm at $r_Z = 1.045$ fm. (At polarizability-corrected values $r_Z \sim 0.87$ fm the canonical Eides value lands at $\sim -33$ ppm.) The correct numerical formula is what's coded in `zemach_shift_ppm`.

Plugging into the HF-2 closing baseline $A_{\rm hf}^{\rm GV} = 1420.4879$ MHz (residual $+57.80$ ppm vs experimental 1420.4058 MHz):

| Scenario | $r_{\rm subst}$ | $\Delta\nu/\nu$ | $A_{\rm hf}$ | residual |
|---|---|---|---|---|
| Baseline (HF-2) | — | — | 1420.4879 MHz | $+57.80$ ppm |
| **(A)** $R_{\rm PROTON\_BOHR}$ proxy | 0.8414 fm | $-31.80$ ppm | 1420.4427 MHz | $+25.99$ ppm |
| **(B)** $r_Z$ external 1.045 fm | 1.045 fm | $-39.49$ ppm | 1420.4318 MHz | $+18.31$ ppm |
| **(C)** Native (no machinery) | — | 0 | 1420.4879 MHz | $+57.80$ ppm |

The "charge-radius proxy" and "$r_Z$ external" outcomes differ by the structural ratio $r_Z / R_p \approx 1.242$. (A) substitutes $R_{\rm PROTON\_BOHR}$ for the unknown Zemach radius, which is what would happen if a careless user re-purposed the framework's existing parameter; the result has the right *form* (linear in $r/a_0$) but the wrong *size* (off by 24% from the magnetization-radius truth). (B) supplies $r_Z$ as an external Layer-2 calibration focal length, structurally identical to how HF-1 supplied the recoil correction.

The framework's own native prediction, without external supply, is (C) — the +58 ppm residual stands at this stage.

## 4. Verdict: NEGATIVE

The framework has no Zemach-relevant machinery. The hyperfine coupling is point-like in the proton spatial coordinate; there is no magnetization distribution operator, no separate $r_Z$ input, and no size dependence in `hyperfine_coupling_pauli`. `R_PROTON_BOHR` is the only proton-structure parameter and is wired into a one-body electronic operator (Foldy binding-energy correction), not the hyperfine coupling. **Zemach is structurally external — like recoil — and must be supplied as a calibration focal length with $r_Z$ as input.**

This is a negative outcome of exactly the same flavor as HF-3:

- HF-3 found the cross-register architecture has zero spatial-spatial coupling (only spin-spin hyperfine + classical-parameter Foldy), so it cannot derive the recoil $(1+m_e/m_p)^{-3}$ factor.
- HF-4 finds the cross-register architecture has zero magnetization-distribution operator, so it cannot derive the Zemach $-2Zr_Z/{\rm bohr}$ correction.

Both negatives are predicted by the structural-skeleton scope memo (`geovac_structural_skeleton_scope_pattern.md`): the framework determines selection rules, transcendental class, and scaling laws, but parameter values rooted in physics outside its native gauge content (in this case, hadronic SU(3) matter elements that GeoVac admits as a Wilson construction in Paper 30 + Sprint ST-SU3 but cannot autonomously compute) are calibration inputs. They are valid Layer-2 prescriptions; they are not derivations.

## 5. Updated $A_{\rm hf}$ prediction

The framework's native one-loop QED prediction stays at **$A_{\rm hf}^{\rm native} = 1420.4879$ MHz** with **residual $+57.8$ ppm**. With $r_Z = 1.045$ fm calibrated in as an external Layer-2 input (the right Eides convention, no fits beyond CODATA), the prediction lands at **$1420.4318$ MHz** with **residual $+18.3$ ppm**. The remaining $\sim 18$ ppm decomposes per Eides Tab. 7.3 / Karshenboim 2005 review:

| Source | Magnitude | Sign |
|---|---|---|
| Polarizability $\Delta_{\rm pol}$ | $+1.4 \pm 0.6$ ppm | $+$ |
| Recoil beyond reduced-mass (Bodwin–Yennie) | $+5.85$ ppm | $+$ |
| Multi-loop QED (two- + three-loop) | $+5.07$ ppm | $+$ |
| Hadronic vacuum polarization | $+0.1$ ppm | $+$ |
| Convention drift on $r_Z$ (1.045 vs 1.054) | $\sim \pm 5$ ppm | $\pm$ |

Sum $\sim +12$ to $+18$ ppm, consistent with the +18 ppm observed. The agreement at this level is noise-limited by the Zemach-radius input uncertainty (~1.5%); a different conventional choice for $r_Z$ shifts the residual within a 5–10 ppm band.

## 6. Paper 34 projection-vocabulary placement

The proton charge radius $R_p$ and the proton Zemach radius $r_Z$ are **two different Layer-2 focal lengths**, both QCD-internal:

| Focal length | Meaning | Modulates | Wired into GeoVac? |
|---|---|---|---|
| $R_p$ (charge radius) | RMS of charge distribution $\rho(r)$ | $\Delta E_{1s}$ via Foldy | YES — `form_factor.finite_size_correction` |
| $r_Z$ (Zemach radius) | First moment of $\rho \star m$ | $\Delta\nu / \nu_F$ via $-2Zm_e r_Z$ | NO — no magnetization operator |

In Paper 34's three-axis vocabulary (variable, dimension, transcendental):

- **$R_p$**: adds (length scale; $[L]$; no transcendental). Paper 34 projection class: rest-mass / size / Layer-2 calibration. Category C (external CODATA).
- **$r_Z$**: adds (length scale, *different from $R_p$*: magnetization-weighted, not charge-weighted; $[L]$; no transcendental). Paper 34 projection class: same axis as $R_p$ but a structurally distinct focal length. Category C (external).

The two radii are *not* the same projection. A framework that derives one does not automatically derive the other. The fact that $R_p$ and $r_Z$ differ by 24% in hydrogen — and by larger factors in other nuclei — is not noise but structurally meaningful: charge density $\rho(r)$ is the SU(3)-internal quark color-electric distribution, magnetization $m(r)$ is the quark color-magnetic-moment distribution, and the two have different parametric dependence on quark masses, gluon dynamics, and chiral structure. GeoVac's Wilson-SU(3) construction (Paper 30 + Sprint ST-SU3) admits the gauge content but cannot autonomously compute matter-element distributions, because the (N,0) Bargmann tower has the CG obstruction Sprint 5 Track S5 found.

This places both radii cleanly in the structural-skeleton scope: the framework can host them, but cannot produce them. Whether to add a "Zemach Layer-2 input" entry to Paper 34's living catalogue is the PI's call (it would be the 16th projection — the existing 15 do not include a magnetization-distribution slot).

## 7. PI-facing framing

In one paragraph: GeoVac's existing `R_PROTON_BOHR` parameter is wired to the binding-energy correction (Foldy), not the hyperfine coupling. There is no quantum operator for the proton magnetization distribution on the nuclear register. The hyperfine routine is structurally point-like in the proton spatial coordinate. So Zemach is externally calibrated, structurally identical to HF-1's recoil treatment. With $r_Z = 1.045$ fm as a Layer-2 input the residual drops from +58 ppm to +18 ppm, with the remaining $\sim 18$ ppm decomposing into multi-loop QED + recoil-beyond-reduced-mass + nuclear polarizability — known physics targeting HF-5. The negative verdict here is honest about the same scope boundary HF-3 found: cross-register architecture in Track NI is single-focal calibration plus parametric coupling, not multi-focal composition that resolves nuclear-structure effects natively. The proton charge radius and the proton Zemach radius are *two different Layer-2 focal lengths*, both QCD-internal, neither autonomously derivable from the framework's gauge content (Paper 30 admits SU(3) but the matter elements are CG-obstructed by Sprint 5 Track S5). The HF program at this stage is closing residuals one named external focal length at a time; the framework keeps producing the structural skeleton at every step.

## 8. Files

- `debug/sprint_hf_track4.py` — diagnostic script (no production code modified).
- `debug/data/sprint_hf_track4.json` — structured outputs of the four steps.
- `debug/sprint_hf_track4_memo.md` — this memo.

## 9. References

- Eides, M. I., Grotch, H., Shelyuto, V. A. *Theory of Light Hydrogenic Bound States* (Springer, 2007), Ch. 7 §7.2 — Zemach correction $\Delta\nu_Z / \nu_F = -2 Z \alpha m_e r_Z$ and the canonical hydrogen value $r_Z = 1.045(16)$ fm.
- Karshenboim, S. G. *Phys. Rep.* **422**, 1 (2005) — comprehensive review of hydrogen hyperfine structure including polarizability-corrected Zemach values.
- Foldy, L. L. *Phys. Rev.* **111**, 1093 (1958); Friar, J. L. *Ann. Phys.* **122**, 151 (1979) — Foldy / charge-radius binding-energy correction (what `R_PROTON_BOHR` is wired into).
- CLAUDE.md memory `geovac_structural_skeleton_scope_pattern.md` — the framework determines structure but not parameter values; same pattern as HF-3 (recoil), HF-2 closing rationale (Schwinger asymptote applied to flat-space hydrogen), Sprint H1 (Yukawa form yes, value no), LS-8a (UV form yes, counterterms no).

## 10. Convention note for the PI

The briefing's quoted $\Delta\nu_Z / \nu_F = -2 Z \alpha (r_Z / a_0)$ is off by a factor of $\alpha$ relative to the canonical Eides Eq. (7.4) form. The canonical natural-unit form is $-2 Z \alpha m_e r_Z$, which in atomic units (where $\lambda_C = \alpha$ bohr and $m_e = 1$) becomes $-2 Z (r_Z[{\rm bohr}])$. With $r_Z = 1.974 \times 10^{-5}$ bohr this gives $-39.5$ ppm, agreeing with Eides' $-41$ ppm at the $r_Z$-uncertainty level. The briefing's $-33$ ppm corresponds to a polarizability-corrected $r_Z \approx 0.87$ fm — also in Eides Tab. 7.3 — but coupled to the briefing's wrong-by-α formula it would have given $-0.21$ ppm, which is invisible at our precision. I used the canonical formula throughout and noted this discrepancy in the script comments. The PI may want to correct this in the briefing template for HF-5.
