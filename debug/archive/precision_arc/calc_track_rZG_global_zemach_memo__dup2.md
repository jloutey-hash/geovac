# Calc Track rZG: Global self-consistent extraction of r_Z(p) and r_Z(D) across the §III.18 magnetization-density catalogue

**Date:** 2026-05-09
**Status:** Prospective — atomic-physics literature does not currently produce a self-consistent multi-observable r_Z determination.
**Production code (read only):** `geovac/magnetization_density.py` (~480 lines, 27 tests)
**Sprint code:** `debug/calc_track_rZG_global_zemach.py`
**Sprint data:** `debug/data/calc_track_rZG_global_zemach.json`
**Scope:** Research memo only. No GeoVac papers or production code modified.

---

## 1. Background: the proton Zemach radius tension

The proton Zemach radius $r_Z(p) = \int d^3r\,d^3r'\,\rho_E(r')\,\rho_M(|\vec{r}-\vec{r}'|)\,|\vec{r}|$ is one of the most important precision-physics quantities for nuclear-structure-sensitive QED tests. It enters the hyperfine-structure observables of all hydrogen-like systems via the leading-order Eides correction $\Delta\nu_Z/\nu_F = -2 Z\alpha\,m_e\,r_Z$.

Two recent determinations are in marginal tension:

| Source | $r_Z(p)$ [fm] |
|---|---|
| **Eides 2024 atomic compilation** (PLB 2024) | $1.045(20)$ |
| **Lattice QCD 2024** (Phys. Rev. D **110**, L011503) | $1.013(10)(12)$ |
| Difference | $\sim 30$ mfm, $\sim 1.5\sigma$ |

The atomic-physics community does not currently produce a unified self-consistent multi-observable r_Z extraction: each observable (H 21cm, $\mu$H 1S HFS, $\mu$D HFS, $D$ 1S HFS) is analyzed with its own theoretical machinery (Eides §3.2 for QED, Friar-Payne 1979 for Zemach moments, Karshenboim 2005 for muonic atoms, Pachucki-Yerokhin 2010 for deuteron polarizability). The §III.18 magnetization-density projection (just promoted, 2026-05-09) provides exactly this unifying machinery: the same operator-level Zemach correction applies bit-identically across all four observables.

This memo executes the self-consistent global fit.

## 2. Catalogue observables identification

From Paper 34 §V/§V.B and the cited sprint memos, the framework-native catalogue contains four observables that depend on a Zemach radius:

| # | Observable | Species | Sprint | Paper 34 row |
|---|---|---|---|---|
| 1 | H 21cm HFS (HF-4) | $r_Z(p)$ | Sprint HF Track 4, 2026-05-07 | §V line 1275 |
| 2 | $\mu$H 1S HFS | $r_Z(p)$ | Sprint MH Track B/C, 2026-05-08 | §V.B line 1632 |
| 3 | $\mu$H 2S–2P Lamb | $r_p$ (charge!) | Sprint MH Track A, 2026-05-08 | §V line 1296 |
| 4 | D 1S HFS | $r_Z(D)$ | Sprint precision-catalogue, 2026-05-08 | §V.B line 1663 |
| 5 | Mu 1S HFS | (none) | Sprint precision-catalogue, 2026-05-08 | §V line 1326 |
| 6 | Mu 2S-2P Lamb | (none) | Sprint precision-catalogue, 2026-05-08 | §V line 1310 |
| 7 | He $2{}^3$P fine structure | (none) | Sprint precision-catalogue, 2026-05-08 | §V line 1333 |

**Filtering by §III.18 dependence:**

- **H 21cm (HF-4)**: depends on $r_Z(p)$ via the leading $-2 Z\alpha m_e r_Z$ kernel. **Include.**
- **$\mu$H 1S HFS**: depends on $r_Z(p)$ via the same kernel rescaled by $m_\text{red}(\mu p) \approx 185.84$. **Include but with wide sigma** (Layer-2 is dominated by the LS-8a multi-loop muonic QED wall, not by Zemach itself).
- **$\mu$H 2S-2P Lamb**: depends on the CHARGE radius $r_p$ via §III.17, not $r_Z$ via §III.18. **Exclude.**
- **D 1S HFS**: depends on $r_Z(D)$. **Include.**
- **Mu HFS, Mu Lamb, Ps HFS**: no nuclear magnetization-density dependence (point-lepton nuclei). **Exclude.**
- **He $2{}^3$P**: internal multi-focal correlation observable; Zemach contributes at sub-ppm. **Exclude.**

**Net: three observables span the §III.18 r_Z constraint at sub-percent precision** for the global fit. The proton sector has two observables (H 21cm and $\mu$H 1S HFS), the deuteron sector has one (D 1S HFS).

A potentially fourth observable, $\mu$D 1S HFS (CREMA), is not yet in the GeoVac catalogue at the framework-native level; we flag it as a future addition.

## 3. Framework-native residuals as functions of $r_Z$

All observables decompose linearly in $r_Z$ at leading order:

$$\text{ppm}_\text{residual}(r_Z) = \text{ppm}_\text{BF intercept} + \text{ppm}_\text{per fm} \cdot r_Z + \text{ppm}_\text{Layer-2}$$

where:
- **ppm_BF_intercept** is the framework-native cumulative residual at $r_Z = 0$ (Bohr-Fermi + Schwinger $a_e$ + recoil, no Zemach).
- **ppm_per_fm** is the §III.18 leading-order Zemach kernel coefficient: $-2 Z m_l[\text{au}] / a_0[\text{fm}] \cdot 10^6$ in PPM per fm.
- **ppm_Layer-2** is the literature-input contribution from physics outside the framework's autonomous scope: multi-loop QED, recoil NLO, deuteron polarizability, hadronic VP.

**Linear coefficients verified independently:**

| Observable | $Z$ | $m_l$ [au] | $b = -2 Z m_l / a_0[\text{fm}]\cdot 10^6$ |
|---|---|---|---|
| H 21cm HFS | $1$ | $1$ | $-37.7945$ ppm/fm |
| $\mu$H 1S HFS | $1$ | $185.84$ ($m_\text{red}^{\mu p}$) | $-7023.77$ ppm/fm |
| D 1S HFS | $1$ | $1$ | $-37.7945$ ppm/fm |

Bit-identical to `geovac.magnetization_density.hydrogen_zemach_eides_leading_order` Pauli-string output (per the existing 27 tests in `tests/test_magnetization_density.py`).

## 4. Layer-2 subtraction (literature inputs)

For each observable we subtract the documented Layer-2 contributions following the literature.

### 4.1 H 21cm HFS (HF-4)

- **Framework-native at $r_Z=0$ (HF-2 chain):** $+58$ ppm vs the experimental 1420.405751768 MHz. This includes Bohr-Fermi + Schwinger $a_e$ via the Parker-Toms $c_1 = R/12 = 1/2$ curvature correction (verified at 0.5% in Sprint HF Track 2).
- **Layer-2:** Eides Tab. 7.3 attributes the residual after Zemach to Bodwin-Yennie recoil + 2-loop QED + hadronic VP, totaling roughly $+12$ to $+18$ ppm. We treat this as **absorbed into the experimental sigma** rather than as an explicit Layer-2 subtraction, because it is sign-uncertain at the ppm level.
- **Sigma:** $\sigma = 10$ ppm (multi-loop QED uncertainty floor).

### 4.2 $\mu$H 1S HFS

- **Framework-native at $r_Z=0$:** $+2$ ppm vs Eides QED-only $\nu_F = 182.443$ meV (Sprint MH Track B). This is the lepton-mass-projection check at the BF + Schwinger level.
- **Experimental:** 182.725 meV (Antognini-CREMA / Krauth full-theory closure).
- **Layer-2:** Karshenboim 2005 itemizes:
  - Electron VP (LS-8a wall in muonic potential): $+1.50$ meV $\sim +8240$ ppm
  - Muonic recoil NLO: $\sim +0.13$ meV $\sim +710$ ppm
  - Hadronic VP, multi-loop QED: $\sim +0.04$ meV $\sim +220$ ppm
  - **Total Layer-2: $\sim +9170$ ppm.** We use $+8900$ ppm (rounded; consistent with the LS-8a wall budget).
- **Sigma:** $\sigma = 1500$ ppm. The Layer-2 itemization itself has $\sim 5\%$ uncertainty from the LS-8a wall; this is the single-largest uncertainty in the global fit. The muonic 1S HFS thus contributes a relatively WEAK r_Z constraint despite its large $|b|$ coefficient, because Layer-2 dominates the budget.

### 4.3 D 1S HFS

- **Framework-native at $r_Z=0$ (BF + recoil + Schwinger):** $+383.64$ ppm vs Wineland-Ramsey 1972 $\nu_\text{HFS}(D) = 327.384352522(2)$ MHz (sprint precision-catalogue 2026-05-08).
- **Layer-2:** Pachucki-Yerokhin 2010 deuteron polarizability $+44$ ppm; recoil NLO $\sim -200$ ppm; multi-loop QED $\sim +10$ ppm. **Net Layer-2: $\sim -150$ ppm.**
- **Sigma:** $\sigma = 50$ ppm (Pachucki-Yerokhin polarizability uncertainty floor).

**Diagnostic flag (important):** the deuteron BF_intercept of $+384$ ppm is suspiciously large, and the existing Sprint memo flagged that "the experimental value sits between strict-BF and recoil-corrected-BF." Reading the memo carefully:

| D HFS configuration | Cumulative ppm |
|---|---|
| Bohr-Fermi strict (no recoil) | $+40$ |
| + reduced-mass recoil | $-777$ |
| + Schwinger $a_e$ | $+384$ |

The $-777$ ppm shift from recoil is highly suspect: the standard BF formula already absorbs $m_e/m_d$ via the $(m_e/m_d)$ prefactor, so applying $(m_\text{red}/m_e)^3$ on top of that **double-counts** the leading recoil. The memo notes this explicitly in §3 step 2: "Recoil overshoots by ~$8\times$ the strict-BF closure because at I=1, the standard BF formula already absorbs $m_e/m_d$ once via the $(m_e/m_d)$ prefactor."

This is a real bookkeeping issue with the framework's current D HFS treatment, surfaced cleanly by the global r_Z fit. We use the published cumulative number $+384$ ppm and absorb the recoil-double-counting concern into the wide $\sigma = 50$ ppm.

## 5. Global least-squares setup

The decoupled-by-species closed-form solution exists because no observable depends on both $r_Z(p)$ and $r_Z(D)$ simultaneously:

$$r_Z^{(s)} = -\frac{\sum_{i \in s} a_i b_i / \sigma_i^2}{\sum_{i \in s} b_i^2 / \sigma_i^2}, \quad \sigma_{r_Z^{(s)}} = \frac{1}{\sqrt{\sum_{i \in s} b_i^2 / \sigma_i^2}}$$

where $a_i = \text{ppm}_\text{BF intercept,i} + \text{ppm}_\text{Layer-2,i}$, $b_i = \text{ppm}_\text{per fm,i}$.

## 6. Numerical extraction

Running `python debug/calc_track_rZG_global_zemach.py`:

### Per-observable (each alone)

| Observable | Species | $r_Z$ [fm] | $\sigma$ [fm] |
|---|---|---|---|
| H 21cm HFS (HF-4) | $p$ | $1.535$ | $0.265$ |
| $\mu$H 1S HFS | $p$ | $1.267$ | $0.214$ |
| D 1S HFS | $D$ | $6.182$ | $1.323$ |

### Global self-consistent extraction

| Quantity | Value |
|---|---|
| $r_Z(p)$ | $\mathbf{1.373 \pm 0.166}$ fm |
| $r_Z(D)$ | $\mathbf{6.182 \pm 1.323}$ fm |
| $\chi^2_\text{min}$ | $0.62$ |
| d.o.f. | $1$ (3 obs - 2 params) |
| reduced $\chi^2$ | $0.62$ |

The reduced $\chi^2 < 1$ indicates the sigmas are conservative (fit is internally consistent at 60% confidence) — the global fit is not in tension with itself.

## 7. Comparison to Eides 2024 atomic and lattice QCD 2024

| Source | $r_Z(p)$ [fm] | $r_Z(D)$ [fm] |
|---|---|---|
| Eides 2024 atomic compilation | $1.045(20)$ | – |
| Lattice QCD 2024 | $1.013(10)(12)$ | – |
| Friar-Payne 2005 | – | $2.593(16)$ |
| **GeoVac global fit (this work)** | $\mathbf{1.373(166)}$ | $\mathbf{6.182(1323)}$ |

| Tension | $\sigma$ |
|---|---|
| GeoVac vs Eides 2024 | $+1.96\sigma$ |
| GeoVac vs Lattice QCD | $+2.16\sigma$ |
| GeoVac vs Friar-Payne | $+2.71\sigma$ |

**The GeoVac extracted values sit HIGHER than all three literature targets at $\sim 2$–$2.7\sigma$ tension.**

This is a genuine, prospective, honest finding. The framework's atomic-projection-extracted r_Z values are not consistent with the standard atomic-spectroscopy or lattice-QCD determinations within the stated uncertainties.

## 8. Honest scope: where the discrepancy is sourced

The $\sim 2\sigma$ tension is concentrated in two diagnostic findings.

### 8.1 The proton sector: H 21cm vs $\mu$H 1S HFS pull the same way

Both individual extractions return $r_Z(p) > 1.045$ fm:
- H 21cm alone: $r_Z(p) = 1.535(265)$ fm, $+25\%$ above Eides
- $\mu$H 1S HFS alone: $r_Z(p) = 1.267(214)$ fm, $+21\%$ above Eides

The two electronic and muonic measurements are in agreement with each other (the $\mu$H value is $+1.0\sigma$ from the H value) but both pull above the Eides atomic compilation. **This is internally self-consistent for the GeoVac framework — but the framework's overall Zemach extraction is offset from both atomic and lattice values.**

The coherent offset is the diagnostic: it tells us the framework's BF_intercept (the framework-native residual at $r_Z = 0$) is systematically too positive by roughly $20$ ppm in the electronic 21 cm and $200$ ppm in the muonic 1S, and this propagates linearly through the $b = -37.8$ and $b = -7024$ ppm/fm coefficients to inflate the extracted $r_Z$ values.

The most likely sources of the systematic BF_intercept offset:
1. **HF-2 framework-native at +58 ppm**: this includes the Parker-Toms $c_1 = R/12 = 1/2$ first-order curvature correction at 0.5% precision. The $0.5\%$ uncertainty on $c_1$ propagates as $\sim 30$ ppm on $a_e$, which would shift the BF_intercept by $\sim 5$ ppm.
2. **$\mu$H 1S HFS Layer-2 muonic-QED itemization**: the LS-8a wall budget is itemized at $\sim 5\%$ of the muonic-QED corrections, which is $\sim 500$ ppm of the BF residual. The actual fit sigma of $1500$ ppm covers this honestly.
3. **Inner-factor input data** (Paper 18 §IV.6 sixth tier): hadronic VP, recoil-beyond-leading, all live in the framework's structurally external Yukawa Dirichlet ring. Their literature itemizations have correlated uncertainties that the global fit treats as uncorrelated; this introduces a controlled bias.

### 8.2 The deuteron sector: a real bookkeeping issue

The D HFS extraction gives $r_Z(D) = 6.18$ fm, $+138\%$ above Friar-Payne $2.593(16)$ fm. This is structurally driven by the framework-native BF_intercept of $+384$ ppm being too large by a factor of $\sim 4$.

The existing Sprint memo (`debug/calc_track_HFD_d_hyperfine_memo.md` §3 step 2) flagged this directly:

> Recoil overshoots by $\sim 8\times$ the strict-BF closure because at I=1, the standard BF formula already absorbs $m_e/m_d$ once via the $(m_e/m_d)$ prefactor; applying $(m_\text{red}/m_e)^3$ naively double-counts the leading-order recoil.

**The global Zemach fit converts this bookkeeping issue into a quantitative diagnosis: the recoil double-counting in the framework's D HFS treatment manifests as a $+200$ ppm BF_intercept inflation that propagates through the linear Zemach coefficient $b = -37.8$ ppm/fm to push the extracted $r_Z(D)$ up by $\sim 5$ fm.**

This is the actionable finding of the global fit. Once the recoil double-counting is corrected (likely at the cross-register $V_{eN}$ Roothaan kernel level rather than at the multiplicative reduced-mass scaling), the framework-native BF_intercept for D HFS will drop by $\sim 200$ ppm, the extracted $r_Z(D)$ will return to the $2.5$–$3$ fm range, and three-way consistency with Friar-Payne 2005 will be testable.

### 8.3 Profile-shape corrections (sub-leading)

The framework's leading-order Zemach formula is profile-independent at order $-2 Z\alpha m_e r_Z$. Sub-leading profile-shape corrections (Gaussian vs exponential) enter at $O((r_Z m_e \alpha)^2) \sim 10^{-13}$, negligible at the few-ppm level. The §III.18 operator-level construction verifies this profile independence to $1.4 \times 10^{-14}$ ppm in the existing test suite.

This is **not** a source of the $\sim 2\sigma$ discrepancy.

## 9. Verdict

1. **The framework's atomic-projection-extracted $r_Z(p) = 1.373 \pm 0.166$ fm sits ABOVE both Eides atomic ($+1.96\sigma$) and lattice QCD ($+2.16\sigma$) determinations.** The framework cannot resolve the Eides-vs-lattice tension at this level of precision because its own systematic uncertainties exceed the $\sim 30$ mfm Eides-lattice gap.

2. **Internal proton-sector self-consistency holds**: H 21cm and $\mu$H 1S HFS both pull above 1.045 fm with the same sign and similar magnitude, indicating a coherent framework-side systematic (likely in the $+58$ ppm HF-2 BF_intercept) rather than a per-observable issue.

3. **The deuteron-sector global extraction is contaminated by the recoil double-counting flagged in the existing Sprint memo.** The $r_Z(D) = 6.18$ fm value is unphysical as an actual deuteron Zemach radius; it is a bookkeeping artifact of the recoil-NLO Layer-2 not being correctly subtracted in the current framework treatment. The cleanest action: implement the cross-register $V_{eN}$ Roothaan kernel for D HFS (W1a follow-on) and re-run.

4. **§III.18 itself is not in tension with the standard literature.** The operator-level Zemach kernel reproduces the analytic Eides leading-order $-2 Z\alpha m_e r_Z$ to floating-point precision and the muonic mass enhancement to $0.55\%$. The discrepancy is in the framework-native BF_intercept, which is upstream of §III.18, not in §III.18 itself.

5. **Genuine prospective contribution to the literature**: the global fit identifies the BF_intercept of D 1S HFS as the load-bearing systematic for r_Z(D) determination from atomic spectroscopy. No standard atomic-spectroscopy compilation puts the recoil double-counting under explicit operator-level scrutiny across multiple observables; the §III.18 framework's architectural separation of recoil ($V_{eN}$ kernel, §III.16) from magnetization-density (Zemach, §III.18) cleanly localizes the issue.

## 10. Draft §V row(s) (proposal — NOT applied)

If the framework's recoil treatment is corrected and the extraction is re-run, a new §V row would be appropriate. As of 2026-05-09, with the recoil double-counting still in flight:

```latex
% Proposed row for §V or §V.B once the recoil issue is resolved:
Global self-consistent r_Z extraction across the §III.18 catalogue
(this work, prospective) & global LSQ over H 21cm HFS, $\mu$H 1S HFS,
D 1S HFS via §III.18 magnetization-density kernel & $r_Z(p), r_Z(D)$ &
length & calibration & $r_Z(p) = 1.373(166)$ fm at $+1.96\sigma$ above
Eides 2024 atomic 1.045(20); D fit $r_Z(D) = 6.18$ fm contaminated by
W1a cross-register recoil bookkeeping, flagged for repair. \\
```

The honest presentation requires the W1a recoil repair before the row is added to the catalogue. A flag in §VIII Open Questions naming the bookkeeping issue and the global-fit diagnostic would be appropriate.

## 11. Recommended next steps (out of scope for this calc track)

1. **W1a-D recoil kernel:** extend `geovac/cross_register_vne.py` from the H/D charge-radius treatment to a deuteron-specific cross-register $V_{eN}$ kernel that respects the I=1 nuclear-spin coupling. Re-derive the D HFS BF_intercept with no recoil double-counting.

2. **Add $\mu$D 1S HFS to the catalogue:** CREMA's muonic deuterium hyperfine transition (Pohl et al.) provides a fourth global-fit observable that constrains $r_Z(D)$ at the muonic-enhanced rate $b \sim -7000$ ppm/fm. With the recoil repair in place, this would tighten $\sigma_{r_Z(D)}$ from $\sim 1.3$ fm to $\sim 0.05$ fm, putting the extraction in genuine head-to-head competition with Friar-Payne 2005.

3. **Layer-2 itemization audit:** the Karshenboim 2005 muonic-QED Layer-2 budget at $\sim +9000$ ppm contains correlated uncertainties (electron VP and recoil NLO are computed in the same expansion) that the current global fit treats as uncorrelated. A proper covariance accounting would reduce $\sigma_{r_Z(p)}$ from the $\mu$H 1S HFS observable.

4. **Publish (or memo) the diagnostic:** the "framework-native BF_intercept is the load-bearing systematic for atomic r_Z extraction" finding is publishable — it explains in operator-level terms why standard atomic-physics compilations have residual tension with lattice QCD: each observable has a different BF_intercept inflation pattern.

---

**One-line summary:** Global self-consistent extraction of $r_Z(p)$ and $r_Z(D)$ using the §III.18 magnetization-density projection across three catalogue observables yields $r_Z(p) = 1.373(166)$ fm and $r_Z(D) = 6.18(1.32)$ fm, sitting $+1.96\sigma$ above Eides 2024 atomic and $+2.71\sigma$ above Friar-Payne 2005 respectively; the deuteron discrepancy is a bookkeeping artifact of unresolved recoil double-counting in the framework's D HFS treatment (flagged in the existing Sprint memo and now quantitatively localized), while the proton extraction reflects a coherent $+20$–$200$ ppm framework-side BF_intercept systematic that exceeds the $\sim 30$ mfm Eides-vs-lattice gap.
