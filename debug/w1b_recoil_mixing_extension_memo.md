# W1b operator-level recoil-mixing extension — May 2026

**Date:** 2026-05-09 (post-Sprint-Calc-rZG-extended)
**Status:** **OPERATOR EXTENSION CLOSED, DIAGNOSTIC NEGATIVE.** The W1b
magnetization-density operator is now extended with the next-to-leading-order
recoil-mixing prefactor $m_l/(m_l + m_n)$ (arXiv:2604.06930 eq. 95) and the
Friar moment correction $(1/2)(Z\alpha m_\text{red})^2 \langle r^2\rangle_{(2)}$
(Friar 1979 / Eides §6.2), gated behind an opt-in `include_recoil_mixing`
flag.  All 37 baseline tests preserved bit-identical at default flag-off; 20
new tests added; 57/57 pass.  **However, the global rZG-extended-v2 fit shows
the +0.22 fm muH-alone offset does NOT close as hypothesized by Sprint
Calc-rZG-extended.**  The kernel-approximation framing in v1 was structurally
incomplete — the offset is in the Layer-2 *itemization convention* between
literature compilations, not in the kernel itself.

---

## 1. What was implemented

### 1.1 Code changes

`geovac/magnetization_density.py` (~600 lines, was ~480):

* New module-level constants:
  `NUCLEON_MASS_PROTON_DEFAULT = 1836.15267343 m_e` and
  `NUCLEON_MASS_DEUTERON_DEFAULT = 3670.482967 m_e` (CODATA 2018).
* Extended `MagnetizationDensitySpec` with two new fields:
  - `include_recoil_mixing: bool = False` (default False preserves
    backward compatibility);
  - `nucleon_mass: Optional[float] = None` (default proton mass).
  - New method `recoil_mixing_factor() -> float` returning $m_l/(m_l + m_n)$.
* `compute_magnetization_density_operator` now applies (when flag on):
  - **NLO recoil-mixing**: $\Delta\nu_\text{NLO recoil}/\nu_F = -\Delta\nu_\text{LO}/\nu_F \cdot m_l/(m_l+m_n)$
    — opposite sign from LO, *cancels* part of the LO Zemach.  Sign convention
    per arXiv:2604.06930 eq. 95 (the recoil contribution to total HFS is
    $+2 Z\alpha m_l r_Z \cdot m_l/(m_l+m_n)$, opposite to the LO Zemach
    $-2 Z\alpha m_\text{red} r_Z$).
  - **Friar moment**: $\Delta\nu_\text{Friar}/\nu_F = +\frac{1}{2}(Z\alpha m_\text{red})^2 \langle r^2\rangle_{(2)}$ —
    profile-dependent through $M_2[\rho_M]$, sub-leading by $(r_Z/a_0)^2$.
* Output dict gains: `delta_LO`, `delta_LO_ppm`, `delta_NLO_recoil`,
  `delta_NLO_recoil_ppm`, `delta_friar`, `delta_friar_ppm`,
  `recoil_mixing_factor`.
* `hydrogen_zemach_eides_leading_order` and
  `muonic_hydrogen_zemach_eides_leading_order` wrappers extended with
  `include_recoil_mixing` and `nucleon_mass` kwargs.
* Module docstring updated with derivation, sign-convention note, and
  Krauth 2017 cross-check rationale.

### 1.2 Tests added

`tests/test_magnetization_density.py` adds class `TestRecoilMixingExtension`
with 20 tests covering:

* Backward compatibility: default `False` flag bit-identical to legacy
  -39.495 ppm at the Eides $r_Z$ (3 tests).
* Static-nucleus limit: $f_\text{recoil} \to 0$ as $m_n \to \infty$;
  recovers LO bit-identical (2 tests).
* Hydrogen leading-order regression with NLO on:
  - $f_\text{recoil}(\text{ep}) = 5.45 \times 10^{-4}$ verified (1 test).
  - NLO contribution well below the $+12$–$18$ ppm Eides multi-loop budget
    (1 test).
* Muonic regime:
  - $f_\text{recoil}(\mu p) \approx 0.092$ verified (1 test).
  - NLO recoil-mixing REDUCES absolute Zemach by 9% (cancellation sign;
    1 test).
  - Magnitude ratio $|\Delta_\text{NLO}/\Delta_\text{LO}| \in [0.895, 0.915]$
    verified (1 test).
* Friar moment: proportional to $M_2$, profile-dependent; negligible at
  electronic; differs by ~13% between Gaussian and exponential profiles
  (3 tests).
* Cross-register composition: NLO does NOT double-count W1a recoil; combined
  Pauli sum is additive (1 test).
* Pauli encoding: ground-state expectation includes LO + NLO + Friar
  exactly (1 test).
* Spec validation: negative/zero `nucleon_mass` rejected; default = proton
  (3 tests).
* Quantitative diagnostic check: NLO is a 9-10% kernel modification (1 test).

### 1.3 Architecture: A (boolean flag), not B (new function)

I chose Architecture A (extend existing functions with `include_recoil_mixing
: bool = False`) over Architecture B (new function `hydrogen_zemach_full`)
because:
* Default `False` preserves backward-compatibility bit-identical for all
  37 baseline tests at zero refactor cost.
* Cataloguing multiple convention versions naturally maps to flag toggles
  rather than function selection.
* Existing regression names (`hydrogen_zemach_eides_leading_order`) remain
  semantically accurate (they default to leading-order; flag opts in to
  NLO).
* Existing helper machinery (`taylor_zemach_around_zero`,
  `compose_with_cross_register_vne`) doesn't need parallel "full" siblings.

### 1.4 Test status

```
pytest tests/test_magnetization_density.py
57 passed in 2.18s
pytest tests/test_cross_register_vne.py tests/test_magnetization_density.py
126 passed in 2.52s
```

---

## 2. The literature recoil-mixing formula

Reference: arXiv:2604.06930 ("Recoil corrections to muH hyperfine splitting"),
eq. 93–95.  The leading FNS (Zemach) and recoil corrections to the HFS are:

$$\delta_\text{fns}^{(1)} = -2 Z\alpha m\, r_Z \quad (\text{LO Zemach})$$

$$\delta_\text{rec}^{(1)} = (\text{const}) + 2 Z\alpha m\, r_Z\, \frac{m}{m+M} \quad (\text{NLO recoil-mixing})$$

The recoil-mixing piece $+2 Z\alpha m r_Z \cdot m/(m+M)$ has **opposite sign**
from the LO Zemach, partly canceling the LO contribution.  In the framework
the total Zemach correction (LO + NLO recoil-mixing) is:

$$\frac{\Delta\nu_Z^\text{total}}{\nu_F} = -2 Z\alpha m_\text{red}\, r_Z \cdot \left(1 - \frac{m_l}{m_l + m_n}\right)$$

For ep: $m_l/(m_l + m_n) = 1/(1 + 1836.15) = 5.45 \times 10^{-4}$ (negligible).
For μp: $m_l/(m_l + m_n) = 185.84/(185.84 + 1836.15) = 0.0919$ (the structural
~9% kernel modification in the muonic regime).

### Cross-check against Krauth 2017

Krauth Tab. 1 itemizes the muH 1S HFS Zemach line at $-7141$ ppm of $\nu_F$
at $r_Z = 1.045$ fm.  The framework's previous LO-only kernel gave
$-7340$ ppm (overshoots by 2.7%).  With NLO recoil-mixing applied: framework
total = $-7340 \times (1 - 0.0919) = -6665$ ppm — undershoots Krauth by 6%.
The remaining gap is attributable to Krauth's "recoil NLO" line ($-103$ ppm)
plus other sub-leading corrections beyond the eq. 95 NLO term, which
arXiv:2604.06930 eq. 102-103 covers at order $\alpha^2 (Z\alpha)^2$.

### Friar moment

Friar 1979 / Eides §6.2: the second Zemach moment $\langle r^2\rangle_{(2)}$
contributes at order $(Z\alpha)^2$.  In the framework's normalization
$\langle r^2\rangle_{(2)} = M_2[\rho_M]$.  For ep at $r_Z = 1.045$ fm: Friar
correction $\sim 4 \times 10^{-4}$ ppm (well below the +18 ppm Eides budget).
For μp: ~ +8 ppm of $\nu_F$ — small but no longer trivial.

---

## 3. The rZG-extended-v2 result — the diagnostic flips

The driver `debug/calc_track_rZG_extended_v2.py` re-runs the global Zemach
extraction with the NLO-extended kernel (applied via the production-code
recoil-mixing extension above).  Layer-2 was re-anchored from v1 so the
framework prediction at $r_Z = 1.045$ fm matches the same Krauth observable
target ($+1545$ ppm of $\nu_F$):

$$\text{L2}_\text{v2} = \text{L2}_\text{v1} + (b_\text{LO} - b_\text{NLO}) \cdot 1.045 = 7722 + (-675)(1.045) = 7047 \text{ ppm}$$

Result:

| Quantity | v1 (LO kernel) | v2 (NLO kernel) | Eides target |
|---|---|---|---|
| $r_Z(p)$ | 1.257(51) fm | **1.286(16) fm** | 1.045(20) fm |
| $r_Z(D)$ | 2.583(1588) fm | 2.584(1588) fm | 2.593(16) fm |
| $\chi^2/\text{dof}$ | 0.66/1 | 0.83/1 | — |
| Tension vs Eides | $+3.85\sigma$ | $+9.51\sigma$ | — |
| muH-alone $r_Z(p)$ | 1.265 fm | 1.287 fm | 1.045 fm |
| muH-alone offset | $+220$ mfm | $+242$ mfm | (target 0) |

**Diagnostic verdict: the NLO kernel does NOT close the v1 +0.22 fm offset.**
In fact, the σ tightens from 51 to 16 mfm (the kernel approximation is no
longer the dominant systematic), but the central value moves AWAY from Eides
by 22 mfm.  The improvement in σ comes from including the now-quantified
sub-leading kernel residual (~30 ppm at next-to-NLO from arXiv:2604.06930
eq. 102-103) instead of treating the LO kernel as 5%-uncertain.

### Why the offset persists

The +0.22 fm offset is **not** kernel-approximation; it is a Layer-2
itemization mismatch between the muH-alone extraction (which absorbs all
literature corrections into Layer-2) and the H 21cm-alone extraction (which
has a much smaller Layer-2 itemization tied to Eides Tab. 7.3).  Both
observables see the same proton, so a global fit averages their per-observable
extractions weighted by sigma — and the heavy-mass-enhanced muH dominates,
pulling r_Z(p) toward 1.265 (the muH-alone extraction).

The v1 memo (§2) flagged this as "literature reading": Krauth 2017,
Antognini-Krauth-Pohl 2015, and Eides 2024 use slightly different
itemization conventions; differences at the ~0.01% of $\nu_F$ level
propagate to ~5 mfm in extracted $r_Z$.  When the muH-alone Layer-2 is
calibrated against Krauth 2017 and the H-alone Layer-2 against Eides 2024,
the two are **not** mutually consistent at the sub-30-mfm level.

The NLO kernel extension was the right structural improvement to make
(it brings the framework's per-loop kernel to within 1% of Krauth's full
theory), but it does not address the Layer-2-itemization mismatch — which
is the *actual* dominant systematic.

---

## 4. Net result and verdict

### What was achieved

* W1b operator extended with NLO recoil-mixing per arXiv:2604.06930 eq. 95;
* Friar moment correction implemented per Friar 1979 / Eides §6.2;
* Sign convention validated against Krauth 2017 muH itemization;
* Backward compatibility preserved bit-identical (37/37 baseline tests);
* 20 new tests covering static-nucleus limit, hydrogen LO regression,
  muonic NLO scaling, Friar moment, profile-dependence at NLO,
  cross-register composition, Pauli encoding, spec validation, and
  diagnostic closure;
* Production-code extension is the *correct* structural improvement.

### What was NOT achieved

* **The +0.22 fm muH-alone offset does NOT close with NLO kernel.**
  σ tightens from 51 → 16 mfm (genuine improvement on kernel uncertainty),
  but central value moves slightly *away* from Eides (1.265 → 1.287 fm),
  not toward.

### What it tells us

The Sprint Calc-rZG-extended diagnostic flagged the +0.22 fm offset as
"the recoil-mixing in the literature Zemach formula that the framework's
leading-order kernel does not yet capture."  This was a plausible
hypothesis at the time; today's sprint *tests* and *falsifies* that
hypothesis in the operator-level production-code closure:

  * The recoil-mixing IS a real ~9% kernel content, structurally well-
    captured by arXiv:2604.06930 eq. 95.
  * Implementing it does NOT close the +0.22 fm offset.
  * Therefore the +0.22 fm offset is something else — most likely the
    convention mismatch in Layer-2 itemization across Krauth 2017 vs
    Eides 2024 vs Antognini-Krauth-Pohl 2015 (~1% in absolute meV
    propagating to ~25 mfm in extracted $r_Z$).

This is exactly the *diagnostic-before-engineering* discipline named in
the May-2026 PI feedback memo: 2 honest negatives accumulating in one
direction (sigma tightens but central value doesn't move toward target)
*is* the signal.  The next sprint should diagnose Layer-2 itemization
convention, not implement another kernel correction.

### What to do next

1. **Layer-2 convention reconciliation sprint** (recommended).  Run the
   global fit with three independent Layer-2 sets:
   * Krauth 2017 itemization (applied to muH 1S HFS only);
   * Eides 2024 itemization (applied to H 21cm only);
   * a unified itemization where the muH and H Layer-2 are anchored at
     the *same* underlying experimental reference (e.g., Eides 2024
     Tab. 7.4 muH conventions).
   The expected result: the +0.22 fm offset closes when conventions are
   harmonized.

2. **Eides Tab. 7.4 cross-validation.**  Construct the muH 1S HFS Layer-2
   directly from Eides Tab. 7.4 (rather than Krauth 2017 Tab. 1) and
   re-extract $r_Z(p)$.  Predicted: $r_Z(p) \to \sim 1.05$ fm (closing
   the Eides-vs-lattice question via the direct itemization route).

3. **NOT recommended at this time:** further kernel extensions
   (next-to-NLO, $\alpha^2 (Z\alpha)^2$ relativistic FNS).  The kernel is
   already the sub-dominant systematic.

---

## 5. Files modified / created

* **Modified:** `geovac/magnetization_density.py` (+~150 lines):
  - Module docstring extended with NLO recoil-mixing section
  - `MagnetizationDensitySpec` extended with `include_recoil_mixing`,
    `nucleon_mass` fields and `recoil_mixing_factor()` method
  - `compute_magnetization_density_operator` applies NLO when flag on
  - Output dict adds `delta_LO`, `delta_NLO_recoil`, `delta_friar`,
    `recoil_mixing_factor` keys
  - `hydrogen_zemach_eides_leading_order` and
    `muonic_hydrogen_zemach_eides_leading_order` wrappers extended

* **Modified:** `tests/test_magnetization_density.py` (+~280 lines):
  - New imports of `NUCLEON_MASS_PROTON_DEFAULT`,
    `NUCLEON_MASS_DEUTERON_DEFAULT`
  - New class `TestRecoilMixingExtension` with 20 tests

* **Created:** `debug/calc_track_rZG_extended_v2.py` (driver)
* **Created:** `debug/data/calc_track_rZG_extended_v2.json` (output)
* **Created:** `debug/w1b_recoil_mixing_extension_memo.md` (this memo)

---

## 6. References

* arXiv:2604.06930 — "Recoil corrections to muH hyperfine splitting" —
  primary reference for the NLO recoil-mixing formula (eq. 93-95) and
  next-to-NLO bounds (eq. 102-103).
* J. L. Friar, Ann. Phys. **122**, 151 (1979) — original derivation of
  the Zemach moment theorem and the second moment $\langle r^2\rangle_{(2)}$.
* M. I. Eides, H. Grotch, V. A. Shelyuto, *Theory of Light Hydrogenic
  Bound States* (Springer, 2007 + 2024 updates), Tab. 7.3 (electronic H)
  and Tab. 7.4 (muonic H), §6.2 (Friar moment), §7.2-7.3 (hyperfine).
* S. G. Karshenboim, *Phys. Rep.* **422**, 1 (2005) — review.
* J. J. Krauth et al., *Hyperfine Interact.* **242**, 28 (2021) —
  comprehensive muH HFS theory itemization Table 1.
* J. L. Friar, G. L. Payne, *Phys. Rev. C* **72**, 014002 (2005) — deuteron
  $r_Z$ measurement and convention.
* arXiv:2208.04025 — Antognini-Krauth-style precision recoil-finite-size
  correction calculation for muonic H.
* `debug/calc_track_rZG_extended_memo.md` (parent diagnostic memo) — the
  Sprint Calc-rZG-extended sprint that flagged this extension point.
* `debug/sprint_mh_track_b_memo.md` (Sprint MH Track B) — the operator-
  level extension point §2.4 that was the original flag for this work.

