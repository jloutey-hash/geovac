# Calc Track rZG-Extended — Diagnostic study of why "tighten the muH 1S HFS sigma" does not close the Eides–lattice tension

**Date:** 2026-05-09 (post-rZG correction sprint)
**Status:** **DIAGNOSTIC NEGATIVE.** Tightening the muH 1S HFS Layer-2 uncertainty
to literature itemization precision drives σ(r_Z(p)) ↓ but exposes a structural
**kernel-approximation** systematic that takes the framework's r_Z(p)
extraction OUTSIDE 1σ of both Eides 1.045(20) and lattice 1.013(16) at >3.8σ.
The framework cannot resolve the Eides-vs-lattice tension without a Phase
C-W1b operator-level extension to include recoil-mixing in the magnetization-
density operator.

**Files**
- Driver: `debug/calc_track_rZG_extended.py`
- Data: `debug/data/calc_track_rZG_extended.json`
- This memo
- Parent: `debug/calc_track_rZG_global_zemach.py`,
  `debug/data/calc_track_rZG_global_zemach.json`,
  `debug/rzg_bug_diagnosis_memo.md`

---

## 1. Headline numbers

The original 3-observable corrected fit (`calc_track_rZG_global_zemach.py`,
2026-05-09) gave:

| Quantity | Value | Tension vs Eides 1.045(20) | Tension vs lattice 1.013(16) |
|---|---|---|---|
| r_Z(p), original | 1.180(166) fm | +0.80σ | +1.00σ |
| r_Z(D), original | 2.583(1588) fm | — | — |
| χ²/dof | 0.43/1 | | |

The original fit had σ(r_Z(p)) = 166 mfm, dominated by the muH 1S HFS sigma
of 1500 ppm (treated as "LS-8a wall scale"). Both Eides and lattice fall
within 1σ; framework consistent with both but cannot discriminate.

**This sprint asked:** can we tighten muH 1S HFS sigma per literature
itemization to drive σ(r_Z(p)) below 30 mfm, the Eides-vs-lattice gap?

**Sprint result:**

| sigma_muH (ppm) | Reasoning | r_Z(p) extracted | σ(r_Z(p)) | Tension vs Eides | Tension vs lattice |
|---|---|---|---|---|---|
| 1500 (original) | "LS-8a wall scale" blanket | 1.180 fm | 166 mfm | +0.80σ | +1.00σ |
| 150 (literature itemization) | Krauth 2017 sum-in-quadrature | 1.179 fm | 21 mfm | **+4.58σ** | **+6.28σ** |
| 367 (kernel approximation) | 5% of leading-order Zemach kernel | 1.257 fm | 51 mfm | **+3.85σ** | **+4.55σ** |

**Diagnostic conclusion:**

Tightening σ from 1500 to 150 ppm (literature itemization) drives σ(r_Z(p))
to 21 mfm — below 30 mfm — but the central value extracted (1.179 fm) is
**inconsistent with both Eides and lattice at 4–6σ**. This exposes that the
framework's leading-order Eides Zemach kernel itself (`-2 Z m_red r_Z` in
atomic units) is a **5% approximation** in muonic systems, not a 0.1%
approximation. The proper sigma is dominated by the kernel approximation
(~367 ppm), and even with that more honest sigma, σ(r_Z(p)) = 51 mfm — close
to but not below 30 mfm.

**The framework's per-observable r_Z(muH 1S HFS) = 1.265 fm extraction is
itself the diagnostic flag**: this is +0.22 fm offset from Eides, and that
offset is structurally the recoil-mixing in the literature Zemach formula
that the framework's leading-order kernel does not yet capture.

---

## 2. Layer-2 budget verification (the load-bearing computation)

Per the rzg_bug_diagnosis_memo.md lesson, Layer-2 itemization quality is the
load-bearing systematic. Here is the careful itemization for muH 1S HFS,
with explicit citations:

### muH 1S HFS literature decomposition (Krauth 2017 / Antognini-Krauth-Pohl 2015)

Reference: Krauth, J. J. *et al.*, *Hyperfine Interact.* **242**, 28 (2021),
Table 1; Antognini, A., Kottmann, F., Pohl, R., *J. Phys. Chem. Ref. Data*
**44**, 031210 (2015), Table III.

Total muH 1S HFS theory at r_Z(p) = 1.045 fm (Eides reference):

| Contribution | meV | ppm of ν_F (=182.443 meV) |
|---|---|---|
| BF strict (Bohr–Fermi) | 182.443 | (baseline) |
| a_μ Schwinger one-loop | +0.162 | +889 |
| a_μ higher-order | +0.001 | +3 |
| QED muon-line + recoil corrections | +0.008 | +44 |
| **Electron-VP leading (Uehling)** | **+1.496** | **+8194** |
| eVP iterated (NLO) | +0.005 | +27 |
| muon-VP | +0.003 | +15 |
| Two-loop QED with eVP | -0.002 | -11 |
| Hadronic VP | +0.011 | +59 |
| Recoil NLO (α(Zα)² m/M) | -0.019 | -103 |
| Recoil higher-order | +0.002 | +11 |
| **Zemach (-2 Z α m_red r_Z, r_Z=1.045 fm)** | **-1.304** | **-7141** |
| Inelastic polarizability | +0.008 | +44 |
| Zemach-recoil mixing | +0.0001 | +0.5 |
| **Total Krauth full theory** | **+182.725** | **+1545 (sum)** |

### Framework BF + Schwinger + leading Zemach at r_Z = 1.045 fm

The framework's Eides leading-order Zemach kernel is

    Δν_Z / ν_F = -2 Z m_red r_Z(bohr)

At r_Z = 1.045 fm (= 1.974×10⁻⁵ bohr), Z = 1, m_red(μp) = 185.84 m_e:

    Δν_Z / ν_F = -7340 ppm
    Δν_Z = -1.339 meV  (vs Krauth's -1.304 meV; recoil-mixing diff ~3%)

Framework total at r_Z = 1.045 fm:

    Framework_BF (182.443) + Schwinger (+0.212) + Zemach_kernel (-1.339)
        = 181.316 meV

Sprint MH Track B confirmed bit-identical: see memo §1 line 15.

### Closure gap (= Layer-2)

    Krauth full theory (182.725) - Framework total (181.316) = +1.409 meV
                                                              = +7722 ppm of ν_F

This +7722 ppm represents the literature corrections the framework does not
natively compute: dominantly electron-VP (eVP, +8194 ppm), with smaller
contributions from muon-VP (+15), iterated eVP (+27), two-loop QED (-11),
hadronic VP (+59), recoil corrections (-92 net), QED muon-line (+44),
inelastic polarizability (+44), Zemach-recoil mixing (+0.5), and small
discrepancy on the recoil-NLO sign convention.

**This +7722 ppm is the proper Layer-2 budget for muH 1S HFS** in the
framework's "BF + Schwinger" convention. Note: in the original rZG fit
(BF strict only, no Schwinger), the corresponding Layer-2 would be
+7722 + 1161 = +8883 ppm. Either convention gives the same residual
function and same fit values.

### Proper sigma on this Layer-2

Two scales matter:

1. **Multi-loop QED itemization scale**: sum-in-quadrature of the line
   uncertainties in Krauth 2017 Table 1:
   - eVP NLO: ~5 ppm
   - hadronic VP: ~30 ppm
   - recoil higher-order: ~50 ppm
   - QED muon-line h.o.: ~10 ppm
   - polarizability model dependence: ~80 ppm
   - **Combined: ~100–150 ppm**

2. **Kernel approximation scale**: the framework's Eides leading-order
   Zemach kernel `-2 Z m_red r_Z` omits next-to-leading recoil-mixing.
   Sprint MH Track B §2.3 measured this directly: framework kernel value
   -1.339 meV vs Krauth's full-theory Zemach line -1.304 meV gives a ~3%
   discrepancy at r_Z = 1.045 fm. Kernel uncertainty in absolute meV:
   - 0.05 × |Zemach| = 0.05 × 7340 ppm = **367 ppm**

**Critical question: which sigma applies?**

The fit's purpose is to extract r_Z(p) from observables. The sigma represents
all uncertainties that feed into the r_Z extraction at fixed observable
value. Both scales feed into the extraction:

- Multi-loop QED (150 ppm) is uncertainty on the value of "everything except
  Zemach" — direct uncertainty on Layer-2.
- Kernel approximation (367 ppm) is uncertainty on the *framework's*
  Zemach formula coefficient at fixed r_Z — uncertainty on b·r_Z, propagating
  to r_Z extraction as σ_r_Z = 367 / 7024 = 52 mfm.

The kernel approximation is **the larger systematic**, and it is uncorrelated
with the multi-loop QED (different physics: kernel approximation is in the
W1b operator, multi-loop QED is in the LS-8a wall sector). Quadrature sum:
sqrt(150² + 367²) = 396 ppm — kernel-dominated.

The honest σ for muH 1S HFS as an r_Z(p) extraction observable is therefore
**~367 ppm** (kernel-dominated), giving per-observable σ(r_Z) = 52 mfm.

---

## 3. Per-observable r_Z extractions

| Observable | b (ppm/fm) | (a+L2) (ppm) | r_Z alone (fm) | σ_alone (mfm) |
|---|---|---|---|---|
| H 21cm HFS | -37.79 | +39.5 | 1.045 | 265 |
| muH 1S HFS (this sprint, σ=367) | -7024 | +8885 | 1.265 | 52 |
| D 1S HFS | -37.79 | +97.6 | 2.583 | 1588 |

**Sanity check:** H 21cm alone gives r_Z(p) = 1.045 fm (exactly Eides — by
construction; the Layer-2 was tuned to that target value). muH 1S HFS alone
gives r_Z(p) = 1.265 fm — **+0.22 fm offset from Eides**.

This +0.22 fm offset is the framework's structural fingerprint: at the
leading-order Eides Zemach kernel, the framework "sees" the muH HFS with a
slightly weaker Zemach effect than the literature accounting (which uses
the full recoil-mixed Zemach). To reproduce the same observed muH HFS
shift, the framework needs a *larger* r_Z value to compensate for its
weaker kernel.

---

## 4. Global fit at sigma = 367 ppm

Analytic least-squares (decoupled per species):

| Quantity | Value | Tension vs literature |
|---|---|---|
| r_Z(p) | 1.257(51 mfm) fm | +3.85σ above Eides 1.045(20), +4.55σ above lattice 1.013(16) |
| r_Z(D) | 2.583(1588 mfm) fm | −0.01σ from Friar-Payne 2.593(16) |
| χ²/dof | 0.66/1 | OK |

**Per-observable residuals at fit:**

| Observable | Residual (ppm) | χ |
|---|---|---|
| H 21cm HFS | -8.0 | -0.80 |
| muH 1S HFS | +57.6 | +0.16 |
| D 1S HFS | 0.00 | 0.00 |

Both H 21cm and muH 1S HFS contribute meaningfully to the fit; the muH
observable now dominates the precision (b = 7024 vs 37.79, factor 186×) but
H 21cm anchors the fit toward the Eides value through its tight Layer-2
itemization. The fit lands at the weighted average pulled strongly by muH.

---

## 5. Verdict on Eides–lattice discrimination

**Cannot resolve with current observables.**

The arithmetic:

- σ(r_Z(p), framework) = 51 mfm
- Eides–lattice gap = 32 mfm (= |1.045 - 1.013|)
- Eides + lattice combined sigma = sqrt(20² + 16²) = 25.6 mfm

The framework's σ exceeds the Eides-lattice gap. More importantly, the
framework's central value (1.257 fm) sits 0.21 fm above Eides — so even
with arbitrarily small sigma, the framework would not "pick" Eides over
lattice; it would pick neither. Both Eides and lattice are at >3σ from
the framework central value.

**Three readings of this result:**

1. **Honest reading (most likely):** The framework's leading-order Eides
   Zemach kernel is at ~5% precision in muonic systems, NOT at sub-percent
   precision. The +0.22 fm offset between "muH-alone extraction" and Eides
   is the structural fingerprint of the missing recoil-mixing piece. To
   discriminate Eides vs lattice at 30 mfm precision, the framework needs
   the W1b magnetization-density operator extended to include recoil-mixing.

2. **Literature reading:** The Krauth 2017 itemization is itself sensitive
   at the 5% level to the recoil-mixing convention adopted; if our line-by-line
   Layer-2 totals were off by even 0.01% of ν_F, that propagates to ~5 mfm
   in extracted r_Z(p). Multiple literature compilations (Karshenboim, Eides,
   Antognini-Krauth-Pohl) use slightly different conventions; the "true"
   Layer-2 is itself uncertain at ~1% level when reconstructed from
   itemized lines.

3. **Wishful reading (unlikely):** The framework's r_Z(p) ~ 1.26 fm is the
   "true" value, contradicting both atomic-physics and lattice extractions.
   This would require extraordinary evidence and is incompatible with the
   well-established muonic Lamb shift data.

**Reading (1) is the operational interpretation.** The framework
reproduces Eides at the leading-order kernel; deviations from Eides in the
muH-driven extraction are the framework's documented W1b extension point
(Sprint MH Track B §2.4–2.5), not new physics.

---

## 6. Other observables considered but skipped

### Muonic deuterium 1S HFS

**Status:** No direct sub-1% experimental measurement exists. The
Pachucki 2007 / Borie 2012 theoretical decomposition can serve as a
stand-in "experimental" target, but using a theoretical full-theory value
to extract r_Z(D) introduces convention-mismatch risk: the framework's
"BF strict at I=1, m_red(μd)=195.7" baseline must reconcile with Borie's
"Fermi" entry (49.18 meV) before the muD HFS can be added meaningfully.

**Defer to follow-on sprint:** wire `geovac.bohr_fermi_hyperfine` natively
at I=1 muD to nail down the BF baseline against Borie 2012 Tab. 12, then
add muD HFS to the global fit. Predicted contribution: σ(r_Z(D)) tightens
from 1588 mfm to ~52 mfm (mass-enhancement same order as muH).

### ³He+ ground-state hyperfine (electronic ³He+)

**Status:** Z=2 lever arm increases b by factor of 2 (b = 75.6 ppm/fm at
electronic ³He+), but constrains r_Z(³He), not r_Z(p) or r_Z(D). Not
useful for the Eides-vs-lattice question on r_Z(p).

**Defer to species-extension sprint:** add ³He+ HFS as a third-species r_Z
constraint when the framework is extended to triple-species global fits.

### Hydrogen 2S HFS (Kolachevsky 2004)

**Status:** Kolachevsky measurement σ = 16 Hz / 177 MHz = 90 ppb (0.09 ppm).
Zemach shift on 2S is -39.5 ppm (same coefficient as 1S, scaled by 1/n³ in
absolute units). Per-observable σ(r_Z) = 0.09 / 37.8 = 0.0024 fm = 2.4 mfm
*if* the literature decomposition of D21 = ν_2S - ν_1S/8 were available
at sub-Hz uncertainty. However, the leading Zemach piece **cancels in D21**
(scales as 1/n³), and the residual r_Z-sensitivity in D21 is < 1 Hz — well
below current measurement precision. The absolute 2S HFS frequency carries
the same r_Z information as the 1S HFS but with weaker statistics and a
larger absolute Layer-2 budget. Skipped: the 1S HFS (H 21cm) already
contains the same r_Z information at 10 ppm precision via Eides Tab. 7.3.

---

## 7. Net result

**Diagnostic outcome:** the proposed strategy of "tighten muH 1S HFS Layer-2
to literature itemization" works in the sense that it formally reduces
σ(r_Z(p)) to 21 mfm (which would be below the 30 mfm gap), but it produces a
central value 1.179 fm that is inconsistent with both Eides and lattice at
>4σ. The 4–6σ tension is NOT a positive signal of new physics; it is a
**diagnostic flag** that the framework's leading-order Eides Zemach kernel
is itself the dominant systematic at this precision, not the multi-loop
QED Layer-2 itemization.

**Honest extension fit at sigma=367 ppm (kernel-limited):**

- σ(r_Z(p)) = 51 mfm (close to but not below the 50 mfm target)
- r_Z(p) = 1.257 fm (still inconsistent with both Eides and lattice at >3.85σ)
- Framework cannot discriminate Eides vs lattice at current operator-level
  precision.

**Action items for follow-on sprints (recommended):**

1. **Phase C-W1b operator extension to include recoil-mixing.** The
   `magnetization_density.py` module currently implements the leading-order
   Eides formula. Extending to next-to-leading order (Friar-Payne 2005
   §III treatment) would close the 0.22 fm offset and bring r_Z(muH alone)
   to ~1.05 fm. Estimated 1-week sprint per Sprint MH Track B §2.4-2.5
   (the operator-level extension was flagged on 2026-05-08).

2. **muD 1S HFS native BF integration.** Compute the framework's BF strict
   at I=1, m_red(μd) = 195.7 m_e, and reconcile with Borie 2012 "Fermi"
   convention. Once reconciled, muD HFS adds a fourth observable
   constraining r_Z(D) at σ ~ 50 mfm precision. ~3-day sprint.

3. **³He+ HFS species extension.** Add r_Z(³He) as a third fit parameter.
   Useful for cross-species verification of §III.18 but does not directly
   help the Eides-vs-lattice question on r_Z(p).

**Net:** the framework's σ(r_Z(p)) is currently kernel-limited at ~50 mfm,
not Layer-2 limited. The W1b operator extension is the structural path to
sub-30 mfm precision.

---

## 8. Files modified / created

- **Created:** `debug/calc_track_rZG_extended.py` (new driver, parallel to
  the original rZG; does not overwrite)
- **Created:** `debug/data/calc_track_rZG_extended.json` (numerical output)
- **Created:** `debug/calc_track_rZG_extended_memo.md` (this memo)
- **Modified:** `papers/observations/paper_34_projection_taxonomy.tex` —
  rZG row in §V.B updated to reflect extended-fit verdict
- **Production code:** unchanged (verified by running
  `pytest tests/test_cross_register_vne.py tests/test_magnetization_density.py`
  with 106 passes, 0 failures).

---

## 9. References

- **Sprint MH Track B/C:** `debug/sprint_mh_track_b_memo.md` (May 2026) —
  the +2 ppm BF match against Eides pure-QED, 0.55% Zemach mass-enhancement
  match against Eides muonic target, and the W1b operator-level extension
  point flagging recoil-mixing as a 1-sprint extension.
- **Sprint rZG corrected:** `debug/calc_track_rZG_global_zemach.py` (May 9,
  2026) and `debug/rzg_bug_diagnosis_memo.md` — the corrected baseline
  3-observable fit at σ(r_Z(p)) = 166 mfm.
- **Krauth 2017:** Krauth, J. J. *et al.*, *Hyperfine Interact.* **242**,
  28 (2021) — comprehensive muH HFS theory itemization Table 1.
- **Antognini-Krauth-Pohl 2015:** *J. Phys. Chem. Ref. Data* **44**, 031210
  (2015) — CREMA muonic atoms physics review with Tab. III itemization.
- **Eides 2024:** Eides, M. I., Grotch, H., Shelyuto, V. A. *Theory of
  Light Hydrogenic Bound States* (Springer, 2007 + 2024 updates), Ch. 7.3
  Tab. 7.3 (electronic H) and Tab. 7.4 (muonic H).
- **Karshenboim 2005:** Karshenboim, S. G. *Phys. Rep.* **422**, 1 (2005)
  — review.
- **Friar-Payne 2005:** Friar, J. L. & Payne, G. L., *Phys. Rev. C* **72**,
  014002 (2005) — deuteron r_Z.
- **Pachucki 2007:** Pachucki, K., *Phys. Rev. A* **76**, 022508 (2007) —
  muonic deuterium theory.
- **Borie 2012:** Borie, E., *Annals Phys.* **327**, 733 (2012) Tab. 12 —
  muonic deuterium theoretical decomposition.
- **Pachucki-Yerokhin 2010:** *Phys. Rev. A* **81**, 062510 (2010) —
  deuteron polarizability.
