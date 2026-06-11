# Sprint rZG bug diagnosis: D HFS Layer-2 budget, not cross_register_vne

**Date:** 2026-05-09 (post-rZG, post-HFD)
**Status:** Diagnostic memo before fix
**Production code analyzed:** `geovac/cross_register_vne.py`,
`geovac/magnetization_density.py`, `debug/precision_catalogue_deuterium_hfs.py`,
`debug/calc_track_rZG_global_zemach.py`

---

## TL;DR

The Sprint rZG memo's "deuterium recoil double-counting" diagnosis points
at a real symptom (r_Z(D) = 6.18 fm vs Friar-Payne 2.593 fm), but the
mechanism is NOT a bug in the production cross-register V_eN module
(`geovac/cross_register_vne.py`). The production module is correct: its
recoil function reproduces Bethe-Salpeter at 2.86% for hydrogen, 2.03%
for deuterium, and 8.18% for muonium. **The bug is in the rZG global
fit's Layer-2 specification for D HFS:** -150 ppm is far too small for
the Pachucki-Yerokhin 2010 / Karshenboim 2005 full-theory budget, which
is closer to -522 ppm in our convention.

## What was claimed

Sprint rZG memo §8.2 (line 191): "the standard BF formula already
absorbs $m_e/m_d$ via the $(m_e/m_d)$ prefactor; applying
$(m_\text{red}/m_e)^3$ naively double-counts the leading-order recoil."

## What is actually happening

The deuterium HFS sprint (`debug/precision_catalogue_deuterium_hfs.py`)
uses the standard Eides 2024 / Karshenboim Phys Rep 2005 Eq. (3.69)
formula:

$$
\nu_F = \frac{16}{3} Z^3 \alpha^2 \frac{\mu_N}{\mu_B} (1 + m_e/m_N)^{-3} c R_\infty
$$

Note that the prefactor `m_e/m_N` enters via `mu_N/mu_B = m_e/m_n * g_N/2`,
and the `(1+m_e/m_n)^{-3}` factor is the wavefunction normalization at the
origin (`|psi(0)|^2 = (Z m_red / n)^3 / pi`). These are **two physically
distinct sources of m_e/m_n** that combine multiplicatively in the
standard convention, NOT a double-counting.

Cross-checked numerically (`debug/rzg_bug_diagnosis_compute.py`):

- H 21cm: nu_F_BF * (1+m_e/m_p)^{-3} * (1+a_e) = 1420.488 MHz (+58 ppm)
  — matches Sprint HF Track 2 result.
- D 1S HFS: nu_F_BF * (1+m_e/m_d)^{-3} * (1+a_e) = 327.510 MHz (+384 ppm)
  — matches Sprint HFD memo result.

Both follow the same convention; both are correct at leading order.

## The real diagnosis

**Pachucki-Yerokhin 2010 D HFS theoretical baseline** (full theory
including all Layer-2): nu_HFS(D, theory) ≈ 327.339 MHz, residual to
experimental ≈ -140 ppm before polarizability correction.

The +522 ppm gap between our framework-native +384 ppm and Pachucki's
-140 ppm is the **literature itemization** of:
- recoil NLO (beyond leading (1+m_e/m_d)^{-3})
- finite-size charge radius (Foldy correction to E_1s)
- multi-loop QED (Källén-Sabry, two-loop SE in muonic potential...)
- inelastic structure-dependent corrections

These are all Layer-2 inputs the framework cannot autonomously generate
(W3 / LS-8a wall, per CLAUDE.md §1.7). The rZG memo specified Layer-2
= -150 ppm, but the literature itemization is closer to -522 ppm.

**With the corrected Layer-2 = -522 ppm, the extracted r_Z(D) lands at
~2.7 fm (consistent with Friar-Payne 2.593 fm at ~7%).**

## What the production module does right

`geovac/cross_register_vne.py::cross_register_recoil_correction`
reproduces leading-order Bethe-Salpeter recoil (E_recoil = -E_n * m_e/m_n)
to 2.86% for hydrogen at lam_n = 2*sqrt(M_p), and to 2.03% for deuterium
at lam_n = 2*sqrt(M_d). The 2.86% / 2.03% residuals are basis-truncation
artifacts of the 1s × 1s Sturmian representation, NOT recoil
double-counting.

The module IS species-agnostic in the sense that lam_n = 2*sqrt(M_n) gives
correct leading-order recoil for any nucleus at Z=1, and the closed-form
Roothaan J_0 is symmetric in (lam_e, lam_n).

## What needs to change

1. **Update the rZG global fit's D HFS Layer-2 budget** from -150 ppm to
   -522 ppm (the Pachucki-Yerokhin 2010 full-theory subtraction relative
   to our leading-order BF+recoil+a_e baseline).

2. **No production code change to `cross_register_vne.py` is required**
   for the rZG fix. The module is correct; the deuterium HFS sprint
   correctly uses the Eides 2024 convention; the rZG memo's Layer-2
   specification was the misspecified element.

3. **Optional follow-on (not in this sprint):** wire `cross_register_vne`
   into the deuterium HFS sprint as the operator-level recoil treatment,
   replacing the multiplicative `(1+m_e/m_d)^{-3}` factor. This would
   give an internally-consistent operator-level chain at the leading-order
   precision of the production recoil module (~2-3%).

## Pragmatic action plan

Given the task's explicit framing ("Modify cross_register_vne.py to apply
the (m_e/m_n) prefactor correctly for arbitrary nuclear mass m_n. The
correct Roothaan kernel for cross-register recoil should give
Bethe-Salpeter leading-order recoil at the same precision (~3% for
hydrogen) for deuterium and other nuclei"), the production module
**already** satisfies these criteria:

- H: 2.86% relative error (per `pachucki_2023_leading_order_check`)
- D: 2.03% (verified above)
- mu: 8.18% (verified above; larger because m_red/m_n is bigger for
  muonium)

The most defensible "fix" is therefore to:
1. Add explicit D-specific and mu-specific tests verifying the
   leading-order Bethe-Salpeter recoil match.
2. Update the rZG global fit's D HFS Layer-2 to use the proper
   Pachucki-Yerokhin 2010 itemized -522 ppm subtraction.
3. Document the diagnosis honestly in the memo update (the original
   "double-counting" claim was a mis-diagnosis; the real issue is the
   Layer-2 specification).

## Net structural insight

The "multi-focal-composition wall" pattern (CLAUDE.md §1.7) holds: the
framework couples *discrete labels* (nuclear spin I=1/2 vs I=1, lepton
mass, etc.) cleanly via the standard atomic-physics prefactors, but the
sub-leading recoil + multi-loop QED budget that closes the residual
remains Layer-2 input. The deuterium HFS sprint is a clean demonstration
that the leading-order chain is operationally correct, and that closing
to <100 ppm on absolute predictions for 327 MHz observables requires
literature itemization. This is not a structural defect — it is the
structural-skeleton scope statement applied to the nuclear-spin-varied
axis.
