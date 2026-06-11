# Precision Catalogue: Positronium 1S Hyperfine Splitting

**Date:** 2026-05-08.
**Driver:** `debug/precision_catalogue_positronium.py`
**Data:** `debug/data/precision_catalogue_positronium.json`
**Verdict:** POSITIVE STRUCTURAL — framework reproduces Fermi-contact part of Ps HFS exactly (4/7 of LO); annihilation channel and multi-loop QED are cleanly LS-8a-attributable Layer-2 inputs.

## Goal

Test multi-focal architecture in the equal-mass limit lambda_a = lambda_b
(both leptons have m = m_e, m_red = m_e/2). Identify how the Roothaan
multipole machinery behaves when there is no asymptotic acceleration, and
quantify the framework's native vs Layer-2 split for the 1S HFS.

## Architecture

Positronium = e⁻ e⁺. The simplest two-particle bound state in QED.

| Quantity | Value |
|:---------|:------|
| Lepton register | electron (m_e, g_e Dirac=2) |
| "Nucleus" register | positron (m_e, g_e Dirac=2) — equal mass, also lepton |
| m_red(ee) | 0.5 m_e |
| Multi-focal regime | λ_lepton = λ_positron = m_red = 0.5; **NO asymptotic small parameter** |
| Experimental | 203,389,100(740) kHz = 203.3891 GHz (Ishida 2014) |

## Decomposition

The leading-order Ps HFS has two structurally distinct contributions
(Karshenboim 2005 review §4):

| Component | Magnitude | LO fraction | Framework status |
|:----------|:----------|:------------|:-----------------|
| Fermi contact (magnetic dipole-dipole) | (1/3) m_e c² α⁴ /h | 4/12 | **NATIVE** |
| Annihilation (e⁺e⁻ → γ → e⁺e⁻) | (1/4) m_e c² α⁴ /h | 3/12 | Layer-2 (LS-8a) |
| **Total LO** | **(7/12) m_e c² α⁴ /h** | **7/12** | mixed |

Higher-order (α^5 m_e c² / h and beyond) corrections from
Czarnecki-Melnikov-Yelkhovsky 2000 and follow-ups close the gap from LO total
(204.39 GHz) to experimental (203.389 GHz) at ppm level.

## Component results

### Step 1: Framework-native Fermi contact via Bohr-Fermi

ν_F(Ps, Dirac g=2) = (2/3) × g_e² × α² × (m_red(ee))³ / (m_e × m_e) × Hartree(m_e)
                    = (2/3) × 4 × α² × (1/8) × Ha
                    = (1/3) α² × Ha
                    = **116.7924 GHz**

Exact match to canonical (1/3) m_e c² α⁴ /h to machine precision.

With full g_e = 2.00231930... (Schwinger included): ν_F = 117.0634 GHz
(0.232% enhancement vs Dirac, as expected for one-loop a_e).

### Step 2: Annihilation channel (Layer-2)

ν_annih = (1/4) m_e c² α⁴ /h = **87.5943 GHz**.

Mechanism: virtual e⁺e⁻ → γ → e⁺e⁻ with both photon and pair on-shell. This
requires a vertex coupling e⁺e⁻ to a *vector* photon AND the photon back to
e⁺e⁻ — a genuinely second-quantized field-theory effect that the framework's
bare action on the graph does not generate. **Same wall as Källén-Sabry
two-loop VP for muonic H.** Sits in the LS-8a-class wall (virtual photon
loops, not native to the bare graph).

### Step 3: Total LO

(7/12) m_e c² α⁴ /h = **204.3866 GHz**.

vs experimental 203.3891 GHz → residual **+998 MHz / +0.490%**.

The 0.49% residual is multi-loop QED, fully external (Czarnecki et al. 2000).

### Step 4: Multi-focal at equal masses

For Ps, λ_lepton = λ_positron = m_red(ee) = 0.5. The Roothaan multipole
expansion 1/r₁₂ = Σ_L (r_<^L / r_>^(L+1)) P_L(cos θ_12) **still terminates at
L_max = 2 l_max** by Gaunt selection rules — no structural issue. There IS
no asymptotic acceleration (no "small parameter" between the two scales),
but the integral converges at every L.

This is structurally distinct from the muonic-hydrogen case where λ_lepton ≪
λ_nucleus and the heavy-nucleus expansion accelerates. **For Ps, the
cross-register kernel reduces to the standard symmetric two-body Coulomb
integral — exactly the case Roothaan 1951 originally handled.**

The multi-focal machinery extends cleanly to λ_a = λ_b without modification;
it just operates without asymptotic compression.

## Verdict and catalogue entry

**POSITIVE STRUCTURAL.** The framework cleanly partitions Ps HFS:

- Framework-native Fermi-contact:  4/7 of LO (116.79 GHz, machine-precision match)
- Annihilation channel:             3/7 of LO (87.59 GHz, Layer-2 input)
- Multi-loop QED:                  ~0.5% of LO (Layer-2 input)

For Paper 34 §V.B (off-precision matches with error code A — approximation order):

| System | Observable | Framework value | Reference | Residual | Code |
|:-------|:-----------|:----------------|:----------|:---------|:----:|
| Ps | 1S HFS (Fermi-contact only) | 116.79 GHz | 203.389 GHz (Ishida 2014) | −42.58% | A |
| Ps | 1S HFS (framework + Layer-2 LO annihilation) | 204.39 GHz | 203.389 GHz | +0.49% | A |

Both rows are off-precision; error code A (approximation order, framework
provides Fermi-contact native; annihilation + multi-loop are next orders that
the framework's bare action does not generate).

## New finding for the multi-focal architecture

**The Roothaan multipole termination at L_max = 2 l_max is preserved at the
equal-mass limit lambda_a = lambda_b.** This was an open question raised in the
muonic-hydrogen Track A memo (Sprint MH): does the multipole convergence
machinery break when the small-parameter expansion is absent? Answer: No.
The termination is a Gaunt selection rule (angular content), not a
small-parameter expansion. It holds independent of the ratio of focal
lengths.

This is worth noting in CLAUDE.md §2 multi-focal sprint outcome: the
architecture covers the full hierarchy from λ_lepton ≪ λ_nucleus
(hydrogen-like, asymptotic acceleration available) through λ_lepton = λ_nucleus
(positronium, symmetric, no acceleration but Gaunt termination preserved) to
λ_lepton ≫ λ_nucleus (would arise in inverted-mass scenarios; not yet tested
in the catalogue, but no obstruction expected by the same argument).

## Connection to other multi-focal observables

| System | λ_l/λ_n | Framework native? | Off-precision residual |
|:-------|:--------|:------------------|:-----------------------|
| H 1S-2S | ~1/1836 | Yes (Bohr+SE+recoil) | ~10⁻¹¹ (Parthey ref) |
| H 21cm (HFS) | ~1/1836 | Bohr-Fermi yes | +18 ppm (multi-loop) |
| μH 2S-2P (Lamb) | ~1/9.7 | Yes (full Uehling kernel) | −0.10% |
| μH 1S HFS | ~1/9.7 | Yes (Bohr-Fermi + Zemach) | +2 ppm |
| Mu 1S-2S | ~1/207 | Yes (rest-mass + recoil) | −0.25 ppm |
| **Ps 1S HFS** | **1/1 (this)** | **Fermi-contact only** | **+0.49% (multi-loop+annih)** |

Ps fills in the equal-mass corner of the matrix. Together with Mu (electron
on antimuon), μH (muon on proton), and H (electron on proton), the multi-focal
architecture has been validated across four systems spanning the e-μ-p mass
hierarchy.
