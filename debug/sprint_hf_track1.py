"""
Sprint HF Track 1 — Bohr-Fermi A_hf for hydrogen 21 cm
========================================================

Goal
----
Derive A_hf for hydrogen (1s ground state) at the Bohr-Fermi level
(g_e = 2 exactly, no anomalous moment, no recoil, no nuclear structure,
no multi-loop QED) from GeoVac's native machinery, and compare to the
experimental 21 cm line 1 420 405 751.768 Hz.

Formula (atomic units, Hartree)
-------------------------------
The Fermi-contact hyperfine coupling for hydrogen 1s with point-particle
Dirac proton, Pauli electron (g_e = 2), no recoil, in atomic units, is

    H_F = (8 pi / 3) * g_e * g_p * mu_B * mu_N * (S_e . I_p) * delta^3(r)

so that the effective coupling A_hf in the convention H_eff = A_hf I.S is

    A_hf = (8 pi / 3) * g_e * g_p * mu_B * mu_N * |psi_1s(0)|^2

with
    mu_B = alpha / 2          (Bohr magneton, atomic units; e h_bar / 2 m_e c)
    mu_N = (m_e/m_p) * alpha / 2   (nuclear magneton, atomic units)
    |psi_1s(0)|^2 = Z^3 / pi  (hydrogenic 1s wavefunction at origin)

Substituting and simplifying:

    A_hf = (8 pi / 3) * g_e * g_p * (alpha/2) * (m_e/m_p)(alpha/2) * (Z^3/pi)
         = (2/3) * g_e * g_p * alpha^2 * (m_e/m_p) * Z^3   Hartree.

For hydrogen (Z = 1), with g_e = 2 (Dirac tree-level):

    A_hf = (4/3) * g_p * alpha^2 * (m_e/m_p)   Hartree.

This is the canonical Bohr-Fermi result. (Numerical: ~ 2.160e-7 Ha
~ 1421.2 MHz at strict CODATA inputs, no recoil.)

Note on conventions
-------------------
The PI briefing wrote A_hf = (4/3) alpha^2 g_e (g_p/2) (m_e/m_p) Hartree,
which is *literally* a factor of 2 short of the canonical formula above
when interpreted as written.  In context, "(g_p/2)" was almost certainly
shorthand for mu_p/mu_N ~ 2.793 (the proton's magnetic moment in nuclear
magnetons), and the briefing's formula was meant to be read with an
implicit factor of 2 inside the magneton definition (i.e., mu_p = (g_p/2)
* 2 mu_N = g_p mu_N for spin I = 1/2).  The Fermi-contact Hamiltonian
written in terms of operators carries g_e g_p, not g_e (g_p/2), and the
value 1418.84 MHz quoted in the briefing's expected outcome corresponds
to the g_e g_p form (with subsequent reduced-mass correction).  We use
the canonical g_e g_p form here.

External focal lengths (calibration data, supplied here)
--------------------------------------------------------
- alpha = 7.2973525693e-3  (CODATA 2018)
- g_p = 5.585694689        (CODATA 2018, proton g-factor)
- m_e/m_p = 1/1836.15267343 (CODATA 2018, electron-proton mass ratio)
- 1 Ha = 6.579683920502e15 Hz   (Hartree-frequency conversion; only used to
                                  report a final number in MHz)

Forbidden inputs (none used)
----------------------------
- A_hf experimental value (only used in the final residual calculation,
  and only there).
- Numerical |psi_1s(0)|^2 supplied as a constant (we obtain it
  symbolically from the GeoVac hydrogenic 1s wavefunction, then
  evaluate at r = 0).

Internal sources (where each focal length comes from in GeoVac)
---------------------------------------------------------------
- |psi_1s(0)|^2 = Z^3/pi: derived from
  geovac.dirac_matrix_elements._hydrogenic_radial_wavefunction(n=1, l=0, r=0, Z=1)
  evaluated symbolically in sympy.
- g_e = 2: this is the Dirac-Pauli value; on S^3 the Camporesi-Higuchi
  Dirac operator gives the same tree-level g_e = 2, but this is a
  property of the Dirac equation, not GeoVac-specific. We flag this
  honestly.
- (4/3) prefactor: assembled from the (8 pi / 3) Fermi contact prefactor
  times (1/pi) from |psi(0)|^2 and the (1/4) from the two magneton
  factors of (alpha/2). All algebraic.

What this sprint can and cannot claim
-------------------------------------
- We CAN claim: "Given calibration of (alpha, g_p, m_e/m_p), GeoVac's
  algebraic hydrogen 1s machinery (Fock projection / Dirac on S^3
  framework) reproduces A_hf at the Bohr-Fermi level, residual ~ 1100
  ppm against experiment." This residual is the well-known a_e ~
  alpha / (2 pi) anomalous moment correction (Track HF-2 will address).
- We CANNOT claim: "GeoVac derives g_e = 2." The Dirac equation gives
  g_e = 2 at tree level; this is independent of the S^3 projection.
- The 1s wavefunction at the origin |psi_1s(0)|^2 = Z^3/pi is derived
  algebraically from the *continuum* hydrogenic wavefunction, not from
  the discrete graph. The Fock projection makes the n-shell graph
  Laplacian's spectrum exact (n^2 - 1), but the spatial profile of the
  1s wavefunction is the continuum embedding — see the algebraic
  registry, CLAUDE.md Section 12.

Outputs
-------
- debug/sprint_hf_track1.py    (this file)
- debug/data/sprint_hf_track1.json   (numerical results, structured)
- debug/sprint_hf_track1_memo.md     (~1500 word narrative)
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from typing import Any, Dict

import sympy as sp
from sympy import Rational, Symbol, sqrt, simplify, exp, Integer, pi as sp_pi, oo

# Make geovac importable when the script is run directly from the project root.
_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from geovac.dirac_matrix_elements import _hydrogenic_radial_wavefunction


# ---------------------------------------------------------------------------
# 1. CODATA 2018 calibration constants (external focal lengths)
# ---------------------------------------------------------------------------

# Fine-structure constant (dimensionless, CODATA 2018)
ALPHA_CODATA: float = 7.2973525693e-3

# Proton g-factor (dimensionless, CODATA 2018)
G_P_CODATA: float = 5.585694689

# Electron-proton mass ratio (dimensionless, CODATA 2018)
ME_OVER_MP_CODATA: float = 1.0 / 1836.15267343

# Hartree-to-Hz conversion (CODATA 2018; used only for final unit display)
HZ_PER_HARTREE: float = 6.579683920502e15

# Experimental 21 cm transition frequency (CODATA / NIST)
NU_HF_EXPERIMENTAL_HZ: float = 1.420405751768e9
NU_HF_EXPERIMENTAL_MHZ: float = NU_HF_EXPERIMENTAL_HZ / 1.0e6

# Reduced-mass correction factor (1 + m_e/m_p)^{-3}.
# This is a *recoil* effect: if both proton and electron are point particles
# but the proton has finite mass, the wavefunction localizes on the reduced-
# mass coordinate, contracting |psi(0)|^2 by (1 + m_e/m_p)^{-3}. The PI
# briefing says "no recoil", so this factor is computed but NOT applied to
# the headline Bohr-Fermi value. It is reported separately as a known
# correction.
REDUCED_MASS_FACTOR: float = (1.0 + ME_OVER_MP_CODATA) ** (-3)


# ---------------------------------------------------------------------------
# 2. The 1s wavefunction at origin from GeoVac's hydrogenic radial machinery
# ---------------------------------------------------------------------------

def psi_1s_squared_at_origin_symbolic() -> sp.Expr:
    """|psi_1s(r=0)|^2 evaluated symbolically.

    The full 1s wavefunction is psi_1s(r) = R_{1,0}(r) * Y_{0,0}(theta,phi),
    with Y_{0,0} = 1 / sqrt(4 pi). The radial function R_{1,0} comes from
    geovac.dirac_matrix_elements._hydrogenic_radial_wavefunction.

    Evaluating at r = 0:
        R_{1,0}(0)   = 2 * Z^{3/2}        (a textbook value, recovered here)
        Y_{0,0}(.)^2 = 1 / (4 pi)
    so
        |psi_1s(0)|^2 = R_{1,0}(0)^2 * (1/(4 pi)) = 4 Z^3 / (4 pi) = Z^3 / pi.

    Returns
    -------
    sympy.Expr
        Z^3 / pi, exact.
    """
    Z = sp.Symbol("Z", positive=True)
    r = sp.Symbol("r", nonnegative=True)
    R10 = _hydrogenic_radial_wavefunction(n=1, l=0, r=r, Z=Z)
    R10_at_0 = sp.simplify(R10.subs(r, 0))
    Y00_squared = Rational(1, 1) / (4 * sp_pi)  # |Y_{0,0}|^2 = 1/(4 pi)
    psi_sq = sp.simplify(R10_at_0 ** 2 * Y00_squared)
    return psi_sq


# ---------------------------------------------------------------------------
# 3. Symbolic Bohr-Fermi formula assembly
# ---------------------------------------------------------------------------

def bohr_fermi_symbolic() -> Dict[str, sp.Expr]:
    """Assemble A_hf at the Bohr-Fermi level symbolically.

    The Fermi-contact Hamiltonian (Bethe-Salpeter Eq. 22.18, atomic units):

        H_F = (8 pi / 3) * g_e * g_p * mu_B * mu_N * (S_e . I_p) * delta^3(r)

    so the coefficient of I.S after taking <delta^3(r)> = |psi(0)|^2 is

        A_hf = (8 pi / 3) * g_e * g_p * mu_B * mu_N * |psi_1s(0)|^2

    with mu_B = alpha/2, mu_N = (m_e/m_p) alpha/2 in atomic units.

    Substituting:

        A_hf = (8 pi / 3) * g_e * g_p * (alpha/2) * (m_e/m_p)(alpha/2)
                                                       * (Z^3/pi)
             = (8 / 3) * g_e * g_p * (alpha^2/4) * (m_e/m_p) * Z^3
             = (2/3) * g_e * g_p * alpha^2 * (m_e/m_p) * Z^3.

    For Z=1, g_e=2 the canonical hydrogen result is

        A_hf = (4/3) * g_p * alpha^2 * (m_e/m_p)   Hartree.
    """
    alpha_s = sp.Symbol("alpha", positive=True)
    g_e_s = sp.Symbol("g_e", positive=True)
    g_p_s = sp.Symbol("g_p", positive=True)
    eta_s = sp.Symbol("m_e_over_m_p", positive=True)
    Z_s = sp.Symbol("Z", positive=True)

    mu_B = alpha_s / 2
    mu_N = eta_s * alpha_s / 2

    psi0_sq = Z_s ** 3 / sp_pi
    fermi_prefactor = (8 * sp_pi) / 3

    # Canonical Fermi contact: (8 pi/3) g_e g_p mu_B mu_N |psi(0)|^2.
    A_hf = sp.simplify(
        fermi_prefactor * g_e_s * g_p_s * mu_B * mu_N * psi0_sq
    )
    # Simplification check.
    A_hf_canonical = sp.simplify(
        Rational(2, 3) * g_e_s * g_p_s * alpha_s ** 2 * eta_s * Z_s ** 3
    )
    assert sp.simplify(A_hf - A_hf_canonical) == 0, (
        f"Symbolic mismatch: {sp.simplify(A_hf)} vs {A_hf_canonical}"
    )

    # Specialize to Z = 1 and g_e = 2 (Bohr-Fermi).
    A_hf_BF = sp.simplify(A_hf.subs({Z_s: 1, g_e_s: 2}))
    # Should equal (4/3) * g_p * alpha^2 * (m_e/m_p)
    A_hf_BF_simple = sp.simplify(
        Rational(4, 3) * g_p_s * alpha_s ** 2 * eta_s
    )
    assert sp.simplify(A_hf_BF - A_hf_BF_simple) == 0, (
        f"Z=1, g_e=2 specialization mismatch: {A_hf_BF} vs {A_hf_BF_simple}"
    )

    return {
        "alpha": alpha_s,
        "g_e": g_e_s,
        "g_p": g_p_s,
        "m_e_over_m_p": eta_s,
        "Z": Z_s,
        "mu_B_au": mu_B,
        "mu_N_au": mu_N,
        "psi_1s_squared_at_origin": psi0_sq,
        "fermi_prefactor": fermi_prefactor,
        "A_hf_full": A_hf,
        "A_hf_BF_Z1_ge2": A_hf_BF,
        "A_hf_BF_canonical_form": A_hf_BF_simple,
    }


def derive_psi_1s_at_origin_from_geovac() -> Dict[str, Any]:
    """Independently re-derive |psi_1s(0)|^2 = Z^3/pi from GeoVac's machinery.

    Returns the intermediate values (R_{1,0}(0), Y_{0,0}^2, product) so that
    the chain of identities is logged in the JSON output.
    """
    Z = sp.Symbol("Z", positive=True)
    r = sp.Symbol("r", nonnegative=True)
    R10 = _hydrogenic_radial_wavefunction(n=1, l=0, r=r, Z=Z)
    R10_at_0 = sp.simplify(R10.subs(r, 0))
    expected_R10 = 2 * Z ** Rational(3, 2)
    assert sp.simplify(R10_at_0 - expected_R10) == 0, (
        f"R_{{1,0}}(0) mismatch: GeoVac gives {R10_at_0}, expected {expected_R10}"
    )
    Y00_squared = Rational(1, 1) / (4 * sp_pi)
    psi_sq = sp.simplify(R10_at_0 ** 2 * Y00_squared)
    expected_psi_sq = Z ** 3 / sp_pi
    assert sp.simplify(psi_sq - expected_psi_sq) == 0, (
        f"|psi_1s(0)|^2 mismatch: GeoVac gives {psi_sq}, expected {expected_psi_sq}"
    )
    return {
        "R_10_at_0_symbolic": str(R10_at_0),
        "Y_00_squared": str(Y00_squared),
        "psi_1s_squared_at_origin": str(psi_sq),
        "Z_equals_1_value": str(psi_sq.subs(Z, 1)),
    }


# ---------------------------------------------------------------------------
# 4. Numerical evaluation
# ---------------------------------------------------------------------------

def evaluate_numeric() -> Dict[str, Any]:
    """Evaluate the Bohr-Fermi A_hf numerically, log per-factor decomposition.

    Returns a dict suitable for JSON serialization.
    """
    sym = bohr_fermi_symbolic()
    subs = {
        sym["alpha"]: ALPHA_CODATA,
        sym["g_e"]: 2.0,
        sym["g_p"]: G_P_CODATA,
        sym["m_e_over_m_p"]: ME_OVER_MP_CODATA,
        sym["Z"]: 1,
    }

    # Per-factor numerics (all dimensionless or in atomic units).
    alpha_sq = ALPHA_CODATA ** 2
    g_e = 2.0
    g_p = G_P_CODATA
    g_p_half = G_P_CODATA / 2.0  # = mu_p / mu_N (proton mag moment in nuclear magnetons)
    eta = ME_OVER_MP_CODATA
    psi0_sq_au = 1.0 / 3.141592653589793  # |psi_1s(0)|^2 = 1/pi at Z=1, in 1/a_0^3

    # A_hf in Hartree
    A_hf_Ha_strict = float(sym["A_hf_BF_Z1_ge2"].subs(subs))
    # Cross-check via the explicit canonical formula:
    # A_hf = (4/3) g_p alpha^2 (m_e/m_p) at Z=1, g_e=2
    A_hf_Ha_check = (4.0 / 3.0) * g_p * alpha_sq * eta
    assert abs(A_hf_Ha_strict - A_hf_Ha_check) < 1e-22, (
        f"Numerical inconsistency: sympy {A_hf_Ha_strict} vs explicit {A_hf_Ha_check}"
    )
    # Second cross-check: the (8 pi /3) g_e g_p mu_B mu_N |psi(0)|^2 form.
    mu_B = ALPHA_CODATA / 2.0
    mu_N = eta * ALPHA_CODATA / 2.0
    A_hf_Ha_check2 = ((8.0 * 3.141592653589793) / 3.0) * g_e * g_p * mu_B * mu_N * psi0_sq_au
    assert abs(A_hf_Ha_strict - A_hf_Ha_check2) < 1e-22, (
        f"Cross-check 2 failed: sympy {A_hf_Ha_strict} vs (8pi/3) form {A_hf_Ha_check2}"
    )

    # In Hz / MHz
    A_hf_Hz_strict = A_hf_Ha_strict * HZ_PER_HARTREE
    A_hf_MHz_strict = A_hf_Hz_strict / 1.0e6

    # Optional: with reduced-mass correction, for context
    A_hf_Hz_with_recoil = A_hf_Hz_strict * REDUCED_MASS_FACTOR
    A_hf_MHz_with_recoil = A_hf_Hz_with_recoil / 1.0e6

    # Residuals
    residual_Hz_strict = A_hf_Hz_strict - NU_HF_EXPERIMENTAL_HZ
    residual_MHz_strict = residual_Hz_strict / 1.0e6
    residual_ppm_strict = residual_Hz_strict / NU_HF_EXPERIMENTAL_HZ * 1.0e6

    residual_Hz_recoil = A_hf_Hz_with_recoil - NU_HF_EXPERIMENTAL_HZ
    residual_MHz_recoil = residual_Hz_recoil / 1.0e6
    residual_ppm_recoil = residual_Hz_recoil / NU_HF_EXPERIMENTAL_HZ * 1.0e6

    # Electron anomalous moment (for context only — Track HF-2 will derive it)
    a_e_codata = 1.15965218091e-3  # CODATA 2018 electron anomalous moment
    a_e_schwinger = ALPHA_CODATA / (2.0 * 3.141592653589793)
    a_e_residual_factor_strict = residual_ppm_strict / 1.0e6  # in dimensionless
    a_e_residual_factor_recoil = residual_ppm_recoil / 1.0e6  # in dimensionless

    return {
        "inputs": {
            "alpha_CODATA": ALPHA_CODATA,
            "g_p_CODATA": G_P_CODATA,
            "m_e_over_m_p_CODATA": ME_OVER_MP_CODATA,
            "g_e_BF": g_e,
            "Z": 1,
            "n": 1,
            "Hz_per_Hartree_CODATA": HZ_PER_HARTREE,
            "nu_hf_experimental_Hz": NU_HF_EXPERIMENTAL_HZ,
            "reduced_mass_factor": REDUCED_MASS_FACTOR,
            "a_e_CODATA": a_e_codata,
            "a_e_Schwinger_alpha_over_2pi": a_e_schwinger,
        },
        "per_factor_decomposition": {
            "alpha_squared": alpha_sq,
            "g_e": g_e,
            "g_p_half_eq_mu_p_over_mu_N": g_p_half,
            "g_p": g_p,
            "m_e_over_m_p": eta,
            "psi_1s_squared_at_origin_au": psi0_sq_au,
            "psi_1s_squared_origin_symbolic": "Z^3 / pi (Z=1 -> 1/pi)",
            "fermi_contact_prefactor": "(8 pi / 3)",
            "two_magneton_factor": "(alpha/2)*(alpha/2) = alpha^2 / 4",
            "assembled_prefactor_for_Z1_ge2": "(4/3)",
        },
        "outputs_strict_no_recoil": {
            "A_hf_Ha": A_hf_Ha_strict,
            "A_hf_Hz": A_hf_Hz_strict,
            "A_hf_MHz": A_hf_MHz_strict,
            "residual_Hz_vs_experiment": residual_Hz_strict,
            "residual_MHz_vs_experiment": residual_MHz_strict,
            "residual_ppm_vs_experiment": residual_ppm_strict,
            "residual_factor_dimensionless": a_e_residual_factor_strict,
        },
        "outputs_with_reduced_mass_for_context": {
            "A_hf_Ha": A_hf_Ha_strict * REDUCED_MASS_FACTOR,
            "A_hf_Hz": A_hf_Hz_with_recoil,
            "A_hf_MHz": A_hf_MHz_with_recoil,
            "residual_Hz_vs_experiment": residual_Hz_recoil,
            "residual_MHz_vs_experiment": residual_MHz_recoil,
            "residual_ppm_vs_experiment": residual_ppm_recoil,
            "residual_factor_dimensionless": a_e_residual_factor_recoil,
            "note": "Reduced-mass factor (1 + m_e/m_p)^{-3} is a recoil correction "
                    "and is NOT part of strict Bohr-Fermi at g_e=2. Reported here "
                    "for context: it explains roughly half of the strict residual.",
        },
    }


# ---------------------------------------------------------------------------
# 5. Per-factor source attribution
# ---------------------------------------------------------------------------

def per_factor_attribution() -> Dict[str, Dict[str, str]]:
    """Catalogue of which factor came from where.

    Distinguishes:
      - GeoVac-internal: derivable from the framework's algebraic machinery.
      - GeoVac-continuum-embedding: comes from the continuum hydrogenic
        wavefunction in atomic units (i.e., the Fock projection's
        embedding of the discrete graph into continuous Hilbert space).
        Strictly, the spatial wavefunction profile is part of the
        framework's *embedding*, not the bare graph. We flag this honestly.
      - External calibration: numerical CODATA value supplied as input.
      - Standard physics (not GeoVac-specific): a result of standard
        Dirac/Pauli theory that happens to hold on S^3 trivially.
    """
    return {
        "fermi_contact_prefactor_(8 pi / 3)": {
            "value": "(8 pi / 3)",
            "source": "Standard Pauli (Dirac NR limit) Fermi contact term. "
                      "Comes from <delta(r)> matrix element of the spin-spin "
                      "magnetic interaction. Not GeoVac-specific.",
            "category": "Standard physics (Dirac NR limit)",
        },
        "g_e_equals_2": {
            "value": "2 (exact)",
            "source": "Tree-level Dirac equation. The Camporesi-Higuchi Dirac "
                      "operator on S^3 (geovac.dirac_s3) gives the same g_e = 2 "
                      "at tree level via |lambda_n| = n + 3/2, but this is the "
                      "Dirac equation acting (g_e = 2 is forced by the gamma "
                      "matrix algebra), not a GeoVac-specific result.",
            "category": "Standard physics (Dirac equation)",
        },
        "mu_B_equals_alpha_over_2": {
            "value": "alpha/2 (atomic units)",
            "source": "Bohr magneton in atomic units: e h_bar / (2 m_e c) = "
                      "alpha/2 with c = 1/alpha. Standard atomic-unit identity, "
                      "uses CODATA alpha as input.",
            "category": "Standard physics + external calibration (alpha)",
        },
        "mu_N_equals_(m_e/m_p)_alpha_over_2": {
            "value": "(m_e/m_p) * alpha/2 (atomic units)",
            "source": "Nuclear magneton in atomic units: scales mu_B by m_e/m_p. "
                      "Uses CODATA m_e/m_p as external input.",
            "category": "Standard physics + external calibration (m_e/m_p)",
        },
        "g_p_factor": {
            "value": "5.5856946893 (= 2 * mu_p/mu_N = 2 * 2.79285)",
            "source": "Proton g-factor — a measured property of the proton's "
                      "internal QCD structure. NOT computable from GeoVac at "
                      "this stage (would require lattice QCD or the full SM "
                      "with strong-coupling content). External CODATA input. "
                      "The Fermi-contact Hamiltonian carries the full g_p, "
                      "not g_p/2 — the latter equals mu_p/mu_N and only "
                      "appears if one writes mu_p = (g_p/2) mu_N as an "
                      "expectation value rather than the operator g_p mu_N S.",
            "category": "External calibration",
        },
        "psi_1s_squared_at_origin": {
            "value": "Z^3 / pi (= 1/pi for Z=1)",
            "source": "Symbolic evaluation of "
                      "geovac.dirac_matrix_elements._hydrogenic_radial_wavefunction"
                      "(n=1, l=0, r, Z) at r=0, multiplied by |Y_{0,0}|^2 = 1/(4 pi). "
                      "GeoVac uses the *continuum* hydrogenic radial wavefunction "
                      "as the spatial profile of the Fock-projected n-shell. "
                      "The graph spectrum n^2 - 1 is exact; the spatial profile "
                      "psi_1s(r) is the framework's continuum embedding (Section 12, "
                      "algebraic registry; Section 4, dimensionless vacuum principle). "
                      "Honest flag: |psi(0)|^2 is a continuum quantity, not a graph "
                      "quantity. The framework does not yet have a *graph-native* "
                      "(non-continuum) notion of |psi(0)|^2.",
            "category": "GeoVac-continuum-embedding (algebraic from hydrogen 1s)",
        },
        "alpha_squared_overall_factor": {
            "value": "alpha^2",
            "source": "Two factors of alpha/2 from the magneton pair, giving "
                      "alpha^2/4. The other factors absorb into the (4/3) in "
                      "the canonical form. External calibration via alpha.",
            "category": "External calibration (alpha twice)",
        },
        "reduced_mass_factor_(NOT_INCLUDED_IN_STRICT_BF)": {
            "value": "(1 + m_e/m_p)^{-3} ~ 0.998366",
            "source": "Recoil correction; corresponds to the proton not being "
                      "infinitely heavy. Briefing says 'no recoil', so this "
                      "factor is excluded from the headline number but reported "
                      "for context. Standard physics (reduced-mass coordinate "
                      "transformation in two-body bound state).",
            "category": "Standard physics (recoil; deliberately excluded here)",
        },
    }


# ---------------------------------------------------------------------------
# 6. Per-factor "what came from GeoVac vs hand-coded"
# ---------------------------------------------------------------------------

def honest_assessment() -> Dict[str, str]:
    """Honest pass-by-pass tally of which pieces came from GeoVac vs hand-coded.

    Three categories:
      A. From existing GeoVac machinery (cited routine).
      B. Standard physics result (Dirac eq / atomic-unit identities).
         Could in principle be re-derived from GeoVac's Dirac-on-S^3
         infrastructure but at tree level the answer is identical.
      C. External calibration (CODATA numerical input, not derivable
         from current GeoVac).
    """
    return {
        "psi_1s_squared_at_origin_(Z3/pi)":
            "A. From GeoVac: "
            "geovac/dirac_matrix_elements.py::_hydrogenic_radial_wavefunction "
            "evaluated symbolically at r=0. Spherical Y_{0,0} normalization "
            "(1/(4 pi)) is standard. Note: this is a continuum quantity (the "
            "wavefunction profile in the framework's continuum embedding), "
            "NOT a graph-native quantity — see CLAUDE.md Section 4.",
        "g_e = 2":
            "B. Standard physics: tree-level Dirac equation. GeoVac's "
            "Camporesi-Higuchi Dirac operator on S^3 has eigenvalues "
            "|lambda_n| = n + 3/2 and g_e = 2 falls out of the same gamma-matrix "
            "algebra. Not GeoVac-specific.",
        "Fermi contact prefactor (8 pi / 3)":
            "B. Standard physics: comes from the delta-function spin-spin "
            "Hamiltonian of QED in the non-relativistic limit (Bethe-Salpeter "
            "Eq. 22.18). Could in principle be re-derived from a graph-native "
            "Fermi contact integral on S^3 (an open question for future "
            "GeoVac work), but at present we use the standard prefactor.",
        "mu_B = alpha/2 in atomic units":
            "B. Standard physics: atomic unit identity using c = 1/alpha. "
            "External calibration enters via alpha (CODATA).",
        "mu_N = (m_e/m_p) alpha/2":
            "C + B. Standard physics; m_e/m_p is external CODATA input.",
        "g_p (proton g-factor)":
            "C. External calibration (CODATA). NOT computable in current "
            "GeoVac — requires QCD/lattice physics. Treat as input.",
        "alpha (fine-structure constant)":
            "C. External calibration (CODATA). Paper 2 in Observations is "
            "the framework's open conjecture about deriving alpha from the "
            "Fock-projected S^3 spectral data; we do NOT use that here. "
            "alpha is supplied as CODATA input.",
        "1 Hartree = 6.580e15 Hz":
            "C. Unit conversion only, used to render the result in MHz.",
        "Reduced-mass factor (1 + m_e/m_p)^{-3}":
            "Excluded from headline result per briefing. Standard physics "
            "two-body recoil; would require a recoil-aware Fock projection "
            "to be derivable from GeoVac (currently not available).",
    }


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 76)
    print("Sprint HF Track 1 — Bohr-Fermi A_hf for hydrogen 21 cm")
    print("=" * 76)
    print()

    # Step 1: derive |psi_1s(0)|^2 from GeoVac machinery.
    print("Step 1: deriving |psi_1s(r=0)|^2 from "
          "geovac.dirac_matrix_elements._hydrogenic_radial_wavefunction...")
    psi_derivation = derive_psi_1s_at_origin_from_geovac()
    for k, v in psi_derivation.items():
        print(f"    {k}: {v}")
    print()

    # Step 2: assemble the Bohr-Fermi formula symbolically.
    print("Step 2: symbolic assembly of A_hf = (8 pi / 3) g_e g_p mu_B mu_N "
          "|psi_1s(0)|^2 ...")
    sym = bohr_fermi_symbolic()
    print(f"    A_hf (general)     = {sp.simplify(sym['A_hf_full'])}")
    print(f"    A_hf (Z=1, g_e=2)  = {sym['A_hf_BF_Z1_ge2']}")
    print(f"    Canonical form     = {sym['A_hf_BF_canonical_form']}")
    print()

    # Step 3: numeric evaluation.
    print("Step 3: numerical evaluation with CODATA 2018 calibration ...")
    numeric = evaluate_numeric()
    out_strict = numeric["outputs_strict_no_recoil"]
    out_recoil = numeric["outputs_with_reduced_mass_for_context"]
    print()
    print(f"    Strict Bohr-Fermi (g_e=2, no recoil):")
    print(f"      A_hf = {out_strict['A_hf_Ha']:.6e} Ha")
    print(f"           = {out_strict['A_hf_MHz']:.4f} MHz")
    print(f"      Experimental nu_hf = {NU_HF_EXPERIMENTAL_MHZ:.4f} MHz")
    print(f"      Residual = {out_strict['residual_MHz_vs_experiment']:+.4f} MHz")
    print(f"               = {out_strict['residual_ppm_vs_experiment']:+.1f} ppm")
    print()
    print(f"    With reduced-mass correction (for context):")
    print(f"      A_hf = {out_recoil['A_hf_MHz']:.4f} MHz")
    print(f"      Residual = {out_recoil['residual_MHz_vs_experiment']:+.4f} MHz")
    print(f"               = {out_recoil['residual_ppm_vs_experiment']:+.1f} ppm")
    print()

    # Step 4: write JSON.
    out_dir = os.path.join(_PROJECT_ROOT, "debug", "data")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "sprint_hf_track1.json")
    payload = {
        "sprint": "HF-1",
        "date_utc": datetime.utcnow().isoformat() + "Z",
        "system": "hydrogen 1s, 21 cm hyperfine",
        "level": "Bohr-Fermi (g_e = 2 exactly, no anomalous moment, no recoil, "
                 "no nuclear structure, no multi-loop QED)",
        "psi_derivation_from_geovac": psi_derivation,
        "symbolic_A_hf_BF_Z1_ge2_canonical": str(sym["A_hf_BF_canonical_form"]),
        "numeric": numeric,
        "per_factor_attribution": per_factor_attribution(),
        "honest_assessment": honest_assessment(),
        "interpretation": {
            "headline_result_strict_BF_MHz": out_strict["A_hf_MHz"],
            "experimental_MHz": NU_HF_EXPERIMENTAL_MHZ,
            "residual_MHz_strict": out_strict["residual_MHz_vs_experiment"],
            "residual_ppm_strict": out_strict["residual_ppm_vs_experiment"],
            "expected_dominant_correction": (
                "a_e ~ alpha / (2 pi) ~ 1161 ppm (electron anomalous moment, "
                "Schwinger 1948). Track HF-2 will derive this from the "
                "GeoVac vertex correction machinery."
            ),
            "secondary_correction": (
                "(1 + m_e/m_p)^{-3} ~ -1635 ppm (reduced-mass / recoil). "
                "The strict Bohr-Fermi over-shoots experiment because no "
                "recoil and no electron anomalous moment cancel each other "
                "partially; including recoil but not a_e undershoots, "
                "including a_e but not recoil overshoots."
            ),
        },
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, default=str)
    print(f"Wrote {out_path}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
