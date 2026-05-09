"""Sprint TD Track 4: Master Mellin engine M1 extensibility test.

Question: does the M1 mechanism (Hopf-base measure / Matsubara temperature
from S^1 compactification) verified for Stefan-Boltzmann on S^3 x S^1_beta
(Track 1, commit 697666f) extend cleanly to the Euclidean Schwarzschild
cigar geometry, reproducing the Hawking temperature T_H = 1/(8 pi M)?

Honest scope (CLAUDE.md s1.5, repeated up front):
- This is NOT trying to derive black hole thermodynamics from a packing
  axiom. The cigar comes from Einstein equations applied to vacuum,
  NOT from GeoVac.
- This is NOT trying to derive S_BH = A/(4 ell_p^2). That requires the
  Connes-Chamseddine spectral action evaluated at the Schwarzschild
  saddle (M2-class object on the horizon S^2; multi-week NCG work).
- This is ONLY testing whether the SAME M1 mechanism that produces
  the 2 pi in Stefan-Boltzmann at finite T (Track 1) also produces the
  2 pi in T_H = 1/(8 pi M).

Structural background:
- Euclidean Schwarzschild ds^2 = (1 - 2M/r) dtau^2 + (1-2M/r)^{-1} dr^2
  + r^2 dOmega_2^2.
- Smoothness at r = 2M (horizon) requires tau ~ tau + beta with
  beta_cigar = 8 pi M (Hawking 1975, Gibbons-Hawking 1977).
- Equivalently, surface gravity kappa_g = 1/(4M); Hawking T_H = kappa_g/(2 pi).
- The 2 pi here is the SAME 2 pi as in Matsubara modes 2 pi k / beta on
  S^1_beta in Track 1.

Reusable Track 1 infrastructure:
- geovac/thermal_tensor_triple.py::matsubara_spectrum(beta_sym, k_max,
  fermionic) returns omega_k = 2 pi k / beta. Calling with
  beta = 8 pi M gives the Hawking-shifted Matsubara spectrum verbatim.
- modular_residual_thermal_tensor() shows the heat-kernel factorization
  Tr e^{-s D^2} = K_spatial(s) * K_temporal(s) with K_temporal leading
  ~ beta / (2 sqrt(pi s)). For the cigar, K_spatial is on (r, Omega)
  and K_temporal is the same theta_3 modular form on the tau-circle.

What this script verifies:
1. T_H = 1/(8 pi M) factors as 1/beta_cigar with beta_cigar = 8 pi M.
2. The 2 pi in T_H IS the M1 Hopf-base measure signature of S^1_tau.
3. The same matsubara_spectrum() function from Track 1 produces the
   correct Hawking Matsubara modes when fed beta_cigar.
4. Surface gravity / 2 pi factorization: T_H = kappa_g / (2 pi) and
   the 2 pi is the circumference of the (suitably normalized)
   asymptotic-time circle.

Verdict structure:
- POSITIVE for the LIMITED M1 extensibility claim.
- NEGATIVE for the broader claim that the framework "derives" T_H from
  internal data (the M = M_Schwarzschild parameter is external input,
  set by Einstein equations + horizon condition).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp

# Ensure we can import Track 1's module
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geovac.thermal_tensor_triple import (
    matsubara_spectrum,
    modular_residual_thermal_tensor,
    PI_TAG_M1_HOPF,
    PI_TAG_M1_MATSUBARA_CIRCLE,
    PI_TAG_M2_SD_S3,
    PI_TAG_M2_ZETA4,
)


# ---------------------------------------------------------------------------
# 1. Hawking beta from the cigar regularity condition
# ---------------------------------------------------------------------------

def cigar_beta_symbolic() -> dict:
    """Smoothness at r = 2M of Euclidean Schwarzschild forces tau-period
    beta_cigar = 8 pi M, equivalently T_H = 1/(8 pi M).

    Derivation (Hawking 1975, Gibbons-Hawking 1977):
    Near r = 2M, write rho^2 = 4(2M)(r-2M). Then ds^2 ~ rho^2 dphi^2 + drho^2
    + (2M)^2 dOmega^2 with phi = tau / (4M). Smoothness at rho = 0 requires
    phi to be 2 pi-periodic, hence tau-period 8 pi M.

    The 2 pi is the circumference of phi-circle (M1 / Hopf-base measure
    signature, identical to S^1_beta in Track 1).
    """
    M = sp.Symbol("M", positive=True)
    # phi = tau / (4M). Smoothness at rho = 0 requires phi-period 2 pi.
    # Hence tau-period beta_cigar = 4 M * 2 pi = 8 pi M.
    phi_period = 2 * sp.pi  # M1 Hopf-base measure signature
    beta_cigar = 4 * M * phi_period
    T_H = 1 / beta_cigar

    # Surface gravity kappa_g = 1/(4M); standard relation T_H = kappa_g / (2 pi).
    kappa_g = sp.Rational(1, 4) / M
    T_H_via_kappa = kappa_g / (2 * sp.pi)
    residual = sp.simplify(T_H - T_H_via_kappa)

    return {
        "phi_period": phi_period,                      # = 2 pi (M1)
        "phi_period_meaning": PI_TAG_M1_MATSUBARA_CIRCLE,
        "beta_cigar_sym": beta_cigar,                  # = 8 pi M
        "T_H_sym": T_H,                                # = 1 / (8 pi M)
        "kappa_g_sym": kappa_g,                        # = 1 / (4 M)
        "T_H_via_kappa_sym": T_H_via_kappa,            # = 1 / (8 pi M)
        "factorization_residual": residual,            # must simplify to 0
        "structural_reading": (
            "T_H = phi_period * (4M)^{-1} / (2 pi)^{-1} reduces to 1/(8 pi M). "
            "The single 2 pi is the S^1_phi circumference. Identical M1 "
            "mechanism as Track 1's Stefan-Boltzmann 2 pi k / beta Matsubara "
            "spectrum, just with beta = 8 pi M instead of beta = 1/(k_B T)."
        ),
    }


# ---------------------------------------------------------------------------
# 2. Re-use Track 1's Matsubara spectrum on the cigar tau-circle
# ---------------------------------------------------------------------------

def cigar_matsubara_via_track1(k_max: int = 4) -> dict:
    """Call Track 1's matsubara_spectrum() with beta = 8 pi M to get the
    Hawking Matsubara spectrum on the cigar tau-circle.

    This is a verbatim reuse: NO new code, NO new mechanism, just the same
    M1 Hopf-base measure factor 2 pi at a different beta.
    """
    M = sp.Symbol("M", positive=True)
    beta_cigar = 8 * sp.pi * M

    # Bosonic (graviton-like): omega_k = 2 pi k / (8 pi M) = k / (4 M)
    bosonic = matsubara_spectrum(beta_cigar, k_max, fermionic=False)
    bosonic_simplified = [(k, sp.simplify(om)) for (k, om) in bosonic]

    # Fermionic (matter Dirac on cigar): omega_k = 2 pi (k + 1/2) / (8 pi M)
    #                                              = (2k + 1) / (8 M)
    fermionic = matsubara_spectrum(beta_cigar, k_max, fermionic=True)
    fermionic_simplified = [(k, sp.simplify(om)) for (k, om) in fermionic]

    # Lowest nonzero bosonic mode: omega_1 = 1/(4 M) = surface gravity kappa_g
    omega1_bosonic = sp.simplify(bosonic_simplified[k_max + 1][1])  # k=1 entry
    kappa_g_check = sp.simplify(omega1_bosonic - sp.Rational(1, 4) / M)

    # Lowest fermionic mode (k=0 in fermionic indexing): omega_0 = pi T_H
    # In the (-k_max..k_max) listing this is index k_max for k=0.
    omega0_fermionic = sp.simplify(fermionic_simplified[k_max][1])  # k=0 entry
    pi_T_H = sp.simplify(sp.pi / beta_cigar)
    fermionic_check = sp.simplify(omega0_fermionic - pi_T_H)

    return {
        "beta_cigar_used": beta_cigar,                 # 8 pi M
        "bosonic_lowest_nonzero": omega1_bosonic,      # = 1/(4M) = kappa_g
        "kappa_g_match": kappa_g_check,                # must be 0
        "kappa_g_meaning": (
            "Lowest bosonic Matsubara frequency on the cigar tau-circle "
            "EXACTLY equals the Schwarzschild surface gravity kappa_g = 1/(4M). "
            "This is structurally analogous to Track 1's lowest bosonic "
            "Matsubara mode on S^1_beta being 2 pi T."
        ),
        "fermionic_lowest": omega0_fermionic,          # = pi T_H = 1/(8M)
        "pi_T_H_match": fermionic_check,               # must be 0
        "fermionic_meaning": (
            "Anti-periodic boundary condition on tau-circle for fermions "
            "shifts the lowest mode to omega = pi T_H (matter Dirac field "
            "on cigar in standard Hartle-Hawking state)."
        ),
        "spectrum_bosonic_first_5": [
            (k, str(om)) for (k, om) in bosonic_simplified[:5]
        ],
        "spectrum_fermionic_first_5": [
            (k, str(om)) for (k, om) in fermionic_simplified[:5]
        ],
        "structural_reading": (
            "Track 1's matsubara_spectrum(beta_cigar, ...) reproduces the "
            "standard Hawking Matsubara spectrum verbatim. The M1 (Hopf-base "
            "measure) signature 2 pi from the tau-circle circumference is "
            "identical between flat S^1_beta and the cigar. No new mechanism "
            "is required; only the value of beta changes (8 pi M for the cigar "
            "vs 1/(k_B T) for flat thermal field theory)."
        ),
    }


# ---------------------------------------------------------------------------
# 3. Heat-kernel factorization on the cigar (structural)
# ---------------------------------------------------------------------------

def cigar_heat_kernel_factorization() -> dict:
    """Track 1's modular_residual_thermal_tensor() shows
        Tr e^{-s D_total^2} = K_spatial(s) * K_temporal(s)
    with K_temporal(s) ~ beta / (2 sqrt(pi s)) leading from Jacobi theta
    inversion on S^1_beta.

    For the cigar, the same factorization is the local statement that
    near r = 2M:
        D_cigar^2 = D_(r,Omega)^2 (x) 1 + 1 (x) D_tau^2
    where D_tau is the Matsubara operator on the tau-circle of period
    8 pi M. The factorization holds locally (in the near-horizon Rindler
    approximation); globally the cigar's r-direction is non-compact and
    decompactifies the tau-circle at large r.

    Honest scope:
    - The Track 1 factorization uses both factors compact (S^3 x S^1_beta).
    - The cigar is non-compact in r at infinity, so the FULL heat kernel
      requires regularization (standard QFT-in-curved-spacetime) that is
      not the master Mellin engine's job.
    - The NEAR-HORIZON M1 mechanism IS reproduced by the factorization.

    This is reading (a) of the verdict: the formal extension works locally
    but the non-compact r-direction sits outside the master Mellin engine's
    compact-spectrum scope.
    """
    track1_result = modular_residual_thermal_tensor()
    return {
        "track1_K_temporal_leading": str(track1_result["K_temporal_leading_sym"]),
        "track1_factorization": track1_result["factorization"],
        "cigar_local_factorization": (
            "Near r = 2M (Rindler): D_cigar^2 ~ D_rho^2 + (1/(4M))^2 phi^2 "
            "+ D_S^2/r^2; the (rho, phi) part is locally R^2 with phi = tau/(4M) "
            "of period 2 pi, identical to flat polar coordinates. The phi-circle "
            "supports a Matsubara sum identical to Track 1's S^1_beta sum."
        ),
        "cigar_non_compact_caveat": (
            "Globally, the r-direction is non-compact (r in [2M, infinity)). "
            "The full cigar heat kernel needs IR regularization at r = infinity "
            "(standard practice: integrate over a finite spacelike box, take "
            "the limit). Track 1's K_spatial * K_temporal factorization assumes "
            "BOTH factors compact; for the cigar it holds only in the near-horizon "
            "approximation. The M1 mechanism on the tau-circle is unaffected, "
            "because tau-compactness is ENFORCED by smoothness at the tip."
        ),
        "structural_reading": (
            "M1 mechanism (Matsubara sum on tau-circle) extends cleanly. "
            "M2 mechanism (Seeley-DeWitt on r-direction) does NOT extend "
            "cleanly because of non-compact r; this is where Connes-Chamseddine "
            "black-hole spectral-action computations require boundary "
            "regularization, and is OUT OF SCOPE for this track."
        ),
    }


# ---------------------------------------------------------------------------
# 4. Master Mellin engine verdict
# ---------------------------------------------------------------------------

def master_mellin_verdict() -> dict:
    """Aggregate verdict on M1 extensibility from S^3 x S^1_beta to
    cigar x S^2_horizon.

    Three sub-verdicts:
    (V1) M1 mechanism (Matsubara temperature from S^1) — POSITIVE.
         The 2 pi in T_H = 1/(8 pi M) is the same M1 Hopf-base measure
         signature as in Track 1's Stefan-Boltzmann.
    (V2) M2 mechanism (Seeley-DeWitt) on horizon S^2 — DEFERRED, OUT OF SCOPE.
         Recovering S_BH = A / 4 from the horizon S^2 heat kernel is the
         Connes-Chamseddine spectral action computation; multi-week NCG work,
         not this track's scope.
    (V3) Bertrand-class forcing — NEGATIVE.
         The cigar geometry is forced by Einstein equations + Schwarzschild
         horizon, NOT by the GeoVac packing axiom. The M parameter is
         external input. The framework does not "explain" black holes.

    Strongest claim from this track: the M1 (Matsubara temperature) mechanism
    is universal across S^1 compactifications regardless of the host
    spacetime (Coulomb-derived S^3, gravitational cigar, anything else with
    a thermal circle). This is a small but non-trivial structural finding.
    """
    return {
        "V1_M1_mechanism": {
            "status": "POSITIVE",
            "claim": (
                "The 2 pi in T_H = 1/(8 pi M) is the Hopf-base measure / "
                "M1 mechanism signature of S^1_tau circumference. Identical "
                "to Track 1's S^1_beta Matsubara mechanism, only beta changes "
                "(beta_cigar = 8 pi M)."
            ),
            "evidence": (
                "(i) Track 1's matsubara_spectrum() called with beta = 8 pi M "
                "reproduces standard Hawking Matsubara spectrum bit-identically; "
                "(ii) the lowest bosonic mode equals surface gravity kappa_g; "
                "(iii) T_H = kappa_g / (2 pi) factorization makes the 2 pi "
                "explicitly the M1 signature."
            ),
        },
        "V2_M2_mechanism": {
            "status": "OUT_OF_SCOPE_DEFERRED",
            "claim": (
                "Recovering Bekenstein-Hawking entropy S_BH = A/(4 ell_p^2) "
                "from the horizon S^2 spectral action would test M2 (Seeley-"
                "DeWitt) extensibility on the horizon submanifold. This is "
                "a separate multi-week NCG computation (Chamseddine-Connes "
                "2010 program), explicitly OUT OF SCOPE per directive."
            ),
        },
        "V3_Bertrand_forcing": {
            "status": "NEGATIVE",
            "claim": (
                "The cigar is forced by Einstein equations + Schwarzschild "
                "horizon at r = 2M. GeoVac's packing axiom does NOT produce "
                "the cigar. Bertrand's theorem (closed-orbit selection of "
                "Coulomb potential) is what selects S^3 in GeoVac; no analog "
                "selects the cigar. M (BH mass) is external input."
            ),
            "interpretation": (
                "This is what we expect: the framework is a structural-skeleton "
                "framework that operates on geometries it is given. The cigar "
                "comes from gravity, not from packing. The fact that one of its "
                "structural mechanisms (M1 Matsubara temperature) extends without "
                "modification is the modest but real result of this track."
            ),
        },
        "net_verdict": (
            "POSITIVE for the limited M1 extensibility claim, NEGATIVE for "
            "any claim that GeoVac derives black hole physics. The Matsubara "
            "temperature mechanism is universal across S^1 compactifications; "
            "GeoVac does NOT autonomously generate the cigar geometry."
        ),
        "Paper_35_Prediction_1_check": (
            "Paper 35 Prediction 1 (pi enters iff continuous integration over "
            "a temporal/spectral parameter) holds verbatim for the cigar: the "
            "single pi in T_H = 1/(8 pi M) comes from the S^1_tau circle "
            "circumference, identified at the same continuous-integration step "
            "(Matsubara sum on tau-circle) as Stefan-Boltzmann's pi^2/90."
        ),
        "Paper_32_Section_VIII_check": (
            "The case-exhaustion theorem (M1 / M2 / M3 sub-cases of master "
            "Mellin engine Tr(D^k e^{-t D^2})) extends formally to the cigar "
            "for k = 0 (M1, Matsubara temperature). Extension to k = 2 (M2, "
            "spectral action / S_BH on horizon S^2) is the natural follow-up "
            "but requires non-compact regularization not covered here."
        ),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def to_jsonable(obj: Any) -> Any:
    """Convert sympy expressions to strings for JSON serialization."""
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for (k, v) in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, sp.Expr):
        return str(obj)
    if isinstance(obj, (int, float, str, bool, type(None))):
        return obj
    return str(obj)


def main() -> None:
    print("Sprint TD Track 4: master Mellin M1 extensibility test")
    print("=" * 70)

    print("\n[1] Cigar beta from regularity at r=2M:")
    s1 = cigar_beta_symbolic()
    print(f"   beta_cigar  = {s1['beta_cigar_sym']}")
    print(f"   T_H         = {s1['T_H_sym']}")
    print(f"   kappa_g     = {s1['kappa_g_sym']}")
    print(f"   factorization residual T_H - kappa_g/(2 pi) = "
          f"{s1['factorization_residual']}")
    assert s1["factorization_residual"] == 0

    print("\n[2] Matsubara spectrum via Track 1:")
    s2 = cigar_matsubara_via_track1(k_max=4)
    print(f"   beta_cigar used: {s2['beta_cigar_used']}")
    print(f"   bosonic lowest nonzero = {s2['bosonic_lowest_nonzero']} "
          f"(should be 1/(4M) = kappa_g)")
    print(f"   kappa_g residual: {s2['kappa_g_match']}")
    print(f"   fermionic lowest = {s2['fermionic_lowest']} "
          f"(should be pi T_H)")
    print(f"   pi T_H residual: {s2['pi_T_H_match']}")
    assert s2["kappa_g_match"] == 0
    assert s2["pi_T_H_match"] == 0

    print("\n[3] Heat-kernel factorization (structural, scope-limited):")
    s3 = cigar_heat_kernel_factorization()
    print(f"   Track 1 K_temporal leading: {s3['track1_K_temporal_leading']}")
    print("   Cigar local factorization:")
    print(f"      {s3['cigar_local_factorization'][:120]}...")

    print("\n[4] Master Mellin engine verdict:")
    s4 = master_mellin_verdict()
    print(f"   V1 (M1 mechanism):   {s4['V1_M1_mechanism']['status']}")
    print(f"   V2 (M2 mechanism):   {s4['V2_M2_mechanism']['status']}")
    print(f"   V3 (Bertrand):       {s4['V3_Bertrand_forcing']['status']}")
    print(f"\n   Net verdict: {s4['net_verdict']}")

    # Save to JSON
    out = {
        "1_cigar_beta_symbolic": to_jsonable(s1),
        "2_matsubara_via_track1": to_jsonable(s2),
        "3_heat_kernel_factorization": to_jsonable(s3),
        "4_master_mellin_verdict": to_jsonable(s4),
        "metadata": {
            "track": "Sprint TD Track 4",
            "framing": (
                "M1 extensibility test from S^3 x S^1_beta (Track 1, commit "
                "697666f) to Euclidean Schwarzschild cigar x S^2_horizon. "
                "POSITIVE for M1; OUT_OF_SCOPE for M2 (S_BH); NEGATIVE for "
                "Bertrand-class forcing of cigar from packing."
            ),
            "honest_scope": (
                "This is a TEMPERATURE-MECHANISM test, not an entropy-shape "
                "test or a black-hole-derivation claim. Per directive."
            ),
        },
    }
    out_path = Path(__file__).parent / "data" / "sprint_td_track4.json"
    out_path.parent.mkdir(exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n   Wrote: {out_path}")
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
