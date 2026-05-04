"""TX-B Observable 2: Spatial Casimir energy on Bargmann-Segal S^5 for the
free 3D HO field.

POSITIVE PREDICTION (KG-3 ANALOG ON S^5).

The Bargmann-Segal lattice (Paper 24) encodes the 3D isotropic harmonic
oscillator on the (N, 0) symmetric SU(3) tower with spectrum
    E_N = hbar omega (N + 3/2),  degeneracy d_N = (N+1)(N+2)/2.

We compute the zeta-regularized Casimir energy

    E_Cas = (1/2) sum_{N >= 0} d_N * E_N
          = (1/2) hbar omega * sum_{N>=0} (N+1)(N+2)/2 * (N + 3/2)

via spectral zeta evaluated at s = -1.  The spectral zeta is

    zeta_X(s) = sum_{N >= 0} d_N * E_N^{-s} / (hbar omega)^{-s}.

Setting hbar omega = 1 for the unit-frequency case (and absorbing the factor
1/2 of d_N = (N+1)(N+2)/2):

    zeta_X(s) = (1/2) sum_{N >= 0} (N+1)(N+2) * (N + 3/2)^{-s}.

Use the Hurwitz form: shift m = N + 3/2.  Then N+1 = m - 1/2, N+2 = m + 1/2,
so (N+1)(N+2) = m^2 - 1/4, giving

    zeta_X(s) = (1/2) sum_{m=3/2,5/2,...} (m^2 - 1/4) m^{-s}
              = (1/2) [zeta_R(s-2, 3/2) - (1/4) zeta_R(s, 3/2)]

(where zeta_R(s, a) is the Hurwitz zeta function).  This is structurally
identical to the Dirac case on S^3 (KG-5)!  In fact:

    zeta_X^HO_S5 (s) = (1/8) zeta_|D|^Dirac_S3 (s)

at s = -1, since zeta_|D|^Dirac_S3 = 4 [zeta_R(s-2, 3/2) - (1/4) zeta_R(s, 3/2)]
and zeta_X^HO_S5 = (1/2) of the same bracket = (1/8) * 4 [bracket] = (1/2) bracket.

So:
    zeta_X^HO_S5(-1) = (1/2) [zeta_R(-3, 3/2) - (1/4) zeta_R(-1, 3/2)]
                    = (1/2) * (-17/240) / 4   (from KG-5: bracket = -17/240/4)
                    = (1/2) * (-17/960)
                    = -17/1920

Casimir BOSON sign (no minus): E_Cas = (1/2) zeta_X(-1) = -17/3840.

Wait, let me redo this carefully.

zeta_R(-1, 3/2) = -B_2(3/2)/2 = -(11/12)/2 = -11/24.
zeta_R(-3, 3/2) = -B_4(3/2)/4 = -(127/240)/4 = -127/960.

Bracket = zeta_R(-3, 3/2) - (1/4) zeta_R(-1, 3/2)
       = -127/960 - (1/4)(-11/24)
       = -127/960 + 11/96
       = -127/960 + 110/960
       = -17/960

So zeta_X^HO_S5(-1) = (1/2) * (-17/960) = -17/1920.

E_Cas (boson) = (1/2) zeta_X(-1) = -17/3840.

This is an EXACT RATIONAL.  No pi appears at any step.

The Paper 35 prediction is CONFIRMED in this sector if the result is rational.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
import mpmath as mp

OUT = Path(__file__).parent / "data" / "tx_b_obs_2.json"
mp.mp.dps = 50


def bargmann_s5_HO_casimir():
    """Compute E_Cas for 3D HO on Bargmann-Segal S^5, zeta-regularized.

    Spectrum E_N = N + 3/2, degeneracy d_N = (N+1)(N+2)/2.
    """
    a = sp.Rational(3, 2)

    # Hurwitz zeta values at negative integers via Bernoulli polynomials.
    B2_a = sp.bernoulli(2, a)             # 11/12
    B4_a = sp.bernoulli(4, a)             # 127/240
    zeta_H_m1 = -B2_a / 2                 # -11/24
    zeta_H_m3 = -B4_a / 4                 # -127/960

    # Spectral zeta at s = -1:
    # zeta_X(s) = (1/2) [zeta_R(s-2, 3/2) - (1/4) zeta_R(s, 3/2)]
    bracket_at_m1 = zeta_H_m3 - sp.Rational(1, 4) * zeta_H_m1
    bracket_at_m1 = sp.simplify(bracket_at_m1)
    zeta_X_m1 = sp.Rational(1, 2) * bracket_at_m1

    # Casimir energy (boson sign): E_Cas = (1/2) zeta_X(-1)
    E_Cas = sp.Rational(1, 2) * zeta_X_m1
    E_Cas_simplified = sp.simplify(E_Cas)

    # Numeric cross-check via mpmath Hurwitz at half-integer shift.
    val_m3_num = mp.zeta(-3, mp.mpf("1.5"))
    val_m1_num = mp.zeta(-1, mp.mpf("1.5"))
    bracket_num = val_m3_num - mp.mpf("0.25") * val_m1_num
    zeta_X_num = mp.mpf("0.5") * bracket_num
    E_Cas_num = mp.mpf("0.5") * zeta_X_num

    # Direct sum cross-check (no zeta-reg) for the partial sum of EXTRA terms,
    # showing that the full sum diverges as expected; we need zeta reg.
    # Just symbolically demonstrate that the per-N contribution lives in Q.
    direct_terms = []
    for N in range(0, 6):
        EN = sp.Rational(2 * N + 3, 2)        # N + 3/2
        dN = sp.Rational((N + 1) * (N + 2), 2)
        contrib = sp.Rational(1, 2) * dN * EN
        direct_terms.append({
            "N": N,
            "E_N": str(EN),
            "d_N": str(dN),
            "contribution_to_E_Cas_no_reg": str(contrib),
            "contribution_in_Q": isinstance(contrib, sp.Rational),
        })

    transcendental_ledger = [
        {"step": "(1) Bargmann-Segal HO spectrum E_N = N + 3/2 (Paper 24, Eq 12)",
         "value": "E_N rational (half-integer); d_N = (N+1)(N+2)/2 rational",
         "transcendentals": "NONE"},
        {"step": "(2) Spectral zeta zeta_X(s) = (1/2) sum_N (N+1)(N+2) (N+3/2)^{-s}",
         "value": "= (1/2) [zeta_R(s-2, 3/2) - (1/4) zeta_R(s, 3/2)] (Hurwitz form)",
         "transcendentals": "NONE in symbolic form"},
        {"step": "(3) Bernoulli polynomial values at 3/2",
         "value": f"B_2(3/2) = {B2_a}, B_4(3/2) = {B4_a}",
         "transcendentals": "NONE -- exact rationals"},
        {"step": "(4) Hurwitz zeta values via -B_{n+1}(a)/(n+1)",
         "value": f"zeta_R(-1, 3/2) = {zeta_H_m1}, zeta_R(-3, 3/2) = {zeta_H_m3}",
         "transcendentals": "NONE -- exact rationals"},
        {"step": "(5) Bracket [zeta_R(-3, 3/2) - (1/4) zeta_R(-1, 3/2)]",
         "value": f"= {bracket_at_m1}",
         "transcendentals": "NONE"},
        {"step": "(6) Spectral zeta at -1: zeta_X(-1) = (1/2) bracket",
         "value": f"= {zeta_X_m1}",
         "transcendentals": "NONE"},
        {"step": "(7) Casimir energy E_Cas = (1/2) zeta_X(-1)",
         "value": f"= {E_Cas_simplified}",
         "transcendentals": "NONE -- pure rational"},
    ]

    pi_in_result = "pi" in str(E_Cas_simplified) or sp.pi in E_Cas_simplified.atoms()

    return {
        "spectrum_source": "Paper 24 Eq (12), HO on (N,0) Bargmann-Segal tower",
        "spectrum": "E_N = N + 3/2, deg d_N = (N+1)(N+2)/2",
        "spectral_zeta_form": "zeta_X(s) = (1/2)[zeta_R(s-2, 3/2) - (1/4) zeta_R(s, 3/2)]",
        "intermediate_zeta_values": {
            "B_2(3/2)": str(B2_a),
            "B_4(3/2)": str(B4_a),
            "zeta_R(-1, 3/2)": str(zeta_H_m1),
            "zeta_R(-3, 3/2)": str(zeta_H_m3),
            "bracket": str(bracket_at_m1),
            "zeta_X(-1)": str(zeta_X_m1),
        },
        "E_Cas_symbolic": str(E_Cas_simplified),
        "E_Cas_numeric_50dps_mpmath": mp.nstr(E_Cas_num, 30),
        "is_rational": isinstance(E_Cas_simplified, sp.Rational),
        "contains_pi": pi_in_result,
        "matches_KG5_Dirac_S3_relation": "zeta_X^HO_S5(-1) = (1/8) zeta_|D|^Dirac_S3(-1) "
                                          "structurally; both arise from the same Hurwitz "
                                          "shift a=3/2 mechanism",
        "direct_partial_sum_first_6_terms": direct_terms,
        "transcendental_ledger": transcendental_ledger,
    }


def main():
    res = bargmann_s5_HO_casimir()

    pi_predicted = False
    pi_observed = res["contains_pi"]
    is_rational = res["is_rational"]

    out = {
        "observable_id": 2,
        "observable_name": "Spatial Casimir on Bargmann-Segal S^5 for 3D HO",
        "predicted_class": "exact rational, no pi",
        "pi_predicted": pi_predicted,
        "computed_class": ("exact rational" if is_rational else "non-rational")
                          + " (E_Cas = " + res["E_Cas_symbolic"] + ")",
        "pi_observed": pi_observed,
        "match_with_prediction": (pi_predicted == pi_observed) and is_rational,
        "result": res,
        "verdict": ("CONFIRMS Paper 35 Prediction 1" if (not pi_observed and is_rational)
                     else "REFUTES Paper 35 Prediction 1"),
        "paper_34_projections": [
            "Bargmann-Segal projection (Paper 24)",
            "Casimir/zeta-regularized integration over discrete spatial spectrum"
        ],
        "temporal_integration_present": False,
        "specific_integration": "None -- pure spatial mode-counting; no Matsubara, "
                                "no proper-time, no Mellin",
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps({
        "observable": "Bargmann-Segal S^5 HO Casimir",
        "predicted_pi": pi_predicted,
        "observed_pi": pi_observed,
        "is_rational": is_rational,
        "E_Cas_symbolic": res["E_Cas_symbolic"],
        "match": (pi_predicted == pi_observed) and is_rational,
    }, indent=2))


if __name__ == "__main__":
    main()
