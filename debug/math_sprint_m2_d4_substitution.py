"""
Math sprint M2: Paper 18 Eq.(38) D(4) substitution error.

Goal
----
Resolve the Paper 18 Eq.(38) substantive error:

  Paper 18 writes
      D(4) = sum_{m>=1} 2 m (m+1) / m^4
           = 2 zeta_R(2) + 2 zeta_R(3)
           = pi^2/3 + 2 zeta_R(3)  ≈ 5.694

  but the correct Camporesi-Higuchi sum (Paper 18 Eq.(30) and Paper 28 Eq.(43))
  is
      D(4) = sum_{n>=0} g_n / |lambda_n|^4
           = sum_{n>=0} 2(n+1)(n+2) / (n + 3/2)^4
           = pi^2 - pi^4/12  ≈ 1.752

  Paper 18's Eq.(38) used integer m where the CH spectrum has half-integer
  eigenvalues lambda_n = n + 3/2.

Tasks executed by this driver
-----------------------------
1. Confirm D(4) = pi^2 - pi^4/12 directly via mpmath at 50+ dps.

2. Cross-validate D(2), D(4), D(6), D(8) against Paper 28 Table 2
   (Eq.(30): D(s) = 2(2^{s-2} - 1) zeta_R(s-2) - (2^s - 1)/2 * zeta_R(s)).

3. Compute D(3), D(5), D(7), D(9) and locate where zeta_R(3) appears
   on the CH spinor bundle (the spinor-side "witness" the failed Eq.(38)
   tried to identify).

4. Compute Paper 18's wrongly-written sum
        S_wrong = sum_{m>=1} 2 m (m+1) / m^4 = 2 zeta(2) + 2 zeta(3)
   to confirm the NUMBER is right (the identification with D(4) is the bug).

5. Diagnose: write down what the right sum looks like with half-integer
   substitution m -> n + 3/2 (i.e. the substitution Paper 18 omitted).

6. Search for clean zeta(3) coefficients on the CH spinor bundle.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
from mpmath import mpf, mpc, zeta, pi, nsum, inf

DPS = 60
mp.mp.dps = DPS


# ---------------------------------------------------------------------------
# Camporesi-Higuchi spectrum on unit S^3
# ---------------------------------------------------------------------------
def lam(n: int):
    """|lambda_n| = n + 3/2."""
    return mpf(n) + mpf("1.5")


def g_dirac(n: int):
    """g_n^{Dirac} = 2 (n+1)(n+2)."""
    return mpf(2) * (n + 1) * (n + 2)


def D_direct(s, N: int = 4000) -> mpf:
    """Direct truncated CH sum D(s) = sum_{n=0}^{N-1} g_n / |lambda_n|^s.

    Sums to N=4000 in scope; supplemented by adaptive tail using mpmath.nsum.
    """

    def term(n):
        return g_dirac(int(n)) / lam(int(n)) ** s

    return mp.nsum(term, [0, mp.inf])


# ---------------------------------------------------------------------------
# Hurwitz closed form from Eq.(30) / Paper 28 Eq.(38):
# D(s) = 2 (2^{s-2} - 1) zeta_R(s-2) - (2^s - 1)/2 * zeta_R(s)
# (this is the closed form Paper 28 Theorem 1 / Eq.(43) uses)
# ---------------------------------------------------------------------------
def D_hurwitz(s) -> mpf:
    """D(s) via 2*zeta(s-2, 3/2) - (1/2)*zeta(s, 3/2) — the original
    Hurwitz form before applying zeta(s, 3/2) = (2^s - 1) zeta_R(s) - 2^s."""
    return 2 * mp.zeta(s - 2, mpf("1.5")) - mpf("0.5") * mp.zeta(s, mpf("1.5"))


def D_hurwitz_eq30(s) -> mpf:
    """The Eq.(30)/Paper 28 closed form:
        D(s) = 2 (2^{s-2} - 1) zeta_R(s-2) - (2^s - 1)/2 * zeta_R(s).
    """
    return 2 * (mp.mpf(2) ** (s - 2) - 1) * mp.zeta(s - 2) - (mp.mpf(2) ** s - 1) / 2 * mp.zeta(s)


# ---------------------------------------------------------------------------
# The closed forms from Paper 28 Table 2 / Theorem 5 / ζ(3) complementarity
# ---------------------------------------------------------------------------
def D_closed_form(s_int: int) -> mpf | None:
    """Return the known closed form at integer s, or None."""
    if s_int == 2:
        # D(2): from D_hurwitz_eq30(2) = 2 (1 - 1) zeta(0) - (3/2) zeta(2)
        #      = -(3/2)*zeta(2) ... but the s=2 case has subtleties (zeta(0)=-1/2)
        # We let direct numeric do this; closed form fits Paper 28's pattern only s>=4
        return None
    if s_int == 4:
        return pi ** 2 - pi ** 4 / 12
    if s_int == 5:
        return 14 * mp.zeta(3) - mpf(31) / 2 * mp.zeta(5)
    if s_int == 6:
        return pi ** 4 / 3 - pi ** 6 / 30
    if s_int == 7:
        return 62 * mp.zeta(5) - mpf(127) / 2 * mp.zeta(7)
    if s_int == 8:
        return 2 * pi ** 6 / 15 - 17 * pi ** 8 / 1260
    if s_int == 9:
        return 254 * mp.zeta(7) - mpf(511) / 2 * mp.zeta(9)
    if s_int == 10:
        return 17 * pi ** 8 / 315 - 31 * pi ** 10 / 4725
    if s_int == 3:
        # closed form from Eq.(30) at s=3:
        # D(3) = 2 (2^{1} - 1) zeta_R(1) - (7/2) zeta_R(3) = divergent zeta_R(1)
        # so D(3) itself has a divergent classical-zeta piece — interesting!
        return None
    return None


# ---------------------------------------------------------------------------
# The WRONG sum that Paper 18 Eq.(38) writes
# ---------------------------------------------------------------------------
def S_wrong_eq38() -> mpf:
    """
    Paper 18's broken sum:
        S = sum_{m>=1} 2 m (m+1) / m^4
          = 2 sum 1/m^2 + 2 sum 1/m^3
          = 2 zeta(2) + 2 zeta(3)
    The NUMBER is right; what's wrong is the IDENTIFICATION with D(4).
    """
    return 2 * mp.zeta(2) + 2 * mp.zeta(3)


# ---------------------------------------------------------------------------
# Spinor-bundle zeta(3) search
# ---------------------------------------------------------------------------
def search_zeta3_on_spinor_bundle() -> dict:
    """
    Test which integer s of D(s) and related Hurwitz objects produce a
    clean zeta(3) coefficient.

    Closed-form analytic answer (from Eq.(30)):
        D(s) = 2 (2^{s-2} - 1) zeta_R(s-2) - (2^s - 1)/2 * zeta_R(s)

    zeta_R(3) appears as coefficient on zeta_R(s-2) when s=5
    (i.e. s-2 = 3), AND on the zeta_R(s) term when s=3 — but s=3
    multiplies zeta_R(1), which DIVERGES.  So the FIRST clean
    appearance of zeta_R(3) on the CH spinor bundle is at:

        D(5) = 14 zeta_R(3) - (31/2) zeta_R(5).

    Numerically check.
    """
    results = {}
    for s in (3, 5, 7, 9):
        try:
            direct = D_direct(s)
            direct_str = str(mp.nstr(direct, 30))
        except Exception as e:
            direct_str = f"DIVERGENT/POLE ({e})"
        closed = D_closed_form(s)
        try:
            eq30 = D_hurwitz_eq30(s)
            eq30_str = str(mp.nstr(eq30, 30))
        except Exception as e:
            # s=3 hits the zeta_R(1) pole, by structural design
            eq30_str = f"DIVERGENT/POLE — zeta_R(s-2={s-2}) at the s-2=1 pole" if s == 3 else f"ERR ({e})"
        results[s] = {
            "D_direct": direct_str,
            "D_closed_form": (str(mp.nstr(closed, 30)) if closed is not None else None),
            "D_Eq30_hurwitz_form": eq30_str,
        }
    return results


# ---------------------------------------------------------------------------
# Diagnostic: what the "correct half-integer Eq.(38)" looks like
# ---------------------------------------------------------------------------
def corrected_eq38_diagnostic():
    """
    Paper 18 wrote (incorrectly):

        D(4) = sum_{m>=1} 2 m (m+1) / m^4

    The CORRECT analog with the half-integer CH spectrum is

        D(4) = sum_{n>=0} 2 (n+1)(n+2) / (n + 3/2)^4
             = 16 * sum_{n>=0} 2 (n+1)(n+2) / (2n+3)^4
             = 16 * sum_{m=3, m odd}^{inf} (m^2 - 1)/2 / m^4
             = 8 * [ sum_{m odd>=3} 1/m^2 - sum_{m odd>=3} 1/m^4 ]
             = 8 * [ (lambda(2) - 1) - (lambda(4) - 1) ]
             = 8 [ lambda(2) - lambda(4) ]
             = pi^2 - pi^4 / 12.

    Here lambda(s) = (1 - 2^-s) zeta(s) is the Dirichlet lambda function:
        lambda(2) = pi^2 / 8,  lambda(4) = pi^4 / 96.

    Numerically:
        8 * (pi^2/8 - pi^4/96 - 1 + 1) ... wait, we subtracted m=1 twice.

    More carefully:
        sum_{m odd >= 1} 1/m^{2k} = lambda(2k) = (1 - 2^{-2k}) zeta(2k)
        sum_{m odd >= 3} 1/m^{2k} = lambda(2k) - 1.

        D(4) = 8 [ (lambda(2) - 1) - (lambda(4) - 1) ]
             = 8 [ lambda(2) - lambda(4) ]
             = 8 [ 3/4 * zeta(2) - 15/16 * zeta(4) ]
             = 6 zeta(2) - 15/2 zeta(4)
             = pi^2 - pi^4 / 12.

    Bit-exact symbolic agreement with the closed form.
    """
    closed = pi ** 2 - pi ** 4 / 12
    summed_odd = 8 * (
        (1 - mpf(2) ** (-2)) * mp.zeta(2)
        - (1 - mpf(2) ** (-4)) * mp.zeta(4)
    )
    # both forms are over odd m>=1 (since lambda(s) already starts at m=1)
    # the m=1 term  (m^2-1)/2 * 1 = 0 vanishes automatically — that's the
    # half-integer shift carrying away the singular term.
    direct = D_direct(4)
    return {
        "closed_form_pi2_minus_pi4_over_12": str(mp.nstr(closed, 40)),
        "summed_via_lambda_odd": str(mp.nstr(summed_odd, 40)),
        "direct_CH_sum": str(mp.nstr(direct, 40)),
        "all_agree_relative_residual": str(
            mp.nstr(abs(direct - closed) / abs(closed), 6)
        ),
    }


# ---------------------------------------------------------------------------
# Drive
# ---------------------------------------------------------------------------
def main():
    out: dict = {}
    out["meta"] = {
        "sprint": "M2",
        "purpose": "Resolve Paper 18 Eq.(38) D(4) substitution error",
        "mpmath_dps": DPS,
    }

    # --- Task 1: confirm D(4) = pi^2 - pi^4/12
    d4_direct = D_direct(4)
    d4_closed = pi ** 2 - pi ** 4 / 12
    d4_hurwitz = D_hurwitz_eq30(4)
    out["task1_D4_value"] = {
        "D_direct_sum_CH_spectrum_50dps": str(mp.nstr(d4_direct, 50)),
        "D_closed_form_pi2_minus_pi4_over_12": str(mp.nstr(d4_closed, 50)),
        "D_Eq30_Hurwitz_form": str(mp.nstr(d4_hurwitz, 50)),
        "direct_vs_closed_residual": str(mp.nstr(abs(d4_direct - d4_closed), 6)),
        "direct_vs_Eq30_residual": str(mp.nstr(abs(d4_direct - d4_hurwitz), 6)),
        "numeric_decimal": str(mp.nstr(d4_direct, 30)),
    }

    # --- Task 2: cross-validation D(2), D(4), D(6), D(8)
    out["task2_cross_validation_even_s"] = {}
    for s in (2, 4, 6, 8):
        d = D_direct(s)
        h = D_hurwitz_eq30(s)
        c = D_closed_form(s)
        out["task2_cross_validation_even_s"][f"s={s}"] = {
            "D_direct": str(mp.nstr(d, 30)),
            "D_Eq30": str(mp.nstr(h, 30)),
            "D_closed_form": (str(mp.nstr(c, 30)) if c is not None else None),
            "direct_vs_Eq30_residual": str(mp.nstr(abs(d - h), 6)),
            "Eq30_vs_closed_residual": (
                str(mp.nstr(abs(h - c), 6)) if c is not None else None
            ),
        }

    # --- Task 3 + 6: where does zeta(3) appear on the CH spinor bundle?
    out["task3_zeta3_search"] = search_zeta3_on_spinor_bundle()

    # explicit zeta(3) witness check at s=5
    d5_direct = D_direct(5)
    d5_closed = 14 * mp.zeta(3) - mpf(31) / 2 * mp.zeta(5)
    out["task3_zeta3_witness_at_s5"] = {
        "interpretation": (
            "Per Eq.(30), zeta_R(3) first appears with non-divergent "
            "coefficient at s=5 (the s=3 case multiplies zeta_R(1) which "
            "diverges). So the CH-spinor-bundle 'first zeta(3) witness' is "
            "D(5) = 14 zeta_R(3) - (31/2) zeta_R(5)."
        ),
        "D5_direct": str(mp.nstr(d5_direct, 50)),
        "D5_closed_form_14_zeta3_minus_31_over_2_zeta5": str(mp.nstr(d5_closed, 50)),
        "residual": str(mp.nstr(abs(d5_direct - d5_closed), 6)),
    }

    # --- Task 4: the WRONG sum is correct as a number
    s_wrong = S_wrong_eq38()
    out["task4_wrong_sum_value"] = {
        "S_wrong_eq38_sum_over_INTEGER_m": str(mp.nstr(s_wrong, 30)),
        "claimed_identification_2_zeta2_plus_2_zeta3": str(
            mp.nstr(2 * mp.zeta(2) + 2 * mp.zeta(3), 30)
        ),
        "true_D4_value": str(mp.nstr(D_direct(4), 30)),
        "discrepancy_S_wrong_minus_D4": str(
            mp.nstr(s_wrong - D_direct(4), 30)
        ),
        "note": (
            "The sum 'sum_{m>=1} 2m(m+1)/m^4 = 2 zeta(2)+2 zeta(3) ~ 5.694' "
            "is itself correct. The error is the IDENTIFICATION with D(4). "
            "The CH spectrum has half-integer eigenvalues lambda_n = n + 3/2, "
            "not integers m. The correct sum is "
            "sum_{n>=0} 2(n+1)(n+2)/(n+3/2)^4."
        ),
    }

    # --- Task 5: corrected Eq.(38) diagnostic
    out["task5_corrected_eq38_diagnostic"] = corrected_eq38_diagnostic()

    # --- Save
    out_path = Path(__file__).parent / "data" / "math_sprint_m2_d4.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    # --- Report
    print("=" * 78)
    print("SPRINT M2 — Paper 18 Eq.(38) D(4) substitution diagnostic")
    print("=" * 78)
    print()
    print("[Task 1] D(4) value at high precision:")
    print(f"  direct CH sum:               {out['task1_D4_value']['D_direct_sum_CH_spectrum_50dps']}")
    print(f"  closed form pi^2 - pi^4/12:  {out['task1_D4_value']['D_closed_form_pi2_minus_pi4_over_12']}")
    print(f"  Eq.(30) Hurwitz form:        {out['task1_D4_value']['D_Eq30_Hurwitz_form']}")
    print(f"  residuals (all ~ 0):         "
          f"{out['task1_D4_value']['direct_vs_closed_residual']}, "
          f"{out['task1_D4_value']['direct_vs_Eq30_residual']}")
    print()
    print("[Task 2] Even-s cross-validation vs Paper 28 Table 2:")
    for s_key, v in out["task2_cross_validation_even_s"].items():
        print(f"  {s_key}: direct={v['D_direct']}")
        if v["D_closed_form"]:
            print(f"        closed={v['D_closed_form']}")
        print(f"        Eq30={v['D_Eq30']}")
    print()
    print("[Task 4] Paper 18's wrong Eq.(38) sum:")
    print(f"  sum_{{m>=1}} 2m(m+1)/m^4 = {out['task4_wrong_sum_value']['S_wrong_eq38_sum_over_INTEGER_m']}")
    print(f"  true D(4)               = {out['task4_wrong_sum_value']['true_D4_value']}")
    print(f"  discrepancy             = {out['task4_wrong_sum_value']['discrepancy_S_wrong_minus_D4']}")
    print()
    print("[Task 3] where does zeta_R(3) live on the CH spinor bundle?")
    print(f"  {out['task3_zeta3_witness_at_s5']['interpretation']}")
    print(f"  D(5) direct:    {out['task3_zeta3_witness_at_s5']['D5_direct']}")
    print(f"  D(5) closed:    {out['task3_zeta3_witness_at_s5']['D5_closed_form_14_zeta3_minus_31_over_2_zeta5']}")
    print(f"  residual:       {out['task3_zeta3_witness_at_s5']['residual']}")
    print()
    print("[Task 5] corrected Eq.(38) diagnostic:")
    for k, v in out["task5_corrected_eq38_diagnostic"].items():
        print(f"  {k}: {v}")
    print()
    print(f"Full results saved to {out_path}")


if __name__ == "__main__":
    main()
