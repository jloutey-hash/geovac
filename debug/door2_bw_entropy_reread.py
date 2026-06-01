"""Door 2 (forcing-catalogue): re-read the BW wedge entropy.

The BW wedge entropy was filed in CLAUDE.md s3 as a NEGATIVE
("area-law REJECTED, R^2=0.83").  The forcing catalogue (docs/forcing_catalogue.md,
Door 2) suspects a positive Cardy-Calabrese boundary-log was shelved underneath.

This driver does the DISCIPLINED re-read:

  (1) The entropy is a CLOSED-FORM analytical sum (verified bit-exact vs the
      modular_hamiltonian numerics in the original sprint).  We therefore work
      with the closed form directly at HIGH PRECISION (mpmath), and extend the
      panel as far as tractable (n_max up to 2000+).

  (2) AUDIT the slope.  We do NOT just fit a + b*log(n).  We:
        - compute the EXACT analytical asymptotic limit of dS/d(log n)
          (the true slope) using the closed form,
        - run windowed slope fits to expose drift,
        - report the local logarithmic derivative n * dS/dn directly.
      Question: does the slope go to EXACTLY 2, to 1.963 (the original fit),
      or drift?

  (3) PSLQ the additive constant (the n->infinity intercept of S - 2*log n)
      against {log 2, log of small integers, gamma_E, log pi, ...}.

  (4) State / rule out the F-theorem (Paper 50) link precisely.

Closed form
===========
Eigenvalue k = 2m+1, m = 0..n-1.  Degeneracy g_m = (n-m)(n-m+1).
Per-state Boltzmann weight w_m = e^{-(2m+1)}.
Z = sum_m g_m w_m.
Per-state probability p_m = w_m / Z.
S = -sum_m g_m p_m log p_m = log Z + <K>,   <K> = sum_m g_m (2m+1) w_m / Z.

Asymptotic (n -> inf): the *dominant* shell is m=0 (k=1, equator),
g_0 = n(n+1) ~ n^2, but the geometric tail e^{-2} per shell is O(1), NOT
vanishing -- so the constant is a genuine series limit, not just log(n^2).
We extract it exactly below.
"""
from __future__ import annotations

import json
import os
import sys
from typing import Dict, List

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

mp.mp.dps = 80  # 80 decimal digits


# ---------------------------------------------------------------------------
# Closed-form entropy at arbitrary n_max, high precision
# ---------------------------------------------------------------------------

def entropy_closed_form(n_max: int) -> mp.mpf:
    """S(rho_W) at canonical BW (beta=2pi baked into K_alpha) via closed form.

    S = log Z + <K>, all at mp precision.
    """
    n = n_max
    Z = mp.mpf(0)
    K_weighted = mp.mpf(0)
    for m in range(n):
        k = 2 * m + 1
        g = mp.mpf((n - m) * (n - m + 1))
        w = mp.e ** (-k)
        gw = g * w
        Z += gw
        K_weighted += gw * k
    avgK = K_weighted / Z
    S = mp.log(Z) + avgK
    return S


def entropy_components(n_max: int) -> Dict[str, mp.mpf]:
    n = n_max
    Z = mp.mpf(0)
    K_weighted = mp.mpf(0)
    for m in range(n):
        k = 2 * m + 1
        g = mp.mpf((n - m) * (n - m + 1))
        w = mp.e ** (-k)
        gw = g * w
        Z += gw
        K_weighted += gw * k
    avgK = K_weighted / Z
    return {"Z": Z, "avgK": avgK, "logZ": mp.log(Z), "S": mp.log(Z) + avgK}


# ---------------------------------------------------------------------------
# Exact asymptotic limit of (S - 2*log n)  and  of the slope dS/dlog n
# ---------------------------------------------------------------------------

def constant_limit_asymptotic(n_terms: int = 400, n_anchor: int = 20000) -> Dict:
    """Extract the exact n->inf constant  C_inf = lim (S(n) - 2 log n).

    Strategy: write the shell index as j = n - m  (j = 1..n), so that the
    DOMINANT shell is j=n? No -- dominant is m=0 i.e. j=n.  Re-anchor on m
    (the small index), which is where the weight concentrates.

    g_m = (n-m)(n-m+1).  For fixed m and n->inf:
        g_m / n^2 -> 1   (since (n-m)(n-m+1)/n^2 -> 1).
    More precisely g_m = n^2 + (1-2m) n + m(m-1).

    Z = sum_m g_m e^{-(2m+1)}
      = n^2 * A0 + n * A1 + A2,    where
        A0 = sum_{m>=0} e^{-(2m+1)}                = e^{-1}/(1-e^{-2})
        A1 = sum_{m>=0} (1-2m) e^{-(2m+1)}
        A2 = sum_{m>=0} m(m-1) e^{-(2m+1)}
    (the upper limit n-1 -> inf is harmless; tail is e^{-2n}).

    <K> = sum_m g_m (2m+1) e^{-(2m+1)} / Z
        = (n^2 B0 + n B1 + B2) / Z,   B's analogous with extra (2m+1) factor.

    Then S = log Z + <K>.  As n->inf:
        log Z = 2 log n + log A0 + (A1/A0)/n + O(1/n^2)
        <K>   = B0/A0 + O(1/n)
    so
        C_inf = log(A0) + B0/A0      [exact closed form in e]
        slope = dS/d(log n) -> 2     [exactly, from the n^2 leading term]

    We compute A0,B0 in closed form and C_inf exactly.
    """
    x = mp.e ** (-2)  # ratio between successive shells in the per-state weight squared sense

    # A0 = sum_{m>=0} e^{-(2m+1)} = e^{-1} sum x^m = e^{-1}/(1-x)
    A0 = mp.e ** (-1) / (1 - x)
    # B0 = sum_{m>=0} (2m+1) e^{-(2m+1)} = e^{-1} sum (2m+1) x^m
    #    sum (2m+1) x^m = (1+x)/(1-x)^2
    B0 = mp.e ** (-1) * (1 + x) / (1 - x) ** 2

    C_inf = mp.log(A0) + B0 / A0

    # Cross-check B0/A0:  this is the limiting <K>.
    avgK_inf = B0 / A0

    # Also compute next-order so we can predict S(n) - (2 log n + C_inf) ~ c1/n
    # A1 = sum (1-2m) e^{-(2m+1)} = e^{-1}[ sum x^m - 2 sum m x^m ]
    #    sum x^m = 1/(1-x);  sum m x^m = x/(1-x)^2
    A1 = mp.e ** (-1) * (1 / (1 - x) - 2 * x / (1 - x) ** 2)
    # slope-correction term in log Z is (A1/A0)/n
    logZ_subleading_coeff = A1 / A0  # coefficient of 1/n in log Z expansion

    return {
        "x_=_e^-2": x,
        "A0": A0,
        "B0": B0,
        "A1": A1,
        "avgK_inf_=_B0/A0": avgK_inf,
        "C_inf_=_log(A0)+B0/A0": C_inf,
        "logZ_subleading_coeff_(1/n)": logZ_subleading_coeff,
    }


# ---------------------------------------------------------------------------
# Windowed slope audit (the curve-fit discipline)
# ---------------------------------------------------------------------------

def windowed_slopes(panel: List[Dict], window: int = 5) -> List[Dict]:
    """Local 2-point and windowed slope dS/dlog(n) across sliding windows."""
    out = []
    for i in range(len(panel) - 1):
        n0, n1 = panel[i]["n_max"], panel[i + 1]["n_max"]
        S0, S1 = mp.mpf(panel[i]["S"]), mp.mpf(panel[i + 1]["S"])
        slope2 = (S1 - S0) / (mp.log(n1) - mp.log(n0))
        out.append({
            "n_lo": n0, "n_hi": n1,
            "two_point_slope": float(slope2),
        })
    return out


def least_squares_loglinear(ns: List[int], Ss: List[mp.mpf]) -> Dict:
    """Fit S = a*log(n) + b at mp precision; return a, b, R2."""
    L = [mp.log(n) for n in ns]
    y = list(Ss)
    N = len(ns)
    sumL = sum(L)
    sumy = sum(y)
    sumLL = sum(l * l for l in L)
    sumLy = sum(l * yy for l, yy in zip(L, y))
    denom = N * sumLL - sumL * sumL
    a = (N * sumLy - sumL * sumy) / denom
    b = (sumy * sumLL - sumL * sumLy) / denom
    ybar = sumy / N
    ss_tot = sum((yy - ybar) ** 2 for yy in y)
    ss_res = sum((yy - (a * l + b)) ** 2 for yy, l in zip(y, L))
    R2 = 1 - ss_res / ss_tot
    return {"slope": float(a), "intercept": float(b), "R2": float(R2)}


# ---------------------------------------------------------------------------
# PSLQ the constant
# ---------------------------------------------------------------------------

def pslq_constant(C_inf: mp.mpf) -> Dict:
    """Try to identify C_inf against a basis of natural constants."""
    results = {}

    # Basis of candidate transcendentals / logs
    basis = {
        "1": mp.mpf(1),
        "log2": mp.log(2),
        "log3": mp.log(3),
        "log5": mp.log(5),
        "logpi": mp.log(mp.pi),
        "log(1-e^-2)": mp.log(1 - mp.e ** -2),
        "log(e^2-1)": mp.log(mp.e ** 2 - 1),
        "gamma_E": mp.euler,
        "1/(e^2-1)": 1 / (mp.e ** 2 - 1),
        "(e^2+1)/(e^2-1)": (mp.e ** 2 + 1) / (mp.e ** 2 - 1),
        "C_inf": C_inf,
    }

    # First: report C_inf to high precision and its naive value
    results["C_inf_value"] = mp.nstr(C_inf, 50)

    # The analytic form is C_inf = log(A0) + B0/A0.
    # A0 = e^{-1}/(1-e^-2) = 1/(e - e^{-1}) = 1/(2 sinh 1).
    A0 = mp.e ** (-1) / (1 - mp.e ** -2)
    results["A0_value"] = mp.nstr(A0, 50)
    results["A0_=_1/(2 sinh 1)?"] = mp.nstr(A0 - 1 / (2 * mp.sinh(1)), 5)
    results["log_A0"] = mp.nstr(mp.log(A0), 50)
    results["-log(2 sinh 1)"] = mp.nstr(-mp.log(2 * mp.sinh(1)), 50)

    # B0/A0 = (1+x)/(1-x) with x=e^-2 = coth(1)
    x = mp.e ** -2
    B0_over_A0 = (1 + x) / (1 - x)
    results["B0/A0_value"] = mp.nstr(B0_over_A0, 50)
    results["B0/A0_=_coth(1)?"] = mp.nstr(B0_over_A0 - mp.coth(1), 5)

    # So C_inf = -log(2 sinh 1) + coth(1).  Verify.
    C_inf_analytic = -mp.log(2 * mp.sinh(1)) + mp.coth(1)
    results["C_inf_analytic_=_coth(1)-log(2 sinh 1)"] = mp.nstr(C_inf_analytic, 50)
    results["residual_C_inf_vs_analytic"] = mp.nstr(C_inf - C_inf_analytic, 5)

    # Now try to see whether this closed form reduces to anything simpler
    # via PSLQ against {1, log2, coth1, log(sinh1), log(e^2-1), ...}
    pslq_vec = [
        C_inf,
        mp.mpf(1),
        mp.log(2),
        mp.coth(1),
        mp.log(mp.sinh(1)),
        mp.log(mp.e ** 2 - 1),
    ]
    rel = mp.pslq(pslq_vec, maxcoeff=10 ** 6, maxsteps=10 ** 5)
    results["pslq_basis"] = ["C_inf", "1", "log2", "coth1", "log(sinh1)", "log(e^2-1)"]
    results["pslq_relation"] = rel

    return results


# ---------------------------------------------------------------------------
# Cardy-Calabrese reading
# ---------------------------------------------------------------------------

def cardy_calabrese_reading(slope_limit: float) -> Dict:
    """A Cardy-Calabrese / boundary-CFT entanglement log has the form
       S = (c_eff / 6) * log(L)  (periodic)  or  (c_eff/6) log(L) (one interval),
    with c_eff the (effective) central charge and L the subsystem size.

    Here the "size" proxy is n_max (the cutoff = number of shells), and the
    slope we measure is dS/d(log n_max).  A literal Cardy-Calabrese identity
    would require this slope to equal c_eff/6 (open) or c_eff/3 (one interval
    in an infinite system) for some physically meaningful c_eff.

    But the discrete construction here is NOT a 1+1D CFT interval; it is a
    degeneracy-log of the lowest K_alpha shell on a (2+1)-D wedge.  The slope
    is the DIMENSION of the equator boundary, not a central charge.  We state
    both readings and let the number decide.
    """
    return {
        "measured_slope_limit": slope_limit,
        "interpretation_A_dim_boundary": (
            "slope = dim(partial W) = dim(S^2 equator) = 2; the entropy is "
            "log of the equator-shell degeneracy n_eq ~ n_max^2, i.e. "
            "log(n_max^2) = 2 log n_max.  This is a DEGENERACY-LOG, not a "
            "1+1D Cardy-Calabrese central-charge log."
        ),
        "interpretation_B_cardy_calabrese": (
            "If forced into Cardy-Calabrese S=(c/6) log L with L=n_max, then "
            f"c_eff = 6 * slope = {6*slope_limit:.4f}.  This c_eff has no "
            "independent derivation and is not the F-theorem F-coefficient "
            "(which is a 3D free energy, not a 2D central charge).  The "
            "Cardy-Calabrese reading is therefore a FORCED analogy, not an "
            "identity."
        ),
        "verdict": (
            "The slope is the equator-boundary DIMENSION (=2), arising as a "
            "degeneracy log.  It is NOT a Cardy-Calabrese central-charge "
            "coefficient.  The '2' is dim(S^2), exact, and provable from the "
            "n^2 leading degeneracy -- but it is NOT c_eff/6 of a boundary CFT."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 74)
    print("Door 2: BW wedge entropy re-read (disciplined audit)")
    print("=" * 74)

    # ---- (1) Extended panel at high precision via closed form ----
    panel_ns = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 80,
                120, 200, 350, 600, 1000, 2000]
    panel = []
    print(f"\n{'n_max':>6} {'S':>22} {'S/log n':>14} {'S-2log n':>20}")
    for n in panel_ns:
        comp = entropy_components(n)
        S = comp["S"]
        s_over_logn = S / mp.log(n)
        s_minus_2logn = S - 2 * mp.log(n)
        panel.append({
            "n_max": n,
            "S": mp.nstr(S, 30),
            "S_float": float(S),
            "logZ": float(comp["logZ"]),
            "avgK": float(comp["avgK"]),
            "S_over_logn": float(s_over_logn),
            "S_minus_2logn": mp.nstr(s_minus_2logn, 30),
            "S_minus_2logn_float": float(s_minus_2logn),
        })
        print(f"{n:6d} {mp.nstr(S,18):>22} {float(s_over_logn):14.8f} "
              f"{mp.nstr(s_minus_2logn,16):>20}")

    # ---- (2) Exact asymptotic ----
    print("\n" + "=" * 74)
    print("(2) EXACT analytical asymptotic (NOT a fit)")
    print("=" * 74)
    asy = constant_limit_asymptotic()
    for k, v in asy.items():
        print(f"  {k:32s} = {mp.nstr(v, 30) if isinstance(v, mp.mpf) else v}")

    C_inf = asy["C_inf_=_log(A0)+B0/A0"]
    print(f"\n  ==> S(n_max) -> 2*log(n_max) + C_inf")
    print(f"      C_inf = {mp.nstr(C_inf, 40)}")
    print(f"      (original fit intercept was 0.540 / 0.56 over small n)")

    # Verify S(n) - 2 log n approaches C_inf
    print("\n  Convergence of [S(n) - 2 log n] -> C_inf:")
    for n in [10, 50, 200, 1000, 2000]:
        comp = entropy_components(n)
        diff = comp["S"] - 2 * mp.log(n) - C_inf
        print(f"    n={n:5d}: S-2logn-C_inf = {mp.nstr(diff, 6):>14} "
              f"(predicted ~ {float(asy['logZ_subleading_coeff_(1/n)'])/n:+.6e}/... )")

    # ---- (3) Slope audit ----
    print("\n" + "=" * 74)
    print("(3) SLOPE AUDIT (windowed + local log-derivative)")
    print("=" * 74)
    Ss = [mp.mpf(p["S"]) for p in panel]
    ns = [p["n_max"] for p in panel]

    print("\n  Two-point local slopes dS/dlog(n):")
    wslopes = windowed_slopes(panel)
    for w in wslopes:
        print(f"    n in [{w['n_lo']:5d},{w['n_hi']:5d}]: slope = {w['two_point_slope']:.6f}")

    # Local logarithmic derivative n dS/dn via finite difference of closed form
    print("\n  Local log-derivative  n * dS/dn  (closed-form FD, h=1):")
    for n in [3, 5, 10, 20, 50, 100, 300, 1000]:
        Sp = entropy_closed_form(n + 1)
        Sm = entropy_closed_form(n - 1)
        dSdlogn = (Sp - Sm) / (mp.log(n + 1) - mp.log(n - 1))
        print(f"    n={n:5d}: dS/dlog(n) = {mp.nstr(dSdlogn, 12)}")

    # Windowed least-squares slope on disjoint windows (drift exposure)
    print("\n  Windowed least-squares slope a in S=a log n + b:")
    windows = [
        ("orig [2..7]", [2, 3, 4, 5, 6, 7]),
        ("[8..20]", [8, 9, 10, 12, 15, 20]),
        ("[30..200]", [30, 50, 80, 120, 200]),
        ("[350..2000]", [350, 600, 1000, 2000]),
        ("[200..2000]", [200, 350, 600, 1000, 2000]),
    ]
    window_fits = {}
    for label, wns in windows:
        wSs = [entropy_closed_form(n) for n in wns]
        fit = least_squares_loglinear(wns, wSs)
        window_fits[label] = fit
        print(f"    {label:14s}: slope={fit['slope']:.6f}  "
              f"intercept={fit['intercept']:.6f}  R2={fit['R2']:.8f}")

    print("\n  ==> The fitted slope DRIFTS toward 2 as the window moves up.")
    print("      The original 1.963 is a SMALL-n artifact; the true slope is 2.")

    # ---- (4) Constant identification ----
    print("\n" + "=" * 74)
    print("(4) CONSTANT IDENTIFICATION (PSLQ + closed form)")
    print("=" * 74)
    cid = pslq_constant(C_inf)
    for k, v in cid.items():
        print(f"  {k:42s} : {v}")

    # ---- (5) Cardy-Calabrese reading ----
    print("\n" + "=" * 74)
    print("(5) CARDY-CALABRESE vs DEGENERACY-LOG reading")
    print("=" * 74)
    cc = cardy_calabrese_reading(2.0)
    for k, v in cc.items():
        print(f"\n  [{k}]\n    {v}")

    # ---- Save ----
    out = {
        "description": (
            "Door 2 disciplined re-read of the BW wedge entropy. "
            "Closed-form S(n_max) at mp precision; exact asymptotic; "
            "slope audit; constant PSLQ; Cardy-Calabrese vs degeneracy-log."
        ),
        "panel": panel,
        "asymptotic": {k: (mp.nstr(v, 40) if isinstance(v, mp.mpf) else v)
                       for k, v in asy.items()},
        "C_inf": mp.nstr(C_inf, 50),
        "C_inf_closed_form": "coth(1) - log(2*sinh(1))",
        "windowed_two_point_slopes": wslopes,
        "windowed_lstsq_fits": window_fits,
        "constant_identification": {k: (mp.nstr(v, 40) if isinstance(v, mp.mpf) else v)
                                    for k, v in cid.items()},
        "cardy_calabrese_reading": cc,
    }
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data", "door2_bw_entropy_reread.json",
    )
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
