"""Sprint G1 — Path 1: Spectral action on S^3_R parameterized by radius.

Tests whether the framework's Connes-Chamseddine-style spectral action
Tr f(D^2/Lambda^2) on the family of round S^3 has gravity-like structure
analogous to CC's "spectral action principle picks out a preferred metric".

Setup
-----
Camporesi-Higuchi Dirac spectrum on S^3 of radius R:
  |lambda_n(R)| = (n + 3/2) / R
  g_n = 2 (n+1) (n+2)              (full 4-component Dirac degeneracy)

Spectral zeta: zeta_R(s) = R^{2s} zeta_unit(s) where
  zeta_unit(s) = sum_n 2(n+1)(n+2) (n+3/2)^{-2s}

Paper 28's two-term exactness theorem states:
  Tr e^{-tD^2}|_{unit} = (sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2} + O(exp(-pi^2/t))

The spectral-zeta translation:
  zeta_unit(s) has SIMPLE POLES at s = 3/2 (residue 1) and s = 1/2 (residue -1/4),
  and is HOLOMORPHIC elsewhere with zeta_unit(-k) = 0 for all k = 0, 1, 2, ...
  (substantive identity — verified numerically below).

Spectral action via Mellin/residues
------------------------------------
For cutoff f with Mellin transform phi(s) = int_0^infty f(x) x^{s-1} dx:

  S(R, Lambda) = Tr f(D^2/Lambda^2) = (1/2 pi i) int phi(s) (Lambda R)^{2s} zeta_unit(s) ds
              = phi(3/2) * (Lambda R)^3 * 1 + phi(1/2) * (Lambda R) * (-1/4)
              = phi(3/2) u^3 - (1/4) phi(1/2) u

(no corrections from phi-poles at non-positive integers because zeta_unit(-k) = 0).

Three cutoffs tested:
  Gaussian:   f(x) = e^{-x},    phi(s) = Gamma(s)
  Polynomial: f(x) = e^{-x^2},  phi(s) = (1/2) Gamma(s/2)
  Sharp:      f(x) = Theta(1-x), phi(s) = 1/s

Findings (verified below)
-------------------------
1. zeta_unit(-k) = 0 for k = 0..5 (exact rational vanishing).
2. S(R, Lambda) = A(f) (Lambda R)^3 - B(f) (Lambda R) is EXACT modulo
   exp-small corrections in (Lambda R) for any "nice" cutoff f.
3. A, B > 0 for all three cutoffs (opposite signs in S = A u^3 - B u).
4. Formal extremum u_crit = sqrt(B/(3A)) = O(1) for all three cutoffs.
5. At large u (UV regime), S_asymp matches S_QM_exact to machine precision.
6. At small u (IR regime), S_QM_exact is tiny (modes Boltzmann-suppressed)
   while S_asymp continues its A u^3 - B u shape into negative values.
7. The "formal extremum at u ~ 0.4-0.5" lies in the IR regime where the
   asymptotic does NOT faithfully represent the QM spectral action.

Reading
-------
The CC spectral action principle takes the *asymptotic expansion* as the
effective classical action for the geometry, regardless of whether the
extremum sits at finite u or in the deep IR. Under that reading, this
sprint identifies the GeoVac analog of CC's Einstein-Hilbert + cosmological
constant structure on S^3_R, with two-term EXACTNESS (no higher-curvature
corrections, ever).

The literal QM exact spectral action has no extremum — it is monotonically
increasing in u. So GeoVac reproduces CC's structural prediction but inherits
the same gap between "spectral action principle" and "literal extremum of the
QM exact sum" that CC has always lived with.
"""

import json
from pathlib import Path

import sympy as sp
from mpmath import mp, mpf, exp, log, gamma, fmul, fdiv, sqrt as msqrt, pi as mp_pi
from sympy import Rational, Integer, zeta as sym_zeta, bernoulli, simplify

mp.dps = 50

OUT_JSON = Path(__file__).parent / "data" / "g1_spectral_action_S3_radius.json"
OUT_JSON.parent.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Three cutoff functions and their Mellin transforms
# ---------------------------------------------------------------------------

def phi_gaussian(s):
    """f(x) = e^{-x},  phi(s) = Gamma(s)."""
    return gamma(s)


def phi_polynomial(s):
    """f(x) = e^{-x^2},  phi(s) = (1/2) Gamma(s/2)."""
    return mpf("0.5") * gamma(s / mpf("2"))


def phi_sharp(s):
    """f(x) = Theta(1 - x),  phi(s) = 1/s."""
    return mpf("1") / mpf(s)


def f_gaussian(x):
    return exp(-x)


def f_polynomial(x):
    return exp(-x * x)


def f_sharp(x):
    return mpf("1") if x <= mpf("1") else mpf("0")


CUTOFFS = {
    "gaussian":   {"phi": phi_gaussian,   "f": f_gaussian,   "label": "f(x) = e^{-x}"},
    "polynomial": {"phi": phi_polynomial, "f": f_polynomial, "label": "f(x) = e^{-x^2}"},
    "sharp":      {"phi": phi_sharp,      "f": f_sharp,      "label": "f(x) = Theta(1-x)"},
}

# ---------------------------------------------------------------------------
# Spectral zeta of Dirac on unit S^3 via Hurwitz zeta
# ---------------------------------------------------------------------------
# zeta_unit(s) = sum_n 2(n+1)(n+2) (n+3/2)^{-2s}
#             = sum_n [2(n+3/2)^2 - 1/2] (n+3/2)^{-2s}
#             = 2 zeta_H(2s-2, 3/2) - (1/2) zeta_H(2s, 3/2)
#
# At s = -k (k = 0, 1, 2, ...):
#   zeta_H(-n, a) = -B_{n+1}(a) / (n+1)   for non-negative integer n
#
# So zeta_unit(-k) = 2 * [-B_{2k+3-1+1}(3/2)/(2k+3-1+1)] - (1/2) * [-B_{2k+1}(3/2)/(2k+1)]
#                 = -2 B_{-2k-1}(3/2)/(?)... let me redo via the standard relation.
#
# Standard: zeta_H(s, a) at s = -n is -B_{n+1}(a)/(n+1).
# Here 2s-2 at s=-k is -2k-2, so zeta_H(-2k-2, 3/2) = -B_{2k+3}(3/2)/(2k+3).
# And 2s at s=-k is -2k, so zeta_H(-2k, 3/2) = -B_{2k+1}(3/2)/(2k+1).
#
# zeta_unit(-k) = 2 * [-B_{2k+3}(3/2)/(2k+3)] - (1/2) * [-B_{2k+1}(3/2)/(2k+1)]
#              = -2 B_{2k+3}(3/2)/(2k+3) + B_{2k+1}(3/2)/(2 (2k+1))
#
# Claim: this is identically zero, i.e.
#   4 (2k+1) B_{2k+3}(3/2) = (2k+3) B_{2k+1}(3/2).
#
# Verify symbolically for k = 0..5.

def zeta_unit_at_neg_k(k):
    """Symbolic value of zeta_{D^2,unit}(-k) for non-negative integer k.

    Uses zeta_unit(s) = 2 zeta_H(2s-2, 3/2) - (1/2) zeta_H(2s, 3/2)
    and zeta_H(-n, a) = -B_{n+1}(a)/(n+1) for non-negative integer n.
    """
    from sympy import bernoulli, Rational, simplify, symbols
    x = symbols("x")
    # B_{2k+3}(3/2) and B_{2k+1}(3/2)
    # bernoulli(n, x) returns the Bernoulli polynomial B_n(x)
    B_top  = bernoulli(2*k + 3, Rational(3, 2))
    B_bot  = bernoulli(2*k + 1, Rational(3, 2))
    val = -2 * B_top / (2*k + 3) + B_bot / (2 * (2*k + 1))
    return simplify(val)


def zeta_unit_at_pole_residues():
    """Residues of zeta_unit at s = 3/2 and s = 1/2.

    Returns symbolic (res_3/2, res_1/2).

    From two-term exactness:
      Tr e^{-tD^2}|_unit = (sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2} + exp small

    Mellin connection: pole of zeta at s = a contributes to Tr e^{-tD^2}
    a term Gamma(a) * Res(zeta, a) * t^{-a}.

    So:
      Gamma(3/2) * Res(zeta, 3/2) = sqrt(pi)/2  => Res = (sqrt(pi)/2) / (sqrt(pi)/2) = 1
      Gamma(1/2) * Res(zeta, 1/2) = -sqrt(pi)/4 => Res = -sqrt(pi)/4 / sqrt(pi) = -1/4
    """
    return Rational(1, 1), Rational(-1, 4)


# ---------------------------------------------------------------------------
# Asymptotic spectral action
# ---------------------------------------------------------------------------

def asymptotic_action(u, phi_fn):
    """S_asymp(u) = phi(3/2) u^3 - (1/4) phi(1/2) u"""
    A = phi_fn(mpf("1.5"))                  # coefficient of u^3
    B = mpf("0.25") * phi_fn(mpf("0.5"))    # coefficient of u
    return A * u**3 - B * u, A, B


def asymptotic_extremum_u(phi_fn):
    """u_crit = sqrt(B / (3A)) where A = phi(3/2), B = (1/4) phi(1/2)."""
    A = phi_fn(mpf("1.5"))
    B = mpf("0.25") * phi_fn(mpf("0.5"))
    return msqrt(B / (3 * A)), A, B


# ---------------------------------------------------------------------------
# QM exact spectral action (truncated discrete sum)
# ---------------------------------------------------------------------------

def qm_exact_action(u, f_fn, n_max=200):
    """S_QM(u) = sum_{n=0}^{n_max} 2(n+1)(n+2) f((n+3/2)^2 / u^2)"""
    total = mpf("0")
    u_sq = u * u
    for n in range(n_max + 1):
        lam = mpf(n) + mpf("1.5")
        g_n = 2 * (n + 1) * (n + 2)
        x = (lam * lam) / u_sq
        total += g_n * f_fn(x)
    return total


def qm_exact_convergence(u, f_fn, n_levels=(50, 100, 200, 400)):
    """Verify convergence in n_max."""
    vals = []
    for nm in n_levels:
        vals.append(qm_exact_action(u, f_fn, n_max=nm))
    return vals


# ---------------------------------------------------------------------------
# Main sprint
# ---------------------------------------------------------------------------

def main():
    results = {}

    print("=" * 72)
    print("Sprint G1 -- Path 1: Spectral action on S^3_R")
    print("Tr f(D^2/Lambda^2) on the family of round S^3 parameterized by R")
    print("=" * 72)

    # -----------------------------------------------------------------------
    # Step 1: Verify zeta_unit(-k) = 0 for k = 0..5 (symbolic).
    # This is the structural identity making the two-term spectral action EXACT.
    # -----------------------------------------------------------------------
    print("\n[Step 1] Verifying zeta_unit(-k) = 0 for k = 0..5 (symbolic exact):")
    zeta_neg_results = {}
    for k in range(6):
        val = zeta_unit_at_neg_k(k)
        zeta_neg_results[k] = str(val)
        status = "OK" if val == 0 else "FAIL"
        print(f"  k={k}:  zeta_unit({-k}) = {val}   [{status}]")
    results["zeta_neg_k"] = zeta_neg_results

    # -----------------------------------------------------------------------
    # Step 2: Pole residues of zeta_unit.
    # -----------------------------------------------------------------------
    print("\n[Step 2] Residues of zeta_unit at the two poles:")
    r32, r12 = zeta_unit_at_pole_residues()
    print(f"  Res(zeta_unit, 3/2) = {r32}")
    print(f"  Res(zeta_unit, 1/2) = {r12}")
    results["pole_residues"] = {"3/2": str(r32), "1/2": str(r12)}

    # -----------------------------------------------------------------------
    # Step 3: For each cutoff, compute asymptotic coefficients and extremum.
    # -----------------------------------------------------------------------
    print("\n[Step 3] Asymptotic spectral action coefficients per cutoff:")
    print(f"  S_asymp(u) = A u^3 - B u  where  A = phi(3/2), B = (1/4) phi(1/2)")
    cutoff_data = {}
    for name, info in CUTOFFS.items():
        u_crit, A, B = asymptotic_extremum_u(info["phi"])
        S_crit, _, _ = asymptotic_action(u_crit, info["phi"])
        print(f"\n  {name} ({info['label']}):")
        print(f"    A = phi(3/2) = {sp.nsimplify(float(A), rational=False)}  ~ {float(A):.6f}")
        print(f"    B = (1/4) phi(1/2) = {sp.nsimplify(float(B), rational=False)}  ~ {float(B):.6f}")
        print(f"    u_crit = sqrt(B/(3A)) = {float(u_crit):.6f}")
        print(f"    S_asymp(u_crit) = {float(S_crit):.6f}  (the formal minimum value)")
        cutoff_data[name] = {
            "A": float(A),
            "B": float(B),
            "u_crit": float(u_crit),
            "S_asymp_at_u_crit": float(S_crit),
        }

    # -----------------------------------------------------------------------
    # Step 4: Compare S_QM_exact vs S_asymp across u range.
    # -----------------------------------------------------------------------
    print("\n[Step 4] S_QM_exact vs S_asymp across u in [0.3, 20] for each cutoff:")
    u_grid = [mpf(x) for x in
              ["0.3", "0.5", "0.7", "1.0", "1.5", "2.0", "3.0", "5.0",
               "7.0", "10.0", "15.0", "20.0"]]
    comparison = {}
    for name, info in CUTOFFS.items():
        phi_fn = info["phi"]
        f_fn = info["f"]
        rows = []
        print(f"\n  {name}:")
        print(f"    {'u':>7s}  {'S_QM':>15s}  {'S_asymp':>15s}  {'rel diff':>10s}")
        for u in u_grid:
            S_qm = qm_exact_action(u, f_fn, n_max=300)
            S_as, _, _ = asymptotic_action(u, phi_fn)
            denom = max(abs(S_qm), mpf("1e-30"))
            rel = abs(S_qm - S_as) / denom
            print(f"    {float(u):>7.3f}  {float(S_qm):>+15.6e}  {float(S_as):>+15.6e}  {float(rel):>10.3e}")
            rows.append({
                "u": float(u),
                "S_QM": float(S_qm),
                "S_asymp": float(S_as),
                "rel_diff": float(rel),
            })
        comparison[name] = rows

    results["cutoff_data"] = cutoff_data
    results["comparison"] = comparison

    # -----------------------------------------------------------------------
    # Step 5: Convergence of QM sum in n_max at u = 3 and u = 10.
    # -----------------------------------------------------------------------
    print("\n[Step 5] Convergence of QM sum in n_max (Gaussian cutoff):")
    convergence = {}
    for u_val in [mpf("3.0"), mpf("10.0")]:
        vals = qm_exact_convergence(u_val, f_gaussian)
        n_levels = (50, 100, 200, 400)
        print(f"\n  u = {float(u_val)}:")
        for nm, v in zip(n_levels, vals):
            print(f"    n_max = {nm}:  S_QM = {float(v):.10e}")
        convergence[f"u={float(u_val)}"] = {
            f"n_max={nm}": float(v) for nm, v in zip(n_levels, vals)
        }
    results["convergence"] = convergence

    # -----------------------------------------------------------------------
    # Step 6: Locate the asymptotic-vs-exact regime boundary.
    # At small u, S_QM -> 0 (modes exponentially suppressed); S_asymp -> 0 too
    # but only linearly (since S_asymp = A u^3 - B u -> -B u).
    # The "agreement" criterion: |S_QM - S_asymp| / |S_asymp| < 0.01.
    # -----------------------------------------------------------------------
    print("\n[Step 6] Regime boundary -- where does S_asymp become a good approximation?")
    regime = {}
    for name, info in CUTOFFS.items():
        phi_fn = info["phi"]
        f_fn = info["f"]
        # Scan u from 0.3 to 20, find smallest u where |S_QM - S_asymp|/|S_asymp| < 0.01
        boundary = None
        for u in [mpf(x) for x in ["0.3", "0.5", "0.7", "1.0", "1.5", "2.0", "2.5",
                                    "3.0", "4.0", "5.0", "7.0", "10.0"]]:
            S_qm = qm_exact_action(u, f_fn, n_max=300)
            S_as, _, _ = asymptotic_action(u, phi_fn)
            if abs(S_as) < mpf("1e-10"):
                continue
            rel = abs(S_qm - S_as) / abs(S_as)
            if rel < mpf("0.01"):
                boundary = float(u)
                break
        regime[name] = {
            "u_agreement": boundary,
            "u_crit_extremum": cutoff_data[name]["u_crit"],
            "extremum_in_agreement_regime": (boundary is not None and
                                              boundary <= cutoff_data[name]["u_crit"]),
        }
        print(f"  {name}: asymptotic agrees with QM for u >= {boundary}")
        print(f"           u_crit (formal extremum) = {cutoff_data[name]['u_crit']:.4f}")
        if boundary is not None and boundary > cutoff_data[name]["u_crit"]:
            print(f"           -> extremum is OUTSIDE asymptotic regime")
        else:
            print(f"           -> extremum is INSIDE asymptotic regime")
    results["regime"] = regime

    # -----------------------------------------------------------------------
    # Save JSON
    # -----------------------------------------------------------------------
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)
    return results


if __name__ == "__main__":
    main()
