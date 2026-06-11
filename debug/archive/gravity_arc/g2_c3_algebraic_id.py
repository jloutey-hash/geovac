"""
g2_c3_algebraic_id.py -- Algebraic identification of curvature coefficient c3.

Three-pronged attack:
  1. High-precision mpmath re-extraction from existing B(n_int) data
  2. Structural formula search: c3 as polynomial in Paper 2 invariants
  3. Extended PSLQ against carefully chosen bases

The key insight: c2 was identified structurally (c2 = (2 - BD - FD - F/B)/5),
not by PSLQ. We should try the same approach for c3.

T9 theorem constraint: c3 = r0 + r2*pi^2 + r4*pi^4  (no odd zeta at one loop)
"""

import json
import sys
import time
from pathlib import Path
from fractions import Fraction

import mpmath
import numpy as np

mpmath.mp.dps = 80

DATA_DIR = Path(__file__).parent / "data"

# --------------------------------------------------------------------------
# Constants at high precision
# --------------------------------------------------------------------------

ALPHA_CODATA = mpmath.mpf("7.2973525693e-3")
ALPHA_UNC = mpmath.mpf("1.1e-12")

PI = mpmath.pi
PI2 = PI ** 2
PI4 = PI ** 4
SCHWINGER = ALPHA_CODATA / (2 * PI)

B = mpmath.mpf(42)
F = PI2 / 6
DELTA = mpmath.mpf(1) / 40

BD = B * DELTA   # 21/20
FD = F * DELTA   # pi^2/240
FB = F / B       # pi^2/252

C1 = mpmath.mpf("0.5")

C2 = (2 - BD - FD - FB) / 5
C2_check = mpmath.mpf(19) / 100 - 41 * PI2 / 25200
assert abs(C2 - C2_check) < mpmath.mpf("1e-50")

LAM = mpmath.mpf(5) / 2
LAM2 = LAM ** 2
X = 1 / LAM2
X2 = X ** 2
X3 = X ** 3


def load_data():
    path = DATA_DIR / "g2_c3_investigation.json"
    with open(path) as f:
        data = json.load(f)
    all_B = {int(k): mpmath.mpf(str(v)) for k, v in data["all_B"].items()}
    V_mag = mpmath.mpf(str(data["V_mag"]))
    return all_B, V_mag


def fit_power_law_mp(all_B, n_min=15):
    keys = sorted(k for k in all_B if k > 0)
    even_n = [n for n in keys if n % 2 == 0 and n >= n_min]
    odd_n = [n for n in keys if n % 2 == 1 and n >= n_min]

    def logfit(ns):
        if len(ns) < 2:
            return mpmath.mpf(1), mpmath.mpf(7)
        n = len(ns)
        sum_x = sum(mpmath.log(k) for k in ns)
        sum_y = sum(mpmath.log(abs(all_B[k])) for k in ns)
        sum_xx = sum(mpmath.log(k) ** 2 for k in ns)
        sum_xy = sum(mpmath.log(k) * mpmath.log(abs(all_B[k])) for k in ns)
        denom = n * sum_xx - sum_x ** 2
        if abs(denom) < mpmath.mpf("1e-60"):
            return mpmath.mpf(1), mpmath.mpf(7)
        slope = (n * sum_xy - sum_x * sum_y) / denom
        intercept = (sum_y - slope * sum_x) / n
        return mpmath.exp(intercept), -slope

    C_e, p_e = logfit(even_n)
    C_o, p_o = logfit(odd_n)
    return C_e, p_e, C_o, p_o


def hurwitz_tail_mp(n_start, C_e, p_e, C_o, p_o):
    n_even = n_start if n_start % 2 == 0 else n_start + 1
    n_odd = n_start if n_start % 2 == 1 else n_start + 1
    t_e = C_e * 2 ** (-p_e) * mpmath.zeta(p_e, mpmath.mpf(n_even) / 2)
    t_o = C_o * 2 ** (-p_o) * mpmath.zeta(p_o, mpmath.mpf(n_odd) / 2)
    return t_e + t_o


def extract_c3_mp(all_B, V_mag, n_cutoff=None):
    if n_cutoff is None:
        n_cutoff = max(all_B.keys())

    cum_B = sum(all_B[k] for k in sorted(all_B) if k <= n_cutoff)

    n_fit = max(5, min(15, n_cutoff // 3))
    C_e, p_e, C_o, p_o = fit_power_law_mp(
        {k: all_B[k] for k in all_B if 0 < k <= n_cutoff}, n_min=n_fit)

    tail = hurwitz_tail_mp(n_cutoff + 1, C_e, p_e, C_o, p_o)

    dp = mpmath.mpf("0.05")
    tail_hi = hurwitz_tail_mp(n_cutoff + 1, C_e, p_e - dp, C_o, p_o - dp)
    tail_lo = hurwitz_tail_mp(n_cutoff + 1, C_e, p_e + dp, C_o, p_o + dp)
    tail_unc = max(abs(tail_hi - tail), abs(tail_lo - tail))

    total_B = cum_B + tail
    F2_S = total_B / V_mag / SCHWINGER
    delta = F2_S - 1

    # F2/S = total_B / (V_mag * alpha/(2pi)), so d(F2/S)/dalpha = -F2/S / alpha
    # dc3/dalpha = -F2/S / (alpha * x^3), propagated uncertainty = |dc3/dalpha| * ALPHA_UNC
    delta_alpha_unc = F2_S / ALPHA_CODATA * ALPHA_UNC / X3

    residual = delta - C1 * X - C2 * X2
    c3 = residual / X3

    c3_tail_unc = tail_unc / V_mag / SCHWINGER / X3
    c3_alpha_unc = delta_alpha_unc / X3
    c3_total_unc = mpmath.sqrt(c3_tail_unc ** 2 + c3_alpha_unc ** 2)

    return {
        "c3": c3,
        "c3_tail_unc": c3_tail_unc,
        "c3_alpha_unc": c3_alpha_unc,
        "c3_total_unc": c3_total_unc,
        "cum_B": cum_B,
        "tail": tail,
        "tail_unc": tail_unc,
        "delta": delta,
        "residual": residual,
        "C_e": C_e, "p_e": p_e,
        "C_o": C_o, "p_o": p_o,
        "n_cutoff": n_cutoff,
    }


# --------------------------------------------------------------------------
# Structural formula search (fast)
# --------------------------------------------------------------------------

def best_rational_approx(x, max_denom=100000):
    """Find best rational approximation to x with denominator <= max_denom."""
    f = Fraction(x).limit_denominator(max_denom)
    return f.numerator, f.denominator


def structural_search(c3_val, c3_unc, verbose=True):
    if verbose:
        print("\n" + "=" * 70)
        print("  STRUCTURAL FORMULA SEARCH")
        print("=" * 70)
        print(f"  Target c3 = {mpmath.nstr(c3_val, 8)}")
        print(f"  Uncertainty = {mpmath.nstr(c3_unc, 3)}")

    bd = float(BD)
    fd = float(FD)
    fb = float(FB)
    d = float(DELTA)
    c2f = float(C2)
    c3f = float(c3_val)
    c3_unc_f = float(c3_unc)

    terms = [
        ("1", 1.0), ("bd", bd), ("fd", fd), ("fb", fb), ("d", d),
        ("bd^2", bd**2), ("fd^2", fd**2), ("fb^2", fb**2),
        ("bd*fd", bd*fd), ("bd*fb", bd*fb), ("fd*fb", fd*fb),
        ("c2", c2f), ("d^2", d**2), ("bd*d", bd*d),
    ]

    best_hits = []

    # Strategy: for each pair of terms with coefficients a, b,
    # compute numer = a*v1 + b*v2, then N = numer / c3.
    # If N is close to integer, we have a hit.
    if verbose:
        print(f"\n  Searching 2-term combinations...")

    t0 = time.time()
    coeff_range = list(range(-50, 51))

    for i, (n1, v1) in enumerate(terms):
        for j, (n2, v2) in enumerate(terms):
            if j < i:
                continue
            for a in coeff_range:
                if a == 0:
                    continue
                for b in coeff_range:
                    if b == 0:
                        continue
                    if i == j and b < a:
                        continue
                    numer = a * v1 + b * v2
                    if abs(numer) < 1e-15:
                        continue
                    N_needed = numer / c3f
                    N_int = round(N_needed)
                    if N_int == 0 or abs(N_int) > 500000:
                        continue
                    candidate = numer / N_int
                    rel_err = abs(candidate - c3f) / abs(c3f)
                    if rel_err < 0.005:
                        desc = f"({a}*{n1} + {b}*{n2}) / {N_int}"
                        best_hits.append((desc, candidate, rel_err, N_int))

    if verbose:
        print(f"    2-term search done ({time.time()-t0:.1f}s), {len(best_hits)} hits within 0.5%")

    # 3-term with small coefficients
    if verbose:
        print(f"  Searching 3-term combinations...")

    t0 = time.time()
    small_range = list(range(-10, 11))
    core_terms = terms[:8]  # first 8 terms

    for i, (n1, v1) in enumerate(core_terms):
        for j, (n2, v2) in enumerate(core_terms):
            if j < i:
                continue
            for k, (n3, v3) in enumerate(core_terms):
                if k < j:
                    continue
                for a in small_range:
                    if a == 0:
                        continue
                    for b in small_range:
                        if b == 0:
                            continue
                        for c in small_range:
                            if c == 0:
                                continue
                            numer = a * v1 + b * v2 + c * v3
                            if abs(numer) < 1e-15:
                                continue
                            N_needed = numer / c3f
                            N_int = round(N_needed)
                            if N_int == 0 or abs(N_int) > 500000:
                                continue
                            candidate = numer / N_int
                            rel_err = abs(candidate - c3f) / abs(c3f)
                            if rel_err < 0.005:
                                desc = f"({a}*{n1} + {b}*{n2} + {c}*{n3}) / {N_int}"
                                best_hits.append((desc, candidate, rel_err, N_int))

    if verbose:
        print(f"    3-term search done ({time.time()-t0:.1f}s)")

    # Sort by relative error
    best_hits.sort(key=lambda x: x[2])

    # Deduplicate: keep only the best per N value
    seen_N = set()
    unique_hits = []
    for h in best_hits:
        if h[3] not in seen_N:
            seen_N.add(h[3])
            unique_hits.append(h)

    if verbose:
        print(f"\n  Top 30 unique hits (by rel_err):")
        for desc, val, rel, N in unique_hits[:30]:
            # Check if within uncertainty
            within_unc = abs(val - c3f) < 3 * c3_unc_f
            mark = " <<<" if within_unc else ""
            print(f"    {desc}")
            print(f"      value={val:.10e}  rel_err={rel:.6e}  N={N}{mark}")

    return unique_hits


# --------------------------------------------------------------------------
# PSLQ with extended bases
# --------------------------------------------------------------------------

def pslq_extended(c3_val, verbose=True):
    if verbose:
        print("\n" + "=" * 70)
        print("  EXTENDED PSLQ IDENTIFICATION")
        print("=" * 70)

    results = {}

    bd_mp = B * DELTA
    fd_mp = F * DELTA
    fb_mp = F / B
    d_mp = DELTA
    c2_mp = C2

    bases = {
        "T9_pure": (
            ["c3", "1", "pi2", "pi4"],
            [c3_val, mpmath.mpf(1), PI2, PI4]
        ),
        "T9_paper2_deg1": (
            ["c3", "1", "pi2", "pi4", "BD", "FD", "FB"],
            [c3_val, mpmath.mpf(1), PI2, PI4, bd_mp, fd_mp, fb_mp]
        ),
        "T9_paper2_deg2": (
            ["c3", "1", "pi2", "pi4", "BD", "FD", "FB",
             "BD^2", "FD^2", "FB^2", "BD*FD", "BD*FB", "FD*FB"],
            [c3_val, mpmath.mpf(1), PI2, PI4, bd_mp, fd_mp, fb_mp,
             bd_mp**2, fd_mp**2, fb_mp**2, bd_mp*fd_mp, bd_mp*fb_mp, fd_mp*fb_mp]
        ),
        "T9_with_c2": (
            ["c3", "1", "pi2", "pi4", "BD", "FD", "FB", "c2", "c2^2", "c1*c2"],
            [c3_val, mpmath.mpf(1), PI2, PI4, bd_mp, fd_mp, fb_mp, c2_mp, c2_mp**2, C1*c2_mp]
        ),
        "T9_with_delta_products": (
            ["c3", "1", "pi2", "pi4", "BD", "FD", "FB", "D", "D^2",
             "BD^2", "FD^2", "FB^2", "BD*D", "FD*D", "FB*D"],
            [c3_val, mpmath.mpf(1), PI2, PI4, bd_mp, fd_mp, fb_mp, d_mp, d_mp**2,
             bd_mp**2, fd_mp**2, fb_mp**2, bd_mp*d_mp, fd_mp*d_mp, fb_mp*d_mp]
        ),
        "T9_c2_cross": (
            ["c3", "c2*BD", "c2*FD", "c2*FB", "c2*D", "c2^2",
             "BD*FD", "BD*FB", "FD*FB", "D^2", "1"],
            [c3_val, c2_mp*bd_mp, c2_mp*fd_mp, c2_mp*fb_mp, c2_mp*d_mp, c2_mp**2,
             bd_mp*fd_mp, bd_mp*fb_mp, fd_mp*fb_mp, d_mp**2, mpmath.mpf(1)]
        ),
        "c3_times_LAM6": (
            ["c3*L6", "1", "pi2", "pi4", "BD", "FD", "FB",
             "BD^2", "FD^2", "FB^2"],
            [c3_val * LAM**6, mpmath.mpf(1), PI2, PI4, bd_mp, fd_mp, fb_mp,
             bd_mp**2, fd_mp**2, fb_mp**2]
        ),
        "residual_direct": (
            ["resid", "1", "pi2", "pi4", "BD", "FD", "FB"],
            [c3_val * X3, mpmath.mpf(1), PI2, PI4, bd_mp, fd_mp, fb_mp]
        ),
    }

    for basis_name, (labels, vals) in bases.items():
        vec = mpmath.matrix(vals)
        try:
            rel = mpmath.pslq(vec, maxcoeff=10000, tol=mpmath.mpf("1e-8"))
            if rel is not None:
                check = sum(float(rel[j]) * float(vals[j]) for j in range(len(rel)))
                if verbose:
                    print(f"\n  {basis_name}:")
                    print(f"    relation: {list(rel)}")
                    print(f"    labels:   {labels}")
                    print(f"    check:    {check:.3e}")
                    # Interpret: if rel[0] != 0, solve for c3
                    if rel[0] != 0:
                        formula_parts = []
                        for ii in range(1, len(rel)):
                            if rel[ii] != 0:
                                formula_parts.append(f"({int(-rel[ii])}/{int(rel[0])})*{labels[ii]}")
                        print(f"    c3 = {' + '.join(formula_parts)}")
                    else:
                        print(f"    NOTE: c3 coefficient is 0 (relation among basis only)")
                results[basis_name] = {
                    "relation": [int(r) for r in rel],
                    "labels": labels,
                    "check": float(check),
                    "c3_coeff": int(rel[0]),
                }
            else:
                if verbose:
                    print(f"  {basis_name}: no relation found")
                results[basis_name] = None
        except Exception as e:
            if verbose:
                print(f"  {basis_name}: error: {e}")
            results[basis_name] = None

    return results


# --------------------------------------------------------------------------
# T9 decomposition probe (using continued fractions)
# --------------------------------------------------------------------------

def t9_decomposition(c3_val, c3_unc, verbose=True):
    """
    Probe c3 = r0 + r2*pi^2 + r4*pi^4 decomposition.
    Use continued fractions for efficient rational approximation.
    """
    if verbose:
        print("\n" + "=" * 70)
        print("  T9 DECOMPOSITION PROBE")
        print("=" * 70)

    c3f = float(c3_val)
    pi2f = float(PI2)
    pi4f = float(PI4)
    c3_unc_f = float(c3_unc)

    hits = []

    # Strategy: for small rationals r0 = p/q, check if (c3 - r0)/pi^2 is rational
    # and if (c3 - r0)/pi^4 is rational
    # Since |c3| ~ 6e-7, need |r0| ~ 6e-7 or r0 ~ 0 with r2*pi^2 ~ 6e-7

    # Case A: r0 = 0, c3 = r2*pi^2 + r4*pi^4
    r2_only = c3f / pi2f  # ~ -6.03e-8
    frac_r2 = Fraction(r2_only).limit_denominator(10000000)
    residual_A = abs(c3f - float(frac_r2) * pi2f)
    if verbose:
        print(f"\n  Case A: c3 = r2*pi^2")
        print(f"    r2 = c3/pi^2 = {r2_only:.12e}")
        print(f"    best rational: {frac_r2} = {float(frac_r2):.12e}")
        print(f"    residual: {residual_A:.3e} (unc: {c3_unc_f:.3e})")

    r4_only = c3f / pi4f  # ~ -6.11e-9
    frac_r4 = Fraction(r4_only).limit_denominator(10000000)
    residual_B = abs(c3f - float(frac_r4) * pi4f)
    if verbose:
        print(f"\n  Case B: c3 = r4*pi^4")
        print(f"    r4 = c3/pi^4 = {r4_only:.12e}")
        print(f"    best rational: {frac_r4} = {float(frac_r4):.12e}")
        print(f"    residual: {residual_B:.3e} (unc: {c3_unc_f:.3e})")

    # Case C: c3 = r0 + r2*pi^2, scan r0 = n/d for small n, large d
    if verbose:
        print(f"\n  Case C: c3 = n/d + r2*pi^2")
        print(f"    Scanning n in [-5,5], d in [1, 10^7]...")

    t0 = time.time()
    for n0 in range(-5, 6):
        if n0 == 0:
            continue
        # c3 = n0/d0 + r2*pi^2
        # For this to work: n0/d0 must be close to c3, i.e. d0 ~ n0/c3
        # Or: r2*pi^2 must be close to c3, meaning n0/d0 ~ 0 (large d0)
        # Let's scan d0 such that |n0/d0| < 10*|c3|
        d_min = max(1, int(abs(n0) / (10 * abs(c3f))))
        d_max = min(10000000, int(abs(n0) / max(abs(c3f) / 100, 1e-20)))
        for d0 in range(d_min, min(d_max + 1, 10000001)):
            r0 = n0 / d0
            remainder = c3f - r0
            r2_cand = remainder / pi2f
            frac = Fraction(r2_cand).limit_denominator(100000)
            r2_approx = float(frac)
            residual = abs(c3f - r0 - r2_approx * pi2f)
            if residual < 3 * c3_unc_f:
                hits.append({
                    "type": "r0+r2*pi2",
                    "n0": n0, "d0": d0,
                    "r2_num": frac.numerator, "r2_den": frac.denominator,
                    "residual": residual,
                    "formula": f"{n0}/{d0} + ({frac.numerator}/{frac.denominator})*pi^2",
                })

    if verbose:
        print(f"    Case C done ({time.time()-t0:.1f}s), {len(hits)} hits within 3*unc")

    # Case D: c3 = r2*pi^2 + r4*pi^4 (no rational part)
    if verbose:
        print(f"\n  Case D: c3 = r2*pi^2 + r4*pi^4")
        print(f"    Scanning r2 = n2/d2 with small denominators...")

    t0 = time.time()
    for d2 in range(1, 100001):
        for n2 in range(-100, 101):
            r2 = n2 / d2
            remainder = c3f - r2 * pi2f
            r4_cand = remainder / pi4f
            frac4 = Fraction(r4_cand).limit_denominator(100000)
            r4_approx = float(frac4)
            residual = abs(c3f - r2 * pi2f - r4_approx * pi4f)
            if residual < 3 * c3_unc_f:
                hits.append({
                    "type": "r2*pi2+r4*pi4",
                    "r2_num": n2, "r2_den": d2,
                    "r4_num": frac4.numerator, "r4_den": frac4.denominator,
                    "residual": residual,
                    "formula": f"({n2}/{d2})*pi^2 + ({frac4.numerator}/{frac4.denominator})*pi^4",
                })

    if verbose:
        print(f"    Case D done ({time.time()-t0:.1f}s), {len([h for h in hits if h['type']=='r2*pi2+r4*pi4'])} hits")

    hits.sort(key=lambda x: x["residual"])

    if verbose:
        print(f"\n  Top 20 T9 decomposition hits:")
        for h in hits[:20]:
            print(f"    {h['formula']}  resid={h['residual']:.3e}")

    return hits


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("  ALGEBRAIC IDENTIFICATION OF c3")
    print("=" * 70)

    all_B, V_mag = load_data()
    n_max = max(all_B.keys())
    print(f"  Loaded B(n_int) for n_int=0..{n_max}")
    print(f"  V_mag = {mpmath.nstr(V_mag, 15)}")

    # Step 1: High-precision extraction
    print("\n" + "=" * 70)
    print("  STEP 1: HIGH-PRECISION EXTRACTION")
    print("=" * 70)

    result = extract_c3_mp(all_B, V_mag)
    c3 = result["c3"]

    print(f"  n_cutoff     = {result['n_cutoff']}")
    print(f"  cum_B        = {mpmath.nstr(result['cum_B'], 20)}")
    print(f"  tail         = {mpmath.nstr(result['tail'], 10)}")
    print(f"  tail_unc     = {mpmath.nstr(result['tail_unc'], 4)}")
    print(f"  delta        = {mpmath.nstr(result['delta'], 20)}")
    print(f"  C1*x         = {mpmath.nstr(C1 * X, 20)}")
    print(f"  C2*x^2       = {mpmath.nstr(C2 * X2, 20)}")
    print(f"  residual     = {mpmath.nstr(result['residual'], 10)}")
    print(f"  c3           = {mpmath.nstr(c3, 10)}")
    print(f"  c3 tail unc  = {mpmath.nstr(result['c3_tail_unc'], 4)}")
    print(f"  c3 alpha unc = {mpmath.nstr(result['c3_alpha_unc'], 4)}")
    print(f"  c3 total unc = {mpmath.nstr(result['c3_total_unc'], 4)}")
    sig = abs(c3) / result['c3_total_unc']
    print(f"  significance = {mpmath.nstr(sig, 4)} sigma")

    # Convergence study
    print("\n  Convergence study:")
    for cutoff in [15, 20, 25, 30, 35, 40, 45, 50]:
        if cutoff > n_max:
            break
        r = extract_c3_mp(all_B, V_mag, n_cutoff=cutoff)
        print(f"    n={cutoff:3d}: c3 = {mpmath.nstr(r['c3'], 8)}  "
              f"unc = {mpmath.nstr(r['c3_total_unc'], 3)}  "
              f"p_even = {mpmath.nstr(r['p_e'], 5)}  "
              f"p_odd = {mpmath.nstr(r['p_o'], 5)}")

    # Step 2: Structural formula search
    best_hits = structural_search(c3, result["c3_total_unc"])

    # Step 3: Extended PSLQ
    pslq_results = pslq_extended(c3)

    # Step 4: T9 decomposition
    t9_hits = t9_decomposition(c3, result["c3_total_unc"])

    # Save results
    output = {
        "c3": float(c3),
        "c3_tail_unc": float(result["c3_tail_unc"]),
        "c3_alpha_unc": float(result["c3_alpha_unc"]),
        "c3_total_unc": float(result["c3_total_unc"]),
        "c3_significance_sigma": float(sig),
        "delta": float(result["delta"]),
        "residual": float(result["residual"]),
        "n_cutoff": result["n_cutoff"],
        "power_law": {
            "C_even": float(result["C_e"]),
            "p_even": float(result["p_e"]),
            "C_odd": float(result["C_o"]),
            "p_odd": float(result["p_o"]),
        },
        "structural_hits_top20": [
            {"desc": d, "value": v, "rel_err": r, "N": n}
            for d, v, r, n in (best_hits[:20] if best_hits else [])
        ],
        "pslq_results": {k: (v if v is not None else "null") for k, v in pslq_results.items()},
        "t9_hits_top20": [h for h in (t9_hits[:20] if t9_hits else [])],
    }

    outfile = DATA_DIR / "g2_c3_algebraic_id.json"
    with open(outfile, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to {outfile}")

    # Summary
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  c3 = {mpmath.nstr(c3, 8)} +/- {mpmath.nstr(result['c3_total_unc'], 3)}")
    dom = 'alpha' if result['c3_alpha_unc'] > result['c3_tail_unc'] else 'tail'
    print(f"  Dominant uncertainty: {dom}")
    if best_hits:
        print(f"  Best structural hit: {best_hits[0][0]} (rel_err={best_hits[0][2]:.3e})")
    if t9_hits:
        print(f"  Best T9 decomposition: {t9_hits[0]['formula']} (residual={t9_hits[0]['residual']:.3e})")


if __name__ == "__main__":
    main()
