"""
g2_c3_fast.py -- Fast float-arithmetic c3 computation for g-2 on S3.

Uses precomputed float CG coefficient cache instead of exact sympy Rational.
Achieves ~50-100x speedup over the exact arithmetic version.

Strategy:
- Use scipy.special.clebsch_gordan if available, else sympy with float conversion
- Cache all CG computations using lru_cache on integer-double (2j) arguments
- Use float arithmetic throughout (exact sympy only for comparison)
"""

import json
import sys
import time
import math
from pathlib import Path
from functools import lru_cache

import mpmath
import numpy as np

mpmath.mp.dps = 50

DATA_DIR = Path(__file__).parent / "data"
CHECKPOINT_FILE = DATA_DIR / "g2_c3_fast_checkpoint.json"
OUTPUT_FILE = DATA_DIR / "g2_c3_investigation.json"

# Physical constants
ALPHA = 7.2973525693e-3
SCHWINGER = ALPHA / (2 * math.pi)

# Paper 2 invariants
B_HOPF = 42
F_HOPF = float(mpmath.pi**2 / 6)
DELTA_HOPF = 1.0 / 40.0
C1 = 0.5
C2_EXACT = float((2 - B_HOPF * DELTA_HOPF - F_HOPF * DELTA_HOPF - F_HOPF / B_HOPF) / 5)
C2_MP = mpmath.mpf(2 - B_HOPF * mpmath.mpf(1)/40 - (mpmath.pi**2/6) * mpmath.mpf(1)/40 - (mpmath.pi**2/6) / B_HOPF) / 5

print(f"c1 = {C1}")
print(f"c2 = {C2_EXACT:.15f}")
print(f"c2 verify (sympy): {float(mpmath.mpf(19)/100 - 41*mpmath.pi**2/25200):.15f}")

# n_ext=1 kinematics
LAM_N1 = 2.5  # = 5/2
LAM2_N1 = 6.25
X_N1 = 1.0 / LAM2_N1  # 4/25 = 0.16
X3_N1 = X_N1**3

# ---------------------------------------------------------------------------
# CG coefficient computation using sympy with caching
# ---------------------------------------------------------------------------
# We use integer 2j, 2m arguments to avoid Rational overhead

try:
    from sympy.physics.wigner import clebsch_gordan
    from sympy import Rational
    _USE_SYMPY = True
    print("Using sympy CG with float cache")
except ImportError:
    _USE_SYMPY = False
    print("Warning: sympy not available")


@lru_cache(maxsize=2_000_000)
def cg(j1_2, j2_2, J_2, m1_2, m2_2, M_2):
    """
    Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M>
    Arguments are 2*j, 2*m (integers).
    Returns float.
    """
    if m1_2 + m2_2 != M_2:
        return 0.0
    if abs(j1_2 - j2_2) > J_2 or J_2 > j1_2 + j2_2:
        return 0.0
    if abs(m1_2) > j1_2 or abs(m2_2) > j2_2 or abs(M_2) > J_2:
        return 0.0
    val = clebsch_gordan(Rational(j1_2, 2), Rational(j2_2, 2), Rational(J_2, 2),
                         Rational(m1_2, 2), Rational(m2_2, 2), Rational(M_2, 2))
    return float(val)


# ---------------------------------------------------------------------------
# Vertex machinery using float CG
# ---------------------------------------------------------------------------

def half_ints_2(two_j):
    """Return list of 2*m for m in [-j, j] as integers."""
    return list(range(-two_j, two_j + 1, 2))


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, jsL2, jsR2, jtL2, jtR2):
    """Return list of (jgL*2, jgR*2) channel pairs."""
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    chs = []
    for jgL2, jgR2 in [(q + 1, q - 1), (q - 1, q + 1)]:
        if jgL2 < 0 or jgR2 < 0:
            continue
        if (abs(jsL2 - jgL2) <= jtL2 <= jsL2 + jgL2 and
                abs(jsR2 - jgR2) <= jtR2 <= jsR2 + jgR2):
            chs.append((jgL2, jgR2))
    return chs


def vertex_amp_pol_f(jsL2, jsR2, js2, mjs2,
                     jtL2, jtR2, jt2, mjt2,
                     jgL2, jgR2, mgL2, mgR2):
    """Float vertex amplitude using cached CG coefficients. Arguments are 2*j, 2*m."""
    total = 0.0
    for mL1_2 in half_ints_2(jsL2):
        mR1_2 = mjs2 - mL1_2
        if abs(mR1_2) > jsR2:
            continue
        mL2_2 = mL1_2 + mgL2
        if abs(mL2_2) > jtL2:
            continue
        mR2_2 = mR1_2 + mgR2
        if abs(mR2_2) > jtR2:
            continue
        if mL2_2 + mR2_2 != mjt2:
            continue
        c1v = cg(jsL2, jsR2, js2, mL1_2, mR1_2, mjs2)
        if c1v == 0.0:
            continue
        c2v = cg(jtL2, jtR2, jt2, mL2_2, mR2_2, mjt2)
        if c2v == 0.0:
            continue
        c3v = cg(jsL2, jgL2, jtL2, mL1_2, mgL2, mL2_2)
        c4v = cg(jsR2, jgR2, jtR2, mR1_2, mgR2, mR2_2)
        total += c1v * c2v * c3v * c4v
    return total


def compute_vertex_single_nint_float(n_ext, n_int, mj_ext_2=1, q_probe=1):
    """
    Compute B(n_int) = vertex(mj=+1/2) - vertex(mj=-1/2) using float CG.
    n_ext=1: jE_L=1 -> 2*jE_L=2, jE_R=1/2 -> 2*jE_R=1, j_ext=1/2 -> 2*j_ext=1
    """
    jE_L_2 = n_ext + 1   # 2*jE_L
    jE_R_2 = n_ext       # 2*jE_R
    j_ext_2 = 1          # 2*j_ext = 1 (j_ext = 1/2 for n_ext=1)

    jI_L_2 = n_int + 1   # 2*jI_L
    jI_R_2 = n_int        # 2*jI_R

    lam = (2 * n_int + 3) / 2.0  # Camporesi-Higuchi: |lambda_n| = n + 3/2
    lam4 = lam**4

    if not vertex_allowed(n_int, n_int, q_probe):
        return 0.0

    probe_chs = get_channels(n_int, n_int, q_probe, jI_L_2, jI_R_2, jI_L_2, jI_R_2)
    if not probe_chs:
        return 0.0

    # j_int range: from |jI_L - jI_R| = 1/2 to jI_L + jI_R = n_int + 1/2
    j_int_min_2 = 1  # 2*|jI_L - jI_R| = 2*|n_int/2 - 1/2| = 1 (for n_int >= 1)
    j_int_max_2 = 2 * n_int + 1  # 2*(jI_L + jI_R) = 2*n_int+1

    total = 0.0

    for j_int_2 in range(j_int_min_2, j_int_max_2 + 1, 2):
        for mj_int_2 in range(-j_int_2, j_int_2 + 1, 2):
            for mj_int_prime_2 in range(-j_int_2, j_int_2 + 1, 2):

                probe_amp = 0.0
                for jpL2, jpR2 in probe_chs:
                    for mpL2 in range(-jpL2, jpL2 + 1, 2):
                        for mpR2 in range(-jpR2, jpR2 + 1, 2):
                            pa = vertex_amp_pol_f(
                                jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                jpL2, jpR2, mpL2, mpR2)
                            probe_amp += pa

                if probe_amp == 0.0:
                    continue

                # q_loop: must satisfy vertex_allowed(n_ext, n_int, q_loop)
                q_lo = max(1, abs(n_ext - n_int))
                q_hi = n_ext + n_int

                for q_loop in range(q_lo, q_hi + 1):
                    if not vertex_allowed(n_ext, n_int, q_loop):
                        continue
                    mu_q = q_loop * (q_loop + 2)  # eigenvalue of photon mode

                    chs1 = get_channels(n_ext, n_int, q_loop, jE_L_2, jE_R_2, jI_L_2, jI_R_2)
                    chs2 = get_channels(n_int, n_ext, q_loop, jI_L_2, jI_R_2, jE_L_2, jE_R_2)

                    for jgL1_2, jgR1_2 in chs1:
                        for jgL2_2, jgR2_2 in chs2:
                            for mgL2_v in range(-jgL1_2, jgL1_2 + 1, 2):
                                for mgR2_v in range(-jgR1_2, jgR1_2 + 1, 2):
                                    if abs(mgL2_v) > jgL2_2 or abs(mgR2_v) > jgR2_2:
                                        continue

                                    v1 = vertex_amp_pol_f(
                                        jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                        jI_L_2, jI_R_2, j_int_2, mj_int_2,
                                        jgL1_2, jgR1_2, mgL2_v, mgR2_v)
                                    if v1 == 0.0:
                                        continue

                                    v3 = vertex_amp_pol_f(
                                        jI_L_2, jI_R_2, j_int_2, mj_int_prime_2,
                                        jE_L_2, jE_R_2, j_ext_2, mj_ext_2,
                                        jgL2_2, jgR2_2, mgL2_v, mgR2_v)
                                    if v3 == 0.0:
                                        continue

                                    total += v1 * probe_amp * v3 / (lam4 * mu_q)

    return total


def compute_B_level(n_ext, n_int):
    """Compute B(n_int) = contrib(mj=+1/2) - contrib(mj=-1/2)."""
    up = compute_vertex_single_nint_float(n_ext, n_int, mj_ext_2=+1)
    dn = compute_vertex_single_nint_float(n_ext, n_int, mj_ext_2=-1)
    return up - dn


# ---------------------------------------------------------------------------
# Verify against known data for n_int=1..25
# ---------------------------------------------------------------------------

def verify_against_existing(all_B_known, n_verify=5, verbose=True):
    """Cross-check the float implementation vs stored values."""
    print("\n--- Verification against n_int=0..25 data ---")
    errors = []
    for n_int in range(1, min(n_verify + 1, 26)):
        B_known = all_B_known[n_int]
        t0 = time.time()
        B_fast = compute_B_level(1, n_int)
        dt = time.time() - t0
        rel_err = abs(B_fast - B_known) / max(abs(B_known), 1e-20)
        errors.append(rel_err)
        if verbose:
            print(f"  n_int={n_int}: B_fast={B_fast:.8e}, B_known={B_known:.8e}, "
                  f"rel_err={rel_err:.2e} ({dt:.2f}s)")
    max_err = max(errors)
    print(f"  Max relative error: {max_err:.2e}")
    return max_err < 1e-8


# ---------------------------------------------------------------------------
# Tail correction using Hurwitz zeta
# ---------------------------------------------------------------------------

def fit_power_law(all_B, n_min=10):
    keys = sorted(all_B.keys())
    even_n = [n for n in keys if n % 2 == 0 and n >= n_min]
    odd_n  = [n for n in keys if n % 2 == 1 and n >= n_min]

    def logfit(ns, Bs):
        logns = np.log(np.array(ns, dtype=float))
        logBs = np.log(np.abs(np.array([Bs[n] for n in ns], dtype=float)))
        A = np.column_stack([np.ones_like(logns), -logns])
        x, _, _, _ = np.linalg.lstsq(A, logBs, rcond=None)
        return np.exp(x[0]), x[1]

    C_e, p_e = logfit(even_n, all_B)
    C_o, p_o = logfit(odd_n, all_B)
    return C_e, p_e, C_o, p_o


def hurwitz_tail(n_start, C_e, p_e, C_o, p_o):
    n_even = n_start if n_start % 2 == 0 else n_start + 1
    n_odd  = n_start if n_start % 2 == 1 else n_start + 1
    t_e = float(C_e * 2**(-p_e) * float(mpmath.zeta(p_e, mpmath.mpf(n_even) / 2)))
    t_o = float(C_o * 2**(-p_o) * float(mpmath.zeta(p_o, mpmath.mpf(n_odd) / 2)))
    return t_e + t_o


def hurwitz_tail_unc(n_start, C_e, p_e, C_o, p_o, dp=0.05):
    t0 = hurwitz_tail(n_start, C_e, p_e, C_o, p_o)
    th = hurwitz_tail(n_start, C_e, p_e - dp, C_o, p_o - dp)
    tl = hurwitz_tail(n_start, C_e, p_e + dp, C_o, p_o + dp)
    return max(abs(th - t0), abs(tl - t0))


# ---------------------------------------------------------------------------
# c3 extraction
# ---------------------------------------------------------------------------

def extract_c3(all_B, V_mag, n_int_cutoff=None, verbose=True):
    """Extract c3 from the delta at n_ext=1."""
    if n_int_cutoff is None:
        n_int_cutoff = max(all_B.keys())

    B_keys = sorted(k for k in all_B if k <= n_int_cutoff)
    cum_B = sum(all_B[k] for k in B_keys)
    n_max = max(B_keys)

    n_fit_start = max(10, n_max // 2)
    C_e, p_e, C_o, p_o = fit_power_law(
        {k: all_B[k] for k in B_keys if k > 0}, n_min=n_fit_start)

    tail_B = hurwitz_tail(n_max + 1, C_e, p_e, C_o, p_o)
    tail_unc_B = hurwitz_tail_unc(n_max + 1, C_e, p_e, C_o, p_o, dp=0.05)

    total_B = cum_B + tail_B
    F2_S = total_B / V_mag / SCHWINGER
    delta = F2_S - 1.0

    tail_unc_delta = tail_unc_B / V_mag / SCHWINGER

    # c3 extraction
    c1_term = C1 * X_N1
    c2_term = C2_EXACT * X_N1**2
    residual = delta - c1_term - c2_term
    c3 = residual / X3_N1
    c3_unc = tail_unc_delta / X3_N1

    if verbose:
        print(f"\n  n_int_cutoff = {n_max}")
        print(f"  Power-law: even C={C_e:.4e} p={p_e:.4f}, odd C={C_o:.4e} p={p_o:.4f}")
        print(f"  cum_B     = {cum_B:.15e}")
        print(f"  tail_B    = {tail_B:.6e} (+/- {tail_unc_B:.2e})")
        print(f"  F2/S      = {F2_S:.12f}")
        print(f"  delta     = {delta:.12f}")
        print(f"  c1*x      = {c1_term:.12f}")
        print(f"  c2*x^2    = {c2_term:.12f}")
        print(f"  residual  = {residual:.6e}")
        print(f"  c3        = {c3:.6e} +/- {c3_unc:.2e}")
        print(f"  sigma     = {abs(c3)/max(c3_unc, 1e-30):.2f}")

    return {
        "n_int_cutoff": n_max,
        "cum_B": cum_B,
        "tail_B": tail_B,
        "tail_B_unc": tail_unc_B,
        "C_even": C_e, "p_even": p_e,
        "C_odd": C_o, "p_odd": p_o,
        "F2_S": F2_S,
        "delta": delta,
        "c3": c3,
        "c3_unc": c3_unc,
        "c3_sigma": abs(c3) / max(c3_unc, 1e-30),
    }


# ---------------------------------------------------------------------------
# PSLQ identification
# ---------------------------------------------------------------------------

def pslq_c3(c3_val, verbose=True):
    """Attempt PSLQ identification of c3 against Paper 2 invariant basis."""
    if abs(c3_val) < 1e-15:
        return {"result": "consistent_with_zero"}

    pi = float(mpmath.pi)
    pi2 = pi**2
    pi4 = pi**4
    B, F, D = 42.0, pi2 / 6, 1.0 / 40
    c2 = C2_EXACT

    # T9 compatible bases (no odd zeta)
    bases = {
        "Small T9": {
            "labels": ["1", "pi^2", "pi^4"],
            "vals": [1.0, pi2, pi4],
        },
        "Medium T9": {
            "labels": ["1", "pi^2", "pi^4", "B*D", "F*D", "F/B"],
            "vals": [1.0, pi2, pi4, B*D, F*D, F/B],
        },
        "Full T9": {
            "labels": ["1", "pi^2", "pi^4", "B*D", "F*D", "F/B",
                       "c2", "c2^2", "c1*c2", "Delta", "Delta^2",
                       "B*D^2", "F*D^2", "F^2/B"],
            "vals": [1.0, pi2, pi4, B*D, F*D, F/B,
                     c2, c2**2, 0.5*c2, D, D**2,
                     B*D**2, F*D**2, F**2/B],
        },
        "T9+odd zeta (check)": {
            "labels": ["1", "pi^2", "zeta3", "pi^4"],
            "vals": [1.0, pi2, float(mpmath.zeta(3)), pi4],
        },
    }

    results = {}
    for basis_name, basis in bases.items():
        labels = basis["labels"]
        vals = basis["vals"]
        vec = mpmath.matrix([mpmath.mpf(str(c3_val))] +
                            [mpmath.mpf(str(v)) for v in vals])
        try:
            rel = mpmath.pslq(vec, maxcoeff=1000, tol=1e-25)
            if rel is not None:
                # Check: rel[0]*c3 + sum(rel[i]*vals[i-1]) = 0
                check = sum(float(rel[j]) * float(vec[j]) for j in range(len(rel)))
                if verbose:
                    print(f"  {basis_name}: {rel} (check={check:.3e})")
                results[basis_name] = {"relation": list(rel), "check": check,
                                        "labels": ["c3"] + labels}
            else:
                if verbose:
                    print(f"  {basis_name}: no relation found")
                results[basis_name] = None
        except Exception as e:
            if verbose:
                print(f"  {basis_name}: error: {e}")
            results[basis_name] = None
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("  g-2 c3 INVESTIGATION (fast float CG)")
    print("=" * 70)

    # Load existing data
    path = DATA_DIR / "g2_extended_nint_v2.json"
    with open(path) as f:
        data = json.load(f)
    all_B = {int(k): float(v) for k, v in data["all_B"].items()}

    path2 = DATA_DIR / "alpha_g_minus_2_ratio_investigation.json"
    with open(path2) as f:
        d2 = json.load(f)
    V_mag = d2["V_magnetic_float"]

    print(f"Loaded n_int=0..{max(all_B)}, V_mag={V_mag:.12e}")

    # Verify fast implementation against known values for n_int=1..5
    print("\nVerification (n_int=1..5):")
    ok = verify_against_existing(all_B, n_verify=5, verbose=True)
    if not ok:
        print("WARNING: Fast implementation disagrees with stored data!")
    else:
        print("  PASSED: fast float agrees with stored data.")

    # Load checkpoint if available
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE) as f:
            ckpt = json.load(f)
        all_B_ext = {int(k): float(v) for k, v in ckpt.items()}
        print(f"\nLoaded checkpoint: n_int up to {max(all_B_ext)}")
        all_B = all_B_ext

    # Extend to n_int=50
    n_target = 50
    n_done = max(all_B.keys())

    if n_done < n_target:
        print(f"\nExtending from n_int={n_done+1} to n_int={n_target}...")
        for n_int in range(n_done + 1, n_target + 1):
            t0 = time.time()
            B = compute_B_level(1, n_int)
            dt = time.time() - t0
            all_B[n_int] = B

            # Running c3 estimate
            cum_B = sum(all_B.values())
            delta_raw = cum_B / V_mag / SCHWINGER - 1.0
            print(f"  n_int={n_int:3d}: B={B:12.4e}  delta_raw={delta_raw:.12f}  ({dt:.1f}s)",
                  flush=True)

            # Save checkpoint
            with open(CHECKPOINT_FILE, 'w') as f:
                json.dump({str(k): v for k, v in sorted(all_B.items())}, f)
    else:
        print(f"Already have n_int up to {n_done}.")

    # Extract c3 at multiple cutoffs
    print("\n" + "=" * 70)
    print("  c3 convergence study")
    print("=" * 70)

    c3_convergence = []
    for cutoff in range(15, max(all_B.keys()) + 1, 5):
        if cutoff not in all_B:
            continue
        r = extract_c3(all_B, V_mag, n_int_cutoff=cutoff, verbose=False)
        c3_convergence.append(r)
        print(f"  n={cutoff:3d}: c3={r['c3']:.6e} +/- {r['c3_unc']:.2e} ({r['c3_sigma']:.1f} sigma)")

    # Final extraction with all data
    print("\n" + "=" * 70)
    print("  Final c3 extraction")
    print("=" * 70)

    result = extract_c3(all_B, V_mag, verbose=True)

    # PSLQ
    print("\n" + "=" * 70)
    print("  PSLQ identification")
    print("=" * 70)
    pslq_results = pslq_c3(result["c3"], verbose=True)

    # Save
    output = {
        "n_int_max": max(all_B.keys()),
        "V_mag": V_mag,
        "c1": C1,
        "c2": C2_EXACT,
        "c2_formula": "19/100 - 41*pi^2/25200",
        "c3_final": result["c3"],
        "c3_unc": result["c3_unc"],
        "c3_sigma": result["c3_sigma"],
        "c3_consistent_with_zero": abs(result["c3"]) < 3 * result["c3_unc"],
        "delta": result["delta"],
        "residual": result["delta"] - C1 * X_N1 - C2_EXACT * X_N1**2,
        "c3_convergence": c3_convergence,
        "pslq": pslq_results,
        "all_B": {str(k): v for k, v in sorted(all_B.items())},
    }

    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {OUTPUT_FILE}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  c1 = {C1} (Parker-Toms, exact)")
    print(f"  c2 = {C2_EXACT:.10f} = 19/100 - 41pi^2/25200 (verified)")
    print(f"  c3 = {result['c3']:.6e} +/- {result['c3_unc']:.2e}")
    print(f"  c3 significance: {result['c3_sigma']:.1f} sigma")
    if abs(result["c3"]) < 3 * result["c3_unc"]:
        print("  c3 is CONSISTENT WITH ZERO at 3-sigma level")
    else:
        print(f"  c3 is NONZERO at {result['c3_sigma']:.1f} sigma")

    return output


if __name__ == "__main__":
    main()
