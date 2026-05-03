"""
g2_c3_investigation.py -- High-precision investigation of c_3 in the S^3 g-2 curvature expansion.

The anomalous magnetic moment F_2 on S^3 has a curvature expansion:

    F_2/Schwinger = 1 + c_1/lambda^2 + c_2/lambda^4 + c_3/lambda^6 + ...

where lambda = n_ext + 3/2 is the Dirac eigenvalue on S^3 (Camporesi-Higuchi).

Known exactly:
    c_1 = 1/2          (Parker-Toms, Riemann curvature coupling)
    c_2 = 19/100 - 41pi^2/25200 = (2 - BDelta - FDelta - F/B)/5 ~= 0.17394231030

    where B=42, F=pi^2/6, Delta=1/40 are the Paper 2 invariants.
    c_2 verified to 8 digits from n_ext=1 spectral sum at n_int=0..25.

This script:
1. Loads the existing n_int=0..25 data (all_B table)
2. Extends to n_int=50 by computing new B levels (dominant computation)
3. Uses Euler-Maclaurin tail correction to handle n_int>50
4. Extracts c_3 = (delta - c_1/lambda^2 - c_2/lambda^4) * lambda^6 with error bounds
5. Runs PSLQ against Paper 2 invariant basis to find a closed form

Strategy for c_3 precision:
- At n_ext=1, lambda=5/2, x=1/lambda^2=4/25=0.16
- delta ~= 0.08445 so c_3 correction is c_3*x^3 = c_3*(0.16)^3 = c_3*0.004096
- To see c_3 at 10^-^7 level we need delta at ~4x10^-^1^0 level
- Current tail at n_int=25 is ~6.7e-9 (too large by factor ~17)
- Running to n_int=50: tail drops to ~6.7e-9 * (25/50)^6.85 ~= 5e-11
- This gives c_3 at ~1e-8 level = 1 digit on the 10^-^7 scale

IMPORTANT: The n_int=25 data shows B(n) follows two interleaved power laws:
  Even n: B ~ 2.09e-3 * n^(-6.74)
  Odd n:  B ~ 8.49e-3 * n^(-6.95)
These exponents give a well-converged tail after n=50.
"""

import json
import sys
import time
from pathlib import Path
from fractions import Fraction

import mpmath
import numpy as np
from sympy import Rational, S
from sympy.physics.wigner import clebsch_gordan

mpmath.mp.dps = 50

# Paths
DATA_DIR = Path(__file__).parent / "data"
OUTPUT_FILE = DATA_DIR / "g2_c3_investigation.json"

# Physical constants
ALPHA = 7.2973525693e-3
SCHWINGER = ALPHA / (2 * np.pi)

# Paper 2 invariants (exact)
B_HOPF = 42
F_HOPF = mpmath.pi**2 / 6
DELTA_HOPF = mpmath.mpf(1) / 40

# Known curvature coefficients
C1 = mpmath.mpf(1) / 2
# c_2 = (2 - B*Delta - F*Delta - F/B) / 5 = 19/100 - 41pi^2/25200
C2_EXACT = (2 - B_HOPF * DELTA_HOPF - F_HOPF * DELTA_HOPF - F_HOPF / B_HOPF) / 5
# Verify: should be 19/100 - 41*pi^2/25200
C2_RATIONAL = mpmath.mpf(19) / 100 - 41 * mpmath.pi**2 / 25200

# n_ext=1: lambda = 5/2
LAM_N1 = mpmath.mpf(5) / 2
LAM2_N1 = LAM_N1**2  # 25/4
LAM4_N1 = LAM_N1**4  # 625/16
LAM6_N1 = LAM_N1**6  # 15625/64


# ---------------------------------------------------------------------------
# CG vertex machinery (copied from converged scripts, same logic)
# ---------------------------------------------------------------------------

def half_ints(j):
    result = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


def vertex_allowed(n1, n2, q):
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_channels(n_src, n_tgt, q, j_sL, j_sR, j_tL, j_tR):
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    chs = []
    for jgL, jgR in [(Rational(q+1, 2), Rational(q-1, 2)),
                     (Rational(q-1, 2), Rational(q+1, 2))]:
        if jgL < 0 or jgR < 0:
            continue
        if (abs(j_sL - jgL) <= j_tL <= j_sL + jgL and
                abs(j_sR - jgR) <= j_tR <= j_sR + jgR):
            chs.append((jgL, jgR))
    return chs


def vertex_amp_pol(j_sL, j_sR, j_s, mj_s,
                   j_tL, j_tR, j_t, mj_t,
                   jgL, jgR, mgL, mgR):
    total = S.Zero
    for mL1 in half_ints(j_sL):
        mR1 = mj_s - mL1
        if abs(mR1) > j_sR:
            continue
        mL2 = mL1 + mgL
        if abs(mL2) > j_tL:
            continue
        mR2 = mR1 + mgR
        if abs(mR2) > j_tR:
            continue
        if mL2 + mR2 != mj_t:
            continue
        c1 = clebsch_gordan(j_sL, j_sR, j_s, mL1, mR1, mj_s)
        if c1 == 0:
            continue
        c2 = clebsch_gordan(j_tL, j_tR, j_t, mL2, mR2, mj_t)
        if c2 == 0:
            continue
        c3 = clebsch_gordan(j_sL, jgL, j_tL, mL1, mgL, mL2)
        c4 = clebsch_gordan(j_sR, jgR, j_tR, mR1, mgR, mR2)
        total += c1 * c2 * c3 * c4
    return total


def compute_vertex_3pt_single_nint(n_ext, j_ext, mj_ext, n_int, q_probe=1):
    """Compute the 3-pt vertex correction contribution from a single n_int level."""
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    jI_L = Rational(n_int + 1, 2)
    jI_R = Rational(n_int, 2)

    lam = Rational(2*n_int + 3, 2)
    lam4 = lam**4

    j_int_max = jI_L + jI_R
    j_int_min = abs(jI_L - jI_R)

    if not vertex_allowed(n_int, n_int, q_probe):
        return S.Zero

    probe_chs = get_channels(n_int, n_int, q_probe, jI_L, jI_R, jI_L, jI_R)
    if not probe_chs:
        return S.Zero

    total = S.Zero
    for j_int in half_ints(j_int_max):
        if j_int < j_int_min:
            continue
        for mj_int in half_ints(j_int):
            for mj_int_prime in half_ints(j_int):
                probe_amp = S.Zero
                for jpL, jpR in probe_chs:
                    for mpL in half_ints(jpL):
                        for mpR in half_ints(jpR):
                            pa = vertex_amp_pol(
                                jI_L, jI_R, j_int, mj_int,
                                jI_L, jI_R, j_int, mj_int_prime,
                                jpL, jpR, mpL, mpR)
                            probe_amp += pa
                if probe_amp == 0:
                    continue

                q_lo = max(1, abs(n_ext - n_int))
                q_hi = n_ext + n_int

                for q_loop in range(q_lo, q_hi + 1):
                    if not vertex_allowed(n_ext, n_int, q_loop):
                        continue
                    mu_q = Rational(q_loop * (q_loop + 2))

                    chs1 = get_channels(n_ext, n_int, q_loop, jE_L, jE_R, jI_L, jI_R)
                    chs2 = get_channels(n_int, n_ext, q_loop, jI_L, jI_R, jE_L, jE_R)

                    for jgL1, jgR1 in chs1:
                        for jgL2, jgR2 in chs2:
                            for mgL in half_ints(jgL1):
                                for mgR in half_ints(jgR1):
                                    if abs(mgL) > jgL2 or abs(mgR) > jgR2:
                                        continue
                                    v1 = vertex_amp_pol(
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jI_L, jI_R, j_int, mj_int,
                                        jgL1, jgR1, mgL, mgR)
                                    if v1 == 0:
                                        continue
                                    v3 = vertex_amp_pol(
                                        jI_L, jI_R, j_int, mj_int_prime,
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jgL2, jgR2, mgL, mgR)
                                    if v3 == 0:
                                        continue
                                    total += v1 * probe_amp * v3 / (lam4 * mu_q)

    return total


# ---------------------------------------------------------------------------
# Power-law tail correction (Euler-Maclaurin improved)
# ---------------------------------------------------------------------------

def fit_power_law(all_B: dict, n_min: int = 10):
    """
    Fit B(n) ~ C * n^(-p) for even and odd n separately.
    Returns (C_even, p_even, C_odd, p_odd).
    """
    even_n = np.array(sorted(n for n in all_B if n % 2 == 0 and n >= n_min), dtype=float)
    odd_n  = np.array(sorted(n for n in all_B if n % 2 == 1 and n >= n_min), dtype=float)

    if len(even_n) < 3 or len(odd_n) < 3:
        raise ValueError(f"Need >= 3 even and >= 3 odd data points above n_min={n_min}")

    even_B = np.array([all_B[int(n)] for n in even_n])
    odd_B  = np.array([all_B[int(n)] for n in odd_n])

    # log-log least squares: log|B| = log(C) - p*log(n)
    def logfit(ns, Bs):
        A = np.column_stack([np.ones_like(ns), -ns])  # [1, -log(n)]
        b = np.log(np.abs(Bs))
        x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        C = np.exp(x[0])
        p = x[1]
        return C, p

    C_e, p_e = logfit(np.log(even_n), even_B)
    C_o, p_o = logfit(np.log(odd_n), odd_B)

    return C_e, p_e, C_o, p_o


def power_law_tail_sum(n_start: int, C_e: float, p_e: float,
                       C_o: float, p_o: float, n_end: int = 2_000_000) -> float:
    """
    Approximate sum_{n=n_start}^{n_end} B(n) where B(n) ~ C_e*n^(-p_e) for even
    and B(n) ~ C_o*n^(-p_o) for odd.

    Uses the Hurwitz zeta function for exact summation of the power-law tail:
        sum_{k=0}^{inf} C * (n0 + k*2)^{-p} = C * 2^{-p} * zeta(p, n0/2)

    where the second argument is the Hurwitz zeta parameter.
    """
    # First even n >= n_start
    n_even_start = n_start if n_start % 2 == 0 else n_start + 1
    # First odd n >= n_start
    n_odd_start = n_start if n_start % 2 == 1 else n_start + 1

    # Even tail: sum_{k=0}^{inf} C_e * (n_even_start + 2k)^{-p_e}
    # = C_e * sum_{k=0}^{inf} (2*(n_even_start/2 + k))^{-p_e}
    # = C_e * 2^{-p_e} * zeta(p_e, n_even_start/2)
    if n_even_start > 0:
        tail_even = float(C_e * 2**(-p_e) *
                          float(mpmath.zeta(p_e, mpmath.mpf(n_even_start) / 2)))
    else:
        tail_even = 0.0

    # Odd tail: sum_{k=0}^{inf} C_o * (n_odd_start + 2k)^{-p_o}
    # = C_o * 2^{-p_o} * zeta(p_o, n_odd_start/2)
    if n_odd_start > 0:
        tail_odd = float(C_o * 2**(-p_o) *
                         float(mpmath.zeta(p_o, mpmath.mpf(n_odd_start) / 2)))
    else:
        tail_odd = 0.0

    return tail_even + tail_odd


def tail_uncertainty(n_start: int, C_e: float, p_e: float,
                     C_o: float, p_o: float, dp: float = 0.1) -> float:
    """
    Estimate uncertainty in the tail sum from +/-dp variation in the power law
    exponent. This is the main source of error when n_start is moderate.
    """
    tail_nominal = power_law_tail_sum(n_start, C_e, p_e, C_o, p_o)
    tail_hi = power_law_tail_sum(n_start, C_e, p_e - dp, C_o, p_o - dp)
    tail_lo = power_law_tail_sum(n_start, C_e, p_e + dp, C_o, p_o + dp)
    return max(abs(tail_hi - tail_nominal), abs(tail_lo - tail_nominal))


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def load_existing_data():
    """Load existing n_int=0..25 data for n_ext=1."""
    path = DATA_DIR / "g2_extended_nint_v2.json"
    with open(path) as f:
        data = json.load(f)

    all_B = {int(k): float(v) for k, v in data["all_B"].items()}
    V_mag_path = DATA_DIR / "alpha_g_minus_2_ratio_investigation.json"
    with open(V_mag_path) as f:
        ratio_data = json.load(f)
    V_mag = ratio_data["V_magnetic_float"]

    return all_B, V_mag


def extend_computation(all_B: dict, V_mag: float,
                       n_int_target: int = 50,
                       checkpoint_file: Path = None):
    """
    Extend the n_int sum from max(all_B) to n_int_target.
    Saves checkpoint after each new level.
    Returns updated all_B dict.
    """
    n_int_done = max(all_B.keys())
    if n_int_done >= n_int_target:
        print(f"  Already have n_int up to {n_int_done} (target={n_int_target}). Skipping.")
        return all_B

    n_ext = 1
    j_ext = Rational(1, 2)

    print(f"  Extending n_int from {n_int_done + 1} to {n_int_target}...", flush=True)

    for n_int in range(n_int_done + 1, n_int_target + 1):
        t0 = time.time()
        contrib_up = compute_vertex_3pt_single_nint(n_ext, j_ext, Rational(1, 2), n_int)
        contrib_dn = compute_vertex_3pt_single_nint(n_ext, j_ext, Rational(-1, 2), n_int)
        B_level = float(contrib_up - contrib_dn)
        elapsed = time.time() - t0

        all_B[n_int] = B_level
        cum_B = sum(all_B.values())
        F2_S = cum_B / V_mag / SCHWINGER
        delta = F2_S - 1.0
        print(f"  n_int={n_int:3d}: B={B_level:12.4e}  delta(raw)={delta:.12f}  "
              f"({elapsed:.1f}s)", flush=True)

        if checkpoint_file is not None:
            _save_checkpoint_data(all_B, checkpoint_file)

    return all_B


def _save_checkpoint_data(all_B: dict, path: Path):
    """Save all_B to a JSON checkpoint."""
    with open(path, 'w') as f:
        json.dump({str(k): v for k, v in sorted(all_B.items())}, f, indent=2)


def load_checkpoint(path: Path) -> dict:
    """Load all_B from a JSON checkpoint."""
    with open(path) as f:
        data = json.load(f)
    return {int(k): float(v) for k, v in data.items()}


# ---------------------------------------------------------------------------
# c_3 extraction
# ---------------------------------------------------------------------------

def extract_c3(all_B: dict, V_mag: float, n_int_cutoff: int = None):
    """
    Extract c_3 from the delta at n_ext=1, given the all_B table.

    Steps:
    1. Sum all_B up to n_int_cutoff (or max available)
    2. Apply Hurwitz-zeta tail correction
    3. Compute delta = F_2/Schwinger - 1
    4. c_3 = (delta - c_1/lambda^2 - c_2/lambda^4) * lambda^6

    Returns dict with all intermediate quantities and error estimates.
    """
    mpmath.mp.dps = 50

    if n_int_cutoff is None:
        n_int_cutoff = max(all_B.keys())

    # Sum up to cutoff
    B_keys = sorted(k for k in all_B if k <= n_int_cutoff)
    cum_B = sum(all_B[k] for k in B_keys)
    n_max_used = max(B_keys)

    # Power-law fit from the upper half of data (better fit in asymptotic regime)
    # Use data from n >= n_max_used//2 to n_max_used
    n_fit_start = max(10, n_max_used // 2)
    C_e, p_e, C_o, p_o = fit_power_law(
        {k: all_B[k] for k in B_keys}, n_min=n_fit_start)

    print(f"\n  Power-law fit (n >= {n_fit_start}):")
    print(f"    Even: C={C_e:.6e}, p={p_e:.6f}")
    print(f"    Odd:  C={C_o:.6e}, p={p_o:.6f}")

    # Hurwitz-zeta tail correction (n_int > n_max_used)
    tail_B = power_law_tail_sum(n_max_used + 1, C_e, p_e, C_o, p_o)
    tail_unc = tail_uncertainty(n_max_used + 1, C_e, p_e, C_o, p_o, dp=0.05)

    print(f"  Tail correction (n > {n_max_used}):")
    print(f"    tail_B = {tail_B:.6e}")
    print(f"    tail_unc = {tail_unc:.6e}")

    # Convert to F_2/Schwinger units
    total_B = cum_B + tail_B
    F2_S = total_B / V_mag / SCHWINGER
    delta = mpmath.mpf(str(F2_S - 1.0))

    # Tail uncertainty in delta units
    tail_unc_delta = mpmath.mpf(str(tail_unc / V_mag / SCHWINGER))

    print(f"\n  F_2/Schwinger:")
    print(f"    cum_B / (V_mag * Schwinger) = {float(cum_B / V_mag / SCHWINGER):.15f}")
    print(f"    tail correction = {float(tail_B / V_mag / SCHWINGER):.6e}")
    print(f"    F_2/S = {float(F2_S):.15f}")
    print(f"    delta = {float(delta):.15f}")

    # Subtract known c_1 and c_2 terms
    x = 1 / LAM2_N1  # = 4/25 = 0.16
    c1_term = C1 * x
    c2_term = C2_EXACT * x**2
    residual = delta - c1_term - c2_term
    c3_estimate = residual / x**3

    print(f"\n  Curvature extraction:")
    print(f"    lambda = {float(LAM_N1)}, x = 1/lambda^2 = {float(x)}")
    print(f"    c_1/lambda^2 = {float(c1_term):.15f}")
    print(f"    c_2/lambda^4 = {float(c2_term):.15f}")
    print(f"    residual = delta - c_1/lambda^2 - c_2/lambda^4 = {float(residual):.12e}")
    print(f"    c_3 = residual / x^3 = {float(c3_estimate):.8e}")
    print(f"    |c_3| = {abs(float(c3_estimate)):.6e}")

    # Error budget
    tail_error_in_c3 = tail_unc_delta / x**3
    c2_uncertainty_in_c3 = mpmath.mpf(0)  # c_2 is exact, no uncertainty
    print(f"\n  Error budget:")
    print(f"    Tail uncertainty -> c_3 error: {float(tail_error_in_c3):.4e}")
    print(f"    Total c_3 error estimate: {float(tail_error_in_c3):.4e}")
    print(f"    |c_3| / error: {abs(float(c3_estimate)) / abs(float(tail_error_in_c3)):.2f}")

    # Is c_3 consistent with zero?
    snr = abs(float(c3_estimate)) / max(float(tail_error_in_c3), 1e-50)
    is_consistent_with_zero = snr < 3.0
    print(f"\n  Is c_3 consistent with zero? {'YES (SNR < 3)' if is_consistent_with_zero else 'NO (SNR >= 3)'}")
    print(f"    Signal/noise = {snr:.2f}")

    return {
        "n_max_used": n_max_used,
        "n_fit_start": n_fit_start,
        "C_even": float(C_e), "p_even": float(p_e),
        "C_odd": float(C_o), "p_odd": float(p_o),
        "cum_B": float(cum_B),
        "tail_B": float(tail_B),
        "tail_unc_B": float(tail_unc),
        "F2_over_S": float(F2_S),
        "delta": float(delta),
        "delta_mpmath": str(delta),
        "c1_term": float(c1_term),
        "c2_term": float(c2_term),
        "residual": float(residual),
        "c3_estimate": float(c3_estimate),
        "c3_abs": abs(float(c3_estimate)),
        "tail_error_in_c3": float(tail_error_in_c3),
        "snr": snr,
        "consistent_with_zero": is_consistent_with_zero,
    }


# ---------------------------------------------------------------------------
# PSLQ identification of c_3 (if nonzero)
# ---------------------------------------------------------------------------

def pslq_c3(c3_val: float, c3_error: float):
    """
    Attempt PSLQ identification of c_3 against a basis of Paper 2 invariants
    and transcendental constants.

    Basis includes:
    - Rational constants: 1, 1/n for small n
    - pi-even: pi^2, pi^4, pi^6
    - Paper 2 invariants: BDelta=21/20, FDelta=pi^2/240, F/B=pi^2/252
    - Products: B*Delta^2, F*Delta^2, F^2/B, BDelta*FDelta, etc.
    - Odd zeta: zeta(3), zeta(5), pi^2*zeta(3) (T9: these should NOT appear at one loop)
    - Catalan G, beta(4) (two-loop, should not appear)
    """
    mpmath.mp.dps = 50
    c3_mp = mpmath.mpf(str(c3_val))

    print(f"\n  PSLQ identification of c_3 = {c3_val:.8e}")
    print(f"  (error ~ {c3_error:.2e})")

    if abs(c3_val) < 3 * c3_error:
        print("  c_3 is consistent with zero (SNR < 3) -- PSLQ skipped as unreliable")
        return {"identified": False, "reason": "consistent_with_zero"}

    # T9-consistent basis: rational + pi^{even} only
    B_mp = mpmath.mpf(42)
    F_mp = mpmath.pi**2 / 6
    D_mp = mpmath.mpf(1) / 40

    BD = B_mp * D_mp  # 21/20
    FD = F_mp * D_mp  # pi^2/240
    FB = F_mp / B_mp  # pi^2/252

    # Build rich basis over T9-consistent content
    basis_labels = [
        "1", "pi^2", "pi^4", "pi^6",
        "B*Delta", "F*Delta", "F/B",
        "B*Delta^2", "F*Delta^2", "F^2/B", "F/B^2",
        "(B*Delta)^2", "(F*Delta)^2",
        "B*Delta*F*Delta", "B*Delta*F/B", "F*Delta*F/B",
        "c2", "c2^2", "c1*c2",
        "1/B", "1/B^2", "Delta", "Delta^2",
    ]
    basis_vals = [
        mpmath.mpf(1), mpmath.pi**2, mpmath.pi**4, mpmath.pi**6,
        BD, FD, FB,
        B_mp * D_mp**2, F_mp * D_mp**2, F_mp**2 / B_mp, F_mp / B_mp**2,
        BD**2, FD**2,
        BD * FD, BD * FB, FD * FB,
        C2_EXACT, C2_EXACT**2, C1 * C2_EXACT,
        mpmath.mpf(1) / B_mp, mpmath.mpf(1) / B_mp**2, D_mp, D_mp**2,
    ]

    # Also add odd-zeta to detect violations of T9
    labels_with_odd = basis_labels + ["zeta(3)", "pi^2*zeta(3)", "zeta(5)", "Catalan_G"]
    vals_with_odd = basis_vals + [
        mpmath.zeta(3), mpmath.pi**2 * mpmath.zeta(3),
        mpmath.zeta(5), mpmath.catalan
    ]

    results = {}

    # Try increasingly large bases
    test_sets = [
        ("Small T9 basis", ["1", "pi^2", "pi^4", "c2", "B*Delta", "F*Delta", "F/B"]),
        ("Medium T9 basis", ["1", "pi^2", "pi^4", "pi^6", "c2", "c2^2",
                              "B*Delta", "F*Delta", "F/B", "B*Delta^2", "F*Delta^2", "F/B^2"]),
        ("Full T9 basis", basis_labels),
        ("T9 + odd-zeta", labels_with_odd),
    ]

    for set_name, label_set in test_sets:
        print(f"\n  --- {set_name} ---")

        # Filter to available labels
        idx = [i for i, l in enumerate(labels_with_odd) if l in label_set]
        sel_labels = [labels_with_odd[i] for i in idx]
        sel_vals = [vals_with_odd[i] for i in idx]

        basis = [c3_mp] + sel_vals

        try:
            # Use tight tolerance relative to c3_error
            tol = max(1e-10, c3_error / abs(c3_val) * 0.1)
            rel = mpmath.pslq(basis, tol=tol, maxcoeff=10000)
        except Exception as e:
            print(f"    PSLQ error: {e}")
            rel = None

        if rel is None or rel[0] == 0:
            print(f"    No relation found (basis size {len(basis)})")
            results[set_name] = {"identified": False}
        else:
            # Reconstruct: c3 = -sum(rel[i+1] * val[i]) / rel[0]
            c3_reconstructed = mpmath.mpf(0)
            terms = []
            for i, (r, l, v) in enumerate(zip(rel[1:], sel_labels, sel_vals)):
                if r != 0:
                    coeff = mpmath.mpf(-r) / mpmath.mpf(rel[0])
                    c3_reconstructed += coeff * v
                    terms.append(f"({-r}/{rel[0]})*{l}")
            residual = abs(c3_mp - c3_reconstructed)
            print(f"    FOUND: c_3 = {' + '.join(terms)}")
            print(f"    Residual: {float(residual):.4e}")
            print(f"    Match digits: {float(-mpmath.log10(residual/abs(c3_mp))):.1f}")
            results[set_name] = {
                "identified": True,
                "terms": terms,
                "residual": float(residual),
                "c3_reconstructed": float(c3_reconstructed),
            }

    return results


# ---------------------------------------------------------------------------
# Convergence study: c_3 vs n_int_cutoff
# ---------------------------------------------------------------------------

def c3_convergence_study(all_B: dict, V_mag: float):
    """
    Track c_3 as a function of the n_int cutoff to see if it converges.
    Uses the same tail correction formula at each cutoff.
    """
    print("\n" + "=" * 70)
    print("  c_3 CONVERGENCE STUDY")
    print("=" * 70)

    n_max = max(all_B.keys())
    cutoffs = [n for n in range(15, n_max + 1, 5)]
    if n_max not in cutoffs:
        cutoffs.append(n_max)

    x = float(1 / LAM2_N1)  # = 0.16

    results = []
    for cutoff in cutoffs:
        B_keys = sorted(k for k in all_B if k <= cutoff)
        cum_B = sum(all_B[k] for k in B_keys)
        n_fit_start = max(10, cutoff // 2)

        try:
            C_e, p_e, C_o, p_o = fit_power_law(
                {k: all_B[k] for k in B_keys}, n_min=n_fit_start)
            tail_B = power_law_tail_sum(cutoff + 1, C_e, p_e, C_o, p_o)
            tail_unc = tail_uncertainty(cutoff + 1, C_e, p_e, C_o, p_o, dp=0.05)
        except Exception:
            tail_B = 0.0
            tail_unc = float('nan')
            C_e = p_e = C_o = p_o = 0.0

        total_B = cum_B + tail_B
        F2_S = total_B / V_mag / SCHWINGER
        delta = F2_S - 1.0

        c1_term = float(C1) * x
        c2_term = float(C2_EXACT) * x**2
        residual = delta - c1_term - c2_term
        c3_est = residual / x**3

        tail_unc_delta = tail_unc / V_mag / SCHWINGER
        c3_error = tail_unc_delta / x**3

        print(f"  n_max={cutoff:3d}: delta={delta:.15f}  c_3={c3_est:.6e}  "
              f"+/-{c3_error:.2e}  (p_e={p_e:.3f}, p_o={p_o:.3f})")

        results.append({
            "n_int_cutoff": cutoff,
            "cum_B": float(cum_B),
            "tail_B": float(tail_B),
            "tail_unc": float(tail_unc),
            "delta": float(delta),
            "c3_estimate": float(c3_est),
            "c3_error": float(c3_error),
            "C_even": float(C_e), "p_even": float(p_e),
            "C_odd": float(C_o), "p_odd": float(p_o),
        })

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("  g-2 c_3 INVESTIGATION: S^3 curvature expansion at n_ext=1")
    print("=" * 70)

    mpmath.mp.dps = 50

    # --- 0. Verify c_2 ---
    print(f"\n  Known constants:")
    print(f"    c_1 = 1/2 = {float(C1):.15f}")
    print(f"    c_2 = (2 - BDelta - FDelta - F/B)/5 = 19/100 - 41pi^2/25200")
    print(f"       = {float(C2_EXACT):.15f}")
    c2_diff = abs(C2_EXACT - C2_RATIONAL)
    print(f"    Verify c_2 = 19/100 - 41pi^2/25200: diff = {float(c2_diff):.4e}")
    assert float(c2_diff) < 1e-40, "c_2 formula mismatch!"
    print(f"    VERIFIED.")
    print(f"    lambda(n_ext=1) = 5/2, x = 1/lambda^2 = 4/25 = 0.16")
    print(f"    x^3 = (4/25)^3 = 64/15625 = {float(1/LAM6_N1):.10f}")

    # --- 1. Load existing data ---
    print(f"\n{'=' * 70}")
    print(f"  Section 1: Loading existing data (n_int=0..25)")
    print(f"{'=' * 70}")
    all_B, V_mag = load_existing_data()
    n_existing = max(all_B.keys())
    print(f"  Loaded n_int=0..{n_existing}, V_mag = {V_mag:.12e}")
    print(f"  cum_B = {sum(all_B.values()):.12e}")

    # Checkpoint file for new data
    checkpoint_file = DATA_DIR / "g2_c3_all_B_checkpoint.json"

    # If checkpoint exists, load it (in case we're resuming)
    if checkpoint_file.exists():
        checkpoint_B = load_checkpoint(checkpoint_file)
        if max(checkpoint_B.keys()) > n_existing:
            print(f"  Resuming from checkpoint: n_int up to {max(checkpoint_B.keys())}")
            all_B.update(checkpoint_B)
            n_existing = max(all_B.keys())

    # --- 2. Extend computation ---
    print(f"\n{'=' * 70}")
    print(f"  Section 2: Extending to n_int=50")
    print(f"{'=' * 70}")

    N_INT_TARGET = 50
    if n_existing < N_INT_TARGET:
        print(f"  Need to compute n_int={n_existing+1}..{N_INT_TARGET}")
        print(f"  Estimated time: ~{(N_INT_TARGET - n_existing) * 30:.0f}s "
              f"(~30s per level at n_int~30-50)", flush=True)
        all_B = extend_computation(all_B, V_mag, N_INT_TARGET, checkpoint_file)
    else:
        print(f"  Already have n_int up to {n_existing}. Proceeding.")

    n_max = max(all_B.keys())
    print(f"\n  Final dataset: n_int=0..{n_max}")

    # --- 3. c_3 extraction at full precision ---
    print(f"\n{'=' * 70}")
    print(f"  Section 3: c_3 extraction at n_max={n_max}")
    print(f"{'=' * 70}")
    extraction = extract_c3(all_B, V_mag, n_int_cutoff=n_max)

    # --- 4. Convergence study ---
    conv_results = c3_convergence_study(all_B, V_mag)

    # --- 5. PSLQ identification ---
    print(f"\n{'=' * 70}")
    print(f"  Section 4: PSLQ identification of c_3")
    print(f"{'=' * 70}")
    c3_val = extraction["c3_estimate"]
    c3_err = extraction["tail_error_in_c3"]
    pslq_results = pslq_c3(c3_val, c3_err)

    # --- 6. Summary ---
    print(f"\n{'=' * 70}")
    print(f"  SUMMARY")
    print(f"{'=' * 70}")
    print(f"  n_int summed: 0..{n_max}")
    print(f"  delta (tail-corrected) = {extraction['delta']:.15f}")
    print(f"  c_1 = 1/2 (exact)")
    print(f"  c_2 = {float(C2_EXACT):.10f} (exact)")
    print(f"  c_3 = {c3_val:.8e} +/- {c3_err:.2e}")
    print(f"  |c_3| consistent with zero: {extraction['consistent_with_zero']}")
    if not extraction['consistent_with_zero']:
        print(f"  Signal/noise = {extraction['snr']:.2f}")

    # Provide context: c_3*x^3 contribution to delta
    x = float(1 / LAM2_N1)
    c3_contribution_to_delta = c3_val * x**3
    print(f"\n  c_3 * x^3 = {c3_contribution_to_delta:.4e} (contribution to delta)")
    print(f"  delta ~= {extraction['delta']:.8f}")
    print(f"  Relative contribution: {abs(c3_contribution_to_delta/extraction['delta']):.4e}")

    # --- 7. Save results ---
    print(f"\n{'=' * 70}")
    print(f"  Section 5: Saving results")
    print(f"{'=' * 70}")

    output = {
        "description": "c3 curvature coefficient of S3 g-2 expansion",
        "n_ext": 1,
        "lambda_n1": 2.5,
        "x_n1": float(1 / LAM2_N1),
        "c1_exact": "1/2",
        "c2_exact": "19/100 - 41*pi^2/25200",
        "c2_numerical": float(C2_EXACT),
        "B_hopf": B_HOPF,
        "F_hopf": float(F_HOPF),
        "Delta_hopf": float(DELTA_HOPF),
        "n_int_max": n_max,
        "V_mag": V_mag,
        "extraction": extraction,
        "convergence_study": conv_results,
        "pslq_results": str(pslq_results),
        "all_B": {str(k): v for k, v in sorted(all_B.items())},
        "summary": {
            "c3_estimate": c3_val,
            "c3_error": c3_err,
            "c3_consistent_with_zero": extraction["consistent_with_zero"],
            "c3_snr": extraction["snr"],
            "delta_corrected": extraction["delta"],
        }
    }

    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"  Saved to {OUTPUT_FILE}")

    return output


if __name__ == "__main__":
    result = main()
    print("\nDone.")
