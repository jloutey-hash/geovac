"""Track RH-M: Spectral-side explicit-formula / zero statistics.

Sprint 3 of the RH-directed series. Sprint 2 Track RH-G ruled out the
Ihara-side Hilbert-Polya candidate (Hashimoto eigenvalue spacings are
Poisson, CV ~ 1, not GUE CV ~ 0.42). This track asks: do the ZEROS of
the SPECTRAL Dirichlet series D(s), D_even(s), D_odd(s) -- as functions
of a complex variable s -- have RH-like statistics?

Functions
---------
D(s)      = 2 * zeta(s-2, 3/2) - (1/2) * zeta(s, 3/2)
D_even(s) = 2^{-s} [ 8 * zeta(s-2, 3/4) - (1/2) * zeta(s, 3/4) ]
D_odd(s)  = 2^{-s} [ 8 * zeta(s-2, 5/4) - (1/2) * zeta(s, 5/4) ]

where zeta(s, a) is the Hurwitz zeta function (analytic continuation in s).

The series converge for Re(s) > 3 (from zeta(s-2, a)); the analytic
continuation has simple poles at s = 1 and s = 3 where the Hurwitz
zetas have their pole.

Method
------
1. Verify pole structure and known special values.
2. Use the argument principle on a rectangular contour to count zeros in
   horizontal strips.
3. For each strip that contains zeros, seed mpmath.findroot and refine.
4. Accumulate zeros, compute normalized spacings (with local unfolding),
   and test the spacing distribution against GUE, Poisson, and picket-fence
   predictions.

References
----------
- Paper 28 (papers/observations/paper_28_qed_s3.tex): vertex parity, D_even/D_odd.
- debug/riemann_limit_memo.md: Sprint 2 RH-G structural redirection.
- geovac/qed_vertex.py: existing Hurwitz implementations.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import mpmath


# ---------------------------------------------------------------------------
# Dirichlet series on the Dirac spectrum of S^3
# ---------------------------------------------------------------------------

def D_full(s) -> mpmath.mpc:
    """D(s) = 2 * zeta(s-2, 3/2) - (1/2) * zeta(s, 3/2).

    Dirac-on-S^3 Dirichlet series, analytic continuation in s via Hurwitz.
    Poles: simple poles at s=1 (from zeta(s, 3/2)) and s=3 (from zeta(s-2, 3/2)).
    """
    s = mpmath.mpc(s)
    return 2 * mpmath.hurwitz(s - 2, mpmath.mpf(3) / 2) \
        - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(3) / 2)


def D_even(s) -> mpmath.mpc:
    """D_even(s) = 2^{-s} * [ 8 * zeta(s-2, 3/4) - (1/2) * zeta(s, 3/4) ].

    Even-n Dirac modes, analytic continuation via Hurwitz zeta at shift 3/4.
    Poles at s=1 and s=3. The factor 2^{-s} is entire and non-zero.
    """
    s = mpmath.mpc(s)
    two_s = mpmath.power(2, -s)
    return two_s * (8 * mpmath.hurwitz(s - 2, mpmath.mpf(3) / 4)
                    - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(3) / 4))


def D_odd(s) -> mpmath.mpc:
    """D_odd(s) = 2^{-s} * [ 8 * zeta(s-2, 5/4) - (1/2) * zeta(s, 5/4) ].

    Odd-n Dirac modes, analytic continuation via Hurwitz zeta at shift 5/4.
    Poles at s=1 and s=3. The factor 2^{-s} is entire and non-zero.
    """
    s = mpmath.mpc(s)
    two_s = mpmath.power(2, -s)
    return two_s * (8 * mpmath.hurwitz(s - 2, mpmath.mpf(5) / 4)
                    - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(5) / 4))


# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

def sanity_checks(dps: int = 50) -> Dict[str, object]:
    """Verify the three functions evaluate correctly on known points.

    - D(4) = pi^2 - pi^4 / 12 (Paper 28, Table 1).
    - D_even(4) + D_odd(4) = D(4).
    - D_even(4) = pi^2/2 - pi^4/24 - 4G + 4*beta(4) (Paper 28 Eq. 13).
    - Residue at simple poles s=1, s=3.
    """
    mpmath.mp.dps = dps
    checks: Dict[str, object] = {}

    pi2 = mpmath.pi ** 2
    pi4 = mpmath.pi ** 4

    # D(4)
    ref_D4 = pi2 - pi4 / 12
    val_D4 = D_full(4)
    checks["D(4)"] = {
        "computed": str(mpmath.nstr(val_D4, 30)),
        "reference (pi^2 - pi^4/12)": str(mpmath.nstr(ref_D4, 30)),
        "abs error": float(abs(val_D4 - ref_D4)),
    }

    # Sum check
    sum_eo = D_even(4) + D_odd(4)
    checks["D_even(4) + D_odd(4) vs D(4)"] = {
        "sum": str(mpmath.nstr(sum_eo, 30)),
        "D(4)": str(mpmath.nstr(val_D4, 30)),
        "abs error": float(abs(sum_eo - val_D4)),
    }

    # D_even(4) per Paper 28 Eq. 13
    G = mpmath.catalan
    beta4 = (mpmath.hurwitz(4, mpmath.mpf(1) / 4)
             - mpmath.hurwitz(4, mpmath.mpf(3) / 4)) / mpmath.power(4, 4)
    ref_Deven4 = pi2 / 2 - pi4 / 24 - 4 * G + 4 * beta4
    val_Deven4 = D_even(4)
    checks["D_even(4) (Paper 28 Eq. 13)"] = {
        "computed": str(mpmath.nstr(val_Deven4, 30)),
        "reference": str(mpmath.nstr(ref_Deven4, 30)),
        "abs error": float(abs(val_Deven4 - ref_Deven4)),
    }

    # Residues (pole at s=3, residue 2 for D_full)
    for fn_name, fn in [("D_full", D_full), ("D_even", D_even), ("D_odd", D_odd)]:
        s_near_3 = mpmath.mpc(3, mpmath.mpf("1e-12"))
        val = fn(s_near_3)
        res_est = val * (s_near_3 - 3)
        checks[f"{fn_name} residue at s=3 (limit)"] = {
            "(s-3)*fn(s) at s=3+1e-12i": str(mpmath.nstr(res_est, 12)),
        }
        s_near_1 = mpmath.mpc(1, mpmath.mpf("1e-10"))
        val = fn(s_near_1)
        res_est = val * (s_near_1 - 1)
        checks[f"{fn_name} residue at s=1 (limit)"] = {
            "(s-1)*fn(s) at s=1+1e-10i": str(mpmath.nstr(res_est, 12)),
        }

    return checks


# ---------------------------------------------------------------------------
# Argument principle: count zeros in a box
# ---------------------------------------------------------------------------

def count_zeros_in_box(
    fn: Callable,
    re_lo: float, re_hi: float,
    im_lo: float, im_hi: float,
    n_per_edge: int = 60,
    pole_shift: float = 0.0,
) -> int:
    """Count zeros minus poles of `fn` inside a rectangle via argument principle.

    Uses discrete arg-unwrapped contour integration. The box must not pass
    through the poles at s=1 or s=3.
    """
    pts = []
    for i in range(n_per_edge):
        pts.append(mpmath.mpc(re_lo + (re_hi - re_lo) * i / n_per_edge, im_lo))
    for i in range(n_per_edge):
        pts.append(mpmath.mpc(re_hi, im_lo + (im_hi - im_lo) * i / n_per_edge))
    for i in range(n_per_edge):
        pts.append(mpmath.mpc(re_hi - (re_hi - re_lo) * i / n_per_edge, im_hi))
    for i in range(n_per_edge):
        pts.append(mpmath.mpc(re_lo, im_hi - (im_hi - im_lo) * i / n_per_edge))
    pts.append(pts[0])

    total = mpmath.mpf(0)
    for i in range(len(pts) - 1):
        v1 = fn(pts[i])
        v2 = fn(pts[i + 1])
        da = mpmath.arg(v2) - mpmath.arg(v1)
        while da > mpmath.pi:
            da -= 2 * mpmath.pi
        while da < -mpmath.pi:
            da += 2 * mpmath.pi
        total += da

    n_zeros = float(total / (2 * mpmath.pi))
    return int(round(n_zeros))


# ---------------------------------------------------------------------------
# Zero refinement from a seed
# ---------------------------------------------------------------------------

def _try_findroot(fn, seed, tol, maxsteps) -> complex:
    """Wrapper that returns None if findroot fails."""
    try:
        z = mpmath.findroot(fn, seed, tol=tol, maxsteps=maxsteps)
        val = fn(z)
        if abs(val) > float(tol) * 1e3:
            return None
        return z
    except Exception:
        return None


def find_zeros_in_strip(
    fn: Callable,
    re_range: Tuple[float, float],
    im_range: Tuple[float, float],
    *,
    n_seeds_re: int = 9,
    n_seeds_im: int = 100,
    pole_near: List[Tuple[float, float]] = None,
    dedup_tol: float = 1e-6,
    tol_zero: float = 1e-25,
    verbose: bool = False,
) -> List[mpmath.mpc]:
    """Find zeros of `fn` in a rectangle by seeding findroot.

    For each seed point in a (n_seeds_re x n_seeds_im) grid, run
    mpmath.findroot; keep unique converged zeros. This differs from
    purely relying on the argument principle -- it provides the actual
    zero locations.

    Parameters
    ----------
    fn : callable
        Function to zero-find.
    re_range, im_range : tuple of float
        Rectangle boundaries.
    n_seeds_re, n_seeds_im : int
        Seed grid density.
    pole_near : list of (re, im) float tuples
        Avoid seeding too close to these points.
    dedup_tol : float
        Distance for treating two zeros as identical.
    tol_zero : float
        mpmath.findroot tolerance.
    """
    if pole_near is None:
        pole_near = [(1.0, 0.0), (3.0, 0.0)]

    re_lo, re_hi = re_range
    im_lo, im_hi = im_range

    # Seed grid
    re_seeds = [re_lo + (re_hi - re_lo) * i / max(n_seeds_re - 1, 1)
                for i in range(n_seeds_re)]
    im_seeds = [im_lo + (im_hi - im_lo) * j / max(n_seeds_im - 1, 1)
                for j in range(n_seeds_im)]

    def too_close_to_pole(s):
        for pr, pi in pole_near:
            if abs(mpmath.re(s) - pr) < 0.1 and abs(mpmath.im(s) - pi) < 0.1:
                return True
        return False

    found: List[mpmath.mpc] = []

    for i, re_seed in enumerate(re_seeds):
        for j, im_seed in enumerate(im_seeds):
            seed = mpmath.mpc(re_seed, im_seed)
            if too_close_to_pole(seed):
                continue

            z = _try_findroot(fn, seed, tol_zero, maxsteps=60)
            if z is None:
                continue

            # Check bounds (loose)
            if not (re_lo - 1.0 <= float(mpmath.re(z)) <= re_hi + 1.0):
                continue
            if not (im_lo - 1.0 <= float(mpmath.im(z)) <= im_hi + 1.0):
                continue

            if too_close_to_pole(z):
                continue

            # Skip the already-known trivial integer zeros at s = 0, -2, -4, ...
            # (integer real with |Im|<1e-8)
            if abs(float(mpmath.im(z))) < 1e-5:
                # Treat as real-axis zero
                pass  # we keep real-axis zeros for now (they count)

            # Dedup
            is_dup = False
            for zf in found:
                if abs(z - zf) < dedup_tol:
                    is_dup = True
                    break
            if not is_dup:
                found.append(z)
                if verbose:
                    print(f"    seed ({re_seed:+.3f},{im_seed:+.3f}) -> "
                          f"z=({float(mpmath.re(z)):+.6f},{float(mpmath.im(z)):+.6f})")

    return found


# ---------------------------------------------------------------------------
# Zero hunting strategy using argument principle + refinement
# ---------------------------------------------------------------------------

def hunt_zeros(
    fn: Callable,
    re_range: Tuple[float, float] = (-4.0, 5.0),
    im_max: float = 60.0,
    strip_height: float = 2.0,
    *,
    verbose: bool = False,
) -> List[mpmath.mpc]:
    """Hunt zeros by: (1) counting per strip via arg principle,
    (2) seeding findroot within strips that contain zeros.

    Returns list of zeros with Im > 0 (upper half plane). The real-axis
    zeros and their conjugates are handled separately.
    """
    re_lo, re_hi = re_range

    # Skim trivial real-axis zeros first (s = 0, -2, -4, ..., -10 in range)
    # and the simple zero at s=0.
    real_zeros = []
    for r in range(-20, 2, 2):
        # Zero at negative even integer
        if re_lo - 0.5 <= r <= re_hi + 0.5:
            val = fn(r)
            if abs(val) < 1e-8:
                real_zeros.append(mpmath.mpc(r, 0))
    # s=0 special zero
    if re_lo <= 0 <= re_hi:
        val = fn(0)
        if abs(val) < 1e-8:
            real_zeros.append(mpmath.mpc(0, 0))

    zeros: List[mpmath.mpc] = []

    im_lo = 0.2  # keep off the real axis; real zeros are counted separately
    while im_lo < im_max:
        im_hi = min(im_lo + strip_height, im_max)
        try:
            n_exp = count_zeros_in_box(fn, re_lo, re_hi, im_lo, im_hi,
                                       n_per_edge=40)
        except Exception:
            im_lo = im_hi
            continue

        if verbose:
            print(f"    Im in [{im_lo:5.2f}, {im_hi:5.2f}]: expected zeros = {n_exp}")

        if n_exp > 0:
            # Seed findroot densely in this strip
            seeds_im = 5 * max(n_exp, 1) + 5
            strip_zeros = find_zeros_in_strip(
                fn,
                re_range=(re_lo, re_hi),
                im_range=(im_lo, im_hi),
                n_seeds_re=5,
                n_seeds_im=seeds_im,
                verbose=False,
            )
            # Keep only those in upper half with Im > 0
            kept = [z for z in strip_zeros
                    if float(mpmath.im(z)) > 0.05
                    and im_lo - 0.5 <= float(mpmath.im(z)) <= im_hi + 0.5]
            for z in kept:
                is_dup = any(abs(z - zf) < 1e-4 for zf in zeros)
                if not is_dup:
                    zeros.append(z)
            if verbose:
                print(f"      Found {len(kept)} zeros, cumulative total {len(zeros)}")

        im_lo = im_hi

    zeros.sort(key=lambda z: float(mpmath.im(z)))
    return zeros, real_zeros


# ---------------------------------------------------------------------------
# Spacing statistics
# ---------------------------------------------------------------------------

def compute_spacing_stats(
    zeros: List[mpmath.mpc],
    *,
    unfold_window: int = 15,
) -> Dict[str, object]:
    """Compute pair-correlation and spacing CV for a sorted list of zeros.

    Sort by Im(s) ascending; compute consecutive spacings. Unfold by dividing
    each spacing by the local mean (window-sized moving average).
    """
    if len(zeros) < 2:
        return {"n_zeros": len(zeros), "cv": None,
                "error": "insufficient zeros"}

    imag_parts = sorted([float(mpmath.im(z)) for z in zeros
                         if float(mpmath.im(z)) > 0.05])
    # Dedup very close zeros
    cleaned = []
    for y in imag_parts:
        if not cleaned or abs(y - cleaned[-1]) > 1e-3:
            cleaned.append(y)
    imag_parts = cleaned

    if len(imag_parts) < 2:
        return {"n_zeros": len(imag_parts), "cv": None,
                "error": "insufficient imag-positive zeros"}

    raw_spacings = [imag_parts[i + 1] - imag_parts[i]
                    for i in range(len(imag_parts) - 1)]

    # Unfold: divide each spacing by local mean
    if len(raw_spacings) >= unfold_window:
        normalized = []
        for i in range(len(raw_spacings)):
            lo = max(0, i - unfold_window // 2)
            hi = min(len(raw_spacings), i + unfold_window // 2 + 1)
            local_mean = sum(raw_spacings[lo:hi]) / (hi - lo)
            if local_mean > 0:
                normalized.append(raw_spacings[i] / local_mean)
        spacings = normalized
    else:
        mean_raw = sum(raw_spacings) / len(raw_spacings)
        spacings = [s / mean_raw for s in raw_spacings] if mean_raw > 0 else []

    if not spacings:
        return {"n_zeros": len(imag_parts), "cv": None,
                "error": "spacings all non-positive"}

    mean_s = sum(spacings) / len(spacings)
    var_s = sum((s - mean_s) ** 2 for s in spacings) / len(spacings)
    std_s = math.sqrt(var_s)
    cv = std_s / mean_s if mean_s > 0 else float("inf")

    # Histogram
    bin_edges = [i * 0.2 for i in range(16)]
    counts = [0] * (len(bin_edges) - 1)
    for s in spacings:
        if s < 0:
            continue
        for b in range(len(bin_edges) - 1):
            if bin_edges[b] <= s < bin_edges[b + 1]:
                counts[b] += 1
                break
    n_total = sum(counts)
    densities = [c / (n_total * 0.2) if n_total > 0 else 0 for c in counts]

    # KS tests
    def poisson_cdf(s):
        return 1.0 - math.exp(-s) if s >= 0 else 0.0

    def gue_cdf(s):
        if s < 0:
            return 0.0
        n_steps = 300
        ds = s / n_steps
        total = 0.0
        for k in range(n_steps):
            x = (k + 0.5) * ds
            total += (32.0 / (math.pi ** 2)) * x ** 2 * \
                math.exp(-4.0 * x ** 2 / math.pi) * ds
        return total

    sorted_s = sorted([s for s in spacings if s >= 0])
    n = len(sorted_s)
    ks_poisson = 0.0
    ks_gue = 0.0
    for i, s in enumerate(sorted_s):
        emp_cdf = (i + 1) / n
        ks_poisson = max(ks_poisson, abs(emp_cdf - poisson_cdf(s)))
        ks_gue = max(ks_gue, abs(emp_cdf - gue_cdf(s)))

    def ks_pvalue(d, n):
        sqrt_n = math.sqrt(n)
        lam = (sqrt_n + 0.12 + 0.11 / sqrt_n) * d
        if lam < 0.05:
            return 1.0
        total = 0.0
        for k in range(1, 100):
            total += 2 * (-1) ** (k - 1) * math.exp(-2 * k ** 2 * lam ** 2)
        return max(0.0, min(1.0, total))

    return {
        "n_zeros": len(imag_parts),
        "imag_parts_head": imag_parts[:20],
        "raw_mean_spacing": sum(raw_spacings) / len(raw_spacings),
        "normalized_mean_spacing": mean_s,
        "normalized_std_spacing": std_s,
        "cv": cv,
        "bin_edges": bin_edges,
        "histogram_counts": counts,
        "histogram_density": densities,
        "ks_poisson_stat": ks_poisson,
        "ks_poisson_pvalue": ks_pvalue(ks_poisson, n),
        "ks_gue_stat": ks_gue,
        "ks_gue_pvalue": ks_pvalue(ks_gue, n),
        "n_spacings": len(spacings),
    }


# ---------------------------------------------------------------------------
# Synthetic calibration
# ---------------------------------------------------------------------------

def synthetic_poisson_spacings(n: int, seed: int = 42) -> List[float]:
    import random
    random.seed(seed)
    return [-math.log(random.random()) for _ in range(n)]


def synthetic_gue_spacings(n: int, seed: int = 42) -> List[float]:
    import random
    random.seed(seed)

    def gue_cdf(s):
        if s < 0:
            return 0.0
        n_steps = 500
        ds = s / n_steps
        total = 0.0
        for k in range(n_steps):
            x = (k + 0.5) * ds
            total += (32.0 / (math.pi ** 2)) * x ** 2 * \
                math.exp(-4.0 * x ** 2 / math.pi) * ds
        return total

    spacings = []
    for _ in range(n):
        u = random.random()
        lo, hi = 0.0, 6.0
        for _ in range(80):
            mid = (lo + hi) / 2
            if gue_cdf(mid) < u:
                lo = mid
            else:
                hi = mid
        spacings.append((lo + hi) / 2)
    return spacings


def _cv_of_list(lst: List[float]) -> float:
    mean = sum(lst) / len(lst)
    var = sum((x - mean) ** 2 for x in lst) / len(lst)
    return math.sqrt(var) / mean if mean > 0 else float("inf")


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main(dps: int = 35,
         im_max: float = 70.0,
         strip_height: float = 3.0,
         re_range: Tuple[float, float] = (-4.0, 5.0)):
    print("=" * 78)
    print("Track RH-M: Spectral-side Dirichlet zero statistics")
    print("=" * 78)
    t0 = time.time()

    mpmath.mp.dps = dps

    results: Dict[str, object] = {
        "mpmath_dps": dps,
        "im_max": im_max,
        "strip_height": strip_height,
        "re_range": re_range,
    }

    # ===== SANITY =====
    print("\n[1/4] Sanity checks")
    checks = sanity_checks(dps=dps)
    results["sanity_checks"] = checks
    # Print key ones
    for k in ["D(4)", "D_even(4) + D_odd(4) vs D(4)", "D_even(4) (Paper 28 Eq. 13)"]:
        v = checks[k]
        print(f"  {k}: abs_error = {v['abs error']:.3e}")

    # ===== SYNTHETIC CALIBRATION =====
    print("\n[2/4] Synthetic CV calibration")
    poi = synthetic_poisson_spacings(500)
    gue = synthetic_gue_spacings(300)
    cv_poi = _cv_of_list(poi)
    cv_gue = _cv_of_list(gue)
    print(f"  Synthetic Poisson (n=500): CV = {cv_poi:.4f}  (target ~1.00)")
    print(f"  Synthetic GUE     (n=300): CV = {cv_gue:.4f}  (target ~0.42)")
    results["synthetic_calibration"] = {
        "poisson_cv": cv_poi,
        "gue_cv": cv_gue,
        "poisson_target": 1.00,
        "gue_target": 0.42,
        "picket_fence_target": 0.0,
    }

    # ===== ZERO FINDING =====
    print(f"\n[3/4] Zero hunting in box [Re={re_range}, Im=0..{im_max}]")
    print(f"       using argument principle with strip height {strip_height}")

    fn_records = [
        ("D_full", D_full),
        ("D_even", D_even),
        ("D_odd", D_odd),
    ]

    fn_zeros_upper: Dict[str, List[mpmath.mpc]] = {}
    fn_real_zeros: Dict[str, List[mpmath.mpc]] = {}
    fn_box_counts: Dict[str, int] = {}

    for fn_name, fn in fn_records:
        print(f"\n  {fn_name}:")
        t_fn = time.time()
        # Full-box winding number first (sanity)
        full_n = count_zeros_in_box(fn, re_range[0], re_range[1], 0.1, im_max,
                                    n_per_edge=100)
        print(f"    Full-box zero count [Re={re_range}, Im=0.1..{im_max}]: {full_n}")
        zeros_upper, zeros_real = hunt_zeros(
            fn, re_range=re_range, im_max=im_max,
            strip_height=strip_height, verbose=False,
        )
        fn_zeros_upper[fn_name] = zeros_upper
        fn_real_zeros[fn_name] = zeros_real
        fn_box_counts[fn_name] = full_n
        print(f"    Upper zeros found: {len(zeros_upper)},  real zeros: {len(zeros_real)}")
        print(f"    First 5 upper zeros (sorted by Im):")
        for z in zeros_upper[:5]:
            print(f"      ({float(mpmath.re(z)):+.5f}, {float(mpmath.im(z)):+.5f})")
        print(f"    Last 3 upper zeros:")
        for z in zeros_upper[-3:]:
            print(f"      ({float(mpmath.re(z)):+.5f}, {float(mpmath.im(z)):+.5f})")
        print(f"    Time: {time.time() - t_fn:.1f}s")

    results["zeros"] = {
        fn_name: [(float(mpmath.re(z)), float(mpmath.im(z))) for z in zeros]
        for fn_name, zeros in fn_zeros_upper.items()
    }
    results["real_zeros"] = {
        fn_name: [(float(mpmath.re(z)), float(mpmath.im(z))) for z in zeros]
        for fn_name, zeros in fn_real_zeros.items()
    }
    results["box_winding_counts"] = fn_box_counts

    # ===== STATISTICS =====
    print("\n[4/4] Spacing statistics")
    stats_by_fn = {}
    for fn_name in ["D_full", "D_even", "D_odd"]:
        zeros = fn_zeros_upper[fn_name]
        stats = compute_spacing_stats(zeros, unfold_window=15)
        stats_by_fn[fn_name] = stats
        cv = stats.get("cv")
        n_z = stats.get("n_zeros")
        print(f"  {fn_name}: {n_z} zeros, CV = {cv if cv is None else f'{cv:.4f}'}")
        if cv is not None:
            print(f"    KS vs Poisson: D={stats['ks_poisson_stat']:.4f}, "
                  f"p={stats['ks_poisson_pvalue']:.4g}")
            print(f"    KS vs GUE:     D={stats['ks_gue_stat']:.4f}, "
                  f"p={stats['ks_gue_pvalue']:.4g}")

    results["statistics"] = stats_by_fn

    # ===== VERDICT =====
    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    verdict_rows = []
    for fn_name, stats in stats_by_fn.items():
        cv = stats.get("cv")
        nz = stats.get("n_zeros")
        if cv is None:
            cat = "insufficient data"
        elif cv < 0.25:
            cat = "picket-fence / rigid"
        elif cv < 0.55:
            cat = "GUE-like (RH candidate)"
        elif cv < 0.85:
            cat = "semi-Poisson / intermediate"
        else:
            cat = "Poisson / integrable"
        verdict_rows.append((fn_name, nz, cv, cat))
        cv_str = "N/A" if cv is None else f"{cv:.4f}"
        print(f"  {fn_name:10s} | n={nz or 0:4d} | CV={cv_str:>8s} | {cat}")

    spectral_hp_viable = any(
        v[2] is not None and 0.30 < v[2] < 0.55
        for v in verdict_rows
    )
    print(f"\nSpectral-side Hilbert-Polya viable: {spectral_hp_viable}")

    results["verdict"] = {
        "rows": [
            {"function": r[0], "n_zeros": r[1], "cv": r[2], "category": r[3]}
            for r in verdict_rows
        ],
        "spectral_hp_viable": spectral_hp_viable,
    }

    elapsed = time.time() - t0
    results["elapsed_seconds"] = elapsed
    print(f"\nTotal elapsed: {elapsed:.1f} s")

    # Write JSON
    out_path = Path(__file__).parent / "data" / "spectral_zero_stats.json"
    out_path.parent.mkdir(exist_ok=True, parents=True)
    with out_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nData written to: {out_path}")

    return results


if __name__ == "__main__":
    main()
