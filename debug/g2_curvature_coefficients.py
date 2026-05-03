"""
Phase 1: Algebraic extraction of curvature expansion coefficients for g-2 on S^3.

Strategy: compute_vertex_3pt_single_level() returns exact sympy Rationals for
each (n_ext, n_int) pair. The per-level magnetic contribution

    f(n_ext, n_int) = [B_up(n_ext, n_int) - B_dn(n_ext, n_int)] / V_tree_magnetic(n_ext)

is an exact Rational for each n_ext. Since the Dirac eigenvalue is
lambda = n_ext + 3/2 = (2*n_ext + 3)/2, and the propagator denominators
are polynomials in lambda, f must be a rational function of lambda.

Phase 1: Compute f(n_ext, n_int) for n_int=0..4, n_ext=1..10.
Phase 2: Identify the rational function R(lambda) for each n_int via
         interpolation on the exact data points.
Phase 3: Laurent-expand at lambda -> infinity to get per-n_int c_k contributions.
"""

import json
import sys
import time

sys.path.insert(0, '.')

from sympy import (Rational, S, cancel, Symbol, series, oo, Poly,
                   factorial, simplify, nsimplify, together, apart)
from sympy import pi as sym_pi
from geovac.qed_anomalous_moment import (
    compute_vertex_3pt_single_level,
    tree_level_probe_magnetic,
)

lam = Symbol('lam')

N_EXT_MAX = 10
N_INT_MAX = 4


def compute_per_level_data():
    """Compute exact Rational f(n_ext, n_int) for all pairs."""
    j_ext = Rational(1, 2)
    data = {}

    for n_ext in range(1, N_EXT_MAX + 1):
        lam_val = Rational(2 * n_ext + 3, 2)
        v_tree = tree_level_probe_magnetic(n_ext, j_ext, q_probe=1)
        if v_tree == 0:
            print(f"  n_ext={n_ext}: V_tree_magnetic = 0 (skip)")
            continue

        data[n_ext] = {'lambda': lam_val, 'V_tree': v_tree, 'levels': {}}
        print(f"\nn_ext={n_ext}  lambda={lam_val}  V_tree={v_tree}")

        for n_int in range(0, N_INT_MAX + 1):
            t0 = time.time()
            b_up = compute_vertex_3pt_single_level(
                n_ext, j_ext, Rational(1, 2), n_int, q_probe=1)
            b_dn = compute_vertex_3pt_single_level(
                n_ext, j_ext, Rational(-1, 2), n_int, q_probe=1)
            b_mag = b_up - b_dn
            f_lev = b_mag / v_tree
            dt = time.time() - t0
            data[n_ext]['levels'][n_int] = f_lev

            if f_lev != 0:
                print(f"  n_int={n_int}: f = {f_lev}  ({float(f_lev):.6e})  [{dt:.1f}s]")
            else:
                print(f"  n_int={n_int}: f = 0  [{dt:.1f}s]")

    return data


def interpolate_rational_function(points, max_total_deg=18):
    """Find minimal-degree rational function through exact (lambda, f) pairs.

    Uses the Cauchy interpolation method: for P(x)/Q(x) = f_i at x_i,
    multiply through to get P(x_i) - f_i * Q(x_i) = 0.
    With Q monic of degree d_q, this is a homogeneous linear system
    in the coefficients of P (degree d_p) and Q (degree d_q - 1 free coeffs).
    """
    from sympy import Matrix, zeros

    n_pts = len(points)
    if n_pts == 0:
        return None

    for total_deg in range(1, max_total_deg + 1):
        for d_p in range(0, total_deg + 1):
            d_q = total_deg - d_p
            if d_q < 0:
                continue
            n_params = (d_p + 1) + d_q  # P has d_p+1 coeffs, Q has d_q free (leading=1)
            if n_params > n_pts:
                continue

            # Build system: for each (x_i, f_i):
            # a_0 + a_1*x + ... + a_{d_p}*x^{d_p} - f_i*(b_0 + b_1*x + ... + b_{d_q-1}*x^{d_q-1}) = f_i * x^{d_q}
            rows = []
            rhs = []
            for x_i, f_i in points:
                row = []
                for k in range(d_p + 1):
                    row.append(x_i ** k)
                for k in range(d_q):
                    row.append(-f_i * x_i ** k)
                rows.append(row)
                rhs.append(f_i * x_i ** d_q)

            A = Matrix(rows[:n_params])
            b_vec = Matrix(rhs[:n_params])

            try:
                sol = A.solve(b_vec)
            except Exception:
                continue

            # Build P and Q
            P = S.Zero
            for k in range(d_p + 1):
                P += sol[k] * lam ** k
            Q = lam ** d_q
            for k in range(d_q):
                Q += sol[d_p + 1 + k] * lam ** k

            # Verify against ALL points
            ok = True
            for x_i, f_i in points:
                pred = cancel(P.subs(lam, x_i) / Q.subs(lam, x_i) - f_i)
                if pred != 0:
                    ok = False
                    break

            if ok:
                R = cancel(P / Q)
                return R, d_p, d_q

    return None


def laurent_at_infinity(R_func, n_terms=8):
    """Laurent-expand R(lambda) at lambda -> infinity.

    Returns dict {power: coeff} where R ~ sum coeff * lam^(-power).
    """
    u = Symbol('u')
    R_u = R_func.subs(lam, 1 / u)
    R_u = cancel(R_u)
    taylor = series(R_u, u, 0, n=n_terms)

    coeffs = {}
    for k in range(n_terms):
        c = taylor.coeff(u, k)
        if c != 0:
            coeffs[k] = c
    return coeffs


def main():
    print("=" * 70)
    print("  PHASE 1: ALGEBRAIC g-2 CURVATURE COEFFICIENTS")
    print("=" * 70)
    print(f"  n_ext = 1..{N_EXT_MAX}, n_int = 0..{N_INT_MAX}")
    print(f"  All in exact sympy Rational arithmetic")
    print()

    t_start = time.time()

    # Step 1: exact per-level data
    print("--- Step 1: Per-level f(n_ext, n_int) ---")
    data = compute_per_level_data()

    t1 = time.time()
    print(f"\nStep 1 took {t1 - t_start:.1f}s")

    # Step 2: rational function identification per n_int
    print("\n--- Step 2: Rational function identification ---")
    rational_functions = {}

    for n_int in range(0, N_INT_MAX + 1):
        pts = []
        for n_ext in sorted(data.keys()):
            f_val = data[n_ext]['levels'].get(n_int, S.Zero)
            if f_val != 0:
                pts.append((data[n_ext]['lambda'], f_val))

        print(f"\nn_int={n_int}: {len(pts)} nonzero points")
        if not pts:
            rational_functions[n_int] = {'status': 'zero'}
            continue

        result = interpolate_rational_function(pts)
        if result is None:
            print(f"  Could not identify rational function")
            rational_functions[n_int] = {
                'status': 'unidentified',
                'points': [(str(x), str(f)) for x, f in pts],
            }
            continue

        R_func, d_p, d_q = result
        print(f"  R(lambda) = {R_func}  [deg_num={d_p}, deg_den={d_q}]")

        # Step 3: Laurent expand
        coeffs = laurent_at_infinity(R_func, n_terms=10)
        print(f"  Laurent at infinity:")
        for power in sorted(coeffs.keys()):
            c = coeffs[power]
            print(f"    lam^(-{power}): {c}  ({float(c):.10e})")

        rational_functions[n_int] = {
            'status': 'found',
            'R_func': str(R_func),
            'deg_num': d_p,
            'deg_den': d_q,
            'laurent': {str(k): str(v) for k, v in coeffs.items()},
        }

    # Step 4: sum Laurent coefficients over n_int
    print("\n" + "=" * 70)
    print("  PARTIAL c_k SUMS (n_int = 0..{})".format(N_INT_MAX))
    print("=" * 70)

    # The curvature expansion of F2/Schwinger:
    #   F2 = sum_{n_int} f(lambda, n_int)
    #   F2/Schwinger = sum_{n_int} f(lambda, n_int) / Schwinger
    #
    # Each f(lambda, n_int) ~ sum_k a_k(n_int) * lam^(-k)
    #
    # The FULL sum: F2(lambda) = sum_{n_int=0}^{inf} R_{n_int}(lambda)
    # has Laurent expansion F2 ~ sum_k [sum_{n_int} a_k(n_int)] * lam^{-k}
    #
    # F2/Schwinger = 1 + c1/lam^2 + c2/lam^4 + ...
    # means the leading term sum_{n_int} a_0(n_int) = Schwinger
    # and c_k = sum_{n_int} a_{2k}(n_int) / Schwinger
    #
    # But we don't need alpha here. We just sum the per-n_int Laurent coefficients
    # at each power and report them. The Schwinger normalization is handled
    # by the fact that F2/Schwinger - 1 = [sum f(lam,n_int) - Schwinger] / Schwinger.

    total = {}
    for n_int in range(0, N_INT_MAX + 1):
        rf = rational_functions[n_int]
        if rf['status'] != 'found':
            continue
        for k_str, v_str in rf['laurent'].items():
            k = int(k_str)
            v = Rational(v_str)
            if k not in total:
                total[k] = S.Zero
            total[k] += v

    print("\nSummed Laurent coefficients:")
    for k in sorted(total.keys()):
        v = total[k]
        print(f"  lam^(-{k}):  {v}  =  {float(v):.12e}")

    # The key question: does this partial sum already reveal the structure?
    # For the curvature expansion c_k, we need the ratio:
    # c_k = a_{2k} / a_0 where a_0 is the flat-space Schwinger limit.
    #
    # a_0 = sum_{n_int} [lam^0 coefficient of R_{n_int}(lam)]
    # This is the Schwinger limit contribution from n_int=0..N_INT_MAX only.

    if 0 in total and total[0] != 0:
        a0_partial = total[0]
        print(f"\n  a_0 (partial, n_int=0..{N_INT_MAX}) = {a0_partial}  ({float(a0_partial):.12e})")
        print(f"  Schwinger = alpha/(2*pi) = 1.1614e-3")
        print(f"  Ratio a_0_partial / Schwinger ~ {float(a0_partial) / 1.1614e-3:.6f}")

        for k in [2, 4, 6, 8]:
            if k in total:
                ck = total[k] / a0_partial
                print(f"\n  c_{k//2} (partial) = a_{k} / a_0 = {ck}  ({float(ck):.10e})")
                if k == 2:
                    print(f"    Parker-Toms: c_1 = R/12 = 1/2")
                    print(f"    Partial c_1 = {float(ck):.8f}")
                if k == 4:
                    print(f"    Known c_2 = 19/100 - 41*pi^2/25200 = 0.17394...")
                    print(f"    Partial c_2 = {float(ck):.8f}")

    # Also: examine the n_int-dependence of each Laurent coefficient
    # to predict the infinite sum
    print("\n" + "=" * 70)
    print("  PER-n_int STRUCTURE OF LAURENT COEFFICIENTS")
    print("=" * 70)

    for power in sorted(total.keys()):
        print(f"\n  Coefficient of lam^(-{power}):")
        vals = []
        for n_int in range(0, N_INT_MAX + 1):
            rf = rational_functions[n_int]
            if rf['status'] != 'found':
                continue
            laurent = rf['laurent']
            v = Rational(laurent.get(str(power), '0'))
            vals.append((n_int, v))
            lambda_int = Rational(2 * n_int + 3, 2)
            print(f"    n_int={n_int} (lam_int={lambda_int}):  {v}  ({float(v):.8e})")

        # Check if these form a Dirichlet-like pattern: a_k(n_int) ~ C / lam_int^p
        if len(vals) >= 3:
            nonzero = [(ni, float(v)) for ni, v in vals if v != 0]
            if len(nonzero) >= 2:
                import math
                lam_ints = [Rational(2 * ni + 3, 2) for ni, _ in nonzero]
                ratios = []
                for i in range(1, len(nonzero)):
                    r = nonzero[i][1] / nonzero[i-1][1]
                    lr = float(lam_ints[i]) / float(lam_ints[i-1])
                    if abs(r) > 1e-30 and abs(lr) > 0:
                        p = -math.log(abs(r)) / math.log(lr)
                        ratios.append(p)
                if ratios:
                    print(f"    Estimated power-law exponent: p ~ {sum(ratios)/len(ratios):.2f}")

    t_total = time.time()
    print(f"\nTotal time: {t_total - t_start:.1f}s")

    # Save
    output = {
        'n_ext_max': N_EXT_MAX,
        'n_int_max': N_INT_MAX,
        'rational_functions': {},
        'partial_laurent_sums': {},
        'per_level': {},
    }

    for n_int, rf in rational_functions.items():
        output['rational_functions'][str(n_int)] = rf

    for k, v in total.items():
        output['partial_laurent_sums'][str(k)] = {
            'exact': str(v),
            'float': float(v),
        }

    for n_ext in sorted(data.keys()):
        d = data[n_ext]
        output['per_level'][str(n_ext)] = {
            'lambda': str(d['lambda']),
            'V_tree': str(d['V_tree']),
            'levels': {str(ni): str(fv) for ni, fv in d['levels'].items()},
        }

    with open('debug/data/g2_curvature_coefficients.json', 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to debug/data/g2_curvature_coefficients.json")


if __name__ == '__main__':
    main()
