"""B1+B2 -- screened-interval rate lemma + analytic continuum target 1/6.

B1: the discrete radial operator L(m) = -d^2/drho^2 + (m^2-1/4)/rho^2 on (0,R],
Dirichlet at rho=R, apex screened (u ~ rho^{m+1/2} -> 0 at 0). Two questions:
  (i) per-eigenvalue rate: FD lambda_0(m;a) vs EXACT continuum lambda_0(m) =
      (j_{m,1}/R)^2 (Bessel zero). Expect O(a^2) (smooth interior + screened apex).
  (ii) heat-trace / tip rate: the aggregate tip(a) -> tip_0. Measured ratio -> 2
      (a halves) => O(a^1). The gap between (i) O(a^2) and (ii) O(a^1) is the
      apex-edge effect (the grid rho_i = i*a misses the strip [0,a]; the
      heat-trace density there contributes O(a)).

B2: with the CORRECT O(a) order, linear-in-a Richardson should hit
    tip_0 = 1/6, the G4-4c Sommerfeld-Cheeger spinor-cone value
    -(1/12)(1/alpha-alpha) -> dK/dalpha|_1 = 1/6. (The order-2 Richardson in the
    L6 driver gave 0.1643, -1.4%, because it assumed the wrong order.)
"""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh_tridiagonal, eigh_tridiagonal
from scipy.special import jn_zeros

OUT = Path(__file__).parent / "data" / "l6_b1b2_rate_and_target.json"
OUT.parent.mkdir(exist_ok=True)


def radial_eigs(m, R, N_rho, n_low=None):
    a = R / N_rho
    k = np.arange(1, N_rho + 1)
    rho = k * a
    diag = 2.0 / a**2 + (m**2 - 0.25) / rho**2
    off = -np.ones(N_rho - 1) / a**2
    if n_low is None:
        return eigvalsh_tridiagonal(diag, off)
    w = eigh_tridiagonal(diag, off, select="i", select_range=(0, n_low - 1),
                         eigvals_only=True)
    return w


def bessel_lambda0(m, R):
    """Exact continuum lowest eigenvalue: lambda_0 = (j_{m,1}/R)^2.

    m here is the centrifugal index of u'' problem; u = sqrt(rho) J_m(sqrt(lam) rho),
    so the Bessel order is exactly m (the half-integer m_eff).
    """
    # jn_zeros needs integer order; for half-integer use mpmath-free approx via
    # scipy.special.jnp? Use a robust root find of J_m via mpmath if available.
    from scipy.special import jv
    from scipy.optimize import brentq
    # first zero of J_m: bracket near m + 1.86 m^{1/3} + ... ; scan
    lo = max(m, 0.1)
    xs = np.linspace(lo, lo + 10 + 2 * m, 4000)
    vals = jv(m, xs)
    sign = np.sign(vals)
    idx = np.where(np.diff(sign) != 0)[0]
    if len(idx) == 0:
        return None
    x0 = brentq(lambda x: jv(m, x), xs[idx[0]], xs[idx[0] + 1])
    return (x0 / R) ** 2


def main():
    res = {}
    R = 10.0
    print("=" * 72)
    print("B1+B2 -- screened-interval rate + continuum target 1/6")
    print("=" * 72)

    # ---------- B1(i): per-eigenvalue rate vs exact Bessel ----------
    print("\n[B1(i)] per-eigenvalue rate: FD lambda_0(m;a) vs exact (j_{m,1}/R)^2")
    per_eig = {}
    for m in [0.5, 1.5, 2.5]:
        exact = bessel_lambda0(m, R)
        print(f"\n  m={m}  exact lambda_0 = {exact:.8f}")
        print(f"  {'N_rho':>6} {'a':>9} {'lam0_FD':>12} {'err':>12} {'ratio':>8}")
        errs = []
        prev = None
        rows = []
        for Nr in [100, 200, 400, 800, 1600]:
            l0 = float(radial_eigs(m, R, Nr, n_low=1)[0])
            err = abs(l0 - exact)
            ratio = (prev / err) if prev else float("nan")
            print(f"  {Nr:>6} {R/Nr:>9.5f} {l0:>12.8f} {err:>12.3e} {ratio:>8.3f}")
            rows.append(dict(N_rho=Nr, a=R / Nr, lam0=l0, err=err, ratio=ratio))
            errs.append((R / Nr, err))
            prev = err
        # fit order: err ~ C a^q ; log-log slope
        aa = np.log([e[0] for e in errs])
        ee = np.log([e[1] for e in errs])
        q = float(np.polyfit(aa, ee, 1)[0])
        print(f"  => per-eigenvalue order q = {q:.3f}")
        per_eig[str(m)] = dict(exact=exact, rows=rows, order_q=q)
    res["B1_per_eigenvalue"] = per_eig

    # ---------- B1(ii)+B2: tip rate + linear extrapolation to 1/6 ----------
    print("\n[B1(ii)/B2] tip(a) rate and linear-in-a extrapolation to 1/6")

    def h(m, t, R, N_rho):
        return float(np.sum(np.exp(-t * radial_eigs(abs(m), R, N_rho))))

    def hprime(m, t, R, N_rho, eps=1e-4):
        return (h(m + eps, t, R, N_rho) - h(m - eps, t, R, N_rho)) / (2 * eps)

    def tip_at(R, N_rho, K_az, t):
        Kd = 4.0 * sum(h(j + 0.5, t, R, N_rho) for j in range(K_az))
        dK = -4.0 * sum((j + 0.5) * hprime(j + 0.5, t, R, N_rho)
                        for j in range(K_az))
        return dK - Kd

    t = 1.0
    K_az = 90   # azimuthal sum converged to machine precision (L6 (3a))
    print(f"\n  R={R}, t={t}, K_az={K_az}")
    print(f"  {'N_rho':>6} {'a':>9} {'tip':>12} {'1/6 - tip':>12} {'ratio':>8}")
    a_rows = []
    prev_d = None
    Nrs = [50, 100, 200, 400, 800, 1600]
    tips = []
    for Nr in Nrs:
        tp = tip_at(R, Nr, K_az, t)
        d = 1 / 6 - tp
        ratio = (prev_d / d) if prev_d else float("nan")
        print(f"  {Nr:>6} {R/Nr:>9.5f} {tp:>12.7f} {d:>12.3e} {ratio:>8.3f}")
        a_rows.append(dict(N_rho=Nr, a=R / Nr, tip=tp, gap_to_sixth=d, ratio=ratio))
        tips.append((R / Nr, tp))
        prev_d = d

    # order of (1/6 - tip) ~ C a^p
    aa = np.log([x[0] for x in tips])
    dd = np.log([abs(1 / 6 - x[1]) for x in tips])
    p_tip = float(np.polyfit(aa, dd, 1)[0])
    print(f"\n  => tip-gap-to-1/6 order p = {p_tip:.3f}  (O(a) signature: p~1)")

    # linear-in-a Richardson (correct order p=1): tip_0 = tip(a/2)+(tip(a/2)-tip(a))
    a_fine, t_fine = tips[-1]
    a_coarse, t_coarse = tips[-2]
    # general 2-point linear extrapolation to a=0 along a-axis
    tip0_lin = t_fine + (t_fine - t_coarse) * (a_fine - 0) / (a_coarse - a_fine)
    # since a_fine = a_coarse/2: tip0 = 2*t_fine - t_coarse
    tip0_simple = 2 * t_fine - t_coarse
    print(f"\n  linear-in-a extrapolation (correct O(a) order):")
    print(f"    tip_0 = 2*tip(a) - tip(2a) = {tip0_simple:.7f}")
    print(f"    1/6 = {1/6:.7f}   rel err = {(tip0_simple-1/6)/(1/6):+.5f}")
    # also a least-squares linear fit tip = tip0 + c1*a over the 3 finest points
    A = np.array([[x[0], 1.0] for x in tips[-4:]])
    b = np.array([x[1] for x in tips[-4:]])
    c1, tip0_fit = np.linalg.lstsq(A, b, rcond=None)[0]
    print(f"    least-sq linear fit (4 finest): tip_0 = {tip0_fit:.7f}, "
          f"c1 = {c1:.5f}  (apex-edge slope)")
    print(f"    rel err vs 1/6 = {(tip0_fit-1/6)/(1/6):+.5f}")

    res["B1B2_tip"] = {
        "rows": a_rows, "order_p": p_tip,
        "tip0_2pt": float(tip0_simple), "tip0_lstsq": float(tip0_fit),
        "c1_apex_edge_slope": float(c1), "target": 1 / 6,
        "rel_err_2pt": float((tip0_simple - 1 / 6) / (1 / 6)),
        "rel_err_lstsq": float((tip0_fit - 1 / 6) / (1 / 6)),
    }

    # ---------- verdict ----------
    print("\n" + "=" * 72)
    q_mean = np.mean([v["order_q"] for v in per_eig.values()])
    b1_ok = abs(p_tip - 1.0) < 0.25 and q_mean > 1.5   # per-eig O(a^2), tip O(a)
    b2_ok = abs((tip0_fit - 1 / 6) / (1 / 6)) < 0.01
    print(f"[B1] per-eigenvalue order q~{q_mean:.2f} (>=~2, smooth+screened apex); "
          f"tip order p~{p_tip:.2f} (~1, apex-edge O(a)). Gap q-p ~ 1 = edge effect.")
    print(f"     => {'CONFIRMED' if b1_ok else 'CHECK'}: rate is O(a), mechanism = apex-edge")
    print(f"[B2] linear extrapolation tip_0 = {tip0_fit:.6f} vs 1/6 = {1/6:.6f} "
          f"({(tip0_fit-1/6)/(1/6)*100:+.3f}%)")
    print(f"     => {'CONFIRMED' if b2_ok else 'CHECK'}: continuum target is 1/6")
    res["verdict"] = {
        "B1_rate_Oa_apex_edge": bool(b1_ok),
        "B2_target_one_sixth": bool(b2_ok),
        "per_eig_order_mean": float(q_mean),
        "tip_order": float(p_tip),
    }
    with OUT.open("w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"\nsaved {OUT}")
    print("=" * 72)


if __name__ == "__main__":
    main()
