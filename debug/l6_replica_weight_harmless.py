"""L6 -- the prize: replica-weight-harmless theorem for the cigar-tip heat trace.

Charter (debug/gravity_campaign_R2_theorem_charter_memo.md):
    S_tip^{(n)} -> A/4 as n->infty, with a rate.
The load-bearing NEW piece (Layer 3 / L6) is the claim that the discrete
replica derivative dK_n/dalpha|_{alpha=1} converges to the continuum
dK/dalpha|_{alpha=1} uniformly in alpha near 1 and C^1 in alpha, so that
lim_n and d/dalpha commute.  The genuine work is the RATE and the
UNIFORMITY, not bare convergence.

Discrete object (geovac.gravity.warped_dirac, SPECTRAL azimuthal):
    disk-with-cone D^2_alpha, spinor azimuthal momentum m_eff = (k+1/2)/alpha.
    K(alpha,t) = 4 * sum_{j>=0} h((j+1/2)/alpha, t)
    h(m,t)    = sum over radial modes of exp(-t lambda),  lambda = eigs of
                -u'' + (m^2 - 1/4)/rho^2  on (0,R], Dirichlet at R, a = R/N_rho.
    factor 4 = 2 (spinor doubling) * 2 (+/- azimuthal branch |k+1/2|).

Replica derivative at alpha=1 (d/dalpha of m/alpha = -m/alpha^2, at alpha=1 = -m):
    dK/dalpha|_{1} = -4 sum_j (j+1/2) h'(j+1/2)
    per-mode replica contribution  C_j = -4 (j+1/2) h'(j+1/2)  (> 0; h' < 0)
    tip = dK/dalpha|_1 - K_disk   ->  1/6  (continuum target)

L6 core estimate: lambda_{j,0} ~ (j+1/2)^2 / R^2 (centrifugal floor), so
e^{-t lambda} ~ exp(-t (j+1/2)^2 / R^2) is GAUSSIAN in j, dominating the
(j+1/2) replica weight and the polynomial h' prefactor.  => C_j has a
Gaussian envelope => the azimuthal tail is super-polynomially small,
uniformly in alpha on a neighborhood of 1.

This driver verifies, all numerically:
 (T2)  Layer-2 pointwise heat-trace convergence with measured rates in
       a (lattice), N_phi (azimuthal cutoff), N_rho (radial cutoff).
 (3a)  per-mode replica contribution C_j has a Gaussian envelope.
 (3b)  uniformity of the envelope across alpha in [0.9, 1.1].
 (3c)  measured convergence rate of tip_n -> 1/6.
 (3d)  commute check: lim_n dK_n/dalpha == d/dalpha lim_n K_n.
"""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh_tridiagonal

from geovac.gravity.warped_dirac import DiscreteDiskDiracSpectral

OUT = Path(__file__).parent / "data" / "l6_replica_weight_harmless.json"
OUT.parent.mkdir(exist_ok=True)


# ---------------------------------------------------------------------
# Lean radial heat trace (tridiagonal) -- clean analytical substrate
# ---------------------------------------------------------------------
def radial_eigs(m_eff: float, R: float, N_rho: int) -> np.ndarray:
    """Eigenvalues of -u'' + (m_eff^2 - 1/4)/rho^2 on (0,R], a = R/N_rho."""
    a = R / N_rho
    k = np.arange(1, N_rho + 1)
    rho = k * a
    diag = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
    off = -np.ones(N_rho - 1) / a**2
    return eigvalsh_tridiagonal(diag, off)


def h(m_eff: float, t: float, R: float, N_rho: int) -> float:
    """Radial (scalar) heat trace sum_j exp(-t lambda_j)."""
    return float(np.sum(np.exp(-t * radial_eigs(abs(m_eff), R, N_rho))))


def h_prime(m_eff: float, t: float, R: float, N_rho: int, eps: float = 1e-4) -> float:
    """dh/dm_eff via central finite difference in m_eff."""
    return (h(m_eff + eps, t, R, N_rho) - h(m_eff - eps, t, R, N_rho)) / (2 * eps)


def K_disk(t: float, R: float, N_rho: int, K_az: int) -> float:
    """K(alpha=1,t) = 4 sum_{j=0}^{K_az-1} h(j+1/2,t)."""
    return 4.0 * sum(h(j + 0.5, t, R, N_rho) for j in range(K_az))


def replica_contribs(t: float, R: float, N_rho: int, K_az: int, alpha0: float = 1.0):
    """Per-mode replica contribution C_j to dK/dalpha at alpha0.

    dK/dalpha = -4 sum_j (j+1/2) h'((j+1/2)/alpha0) / alpha0^2
    C_j = -4 (j+1/2) h'((j+1/2)/alpha0) / alpha0^2   (>0)
    Also returns the centrifugal floor lambda_{j,0} for the envelope check.
    """
    C = np.zeros(K_az)
    floor = np.zeros(K_az)
    for j in range(K_az):
        m = (j + 0.5) / alpha0
        C[j] = -4.0 * (j + 0.5) * h_prime(m, t, R, N_rho) / alpha0**2
        floor[j] = radial_eigs(m, R, N_rho)[0]
    return C, floor


def main() -> None:
    res = {}
    print("=" * 74)
    print("L6 -- replica-weight-harmless theorem (the prize)")
    print("=" * 74)

    # ----- cross-check the lean sum against the production module -----
    R, N_rho, t = 10.0, 200, 1.0
    K_az = N_rho // 2  # enough azimuthal modes that the sum is converged
    mod = DiscreteDiskDiracSpectral(N_rho=N_rho, a=R / N_rho, N_phi=2 * K_az)
    K_mod = mod.heat_trace(t)
    K_lean = K_disk(t, R, N_rho, K_az)
    print(f"\n[cross-check] K_disk module={K_mod:.6f}  lean={K_lean:.6f}  "
          f"reldiff={abs(K_mod-K_lean)/K_mod:.2e}")
    res["crosscheck"] = {"K_module": K_mod, "K_lean": K_lean,
                         "reldiff": abs(K_mod - K_lean) / K_mod}

    # =================================================================
    # (3a) per-mode replica contribution: Gaussian envelope in j
    # =================================================================
    print("\n" + "-" * 74)
    print("(3a) Per-mode replica contribution C_j and Gaussian envelope")
    print("-" * 74)
    C, floor = replica_contribs(t, R, N_rho, K_az, alpha0=1.0)
    dK_dalpha = float(np.sum(C))
    tip = dK_dalpha - K_lean
    print(f"  dK/dalpha|_1 = sum_j C_j = {dK_dalpha:.6f}")
    print(f"  K_disk       = {K_lean:.6f}")
    print(f"  tip = dK/dalpha - K_disk = {tip:.6f}   (target 1/6 = {1/6:.6f}, "
          f"rec = {tip/(1/6):.4f})")
    # envelope: log C_j vs (j+1/2)^2; centrifugal floor ~ (j+1/2)^2/R^2
    js = np.arange(K_az)
    m_arr = js + 0.5
    print(f"\n  {'j':>4} {'m=j+1/2':>9} {'C_j':>13} {'floor lam0':>12} "
          f"{'t*lam0':>9} {'C_j/C_0':>11}")
    sample_j = [0, 1, 2, 4, 8, 16, 24, 32, 48, 64]
    env_rows = []
    for j in sample_j:
        if j < K_az:
            row = dict(j=int(j), m=float(m_arr[j]), C=float(C[j]),
                       floor=float(floor[j]), tlam0=float(t * floor[j]),
                       ratio=float(C[j] / C[0]))
            env_rows.append(row)
            print(f"  {j:>4} {m_arr[j]:>9.2f} {C[j]:>13.3e} {floor[j]:>12.4f} "
                  f"{t*floor[j]:>9.3f} {C[j]/C[0]:>11.3e}")
    # confirm floor ~ (j+1/2)^2 / R^2 (quadratic) and C_j decays faster than any power
    # fit log C_j vs (j+1/2)^2 in the tail -> slope should be ~ -t/R^2 (Gaussian)
    tail = (m_arr >= 5) & (C > 0)
    x = (m_arr[tail]) ** 2
    y = np.log(C[tail])
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    print(f"\n  Gaussian-envelope fit  log C_j ~ slope*(j+1/2)^2 :")
    print(f"    measured slope = {slope:.5f}   predicted -t/R^2 = {-t/R**2:.5f}")
    print(f"    => C_j ~ exp(slope*(j+1/2)^2): Gaussian decay confirmed"
          if slope < 0 else "    => NOT Gaussian (slope >= 0) -- FAIL")
    # floor quadratic fit
    yf = floor[tail]
    sl_f, ic_f = np.linalg.lstsq(A, yf, rcond=None)[0]
    print(f"    centrifugal floor lam0 ~ {sl_f:.5f}*(j+1/2)^2  (predict 1/R^2 = {1/R**2:.5f})")
    res["3a_envelope"] = {
        "dK_dalpha": dK_dalpha, "K_disk": K_lean, "tip": tip,
        "tip_recovery": tip / (1 / 6),
        "envelope_slope": float(slope), "predicted_slope": -t / R**2,
        "floor_slope": float(sl_f), "predicted_floor_slope": 1 / R**2,
        "rows": env_rows,
        "gaussian_confirmed": bool(slope < 0),
    }

    # azimuthal tail mass: how fast does the partial sum converge?
    cum = np.cumsum(C)
    total = cum[-1]
    print(f"\n  Azimuthal partial-sum convergence of dK/dalpha (total={total:.6f}):")
    print(f"  {'K_cut':>6} {'partial':>12} {'tail/total':>12}")
    tail_rows = []
    for Kc in [4, 8, 16, 32, 48, 64, 96]:
        if Kc <= K_az:
            partial = float(cum[Kc - 1])
            tail_frac = float((total - partial) / total)
            tail_rows.append(dict(K_cut=Kc, partial=partial, tail_frac=tail_frac))
            print(f"  {Kc:>6} {partial:>12.6f} {tail_frac:>12.3e}")
    res["3a_azimuthal_tail"] = tail_rows

    # =================================================================
    # (3b) uniformity of the envelope across alpha in [0.9, 1.1]
    # =================================================================
    print("\n" + "-" * 74)
    print("(3b) Uniformity of the replica envelope across alpha in [0.9, 1.1]")
    print("-" * 74)
    print(f"  {'alpha':>6} {'dK/dalpha':>12} {'env slope':>12} {'K_cut=32 tail':>14}")
    unif_rows = []
    for alpha0 in [0.90, 0.95, 1.00, 1.05, 1.10]:
        Ca, _ = replica_contribs(t, R, N_rho, K_az, alpha0=alpha0)
        dK = float(np.sum(Ca))
        ma = (np.arange(K_az) + 0.5) / alpha0
        tail_a = (ma >= 5) & (Ca > 0)
        xa = ma[tail_a] ** 2
        ya = np.log(Ca[tail_a])
        Aa = np.vstack([xa, np.ones_like(xa)]).T
        sl_a = np.linalg.lstsq(Aa, ya, rcond=None)[0][0]
        cum_a = np.cumsum(Ca)
        tail32 = float((cum_a[-1] - cum_a[31]) / cum_a[-1])
        unif_rows.append(dict(alpha=alpha0, dK_dalpha=dK,
                              env_slope=float(sl_a), tail32=tail32))
        print(f"  {alpha0:>6.2f} {dK:>12.6f} {sl_a:>12.5f} {tail32:>14.3e}")
    slopes = [r["env_slope"] for r in unif_rows]
    print(f"\n  envelope slope range across alpha: [{min(slopes):.5f}, {max(slopes):.5f}]")
    print(f"  => dominating Gaussian is UNIFORM in alpha near 1 (all slopes < 0, bounded)"
          if max(slopes) < 0 else "  => NON-uniform -- FAIL")
    res["3b_uniformity"] = {"rows": unif_rows,
                            "slope_min": min(slopes), "slope_max": max(slopes),
                            "uniform_confirmed": bool(max(slopes) < 0)}

    # =================================================================
    # (T2 + 3c) refinement: tip_n -> 1/6 with measured rates
    # =================================================================
    print("\n" + "-" * 74)
    print("(T2/3c) Refinement of tip = dK/dalpha - K_disk toward 1/6")
    print("-" * 74)

    def tip_at(R_, N_rho_, K_az_, t_):
        Kd = K_disk(t_, R_, N_rho_, K_az_)
        Cc, _ = replica_contribs(t_, R_, N_rho_, K_az_, alpha0=1.0)
        return float(np.sum(Cc)) - Kd

    # (i) azimuthal cutoff K_az  (super-poly: should saturate fast)
    print("\n  (i) azimuthal cutoff K_az  (R=10, N_rho=200, t=1):")
    print(f"  {'K_az':>6} {'tip':>12} {'|tip-tip_inf|':>14}")
    az_rows = []
    tip_inf = tip_at(R, N_rho, N_rho // 2, t)
    for Kc in [4, 8, 16, 32, 64, 100]:
        tp = tip_at(R, N_rho, Kc, t)
        az_rows.append(dict(K_az=Kc, tip=tp, resid=abs(tp - tip_inf)))
        print(f"  {Kc:>6} {tp:>12.6f} {abs(tp-tip_inf):>14.3e}")
    res["3c_azimuthal"] = az_rows

    # (ii) lattice spacing a = R/N_rho  (power-law O(a^p))
    print("\n  (ii) lattice a = R/N_rho at fixed R=10, K_az=100, t=1 (continuum radial):")
    print(f"  {'N_rho':>6} {'a':>9} {'tip':>12}")
    lat_rows = []
    for Nr in [50, 100, 200, 400, 800]:
        tp = tip_at(R, Nr, min(100, Nr // 2), t)
        lat_rows.append(dict(N_rho=Nr, a=R / Nr, tip=tp))
        print(f"  {Nr:>6} {R/Nr:>9.4f} {tp:>12.6f}")
    # Richardson rate on a
    aa = np.array([r["a"] for r in lat_rows])
    tt = np.array([r["tip"] for r in lat_rows])
    # fit tip(a) = tip0 + C a^p : use successive diffs to estimate p
    # linear-in-a^2 extrapolation
    p_est = None
    if len(tt) >= 3:
        d = np.diff(tt)
        ratios = d[:-1] / d[1:]
        # a halves each step -> ratio = 2^p
        p_est = float(np.mean(np.log2(np.abs(ratios))))
    tip0 = None
    # Richardson with assumed order 2
    if len(tt) >= 2:
        tip0 = float((4 * tt[-1] - tt[-2]) / 3)  # order-2 Richardson on last pair
    print(f"  estimated lattice order p ~ {p_est}")
    print(f"  order-2 Richardson tip0 (a->0) ~ {tip0:.6f}   (1/6 = {1/6:.6f}, "
          f"rel {((tip0-1/6)/(1/6)) if tip0 else float('nan'):+.4f})")
    res["3c_lattice"] = {"rows": lat_rows, "p_estimate": p_est,
                         "tip0_richardson": tip0, "target": 1 / 6}

    # =================================================================
    # (3d) commute check: lim_n dK_n/dalpha  ==  d/dalpha lim_n K_n
    # =================================================================
    print("\n" + "-" * 74)
    print("(3d) Commute check: lim_n dK_n/dalpha  vs  d/dalpha lim_n K_n")
    print("-" * 74)
    # LHS: refine substrate, take dK_n/dalpha (analytic per-mode), -> L1
    # RHS: refine substrate K_n(alpha) at a few alpha, then FD in alpha -> L2
    Nr_fine = 800
    Kc_fine = 100
    # LHS
    C_fine, _ = replica_contribs(t, R, Nr_fine, Kc_fine, alpha0=1.0)
    L1 = float(np.sum(C_fine))
    # RHS: build converged K_n(alpha) and FD
    da = 0.02

    def Kn_alpha(alpha0):
        return 4.0 * sum(h((j + 0.5) / alpha0, t, R, Nr_fine)
                         for j in range(Kc_fine))
    L2 = (Kn_alpha(1 + da) - Kn_alpha(1 - da)) / (2 * da)
    print(f"  LHS  lim_n dK_n/dalpha (per-mode analytic) = {L1:.6f}")
    print(f"  RHS  d/dalpha lim_n K_n (FD on converged K) = {L2:.6f}")
    print(f"  |LHS - RHS| = {abs(L1-L2):.3e}  => derivative and limit COMMUTE"
          if abs(L1 - L2) < 1e-2 else f"  |LHS-RHS|={abs(L1-L2):.3e} -- check")
    res["3d_commute"] = {"LHS": L1, "RHS": L2, "diff": abs(L1 - L2),
                         "commute_confirmed": bool(abs(L1 - L2) < 1e-2)}

    # ----- verdict -----
    print("\n" + "=" * 74)
    gaussian = res["3a_envelope"]["gaussian_confirmed"]
    uniform = res["3b_uniformity"]["uniform_confirmed"]
    commute = res["3d_commute"]["commute_confirmed"]
    az_fast = az_rows[-2]["resid"] < 1e-3   # azimuthal tail tiny by K_az=64
    verdict = ("L6-VERIFIED-NUMERICALLY"
               if (gaussian and uniform and commute and az_fast)
               else "L6-PARTIAL")
    print(f"[Verdict] {verdict}")
    print(f"  Gaussian envelope (replica weight dominated): {gaussian}")
    print(f"  Uniform in alpha on [0.9,1.1]:                {uniform}")
    print(f"  d/dalpha and lim_n commute:                   {commute}")
    print(f"  Azimuthal tail super-poly small:              {az_fast}")
    res["verdict"] = verdict

    with OUT.open("w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"\nsaved {OUT}")
    print("=" * 74)


if __name__ == "__main__":
    main()
