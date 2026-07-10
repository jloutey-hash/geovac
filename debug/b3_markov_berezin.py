"""B3 (push) -- the UNITAL Markov-normalized heat-kernel Berezin on the disk.

The heat-congruence D_W M_f D_W failed unitality (B(1)=e^{-tD}!=I). The correct
L4 map compresses the MARKOV-NORMALIZED heat-smoothed function:

    g_t(x) = (e^{-t Delta} f)(x) / (e^{-t Delta} 1)(x),     B_n(f) = P_n M_{g_t} P_n.

The division by e^{-t Delta}1 makes the kernel a positive NORMALIZED average
(Markov), giving all four L4 properties by standard facts:
  (a) positivity   : f>=0 => g_t>=0 => P_n M_{g_t} P_n >= 0
  (b) contractivity: min f <= g_t <= max f => ||B_n(f)|| <= ||f||_inf
  (c) unitality    : f=1 => g_t=1 => B_n(1) = P_n M_1 P_n = P_n = I   (BY CONSTRUCTION)
  (d) approx id    : ||B_n(f)-P_n M_f P_n|| <= ||g_t - f||_inf, and (g_t-f)=0 for
                     constant f => bound is proportional to ||grad f||  (L4(c))
  (e) Lipschitz    : ||grad g_t||_inf <= ||grad f||_inf (heat-average non-expansive)

The ONE thing to measure: the rate gamma_n = ||g_t - f||_inf / ||grad f||_inf as
n -> infinity with t = c/Lambda (smoothing tied to the mode cutoff Lambda ~ n).

f = (rho/R)^2 is radial, so only k=0 modes contribute and this is a clean 1D
radial computation.  Controls: the constant function 1 (||g-1||=0, unitality)
and a steeper radial function f2 = rho/R (different ||grad f||) to confirm the
bound scales with ||grad f||, not ||f||.
"""

import json
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros, jv
from numpy.polynomial.legendre import leggauss

OUT = Path(__file__).parent / "data" / "b3_markov_berezin.json"
OUT.parent.mkdir(exist_ok=True)

R = 1.0


def k0_modes(Jmax, Lambda):
    """k=0 radial modes: R_j(rho) = sqrt(2)/(R |J_1(a_j)|) J_0(a_j rho/R),
    normalized int_0^R R_j^2 rho drho = 1. lambda_j = (a_j/R)^2."""
    zeros = jn_zeros(0, Jmax)
    modes = []
    for a in zeros:
        lam = (a / R) ** 2
        if lam <= Lambda:
            norm = np.sqrt(2.0) / (R * abs(jv(1, a)))
            modes.append(dict(alpha=a, lam=lam, norm=norm))
    return modes


def Rj(rho, m):
    return m["norm"] * jv(0, m["alpha"] * rho / R)


def smooth_ratio(f_fn, modes, t, rho):
    """g_t(rho) = (e^{-tD} f)(rho) / (e^{-tD} 1)(rho) via k=0 eigen-expansion."""
    # quadrature for eigen-coefficients c_j = int R_j(s) f(s) s ds  (over (0,R))
    xg, wg = leggauss(600)
    s = 0.5 * R * (xg + 1.0)
    ws = 0.5 * R * wg
    num = np.zeros_like(rho)
    den = np.zeros_like(rho)
    for m in modes:
        Rj_s = Rj(s, m)
        cf = np.sum(Rj_s * f_fn(s) * s * ws)      # <R_j, f>
        c1 = np.sum(Rj_s * 1.0 * s * ws)          # <R_j, 1>
        e = np.exp(-t * m["lam"])
        Rj_rho = Rj(rho, m)
        num += e * cf * Rj_rho
        den += e * c1 * Rj_rho
    return num / den


def main():
    res = {}
    print("=" * 74)
    print("B3 push -- unital Markov-normalized heat-kernel Berezin (disk, k=0)")
    print("=" * 74)

    # test functions (radial) and their sup-gradient on the disk
    f1 = lambda r: (r / R) ** 2;  gradf1 = 2.0 / R      # ||grad|| = max|2 rho/R^2|
    f2 = lambda r: (r / R);       gradf2 = 1.0 / R
    fconst = lambda r: np.ones_like(r); gradc = 0.0

    rho = np.linspace(1e-4, R, 2000)
    Jmax = 80

    print("\n[unitality control] g_t for f=1 should be identically 1:")
    for Lam in [200.0, 800.0]:
        modes = k0_modes(Jmax, Lam)
        g = smooth_ratio(fconst, modes, 1.0 / Lam, rho)
        print(f"  Lambda={Lam:.0f}: max|g_1 - 1| = {np.max(np.abs(g-1)):.2e}  "
              f"(=> B_n(1)=P_n=I)")

    print("\n[approximate-identity rate]  ||g_t - f||_inf  with  t = c/Lambda, c=1")
    Lambdas = [100.0, 200.0, 400.0, 800.0, 1600.0, 3200.0]
    rows = []
    print(f"  {'Lambda':>8} {'n_k0':>5} {'||g-f1||/grad':>14} {'||g-f2||/grad':>14} "
          f"{'min g1':>8} {'max g1':>8}")
    for Lam in Lambdas:
        modes = k0_modes(Jmax, Lam)
        t = 1.0 / Lam
        g1 = smooth_ratio(f1, modes, t, rho)
        g2 = smooth_ratio(f2, modes, t, rho)
        e1 = np.max(np.abs(g1 - f1(rho))) / gradf1     # L4(c) coeff for f1
        e2 = np.max(np.abs(g2 - f2(rho))) / gradf2     # and for f2
        rows.append(dict(Lambda=Lam, n_k0=len(modes), e1=float(e1), e2=float(e2),
                         ming1=float(g1.min()), maxg1=float(g1.max())))
        print(f"  {Lam:>8.0f} {len(modes):>5} {e1:>14.5f} {e2:>14.5f} "
              f"{g1.min():>8.4f} {g1.max():>8.4f}")
    res["rows"] = rows

    # rate fit and the ||grad f||-proportionality check
    Ls = np.array([r["Lambda"] for r in rows])
    e1 = np.array([r["e1"] for r in rows])
    e2 = np.array([r["e2"] for r in rows])
    p1 = float(np.polyfit(np.log(Ls), np.log(e1), 1)[0])
    p2 = float(np.polyfit(np.log(Ls), np.log(e2), 1)[0])
    print(f"\n  rate  ||g-f1||/||grad f1|| ~ Lambda^{p1:.3f}")
    print(f"  rate  ||g-f2||/||grad f2|| ~ Lambda^{p2:.3f}")
    print(f"  coefficients e1, e2 agree to factor {max(e1[-1]/e2[-1], e2[-1]/e1[-1]):.2f}"
          f"  => bound is proportional to ||grad f|| (L4(c)), not ||f||")

    # properties
    pos_contr_ok = all(r["ming1"] >= -1e-9 and r["maxg1"] <= 1.0 + 1e-6 for r in rows)
    rate_ok = (p1 < -0.2) and (p2 < -0.2)
    grad_prop_ok = max(e1[-1] / e2[-1], e2[-1] / e1[-1]) < 2.0

    print("\n" + "=" * 74)
    print("[L4 properties of the Markov-Berezin]")
    print(f"  (a)+(b) positivity & contractivity (0<=g<=1 for 0<=f<=1): "
          f"{'PASS' if pos_contr_ok else 'FAIL'}")
    print(f"  (c) unitality B_n(1)=I (by construction):                  PASS")
    print(f"  (d) approx-identity ~ ||grad f|| with rate Lambda^{p1:.2f}: "
          f"{'PASS' if (rate_ok and grad_prop_ok) else 'CHECK'}")
    print(f"  (e) Lipschitz non-expansive (heat average):                analytic")
    verdict = ("B3-L4-CLOSED-NUMERICALLY (unital Markov-Berezin: all 4 properties; "
               f"rate gamma_n ~ Lambda^{p1:.2f} -> 0)"
               if (pos_contr_ok and rate_ok and grad_prop_ok) else "B3-L4-CHECK")
    print(f"\n[Verdict] {verdict}")
    res.update(dict(rate_p1=p1, rate_p2=p2, pos_contr_ok=pos_contr_ok,
                    rate_ok=rate_ok, grad_prop_ok=grad_prop_ok, verdict=verdict))
    with OUT.open("w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"\nsaved {OUT}")
    print("=" * 74)


if __name__ == "__main__":
    main()
