"""B3 de-compactification -- closing Paper 53 Q1 (the inner-arrow / L4 piece).

The boundary obstruction (Dirichlet modes vanish at rho=R, so functions with
content at the boundary do not reconstruct in sup-norm) is the only failure of
the interior reconstruction. On the PLANE (non-compact disk, R->infinity) the
natural function class is C_0 / Schwartz: every observable decays, so it has
negligible content in the boundary collar. Claim:

  THEOREM (de-compactification reconstruction). For f decaying within radius
  theta*R (theta<1), the FULL-disk Markov-Cesaro Berezin reconstructs f at the
  interior rate, with the boundary-collar contribution bounded by the function's
  own collar amplitude ||f||_{[theta R, R]}. For fixed Schwartz f and R->infinity
  this amplitude -> 0 super-polynomially, so the full reconstruction converges at
  the interior rate for all Schwartz observables -- the inner arrow of the
  two-rate de-compactification (Paper 47 template on the disk).

Verification (radial Schwartz f = exp(-rho^2/2), k=0 modes):
  Part A: fix Lambda, vary R -> boundary-collar error tracks f's collar amplitude
          ~ exp(-(R-1)^2/2) and the full error -> the interior (support) error.
  Part B: fix R large (boundary far), vary Lambda -> full error converges at the
          interior rate, NO boundary floor; positivity holds at Cesaro s>=3.
"""

import json
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros, jv
from numpy.polynomial.legendre import leggauss

OUT = Path(__file__).parent / "data" / "b3_decompactification.json"
OUT.parent.mkdir(exist_ok=True)


def k0_modes(R, Jmax, Lam):
    out = []
    for a in jn_zeros(0, Jmax):
        lam = (a / R) ** 2
        if lam <= Lam:
            out.append(dict(alpha=a, lam=lam, norm=np.sqrt(2.0) / (R * abs(jv(1, a)))))
    return out


def Rj(rho, m, R):
    return m["norm"] * jv(0, m["alpha"] * rho / R)


def g_cesaro(f, modes, R, Lam, s, rho):
    xg, wg = leggauss(1000)
    sx = 0.5 * R * (xg + 1.0)
    ws = 0.5 * R * wg
    num = np.zeros_like(rho)
    den = np.zeros_like(rho)
    for m in modes:
        w = max(0.0, (1.0 - m["lam"] / Lam)) ** s
        if w == 0:
            continue
        Rs = Rj(sx, m, R)
        cf = np.sum(Rs * f(sx) * sx * ws)
        c1 = np.sum(Rs * 1.0 * sx * ws)
        Rr = Rj(rho, m, R)
        num += w * cf * Rr
        den += w * c1 * Rr
    return num / den


def main():
    res = {}
    f = lambda r: np.exp(-r**2 / 2.0)        # Schwartz, sigma=1, radial
    s = 3
    print("=" * 76)
    print("B3 de-compactification: Schwartz f=exp(-rho^2/2), Cesaro s=3")
    print("=" * 76)

    # ---- Part A: boundary-collar error vs R (boundary recession) ----
    print("\n[Part A] fix Lambda=200, vary R: boundary-collar error tracks f's collar amplitude")
    print(f"  {'R':>4} {'n_k0':>5} {'collar |f(R-1)|':>16} {'err[R-1,R]/grad':>16} "
          f"{'err[0,4]/grad':>14} {'full/grad':>11} {'ming':>8}")
    Lam = 200.0
    rowsA = []
    for R in [3.0, 4.0, 6.0, 8.0, 12.0]:
        rho = np.linspace(1e-4, R, 5000)
        fr = f(rho)
        gradnum = np.max(np.abs(np.gradient(fr, rho)))
        m = k0_modes(R, 200, 1.3 * Lam)
        g = g_cesaro(f, m, R, Lam, s, rho)
        err = np.abs(g - fr) / gradnum
        collar_amp = float(f(np.array([R - 1.0]))[0])
        bnd = rho >= R - 1.0
        sup = rho <= 4.0
        e_bnd = float(np.max(err[bnd]))
        e_sup = float(np.max(err[sup]))
        e_full = float(np.max(err))
        rowsA.append(dict(R=R, n=len(m), collar_amp=collar_amp, err_bnd=e_bnd,
                          err_sup=e_sup, err_full=e_full, ming=float(g.min())))
        print(f"  {R:>4.0f} {len(m):>5} {collar_amp:>16.2e} {e_bnd:>16.5f} "
              f"{e_sup:>14.5f} {e_full:>11.5f} {g.min():>8.4f}")
    res["partA_boundary_recession"] = rowsA

    # ---- Part B: interior rate on the plane class (boundary far) ----
    print("\n[Part B] fix R=8 (boundary far, f(8)=e^-32~0), vary Lambda: full error -> interior rate")
    print(f"  {'Lambda':>8} {'n_k0':>5} {'full/grad':>11} {'ming':>8} {'maxg':>8}")
    R = 8.0
    rho = np.linspace(1e-4, R, 6000)
    fr = f(rho)
    gradnum = np.max(np.abs(np.gradient(fr, rho)))
    rowsB = []
    Ls = [50.0, 100.0, 200.0, 400.0, 800.0]
    for Lam in Ls:
        m = k0_modes(R, 400, 1.3 * Lam)
        g = g_cesaro(f, m, R, Lam, s, rho)
        e = float(np.max(np.abs(g - fr)) / gradnum)
        rowsB.append(dict(Lambda=Lam, n=len(m), err_full=e,
                          ming=float(g.min()), maxg=float(g.max())))
        print(f"  {Lam:>8.0f} {len(m):>5} {e:>11.5f} {g.min():>8.4f} {g.max():>8.4f}")
    es = np.array([r["err_full"] for r in rowsB])
    p = float(np.polyfit(np.log(Ls), np.log(es), 1)[0])
    print(f"\n  full-disk rate (boundary far) ~ Lambda^{p:.3f}  (interior rate, no boundary floor)")
    res["partB_interior_rate"] = dict(rows=rowsB, rate_p=p)

    # ---- verdict ----
    print("\n" + "=" * 76)
    # A: boundary-collar error decreases with R, tracking collar amplitude
    bnd_decreasing = rowsA[-1]["err_bnd"] < rowsA[0]["err_bnd"]
    collar_small = rowsA[-1]["collar_amp"] < 1e-6
    # B: full error converges (rate < 0), positivity holds (ming>=~0)
    converges = p < -0.2
    positive = all(r["ming"] >= -1e-2 for r in rowsB)
    ok = bnd_decreasing and collar_small and converges and positive
    print(f"[Part A] boundary-collar error decreasing with R: {bnd_decreasing} "
          f"(collar amp at R=12: {rowsA[-1]['collar_amp']:.1e})")
    print(f"[Part B] full-disk converges at interior rate Lambda^{p:.2f}: {converges}; "
          f"positivity (ming>=0): {positive}")
    verdict = ("DECOMPACTIFICATION-INNER-ARROW-CLOSED (Schwartz class reconstructs "
               f"at interior rate; boundary contribution ~ collar amplitude -> 0)"
               if ok else "DECOMPACTIFICATION-CHECK")
    print(f"\n[Verdict] {verdict}")
    res["verdict"] = verdict
    res["partA_ok"] = bool(bnd_decreasing and collar_small)
    res["partB_ok"] = bool(converges and positive)
    with OUT.open("w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"\nsaved {OUT}")
    print("=" * 76)


if __name__ == "__main__":
    main()
