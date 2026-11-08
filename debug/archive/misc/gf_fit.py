"""GF1-GF3 (record 2026-10-25_gaussian_fit_geminal.md): n-term Gaussian least-squares fits of exp(-r).
  W1: min int_0^oo r^2 [e^{-r} - sum c_i e^{-alpha_i r^2}]^2 dr
  W2: the same with the extra weight e^{-r^2}.
Variable projection: c = G^{-1} b solved exactly (mpmath) at each alpha; the exponents by BFGS in log alpha with the
analytic gradient dF/dalpha_k = 2 c_k [<f, r^2 g_k> - sum_j c_j <g_j, r^2 g_k>]; several even-tempered starts and a warm
start from the (n-1)-term fit.  Control: W1, n = 3 against STO-3G (Hehre-Stewart-Pople 1969, zeta = 1).
usage: python gf_fit.py <nmax>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

mp.mp.dps = 60
DATA = Path(__file__).parent / "data"


def Ik(a, kmax):
    """I_k(a) = int_0^oo r^k e^{-a r^2 - r} dr, k = 0..kmax (upward recurrence; stable for these k)."""
    a = mp.mpf(a)
    I = [mp.sqrt(mp.pi / (4 * a)) * mp.exp(1 / (4 * a)) * mp.erfc(1 / (2 * mp.sqrt(a)))]
    I.append((1 - I[0]) / (2 * a))
    for k in range(1, kmax):
        I.append((k * I[k - 1] - I[k]) / (2 * a))
    return I


def parts(al, shift):
    n = len(al)
    G = mp.matrix(n, n); G4 = mp.matrix(n, n)
    b = mp.matrix(n, 1); b4 = mp.matrix(n, 1)
    for i in range(n):
        I = Ik(al[i] + shift, 4)
        b[i] = I[2]; b4[i] = I[4]
        for j in range(n):
            s = al[i] + al[j] + shift
            G[i, j] = mp.sqrt(mp.pi) / (4 * s ** mp.mpf(1.5))
            G4[i, j] = 3 * mp.sqrt(mp.pi) / (8 * s ** mp.mpf(2.5))
    return G, G4, b, b4


def ff(shift):
    return mp.quad(lambda r: r * r * mp.exp(-2 * r - shift * r * r), [0, 1, 5, mp.inf])


def F_grad(la, shift, f2):
    al = [mp.e ** mp.mpf(x) for x in la]
    G, G4, b, b4 = parts(al, shift)
    try:
        c = mp.lu_solve(G, b)
    except ZeroDivisionError:
        # coalescing exponents: a bad point for the line search, not a fit
        return float(f2), np.zeros(len(al)), None, al
    F = f2 - (b.T * c)[0]
    g = []
    for k in range(len(al)):
        s = b4[k] - sum(c[j] * G4[j, k] for j in range(len(al)))
        g.append(float(2 * c[k] * s * al[k]))  # d/d log alpha
    return float(F), np.array(g), c, al


def fit(n, shift, f2, starts):
    best = None
    for la0 in starts:
        res = minimize(lambda x: F_grad(x, shift, f2)[0], np.array(la0, float), jac=lambda x: F_grad(x, shift, f2)[1],
                       method="BFGS", options={"gtol": 1e-30, "maxiter": 4000})
        F, g, c, al = F_grad(res.x, shift, f2)
        if c is None:
            continue
        if best is None or F < best[0]:
            best = (F, res.x, c, al, float(np.max(np.abs(g))))
    return best


def main():
    nmax = int(sys.argv[1])
    out = {}
    for wname, shift in (("W1", 0), ("W2", 1)):
        f2 = ff(shift)
        out[wname] = {}
        prev = None
        for n in range(1, nmax + 1):
            starts = [np.log(a0) + np.log(beta) * np.arange(n) for a0 in (0.05, 0.1, 0.3) for beta in (2.5, 4.0, 7.0)]
            if prev is not None:
                starts.append(np.sort(np.append(prev, prev.max() + 2.0)))
                starts.append(np.sort(np.append(prev, prev.min() - 1.5)))
            F, la, c, al, gmax = fit(n, shift, f2, starts)
            prev = np.array(la)
            order = sorted(range(n), key=lambda i: al[i])
            rel = mp.sqrt(F / f2) if F > 0 else mp.mpf(0)
            rr = [mp.mpf(k) / 100 for k in range(0, 1001)]
            fitv = lambda r: sum(c[i] * mp.exp(-al[i] * r * r) for i in range(n))
            maxabs = max(abs(mp.exp(-r) - fitv(r)) for r in rr)
            out[wname][str(n)] = {"alpha": [mp.nstr(al[i], 30) for i in order], "c": [mp.nstr(c[i], 30) for i in order],
                                  "rel_norm_error": mp.nstr(rel, 6), "max_abs_error_0_10": mp.nstr(maxabs, 6),
                                  "grad_max": gmax}
            print(wname, n, "rel norm err", mp.nstr(rel, 4), "max|err| on [0,10]", mp.nstr(maxabs, 4),
                  "alpha", [mp.nstr(al[i], 6) for i in order], "grad", f"{gmax:.1e}", flush=True)
    (DATA / "gf_fits.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
