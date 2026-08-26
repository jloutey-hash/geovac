"""aha Track 1 / DERIVATION -- where N^1.85 and N^1.97 come from.

Key structural observation
--------------------------
With the sine basis  e_a(chi) = sqrt(2/pi) sin(a chi)  on L^2[0,pi], the SW cross-center
block is EXACTLY the finite section of a MULTIPLICATION operator:

    C_ab(s) = (2/pi) int_0^pi sin(a chi) sin(b chi) W_s(chi) dchi = <e_a, M_{W_s} e_b>,
    W_s(chi) = j_0(s cot(chi/2)),        s = k R.

So  C_n = P_n M_{W_s} P_n  and  sigma_max(n) = ||P_n M_W P_n|| -> ||M_W|| = sup |W| = 1,
attained at chi = pi (cot -> 0, j_0(0) = 1).  Everything follows:

  * sigma_max -> 1 for ANY R -- the near-linear dependence is the symbol touching 1.
  * The RATE is a band-limited concentration problem.  Near chi = pi (delta = pi - chi),
        1 - W_s = 1 - j_0(s tan(delta/2)) = (s^2/24) delta^2 + O(delta^4),
    and  1 - sigma_max(n) = min { <f,(1-W)f> : f in span{sin(a delta)}_{a<=n}, ||f||=1 }.
    Rescaling delta = x/n turns this into  min <x^2> for an odd band-limited profile,
    whose Fourier-dual is the Dirichlet-Dirichlet problem -d^2/dt^2 on [0,1]
    (ct(0) = 0 because a >= 1, ct(1) = 0 because a finite second moment forbids a hard
    band edge), with lowest eigenvalue pi^2.  Hence

        min <delta^2> = pi^2 / n^2      =>      1 - sigma_max  ~=  c_sym * pi^2 / n^2,

    with  c_sym = coefficient of delta^2 in (1 - symbol) at its maximum  ( = s^2/24 for SW ),
    i.e.  1 - sigma_max ~= pi^2 s^2 / (24 n^2) = pi^2 s^2 / (6 N^2),   cond ~= 12 N^2/(pi^2 s^2).

  ==> the exponent is EXACTLY 2, not 1.85 / 1.97; the fitted values are pre-asymptotic
      slopes over the paper's window N = 4..20 (resp. 6..36).

  * Same machinery for water: the A1 whitened coupling is (in the limit) multiplication by
        g(chi) = sqrt2 * W_{R_OH}(chi) / sqrt(1 + W_{R_HH}(chi)),   sup g = g(pi) = 1,
    so water inherits the SAME exponent 2 with a different curvature c_sym.

  * It also explains the sigma-SPECTRUM: by Szego/finite-section, the sigma_k are a
    discretization of the RANGE of |symbol|, so #{sigma_k > tau} grows LINEARLY in n.
    The O<->H coupling therefore cannot be low-rank -- which is why a Woodbury/low-rank
    whitening buys a constant factor but no change of exponent.

Run:  python debug/aha_t1_asymptote.py
"""
from __future__ import annotations

import sys, os, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from debug.aha_t1_core import sw_block_trapz, loglog_fit

R_OH = 1.809
R_HH = 2.0 * R_OH * np.sin(np.deg2rad(104.5) / 2.0)

# ------------------------------------------------------- high-n sine-basis machinery
def _grid(M: int):
    chi = np.linspace(1e-9, np.pi, M)
    w = np.full(M, (chi[1] - chi[0]))
    w[0] = w[-1] = (chi[1] - chi[0]) / 2.0          # trapezoid weights
    return chi, w


def _sinc(x):
    out = np.ones_like(x)
    nz = x != 0.0
    out[nz] = np.sin(x[nz]) / x[nz]
    return out


def sym_block(symbol, nmax: int, M: int = 300001, chunk: int = 60000) -> np.ndarray:
    """<e_a, M_f e_b> for a = 1..nmax, f = symbol(chi); chunked so memory stays bounded."""
    chi, w = _grid(M)
    F = symbol(chi)
    out = np.zeros((nmax, nmax))
    a = np.arange(1, nmax + 1)[:, None]
    for i0 in range(0, M, chunk):
        c = chi[None, i0:i0 + chunk]
        Sm = np.sin(a * c)
        out += (Sm * (w[i0:i0 + chunk] * F[i0:i0 + chunk])) @ Sm.T
    return (2.0 / np.pi) * out


def W_sw(s):
    return lambda chi: _sinc(s / np.tan(chi / 2.0))


def g_water(chi):
    t = 1.0 / np.tan(chi / 2.0)
    return np.sqrt(2.0) * _sinc(R_OH * t) / np.sqrt(1.0 + _sinc(R_HH * t))


def main() -> None:
    print("=" * 100)
    print("DERIVATION -- the cross-center block is a finite section of a MULTIPLICATION operator")
    print("=" * 100)

    # -------------------------------------------------- 1. validate the symbol picture
    print("\n[1] Validate: sym_block(W_s) reproduces the Paper-60 SW cross-block")
    for s in (1.4, 2.0, 4.0):
        Cref = sw_block_trapz(s, 8, "S")
        Cnew = sym_block(W_sw(s), 8, M=300001)
        print(f"    s={s:>4.1f}:  max|C_ref - C_symbol| = {np.abs(Cref-Cnew).max():.2e}"
              f"   (max|C| = {np.abs(Cref).max():.4f})")

    # ----------------------------------- 2. the pure band-limited concentration constant
    print("\n[2] The band-limited concentration constant:  n^2 * min <delta^2>  ->  pi^2 ?")
    print("    (Dirichlet-Dirichlet -d^2/dt^2 on [0,1];  pi^2 = 9.8696)")
    print(f"  {'n':>5} | {'n^2 min<delta^2>':>17} {'/ pi^2':>9}")
    print("  " + "-" * 38)
    for n in (10, 20, 40, 80, 160):
        D = sym_block(lambda c: (np.pi - c) ** 2, n, M=300001)
        lm = float(np.linalg.eigvalsh(D).min())
        print(f"  {n:>5} | {n*n*lm:>17.5f} {n*n*lm/np.pi**2:>9.5f}")

    # ---------------------------------------- 3. SW: local exponent -> 2, prefactor check
    print("\n[3] SW two-center metric: local exponent of  1 - sigma_max  vs n")
    print("    predicted:  1 - sigma_max = (pi^2 s^2 / 24) n^-2 ,  exponent -> 2 exactly")
    ns = [5, 10, 20, 40, 80, 160]
    res = {}
    for s in (1.4, 2.0, 4.0):
        vals = []
        for n in ns:
            C = sym_block(W_sw(s), n, M=300001)
            lam = float(np.linalg.eigvalsh(C).max())
            vals.append(1.0 - lam)
        res[s] = vals
        pred = [np.pi ** 2 * s ** 2 / (24.0 * n * n) for n in ns]
        loc = [np.log(vals[i] / vals[i + 1]) / np.log(ns[i + 1] / ns[i])
               for i in range(len(ns) - 1)]
        print(f"\n    s = kR = {s}")
        print(f"      {'n':>5} " + " ".join(f"{n:>11}" for n in ns))
        print(f"      {'1-sig':>5} " + " ".join(f"{v:>11.3e}" for v in vals))
        print(f"      {'pred':>5} " + " ".join(f"{p:>11.3e}" for p in pred))
        print(f"      {'ratio':>5} " + " ".join(f"{v/p:>11.4f}" for v, p in zip(vals, pred)))
        print(f"      local exponents (successive n-pairs): "
              + "  ".join(f"{x:.4f}" for x in loc))

    print("\n    Paper-60 window fit for comparison (N = 2n = 4..20, i.e. n = 2..10):")
    for s in (1.4, 2.0, 4.0):
        Nw = [2 * n for n in range(2, 11)]
        vv = []
        for n in range(2, 11):
            C = sym_block(W_sw(s), n, M=300001)
            vv.append(1.0 - float(np.linalg.eigvalsh(C).max()))
        a, r2, mr = loglog_fit(Nw, vv)
        cd = [(2 - v) / v for v in vv]
        b, r2b, mrb = loglog_fit(Nw, cd)
        print(f"      s={s:>4.1f}:  1-sigma ~ N^{a:+.3f} (R2={r2:.4f})    "
              f"cond ~ N^{b:+.3f} (R2={r2b:.4f})   <-- the 'fit-only' exponent")

    # -------------------------------------------- 3c. the (N, R) collapse, asymptotic
    print("\n[3c] COLLAPSE.  The derived law makes 1 - sigma_max a function of the SINGLE")
    print("     variable x = n/s = n/(kR) alone:   (1 - sigma_max) * x^2  ->  pi^2/24 ="
          f" {np.pi**2/24:.6f}")
    print(f"  {'n \\ s':>7} | " + "  ".join(f"{s:>9.1f}" for s in (1.4, 2.0, 3.0, 4.0, 6.0, 10.0)))
    print("  " + "-" * 76)
    for n in (10, 20, 40, 80, 160):
        row = []
        for s in (1.4, 2.0, 3.0, 4.0, 6.0, 10.0):
            C = sym_block(W_sw(s), n, M=300001)
            v = 1.0 - float(np.linalg.eigvalsh(C).max())
            row.append(v * (n / s) ** 2)
        print(f"  {n:>7} | " + "  ".join(f"{x:>9.5f}" for x in row))
    print(f"  {'limit':>7} | " + "  ".join(f"{np.pi**2/24:>9.5f}" for _ in range(6)))
    print("     => the collapse IS clean; the 'no clean collapse' seen on the paper's")
    print("        N = 4..20 window is purely a pre-asymptotic (n ~ s) artifact.")

    # ------------------------------------------------------------- 4. water symbol
    print("\n[4] Water A1: the whitened O<->H+ coupling is multiplication by")
    print("    g(chi) = sqrt2 j0(R_OH cot(chi/2)) / sqrt(1 + j0(R_HH cot(chi/2)))")
    d = np.linspace(1e-7, 0.35, 9)
    gv = g_water(np.pi - d)
    cfit = np.polyfit(d ** 2, 1 - gv, 1)
    c_sym = float(cfit[0])
    print(f"    g(pi) = {g_water(np.array([np.pi-1e-9]))[0]:.12f}   (sup of the symbol = 1)")
    print(f"    curvature fit  1 - g = c_sym * delta^2 + ...   ->  c_sym = {c_sym:.6f}")
    print(f"    analytic:  c_sym = R_OH^2/24 - R_HH^2/96 = "
          f"{R_OH**2/24 - R_HH**2/96:.6f}")
    print(f"    predicted law:  1 - sigma_max = c_sym * pi^2 / n^2 = "
          f"{c_sym*np.pi**2:.5f} / n^2")
    print(f"\n  {'n':>5} | {'1-sigma_max':>13} {'predicted':>12} {'ratio':>8}"
          f" {'local exp':>10}")
    print("  " + "-" * 56)
    prev = None
    wvals = []
    for n in (6, 12, 24, 48, 96):
        C = sym_block(g_water, n, M=300001)
        lam = float(np.linalg.eigvalsh(C).max())
        v = 1.0 - lam
        wvals.append((n, v))
        p = c_sym * np.pi ** 2 / n ** 2
        le = "" if prev is None else f"{np.log(prev[1]/v)/np.log(n/prev[0]):>10.4f}"
        print(f"  {n:>5} | {v:>13.6e} {p:>12.6e} {v/p:>8.4f} {le:>10}")
        prev = (n, v)

    print("\n    Paper-60 water window (N = 3*nmax = 6..36):")
    Nw, vv = [], []
    for nmax in range(2, 13):
        C = sym_block(g_water, nmax, M=300001)
        v = 1.0 - float(np.linalg.eigvalsh(C).max())
        Nw.append(3 * nmax); vv.append(v)
    a, r2, mr = loglog_fit(Nw, vv)
    print(f"      1 - sigma_max ~ N^{a:+.3f} (R2={r2:.4f})   <-- pre-asymptotic; true exponent = -2")

    # ------------------------------------------ 5. why the coupling is NOT low rank
    print("\n[5] Why the O<->H coupling cannot be low-rank (Szego / finite-section):")
    print("    the sigma_k discretize the RANGE of the symbol, so #{sigma_k > tau} ~ n * mu(tau)/pi")
    print("    with mu(tau) = |{chi : |symbol(chi)| > tau}|.")
    chi = np.linspace(1e-9, np.pi, 400001)
    gg = np.abs(g_water(chi))
    print(f"  {'tau':>6} {'mu/pi (measure frac)':>21} | " +
          " ".join(f"n={n}" for n in (6, 12, 24, 48)))
    print("  " + "-" * 74)
    for tau in (0.9, 0.7, 0.5, 0.3, 0.1):
        frac = float((gg > tau).mean())
        counts = []
        for n in (6, 12, 24, 48):
            C = sym_block(g_water, n, M=300001)
            sv = np.linalg.svd(C, compute_uv=False)
            counts.append(f"{int((sv>tau).sum())}/{n} ({(sv>tau).mean():.2f})")
        print(f"  {tau:>6.2f} {frac:>21.4f} | " + "  ".join(counts))
    print("\n    => the fraction of sigma_k above any threshold converges to the symbol")
    print("       measure: the coupling rank grows LINEARLY with n.  No low-rank shortcut.")

    # -------------------------------- 5b. the WHOLE top of the spectrum, and Woodbury
    print("\n[5b] The k-th canonical correlation obeys the SAME law with the k-th Dirichlet")
    print("     eigenvalue:   1 - sigma_k  ~=  c_sym * (k pi)^2 / n^2.")
    print("     => removing the top r directions exactly (Woodbury) leaves")
    print("            kappa_res(r) ~= 2 n^2 / (c_sym (r+1)^2 pi^2),")
    print("        i.e. a CONSTANT factor (r+1)^2 -- the N^2 exponent is untouched at fixed r.")
    for s, lab in ((2.0, "SW s=2.0"), (None, "water")):
        cs = (s ** 2 / 24.0) if s is not None else c_sym
        print(f"\n     {lab}   (c_sym = {cs:.6f})")
        print(f"       {'n':>5} | " + "  ".join(f"k={k}" for k in range(1, 6))
              + "      (each cell: measured / predicted)")
        for n in (40, 80, 160):
            C = sym_block(W_sw(s) if s is not None else g_water, n, M=300001)
            ev = np.sort(np.linalg.eigvalsh(C))[::-1]
            cells = []
            for k in range(1, 6):
                meas = 1.0 - ev[k - 1]
                pred = cs * (k * np.pi) ** 2 / n ** 2
                cells.append(f"{meas/pred:.4f}")
            print(f"       {n:>5} | " + "  ".join(f"{c:>7}" for c in cells))

    # ------------------------------------------------ 6. the gerade lever, DERIVED
    print("\n[6] The H2+ gerade lever is the SAME symbol statement, at the other end of the range.")
    print("    For equivalent centers C is symmetric and the symmetry exchanges the two")
    print("    subspaces, so the eigenvectors of S are (u, +/-u):  gerade block = I + C,")
    print("    ungerade block = I - C.  Hence")
    print("        cond(gerade)   -> (1 + sup W) / (1 + inf W) = 2 / (1 + min_x j0(x))")
    print("        cond(ungerade) -> (1 - inf W) / (1 - sup W) = (1 - min j0) / (1 - sigma_max) ~ N^2")
    print("    min_x j0(x) = min_x sin x / x = -0.21723362821...  (at tan x = x, x = 4.49341)")
    kg_inf = 2.0 / (1.0 - 0.2172336282112216)
    print(f"    => cond(gerade) -> 2/(1 - 0.2172336) = {kg_inf:.6f}, INDEPENDENT of R and of N.")
    print(f"\n  {'n':>5} | " + "  ".join(f"cond(I+C) s={s}" for s in (1.4, 2.0, 4.0)))
    print("  " + "-" * 62)
    for n in (10, 20, 40, 80, 160):
        row = []
        for s in (1.4, 2.0, 4.0):
            C = sym_block(W_sw(s), n, M=300001)
            row.append(float(np.linalg.cond(np.eye(n) + C)))
        print(f"  {n:>5} | " + "  ".join(f"{x:>14.6f}" for x in row))
    print(f"    limit  = {kg_inf:.6f} for every R  <-- the 'flat kappa ~ 2' of Paper 60 sec:resource,")
    print("             now an exact constant of the sinc symbol rather than an observation.")
    print("\n    Water has no such block because O and H+ are NOT exchanged by any symmetry:")
    print("    the intra-blocks differ (I vs I+Q), so the 1 +/- sigma eigenvectors are")
    print("    (alpha u, +/- beta v) with alpha != beta and BOTH stay inside A1.  That is the")
    print("    exact structural content of 'the lever is an equivalent-center property' (v4.95.0).")

    os.makedirs("debug/data", exist_ok=True)
    with open("debug/data/aha_t1_asymptote.json", "w") as fh:
        json.dump(dict(sw={str(k): v for k, v in res.items()},
                       water=[[n, v] for n, v in wvals], c_sym_water=c_sym,
                       kappa_gerade_limit=kg_inf), fh, indent=1)
    print("\n  wrote debug/data/aha_t1_asymptote.json")
    print("DONE.")


if __name__ == "__main__":
    main()
