"""FD1-FD3 step B: normalization, momentum-space control, and L S by quadrature
(record docs/predictions/2026-10-21_four_density_star.md, commit 9a33fbef).

Normalization (derivation).  Fourier convention rho(r) = (2 pi)^-3 int rho~(k) e^{ik.r} d^3k,
rho~(0) = 1, rho~_v = (c_v k^2 + 1)^(-3/2), c = 1/a^2.  Since int_0^inf k sin(kr) (1+k^2)^(-3/2) dk
= r K0(r), the density is rho_a(r) = a^3 K0(a r)/(2 pi^2).  Its potential is
phi(r) = int rho(r')/|r - r'| d^3r' = (2/(pi r)) int_0^inf sin(kr)/k (1+k^2)^(-3/2) dk
       = (2/pi) chi(r)/r,  since d/dr of the integral is int cos(kr)(1+k^2)^(-3/2) dk = r K1(r).
A star with central density of scale t and n_l unit leaves is
  int d^3r rho_c(r) phi(r)^n_l = (4 pi/(2 pi^2)) (2/pi)^n_l t^3 int chi^n_l K0(t x) x^(2 - n_l) dx
                               = (2/pi)^(n_l + 1) t^3 int chi^n_l K0(t x) x^(2 - n_l) dx.
n_l = 2 (the chain, n = 3): 8/pi^3, as registered.  n_l = 3 (the star): N = 16/pi^4.

Momentum control (independent of chi and K0).  With 1/r = (2 pi)^-3 int (4 pi/k^2) e^{ik.r} d^3k,
  S = (4 pi)^3 (2 pi)^-9 I = I/(8 pi^6),
  I = int d^3k1 d^3k2 d^3k3 f(k1) f(k2) f(k3) g(|k1 + k2 + k3|) = 4 pi int_0^inf P^2 H(P) M(P) dP,
  f(k) = (1+k^2)^(-3/2)/k^2,  g(k) = (1 + k^2/t^2)^(-3/2),  H = f*f,  M = f*g  (3-D convolutions),
  (a*b)(P) = (2 pi/P) int_0^inf k a(k) [B(P+k) - B(|P-k|)] dk,  B' = q b(q):
  B_f(q) = 1/w + (1/2) log((w-1)/(w+1)), w = sqrt(1+q^2);   B_g(q) = -t^2/sqrt(1 + q^2/t^2).
Same code with two leaves (S2 = 4 pi int P^2 g H dP/(4 pi^4)) must give the chain, B.C_mom/C_pos.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import density_loop_B_ode as B  # noqa: E402
import wg_closure_A as W  # noqa: E402

u = B.u
OUT = HERE / "data" / "four_density_star_B.json"
_CHI = {}


def chi3(xx):
    key = (mp.mp.dps, xx)
    if key not in _CHI:
        _CHI[key] = B.chi(xx) ** 3
    return _CHI[key]


def S_pos(tv):
    tv = mp.mpf(tv)
    f = lambda xx: chi3(xx) / xx * mp.besselk(0, tv * xx) if xx != 0 else mp.mpf(0)
    return 16 / mp.pi**4 * tv**3 * mp.quad(f, [0, mp.mpf(1) / 8, 1, 4, 16, 64, mp.inf])


def _conv(kf, Bfun, P):
    def h(k):
        if k == 0:
            return mp.mpf(0)
        return kf(k) * (Bfun(P + k) - Bfun(abs(P - k)))
    return 2 * mp.pi / P * mp.quad(h, [0, P / 2, P, 2 * P + 1, mp.inf])


def S_mom(tv, leaves=3):
    tv = mp.mpf(tv)
    kf = lambda k: (1 + k * k) ** mp.mpf(-1.5) / k

    def Bf(q):
        if q == 0:
            return -mp.inf
        w = mp.sqrt(1 + q * q)
        return 1 / w + mp.log(q / (w + 1))      # (1/2) log((w-1)/(w+1)), w - 1 = q^2/(w + 1)
    Bg = lambda q: -tv**2 / mp.sqrt(1 + q * q / tv**2)
    g = lambda P: (1 + P * P / tv**2) ** mp.mpf(-1.5)
    if leaves == 3:
        integrand = lambda P: P * P * _conv(kf, Bf, P) * _conv(kf, Bg, P) if P != 0 else mp.mpf(0)
        return 4 * mp.pi * mp.quad(integrand, [0, mp.mpf(1) / 4, 1, 3, 10, mp.inf]) / (8 * mp.pi**6)
    integrand = lambda P: P * P * _conv(kf, Bf, P) * g(P) if P != 0 else mp.mpf(0)
    return 4 * mp.pi * mp.quad(integrand, [0, mp.mpf(1) / 4, 1, 3, 10, mp.inf]) / (4 * mp.pi**4)


def LS_num(Q, tv):
    return W.LF_num(Q, chi3, 0, -1, 3, tv, pref=16 / mp.pi**4)


def main() -> None:
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    QS = [W.rebind(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    QS = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in QS]
    # --- control: N and momentum
    if "momentum" not in rec:
        mp.mp.dps = 20
        t0 = time.time()
        rows = {}
        c_pos, c_mom2 = B.C_pos("0.7"), S_mom("0.7", leaves=2)
        rows["chain_t0.7"] = {"C_pos": mp.nstr(c_pos, 18), "C_mom_my_code": mp.nstr(c_mom2, 18),
                              "rel": mp.nstr(abs(c_mom2 / c_pos - 1), 3)}
        print("chain (my momentum code, 2 leaves) vs B.C_pos:", rows["chain_t0.7"], f"{time.time()-t0:.0f}s", flush=True)
        for tv in ["0.7", "1.6"]:
            a, b = S_pos(tv), S_mom(tv)
            rows[f"star_t{tv}"] = {"S_pos": mp.nstr(a, 18), "S_mom": mp.nstr(b, 18), "rel": mp.nstr(abs(a / b - 1), 3)}
            print(f"star t={tv}:", rows[f"star_t{tv}"], f"{time.time()-t0:.0f}s", flush=True)
        rec["momentum"] = rows
        OUT.write_text(json.dumps(rec, indent=1))
    # --- L S by quadrature at 40 digits, and L 1
    mp.mp.dps = 40
    rec.setdefault("LS", {})
    for tv in ["0.5", "1.3", "2.5", "3.5"]:
        if tv in rec["LS"]:
            continue
        t0 = time.time()
        val = LS_num(QS, tv)
        sv = S_pos(tv)
        one = sum(mp.mpf(tv) ** (2 * i) * mp.mpf(sp.Rational(sp.Poly(q, u).eval(0)).p) / sp.Rational(sp.Poly(q, u).eval(0)).q
                  for i, q in enumerate(QS))
        # largest single term t^(2i) [u^j] Q_i theta^j S, for scale
        rec["LS"][tv] = {"S": mp.nstr(sv, 30), "LS": mp.nstr(val, 25), "L1": mp.nstr(one, 25),
                         "LS_over_L1": mp.nstr(val / one, 25)}
        print(f"t={tv}: S = {mp.nstr(sv, 25)}  L S = {mp.nstr(val, 20)}  L 1 = {mp.nstr(one, 20)}"
              f"  ratio {mp.nstr(val / one, 20)}  ({time.time()-t0:.0f}s)", flush=True)
        OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
