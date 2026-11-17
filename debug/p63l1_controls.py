"""Controls for record P63L1 (Paper 63 tree fibre at l = 1).  No evaluation of the l = 1 fibre near t = -2.

C1  s-type reduction: chain of three (t;1,1) recursion (both roots) vs the paper's closed formula; SA's branch
    integrals re-measured (n = 3 eps^-3/2 vs pi/8; n = 5 against the record's numbers; n = 4 control).
C2  the L = 1 leaf potential: direct quadrature (no multipole expansion) at r = 0.5, 2, 7 vs the closed form
    g(r) = (2/pi)[K1 - chi/r^2] and vs the multipole formula, mpmath dps 30; plus an off-axis 3-D check.
C3  the large-r tail of the L = 1 leaf potential (sympy, exact), numeric check at r = 30..60, and the closed-form
    eps^-3/2 amplitude of the branch at t = -2 (derived from the tail only), with its O(eps) correction.
C4  full 3-D Monte Carlo of the l = 1 chain at t = 1/2; root (s) vs leaf rooting; precision self-differences.
Output: debug/data/p63l1_controls.json
"""
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import k0, k1

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "tests"))
import p63l1_engine as E  # noqa: E402
import test_paper63_star_attainment as SAG  # noqa: E402

out = {}
T0 = time.time()

# ================================================================== C1
G = E.Grid()
c1 = {}
for t in (0.5, 1.6):
    r0, r1, cl = E.chain_s_recursion(t, G, 0), E.chain_s_recursion(t, G, 1), E.chain_s_closed(t, G)
    c1[f"chain_s_t{t}"] = {"root_middle": r0, "root_leaf": r1, "closed_8_over_pi3": cl,
                           "rel_root_vs_leaf": (r1 - r0) / r0, "rel_closed_vs_recursion": (cl - r0) / r0}
    print(f"C1 chain_s t={t}: recursion {r0:.12g} leaf-root rel {(r1-r0)/r0:.2e} closed rel {(cl-r0)/r0:.2e}")
sa3 = []
for eps in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
    v = eps**1.5 * SAG._branch(eps, 3)
    sa3.append({"eps": eps, "eps32_I3": v, "ratio_to_pi_over_8": v / (np.pi / 8)})
    print(f"C1 SA n=3 eps={eps:g}: eps^1.5 I3 = {v:.6f}  (pi/8 = {np.pi/8:.6f}, ratio {v/(np.pi/8):.5f})")
c1["SA_n3_branch"] = sa3
sa5 = {f"{e:g}": float(np.sqrt(e) * SAG._branch(e, 5)) for e in (1e-2, 1e-3)}
sa4 = float(1e-3 * SAG._branch(1e-3, 4))
c1["SA_n5_sqrt_eps_I5"] = sa5
c1["SA_n5_record"] = {"1e-2": 1.14957, "0.001": 0.97256}
c1["SA_n4_eps_I4_1e-3"] = sa4
c1["SA_n4_record_1e-3"] = 0.46110
print("C1 SA n=5 re-measured", sa5, "record 1.14957, 0.97256 | n=4 eps I4", sa4, "record 0.46110")
out["C1"] = c1

# ================================================================== C2
mp.mp.dps = 30


def fp_mp(s):
    return -mp.besselk(1, s) / (2 * mp.pi**2)


def V_direct_axis(r):
    """V_p at r zhat by quadrature over (s, mu) of rho_p/|r - r'| (azimuth = 2 pi exactly by symmetry)."""
    r = mp.mpf(r)

    def inner(s):
        return 2 * mp.pi * mp.quad(lambda m: m / mp.sqrt(r * r + s * s - 2 * r * s * m), [-1, 0, 1])
    pts = [0, r / 2, r, r * mp.mpf("1.5"), 2 * r, r + 10, r + 40, mp.inf]
    return mp.quad(lambda s: s * s * fp_mp(s) * inner(s), pts)


def V_multipole(r):
    r = mp.mpf(r)
    a = mp.quad(lambda s: s**3 * fp_mp(s), [0, r])
    b = mp.quad(fp_mp, [r, r + 10, mp.inf])
    return 4 * mp.pi / 3 * (a / r**2 + r * b)


c2 = []
for r in ("0.5", "2", "7"):
    d, m_, c = V_direct_axis(r), V_multipole(r), E.g_closed_mp(r)
    row = {"r": float(r), "direct": mp.nstr(d, 20), "multipole": mp.nstr(m_, 20), "closed": mp.nstr(c, 20),
           "rel_direct_vs_closed": float((d - c) / c), "rel_multipole_vs_closed": float((m_ - c) / c)}
    c2.append(row)
    print(f"C2 r={r}: direct {mp.nstr(d,16)} closed {mp.nstr(c,16)} rel {row['rel_direct_vs_closed']:.1e}"
          f" | multipole rel {row['rel_multipole_vs_closed']:.1e}", flush=True)
# off-axis 3-D quadrature (scipy, all three angles numerically) at r = 1.3, theta = 1.0 rad
rv = 1.3 * np.array([np.sin(1.0), 0.0, np.cos(1.0)])


def f3(phi, mu, s):
    st = math.sqrt(1 - mu * mu)
    d = math.sqrt((s * st * math.cos(phi) - rv[0]) ** 2 + (s * st * math.sin(phi)) ** 2 + (s * mu - rv[2]) ** 2)
    return s * s * (-k1(s) / (2 * np.pi**2)) * mu / d


v3 = 0.0
for lo, hi in ((0, 0.65), (0.65, 1.3), (1.3, 1.95), (1.95, 5), (5, 40)):
    v3 += __import__("scipy.integrate", fromlist=["tplquad"]).tplquad(
        f3, lo, hi, -1, 1, 0, 2 * np.pi, epsabs=1e-11, epsrel=1e-9)[0]
pred3 = float(E.g_closed_mp(1.3)) * math.cos(1.0)
c2.append({"offaxis_r": 1.3, "theta": 1.0, "direct3d": v3, "g_cos": pred3, "rel": (v3 - pred3) / pred3})
print(f"C2 off-axis 3D: {v3:.12g} vs g cos = {pred3:.12g}, rel {(v3-pred3)/pred3:.1e}", flush=True)
out["C2"] = c2

# ================================================================== C3 (symbolic)
x = sp.symbols("x", positive=True)
NT = 7


def hankel(nu, n=NT):
    """K_nu(x) e^x sqrt(2x/pi) = sum a_k x^-k."""
    return [sp.Rational(sp.prod([4 * nu**2 - (2 * j - 1) ** 2 for j in range(1, k + 1)]), sp.factorial(k) * 8**k)
            for k in range(n)]


def tail_int(coeffs, a, n=NT):
    """int_r^oo s^(a - k) e^-s sum_k c_k ds ~ r^a e^-r sum_m d_m r^-m (Gamma(b, r) asymptotic)."""
    d = [sp.Integer(0)] * n
    for k, ck in enumerate(coeffs):
        b = a - k
        for j in range(n - k):
            d[k + j] += ck * sp.ff(b, j)
    return d


def ser(cs):
    return sum(c * x**-k for k, c in enumerate(cs))


s1 = hankel(1)
s0 = hankel(0)
# R = int_r^oo s K1 = sqrt(pi/2) int s^(1/2) e^-s sum a_k s^-k
Rser = tail_int(s1, sp.Rational(1, 2))                 # R ~ sqrt(pi/2) r^(1/2) e^-r Rser
# E e^r = (2/pi)[K1 e^r + R e^r / r^2] = (2/pi) sqrt(pi/2) r^(-1/2) [ s1 + r^-1.5 * r^.5 Rser ... ]
Eser = sp.expand(ser(s1) + x**-1 * ser(Rser))         # E ~ sqrt(2/pi) r^-1/2 e^-r * Eser
Eser = sp.series(Eser, x, sp.oo, NT).removeO()
lead_E = sp.sqrt(2 / sp.pi)
# multipole route (shows the r^(1/2) cancellation): exp part of (4 pi/3)(-1/(2 pi^2)) [ -r^-2 int_r^oo s^3 K1 + r K0 ]
P3 = tail_int(s1, sp.Rational(5, 2))                  # int_r^oo s^3 K1 ~ sqrt(pi/2) r^(5/2) e^-r P3
piece_a = sp.expand(-x**-2 * x**sp.Rational(5, 2) * ser(P3))  # times sqrt(pi/2) e^-r
piece_b = sp.expand(x * x**sp.Rational(-1, 2) * ser(s0))
mult = sp.expand(-sp.Rational(2, 3) / sp.pi * sp.sqrt(sp.pi / 2) * (piece_a + piece_b) * sp.sqrt(x))
mult_lead_half = sp.limit(mult / x, x, sp.oo)          # coefficient of r^(1/2) (times sqrt(x)/x)
mult_ser = sp.series(mult, x, sp.oo, 4).removeO()
print("C3 E tail: E(r) ~ sqrt(2/pi) r^(-1/2) e^-r [", sp.nsimplify(Eser), "]")
coef_mult = sp.simplify(sp.limit(mult, x, sp.oo))
print("C3 multipole route: r^(1/2) e^-r coefficient", sp.simplify(mult_lead_half), "| r^(-1/2) e^-r coefficient", coef_mult)
assert sp.simplify(coef_mult - lead_E) == 0

# numeric check at r = 30..60, mpmath dps 70 through the MULTIPOLE formula (exercises the cancellation)
mp.mp.dps = 70
c3num = []
Ecoef = [sp.nsimplify(sp.expand(Eser).coeff(x, -k)) for k in range(NT)]
for r in (30, 40, 50, 60):
    rr = mp.mpf(r)
    Vm = V_multipole(rr)
    Enum = Vm + 1 / rr**2
    lead = mp.sqrt(2 / mp.pi) * rr ** mp.mpf(-0.5) * mp.exp(-rr)
    asym = lead * sum(mp.mpf(sp.Rational(c).p) / sp.Rational(c).q * rr**-k for k, c in enumerate(Ecoef))
    pa = 4 * mp.pi / 3 * (-1 / (2 * mp.pi**2)) * (-mp.quad(lambda s: s**3 * mp.besselk(1, s), [rr, rr + 20, mp.inf]) / rr**2)
    row = {"r": r, "E_over_lead": float(Enum / lead), "E_over_asym_minus1": float(Enum / asym - 1),
           "multipole_piece_a_over_lead": float(pa / lead)}
    c3num.append(row)
    print(f"C3 r={r}: E/(sqrt(2/pi) r^-1/2 e^-r) = {float(Enum/lead):.12f}; vs {NT}-term series rel"
          f" {float(Enum/asym-1):.1e}; single multipole piece / lead = {float(pa/lead):.3f} (O(r))", flush=True)
mp.mp.dps = 30

# closed-form amplitude at t = -2 (from the tail; no fibre evaluation)
eps = sp.symbols("epsilon", positive=True)


def amplitude(hser_coeffs, pref):
    """H_sing(x) ~ pref * x * e^-2x * sum h_k x^-k.  I0(y) e^-y sqrt(2 pi y) ~ sum b_k y^-k (y = (2-eps) x).
    Singular part of int_1^oo I0((2-eps)x) H_sing dx = sum_k C_k(eps) Gamma(3/2 - k) eps^(k - 3/2)."""
    b = [sp.Rational(sp.prod([(2 * j - 1) ** 2 for j in range(1, k + 1)]), sp.factorial(k) * 8**k) for k in range(4)]
    a2 = 2 - eps
    total = 0
    for k in range(3):
        Ck = 0
        for i in range(k + 1):
            Ck += hser_coeffs[i] * b[k - i] * a2 ** (-(k - i))
        Ck *= pref / sp.sqrt(2 * sp.pi * a2)
        total += Ck * sp.gamma(sp.Rational(3, 2) - k) * eps ** (k - sp.Rational(3, 2))
    s = sp.series(sp.expand(total * eps ** sp.Rational(3, 2)), eps, 0, 2).removeO()
    A = sp.simplify(s.subs(eps, 0))
    kappa = sp.simplify(s.coeff(eps, 1) / A)
    return A, kappa


# l = 1: H_sing = (2/(3pi)) x^2 E^2 = (2/(3pi)) (2/pi) x e^-2x Eser^2
h_l1 = sp.series(sp.expand(Eser**2), x, sp.oo, 4).removeO()
h_l1c = [sp.expand(h_l1).coeff(x, -k) for k in range(3)]
A_l1, kap_l1 = amplitude(h_l1c, sp.Rational(2, 3) / sp.pi * 2 / sp.pi)
# s: H_sing = (2/pi)^3 R^2 = (2/pi)^3 (pi/2) x e^-2x Rser^2
h_s = sp.series(sp.expand(ser(Rser) ** 2), x, sp.oo, 4).removeO()
h_sc = [sp.expand(h_s).coeff(x, -k) for k in range(3)]
A_s, kap_s = amplitude(h_sc, (2 / sp.pi) ** 3 * sp.pi / 2)
print(f"C3 amplitude l=1: eps^1.5 S -> {A_l1} = {float(A_l1):.10f}, (1 + kappa eps), kappa = {kap_l1}")
print(f"C3 amplitude s  : eps^1.5 S -> {A_s} = {float(A_s):.10f}, kappa = {kap_s}; SA-normalised pi/8 check:"
      f" {sp.simplify(A_s / (2 / sp.pi) ** 3)}")
out["C3"] = {
    "E_tail_series": "E(r) = sqrt(2/pi) r^(-1/2) e^(-r) * (" + str(sp.expand(Eser)) + ")",
    "E_lead_coefficient": "sqrt(2/pi)", "E_lead_numeric": float(lead_E),
    "multipole_route_r_half_coefficient": str(sp.simplify(mult_lead_half)),
    "multipole_route_r_minus_half_coefficient": str(coef_mult),
    "numeric": c3num,
    "amplitude_l1": {"exact": str(A_l1), "numeric": float(A_l1), "kappa": str(kap_l1), "kappa_numeric": float(kap_l1)},
    "amplitude_s": {"exact": str(A_s), "numeric": float(A_s), "kappa": str(kap_s), "kappa_numeric": float(kap_s)},
    "J_l1_unprefixed": {"definition": "int_1^oo I0((2-eps)x) x^2 E^2 dx", "exact": str(sp.simplify(A_l1 * 3 * sp.pi / 2)),
                        "numeric": float(A_l1 * 3 * sp.pi / 2)},
}

# ================================================================== C4
c4 = {}
t = 0.5
G2 = E.Grid(400001)
fr, fl, fc = E.fibre_l1_root(t, G), E.fibre_l1_leaf(t, G), E.fibre_l1_closed(t, G)
fr2 = E.fibre_l1_root(t, G2)
c4["grid"] = {"root": fr, "leaf": fl, "closed": fc, "rel_leaf_vs_root": (fl - fr) / fr,
              "rel_closed_vs_root": (fc - fr) / fr, "root_400k": fr2, "rel_grid_doubling": (fr2 - fr) / fr}
print(f"C4 t=1/2 grid: root {fr:.14g} leaf rel {(fl-fr)/fr:.1e} closed rel {(fc-fr)/fr:.1e} 2x grid rel {(fr2-fr)/fr:.1e}")
for tt in (0.5, 1.6):
    a1 = E.fibre_l1_root(tt, G)
    a2 = E.fibre_l1_root(tt, G2)
    c4[f"precision_selfdiff_t{tt}"] = (a2 - a1) / a1

# high-precision closed form at t = 1/2 (mpmath dps 20)
mp.mp.dps = 20
tm = mp.mpf(t)
Fmp = tm**3 * mp.quad(lambda r: mp.besselk(0, tm * r) * 2 / (3 * mp.pi) * r**2 * E.g_closed_mp(r) ** 2,
                      [0, 0.5, 2, 8, 30, 80])
c4["mp_closed_t0.5"] = mp.nstr(Fmp, 18)
c4["rel_grid_vs_mp"] = float((fr - Fmp) / Fmp)
print(f"C4 mp closed t=1/2: {mp.nstr(Fmp,16)}; grid rel {float((fr-Fmp)/Fmp):.1e}", flush=True)

# Monte Carlo of the full 3-D integral, no reduction
rng = np.random.default_rng(11)
Rg = G.R
cdf_root = G.cum0((2 / np.pi) * t**3 * Rg**2 * k0(t * Rg))
cdf_root /= cdf_root[-1]
cdf_leaf = G.cum0(Rg**2 * k1(Rg))
NP = cdf_leaf[-1] / (2 * np.pi**2) * 2 * np.pi        # int |rho_p| d^3r = (int r^2 K1)(2 pi)/(2 pi^2) = 2/pi
cdf_leaf /= cdf_leaf[-1]


def sample_dir(n, signed):
    phi = rng.random(n) * 2 * np.pi
    if signed:
        mu = np.sqrt(rng.random(n)) * np.where(rng.random(n) < 0.5, -1.0, 1.0)   # density |mu|
    else:
        mu = 2 * rng.random(n) - 1
    st = np.sqrt(1 - mu * mu)
    return np.stack([st * np.cos(phi), st * np.sin(phi), mu], 1), mu


acc = []
n = 1_000_000
for _ in range(40):
    r0 = np.interp(rng.random(n), cdf_root, Rg)
    d0, _ = sample_dir(n, False)
    X0 = d0 * r0[:, None]
    r1 = np.interp(rng.random(n), cdf_leaf, Rg)
    d1, m1 = sample_dir(n, True)
    r2 = np.interp(rng.random(n), cdf_leaf, Rg)
    d2, m2 = sample_dir(n, True)
    X1, X2 = d1 * r1[:, None], d2 * r2[:, None]
    # antithetic: |rho_p| is even under z -> -z and rho_p odd, so each leaf contributes
    # sign(mu) [1/|X0 - X| - 1/|X0 - X'|] / 2 with X' the reflected point (unbiased; no reduction used)
    flip = np.array([1.0, 1.0, -1.0])
    a1 = 1 / np.linalg.norm(X0 - X1, axis=1) - 1 / np.linalg.norm(X0 - X1 * flip, axis=1)
    a2 = 1 / np.linalg.norm(X0 - X2, axis=1) - 1 / np.linalg.norm(X0 - X2 * flip, axis=1)
    w = np.sign(m1) * np.sign(m2) * a1 * a2 / 4
    acc.append(NP**2 * w.mean())     # signs: rho_p = -|rho_p| sign(mu); two leaves -> (+)
acc = np.array(acc)
mc, se = acc.mean(), acc.std(ddof=1) / math.sqrt(len(acc))
c4["MC"] = {"t": t, "samples": 40 * n, "mean": mc, "se": se, "recursion": fr, "z": (mc - fr) / se,
            "rel_se": se / abs(fr), "norm_abs_rho_p": NP}
print(f"C4 MC t=1/2: {mc:.6g} +- {se:.2g} vs recursion {fr:.6g}  (z = {(mc-fr)/se:.2f}, rel se {se/fr:.1e})")
out["C4"] = c4
out["runtime_s"] = time.time() - T0
(HERE / "data" / "p63l1_controls.json").write_text(json.dumps(out, indent=1, default=str))
print("wrote debug/data/p63l1_controls.json")
