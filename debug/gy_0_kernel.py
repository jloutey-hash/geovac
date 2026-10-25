"""GY stage 0 (docs/predictions/2026-10-23_gaussian_geminal_triangle.md): derivation of the Gaussian loop and of
the kite kernel, each step checked symbolically (sympy) or numerically (mpmath / float64 3-D quadrature).

 (1) Gaussian pair identity: e^{-a|l-x1|^2} e^{-a|l-x2|^2} = e^{-2a|l-c|^2} e^{-a|x1-x2|^2/2}, c = (x1+x2)/2 [SYMBOLIC].
 (2) Loop J(beta, R) = int d^3l e^{-beta l^2}/|l-R|^2 = pi^(3/2) beta^(-1/2) int_0^1 (1-s)^(-1/2) e^{-beta R^2 s} ds
     (from 1/q^2 = int_0^oo dlam e^{-lam q^2} and the Gaussian l-integral) = 2 pi^(3/2) Daw(sqrt(beta) R)/(beta R)
     [SYMBOLIC for the s-integral; numerical against a direct radial-angular quadrature at random (beta, R)].
     Hence loop(x1,x2,x3) = e^{-s12/(8g)} 4 g pi^(3/2) Daw(R/sqrt(2g))/R, R^2 = (s13+s23)/2 - s12/4 (2a = 1/(2g), beta = 2a).
     Checked against direct 3-D float64 quadrature at five random triples (gate 1e-8).
 (3) I(b) = int_0^oo e^{-g r^2/2 - b r} dr closed form, int r e^{..} = (1 - b I)/g [SYMBOLIC], and
     K_G(M, kappa) = [H(b-) - H(b+)]/(2 i kappa), H = I(b)(1 - M b/g), against direct r-quadrature at five (M, kappa)
     including imaginary kappa (complex z) [NUMERICAL, 1e-30].
 (4) The kite prefactor sqrt(pi/8) g^(-3/2): the chain of constants (8/pi) [density Feynman weight] *
     (pi/g)^3 (2 pi g)^(3/2) (2 pi)^-3 [Gaussian pair in position space] * 4 pi / (32 pi) [angular integral / F_3 constant]
     [SYMBOLIC]; the full normalisation is then tested by route A vs route B in gy_1_checks.py.
 (5) Split I(b) = E P(b) - (1/g) P(b) Qh(b) [SYMBOLIC via the ODE and NUMERICAL], and phi_0 = sqrt(pi/8) (pi/8) E at g = 1,
     i.e. phi_0 = (pi^2/32) e^{1/2} erfc(1/sqrt 2) [NUMERICAL against phi_z(0)].
Output: debug/data/gy_0_kernel.json"""
from __future__ import annotations
import json, random, sys
from pathlib import Path
import numpy as np
import mpmath as mp
import sympy as sp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import gy_common as GY

DATA = Path(__file__).parent / "data"
out = {}

# (1) Gaussian pair identity
l1, l2, l3, a1, a2, a3, b1, b2, b3, al = sp.symbols("l1 l2 l3 a1 a2 a3 b1 b2 b3 alpha", real=True)
L = sp.Matrix([l1, l2, l3]); X1 = sp.Matrix([a1, a2, a3]); X2 = sp.Matrix([b1, b2, b3])
lhs = al * ((L - X1).dot(L - X1) + (L - X2).dot(L - X2))
C = (X1 + X2) / 2
rhs = 2 * al * (L - C).dot(L - C) + al * (X1 - X2).dot(X1 - X2) / 2
out["pair_identity"] = bool(sp.simplify(lhs - rhs) == 0)
print("(1) Gaussian pair identity:", out["pair_identity"], flush=True)

# (2) loop: s-integral symbolically
s, beta, R = sp.symbols("s beta R", positive=True)
yv0, cc = sp.symbols("y0 c", positive=True)
# substitute 1 - s = y^2: int_0^1 (1-s)^(-1/2) e^{-c s} ds = 2 e^{-c} int_0^1 e^{c y^2} dy
Js = 2 * sp.exp(-cc) * sp.integrate(sp.exp(cc * yv0**2), (yv0, 0, 1))
closed = 2 * sp.exp(-cc) * sp.sqrt(sp.pi) * sp.erfi(sp.sqrt(cc)) / (2 * sp.sqrt(cc))   # = 2 Daw(sqrt c)/sqrt c, c = beta R^2
diff = sp.simplify(Js - closed)
out["loop_s_integral_symbolic"] = bool(diff == 0)
print("(2) s-integral = 2 Daw(sqrt(beta) R)/(sqrt(beta) R):", out["loop_s_integral_symbolic"], "  (sympy:", Js, ")", flush=True)
# numerical J(beta, R) against a direct radial-angular quadrature: J = (2 pi/R) int_0^oo l e^{-beta l^2} ln|(l+R)/(l-R)| dl
mp.mp.dps = 30
rows = []
rng = random.Random(20261023)
for _ in range(4):
    be, Rv = mp.mpf(rng.uniform(0.2, 3)), mp.mpf(rng.uniform(0.1, 4))
    direct = 2 * mp.pi / Rv * mp.quad(lambda l: l * mp.exp(-be * l * l) * mp.log(abs((l + Rv) / (l - Rv))), [0, Rv / 2, Rv, 2 * Rv, mp.inf])
    cl = 2 * mp.pi ** mp.mpf(1.5) * GY.dawson(mp.sqrt(be) * Rv) / (be * Rv)
    rows.append({"beta": float(be), "R": float(Rv), "rel": float(abs(direct / cl - 1))})
out["loop_J_vs_radial_quadrature"] = rows
print("(2) J(beta,R) closed vs radial quadrature, worst rel:", max(r["rel"] for r in rows), flush=True)


def direct3d(x1, x2, g, Nang=160, Nr=240):
    """int d^3l g(l-x1) g(l-x2)/|l|^2, g(q) = exp(-q^2/(4g)), float64 Gauss rules (copy of ga3_1_checks.direct)."""
    x1, x2 = np.array(x1), np.array(x2)
    uu, wu = np.polynomial.legendre.leggauss(Nang)
    ph = np.linspace(0, 2 * np.pi, Nang, endpoint=False); wph = 2 * np.pi / Nang
    st = np.sqrt(1 - uu**2)
    n = np.stack([np.outer(st, np.cos(ph)), np.outer(st, np.sin(ph)), np.outer(uu, np.ones_like(ph))], -1)
    d1, d2 = n @ x1, n @ x2
    Rr = max(np.linalg.norm(x1), np.linalg.norm(x2)) + 6 * np.sqrt(g)
    sr, wr = np.polynomial.legendre.leggauss(Nr)
    gfun = lambda q2: np.exp(-q2 / (4 * g))

    def f_at(ri):
        q1 = ri**2 - 2 * ri * d1 + x1 @ x1
        q2 = ri**2 - 2 * ri * d2 + x2 @ x2
        return gfun(q1) * gfun(q2)
    tot = 0.0
    for lo, hi in [(0, Rr / 3), (Rr / 3, 2 * Rr / 3), (2 * Rr / 3, Rr)]:
        r = (hi - lo) / 2 * sr + (hi + lo) / 2; w = (hi - lo) / 2 * wr
        for ri, wi in zip(r, w):
            tot += wi * np.sum(wu[:, None] * f_at(ri)) * wph
    return tot


rows = []
g = 1.0
for k in range(5):
    x1 = [rng.uniform(-2, 2) for _ in range(3)]; x2 = [rng.uniform(-2, 2) for _ in range(3)]
    s12 = sum((mp.mpf(a) - b) ** 2 for a, b in zip(x1, x2)); s13 = sum(mp.mpf(a) ** 2 for a in x1); s23 = sum(mp.mpf(b) ** 2 for b in x2)
    cl = GY.loop_gauss(g, s12, s13, s23)
    dd = direct3d(x1, x2, g)
    rel = abs(dd / float(cl) - 1)
    rows.append({"x1": x1, "x2": x2, "direct": dd, "closed": float(cl), "rel": rel})
    print(f"(2) loop triple {k}: direct {dd:.15e} closed {float(cl):.15e} rel {rel:.2e}", flush=True)
out["loop_vs_direct3d"] = rows
out["loop_vs_direct3d_worst"] = max(r["rel"] for r in rows)

# (3) I(b), the r-moment identity, K_G against direct quadrature
r, b, gs = sp.symbols("r b gamma", positive=True)
I_closed = sp.sqrt(sp.pi / (2 * gs)) * sp.exp(b**2 / (2 * gs)) * sp.erfc(b / sp.sqrt(2 * gs))
# I(b) is fixed by its ODE I' = (b/g) I - 1/g (differentiate under the integral and integrate by parts) and I(0) = sqrt(pi/(2g))
I0 = sp.integrate(sp.exp(-gs * r**2 / 2), (r, 0, sp.oo))
out["I_closed_symbolic"] = bool(sp.simplify(sp.diff(I_closed, b) - (b / gs * I_closed - 1 / gs)) == 0 and sp.simplify(I_closed.subs(b, 0) - I0) == 0)
# int r e^{-g r^2/2 - b r} dr = -I'(b) = (1 - b I)/g by the ODE; checked numerically too
J1_num = mp.quad(lambda rr: rr * mp.exp(-mp.mpf("1.3") * rr * rr / 2 - mp.mpf("0.7") * rr), [0, mp.inf])
I_num = sp.N(I_closed.subs({gs: sp.Rational(13, 10), b: sp.Rational(7, 10)}), 30)
out["r_moment_identity"] = bool(abs(J1_num - (1 - mp.mpf("0.7") * mp.mpf(str(I_num))) / mp.mpf("1.3")) < mp.mpf(10) ** -25)
print("(3) I(b) closed form:", out["I_closed_symbolic"], " int r e^{..} = (1 - b I)/g:", out["r_moment_identity"], flush=True)
ctx.dps = 40
rows = []
for k in range(5):
    Mv = mp.mpf(rng.uniform(0.8, 2.5)); kv = mp.mpf(rng.uniform(0.1, 2.0))
    if k >= 3:
        kappa = mp.mpc(0, kv)         # imaginary kappa: kappa^2 < 0 (complex z)
        kap2 = -kv * kv
    else:
        kappa = kv
        kap2 = kv * kv
    direct = GY.KG_direct(Mv, kappa, 1, dps=40)
    cl = GY.KG_parts(acb(mp.nstr(Mv, 45)), acb(mp.nstr(kap2, 45)), acb(0), arb(1), "full")
    clm = GY.to_mp if hasattr(GY, "to_mp") else None
    cl_mp = mp.mpc(cl.real.mid().str(45, radius=False), cl.imag.mid().str(45, radius=False))
    rel = float(abs(cl_mp / direct - 1))
    rows.append({"M": float(Mv), "kappa2": float(kap2), "rel": rel})
    print(f"(3) K_G M={float(Mv):.4f} kappa^2={float(kap2):+.4f}: closed {mp.nstr(cl_mp, 20)} direct {mp.nstr(direct, 20)} rel {rel:.1e}", flush=True)
# kappa = 0 limit
direct0 = GY.KG_direct(mp.mpf("1.3"), 0, 1, dps=40)
cl0 = GY.KG_parts(acb("1.3"), 0, acb(0), arb(1), "full")
rel0 = float(abs(mp.mpf(cl0.real.mid().str(45, radius=False)) / direct0 - 1))
rows.append({"M": 1.3, "kappa2": 0.0, "rel": rel0})
print(f"(3) K_G at kappa = 0: rel {rel0:.1e}", flush=True)
out["KG_vs_direct"] = rows
out["KG_vs_direct_worst"] = max(r["rel"] for r in rows)

# (4) the prefactor chain
pref = sp.Rational(8) / sp.pi * (sp.pi / gs) ** 3 * (2 * sp.pi * gs) ** sp.Rational(3, 2) * (2 * sp.pi) ** -3 * 4 * sp.pi / (32 * sp.pi)
out["prefactor_symbolic"] = str(sp.simplify(pref))
out["prefactor_equals_sqrt_pi_over_8_g32"] = bool(sp.simplify(pref - sp.sqrt(sp.pi / 8) * gs ** sp.Rational(-3, 2)) == 0)
print("(4) prefactor:", sp.simplify(pref), out["prefactor_equals_sqrt_pi_over_8_g32"], flush=True)

# (5) split: ODE check I' = (b/g) I - 1/g, P' = (b/g) P, (P Qh)' = (b/g) P Qh + 1; numerical identity at random complex b
Pb = sp.exp((b**2 - 1) / (2 * gs))
yv = sp.Symbol("y", positive=True)
Qhb = sp.integrate(sp.exp(-(yv**2 - 1) / (2 * gs)), (yv, 1, b))
R_sym = -Pb * Qhb / gs
out["R_ode"] = bool(sp.simplify(sp.diff(R_sym, b) - (b / gs * R_sym - 1 / gs)) == 0)
out["I_ode"] = bool(sp.simplify(sp.diff(I_closed, b) - (b / gs * I_closed - 1 / gs)) == 0)
out["R_at_1"] = bool(sp.simplify(R_sym.subs(b, 1)) == 0)
print("(5) ODEs: I", out["I_ode"], " R = -(1/g) P Qh", out["R_ode"], " R(1) = 0", out["R_at_1"], flush=True)
E = GY.E_const(arb(1))
rows = []
for _ in range(4):
    bb = acb(rng.uniform(0.3, 2.5), rng.uniform(-1.5, 1.5))
    lhs = GY.I_fun(bb, arb(1))
    rhs = E * GY.P_fun(bb, arb(1)) - GY.P_fun(bb, arb(1)) * GY.Qh_fun(bb, arb(1))
    rows.append(float(abs(lhs - rhs) / abs(lhs)))
out["split_identity_rel"] = rows
# phi_0
phi0 = GY.phi_z(acb(0), arb(1), 400, "full")
phi0_hand = arb.pi() ** 2 / 32 * arb(1) / 2 * 0 + arb.pi() ** 2 / 32 * (arb(1) / 2).exp() * (1 / arb(2).sqrt()).erfc()
rel = float(abs(phi0.real - phi0_hand) / phi0_hand)
out["phi0_vs_hand"] = {"phi0": phi0.real.mid().str(30, radius=False), "hand": phi0_hand.mid().str(30, radius=False), "rel": rel}
out["E"] = E.mid().str(40, radius=False)
print("(5) phi_0 =", phi0.real.mid().str(25, radius=False), " hand (pi^2/32) e^{1/2} erfc(1/sqrt2) =", phi0_hand.mid().str(25, radius=False), "rel", rel, flush=True)
# E phi_a + phi_b = phi at a complex z
zz = acb("0.9", "1.4")
pa, pb_, pf = (GY.phi_z(zz, arb(1), 400, w) for w in ("a", "b", "full"))
rel = float(abs(E * pa + pb_ - pf) / abs(pf))
out["E_phia_plus_phib_equals_phi"] = rel
print("(5) E phi_a + phi_b = phi at z = 0.9+1.4i: rel", rel, flush=True)
worst = max(out["loop_vs_direct3d_worst"], out["KG_vs_direct_worst"], max(rows), out["phi0_vs_hand"]["rel"], rel)
out["worst"] = worst
out["gate_1e-8"] = "PASS" if worst < 1e-8 and all(out[k] for k in ("pair_identity", "loop_s_integral_symbolic", "I_closed_symbolic", "r_moment_identity", "R_ode", "I_ode", "R_at_1", "prefactor_equals_sqrt_pi_over_8_g32")) else "FAIL"
print("WORST", worst, "GATE", out["gate_1e-8"], flush=True)
(DATA / "gy_0_kernel.json").write_text(json.dumps(out, indent=1, default=str))
