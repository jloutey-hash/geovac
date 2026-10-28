"""WG derivation (record docs/predictions/2026-10-28_water_gaussian_contact.md, phase 1).

Derivation only.  Computes NO scored F_l.  What it does:
  1. STO-3G / STO-6G 1s expansions (zeta = 1) from the Psi4 H rows (exponents / 1.24^2); normalization and
     overlap with the Slater 1s (transcription check).
  2. Single-Gaussian case (i): the closed form F_l = Pref (2l+1) int_0^1 e^{-T0 u^2} i_l(T1 u^2) du, its large-l
     expansion R_l = 1 + c1/l + c2/l^2 (symbolic), checked on a TOY (T0, T1) that is not a scored case, against
     a direct 40-digit Legendre projection and the k-space (Weber) representation.
  3. Case (ii) tail-filter model: F_l(n) / S_l ~ [chi_X^(n)(D)/chi_X^S(D)]^2 * Leg_l[V^(n)] / Leg_l[V^S], where V is
     the Coulomb integral of the H orbital with itself as a charge distribution, at separation w = D sqrt(2(1-x)).
     Model-Slater is compared with the stored WA5 Slater S_l (known data, disclosed); the model's predicted
     crossover l* (ratio = 1/2) is printed.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
OUT = Path(__file__).parent / "data"
D = mp.mpf('1.809')
Z2H = mp.mpf('1.24') ** 2
# Psi4 share/psi4/basis/sto-3g.gbs and sto-6g.gbs, element H (zeta_H = 1.24); exponents scaled back to zeta = 1.
STO = {
    3: ([mp.mpf('3.42525091'), mp.mpf('0.62391373'), mp.mpf('0.16885540')],
        [mp.mpf('0.15432897'), mp.mpf('0.53532814'), mp.mpf('0.44463454')]),
    6: ([mp.mpf('0.3552322122E+02'), mp.mpf('0.6513143725E+01'), mp.mpf('0.1822142904E+01'),
         mp.mpf('0.6259552659E+00'), mp.mpf('0.2430767471E+00'), mp.mpf('0.1001124280E+00')],
        [mp.mpf('0.9163596281E-02'), mp.mpf('0.4936149294E-01'), mp.mpf('0.1685383049E+00'),
         mp.mpf('0.3705627997E+00'), mp.mpf('0.4164915298E+00'), mp.mpf('0.1303340841E+00')]),
}
STO = {n: ([a / Z2H for a in al], d) for n, (al, d) in STO.items()}
log = []


def say(s):
    print(s)
    log.append(str(s))


def contraction(n, zeta):
    """Normalized STO-nG of exponent zeta: list of (exponent, coefficient of e^{-a r^2})."""
    al, d = STO[n]
    terms = [(zeta ** 2 * a, di * (2 * zeta ** 2 * a / mp.pi) ** mp.mpf(0.75)) for a, di in zip(al, d)]
    S = sum(ci * cj * (mp.pi / (ai + aj)) ** 1.5 for ai, ci in terms for aj, cj in terms)
    return [(a, c / mp.sqrt(S)) for a, c in terms], S


# ---- 1. transcription check
for n in (3, 6):
    terms, S = contraction(n, 1)
    ov = sum(c * mp.quad(lambda r: 4 * mp.pi * r ** 2 * mp.exp(-a * r ** 2) * mp.exp(-r) / mp.sqrt(mp.pi), [0, 5, mp.inf])
             for a, c in terms)
    say(f"STO-{n}G zeta=1: alpha = {[mp.nstr(a, 7) for a, _ in terms]}; raw norm {mp.nstr(S, 10)}; "
        f"overlap with Slater 1s {mp.nstr(ov, 10)}")


# ---- 2. single-Gaussian closed form and its large-l expansion
def il(l, z):
    return mp.sqrt(mp.pi / (2 * z)) * mp.besseli(l + mp.mpf(1) / 2, z)


def F_closed_quad(T0, T1, l):
    """Quadrature form; loses accuracy at large l (u^(2l) piles up at u = 1).  Toy checks only."""
    return (2 * l + 1) * mp.quad(lambda u: mp.exp(-T0 * u ** 2) * il(l, T1 * u ** 2), [0, 1])


def F_closed(T0, T1, l):
    """Series form, exact: (2l+1) sum_k T1^(l+2k) / (2^k k! (2l+2k+1)!!) * e^(-T0) M(1, l+2k+3/2, T0) / (2l+4k+1)."""
    tot, k = mp.mpf(0), 0
    while True:
        t = (T1 ** (l + 2 * k) / (2 ** k * mp.factorial(k) * mp.fac2(2 * l + 2 * k + 1)) * mp.exp(-T0)
             * mp.hyp1f1(1, l + 2 * k + mp.mpf(3) / 2, T0) / (2 * l + 4 * k + 1))
        tot += t
        k += 1
        if abs(t) < abs(tot) * mp.mpf(10) ** (-(mp.mp.dps - 3)):
            return (2 * l + 1) * tot


def F_direct(T0, T1, l):
    f = lambda x: mp.quad(lambda u: mp.exp(-(T0 - T1 * x) * u ** 2), [0, 1])
    return (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: f(x) * mp.legendre(l, x), [-1, 0, 1])


def F_weber(sP, sQ, a, l):
    """Same coefficient via k-space: (2/pi) int dk e^{-a k^2} j_l(k sP) j_l(k sQ), normalized to T0, T1 form.
    For unit-charge Gaussians at P, Q with 1/(4p)+1/(4q) = a: (rho1|rho2) = erf(sqrt(mu)|P-Q|)/|P-Q|, mu = 1/(4a)."""
    return (2 * (2 * l + 1) / mp.pi) * mp.quad(lambda k: mp.exp(-a * k ** 2) * mp.sqrt(mp.pi / (2 * k * sP))
                                                * mp.besselj(l + 0.5, k * sP) * mp.sqrt(mp.pi / (2 * k * sQ))
                                                * mp.besselj(l + 0.5, k * sQ), [0, 5, 20, mp.inf])


# toy, not a scored case: mu = 1/(4a) = 0.35, |P| = 1.1, |Q| = 0.8
a_toy = mp.mpf(1) / (4 * mp.mpf('0.35'))
sP, sQ = mp.mpf('1.1'), mp.mpf('0.8')
mu = 1 / (4 * a_toy)
T0t, T1t = mu * (sP ** 2 + sQ ** 2), 2 * mu * sP * sQ
pref = 2 * mp.sqrt(mu / mp.pi)
for l in (0, 3, 8):
    c = pref * F_closed(T0t, T1t, l)
    cq = pref * F_closed_quad(T0t, T1t, l)
    dd = pref * F_direct(T0t, T1t, l)
    w = F_weber(sP, sQ, a_toy, l)
    say(f"toy l={l}: series {mp.nstr(c, 15)}  u-quad {mp.nstr(cq, 15)}  direct {mp.nstr(dd, 15)}  weber {mp.nstr(w, 15)}")

# symbolic large-l expansion of R_l = F_l / (e^{-T0} T1^l / (2l+1)!!)
e, T0s, T1s = sp.symbols('epsilon T0 T1', positive=True)
lsym = 1 / e


def Rk(k, J):
    # term k of the i_l series, times its Kummer factor truncated at J terms, relative to the leading form
    b = lsym + 2 * k + sp.Rational(3, 2)
    M = sum(T0s ** j / sp.rf(b, j) for j in range(J))
    dbl = sp.Mul(*[1 / (2 * lsym + 2 * i + 1) for i in range(1, k + 1)])   # (2l+1)!!/(2l+2k+1)!!
    return T1s ** (2 * k) / (2 ** k * sp.factorial(k)) * dbl * (2 * lsym + 1) / (2 * lsym + 4 * k + 1) * M


R = sum(Rk(k, 4) for k in range(3))
ser = sp.series(R, e, 0, 3).removeO()
c1 = sp.simplify(ser.coeff(e, 1))
c2 = sp.simplify(ser.coeff(e, 2))
say(f"R_l = 1 + c1/l + c2/l^2 + O(l^-3):  c1 = {c1};  c2 = {sp.expand(c2)}")
q = sp.series((R.subs(e, e / (1 + e))) / R, e, 0, 4).removeO()   # R_{l+1}/R_l
say(f"q_l = R_(l+1)/R_l = {sp.expand(q)}")
qc2, qc3 = sp.expand(q.coeff(e, 2)), sp.expand(q.coeff(e, 3))

GCASES = {'g1': (mp.mpf(1), mp.mpf(1)), 'g2': (mp.mpf('0.3'), mp.mpf(2))}   # (alpha_X, alpha_H)
gout = {}
for name, (aX, aH) in GCASES.items():
    p = aX + aH
    muq = p / 2
    s = aH / p
    T0, T1 = muq * D ** 2 * 2 * s ** 2, 2 * muq * D ** 2 * s ** 2
    v1, v2 = float(c1.subs({T0s: T0, T1s: T1})), float(c2.subs({T0s: T0, T1s: T1}))
    w2, w3 = float(qc2.subs({T0s: T0, T1s: T1})), float(qc3.subs({T0s: T0, T1s: T1}))
    gout[name] = dict(aX=float(aX), aH=float(aH), T0=float(T0), T1=float(T1), c1=v1, c2=v2, q2=w2, q3=w3)
    say(f"{name}: T0 = T1 = {mp.nstr(T1, 12)}; c1 = {v1:.6g}, c2 = {v2:.6g}; q_l - 1 = {w2:.6g}/l^2 + {w3:.6g}/l^3")


# ---- 3. tail-filter model for STO-nG vs Slater
def V_slater(w, zH):
    N = mp.sqrt(zH ** 3 / mp.pi)
    Q = N * 8 * mp.pi / zH ** 3
    lam = zH
    if w < mp.mpf('1e-6'):
        return Q ** 2 * (5 * lam / 16)
    J = 1 / w - mp.exp(-lam * w) * (1 / w + 11 * lam / 16 + 3 * lam ** 2 * w / 16 + lam ** 3 * w ** 2 / 48)
    return Q ** 2 * J


def leg_slater(zH, l):
    f = lambda x: V_slater(D * mp.sqrt(2 * (1 - x)), zH)
    return (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: f(x) * mp.legendre(l, x), mp.linspace(-1, 1, 9))


def leg_gauss(n, zH, l):
    terms, _ = contraction(n, zH)
    tot = 0
    for ai, ci in terms:
        for aj, cj in terms:
            qi, qj = ci * (mp.pi / ai) ** 1.5, cj * (mp.pi / aj) ** 1.5
            m = ai * aj / (ai + aj)
            T = 2 * m * D ** 2
            tot += qi * qj * 2 * mp.sqrt(m / mp.pi) * F_closed(T, T, l)
    return tot


def chi_slater(z, r):
    return mp.sqrt(z ** 3 / mp.pi) * mp.exp(-z * r)


def chi_gauss(n, z, r):
    terms, _ = contraction(n, z)
    return sum(c * mp.exp(-a * r ** 2) for a, c in terms)


stored = json.load(open(OUT / "wa5_cusp.json"))
SETS = {'E1': (mp.mpf(1), mp.mpf(1)), 'E2': (mp.mpf('2.25'), mp.mpf('1.24'))}
LS = [4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 48, 56, 60]
model = {}
for sname, (zX, zH) in SETS.items():
    S = stored[f"{sname} (1.809, 1.809)"]["F"]
    lsv = {l: leg_slater(zH, l) for l in LS}
    cS = chi_slater(zX, D) ** 2
    say(f"{sname}: model-Slater / stored S_l: " +
        ", ".join(f"l={l}: {mp.nstr(cS * lsv[l] / S[l], 6)}" for l in LS))
    for n in (3, 6):
        terms, _ = contraction(n, zH)
        amax = max(a for a, _ in terms)
        cx = (chi_gauss(n, zX, D) / chi_slater(zX, D)) ** 2
        rat = {l: cx * leg_gauss(n, zH, l) / lsv[l] for l in LS}
        lstar = next((l for l in LS if rat[l] < 0.5), None)
        model[f"{sname}_n{n}"] = dict(cx=float(cx), alpha_max_H=float(amax), ratio={l: float(v) for l, v in rat.items()},
                                     lstar_grid=lstar)
        say(f"{sname} STO-{n}G: chi_X ratio^2 {mp.nstr(cx, 8)}, alpha_max(H) {mp.nstr(amax, 6)}, "
            f"D sqrt(alpha_max) {mp.nstr(D * mp.sqrt(amax), 5)}; model ratio " +
            ", ".join(f"{l}: {mp.nstr(v, 4)}" for l, v in rat.items()))

(OUT / "wg_derive.json").write_text(json.dumps(dict(single=gout, model=model, c1=str(c1), c2=str(sp.expand(c2)),
                                                     q=str(sp.expand(q))), indent=1))
(OUT / "wg_derive.log").write_text("\n".join(log) + "\n")
