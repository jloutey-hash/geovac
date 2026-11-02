"""Record IC, Phase 1: the controls of the order-one peel, measured before registering.  None touches (XY|XZ).

Legendre-series controls (geometry of the unscored case S1: D1 = 9/5, D2 = 13/5, rho = D1/D2 = 9/13; w = |Y - Z|,
x = (D1^2 + D2^2 - w^2) / (2 D1 D2)).  Each sequence -> minimal bond-angle operator (debug/na_modguess.py, the AO3
search, unchanged) -> pullback to w -> peel (debug/intclass_lib.py).
  P1  sum rho^L C(2L,L)/4^L P_L(x) = (1/pi) int_0^rho dt / sqrt(t (rho - t)(1 - 2xt + t^2)): a period of a family of
      elliptic curves                                                      expected STUCK, order-2 residual, log
  P3  P1 + N1 (one component)                                              expected STUCK at order 2 after peeling
  N1  sum rho^L/(L+1) P_L(x) = log of a rational function of w (/rho)      expected TRIANGULAR
  N2  the two-centre Coulomb integral (1s_Y 1s_Y | 1s_Z 1s_Z), zeta = 1, as a function of the angle: its three
      constant components (1, e^{-2(D2-D1)}, e^{-2(D1+D2)}); the exponential components grow factorially
      (formal series), as the (XY|XZ) components may                       expected TRIANGULAR (each)
  N3  the Yukawa function e^{-2w}/w: its two constant components            expected TRIANGULAR (each)
Operator-level controls:
  P4  P1's w-operator twisted by e^{3w}                                    expected STUCK
  P5  Y's minimal operator M13 (Paper 59 sec:density_ring; pieces: K3 R_V order 4, elliptic z1, z2 order 2,
      order-one F's), from tests/data/paper59_y_2026_10_26.json, variable x = t^2   expected STUCK
  P6  the radial core's Picard-Fuchs operator eq:pf at rho = 1/5 (order 4, irreducible, irregular)  expected STUCK
  B1  Bessel K0: D y'' + y' - D y (non-Liouvillian, genus-zero confluent)  expected STUCK (shows STUCK != elliptic)
  B2  a triangular operator of the same shape as N2's, built by hand: (d - 2 + 1/w) o (w^2 d - ...) [see code]
usage: python debug/intclass_controls.py [prime]   ->  debug/data/intclass_controls_<p>.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint

sys.path.insert(0, str(Path(__file__).parent))
import intclass_lib as I
import na_modguess as G
import na_modanalyze as A

ROOT = Path(__file__).resolve().parent.parent
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1073741789
F = lambda q: I.to_fp(Fraction(q), p)
D1, D2 = Fraction(9, 5), Fraction(13, 5)
RHO = D1 / D2
XSTAR = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
A0, A2 = F(XSTAR), F(Fraction(-1) / (2 * D1 * D2))
LMAX = 900


def guess_x_operator(seq, Rmax=16, margin=30):
    m = I.minimal_operator(seq, p, Rmax, margin)
    if m is None:
        return None
    return dict(R=m['order'], Q0=m['degree'], found=m['found'], c=m['c'])


def to_w(c):
    return I.pullback_quadratic(c, A0, A2, p)


def run(name, c, expect, extra=None):
    t0 = time.time()
    r = I.peel(c, p, verbose=False)
    out = dict(name=name, expect=expect, verdict=r['verdict'], order=len(c) - 1,
               maxdeg=max(x.degree() for x in c), steps=r['steps'], time=round(time.time() - t0, 2))
    if r['verdict'] == 'STUCK':
        out['residual_order'] = r['residual_order']
        out['residual_diag'] = r['residual_diag']
        out['adjoint_has_order1'] = r['adjoint_has_order1']
        out['class'] = I.classify_residual(r, p)
        out['residual_diag'] = None
    if extra:
        out.update(extra)
    out['pass'] = (r['verdict'] == expect)
    print(f"{name:4s} order {out['order']:2d} deg {out['maxdeg']:4d}: {r['verdict']:10s} (expect {expect}) "
          f"steps {len(r['steps'])} {'res ' + str(out.get('residual_order')) if 'residual_order' in out else ''} "
          f"{out['time']}s  {'PASS' if out['pass'] else 'FAIL'}", flush=True)
    if 'class' in out:
        c_ = out['class']
        print(f"      class {c_['type']}, clean {c_['clean']}, lam {c_['lam']}, adj-order1 {c_['adjoint_has_order1']}, "
              f"flags {c_['flags']}", flush=True)
    return out


# ------------------------------------------------------------------------------------ sequences
def seq_P1():
    out, cb = [], Fraction(1)
    for L in range(LMAX):
        out.append(F(RHO ** L * cb))
        cb = cb * Fraction((2 * L + 1) * (2 * L + 2), (L + 1) ** 2 * 4)
    return out


def seq_N1():
    return [F(RHO ** L / (L + 1)) for L in range(LMAX)]


class S:
    """truncated eps-series mod p, length 4"""
    n = 4

    def __init__(self, c):
        self.c = [x % p for x in c] + [0] * (S.n - len(c))

    def __add__(s, o):
        return S([(a + b) % p for a, b in zip(s.c, o.c)])

    def __sub__(s, o):
        return S([(a - b) % p for a, b in zip(s.c, o.c)])

    def __mul__(s, o):
        if isinstance(o, int):
            return S([a * o % p for a in s.c])
        r = [0] * S.n
        for i, a in enumerate(s.c):
            if a:
                for j in range(S.n - i):
                    r[i + j] = (r[i + j] + a * o.c[j]) % p
        return S(r)


def ser_inv_lin(a, k0):
    """1/(a (k0 + eps)) as series"""
    inv = pow(a * k0 % p, p - 2, p)
    q = pow(k0, p - 2, p)
    return S([inv * pow(-q % p, i, p) % p for i in range(S.n)])


def ser_exp(a):
    """e^{a eps}"""
    out, t, f = [], 1, 1
    for i in range(S.n):
        out.append(t * pow(f, p - 2, p) % p)
        t = t * a % p
        f = f * (i + 1)
    return S(out)


def yukawa_components(Lmax, kappa=2):
    """E(kappa) = e^{-kappa w}/w = sum_L (2L+1) kappa i_L(kappa a) k~_L(kappa b) e^{-kappa b} P_L, a = D1 < b = D2,
    i_L = e^{z} al_L + e^{-z} be_L, k~_L(z) = e^{z} k_L(z) (2/pi) polynomial in 1/z.  Series in eps (kappa -> kappa +
    eps).  Returns (cA, cB): per L, eps-series of the coefficient of e^{-kappa(b-a)} and of e^{-kappa(a+b)}."""
    a, b = F(D1), F(D2)
    za = ser_inv_lin(a, kappa)        # 1/z, z = (kappa + eps) a
    zb = ser_inv_lin(b, kappa)
    half = pow(2, p - 2, p)
    one = S([1])
    al = [za * half, (za - za * za) * half]
    be = [za * (p - half), (za + za * za) * half]
    kt = [zb, zb + zb * zb]
    for L in range(1, Lmax):
        al.append(al[L - 1] - za * al[L] * (2 * L + 1))
        be.append(be[L - 1] - za * be[L] * (2 * L + 1))
        kt.append(kt[L - 1] + zb * kt[L] * (2 * L + 1))
    kap = S([kappa, 1])
    eA = ser_exp((a - b) % p)         # e^{eps (a - b)}
    eB = ser_exp((-a - b) % p)
    cA = [kap * al[L] * kt[L] * eA * (2 * L + 1) for L in range(Lmax)]
    cB = [kap * be[L] * kt[L] * eB * (2 * L + 1) for L in range(Lmax)]
    return cA, cB


def coulomb_components(Lmax):
    """J(w) = 1/w - [E - (11/8) E' + (3/4) E'' - (1/6) E'''], E^(n) = (-d/dkappa)^n E at kappa = 2 (zeta = 1):
    the 1s-1s two-centre Coulomb integral.  Components: constant 1, e^{-2(b-a)}, e^{-2(a+b)}."""
    cA, cB = yukawa_components(Lmax)
    fac = [1, p - 1, 2, p - 6]        # (-1)^n n!
    w = [1, F(Fraction(11, 8)), F(Fraction(3, 4)), F(Fraction(1, 6))]

    def comb_(c):
        return [(-sum(w[n] * fac[n] % p * c[L].c[n] for n in range(4))) % p for L in range(Lmax)]
    one = [F(D1 ** L / D2 ** (L + 1)) for L in range(Lmax)]
    return one, comb_(cA), comb_(cB)


def e1_components(Lmax, kappa=2):
    """E1(kappa w) = E(kappa)/kappa - int_kappa^inf E(s)/s^2 ds,  E(s) = e^{-s w}/w.  With u = 1/s,
    E_L(s)/s^2 = (2L+1) u [e^{-s(b-a)} al_L k~_L + e^{-s(a+b)} be_L k~_L] (al, be, k~ polynomials in u), and
    int_kappa^inf s^{-n} e^{-cs} ds = kappa^{1-n} (alpha_n e^{-c kappa} + beta_n E1(c kappa)),
    E_n = (e^{-z} - z E_{n-1})/(n-1).  Components: e^{-k(b-a)}, e^{-k(a+b)}, E1(k(b-a)), E1(k(a+b))."""
    a, b = F(D1), F(D2)
    ia, ib = pow(a, p - 2, p), pow(b, p - 2, p)
    P = lambda c: flint.nmod_poly([x % p for x in c], p)
    u = P([0, 1])
    half = pow(2, p - 2, p)
    al = [P([0, ia * half]), P([0, ia * half, -ia * ia * half])]
    be = [P([0, -ia * half]), P([0, ia * half, ia * ia * half])]
    kt = [P([0, ib]), P([0, ib, ib * ib])]
    for L in range(1, Lmax):
        al.append(al[L - 1] - u * al[L] * ((2 * L + 1) * ia % p))
        be.append(be[L - 1] - u * be[L] * ((2 * L + 1) * ia % p))
        kt.append(kt[L - 1] + u * kt[L] * ((2 * L + 1) * ib % p))
    nmax = 2 * Lmax + 8
    out = {}
    for tag, c in (('m', (b - a) % p), ('p', (a + b) % p)):
        z = c * kappa % p
        alpha, beta = [0, 0], [0, 1]
        for n in range(2, nmax):
            inv = pow(n - 1, p - 2, p)
            alpha.append((1 - z * alpha[-1]) * inv % p)
            beta.append((-z * beta[-1]) * inv % p)
        ik = pow(kappa, p - 2, p)
        kp = [1]                       # kappa^{1-n}: n = 1 -> 1, n -> ik^{n-1}
        for n in range(1, nmax):
            kp.append(kp[-1] * ik % p)
        cE, cE1 = [], []
        for L in range(Lmax):
            poly = (al[L] if tag == 'm' else be[L]) * kt[L]
            ev = poly(flint.nmod(ik, p))
            top = (2 * L + 1) * int(ev) % p                      # E_L(kappa)/kappa = (2L+1) i_L k~_L
            cs = (u * poly).coeffs()
            s1 = s2 = 0
            for n, v in enumerate(cs):
                if n == 0 or int(v) == 0:
                    continue
                w_ = (2 * L + 1) * int(v) % p * kp[n - 1] % p
                s1 = (s1 + w_ * alpha[n]) % p
                s2 = (s2 + w_ * beta[n]) % p
            cE.append((top - s1) % p)
            cE1.append((-s2) % p)
        out['exp_' + tag], out['E1_' + tag] = cE, cE1
    return out


# ------------------------------------------------------------------------------------ operator-level controls
def m13_dform():
    d = json.loads((ROOT / 'tests' / 'data' / 'paper59_y_2026_10_26.json').read_text())
    cth = [flint.nmod_poly([int(v) % p for v in c], p) for c in d['M13']]
    return I.from_theta(cth, p, s=2)        # theta = t d/dt = 2 x d/dx, x = t^2


def eqpf(rho=Fraction(1, 5)):
    r = F(rho)
    P = lambda *c: flint.nmod_poly([x % p for x in c], p)
    # D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L
    return [P(0, -(1 - r)), P((1 - 2 * r) % p), P(0, (1 - 2 * r) % p), P(2 * r % p), P(0, r)]


def bessel():
    P = lambda *c: flint.nmod_poly([x % p for x in c], p)
    return [P(0, -1), P(1), P(0, 1)]


def twist(c, lam):
    return I.content_reduce(I.shift_lambda(c, (-lam) % p, p))   # operator for e^{lam w} * (solutions of c)


if __name__ == '__main__':
    T0 = time.time()
    res = dict(prime=p, case='S1 geometry (D1, D2) = (9/5, 13/5)', results=[])
    ops = {}
    for name, seq, expect in [('P1', seq_P1(), 'STUCK'), ('N1', seq_N1(), 'TRIANGULAR')]:
        g = guess_x_operator(seq)
        ops[name] = g
        print(name, 'x-operator', g['R'], g['Q0'], flush=True)
        res['results'].append(run(name, to_w(g['c']), expect, dict(x_order=g['R'], x_degree=g['Q0'])))
    s3 = [(a + b) % p for a, b in zip(seq_P1(), seq_N1())]
    g = guess_x_operator(s3)
    print('P3 x-operator', g['R'], g['Q0'], flush=True)
    res['results'].append(run('P3', to_w(g['c']), 'STUCK', dict(x_order=g['R'], x_degree=g['Q0'])))
    one, cA, cB = coulomb_components(LMAX)
    for nm, sq in [('N2a', one), ('N2b', cA), ('N2c', cB)]:
        g = guess_x_operator(sq)
        print(nm, 'x-operator', g and (g['R'], g['Q0']), flush=True)
        res['results'].append(run(nm, to_w(g['c']), 'TRIANGULAR', dict(x_order=g['R'], x_degree=g['Q0'])))
    yA, yB = yukawa_components(LMAX)
    for nm, sq in [('N3a', [x.c[0] for x in yA]), ('N3b', [x.c[0] for x in yB])]:
        g = guess_x_operator(sq)
        print(nm, 'x-operator', g and (g['R'], g['Q0']), flush=True)
        res['results'].append(run(nm, to_w(g['c']), 'TRIANGULAR', dict(x_order=g['R'], x_degree=g['Q0'])))
    for nm, sq in sorted(e1_components(300).items()):
        g = guess_x_operator(sq)
        print('N5' + nm, 'x-operator', g and (g['R'], g['Q0']), flush=True)
        res['results'].append(run('N5' + nm, to_w(g['c']), 'TRIANGULAR', dict(x_order=g['R'], x_degree=g['Q0'])))
    s4 = [F(RHO ** L / (L + 1) ** 2) for L in range(LMAX)]
    g = guess_x_operator(s4)
    print('N4 x-operator', g['R'], g['Q0'], flush=True)
    res['results'].append(run('N4', to_w(g['c']), 'TRIANGULAR', dict(x_order=g['R'], x_degree=g['Q0'])))
    s6 = [F(sum(Fraction(k) * RHO ** L / (L + k) for k in range(1, 6)) + 3 * RHO ** L / Fraction(L + 2) ** 2)
          for L in range(LMAX)]
    g = guess_x_operator(s6)
    print('N6 x-operator', g['R'], g['Q0'], flush=True)
    res['results'].append(run('N6', to_w(g['c']), 'TRIANGULAR', dict(x_order=g['R'], x_degree=g['Q0'])))
    s7 = [(a_ + b_) % p for a_, b_ in zip(s6, seq_P1())]
    g = guess_x_operator(s7)
    print('P7 x-operator', g['R'], g['Q0'], flush=True)
    res['results'].append(run('P7', to_w(g['c']), 'STUCK', dict(x_order=g['R'], x_degree=g['Q0'])))
    y3A, y3B = yukawa_components(LMAX, kappa=3)
    e1c = e1_components(LMAX)
    s8 = [(cA[L] + 2 * y3A[L].c[0] + 5 * y3B[L].c[0] + 7 * e1c['exp_m'][L] + 11 * s6[L]) % p for L in range(LMAX)]
    g = guess_x_operator(s8, Rmax=24)
    print('N8 x-operator', g['R'], g['Q0'], flush=True)
    res['results'].append(run('N8', to_w(g['c']), 'TRIANGULAR', dict(x_order=g['R'], x_degree=g['Q0'])))
    s9 = [(a_ + b_) % p for a_, b_ in zip(s8, seq_P1())]
    g = guess_x_operator(s9, Rmax=24)
    print('P8 x-operator', g['R'], g['Q0'], flush=True)
    res['results'].append(run('P8', to_w(g['c']), 'STUCK', dict(x_order=g['R'], x_degree=g['Q0'])))
    res['results'].append(run('P4', twist(to_w(ops['P1']['c']), 3), 'STUCK'))
    res['results'].append(run('P5', m13_dform(), 'STUCK'))
    res['results'].append(run('P6', eqpf(), 'STUCK'))
    res['results'].append(run('B1', bessel(), 'STUCK'))
    res['all_pass'] = all(r['pass'] for r in res['results'])
    res['time'] = round(time.time() - T0, 1)
    print('ALL PASS' if res['all_pass'] else 'SOME FAIL', res['time'], 's', flush=True)
    out = Path(__file__).parent / 'data' / f'intclass_controls_{p}.json'
    out.write_text(json.dumps(res, indent=1, default=str))
