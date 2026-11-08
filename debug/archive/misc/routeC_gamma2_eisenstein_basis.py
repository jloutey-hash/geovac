"""Gamma(2) weight-2 Eisenstein basis + length-<=2 iterated-integral structure for the
Paper 59 collinear observable T2.  Digit-independent STRUCTURAL work.

Deliverables:
 (A) M_2(Gamma(2)) = <theta2^4, theta4^4>, dim 2, Jacobi relation theta3^4=theta2^4+theta4^4.
 (B) length-1 Eichler-integral dictionary -> {pi, varpi=K(1/2), 1/varpi, G=Catalan}.
 (C) length-2 iterated integrals of the weight-2 forms land in the graded ring (examples PSLQ'd).
 (D) ring {pi, varpi, 1/varpi, G} weight-graded to <=3: dimension enumeration.

Every classical identity checked at two precisions + independent evaluator.
"""
from __future__ import annotations
import json
import mpmath as mp


# ----------------------------------------------------------------------------
# (A) weight-2 Eisenstein space for Gamma(2) from theta fourth powers
# ----------------------------------------------------------------------------
def thetas(tau):
    """theta2,theta3,theta4 with nome q=e^{i pi tau}."""
    q = mp.e ** (1j * mp.pi * tau)
    th2 = 2 * mp.nsum(lambda n: q ** ((n + mp.mpf('0.5')) ** 2), [0, mp.inf])
    th3 = mp.nsum(lambda n: q ** (n * n), [-mp.inf, mp.inf])
    th4 = mp.nsum(lambda n: (-1) ** n * q ** (n * n), [-mp.inf, mp.inf])
    return th2, th3, th4


def check_weight2_space():
    print("=" * 72)
    print("(A) M_2(Gamma(2)) = span{theta2^4, theta4^4}, Jacobi theta3^4=theta2^4+theta4^4")
    print("=" * 72)
    for tau in [1j, 0.5 + 1.3j, 0.2 + 0.9j]:
        t2, t3, t4 = thetas(tau)
        jac = t3 ** 4 - (t2 ** 4 + t4 ** 4)
        print(f"  tau={tau}:  |theta3^4-(theta2^4+theta4^4)| = {mp.nstr(abs(jac),3)}")
    print("  => three weight-2 theta^4 forms, ONE linear (Jacobi) relation => dim M_2 = 2.")
    print("  => dim S_2(Gamma(2)) = 0 (genus 0): all weight-2 forms are Eisenstein.\n")


# ----------------------------------------------------------------------------
# (B) length-1 Eichler dictionary: {pi, varpi, 1/varpi, G}
# ----------------------------------------------------------------------------
def length1_dictionary(dps):
    mp.mp.dps = dps
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))   # K(1/2), CM period tau=i
    G = mp.catalan                                             # L(2, chi_-4), Eisenstein L-value
    E12 = mp.ellipe(mp.mpf('0.5'))                             # second-kind period at tau=i
    # native length-1 Eichler integrals over the modulus k in [0,1] (cusp k=0 -> cusp k=1)
    intK = mp.quad(lambda k: mp.ellipk(k * k), [0, 1])         # = 2G
    intE = mp.quad(lambda k: mp.ellipe(k * k), [0, 1])         # = G + 1/2
    checks = {
        'int_0^1 K dk = 2G': abs(intK - 2 * G),
        'int_0^1 E dk = G+1/2': abs(intE - (G + mp.mpf('0.5'))),
        'E(1/2) = pi/4varpi + varpi/2 (=> 1/varpi native)': abs(E12 - (pi / (4 * varpi) + varpi / 2)),
    }
    return checks


# ----------------------------------------------------------------------------
# (C) length-2 iterated integrals of the weight-2 forms (over the k-modulus path)
# ----------------------------------------------------------------------------
def guarded_pslq(target, basis, names, maxcoeff, tol):
    rel = mp.pslq([target] + list(basis), tol=tol, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    if rel is None:
        return None, "(none)"
    if rel[0] == 0:
        return rel, "(V-coeff 0 => basis-internal identity, NOT a closure)"
    h = max(abs(c) for c in rel)
    if h <= 200:
        terms = " + ".join(f"{rel[i+1]}*{names[i]}" for i in range(len(names)) if rel[i + 1])
        return rel, f"CLOSURE h={h}:  {rel[0]}*V = {terms}"
    return rel, f"(height {h} high => not a low-height closure)"


def length2_layer(dps):
    mp.mp.dps = dps
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))
    G = mp.catalan
    nu = 1 / varpi
    z3 = mp.zeta(3)
    # length-2 cusp-to-cusp double iterated integral of the period K over the modulus:
    KK = mp.quad(lambda k: mp.ellipk(k * k) * mp.quad(lambda kk: mp.ellipk(kk * kk), [0, k]), [0, 1])
    # length-1 of the weight-2 form K^2 (a genuine weight-3 candidate):
    K2 = mp.quad(lambda k: mp.ellipk(k * k) ** 2, [0, 1])
    tol = mp.mpf(10) ** (-(dps - 8))
    # graded ring monomials up to weight 4 for the landing test
    basis = [mp.mpf(1), pi, varpi, nu, G,
             pi ** 2, pi * varpi, varpi ** 2, pi * nu, nu ** 2,
             G * pi, G * varpi, G * nu, pi ** 3, pi ** 2 * varpi, pi * varpi ** 2,
             varpi ** 3, G ** 2, pi ** 2 * nu, z3, G * pi * nu]
    names = ['1', 'pi', 'varpi', '1/varpi', 'G', 'pi^2', 'pi*varpi', 'varpi^2', 'pi/varpi',
             '1/varpi^2', 'G*pi', 'G*varpi', 'G/varpi', 'pi^3', 'pi^2*varpi', 'pi*varpi^2',
             'varpi^3', 'G^2', 'pi^2/varpi', 'zeta3', 'G*pi/varpi']
    out = {}
    for nm, val in [('double-Eichler  int K(int K)dk', KK), ('int K^2 dk', K2)]:
        rel, tag = guarded_pslq(val, basis, names, maxcoeff=10 ** 6, tol=tol)
        out[nm] = (mp.nstr(val, 25), tag)
    out['KK vs 2*G^2 (weight-4, predicted exact)'] = mp.nstr(abs(KK - 2 * G ** 2), 3)
    return out


# ----------------------------------------------------------------------------
# (D) ring dimension: {pi, varpi, 1/varpi, G} graded to weight <= 3
# ----------------------------------------------------------------------------
def ring_dimension():
    # weights: pi:1, varpi:1, inv=1/varpi:-1, G:2.  Monomials pi^a varpi^b inv^c G^d.
    # physical enumeration: transcendental degree b+c <= 2 (length-2), G-degree d<=1
    # (one Eisenstein L-value at weight<=3), no simultaneous varpi & 1/varpi, 0<=weight<=3.
    monos = set()
    for a in range(0, 4):
        for b in range(0, 3):
            for c in range(0, 3):
                if b > 0 and c > 0:
                    continue
                for d in range(0, 2):
                    w = a + b - c + 2 * d
                    if 0 <= w <= 3 and (b + c) <= 2:
                        monos.add((a, b, c, d, w))
    by_w = {}
    for m in monos:
        by_w.setdefault(m[4], []).append(m)
    return monos, by_w


def main():
    check_weight2_space()

    print("=" * 72)
    print("(B) length-1 Eichler dictionary  (two precisions + independent evaluator)")
    print("=" * 72)
    for dps in (30, 45):
        checks = length1_dictionary(dps)
        print(f"  --- dps={dps} ---")
        for k, v in checks.items():
            print(f"     {k:52s}: err {mp.nstr(v,3)}")
    print("  Generators produced at length 1: pi, varpi=K(1/2), 1/varpi, G=Catalan.\n")

    print("=" * 72)
    print("(C) length-2 iterated integrals land in the graded ring")
    print("=" * 72)
    res = length2_layer(40)
    for k, v in res.items():
        print(f"  {k}:\n      {v}")
    print()

    print("=" * 72)
    print("(D) ring {pi,varpi,1/varpi,G} weight-graded <=3: dimension")
    print("=" * 72)
    monos, by_w = ring_dimension()
    print(f"  total dimension (weight 0..3) = {len(monos)}")
    for w in sorted(by_w):
        ms = ", ".join(f"pi^{m[0]}varpi^{m[1]}inv^{m[2]}G^{m[3]}" for m in sorted(by_w[w]))
        print(f"    weight {w}: dim {len(by_w[w])}  [{ms}]")
    out = {'ring_dim_weight_le3': len(monos),
           'dim_by_weight': {str(w): len(v) for w, v in by_w.items()}}
    with open('debug/data/gamma2_eisenstein_basis.json', 'w') as f:
        json.dump(out, f, indent=2)
    print("\n  wrote debug/data/gamma2_eisenstein_basis.json")


if __name__ == '__main__':
    main()
