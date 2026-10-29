"""Large-L structure of a recurrence sum_k p_k(L) F_{L+k} = 0 mod p (record B2; calibration on BU's one-electron
recurrence).  Characteristic roots lam (leading coefficients), each mapped to x = (lam + 1/lam)/2; at each rational root
the exponents beta of the formal solutions lam^L L^beta (indicial polynomial of sum_j q_j(L) Delta^j, q_j(L) =
sum_k C(k, j) lam^k p_k(L), on the top edge of the Newton polygon).  A singular term (x0 - x)^alpha maps to
beta = -alpha - 1/2.    python debug/b2_asym.py <recurrence.json>
"""
import json
import math
import sys

import flint


def ratrec(a, p):
    r0, r1 = p, int(a) % p
    s0, s1 = 0, 1
    B = math.isqrt(p // 2)
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    from fractions import Fraction
    return Fraction(r1 * (1 if s1 > 0 else -1), abs(s1))


def analyse(polys, p):
    P = [flint.nmod_poly(c, p) for c in polys]
    m = len(P) - 1
    d = max(q.degree() for q in P)
    chi = flint.nmod_poly([int(q.coeffs()[d]) if q.degree() == d else 0 for q in P], p)
    out = []
    for f, mult in chi.factor()[1]:
        if f.degree() != 1:
            out.append(dict(root=f'irreducible deg {f.degree()}', mult=mult))
            continue
        lam = -f.coeffs()[0] / f.coeffs()[1]
        # q_j(L) = sum_k C(k, j) lam^k p_k(L)
        q = []
        for j in range(m + 1):
            s = flint.nmod_poly([], p)
            for k in range(j, m + 1):
                s += P[k] * (flint.nmod(math.comb(k, j), p) * lam ** k)
            q.append(s)
        top = max(qj.degree() - j for j, qj in enumerate(q) if qj.degree() >= 0)
        ind = flint.nmod_poly([], p)
        for j, qj in enumerate(q):
            if qj.degree() >= 0 and qj.degree() - j == top:
                fall = flint.nmod_poly([1], p)
                for t in range(j):
                    fall *= flint.nmod_poly([-t, 1], p)
                ind += fall * flint.nmod_poly([int(qj.coeffs()[qj.degree()])], p)
        betas = []
        for g, mm in ind.factor()[1]:
            betas.append((str(ratrec(-g.coeffs()[0] / g.coeffs()[1], p)) if g.degree() == 1 else f'deg{g.degree()}', mm))
        lr = ratrec(lam, p)
        out.append(dict(root=str(lr), x=str((lr + 1 / lr) / 2) if lr != 0 else None, mult=mult, betas=betas))
    return out


if __name__ == '__main__':
    d = json.load(open(sys.argv[1]))
    for r in analyse(d['polys'], d['p']):
        print(r)
