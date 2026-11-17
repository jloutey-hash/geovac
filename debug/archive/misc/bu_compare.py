"""Record BU, FINAL COMPARISON ONLY: compute NA's guessed operator exactly as debug/na_run.py did
(G.nullity(B, 8, 31, 30) on the generic sequence L = 0..1000) and test whether it right-divides the derived operator
data/bu_angle_v2_<p>_9_30.json, by exact Ore right division over F_p(x) in standard form sum c_k(x) d^k.
Also tests the derived operator against the guessed one's solution count (orders) and reports the left quotient."""
import json

import flint

import na_modanalyze as A
import na_modguess as G

PRIMES = [(1 << 61) - 1, 4611686018427387847]


def dmul(c):
    """d/dx o (sum c_k d^k) = sum (c_k' d^k + c_k d^{k+1})."""
    p = c[0].modulus()
    out = [flint.nmod_poly([], p) for _ in range(len(c) + 1)]
    for k, ck in enumerate(c):
        out[k] += ck.derivative()
        out[k + 1] += ck
    return out


def trim(c):
    while len(c) > 1 and c[-1].degree() < 0:
        c = c[:-1]
    return c


def right_divide(D, Gs):
    """Fraction-free right division of D by Gs (standard forms). Returns (is_exact, quotient pieces)."""
    p = Gs[0].modulus()
    g = Gs[-1]
    rG = len(Gs) - 1
    R = [x for x in D]
    scale = flint.nmod_poly([1], p)
    quot = []                                  # (k, numerator poly, scale) : R_scaled = sum ... + remainder
    while len(trim(R)) - 1 >= rG and any(x.degree() >= 0 for x in R):
        R = trim(R)
        k = len(R) - 1 - rG
        lead = R[-1]
        # g^? * R - lead * d^k o G  removes the top term
        DG = list(Gs)
        for _ in range(k):
            DG = dmul(DG)
        R = [g * x for x in R]
        R = [x - lead * y for x, y in zip(R + [flint.nmod_poly([], p)] * (len(DG) - len(R)), DG)]
        scale = scale * g
        quot.append((k, lead))
        R = trim(R)
        if len(R) - 1 < rG:
            break
    exact = all(x.degree() < 0 for x in R)
    return exact, quot, R


if __name__ == '__main__':
    res = {}
    for p in PRIMES:
        seq = json.load(open(f'data/bu_check_seq_{p}.json'))['gen']
        B = G.Basis([flint.nmod(v, p) for v in seq], p)
        nul, keys, (N, nrow) = G.nullity(B, 8, 31, 30)
        assert nul == 1, nul
        gvec = [N[i, 0] for i in range(N.nrows())]
        Gs = A.standard_form(keys, gvec, p)
        d = json.load(open(f'data/bu_angle_v2_{p}_9_30.json'))
        Ds = A.standard_form([tuple(k) for k in d['keys']], [flint.nmod(v, p) for v in d['vecs'][0]], p)
        exact, quot, rem = right_divide(Ds, Gs)
        # quotient as a first-order operator q1 d + q0 (rational): report degrees of its numerators
        info = dict(guessed_nullity=nul, guessed_orders=len(Gs) - 1, guessed_degrees=[c.degree() for c in Gs],
                    derived_orders=len(Ds) - 1, derived_degrees=[c.degree() for c in Ds],
                    right_division_exact=exact, quotient_steps=[(k, q.degree()) for k, q in quot],
                    remainder_degrees=[c.degree() for c in rem])
        # does the guessed operator already annihilate what the derived one annihilates? (G is order 8 < 9)
        # decisive: exact right division
        print(p, info, flush=True)
        res[str(p)] = info
    json.dump(res, open('data/bu_compare.json', 'w'), indent=1)
