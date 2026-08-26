"""PSLQ fit of V=T2 against the weight-1/weight-2 Eisenstein/CM-Gamma ring,
with a decoy control (same-magnitude meaningless constant that must NOT also
fit for a closure to be trustworthy).

Weight-1 basis:
    1, pi, varpi = Gamma(1/4)^2/(4 sqrt(pi))  [K(1/2), lemniscate constant]
    E12 = E(1/2)   [complete elliptic E at parameter 1/2]
    K2  = (1+sqrt(2))^(1/2) Gamma(1/8) Gamma(3/8) / (2^(13/4) sqrt(pi))  [disc-8 period]

Weight-2 basis: all pairwise products (with repetition) of the weight-1 basis
EXCLUDING the constant 1 itself as a weight-2 generator (1*1=1 is weight-0,
already covered) -- i.e. all degree-2 monomials in {pi,varpi,E12,K2}, plus the
weight-1 basis itself (so a relation V = sum c_i b_i can land at either weight).
"""
from __future__ import annotations
import sys
import mpmath as mp


def build_basis(dps):
    mp.mp.dps = dps + 20
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))
    E12 = mp.ellipe(mp.mpf(1) / 2)
    K12 = mp.ellipk(mp.mpf(1) / 2)
    K2 = (mp.sqrt(1 + mp.sqrt(2))
          * mp.gamma(mp.mpf(1) / 8) * mp.gamma(mp.mpf(3) / 8)
          / (2 ** mp.mpf('3.25') * mp.sqrt(pi)))
    # sanity: varpi should equal K(1/2) in mpmath's parameter convention
    chk = abs(varpi - K12) / abs(varpi)
    w1 = {'1': mp.mpf(1), 'pi': pi, 'varpi': varpi, 'E12': E12, 'K2': K2}
    names1 = ['pi', 'varpi', 'E12', 'K2']  # exclude '1' from generators of products
    w2 = {}
    for i, a in enumerate(names1):
        for j, b in enumerate(names1):
            if j < i:
                continue
            key = f"{a}*{b}"
            w2[key] = w1[a] * w1[b]
    mp.mp.dps = dps
    return w1, w2, chk


def run_pslq(V, basis_dict, dps, maxcoeff=10 ** 6, label=""):
    names = list(basis_dict.keys())
    vec = [-V] + [basis_dict[n] for n in names]
    mp.mp.dps = dps
    try:
        rel = mp.pslq(vec, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    except Exception as e:
        print(f"  [{label}] pslq exception: {e}")
        return None
    if rel is None:
        print(f"  [{label}] PSLQ: no relation found (maxcoeff={maxcoeff})")
        return None
    cV = rel[0]
    coeffs = rel[1:]
    height = max(abs(c) for c in rel)
    print(f"  [{label}] relation found: V-coeff={cV}  coeffs={dict(zip(names, coeffs))}  height={height}")
    return rel


def two_precision_check(V_str, basis_builder, dps_list, maxcoeff, label):
    print(f"--- {label} ---")
    rels = []
    for dps in dps_list:
        V = mp.mpf(V_str)
        w1, w2, chk = basis_builder(dps)
        combined = {**w1, **w2}
        del combined['1']  # '1' handled implicitly via V-coeff / const term below
        rel = run_pslq(V, combined, dps, maxcoeff=maxcoeff, label=f"dps={dps}")
        rels.append(rel)
    return rels


def main():
    dps_list = [int(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else [30, 40]
    maxcoeff = 10 ** 6

    V_str = sys.argv[0]  # placeholder, real value passed via env / edited below
    print("This module is imported by routeC_pslq_run.py with the actual V string.")


if __name__ == '__main__':
    main()
