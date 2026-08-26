"""Guarded weight-graded PSLQ for the Paper 59 collinear observable.

Fit the RAW period  W = V * pi / 8 = int_0^1 ds int_0^1 dt J(s,t)  (the (8/pi)
prefactor stripped, so W is the natural period) against a weight-graded basis of
monomials in the CM-fibre constants the parameter integral actually visits:

  pi                                   (weight 1)
  varpi = K(1/2) = Gamma(1/4)^2/(4 sqrt pi)   -- disc-4 fibre rho=1/2 (tau=i)  (weight 1)
  P8    = disc-8 fibre period          (weight 1)

NOTE (correction to the first-pass memo basis): E(1/2) is NOT an independent
weight-1 generator. The Legendre relation at the self-dual point m=1/2
(K=K', E=E') is  2 E K - K^2 = pi/2, i.e.  E(1/2) = pi/(4 varpi) + varpi/2  --
weight-MIXED (carries a 1/varpi = weight -1 piece). Including it makes the
weight-2 level rank-deficient and manufactures spurious V-coeff-0 "hits"
(exactly the trap the first-pass PSLQ fell into). So E12 is EXCLUDED as a
generator; if the true closed form needs an E, it appears via this identity.

Basis = all monomials pi^a varpi^b P8^c with 1 <= a+b+c <= wmax (plus the
constant 1). Guarded acceptance:
  * run at several precisions (subsets of the trusted digits);
  * a DECOY of matching magnitude (structureless constant) is fit in parallel
    and MUST NOT find a comparable-height relation;
  * report height; a trustworthy closure = same low-height relation stable
    across precisions for V, absent for the decoy.

Usage: python routeC_pslq_v2.py <V_or_W> <mode:V|W> <trusted_digits> [wmax] [use_P8]
"""
from __future__ import annotations
import sys
import itertools
import mpmath as mp


def cm_constants(dps):
    mp.mp.dps = dps + 25
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))          # K(1/2)
    P8 = (mp.sqrt(1 + mp.sqrt(2)) * mp.gamma(mp.mpf(1) / 8)
          * mp.gamma(mp.mpf(3) / 8) / (2 ** mp.mpf('3.25') * mp.sqrt(pi)))  # disc-8
    return {'pi': pi, 'varpi': varpi, 'P8': P8}


def build_basis(dps, wmax, use_P8):
    c = cm_constants(dps)
    gens = ['pi', 'varpi'] + (['P8'] if use_P8 else [])
    basis = {'1': mp.mpf(1)}
    # monomials of total degree 1..wmax
    for w in range(1, wmax + 1):
        for combo in itertools.combinations_with_replacement(gens, w):
            name = '*'.join(combo)
            val = mp.mpf(1)
            for g in combo:
                val *= c[g]
            basis[name] = val
    mp.mp.dps = dps
    return basis


def guarded_pslq(target, basis, dps, maxcoeff, label):
    names = list(basis.keys())
    mp.mp.dps = dps
    vec = [target] + [basis[n] for n in names]
    try:
        rel = mp.pslq(vec, maxcoeff=maxcoeff, maxsteps=2 * 10 ** 6)
    except Exception as e:
        print(f"    [{label} dps={dps}] exception: {e}", flush=True)
        return None
    if rel is None:
        print(f"    [{label} dps={dps}] NO relation (maxcoeff={maxcoeff})", flush=True)
        return None
    cT = rel[0]
    height = max(abs(x) for x in rel)
    terms = {n: rel[i + 1] for i, n in enumerate(names) if rel[i + 1] != 0}
    print(f"    [{label} dps={dps}] REL: target-coeff={cT} height={height}  {terms}", flush=True)
    return rel


def main():
    val = mp.mpf(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else 'W'
    trusted = int(sys.argv[3]) if len(sys.argv) > 3 else 18
    wmax = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    use_P8 = (sys.argv[5].lower() in ('1', 'true', 'yes')) if len(sys.argv) > 5 else True

    mp.mp.dps = trusted + 25
    if mode == 'V':
        W = val * mp.pi / 8
    else:
        W = val
    # decoy: structureless transcendental of COMPARABLE magnitude to W (~0.155),
    # algebraically unrelated to pi and Gamma(1/4). NOT normalized to |W| (that
    # would make it EQUAL W). Must NOT find the same low-height relation as W.
    decoy = mp.log(mp.mpf(11)) / mp.mpf('15.4')   # ~0.15571

    print(f"W = {mp.nstr(W, trusted)}  (mode={mode}, trusted~{trusted} digits, wmax={wmax}, P8={use_P8})", flush=True)
    print(f"decoy = {mp.nstr(decoy, trusted)}", flush=True)

    dps_grid = sorted(set([max(12, trusted - 6), max(14, trusted - 3), trusted]))
    for dps in dps_grid:
        basis = build_basis(dps, wmax, use_P8)
        n = len(basis)
        print(f"\n== dps={dps}  basis size n={n} ==", flush=True)
        # false-positive height scale ~ 10^(dps/(n-1))
        fp = mp.power(10, mp.mpf(dps) / (n - 1))
        print(f"   (spurious-relation height scale ~ 10^(dps/(n-1)) = {mp.nstr(fp,4)})", flush=True)
        guarded_pslq(-W, basis, dps, 10 ** 8, "REAL")
        guarded_pslq(-decoy, basis, dps, 10 ** 8, "DECOY")


if __name__ == '__main__':
    main()
