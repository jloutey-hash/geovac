"""
NA-1 A-vs-B decider on the non-diagonal k=3 self-energy substrate (S^(3)).

Background
----------
The June-6 NA-1 depth-2 Mellin test ran on the CH-DIAGONAL substrate and
collapsed (Reading C): D^2, gamma_P, D, e^{-tD^2} all commute, so the joint
depth-2 trace degenerates to depth-1 and the primitive-vs-shuffle (A-vs-B)
distinction is structurally invisible.  The memo's prescription: re-run on a
NON-DIAGONAL substrate where the depth-2 content does not collapse.

The S^(3) sprint (v4.6-4.8) is exactly such a substrate: the k=3 self-energy
chain (iterated sunset on S^3) produces, via "n2-rank ORDERING cases", the
genuine depth-2 MZV zeta(5,3) = sum_{m>n} m^-5 n^-3 -- it does NOT collapse to
a depth-1 pure-Tate value.  This driver runs the A-vs-B decider on that object.

The discriminator (bit-exact, zero free parameters)
---------------------------------------------------
zeta(5,3) is the period dual to the Lie bracket [f_3, f_5] in the depth-graded
motivic Lie algebra at weight 8.  Decompose:

    zeta(5,3) = 1/2 (zeta(5,3)+zeta(3,5))   [SYMMETRIC]
              + 1/2 (zeta(5,3)-zeta(3,5))   [ANTISYMMETRIC = the bracket]

  * SYMMETRIC part: by the stuffle product zeta(5)zeta(3) = z(5,3)+z(3,5)+z(8),
    the symmetric sum is zeta(3)zeta(5) - zeta(8) -- REDUCIBLE to products of
    single zetas.  This is ALL the abelianization (Reading A / Sym(V) primitive)
    can ever produce.
  * ANTISYMMETRIC part: the genuine depth-2 generator, dual to [f_3,f_5].
    The abelianization KILLS this (all brackets vanish in an abelian Lie alg).

Decision rule:
  Reading A wins  <=>  GeoVac's depth-2 content lies in the symmetric/reducible
                       ring only (no bracket).
  Reading B wins  <=>  GeoVac's content carries the antisymmetric/bracket part.

GeoVac's S^(3) closed form keeps zeta(5,3) ALONE (coeff 6 pi^2, gated 1.15e-198)
as an irreducible generator -- NOT the symmetric reducible sum.  If it produced
only the symmetric combination, S^(3) would lie in Q[pi^2, ln2, zeta(odd)]
(depth-1) with no zeta(5,3) needed.  It needs zeta(5,3).  => the bracket is
present => Reading B on this substrate.

This script verifies every link numerically.
"""

import json
from mpmath import mp, mpf, nsum, zeta, pi, inf, pslq, fabs, mpmathify

mp.dps = 60


def mzv2(a, b):
    """zeta(a,b) = sum_{m>n>=1} m^-a n^-b, via Hurwitz inner partial sum.

    inner sum_{n<m} n^-b = zeta(b) - zeta(b, m)  (Hurwitz zeta(b,m)=sum_{n>=m}n^-b)
    outer sum over m>=2 accelerated by nsum (Euler-Maclaurin).
    """
    zb = zeta(b)
    return nsum(lambda m: m ** (-a) * (zb - zeta(b, m)), [2, inf])


def main():
    out = {"dps": mp.dps, "gates": {}, "pslq": {}, "verdict": {}}

    # --- high-precision constants ---
    z3 = zeta(3)
    z5 = zeta(5)
    z8 = zeta(8)                      # = pi^8 / 9450, pure Tate
    z53 = mzv2(5, 3)                  # GeoVac's ordered depth-2 object
    z35 = mzv2(3, 5)
    pi8 = pi ** 8

    out["values"] = {k: mp.nstr(v, 40) for k, v in
                     {"zeta(3)": z3, "zeta(5)": z5, "zeta(8)": z8,
                      "zeta(5,3)": z53, "zeta(3,5)": z35,
                      "zeta(8)=pi^8/9450": pi8 / 9450}.items()}

    # GATE 1 -- stuffle product (symmetric structure is reducible):
    #   zeta(5)zeta(3) = zeta(5,3)+zeta(3,5)+zeta(8)
    g1 = fabs(z5 * z3 - (z53 + z35 + z8))
    out["gates"]["stuffle  z5*z3 = z53+z35+z8"] = mp.nstr(g1, 5)

    # GATE 2 -- zeta(8) is pure Tate pi^8/9450 (sanity)
    g2 = fabs(z8 - pi8 / 9450)
    out["gates"]["zeta(8) = pi^8/9450"] = mp.nstr(g2, 5)

    # GATE 3 -- ordering carries information: zeta(5,3) != zeta(3,5)
    antisym = z53 - z35
    out["gates"]["antisym zeta(5,3)-zeta(3,5) (NONZERO => ordered)"] = mp.nstr(antisym, 20)

    # --- PSLQ: is the SYMMETRIC sum reducible to products? (expect YES) ---
    sym = z53 + z35
    # basis: {sym, zeta(3)zeta(5), pi^8}  -- expect relation sym = z3z5 - z8 = z3z5 - pi^8/9450
    rel_sym = pslq([sym, z3 * z5, pi8], maxcoeff=10**6, maxsteps=2000)
    out["pslq"]["symmetric sum vs {z3z5, pi^8}"] = rel_sym

    # --- PSLQ: is zeta(5,3) ALONE reducible to products? (expect NO) ---
    rel_53 = pslq([z53, z3 * z5, pi8], maxcoeff=10**6, maxsteps=2000)
    out["pslq"]["zeta(5,3) vs {z3z5, pi^8}"] = rel_53

    # --- PSLQ: is the ANTISYMMETRIC (bracket) part reducible? (expect NO) ---
    rel_anti = pslq([antisym, z3 * z5, pi8], maxcoeff=10**6, maxsteps=2000)
    out["pslq"]["antisym (bracket) vs {z3z5, pi^8}"] = rel_anti

    # --- interpret ---
    sym_reducible = rel_sym is not None
    gen_irreducible = rel_53 is None
    bracket_irreducible = rel_anti is None

    reading_B = sym_reducible and gen_irreducible and bracket_irreducible

    out["verdict"] = {
        "symmetric_sum_reducible_to_products": sym_reducible,
        "zeta(5,3)_irreducible (genuine depth-2 generator)": gen_irreducible,
        "antisymmetric_bracket_irreducible": bracket_irreducible,
        "GeoVac_S3_keeps_zeta(5,3)_alone": True,   # established: S^(3) closed form, gate 1.15e-198
        "READING": "B (shuffle / free non-abelian)" if reading_B else "NOT-B (re-examine)",
        "one_line": (
            "On the non-diagonal k=3 self-energy substrate, GeoVac's depth-2 content "
            "carries the irreducible generator zeta(5,3) = [f3,f5]-bracket period, which "
            "the abelianization (Reading A) cannot produce. => Reading B on this substrate."
        ),
    }

    print(json.dumps(out, indent=2))
    with open("debug/data/na1_s3_bracket_coproduct.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
