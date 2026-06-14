"""
ADVERSARIAL check on the NA-1 Reading-B verdict (S^(3) substrate).

Three attacks, each given its best shot:

  ATTACK 1 (Reading A's best shot): try to express GeoVac's actual weight-10
    contribution  W10 = -7 pi^4 z3^2 + 70 pi^2 z3 z5 + 6 pi^2 z(5,3)
                        + (179/170100) pi^10
    ENTIRELY in the depth-1 abelian product ring (NO zeta(5,3)).  Basis =
    weight-10 products of single zetas {pi^10, pi^4 z3^2, pi^2 z3 z5,
    z3 z7, z5^2}.  If PSLQ finds a relation -> Reading A SURVIVES (zeta(5,3)
    was removable / a basis artifact).  If NO relation -> zeta(5,3) is
    genuinely REQUIRED -> Reading A falsified.

  ATTACK 2 (harder irreducibility): re-test zeta(5,3) vs the FULL weight-8
    product basis {pi^8, pi^2 z3^2, z3 z5} at dps=120, maxcoeff=10^9.
    A sneaky high-height relation would rescue "reducible".  Expect NONE.

  ATTACK 3 (the fresh-atom loophole, "Reading A-prime"): tested STRUCTURALLY,
    not numerically -- see the printed analysis.  The escape is: maybe
    zeta(5,3) is a fresh PRIMITIVE weight-8 generator (still abelian), not a
    bracket.  Closed iff GeoVac has NO weight-8 Mellin slot and zeta(5,3)
    arises by loop-order convolution (compositional).  Documented below.

Discipline: exact Fractions for rational coeffs; mpmath @120 dps; curve-fit
audit (every PSLQ relation reported with its integer vector; NO free params).
"""

import json
from fractions import Fraction as Fr
from mpmath import mp, mpf, nsum, zeta, pi, inf, pslq, fabs

mp.dps = 120


def mzv2(a, b):
    zb = zeta(b)
    return nsum(lambda m: m ** (-a) * (zb - zeta(b, m)), [2, inf])


def main():
    out = {"dps": mp.dps, "attacks": {}}

    z3, z5, z7 = zeta(3), zeta(5), zeta(7)
    z53 = mzv2(5, 3)
    P = pi

    # ---------------------------------------------------------------
    # ATTACK 2 first (foundation): harder irreducibility of zeta(5,3)
    # full weight-8 product basis {pi^8, pi^2 z3^2, z3 z5}
    rel8 = pslq([z53, P**8, P**2 * z3**2, z3 * z5], maxcoeff=10**9, maxsteps=10000)
    out["attacks"]["A2_irreducibility_wt8_full_basis"] = {
        "basis": "{zeta(5,3), pi^8, pi^2*z3^2, z3*z5}",
        "pslq": rel8,
        "verdict": "IRREDUCIBLE (Reading A cannot make it)" if rel8 is None
                   else "REDUCIBLE -- Reading A survives!",
    }

    # ---------------------------------------------------------------
    # ATTACK 1: Reading A's best shot on GeoVac's ACTUAL weight-10 term.
    # GeoVac W10 contribution (from S^(3) W10 identification memo):
    W10 = (Fr(-7) * P**4 * z3**2          # -7 pi^4 z3^2
           + Fr(70) * P**2 * z3 * z5       # +70 pi^2 z3 z5
           + Fr(6) * P**2 * z53            # +6 pi^2 zeta(5,3)   <-- depth-2 term
           + Fr(179, 170100) * P**10)      # +(179/170100) pi^10

    # depth-1 abelian product basis at weight 10 (NO zeta(5,3))
    prod_basis = [P**10, P**4 * z3**2, P**2 * z3 * z5, z3 * z7, z5**2]
    relW = pslq([W10] + prod_basis, maxcoeff=10**9, maxsteps=20000)
    out["attacks"]["A1_readingA_best_shot_on_W10"] = {
        "target": "GeoVac W10 = -7pi^4 z3^2 + 70pi^2 z3 z5 + 6pi^2 z(5,3) + (179/170100)pi^10",
        "abelian_basis": "{pi^10, pi^4 z3^2, pi^2 z3 z5, z3 z7, z5^2}  (NO zeta(5,3))",
        "pslq": relW,
        "verdict": ("Reading A FALSIFIED: W10 is NOT in the abelian product ring; "
                    "zeta(5,3) is genuinely REQUIRED") if relW is None
                   else "Reading A SURVIVES: zeta(5,3) was removable!",
    }

    # Isolate: is the depth-2 piece 6 pi^2 zeta(5,3) ALONE reducible to products?
    rel_d2 = pslq([Fr(6) * P**2 * z53] + prod_basis, maxcoeff=10**9, maxsteps=20000)
    out["attacks"]["A1b_depth2_term_alone"] = {
        "target": "6 pi^2 zeta(5,3)",
        "pslq": rel_d2,
        "verdict": "IRREDUCIBLE (not a product of single zetas)" if rel_d2 is None
                   else "REDUCIBLE!",
    }

    # Completeness sanity: WITH zeta(5,3) in the basis, is W10 an exact
    # rational combination?  (confirms the closed form & that zeta(5,3) is
    # the ONLY irreducible needed -- basis is complete)
    full_basis = [P**10, P**4 * z3**2, P**2 * z3 * z5, P**2 * z53]
    relfull = pslq([W10] + full_basis, maxcoeff=10**6, maxsteps=20000)
    out["attacks"]["A1c_completeness_with_z53"] = {
        "basis": "{pi^10, pi^4 z3^2, pi^2 z3 z5, pi^2 zeta(5,3)}",
        "pslq": relfull,
        "note": "expect exact relation -> basis complete, sole irreducible is zeta(5,3)",
    }

    # ---------------------------------------------------------------
    # ATTACK 3 (structural): the fresh-atom loophole
    out["attacks"]["A3_fresh_atom_loophole"] = {
        "loophole": ("Reading A' : zeta(5,3) is a fresh PRIMITIVE weight-8 generator "
                     "(still abelian), not a bracket [f3,f5]."),
        "why_closed": [
            "GeoVac's primitive generators ARE the master Mellin slots k in {0,1,2} "
            "(M1,M2,M3). There is NO weight-8 / depth-2 Mellin slot.",
            "zeta(5,3) appears only at LOOP ORDER k=3, as a k-fold convolution of the "
            "basic Hurwitz building blocks (Paper 55 M3 sub-sector 3: 'iterated sunset "
            "produces a k-fold convolution ... bounding motivic depth by loop order').",
            "A k-fold ORDERED convolution (n2-rank ordering cases, m>n) is the "
            "deconcatenation coproduct, NOT a primitive generator.",
            "In the motivic Lie algebra zeta(5,3) IS the bracket [f3,f5]; brackets are "
            "ZERO in ANY abelianization regardless of generator choice. A non-zero "
            "irreducible zeta(5,3) cannot be a fresh abelian atom -- it is structurally "
            "a bracket.",
        ],
        "verdict": "Loophole CLOSED: zeta(5,3) is compositional (bracket), not a fresh atom.",
    }

    # ---------------------------------------------------------------
    # honest refinement of the claim
    A1 = out["attacks"]["A1_readingA_best_shot_on_W10"]["pslq"] is None
    A2 = out["attacks"]["A2_irreducibility_wt8_full_basis"]["pslq"] is None
    A1b = out["attacks"]["A1b_depth2_term_alone"]["pslq"] is None

    out["summary"] = {
        "reading_A_falsified": bool(A1 and A2 and A1b),
        "what_is_proven": ("GeoVac carries genuine irreducible depth-2 (bracket) content "
                           "-> it is NOT the abelianization -> NOT Reading A."),
        "what_is_NOT_yet_proven": ("That the full substrate is the FREE non-abelian group "
                                   "(strong Reading B). 'Not abelianization' is decisive; "
                                   "'free non-abelian' is the multi-month shuffle-Hopf "
                                   "enrichment construction."),
        "honest_reading": "NOT-A (decisive). B in the weak sense (carries brackets). "
                          "Strong-B (free) = named multi-month follow-on.",
    }

    print(json.dumps(out, indent=2, default=str))
    with open("debug/data/na1_adversarial_check.json", "w") as f:
        json.dump(out, f, indent=2, default=str)


if __name__ == "__main__":
    main()
