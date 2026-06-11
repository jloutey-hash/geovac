"""S^(3) round-2 decomposition — Step 2: exact decomposition engine.

Decomposes the o-space cores of S^(3) into Hoffman multiple t-values.

ORDERING CONVENTION (stated per charter):
    t3(b1,b2,b3) = sum_{o3 > o2 > o1 >= 1, odd} o3^(-b1) o2^(-b2) o1^(-b3)
    t2(b1,b2)    = sum_{o2 > o1 >= 1, odd}      o2^(-b1) o1^(-b2)
    lam(s)       = sum_{o >= 1 odd} o^(-s) = (1 - 2^(-s)) zeta(s)
i.e. the FIRST argument is carried by the LARGEST summation variable.

Inputs (verified in s3_decomp_setup.py, all bit-exact at N=12):
    S3 = C - 2G - 16 S_min P - 8 P^3 + 8 P Q + R3
    C  = sum_{o_i >= 5 odd} min(o1,o2) min(o2,o3) psi(o1) psi(o2) psi(o3)
    G  = sum_{o_i >= 5 odd} min(o1,o2) psi(o1) psi(o2)^2
    psi(o) = 8(o^2-1)/o^4,  psi(1) = 0,  psi(3) = 64/81.

This driver:
  E1. Builds the exact rational coefficient tables decomposing C and G
      (and the auxiliary Mfull2) over {t3, t2, lam, lam*lam, const}.
      Region split for C over o >= 1 (13 regions: 6 strict + 6 single-tie
      + diagonal), then inclusion-exclusion in the o_i = 3 pins
      (psi(1) = 0 makes the o=1 slices vanish identically).
  E2. Verifies the decomposition BIT-EXACTLY in Fraction arithmetic at
      truncation o <= O for O in (25, 41): every symbol truncated to the
      same box; identity holds exactly at every cutoff.
  E3. Emits the coefficient table debug/data/s3_decomp_table.json.

Region algebra for C3 := C extended to o >= 1 (derivation in memo):
  C3 = 4*S_Ra + 2*S_Rc + 3*S_T12a + 2*S_T12b + S_T13b + S_T123
  S_Ra   (strict, weight o_low*o_mid; x4: R_a,R_b,R_e,R_f)
       = 512 sum_{E in {2+,4-}, m in {1+,3-}, l in {1+,3-}} s t3(E,m,l)
  S_Rc   (strict, weight o_low^2 = o2^2, o2 smallest; x2: R_c,R_d)
       = 512 sum_{E in {2+,4-}, f in {2+,4-}} s [ (1/2)(t2(E,f-1)-t2(E,f))
                                                  - t3(E,f,2) ]
  S_T12a (tie low pair, weight o_low^2; x3: T_12a,T_23a,T_13a)
       = 512 sum_{E in {2+,4-}} s [t2(E,2) - 2 t2(E,4) + t2(E,6)]
  S_T12b (tie high pair, weight o_low*o_high; x2: T_12b,T_23b)
       = 512 sum_{B in {3+,5(-2),7+}, e in {1+,3-}} c s t2(B,e)
  S_T13b (tie outer pair o1=o3 large, weight o2^2; x1)
       = 512 sum_{B in {4+,6(-2),8+}} c [ (1/2)(lam(B-1)-lam(B)) - t2(B,2) ]
  S_T123 = 512 (lam4 - 3 lam6 + 3 lam8 - lam10)

  C  = C3 - 6 psi3 Mfull2 - 9 psi3 pt^2 + 27 psi3^2 pt - 9 psi3^3,
       pt = 8 lam2 - 8 lam4,
       Mfull2 = 128 [t2(2,1)-t2(2,3)-t2(4,1)+t2(4,3)] + 64 (lam3-2lam5+lam7)

  Gfull = 512 [ t2(4,1)-t2(4,3)-2t2(6,1)+2t2(6,3)+t2(8,1)-t2(8,3)
              + t2(2,3)-2t2(2,5)+t2(2,7)-t2(4,3)+2t2(4,5)-t2(4,7) ]
        + 512 (lam5 - 3 lam7 + 3 lam9 - lam11)
  G  = Gfull - 3 psi3 qt - 3 psi3^2 pt + 3 psi3^3,
       qt = 64 (lam4 - 2 lam6 + lam8)

Output: debug/data/s3_decomp_table.json
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

PSI3 = Fraction(64, 81)


def psi(o: int) -> Fraction:
    return Fraction(8 * (o * o - 1), o ** 4)


# ---------------------------------------------------------------------------
# Symbolic linear combinations: dict symbol -> Fraction
#   symbols: ('t3', b1, b2, b3) ('t2', b, c) ('lam', s)
#            ('lam2', s1, s2)  [product lam(s1)*lam(s2)]   ('one',)
# ---------------------------------------------------------------------------

def lc_add(*lcs):
    out = {}
    for lc in lcs:
        for k, v in lc.items():
            out[k] = out.get(k, Fraction(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def lc_scale(lc, s):
    return {k: v * s for k, v in lc.items() if v * s != 0}


def build_C3() -> dict:
    E_ = [(2, 1), (4, -1)]      # psi exponent/sign on a single variable
    M_ = [(1, 1), (3, -1)]      # o*psi
    lc = {}

    def add(sym, c):
        lc[sym] = lc.get(sym, Fraction(0)) + Fraction(c)

    # 4 * S_Ra
    for E, sE in E_:
        for m, sm in M_:
            for l, sl in M_:
                add(('t3', E, m, l), 4 * 512 * sE * sm * sl)
    # 2 * S_Rc
    for E, sE in E_:
        for f, sf in E_:
            s = sE * sf
            add(('t3', E, f, 2), -2 * 512 * s)
            add(('t2', E, f - 1), 2 * 512 * s * Fraction(1, 2))
            add(('t2', E, f), -2 * 512 * s * Fraction(1, 2))
    # 3 * S_T12a
    for E, sE in E_:
        for e, c in [(2, 1), (4, -2), (6, 1)]:
            add(('t2', E, e), 3 * 512 * sE * c)
    # 2 * S_T12b
    for B, c in [(3, 1), (5, -2), (7, 1)]:
        for e, se in M_:
            add(('t2', B, e), 2 * 512 * c * se)
    # S_T13b
    for B, c in [(4, 1), (6, -2), (8, 1)]:
        add(('lam', B - 1), 512 * c * Fraction(1, 2))
        add(('lam', B), -512 * c * Fraction(1, 2))
        add(('t2', B, 2), -512 * c)
    # S_T123
    for s, c in [(4, 1), (6, -3), (8, 3), (10, -1)]:
        add(('lam', s), 512 * c)
    return {k: v for k, v in lc.items() if v != 0}


def build_Mfull2() -> dict:
    lc = {}

    def add(sym, c):
        lc[sym] = lc.get(sym, Fraction(0)) + Fraction(c)

    for E, sE in [(2, 1), (4, -1)]:
        for e, se in [(1, 1), (3, -1)]:
            add(('t2', E, e), 128 * sE * se)
    for s, c in [(3, 1), (5, -2), (7, 1)]:
        add(('lam', s), 64 * c)
    return lc


def build_Gfull() -> dict:
    lc = {}

    def add(sym, c):
        lc[sym] = lc.get(sym, Fraction(0)) + Fraction(c)

    # o1 < o2 : small o*psi {1+,3-}, big psi^2 {4:+1,6:-2,8:+1}
    for B, c in [(4, 1), (6, -2), (8, 1)]:
        for e, se in [(1, 1), (3, -1)]:
            add(('t2', B, e), 512 * c * se)
    # o1 > o2 : big psi {2+,4-}, small o*psi^2 {3:+1,5:-2,7:+1}
    for E, sE in [(2, 1), (4, -1)]:
        for e, c in [(3, 1), (5, -2), (7, 1)]:
            add(('t2', E, e), 512 * sE * c)
    # diagonal o*psi^3
    for s, c in [(5, 1), (7, -3), (9, 3), (11, -1)]:
        add(('lam', s), 512 * c)
    return lc


def build_C() -> dict:
    """C over o >= 5, fully linear in {t3,t2,lam,lam2,one}."""
    pt = {('lam', 2): Fraction(8), ('lam', 4): Fraction(-8)}
    pt2 = {('lam2', 2, 2): Fraction(64), ('lam2', 2, 4): Fraction(-128),
           ('lam2', 4, 4): Fraction(64)}
    return lc_add(
        build_C3(),
        lc_scale(build_Mfull2(), -6 * PSI3),
        lc_scale(pt2, -9 * PSI3),
        lc_scale(pt, 27 * PSI3 ** 2),
        {('one',): -9 * PSI3 ** 3},
    )


def build_G() -> dict:
    pt = {('lam', 2): Fraction(8), ('lam', 4): Fraction(-8)}
    qt = {('lam', 4): Fraction(64), ('lam', 6): Fraction(-128),
          ('lam', 8): Fraction(64)}
    return lc_add(
        build_Gfull(),
        lc_scale(qt, -3 * PSI3),
        lc_scale(pt, -3 * PSI3 ** 2),
        {('one',): 3 * PSI3 ** 3},
    )


# ---------------------------------------------------------------------------
# Exact truncated evaluation (box o <= O, odd)
# ---------------------------------------------------------------------------

def eval_sym_trunc(sym, O: int) -> Fraction:
    odds = list(range(1, O + 1, 2))
    if sym == ('one',):
        return Fraction(1)
    kind = sym[0]
    if kind == 'lam':
        s = sym[1]
        return sum(Fraction(1, o ** s) for o in odds)
    if kind == 'lam2':
        s1, s2 = sym[1], sym[2]
        a = sum(Fraction(1, o ** s1) for o in odds)
        b = sum(Fraction(1, o ** s2) for o in odds)
        return a * b
    if kind == 't2':
        b, c = sym[1], sym[2]
        tot = Fraction(0)
        # iterate small variable, accumulate
        acc = Fraction(0)  # sum over o1 < o2 of o1^-c
        prev = None
        for o2 in odds:
            if prev is not None:
                acc += Fraction(1, prev ** c)
            tot += Fraction(1, o2 ** b) * acc
            prev = o2
        return tot
    if kind == 't3':
        b1, b2, b3 = sym[1], sym[2], sym[3]
        tot = Fraction(0)
        acc1 = {}  # o2 -> sum_{o1<o2} o1^-b3 (built incrementally)
        s1 = Fraction(0)
        prev = None
        inner2 = Fraction(0)  # sum_{o2<o3} o2^-b2 * acc1(o2)
        for o3 in odds:
            if prev is not None:
                # extend: include o2 = prev into inner2; first update s1
                inner2 += Fraction(1, prev ** b2) * acc1[prev]
            tot += Fraction(1, o3 ** b1) * inner2
            # build acc1 for the *next* odd value (= current o3): needs
            # sum_{o1 < o3}:
            if prev is not None:
                s1 += Fraction(1, prev ** b3)
            acc1[o3] = s1
            prev = o3
        return tot
    raise ValueError(sym)


def eval_lc_trunc(lc, O: int) -> Fraction:
    return sum(c * eval_sym_trunc(sym, O) for sym, c in lc.items())


def direct_C_trunc(O: int) -> Fraction:
    odds = list(range(5, O + 1, 2))
    P = {o: psi(o) for o in odds}
    return sum(min(a, b) * min(b, c) * P[a] * P[b] * P[c]
               for a in odds for b in odds for c in odds)


def direct_G_trunc(O: int) -> Fraction:
    odds = list(range(5, O + 1, 2))
    P = {o: psi(o) for o in odds}
    return sum(min(a, b) * P[a] * P[b] ** 2 for a in odds for b in odds)


def main():
    t0 = time.time()
    OUT = {}
    C_lc, G_lc = build_C(), build_G()

    print("E2: bit-exact truncated verification ...")
    ok_all = True
    for O in (25, 41):
        dC, eC = direct_C_trunc(O), eval_lc_trunc(C_lc, O)
        dG, eG = direct_G_trunc(O), eval_lc_trunc(G_lc, O)
        okC, okG = dC == eC, dG == eG
        ok_all &= okC and okG
        print("  O=%2d: C %s   G %s" % (O, "EXACT" if okC else "FAIL",
                                        "EXACT" if okG else "FAIL"))
        OUT["C_exact_O%d" % O] = okC
        OUT["G_exact_O%d" % O] = okG

    # full S3 combination: S3 = C - 2G - 16 S_min P - 8P^3 + 8PQ + R3
    # box-consistency of the product terms checked in setup V3; here record
    # the final LINEAR table for C - 2G:
    S3lin = lc_add(C_lc, lc_scale(G_lc, -2))
    OUT["all_pass"] = ok_all
    print("ALL EXACT: %s  (%.1fs)" % (ok_all, time.time() - t0))

    def serialize(lc):
        return {repr(k): [str(v.numerator), str(v.denominator)]
                for k, v in sorted(lc.items(), key=lambda kv: repr(kv[0]))}

    tab = {
        "convention": "t3(b1,b2,b3)=sum_{o3>o2>o1>=1 odd} "
                      "o3^-b1 o2^-b2 o1^-b3; first arg = largest variable",
        "S3_full_relation": "S3 = [C - 2G linear table below] "
                            "- 16*S_min*P - 8*P^3 + 8*P*Q + R3",
        "P": "pi^2 - pi^4/12 - 64/81",
        "Q": "4 zeta(4,5/2) - 2 zeta(6,5/2) + (1/4) zeta(8,5/2)",
        "R3": "8 zeta(6,5/2) - 6 zeta(8,5/2) + (3/2) zeta(10,5/2)"
              " - (1/8) zeta(12,5/2)",
        "S_min": "round-1 constant 2.47993693803422255441...",
        "C_table": serialize(C_lc),
        "G_table": serialize(G_lc),
        "S3_linear_table": serialize(S3lin),
        "verification": OUT,
    }
    outp = Path(__file__).parent / "data" / "s3_decomp_table.json"
    outp.write_text(json.dumps(tab, indent=2))
    print("Saved:", outp)

    # print the headline triple-t content
    print("\nTriple t-values in S3 (= those of C; G has none):")
    for k, v in sorted(C_lc.items()):
        if k[0] == 't3':
            w = sum(k[1:])
            print("  %+6d * t%s   [weight %d]" % (v, k[1:], w))


if __name__ == "__main__":
    main()
