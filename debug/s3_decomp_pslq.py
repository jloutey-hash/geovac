"""S^(3) round-2 decomposition — Steps 3-5: stuffle reduction + targeted PSLQ.

Stage A (dps 220):
  R. Verify the 8 convergent stuffle relations used to symmetrize the
     12 triple t-values of S^(3) (t(a)t(b,c) = t(a,b,c)+t(b,a,c)+t(b,c,a)
     +t(a+b,c)+t(b,a+c), with t(a)=lam(a)):
       R1: lam2 t(2,2) = 3 t(2,2,2) + t(4,2) + t(2,4)
       R2: lam3 t(2,1) = t(3,2,1)+t(2,3,1)+t(2,1,3) + t(5,1) + t(2,4)
       R4: lam2 t(4,2) = t(2,4,2)+2t(4,2,2) + t(6,2) + t(4,4)
       R5: lam3 t(2,3) = t(3,2,3)+2t(2,3,3) + t(5,3) + t(2,6)
       R6: lam3 t(4,1) = t(3,4,1)+t(4,3,1)+t(4,1,3) + t(7,1) + t(4,4)
       R7: lam3 t(4,3) = t(3,4,3)+2t(4,3,3) + t(7,3) + t(4,6)
       R8: lam4 t(4,2) = 2t(4,4,2)+t(4,2,4) + t(8,2) + t(4,6)
     Post-stuffle depth-3 survivors (one per argument-multiset family):
       w4: t3(2,1,1); w6: t3(3,2,1), t3(4,1,1);
       w8: t3(3,2,3), t3(3,4,1), t3(4,2,2); w10: t3(3,4,3), t3(4,2,4).
  P. Targeted PSLQ inside the level-2 (MT(Z[1/2]) / alternating-MZV) room,
     weight-homogeneous graded bases (conjectural dims 1,1,2,3,5,8,13,21,34):
       w4 (dim 5):  t3(2,1,1); t2(3,1)
       w6 (dim 13 = 12 products + rep t2(5,1)): t3(4,1,1), t3(3,2,1),
            antisym a6 = t2(4,2)-t2(2,4); irreducibility probe t2(5,1)
       odd-weight doubles (depth-2 parity room, restricted product bases):
            w5: t2(4,1), t2(2,3);  w7: t2(4,3), t2(6,1), t2(2,5)
            w9: t2(6,3), t2(8,1), t2(2,7), t2(4,5);  w11: t2(8,3), t2(4,7)

Stage B (dps 340): w8 survivors t3(4,2,2), t3(3,2,3), t3(3,4,1) and the
     even doubles t2(6,2), t2(2,6) against the 30-element w8 basis
     (27 products + Li8(1/2) + depth-2 reps t2(7,1), t2(5,3)).

w10 (graded piece -2048 t3(4,3,3) - 1024 t3(4,4,2)): basis dim 89 > 40 cap
     -> SKIPPED per charter; recorded with required basis size.

Values cached in debug/data/s3_pslq_cache.json.
Results in debug/data/s3_pslq_stageA.json / s3_pslq_stageB.json.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath

HERE = Path(__file__).parent
CACHE_PATH = HERE / "data" / "s3_pslq_cache.json"
CACHE = json.loads(CACHE_PATH.read_text()) if CACHE_PATH.exists() else {}


def cache_get(key, dps, fn):
    k = "%s@%d" % (key, dps)
    if k in CACHE:
        return mpmath.mpf(CACHE[k])
    t0 = time.time()
    v = fn()
    CACHE[k] = mpmath.nstr(v, dps - 5)
    CACHE_PATH.write_text(json.dumps(CACHE, indent=1))
    print("    [%s computed, %.0fs]" % (k, time.time() - t0))
    return v


def lam(s):
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


def t2_raw(b, c):
    def term(j):
        j = int(j)
        return (mpmath.mpf(2 * j + 1) ** -c * mpmath.power(2, -b)
                * mpmath.zeta(b, mpmath.mpf(j) + mpmath.mpf(3) / 2))
    return mpmath.nsum(term, [0, mpmath.inf], method='levin')


def t3_raw(b1, b2, b3):
    lam3 = lam(b3) if b3 >= 2 else None
    psih = mpmath.digamma(mpmath.mpf(1) / 2)

    def term(j):
        j = int(j)
        tau = mpmath.power(2, -b1) * mpmath.zeta(
            b1, mpmath.mpf(j) + mpmath.mpf(3) / 2)
        if b3 == 1:
            pL = (mpmath.digamma(mpmath.mpf(j) + mpmath.mpf(1) / 2)
                  - psih) / 2
        else:
            pL = lam3 - mpmath.power(2, -b3) * mpmath.zeta(
                b3, mpmath.mpf(j) + mpmath.mpf(1) / 2)
        return mpmath.mpf(2 * j + 1) ** -b2 * tau * pL
    return mpmath.nsum(term, [1, mpmath.inf], method='levin')


def T2(b, c, dps):
    return cache_get("t2(%d,%d)" % (b, c), dps, lambda: t2_raw(b, c))


def T3(b1, b2, b3, dps):
    return cache_get("t3(%d,%d,%d)" % (b1, b2, b3), dps,
                     lambda: t3_raw(b1, b2, b3))


def monomials(gens, w):
    """All multiset monomials of total weight w over gens [(name,wt,val)]."""
    out = []

    def rec(i, rem, label, val):
        if rem == 0:
            out.append(("*".join(label) if label else "1", val))
            return
        if i >= len(gens):
            return
        name, wt, v = gens[i]
        kmax = rem // wt
        for k in range(kmax + 1):
            rec(i + 1, rem - k * wt, label + [name] * k, val * v ** k)

    rec(0, w, [], mpmath.mpf(1))
    return out


def identify(name, value, basis, dps, maxcoeff=10 ** 10):
    tol = mpmath.mpf(10) ** -(dps - 40)
    vec = [value] + [v for _, v in basis]
    rel = mpmath.pslq(vec, tol=tol, maxcoeff=maxcoeff, maxsteps=500000)
    if rel is None or rel[0] == 0:
        print("  %s: NO relation (basis dim %d, dps %d)"
              % (name, len(basis), dps))
        return None
    coeffs = {basis[i - 1][0]: Fraction(-rel[i], rel[0])
              for i in range(1, len(rel)) if rel[i] != 0}
    recon = sum(mpmath.mpf(f.numerator) / f.denominator * dict(basis)[nm]
                for nm, f in coeffs.items())
    res = abs(value - recon)
    ok = res < mpmath.mpf(10) ** -(dps - 60)
    print("  %s = %s" % (name, " + ".join(
        "(%s)*%s" % (f, nm) for nm, f in sorted(coeffs.items()))))
    print("    residual %s -> %s" % (mpmath.nstr(res, 4),
                                     "ACCEPT" if ok else "REJECT"))
    return {"coeffs": {nm: str(f) for nm, f in coeffs.items()},
            "residual": mpmath.nstr(res, 4)} if ok else None


# ===================================================================
def stage_A():
    dps = 220
    mpmath.mp.dps = dps
    OUT = {"dps": dps}
    ln2 = mpmath.log(2)
    pi = mpmath.pi
    z3, z5, z7, z9, z11 = (mpmath.zeta(k) for k in (3, 5, 7, 9, 11))
    Li = {k: mpmath.polylog(k, mpmath.mpf(1) / 2) for k in range(4, 8)}

    print("=== Stage A: stuffle relation verification (dps %d) ===" % dps)
    need_t3 = [(2, 1, 1), (2, 1, 3), (2, 2, 2), (2, 3, 1), (2, 3, 3),
               (2, 4, 2), (4, 1, 1), (4, 1, 3), (4, 2, 2), (4, 3, 1),
               (4, 3, 3), (4, 4, 2), (3, 2, 1), (3, 2, 3), (3, 4, 1),
               (3, 4, 3), (4, 2, 4)]
    t3v = {k: T3(*k, dps=dps) for k in need_t3}
    need_t2 = [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7),
               (3, 1), (3, 3), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5),
               (4, 6), (4, 7), (5, 1), (5, 3), (6, 1), (6, 2), (6, 3),
               (7, 1), (7, 3), (8, 1), (8, 2), (8, 3)]
    t2v = {k: T2(*k, dps=dps) for k in need_t2}

    rels = {
        "R1": lam(2) * t2v[(2, 2)] - (3 * t3v[(2, 2, 2)]
                                      + t2v[(4, 2)] + t2v[(2, 4)]),
        "R2": lam(3) * t2v[(2, 1)] - (t3v[(3, 2, 1)] + t3v[(2, 3, 1)]
                                      + t3v[(2, 1, 3)] + t2v[(5, 1)]
                                      + t2v[(2, 4)]),
        "R4": lam(2) * t2v[(4, 2)] - (t3v[(2, 4, 2)] + 2 * t3v[(4, 2, 2)]
                                      + t2v[(6, 2)] + t2v[(4, 4)]),
        "R5": lam(3) * t2v[(2, 3)] - (t3v[(3, 2, 3)] + 2 * t3v[(2, 3, 3)]
                                      + t2v[(5, 3)] + t2v[(2, 6)]),
        "R6": lam(3) * t2v[(4, 1)] - (t3v[(3, 4, 1)] + t3v[(4, 3, 1)]
                                      + t3v[(4, 1, 3)] + t2v[(7, 1)]
                                      + t2v[(4, 4)]),
        "R7": lam(3) * t2v[(4, 3)] - (t3v[(3, 4, 3)] + 2 * t3v[(4, 3, 3)]
                                      + t2v[(7, 3)] + t2v[(4, 6)]),
        "R8": lam(4) * t2v[(4, 2)] - (2 * t3v[(4, 4, 2)] + t3v[(4, 2, 4)]
                                      + t2v[(8, 2)] + t2v[(4, 6)]),
    }
    for k, v in rels.items():
        print("  %s residual: %s" % (k, mpmath.nstr(abs(v), 4)))
    OUT["stuffle_residuals"] = {k: mpmath.nstr(abs(v), 4)
                                for k, v in rels.items()}
    OUT["stuffle_all_pass"] = all(
        abs(v) < mpmath.mpf(10) ** -(dps - 40) for v in rels.values())
    print("  stuffle all pass:", OUT["stuffle_all_pass"])

    # ---------------- PSLQ: w4 ----------------
    print("\n=== Stage A: PSLQ w4 (dim 5) ===")
    gens4 = [("ln2", 1, ln2), ("pi2", 2, pi ** 2), ("z3", 3, z3),
             ("Li4", 4, Li[4])]
    b4 = monomials(gens4, 4)
    OUT["t3(2,1,1)"] = identify("t3(2,1,1)", t3v[(2, 1, 1)], b4, dps)
    OUT["t2(3,1)"] = identify("t2(3,1)", t2v[(3, 1)], b4, dps)

    # ---------------- PSLQ: w6 ----------------
    print("\n=== Stage A: PSLQ w6 (dim 13) ===")
    gens6 = gens4 + [("z5", 5, z5), ("Li5", 5, Li[5]), ("Li6", 6, Li[6]),
                     ("t2(5,1)", 6, t2v[(5, 1)])]
    b6 = monomials(gens6, 6)
    print("  w6 basis dim = %d (conjectural 13)" % len(b6))
    OUT["w6_dim"] = len(b6)
    OUT["t3(4,1,1)"] = identify("t3(4,1,1)", t3v[(4, 1, 1)], b6, dps)
    OUT["t3(3,2,1)"] = identify("t3(3,2,1)", t3v[(3, 2, 1)], b6, dps)
    a6 = t2v[(4, 2)] - t2v[(2, 4)]
    OUT["a6=t2(4,2)-t2(2,4)"] = identify("t2(4,2)-t2(2,4)", a6, b6, dps)
    # irreducibility probe for the rep itself (against products only)
    b6_prod = [x for x in b6 if x[0] != "t2(5,1)"]
    OUT["t2(5,1)_probe"] = identify("t2(5,1) [products only]",
                                    t2v[(5, 1)], b6_prod, dps)

    # ---------------- PSLQ: odd-weight doubles ----------------
    print("\n=== Stage A: odd-weight doubles ===")
    gens5 = gens4 + [("z5", 5, z5), ("Li5", 5, Li[5])]
    b5 = monomials(gens5, 5)
    print("  w5 basis dim = %d (conjectural 8)" % len(b5))
    OUT["t2(4,1)"] = identify("t2(4,1)", t2v[(4, 1)], b5, dps)
    OUT["t2(2,3)"] = identify("t2(2,3)", t2v[(2, 3)], b5, dps)

    gens7 = gens6[:-1] + [("z7", 7, z7), ("Li7", 7, Li[7]),
                          ("t2(5,1)", 6, t2v[(5, 1)])]
    b7 = monomials(gens7, 7)
    print("  w7 basis dim = %d (19 products incl t2(5,1)*ln2)" % len(b7))
    for key in ((4, 3), (6, 1), (2, 5)):
        OUT["t2%s" % (key,)] = identify("t2%s" % (key,), t2v[key], b7, dps)

    b9 = [("z9", z9), ("pi2*z7", pi ** 2 * z7), ("pi4*z5", pi ** 4 * z5),
          ("pi6*z3", pi ** 6 * z3), ("pi8*ln2", pi ** 8 * ln2)]
    print("  w9 restricted basis dim = 5 (depth-1 parity room)")
    for key in ((6, 3), (8, 1), (2, 7), (4, 5)):
        OUT["t2%s" % (key,)] = identify("t2%s" % (key,), t2v[key], b9, dps)

    b11 = [("z11", z11), ("pi2*z9", pi ** 2 * z9), ("pi4*z7", pi ** 4 * z7),
           ("pi6*z5", pi ** 6 * z5), ("pi8*z3", pi ** 8 * z3),
           ("pi10*ln2", pi ** 10 * ln2)]
    print("  w11 restricted basis dim = 6")
    for key in ((8, 3), (4, 7)):
        OUT["t2%s" % (key,)] = identify("t2%s" % (key,), t2v[key], b11, dps)

    OUT["w10_skipped"] = ("graded piece -2048 t3(4,3,3) - 1024 t3(4,4,2): "
                          "level-2 w10 basis dim 89 > 40 cap; "
                          "stuffle pins t3(3,4,3) + 2 t3(4,3,3) and "
                          "2 t3(4,4,2) + t3(4,2,4) (R7, R8); parity theorem "
                          "(level 2, w+d = 13 odd) predicts depth<=2 "
                          "reducibility; NOT executed.")
    (HERE / "data" / "s3_pslq_stageA.json").write_text(
        json.dumps(OUT, indent=2))
    print("Saved stage A.")


# ===================================================================
def stage_B():
    dps = 340
    mpmath.mp.dps = dps
    OUT = {"dps": dps}
    ln2 = mpmath.log(2)
    pi = mpmath.pi
    z3, z5, z7 = mpmath.zeta(3), mpmath.zeta(5), mpmath.zeta(7)
    Li = {k: mpmath.polylog(k, mpmath.mpf(1) / 2) for k in range(4, 9)}

    print("=== Stage B: w8 room (dps %d) ===" % dps)
    t251 = T2(5, 1, dps)
    t271 = T2(7, 1, dps)
    t253 = T2(5, 3, dps)
    gens8 = [("ln2", 1, ln2), ("pi2", 2, pi ** 2), ("z3", 3, z3),
             ("Li4", 4, Li[4]), ("z5", 5, z5), ("Li5", 5, Li[5]),
             ("Li6", 6, Li[6]), ("t2(5,1)", 6, t251),
             ("z7", 7, z7), ("Li7", 7, Li[7]),
             ("Li8", 8, Li[8]), ("t2(7,1)", 8, t271),
             ("t2(5,3)", 8, t253)]
    b8 = monomials(gens8, 8)
    print("  w8 basis dim = %d (27 products + 3 w8 gens; conj. full 34 "
          "incl. depth>2 gens not needed for parity-reduced targets)"
          % len(b8))
    OUT["w8_dim"] = len(b8)

    for args in ((4, 2, 2), (3, 2, 3), (3, 4, 1)):
        v = T3(*args, dps=dps)
        OUT["t3%s" % (args,)] = identify("t3%s" % (args,), v, b8, dps)

    for key in ((6, 2), (2, 6)):
        v = T2(*key, dps=dps)
        basis = [x for x in b8 if x[0] != "t2%s" % (key,)]
        OUT["t2%s" % (key,)] = identify("t2%s" % (key,), v, basis, dps)

    (HERE / "data" / "s3_pslq_stageB.json").write_text(
        json.dumps(OUT, indent=2))
    print("Saved stage B.")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "A"
    (stage_A if stage.upper() == "A" else stage_B)()
