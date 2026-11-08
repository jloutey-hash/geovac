"""P3e-2W basis blocks and spans, exactly as fixed in docs/predictions/2026-10-19_p3e2_widened.md
(Model: Basis blocks, Spans, Threshold spans) and debug/p3e2_scoping_memo.md Sec. 3.

Letters: nonzero values of a, b, c, a+-b, a+-c, b+-c, a+b+c, -a+b+c, a-b+c, a+b-c; Pr = primes
dividing them.  Lambda = residue moduli k_ab, k_ac (|v-u|/(v+u), the R-block curves; the table's
'residue moduli' for T2 and T4) and the three threshold ratios, distinct, 0 < |l| < 1.
Bailey: c = max mass, p <= q the other two, k+- = [sqrt((c+p)^2-q^2) +- sqrt((c-p)^2-q^2)]/(2c),
kappa = pi/(2c), V = kappa [K(k-)K(k+') + K(k+)K(k-')] (checked numerically: assert).
"""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
import sympy as sp

DATA = Path(__file__).parent / "data"
THRESHOLD = {"S1", "S3"}


def letters(m):
    a, b, c = m
    vals = [a, b, c, a + b, a - b, a + c, a - c, b + c, b - c, a + b + c, -a + b + c, a - b + c, a + b - c]
    vals = [abs(v) for v in vals if v != 0]
    Pr = sorted({p for v in vals for p in sp.primefactors(v)})
    lam = []
    for u, v in [(a, b), (a, c)]:
        lam.append(Fraction(abs(v - u), v + u))
    s = a + b + c
    lam += [Fraction(-a + b + c, s), Fraction(a - b + c, s), Fraction(a + b - c, s)]
    L = []
    for l in lam:
        if 0 < abs(l) < 1 and l not in L:
            L.append(l)
    return Pr, L


def squarefree(n: int) -> int:
    f = sp.factorint(n)
    out = 1
    for p, e in f.items():
        if e % 2:
            out *= p
    return out


def blocks(ID: str, masses, vals: dict):
    """vals: dict of mpf grid values (current precision).  Returns ordered dict name -> value per block."""
    a, b, c = masses
    Pr, Lam = letters(masses)
    pi = mp.pi
    K = lambda k: mp.ellipk(k * k)
    E = lambda k: mp.ellipe(k * k)
    cm = max(masses)
    rest = sorted(masses)
    rest.remove(cm)
    p, q = rest
    A_ = (cm + p) ** 2 - q**2
    B_ = (cm - p) ** 2 - q**2
    kap = pi / (2 * cm)
    out = {}
    out["K3"] = {k: vals[k] for k in ("V", "D1", "D2", "S12")}
    FUN = {"K": K, "E": E}
    if ID in THRESHOLD:
        assert B_ == 0
        k = mp.sqrt(A_) / (2 * cm)
        kp = mp.sqrt(1 - k * k)
        out["DegB"] = {f"kap{X}{Y}'": kap * FUN[X](k) * FUN[Y](kp) for X in "KE" for Y in "KE"}
        bailey_V = 2 * kap * K(k) * K(kp)
    elif B_ > 0:
        kpl, kmi = (mp.sqrt(A_) + mp.sqrt(B_)) / (2 * cm), (mp.sqrt(A_) - mp.sqrt(B_)) / (2 * cm)
        kplp, kmip = mp.sqrt(1 - kpl**2), mp.sqrt(1 - kmi**2)
        d = squarefree(B_)
        Bb = {}
        for X in "KE":
            for Y in "KE":
                s_ = kap * (FUN[X](kmi) * FUN[Y](kplp) + FUN[X](kpl) * FUN[Y](kmip))
                a_ = kap * mp.sqrt(d) * (FUN[X](kmi) * FUN[Y](kplp) - FUN[X](kpl) * FUN[Y](kmip))
                if X + Y != "KK":
                    Bb[f"S_{X}{Y}"] = s_
                else:
                    bailey_V = s_
                Bb[f"A_{X}{Y}"] = a_
        out["B"] = {k_: Bb[k_] for k_ in ["S_KE", "S_EK", "S_EE", "A_KK", "A_KE", "A_EK", "A_EE"]}
    else:
        kpl = (mp.sqrt(A_) + mp.mpc(0, 1) * mp.sqrt(-B_)) / (2 * cm)
        kmi = mp.conj(kpl)
        kmip = mp.sqrt(1 - kmi**2)
        Bb = {}
        for X in "KE":
            for Y in "KE":
                z = kap * FUN[X](kpl) * FUN[Y](kmip)
                if X + Y != "KK":
                    Bb[f"S_{X}{Y}"] = 2 * mp.re(z)
                else:
                    bailey_V = 2 * mp.re(z)
                Bb[f"A_{X}{Y}"] = 2 * mp.im(z)
        out["B"] = {k_: Bb[k_] for k_ in ["S_KE", "S_EK", "S_EE", "A_KK", "A_KE", "A_EK", "A_EE"]}
    out["_bailey_check"] = abs(bailey_V - vals["V"])
    out["R"] = {"piPab": pi * vals["Pab"], "piQab": pi * vals["Qab"], "piPac": pi * vals["Pac"], "piQac": pi * vals["Qac"]}
    P1, Q1, P2, Q2 = vals["Pab"], vals["Qab"], vals["Pac"], vals["Qac"]
    out["RR"] = {"PabPac/pi": P1 * P2 / pi, "PabQac/pi": P1 * Q2 / pi, "QabPac/pi": Q1 * P2 / pi, "QabQac/pi": Q1 * Q2 / pi}
    T3 = {"pi^3": pi**3, "zeta3": mp.zeta(3)}
    for pr in Pr:
        T3[f"pi^2log{pr}"] = pi**2 * mp.log(pr)
    out["T3"] = T3
    Li = {}
    for l in Lam:
        lv = mp.mpf(l.numerator) / l.denominator
        Li[f"piLi2({l})"] = pi * mp.polylog(2, lv)
        Li[f"Li3({l})"] = mp.polylog(3, lv)
    out["Li"] = Li
    out["O"] = {k: vals[k] for k in ("U", "Ua", "Ub", "Uc")}
    LR = {}
    for pr in Pr:
        for nm in ("Pab", "Qab", "Pac", "Qac"):
            LR[f"log{pr}*{nm}"] = mp.log(pr) * vals[nm]
    out["LR"] = LR
    out["pi3"] = {"pi^3": pi**3}
    out["_letters"] = {"Pr": Pr, "Lambda": [str(l) for l in Lam]}
    return out


SPANS = {"SA": ["K3", "B", "R", "RR", "pi3"], "SB": ["K3", "T3"], "SC": ["K3", "T3", "Li"],
         "SD": ["K3", "O", "T3", "R"], "SE": ["K3", "R", "LR", "pi3"]}
DEGREE = {"O": 2}   # every other block element has degree 3 under the file's counting


def span(bl: dict, name: str, ID: str):
    """Ordered list of (element name, block, value).  At threshold triples K3 and B -> DegB."""
    blist = SPANS[name]
    if ID in THRESHOLD:
        blist = [x for x in blist if x not in ("K3", "B")]
        blist = ["DegB"] + blist
    out, seen = [], set()
    for blk in blist:
        for k, v in bl[blk].items():
            if k in seen:
                continue
            seen.add(k)
            out.append((k, blk, v))
    return out


def load_values(ID: str, W: int = 520) -> dict:
    rec = json.loads((DATA / f"p3e2w_values_{ID}_{W}.json").read_text())
    return rec
