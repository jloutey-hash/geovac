"""VV stage A driver: the physical vector-photon vertex on round unit S^3.

Executes docs/predictions/2026-11-07_physical_vector_vertex.md (committed
974e3d45 before this driver existed).  Route 1 = group theory
(debug/vector_vertex_s3_route1.py); route 2 = R^4 coordinates
(debug/vector_vertex_s3_route2.py, imports neither route 1 nor geovac).
geovac.qed_anomalous_moment is imported ONLY for control C3.

Writes debug/data/vector_vertex_s3.json.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from fractions import Fraction
from typing import Dict, List, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import vector_vertex_s3_route1 as R1  # noqa: E402

SIGNS = (1, -1)
COMPS = ("A", "B")
N1MAX, Q1MAX = 3, 7
N2MAX, Q2MAX = 2, 5
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def sgn(s: int) -> str:
    return "+" if s > 0 else "-"


# ---------------------------------------------------------------------------
# the record's rule, transcribed literally
# ---------------------------------------------------------------------------
def rule_record(s1: int, s2: int, comp: str, n1: int, n2: int, q: int) -> bool:
    if q < 1:
        return False
    if s1 == s2:
        if (n1 + n2 + q) % 2 != 1:
            return False
        c = comp if s1 > 0 else {"A": "B", "B": "A"}[comp]   # -- is the mirror
        if c == "A":
            return abs(n1 - n2) + 1 <= q <= n1 + n2 + 1
        return abs(n1 - n2) + 1 <= q <= n1 + n2 - 1
    # opposite sign: stated for (+, -); (-, +) is its transpose n1 <-> n2
    if (n1 + n2 + q) % 2 != 0:
        return False
    if s1 < 0:
        n1, n2 = n2, n1
    if comp == "A":
        return (abs(n1 + 1 - n2) <= q + 1 <= n1 + n2 + 1) and (abs(n1 - n2 - 1) <= q - 1 <= n1 + n2 + 1)
    return (abs(n1 + 1 - n2) <= q - 1 <= n1 + n2 + 1) and (abs(n1 - n2 - 1) <= q + 1 <= n1 + n2 + 1)


def content_first_second(s: int, n: int, convention: int = 1) -> Tuple[int, int]:
    """Doubled (first, second) content; convention -1 = the mirror labelling."""
    a, b = (n + 1, n) if s * convention > 0 else (n, n + 1)
    return a, b


def comp_content(comp: str, q: int, convention: int = 1) -> Tuple[int, int]:
    A = (q + 1, q - 1)
    B = (q - 1, q + 1)
    if convention < 0:
        A, B = B, A
    return A if comp == "A" else B


def tri(a: int, b: int, c: int) -> bool:
    return abs(a - b) <= c <= a + b and (a + b + c) % 2 == 0


def rule_reptheory(s1, s2, comp, n1, n2, q, convention=1) -> bool:
    if q < 1:
        return False
    f1, g1 = content_first_second(s1, n1, convention)
    f2, g2 = content_first_second(s2, n2, convention)
    fa, ga = comp_content(comp, q, convention)
    return tri(f1, fa, f2) and tri(g1, ga, g2)


def model_rule(n1, n2, q) -> bool:
    return q >= 1 and abs(n1 - n2) <= q <= n1 + n2 and (n1 + n2 + q) % 2 == 1


def key(s1, s2, comp, n1, n2, q) -> str:
    return f"{sgn(s1)}{sgn(s2)}|{comp}|{n1},{n2},{q}"


OUT: Dict = {"record": "docs/predictions/2026-11-07_physical_vector_vertex.md",
             "record_commit": "974e3d45"}

# ---------------------------------------------------------------------------
# rule bookkeeping: record rule vs rep theory, mirror invariance
# ---------------------------------------------------------------------------
log("rule consistency")
rule_checks = {"record_vs_reptheory_mismatch": [], "mirror_mismatch": []}
for s1 in SIGNS:
    for s2 in SIGNS:
        for comp in COMPS:
            for n1 in range(N1MAX + 1):
                for n2 in range(N1MAX + 1):
                    for q in range(1, Q1MAX + 1):
                        a = rule_record(s1, s2, comp, n1, n2, q)
                        b = rule_reptheory(s1, s2, comp, n1, n2, q)
                        if a != b:
                            rule_checks["record_vs_reptheory_mismatch"].append(key(s1, s2, comp, n1, n2, q))
                        # mirror: swap which SU(2) is "first": sign -> -sign, A <-> B
                        m = rule_record(-s1, -s2, {"A": "B", "B": "A"}[comp], n1, n2, q)
                        if a != m:
                            rule_checks["mirror_mismatch"].append(key(s1, s2, comp, n1, n2, q))
                        c = rule_reptheory(s1, s2, comp, n1, n2, q, convention=-1)
                        d = rule_record(-s1, -s2, {"A": "B", "B": "A"}[comp], n1, n2, q)
                        if c != d:
                            rule_checks["mirror_mismatch"].append("conv-1:" + key(s1, s2, comp, n1, n2, q))
OUT["rule_checks"] = rule_checks
log("rule checks", {k: len(v) for k, v in rule_checks.items()})

# ---------------------------------------------------------------------------
# C1 / C2 route 1
# ---------------------------------------------------------------------------
log("C1/C2 route 1")
c1_r1 = {}
for n in range(N1MAX + 1):
    for s in SIGNS:
        tJ, tj = R1.spinor_content(s, n)
        lam = Fraction(s * (2 * n + 3), 2)
        ok = R1.check_eigen(tj, 1, tJ, Fraction(4), Fraction(3, 2), lam)  # D = 2 sigma.R + 3/2 = 4 S.R + 3/2
        mult = (tJ + 1) * (tj + 1)
        c1_r1[f"{sgn(s)}{n}"] = {"eigen_ok": ok, "mult": mult, "mult_ok": mult == (n + 1) * (n + 2)}
# exhaustion: functions (j) x spin 1/2 = J = j +- 1/2 ; every level accounted for
c1_r1_pass = all(v["eigen_ok"] and v["mult_ok"] for v in c1_r1.values())
c2_r1 = {}
for q in range(1, Q1MAX + 1):
    for comp in COMPS:
        tJ, tj = R1.vector_content(comp, q)
        lam = Fraction(-(q + 1) if comp == "A" else (q + 1))
        ok = R1.check_eigen(tj, 2, tJ, Fraction(-2), Fraction(-2), lam)  # *d = -2 S.R - 2
        mult = (tJ + 1) * (tj + 1)
        c2_r1[f"{comp}{q}"] = {"stard_eigen_ok": ok, "hodge": (q + 1) ** 2, "mult": mult,
                               "mult_ok": mult == q * (q + 2)}
    # exact forms J = j: *d = 0
    tj = q + 1
    c2_r1[f"exact{q}"] = {"stard_zero_ok": R1.check_eigen(tj, 2, tj, Fraction(-2), Fraction(-2), Fraction(0))}
c2_r1_pass = all(all(x for x in v.values() if isinstance(x, bool)) for v in c2_r1.values())
OUT["C1_route1"] = {"pass": c1_r1_pass, "detail": c1_r1}
OUT["C2_route1"] = {"pass": c2_r1_pass, "detail": c2_r1,
                    "coclosed": "by Schur: delta maps (j, J) to functions (j, j); J = j +- 1 has no image"}
log("C1 r1", c1_r1_pass, "C2 r1", c2_r1_pass)

# ---------------------------------------------------------------------------
# route 1 vertex table
# ---------------------------------------------------------------------------
log("route 1 vertex table")
R1TAB: Dict[str, Dict] = {}
for s1 in SIGNS:
    for s2 in SIGNS:
        for comp in COMPS:
            for n1 in range(N1MAX + 1):
                for n2 in range(N1MAX + 1):
                    for q in range(1, Q1MAX + 1):
                        R1TAB[key(s1, s2, comp, n1, n2, q)] = R1.vertex_route1(s1, n1, comp, q, s2, n2)
log("route 1 done", sum(v["open"] for v in R1TAB.values()), "open of", len(R1TAB))
mcons_all = all(v.get("m_conserved", True) for v in R1TAB.values())


# closed form: Sum|V|^2 = (2jA+1)(2j2+1) * 6 (2j1+1)(2JA+1)(2J2+1)(2J1+1) {9j}^2
def S_9j(s1, n1, comp, q, s2, n2) -> Fraction:
    from sympy.physics.wigner import wigner_9j
    from sympy import Rational, nsimplify
    tJ1, tj1 = R1.spinor_content(s1, n1)
    tJA, tjA = R1.vector_content(comp, q)
    tJ2, tj2 = R1.spinor_content(s2, n2)
    h = lambda x: Rational(x, 2)
    if not (tri(tjA, tj2, tj1) and tri(tJA, tJ2, tJ1)):
        return Fraction(0)
    w = wigner_9j(h(tjA), h(tj2), h(tj1), 1, h(1), h(1), h(tJA), h(tJ2), h(tJ1), prec=None)
    w2 = Rational(nsimplify(w ** 2))
    if not R1.left_open(tj1, tjA, tj2):
        return Fraction(0)
    val = (tjA + 1) * (tj2 + 1) * 6 * (tj1 + 1) * (tJA + 1) * (tJ2 + 1) * (tJ1 + 1) * w2
    val = Rational(val)
    return Fraction(int(val.p), int(val.q))


log("9j closed form vs brute force")
nine_mism = []
for k_, v in R1TAB.items():
    sp, comp, nn = k_.split("|")
    s1 = 1 if sp[0] == "+" else -1
    s2 = 1 if sp[1] == "+" else -1
    n1, n2, q = map(int, nn.split(","))
    if v["S"] != S_9j(s1, n1, comp, q, s2, n2):
        nine_mism.append(k_)
OUT["closed_form_9j_vs_bruteforce"] = {"mismatches": nine_mism, "n_checked": len(R1TAB)}
log("9j mismatches", len(nine_mism))

# ---------------------------------------------------------------------------
# C3: model restriction vs geovac get_photon_channels
# ---------------------------------------------------------------------------
log("C3")
from sympy import Rational as SR  # noqa: E402
from geovac.qed_anomalous_moment import get_photon_channels, vertex_allowed  # noqa: E402

c3_rows = []
c3_ok = True
for s in SIGNS:
    for n1 in range(N1MAX + 1):
        for n2 in range(N1MAX + 1):
            for q in range(1, Q1MAX + 1):
                if not model_rule(n1, n2, q):
                    if vertex_allowed(n1, n2, q):
                        c3_ok = False
                    continue
                if s > 0:
                    cont = (SR(n1 + 1, 2), SR(n1, 2), SR(n2 + 1, 2), SR(n2, 2))
                else:
                    cont = (SR(n1, 2), SR(n1 + 1, 2), SR(n2, 2), SR(n2 + 1, 2))
                code = len(get_photon_channels(n1, n2, q, *cont))
                r1 = sum(R1TAB[key(s, s, c, n1, n2, q)]["open"] for c in COMPS)
                c3_rows.append({"sign": sgn(s) + sgn(s), "n1": n1, "n2": n2, "q": q, "code_W": code, "route1_W": r1})
                if code != r1:
                    c3_ok = False
OUT["C3"] = {"pass": c3_ok, "rows": c3_rows,
             "note": "same-sign content passed to get_photon_channels exactly as qed_anomalous_moment does "
                     "(jL=(n+1)/2, jR=n/2 for +); -- uses the mirrored content."}
log("C3", c3_ok)

# ---------------------------------------------------------------------------
# route 2
# ---------------------------------------------------------------------------
import vector_vertex_s3_route2 as R2  # noqa: E402

log("route 2: spinor eigenspaces")
c1_r2 = {}
SP2: Dict[Tuple[int, int], List] = {}   # (Dsign, n) -> fields
CHOICE = {}
for c0 in (Fraction(3, 2), Fraction(-3, 2)):
    good = True
    det = {}
    spaces = {}
    for k in range(0, N2MAX + 2):
        es, dim = R2.spinor_eigenspaces(k, c0)
        tot = sum(len(v) for v in es.values())
        det[f"k{k}"] = {"dim_S+xH_k": dim, "eigs": {str(l): len(v) for l, v in es.items() if v}}
        if tot != dim:
            good = False
        for lam, fields in es.items():
            if not fields:
                continue
            s = 1 if lam > 0 else -1
            n = abs(lam) - Fraction(3, 2)
            if n.denominator != 1 or n < 0 or len(fields) != (int(n) + 1) * (int(n) + 2):
                good = False
                continue
            spaces[(s, int(n))] = fields
    # symmetric spectrum: both signs present at every level n <= N2MAX
    for n in range(N2MAX + 1):
        for s in SIGNS:
            if (s, n) not in spaces:
                good = False
    c1_r2[str(c0)] = {"symmetric_spectrum_pm(n+3/2)": good, "detail": det}
    if good:
        CHOICE["c0"] = c0
        SP2 = spaces
assert "c0" in CHOICE, "C1 route 2 failed for both constants"
c0 = CHOICE["c0"]
# explicit D-check on every basis field
dcheck = all(R2.check_spinor_eigen(f, c0, Fraction(s * (2 * n + 3), 2)) for (s, n), fs in SP2.items() for f in fs)
c1_r2_pass = dcheck and sum(1 for v in c1_r2.values() if v["symmetric_spectrum_pm(n+3/2)"]) == 1
OUT["C1_route2"] = {"pass": c1_r2_pass, "constant_c0": str(c0), "explicit_D_check": dcheck, "detail": c1_r2}
log("C1 r2", c1_r2_pass, "c0 =", c0)

log("route 2: vector eigenspaces")
VEC2: Dict[Tuple[int, int], List] = {}  # (stard sign, q)
c2_r2 = {}
for q in range(1, Q2MAX + 1):
    for sd in (1, -1):
        fs = R2.vector_eigenspace(q, sd * (q + 1))
        tang = all(not R2.x_dot(A) for A in fs)
        cocl = all(not R2.div(A) for A in fs)
        hodge = all(R2.eigen_ratio(A, R2.star_d(R2.star_d(A))) is not None
                    and R2.eigen_ratio(A, R2.star_d(R2.star_d(A))).re == (q + 1) ** 2 for A in fs)
        VEC2[(sd, q)] = fs
        c2_r2[f"q{q}{'+' if sd > 0 else '-'}"] = {"dim": len(fs), "dim_ok": len(fs) == q * (q + 2),
                                                 "tangent": tang, "coclosed": cocl, "hodge_(q+1)^2": hodge}
# exact forms are closed: *d(grad f) = 0 (orientation-free check of the *d formula)
exact_ok = True
for p in (1, 2, 3):
    for h in R2.harmonic_basis(p)[:4]:
        grad = [R2.p_d(h, a) for a in range(4)]
        # tangential part, homogenized: r^2 grad h - p x h
        r2h = [R2.p_add(R2.p_add(R2.p_add(R2.p_x(R2.p_x(g, 0), 0), R2.p_x(R2.p_x(g, 1), 1)),
                                 R2.p_add(R2.p_x(R2.p_x(g, 2), 2), R2.p_x(R2.p_x(g, 3), 3))),
                        R2.p_x(h, a), R2.G(-p)) for a, g in enumerate(grad)]
        if any(R2.star_d(r2h)):
            exact_ok = False
c2_r2_pass = exact_ok and all(all(x for x in v.values() if isinstance(x, bool)) for v in c2_r2.values())
OUT["C2_route2"] = {"pass": c2_r2_pass, "detail": c2_r2, "stard_kills_gradients": exact_ok,
                    "note": "Delta_Hodge = delta d on coclosed forms = (*d)^2 in 3D"}
log("C2 r2", c2_r2_pass)

# su(2) x su(2) content by Casimirs, calibrated on the scalar x1 (content (1/2,1/2))
log("route 2: content")
cal = {}
for sg in (1, -1):
    r = R2.eigen_ratio({(1, 0, 0, 0): R2.ONE}, R2.casimir({(1, 0, 0, 0): R2.ONE}, R2.G_scalar, sg))
    cal[sg] = r.re


def content_r2(fields, Gop) -> Tuple[int, int]:
    """(2 j_+, 2 j_-) measured; asserts the Casimir is scalar on every basis field."""
    res = []
    for sg in (1, -1):
        vals = set()
        for f in fields:
            r = R2.eigen_ratio(f, R2.casimir(f, Gop, sg))
            assert r is not None and r.im == 0, "Casimir not scalar"
            vals.add(r.re)
        assert len(vals) == 1
        jj = vals.pop() / cal[sg] * Fraction(3, 4)    # j(j+1)
        # solve 2j: (2j)(2j+2) = 4 j(j+1)
        t = 0
        while Fraction(t * (t + 2), 4) < jj:
            t += 1
        assert Fraction(t * (t + 2), 4) == jj
        res.append(t)
    return tuple(res)


CONT2 = {}
for (s, n), fs in SP2.items():
    CONT2[("spin", s, n)] = content_r2(fs, R2.G_spinor)
for (sd, q), fs in VEC2.items():
    CONT2[("vec", sd, q)] = content_r2(fs, R2.G_vector)
# first factor := the factor in which D=+ spinors carry the larger spin (record convention)
c_plus0 = CONT2[("spin", 1, 0)]
FIRST = 0 if c_plus0[0] > c_plus0[1] else 1


def first_second(c):
    return (c[FIRST], c[1 - FIRST])


content_ok = True
for (s, n) in SP2:
    want = (n + 1, n) if s > 0 else (n, n + 1)
    if first_second(CONT2[("spin", s, n)]) != want:
        content_ok = False
COMP2 = {}
for (sd, q) in VEC2:
    fsd = first_second(CONT2[("vec", sd, q)])
    if fsd == (q + 1, q - 1):
        COMP2[q, "A"] = sd
    elif fsd == (q - 1, q + 1):
        COMP2[q, "B"] = sd
    else:
        content_ok = False
# equivariance: D and *d commute with the so(4) generators (spot check on every space)
equiv_ok = True
for (s, n), fs in SP2.items():
    lam = Fraction(s * (2 * n + 3), 2)
    for a in range(4):
        for b in range(a + 1, 4):
            if not R2.check_spinor_eigen(R2.G_spinor(fs[0], a, b), c0, lam):
                equiv_ok = False
for (sd, q), fs in VEC2.items():
    for a in range(4):
        for b in range(a + 1, 4):
            g = R2.G_vector(fs[0], a, b)
            sdg = R2.star_d(g)
            if any(R2.p_add(x, y, R2.G(-sd * (q + 1))) for x, y in zip(sdg, g)):
                equiv_ok = False
OUT["route2_content"] = {"first_factor_index(0=J+,1=J-)": FIRST, "content_ok": content_ok,
                         "so4_equivariance_ok": equiv_ok,
                         "component_stard_sign": {f"{c}{q}": v for (q, c), v in COMP2.items()},
                         "spinor_content": {f"{sgn(s)}{n}": first_second(CONT2[('spin', s, n)]) for (s, n) in SP2}}
log("content", content_ok, "equiv", equiv_ok, "FIRST", FIRST, "comp *d signs", COMP2)

# orthogonal integer-scaled bases
log("route 2: Gram-Schmidt")
SPB, VCB = {}, {}
for kk, fs in SP2.items():
    ob, nrm = R2.gram_schmidt(fs)
    lst = []
    for f, n2 in zip(ob, nrm):
        fi, L = R2.to_int_field(f)
        lst.append((fi, n2 * L * L))
    SPB[kk] = lst
for (q, comp), sd in COMP2.items():
    ob, nrm = R2.gram_schmidt(VEC2[(sd, q)])
    lst = []
    for f, n2 in zip(ob, nrm):
        fi, L = R2.to_int_field(f)
        lst.append((fi, n2 * L * L))
    VCB[(q, comp)] = lst


def spin_deg(s, n):
    # degree of the polynomial spinors at (sign, n): D-sign with large multiplicity at k=n is
    # the one with (n+1)(n+2) states at degree k = n; the other sign sits at k = n+1
    f = SP2[(s, n)][0]
    for p in f:
        for k_ in p:
            return sum(k_)


log("route 2: vertex tensors")
R2TAB: Dict[str, Dict] = {}
WM_CACHE = {}
for s1 in SIGNS:
    for n1 in range(N2MAX + 1):
        U = SPB[(s1, n1)]
        d1 = spin_deg(s1, n1)
        for s2 in SIGNS:
            for n2 in range(N2MAX + 1):
                Vv = SPB[(s2, n2)]
                d2 = spin_deg(s2, n2)
                dK = d1 + d2 + 1
                curr = [[R2.current(u, v, dK) for (v, _) in Vv] for (u, _) in U]
                for q in range(1, Q2MAX + 1):
                    for comp in COMPS:
                        W = VCB[(q, comp)]
                        ck = (q, comp, dK)
                        if ck not in WM_CACHE:
                            WM_CACHE[ck] = R2.wm_matrix([w for (w, _) in W], q, dK)
                        WM = WM_CACHE[ck]
                        C = R2.moment_den(dK + q)
                        S = Fraction(0)
                        nz = False
                        for i, (u, nu) in enumerate(U):
                            for kv, (v, nv) in enumerate(Vv):
                                re, im = curr[i][kv]
                                vr = re.dot(WM)
                                vi = im.dot(WM)
                                part = Fraction(0)
                                for jw, (w, nw) in enumerate(W):
                                    a2 = vr[jw] * vr[jw] + vi[jw] * vi[jw]
                                    if a2:
                                        nz = True
                                        part += Fraction(a2, nw)
                                if part:
                                    S += part / (C * C * nu * nv)
                        R2TAB[key(s1, s2, comp, n1, n2, q)] = {"S": S, "open": nz}
        log("  done s1,n1 =", s1, n1)

# SA4 diagonal in route 2: <a| c(A) |a> for every basis state a at n = 0, q = 1
sa4_r2 = {}
for s in SIGNS:
    U = SPB[(s, 0)]
    d = spin_deg(s, 0)
    per = []
    for (u, nu) in U:
        re, im = R2.current(u, u, 2 * d + 1)
        tot = Fraction(0)
        for comp in COMPS:
            W = VCB[(1, comp)]
            ck = (1, comp, 2 * d + 1)
            WM = WM_CACHE[ck] if ck in WM_CACHE else R2.wm_matrix([w for (w, _) in W], 1, 2 * d + 1)
            C = R2.moment_den(2 * d + 2)
            vr, vi = re.dot(WM), im.dot(WM)
            for jw, (w, nw) in enumerate(W):
                tot += Fraction(vr[jw] ** 2 + vi[jw] ** 2, C * C * nu * nu * nw)
        per.append(str(tot))
    sa4_r2[sgn(s)] = per
OUT["SA4_route2_sum_photon_|<a|c(A)|a>|^2_per_basis_state"] = sa4_r2
log("SA4 r2", sa4_r2)

# Ward identity (Paper 33 rule 6, reported): pure-gauge photon A = d f couples as
# <u, c(grad_S f) v> = (lam_u - lam_v) <u, f v>   (operator identity [D, f] = c(df)).
log("route 2: Ward check")
ward_ok, ward_n, ward_nonzero = True, 0, 0
for p in (1, 2):
    for f in R2.harmonic_basis(p)[:3]:
        grad = [R2.p_d(f, a) for a in range(4)]
        r2g = []
        for a, g in enumerate(grad):
            t = {}
            for b in range(4):
                t = R2.p_add(t, R2.p_x(R2.p_x(g, b), b))
            r2g.append(R2.p_add(t, R2.p_x(f, a), R2.G(-p)))     # tangential gradient, homogenized

        def cliff(Avec, psi):
            out = [{}, {}]
            for a in range(4):
                for b in range(4):
                    xa = [R2.p_x(R2.p_mul(Avec[b], psi[0]), a), R2.p_x(R2.p_mul(Avec[b], psi[1]), a)]
                    t = R2.mat_apply(R2.TT[a][b], xa)
                    out = [R2.p_add(out[0], t[0]), R2.p_add(out[1], t[1])]
            return out
        for (su, nu_) in [(1, 0), (-1, 0), (1, 1), (-1, 1)]:
            for (sv, nv_) in [(1, 0), (-1, 0), (1, 1), (-1, 1)]:
                lu = Fraction(su * (2 * nu_ + 3), 2)
                lv = Fraction(sv * (2 * nv_ + 3), 2)
                for u in SP2[(su, nu_)][:2]:
                    for v in SP2[(sv, nv_)][:2]:
                        lhs = R2.inner(u, cliff(r2g, v))
                        rhs = R2.inner(u, [R2.p_mul(f, v[0]), R2.p_mul(f, v[1])]) * (lu - lv)
                        ward_n += 1
                        if not (lhs - rhs).iszero():
                            ward_ok = False
                        if not lhs.iszero():
                            ward_nonzero += 1
OUT["Ward_route2"] = {"pass": ward_ok, "n_checked": ward_n, "n_nonzero_lhs": ward_nonzero,
                      "identity": "<u, c(grad f) v> = (lam_u - lam_v) <u, f v>, f harmonic deg 1,2; n <= 1"}
log("Ward", ward_ok, ward_n, ward_nonzero)

# ---------------------------------------------------------------------------
# C4
# ---------------------------------------------------------------------------
c4_mism = []
n_c4 = 0
for k_, v in R2TAB.items():
    n_c4 += 1
    w = R1TAB[k_]
    if v["open"] != w["open"] or v["S"] != w["S"]:
        c4_mism.append({"key": k_, "r1": [w["open"], str(w["S"])], "r2": [v["open"], str(v["S"])]})
c4_pass = not c4_mism
OUT["C4"] = {"pass": c4_pass, "n_compared": n_c4, "mismatches": c4_mism[:50]}
log("C4", c4_pass, n_c4, len(c4_mism))

controls_pass = c1_r1_pass and c1_r2_pass and c2_r1_pass and c2_r2_pass and c3_ok and c4_pass
OUT["controls_pass"] = controls_pass

# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------
def score_set(pred_fn, keys):
    extras, acc0 = [], []
    for k_ in keys:
        sp, comp, nn = k_.split("|")
        s1 = 1 if sp[0] == "+" else -1
        s2 = 1 if sp[1] == "+" else -1
        n1, n2, q = map(int, nn.split(","))
        pred = pred_fn(s1, s2, comp, n1, n2, q)
        got = R1TAB[k_]["open"]
        if got and not pred:
            extras.append(k_)
        if pred and not got:
            acc0.append(k_)
    return extras, acc0


ALL = list(R1TAB)
same = [k_ for k_ in ALL if k_[0] == k_[1]]
opp = [k_ for k_ in ALL if k_[0] != k_[1]]
SC = {}
if controls_pass:
    ex, az = score_set(rule_record, same)
    top = []
    for s in SIGNS:
        for n1 in range(N1MAX + 1):
            for n2 in range(N1MAX + 1):
                q = n1 + n2 + 1
                nopen = sum(R1TAB[key(s, s, c, n1, n2, q)]["open"] for c in COMPS)
                if nopen != 1:
                    top.append(key(s, s, "*", n1, n2, q))
    SC["SA1"] = {"verdict": "HELD" if not ex and not az and not top else "MISSED",
                 "extras": ex, "accidental_zeros": az, "q=n1+n2+1_not_exactly_one_component": top}
    ex, az = score_set(rule_record, [k_ for k_ in opp if k_.startswith("+-")])
    ex2, az2 = score_set(rule_record, [k_ for k_ in opp if k_.startswith("-+")])
    parity_bad = [k_ for k_ in opp if R1TAB[k_]["open"] and sum(map(int, k_.split("|")[2].split(","))) % 2 != 0]
    nonempty = any(R1TAB[k_]["open"] for k_ in opp)
    witness = R1TAB[key(1, -1, "B", 0, 1, 1)]["open"]
    SC["SA2"] = {"verdict": "HELD" if not ex and not az and not parity_bad and nonempty and witness else "MISSED",
                 "extras(+-)": ex, "accidental_zeros(+-)": az, "odd_parity_open": parity_bad,
                 "nonempty": nonempty, "witness_(+-,B,0,1,1)_open": witness,
                 "reported_-+_(transposed rule, not stated in the record)": {"extras": ex2, "accidental_zeros": az2}}
    # SA3
    sa3_bad = []
    for ext_pos in (1, 2):
        for se in SIGNS:
            for si in SIGNS:
                for comp in COMPS:
                    for ni in range(N1MAX + 1):
                        for q in range(1, Q1MAX + 1):
                            if ext_pos == 1:
                                k_ = key(se, si, comp, 0, ni, q)
                            else:
                                k_ = key(si, se, comp, ni, 0, q)
                            got = R1TAB[k_]["open"]
                            if se == si:
                                want_comp = "A" if se > 0 else "B"
                                pred = (q == ni + 1 and comp == want_comp)
                            else:
                                want_comp = "B" if se > 0 else "A"
                                pred = (q == ni and ni >= 1 and comp == want_comp)
                            if got != pred:
                                sa3_bad.append(k_)
    SC["SA3"] = {"verdict": "HELD" if not sa3_bad else "MISSED", "mismatches": sa3_bad}
    # SA4
    sa4 = {s_: R1.diagonal_ground(s_) for s_ in SIGNS}
    sa4_ok = True
    for s_, d in sa4.items():
        vals = []
        for comp, info in d.items():
            if info.get("left_open"):
                vals.append(info["sum_mA_|V(a,A,a)|^2"])
        per_state_any = {}
        for dct in vals:
            for st, val in dct.items():
                per_state_any[st] = per_state_any.get(st, False) or (val != "0")
        if not per_state_any or not all(per_state_any.values()):
            sa4_ok = False
    sa4_ok = sa4_ok and all(Fraction(x) != 0 for v in sa4_r2.values() for x in v)
    SC["SA4"] = {"verdict": "HELD" if sa4_ok else "MISSED", "route1": sa4, "route2": sa4_r2}
    # SA5
    az5 = [k_ for k_ in ALL if rule_record(*(lambda sp, c, nn: ((1 if sp[0] == '+' else -1), (1 if sp[1] == '+' else -1), c, *map(int, nn.split(","))))(*k_.split("|"))) and not R1TAB[k_]["open"]]
    SC["SA5"] = {"verdict": "HELD" if not az5 else "MISSED", "accidental_zeros": az5}
else:
    SC = {"VOID": True}
OUT["scores"] = SC
log("scores", {k: v.get("verdict") for k, v in SC.items() if isinstance(v, dict)})

# ---------------------------------------------------------------------------
# reported: Sigma_phys(n_ext)
# ---------------------------------------------------------------------------
log("Sigma tables")
from sympy.physics.wigner import wigner_9j  # noqa: E402

S_CACHE: Dict = {}


def S_any(s1, n1, comp, q, s2, n2) -> Fraction:
    k_ = key(s1, s2, comp, n1, n2, q)
    if k_ in R1TAB:
        return R1TAB[k_]["S"]
    if k_ not in S_CACHE:
        S_CACHE[k_] = S_9j(s1, n1, comp, q, s2, n2)
    return S_CACHE[k_]


sig = {}
for n_ext in range(4):
    for nmax in (3, 5, 8, 12):
        cnt = Fraction(0)
        cnt_hodge = Fraction(0)
        amp = Fraction(0)
        amp_hodge = Fraction(0)
        for s_int in SIGNS:
            for n_int in range(nmax + 1):
                lam4 = Fraction(2 * n_int + 3, 2) ** 4
                g = (n_int + 1) * (n_int + 2)
                for q in range(1, n_ext + n_int + 2):
                    for comp in COMPS:
                        S = S_any(1, n_ext, comp, q, s_int, n_int)
                        if S == 0:
                            continue
                        dT = q * (q + 2)
                        cnt += Fraction(g * dT) / (lam4 * dT)        # mu_q = q(q+2) (eq:hodge1_spectrum)
                        cnt_hodge += Fraction(g * dT) / (lam4 * (q + 1) ** 2)
                        d_ext = (n_ext + 1) * (n_ext + 2)
                        amp += S / (d_ext * lam4 * dT)
                        amp_hodge += S / (d_ext * lam4 * (q + 1) ** 2)
        sig[f"n_ext={n_ext},n_max={nmax}"] = {
            "count_muq(q+2)": [str(cnt), float(cnt)], "count_mu(q+1)^2": [str(cnt_hodge), float(cnt_hodge)],
            "amp_muq(q+2)": [str(amp), float(amp)], "amp_mu(q+1)^2": [str(amp_hodge), float(amp_hodge)]}
OUT["Sigma_phys"] = {
    "definition": "ext sign +, psi1 = ext, both internal signs; count: sum_{s_int,n_int,q,comp open} "
                  "g_s(n_int) d_q^T/(|lambda|^4 mu_q), g_s = (n+1)(n+2) per sign; amp: W*g*d replaced by "
                  "Sum_{m_int,m_gamma}|V|^2 at fixed external state = Sum_all|V|^2/d_ext (Haar-normalized "
                  "volume, unnormalized coupling); mu_q = q(q+2) as in eq:hodge1_spectrum, and (q+1)^2 variant.",
    "table": sig}
# mirror: ext sign - gives the same (checked at n_max = 5)
mir_ok = True
for n_ext in range(4):
    a = b = Fraction(0)
    for s_int in SIGNS:
        for n_int in range(6):
            for q in range(1, n_ext + n_int + 2):
                for comp in COMPS:
                    a += S_any(1, n_ext, comp, q, s_int, n_int)
                    b += S_any(-1, n_ext, comp, q, -s_int, n_int)
    mir_ok = mir_ok and a == b
OUT["Sigma_phys"]["ext_sign_mirror_ok"] = mir_ok
log("Sigma done, mirror", mir_ok)

# extra observation: geovac.qed_self_energy's W uses mixed-sign content without integrality
from geovac.qed_self_energy import _so4_channel_count  # noqa: E402
obs = []
for n1 in range(4):
    for n2 in range(4):
        for q in range(1, 8):
            w_code = _so4_channel_count(n1, n2, q)
            w_pp = sum(R1TAB[key(1, 1, c, n1, n2, q)]["open"] for c in COMPS)
            w_pm = sum(R1TAB[key(1, -1, c, n1, n2, q)]["open"] for c in COMPS)
            if w_code != w_pp:
                obs.append({"n1": n1, "n2": n2, "q": q, "in_model_window": model_rule(n1, n2, q), "self_energy_code_W": w_code, "phys_W++": w_pp, "phys_W+-": w_pm})
OUT["observation_qed_self_energy_W_vs_physical"] = obs

OUT["m_conservation_both_factors_route1"] = mcons_all
OUT["route1_open_set"] = {k_: [v["open"], str(v["S"])] for k_, v in R1TAB.items()}
OUT["route2_open_set"] = {k_: [v["open"], str(v["S"])] for k_, v in R2TAB.items()}
with open(__file__, "rb") as fh:
    OUT["driver_sha256"] = hashlib.sha256(fh.read()).hexdigest()
os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
with open(os.path.join(HERE, "data", "vector_vertex_s3.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, default=str)
log("wrote json; sha256", OUT["driver_sha256"])
