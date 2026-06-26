# -*- coding: utf-8 -*-
"""B3 Phase 3, Sprint 3b: the fold-transfer rule in closed form (2026-06-10).

Exact-arithmetic (sympy Clebsch-Gordan) proof-by-computation of the wedge
fold-transfer rule observed numerically in Sprint 3, per the
discrete-for-skeleton rule (Layer-1 claim => exact arithmetic, no PSLQ).

The objects: window compressions of band-b multipliers on the Peter-Weyl
window j <= j_max,

  C^b_{mu',mu}[(j1,m1',m1),(j2,m2',m2)]
      = sqrt((2j2+1)(2b+1)/(2j1+1)) <b mu'; j2 m2' | j1 m1'> <b mu; j2 m2 | j1 m1>,

Hermitized G = C + C^T (real entries), the plain-swap wedge reflection
R: m' -> -m', and the folding G_W = V^dag G V onto the +1 eigenspace.

Checks:
  T0  convention validation: exact entries reproduce the numeric (quadrature)
      compressions of the b1 machinery for all seven Sprint-1 classes (j<=1).
  T1  BLOCK REFLECTION IDENTITY (exact, the lemma):
        [R C^b_{mu',mu} R]_(j1,j2) = (-1)^{b+j2-j1} [C^b_{-mu',mu}]_(j1,j2).
      Pure CG symmetry <j1,-m1'|b,mu',j2,-m2'> = (-1)^{b+j2-j1}<j1 m1'|b,-mu',j2,m2'>.
  T2  PARITY RULE (mu = 0 column, exact): each (j1,j2) block of G is an
      R-conjugation eigenblock with parity eps = (-1)^{b+j2-j1+mu'};
      blocks with eps = -1 are annihilated by the folding, eps = +1 survive.
      Exact Frobenius^2 fold ratios; (2,1) annihilation at j_max = 1 is the
      special case "all blocks odd" (only block (1,1), eps = -1).
  T3  mu', mu both nonzero: R maps the (mu',mu) component to (-mu',mu),
      disjoint from G's span -- no parity eigenstructure, partial fold
      ratios (exact values reported).
  T4  WINDOW-EDGE THEOREM (j_max = 3/2): both Sprint-3 headline facts are
      j_max = 1 edge effects --
        (a) (2,1) is NO LONGER annihilated (half-integer blocks have
            eps = +1 and survive);
        (b) folded (2,2) NO LONGER commutes with K_W (non-mirror
            transitions m' -> m'+2 with m' != -1 enter; at j_max = 1 the
            range forces the single mirror transition -1 -> +1).
      The invariant content is the parity rule + the mirror decomposition
      (mirror part = transitions m1' = -m2', always flow-commuting), not
      the annihilation/conversion themselves. mu' = 0 classes fold to
      weight-diagonal (commuting) generators at EVERY window (invariant).

Frozen falsifier: tests/test_wh7_b3_fold_rule.py.
"""
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

HALF = sp.Rational(1, 2)

# the seven Sprint-1 classes as (b, mu', mu) in exact rationals
CLASSES = {
    "(0.5,0.5) spacelike": (HALF, HALF, HALF),
    "(1.0,0.0) spacelike": (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
    "(2.0,0.0) spacelike": (sp.Integer(2), sp.Integer(0), sp.Integer(0)),
    "(2.0,1.0) spacelike": (sp.Integer(2), sp.Integer(1), sp.Integer(0)),
    "(1.0,1.0) null":      (sp.Integer(1), sp.Integer(1), sp.Integer(0)),
    "(1.5,1.5) timelike":  (sp.Rational(3, 2), sp.Rational(3, 2), HALF),
    "(2.0,2.0) timelike":  (sp.Integer(2), sp.Integer(2), sp.Integer(0)),
}
MU0_CLASSES = [n for n, (b, mp, mu) in CLASSES.items() if mu == 0]


def jlist(jmax):
    out, j = [], sp.Integer(0)
    while j <= jmax:
        out.append(j)
        j += HALF
    return out


def labels(jmax):
    labs = []
    for j in jlist(jmax):
        ms = [j - k for k in range(int(2 * j) + 1)]
        for mp in ms:
            for m in ms:
                labs.append((j, mp, m))
    return labs


_CG_CACHE = {}


def cg(a, al, b, be, c, ga):
    key = (a, al, b, be, c, ga)
    if key not in _CG_CACHE:
        _CG_CACHE[key] = CG(a, al, b, be, c, ga).doit()
    return _CG_CACHE[key]


def build_C(b, mu_p, mu, labs):
    """Sparse dict {(row_label, col_label): exact entry} of C^b_{mu',mu}."""
    out = {}
    for (j2, m2p, m2) in labs:
        m1p, m1 = m2p + mu_p, m2 + mu
        for j1 in jlist(max(j for j, _, _ in labs)):
            if abs(m1p) > j1 or abs(m1) > j1:
                continue
            if not (abs(b - j2) <= j1 <= b + j2):
                continue
            c1 = cg(b, mu_p, j2, m2p, j1, m1p)
            if c1 == 0:
                continue
            c2 = cg(b, mu, j2, m2, j1, m1)
            if c2 == 0:
                continue
            N = sp.sqrt((2 * j2 + 1) * (2 * b + 1) / (2 * j1 + 1))
            out[((j1, m1p, m1), (j2, m2p, m2))] = N * c1 * c2
    return out


def herm(C):
    """G = C + C^T (entries are real)."""
    G = {}
    for (i, k), v in C.items():
        G[(i, k)] = G.get((i, k), 0) + v
        G[(k, i)] = G.get((k, i), 0) + v
    return {ik: sp.simplify(v) for ik, v in G.items() if sp.simplify(v) != 0}


def refl(lab):
    j, mp, m = lab
    return (j, -mp, m)


def conj_R(M):
    """R M R for the plain-swap reflection (permutation conjugation)."""
    return {(refl(i), refl(k)): v for (i, k), v in M.items()}


def blocks_of(M):
    bl = {}
    for (i, k), v in M.items():
        bl.setdefault((i[0], k[0]), {})[(i, k)] = v
    return bl


def is_zero(e):
    s = sp.simplify(sp.expand(e))
    if s == 0:
        return True
    r = s.equals(0)
    return bool(r) if r is not None else False


def fold(G, labs):
    """G_W = V^dag G V onto the R +1 eigenspace. Wedge labels: (j, mp, m)
    with mp > 0 (sym pair, 1/sqrt2 each) or mp = 0 (R-fixed)."""
    wlabs = [l for l in labs if l[1] >= 0]

    def comps(l):
        j, mp, m = l
        if mp == 0:
            return [(l, sp.Integer(1))]
        s = 1 / sp.sqrt(2)
        return [(l, s), ((j, -mp, m), s)]

    GW = {}
    for wi in wlabs:
        for wk in wlabs:
            tot = sp.Integer(0)
            for (li, ci) in comps(wi):
                for (lk, ck) in comps(wk):
                    if (li, lk) in G:
                        tot += ci * ck * G[(li, lk)]
            tot = sp.simplify(tot)
            if tot != 0:
                GW[(wi, wk)] = tot
    return GW, wlabs


def frob2(M):
    return sp.simplify(sum(v ** 2 for v in M.values()))


def kw_commutes(GW):
    """[K_W, G_W] = 0 iff every entry connects equal |2 m'| weights.
    Also return the mirror/non-mirror split (mirror: |mp_i| == |mp_k|)."""
    nonmirror = {ik: v for ik, v in GW.items()
                 if abs(ik[0][1]) != abs(ik[1][1])}
    return len(nonmirror) == 0, nonmirror


def to_float_matrix(M, labs):
    idx = {l: i for i, l in enumerate(labs)}
    A = np.zeros((len(labs), len(labs)))
    for (i, k), v in M.items():
        A[idx[i], idx[k]] = float(sp.N(v))
    return A


def run():
    out = {}
    jmax1 = sp.Integer(1)
    labs1 = labels(jmax1)

    # ---------- T0: convention validation against the b1 numerics ----------
    import wh7_b3_phase3_state_intervals as s1  # noqa: E402 (numeric side)
    import wh7_b1_joint_product_gh as b1  # noqa: E402
    order = [(sp.nsimplify(j), sp.nsimplify(j) - a, sp.nsimplify(j) - c)
             for (j, a, c) in b1.LABELS]
    t0 = 0.0
    SROWCOL = {"(0.5,0.5) spacelike": (0.5, 0, 0),
               "(1.0,0.0) spacelike": (1.0, 1, 1),
               "(2.0,0.0) spacelike": (2.0, 2, 2),
               "(2.0,1.0) spacelike": (2.0, 1, 2),
               "(1.0,1.0) null": (1.0, 0, 1),
               "(1.5,1.5) timelike": (1.5, 0, 1),
               "(2.0,2.0) timelike": (2.0, 0, 2)}
    for name, (b, mp, mu) in CLASSES.items():
        G = herm(build_C(b, mp, mu, labs1))
        A = to_float_matrix(G, order)
        nrm = np.linalg.norm(A, 2)
        bf, rf, cf = SROWCOL[name]
        Gnum = s1.class_generator(bf, rf, cf)
        t0 = max(t0, float(np.max(np.abs(A / nrm - np.real(Gnum)))))
    out["T0_numeric_crosscheck_max_dev"] = t0

    # ---------- T1: block reflection identity (exact) ----------
    t1_ok, t1_checked = True, 0
    for name, (b, mp, mu) in CLASSES.items():
        C = build_C(b, mp, mu, labs1)
        Cr = build_C(b, -mp, mu, labs1)
        RCR = conj_R(C)
        keys = set(RCR) | set(Cr)
        for (i, k) in keys:
            j1, j2 = i[0], k[0]
            tau = sp.Integer(-1) ** sp.nsimplify(b + j2 - j1)
            d = RCR.get((i, k), 0) - tau * Cr.get((i, k), 0)
            t1_checked += 1
            if not is_zero(d):
                t1_ok = False
    out["T1_block_reflection_identity"] = {"holds": t1_ok,
                                           "entries_checked": t1_checked}

    # ---------- T2: parity rule on the mu = 0 column (exact) ----------
    t2 = {}
    for name in MU0_CLASSES:
        b, mp, mu = CLASSES[name]
        G = herm(build_C(b, mp, mu, labs1))
        RGR = conj_R(G)
        bl_G, bl_R = blocks_of(G), blocks_of(RGR)
        entry = {"blocks": {}, "parity_rule_holds": True}
        for blk in sorted(set(bl_G) | set(bl_R), key=str):
            j1, j2 = blk
            eps = sp.Integer(-1) ** sp.nsimplify(b + j2 - j1 + mp)
            ok = all(is_zero(bl_R.get(blk, {}).get(ik, 0)
                             - eps * v)
                     for ik, v in bl_G.get(blk, {}).items())
            ok = ok and all(is_zero(bl_G.get(blk, {}).get(ik, 0)
                                    - v / eps)
                            for ik, v in bl_R.get(blk, {}).items())
            entry["blocks"][str(blk)] = {"eps": int(eps), "eigen_ok": bool(ok)}
            if not ok:
                entry["parity_rule_holds"] = False
        GW, _ = fold(G, labs1)
        f2G, f2W = frob2(G), frob2(GW)
        ratio = sp.simplify(f2W / f2G) if f2G != 0 else sp.Integer(0)
        entry["fold_frob2_ratio"] = str(ratio)
        entry["fold_frob2_ratio_float"] = float(sp.N(ratio))
        # the falsifiable form of the rule: PER BLOCK, the fold annihilates
        # the block iff eps = -1 (an eps = +1 block folds to its P-restriction,
        # nonzero but generally Frobenius-partial -- the complement lives in
        # the antisymmetric wedge complement)
        per_block_ok = True
        for blk, Mb in bl_G.items():
            eps = sp.Integer(-1) ** sp.nsimplify(b + blk[1] - blk[0] + mp)
            GWb, _ = fold(Mb, labs1)
            dead = (len(GWb) == 0)
            if dead != (eps == -1):
                per_block_ok = False
        entry["block_annihilation_iff_eps_odd"] = bool(per_block_ok)
        entry["annihilated"] = (len(GW) == 0)
        comm, nonm = kw_commutes(GW)
        entry["KW_commutes"] = bool(comm)
        t2[name] = entry
    out["T2_parity_rule_mu0"] = t2

    # ---------- T3: mu', mu both nonzero (partial folds, exact) ----------
    t3 = {}
    for name in ("(0.5,0.5) spacelike", "(1.5,1.5) timelike"):
        b, mp, mu = CLASSES[name]
        G = herm(build_C(b, mp, mu, labs1))
        GW, _ = fold(G, labs1)
        ratio = sp.simplify(frob2(GW) / frob2(G))
        t3[name] = {"fold_frob2_ratio": str(ratio),
                    "fold_frob2_ratio_float": float(sp.N(ratio)),
                    "strictly_partial": bool(0 < float(sp.N(ratio)) < 1)}
    out["T3_partial_folds"] = t3

    # ---------- T4: window-edge theorem at j_max = 3/2 ----------
    labs32 = labels(sp.Rational(3, 2))
    t4 = {}
    # (a) (2,1) revives: half-integer blocks have eps = +1
    b, mp, mu = CLASSES["(2.0,1.0) spacelike"]
    G21 = herm(build_C(b, mp, mu, labs32))
    GW21, _ = fold(G21, labs32)
    surv = sorted({(str(i[0]), str(k[0])) for (i, k) in GW21})
    t4["(2,1)_revived_at_3/2"] = {
        "GW_nonzero": len(GW21) > 0,
        "fold_frob2_ratio": str(sp.simplify(frob2(GW21) / frob2(G21))),
        "surviving_blocks": [f"({a},{c})" for a, c in surv]}
    # (b) folded (2,2) no longer commutes: non-mirror transitions enter
    b, mp, mu = CLASSES["(2.0,2.0) timelike"]
    G22 = herm(build_C(b, mp, mu, labs32))
    GW22, _ = fold(G22, labs32)
    comm32, nonm32 = kw_commutes(GW22)
    # mirror decomposition: the mirror part always commutes
    mirror = {ik: v for ik, v in GW22.items()
              if abs(ik[0][1]) == abs(ik[1][1])}
    mfrac = sp.simplify(frob2(mirror) / frob2(GW22)) if GW22 else sp.Integer(0)
    t4["(2,2)_at_3/2"] = {
        "KW_commutes": bool(comm32),
        "n_nonmirror_entries": len(nonm32),
        "mirror_part_nonzero": len(mirror) > 0,
        "mirror_frob2_fraction_exact": str(mfrac),
        "mirror_frob2_fraction": float(sp.N(mfrac))}
    # at j_max = 1 the folded (2,2) IS purely mirror (re-pin)
    G22_1 = herm(build_C(b, mp, mu, labs1))
    GW22_1, _ = fold(G22_1, labs1)
    comm1, _ = kw_commutes(GW22_1)
    t4["(2,2)_at_1_purely_mirror"] = bool(comm1) and all(
        abs(ik[0][1]) == abs(ik[1][1]) for ik in GW22_1)
    # invariant: mu' = 0 classes commute at BOTH windows
    inv = True
    for name in ("(1.0,0.0) spacelike", "(2.0,0.0) spacelike"):
        bb, mmp, mmu = CLASSES[name]
        for ll in (labs1, labs32):
            GWm, _ = fold(herm(build_C(bb, mmp, mmu, ll)), ll)
            inv = inv and kw_commutes(GWm)[0]
    t4["mu0_classes_commute_both_windows"] = bool(inv)
    out["T4_window_edge"] = t4
    return out


if __name__ == "__main__":
    res = run()
    (ROOT / "data" / "wh7_b3_phase3_sprint3b_fold_rule.json").write_text(
        json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(f"T0 numeric cross-check max dev: {res['T0_numeric_crosscheck_max_dev']:.2e}")
    t1 = res["T1_block_reflection_identity"]
    print(f"T1 block reflection identity (exact): holds={t1['holds']} "
          f"({t1['entries_checked']} entries)")
    print("T2 parity rule, mu = 0 column (exact):")
    for name, e in res["T2_parity_rule_mu0"].items():
        bl = " ".join(f"{k}:eps={v['eps']}{'ok' if v['eigen_ok'] else 'FAIL'}"
                      for k, v in e["blocks"].items())
        print(f"   {name:24s}: rule={e['parity_rule_holds']} "
              f"fold_F2={e['fold_frob2_ratio']} "
              f"(={e['fold_frob2_ratio_float']:.4f}, "
              f"blk_annih_iff_odd={e['block_annihilation_iff_eps_odd']}) "
              f"annihilated={e['annihilated']} KW_comm={e['KW_commutes']}")
        print(f"        {bl}")
    print("T3 partial folds (mu', mu != 0):")
    for name, e in res["T3_partial_folds"].items():
        print(f"   {name:24s}: fold_F2={e['fold_frob2_ratio']} "
              f"(={e['fold_frob2_ratio_float']:.4f}) "
              f"partial={e['strictly_partial']}")
    t4 = res["T4_window_edge"]
    print("T4 window-edge theorem at j_max=3/2:")
    a = t4["(2,1)_revived_at_3/2"]
    print(f"   (2,1) revived: GW != 0: {a['GW_nonzero']} "
          f"fold_F2={a['fold_frob2_ratio']} blocks={a['surviving_blocks']}")
    bb = t4["(2,2)_at_3/2"]
    print(f"   (2,2) at 3/2: KW_commutes={bb['KW_commutes']} "
          f"nonmirror={bb['n_nonmirror_entries']} "
          f"mirror_F2_frac={bb['mirror_frob2_fraction']:.4f}")
    print(f"   (2,2) at 1 purely mirror: {t4['(2,2)_at_1_purely_mirror']}")
    print(f"   mu'=0 classes commute at both windows: "
          f"{t4['mu0_classes_commute_both_windows']}")
