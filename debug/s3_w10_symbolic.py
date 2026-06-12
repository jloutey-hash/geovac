"""S^(3) w10 layer — symbolic-only reduction attempt (NO PSLQ, PI directive).

Target: W10 = -2048 t3(4,3,3) - 1024 t3(4,4,2) - 512 t2(8,2) - 1536 t2(4,6)
(coefficients read from debug/data/s3_decomp_table.json), the weight-10
multiple-t layer of the S^(3) decomposition.

Convention (matches s3_decomp_table.json / s3_decomp_numerics.py):
  t(s1,...,sk) = sum_{o1>o2>...>ok>=1 odd} o1^-s1 ... ok^-sk
  (first argument on the LARGEST variable); t(s) = lam(s) = (1-2^-s) zeta(s).
These are Hoffman multiple t-values (arXiv:1612.05232).

Only EXACT relations are used: the quasi-shuffle (stuffle) expansions of
products of convergent multiple t-values.  The stuffle product is an index-
combinatorial identity (interleave the summation chains; equalities merge
exponents), valid verbatim for the odd-restricted sums:
  t(u) * t(w) = sum of quasi-shuffle terms qsh(u,w).
No PSLQ, no integer-relation search, no numerically-tuned coefficient
appears anywhere in this driver.

Three nested exact linear systems over Q (unknowns = convergent weight-10
words of depth >= 2; the depth-1 word (10) = lam(10) is a known):

  Stage 1 (charter core): products lam*lam and lam*t2 -> relations supported
          on depth <= 3 words (8 double + 28 triple unknowns).
  Stage 2: + lam*t3 and t2*t2 -> depth <= 4 words (+56 quadruple unknowns,
          which must cancel exactly in any reduction of the depth-<=3
          target; no elimination shortcuts -- membership is tested in the
          full coordinate space).
  Stage 3 (diagnostic, symbolic only): the FULL stuffle product closure
          at weight 10 -- every pairwise product of convergent words of
          weights (2,8),(3,7),(4,6),(5,5).  Settles definitively whether
          W10 lies in the stuffle-decomposable subspace.

Gaussian elimination over Fraction with known-side tracking: reducing the
target vector v against the pivot rows yields
   W10 = sum_r c_r * (product known)_r + sum_j d_j * (free word)_j
exactly; the d_j-part is the stuffle-irreducible residual ("surviving
generators").

Numerics subcommand verifies every Stage-1 relation and every derived
identity against high-precision values (cached 220-dps values in
s3_pslq_cache.json where available; missing values computed by the
Levin-safe standard evaluator -- all summands here are log-free -- and by
the Abel-summation Trailing1 evaluator for t3(b1,b2,1), per the closure
sprint's corrected design; cf. the mpmath-levin-on-logs failure mode,
CLAUDE.md section 3).

Usage:
  python s3_w10_symbolic.py symbolic           # exact linear algebra
  python s3_w10_symbolic.py numerics [dps]     # verify relations+identities
Output: debug/data/s3_w10_symbolic.json (cumulative, checkpointed)
        debug/data/s3_w10_value_cache.json (new high-precision values)
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / "data"
OUT_PATH = DATA / "s3_w10_symbolic.json"
VCACHE_PATH = DATA / "s3_w10_value_cache.json"
SHARED_CACHE_PATH = DATA / "s3_pslq_cache.json"

# ===================================================================
# quasi-shuffle product on composition words (first index = largest var)
# ===================================================================
_QSH_MEMO: dict = {}


def qsh(u: tuple, w: tuple) -> dict:
    """Quasi-shuffle expansion of t(u)*t(w) as {word: int coefficient}."""
    key = (u, w)
    if key in _QSH_MEMO:
        return _QSH_MEMO[key]
    if not u:
        r = {w: 1}
    elif not w:
        r = {u: 1}
    else:
        r = {}
        for t, c in qsh(u[1:], w).items():
            k = (u[0],) + t
            r[k] = r.get(k, 0) + c
        for t, c in qsh(u, w[1:]).items():
            k = (w[0],) + t
            r[k] = r.get(k, 0) + c
        for t, c in qsh(u[1:], w[1:]).items():
            k = (u[0] + w[0],) + t
            r[k] = r.get(k, 0) + c
    _QSH_MEMO[key] = r
    return r


def comps(total: int, first_min: int = 2) -> list:
    """Compositions (s1,..,sk), s1 >= first_min, si >= 1, k >= 1."""
    out = []

    def rec(rem, acc):
        if rem == 0:
            if acc:
                out.append(tuple(acc))
            return
        lo = first_min if not acc else 1
        for p in range(lo, rem + 1):
            rec(rem - p, acc + [p])

    rec(total, [])
    return out


def wlabel(wd: tuple) -> str:
    if len(wd) == 1:
        return "lam(%d)" % wd[0]
    return "t%d(%s)" % (len(wd), ",".join(map(str, wd)))


def klabel(k) -> str:
    if k[0] == "prod":
        return "%s*%s" % (wlabel(k[1]), wlabel(k[2]))
    if k[0] == "word":
        return wlabel(k[1])
    raise ValueError(k)


# ===================================================================
# relation construction: t(u)*t(w) = sum(qsh) ->  vec . x = known
# ===================================================================
def build_relation(u: tuple, w: tuple):
    """Return (vec, known): vec.x = known, x = depth>=2 weight-10 words.

    Depth-1 expansion terms (the word (10,)) are moved to the known side.
    known is a dict over symbols ('prod',u,w) and ('word',(10,))."""
    exp = qsh(u, w)
    vec, known = {}, {("prod", u, w): Fraction(1)}
    for word, c in exp.items():
        if len(word) == 1:
            known[("word", word)] = known.get(("word", word),
                                              Fraction(0)) - Fraction(c)
        else:
            vec[word] = vec.get(word, Fraction(0)) + Fraction(c)
    return vec, known


def axpy(dst: dict, src: dict, f: Fraction):
    for k, v in src.items():
        nv = dst.get(k, Fraction(0)) + f * v
        if nv:
            dst[k] = nv
        else:
            dst.pop(k, None)


def eliminate(rows, col_pos):
    """Forward elimination; returns pivots {col: (vec, known)}."""
    pivots = {}
    for vec, known in rows:
        vec, known = dict(vec), dict(known)
        while vec:
            lead = min(vec, key=lambda c: col_pos[c])
            if lead in pivots:
                pvec, pknown = pivots[lead]
                f = vec[lead] / pvec[lead]
                axpy(vec, pvec, -f)
                axpy(known, pknown, -f)
            else:
                pivots[lead] = (vec, known)
                break
    return pivots


def reduce_vector(v: dict, pivots: dict, col_pos):
    """v = sum_r c_r A_r + d.  Returns (d, acc_known) with
    v.x = acc_known(values) + d.x  exactly.  Full reduction: repeat until
    no pivot-column entry remains (each step zeroes the earliest one and
    introduces only strictly later entries, so this terminates)."""
    d, acc = dict(v), {}
    while True:
        cands = [c for c in d if c in pivots]
        if not cands:
            break
        lead = min(cands, key=lambda c: col_pos[c])
        pvec, pknown = pivots[lead]
        f = d[lead] / pvec[lead]
        axpy(d, pvec, -f)
        axpy(acc, pknown, f)
    assert not any(c in pivots for c in d)
    return d, acc


# ===================================================================
# stage definitions
# ===================================================================
D2_ORDER = [(8, 2), (4, 6), (6, 4), (5, 5), (2, 8), (3, 7), (7, 3), (9, 1)]


def column_order(max_depth: int):
    """Pivot-preference order: deepest words first, then triples (table
    objects first), then doubles ending with the preferred surviving
    generators t2(7,3), t2(9,1)."""
    cols = []
    for dpt in range(max_depth, 3, -1):
        cols += sorted(w for w in comps(10) if len(w) == dpt)
    front3 = [(4, 3, 3), (4, 4, 2), (3, 4, 3), (4, 2, 4)]
    rest3 = sorted(w for w in comps(10)
                   if len(w) == 3 and w not in front3)
    cols += front3 + rest3
    cols += D2_ORDER
    return cols


def stage_products(stage: int):
    prods = []
    if stage <= 2:
        for a in range(2, 9):           # lam(a) * lam(10-a)
            b = 10 - a
            if 2 <= b and a <= b:
                prods.append(((a,), (b,)))
        for a in range(2, 9):           # lam(a) * t2(b,c)
            for u in comps(10 - a):
                if len(u) == 2:
                    prods.append(((a,), u))
        if stage == 2:
            for a in range(2, 9):       # lam(a) * t3(f,g,h)
                for u in comps(10 - a):
                    if len(u) == 3:
                        prods.append(((a,), u))
            for w1 in range(3, 8):      # t2 * t2 (unordered)
                w2 = 10 - w1
                for u in comps(w1):
                    if len(u) != 2:
                        continue
                    for v_ in comps(w2):
                        if len(v_) != 2:
                            continue
                        if w1 > w2 or (w1 == w2 and u > v_):
                            continue
                        prods.append((u, v_))
    else:                               # stage 3: full closure
        for w1 in range(2, 9):
            w2 = 10 - w1
            if w1 > w2:
                continue
            for u in comps(w1):
                for v_ in comps(w2):
                    if w1 == w2 and u > v_:
                        continue
                    prods.append((u, v_))
    return prods


def run_stage(stage: int, v_target: dict):
    max_depth = {1: 3, 2: 4, 3: 10}[stage]
    prods = stage_products(stage)
    rows = [build_relation(u, w) for u, w in prods]
    # sanity: support depth bound
    for (u, w), (vec, _) in zip(prods, rows):
        assert all(2 <= len(wd) <= max_depth for wd in vec), (u, w)
    cols = column_order(max_depth)
    col_pos = {c: i for i, c in enumerate(cols)}
    pivots = eliminate(rows, col_pos)
    rank = len(pivots)
    free_cols = [c for c in cols if c not in pivots
                 and any(c in vec for vec, _ in rows)]
    never_touched = [c for c in cols if c not in pivots
                     and not any(c in vec for vec, _ in rows)]
    d, acc = reduce_vector(v_target, pivots, col_pos)
    d3_pivots = [c for c in pivots if len(c) == 3]
    d2_pivots = [c for c in pivots if len(c) == 2]
    res = {
        "n_products": len(prods),
        "n_unknowns": len(cols),
        "rank": rank,
        "quotient_dim": len(cols) - rank,
        "pivots_depth3": len(d3_pivots),
        "pivots_depth2": len(d2_pivots),
        "free_columns_touched": [wlabel(c) for c in free_cols],
        "columns_never_in_any_relation": [wlabel(c) for c in never_touched],
        "v_in_rowspace": not d,
        "residual": {wlabel(k): str(val) for k, val in
                     sorted(d.items(), key=lambda kv: col_pos[kv[0]])},
        "known_combination": {klabel(k): str(val) for k, val in
                              sorted(acc.items(), key=lambda kv: klabel(kv[0]))},
    }
    # machine-readable copies for the numerics pass
    res["_residual_raw"] = [[list(k), str(val)] for k, val in d.items()]
    res["_known_raw"] = [[k[0], list(k[1]), list(k[2]) if k[0] == "prod"
                          else None, str(val)] for k, val in acc.items()]
    return res, pivots, col_pos


# ===================================================================
# symbolic main
# ===================================================================
def symbolic_main():
    out = json.loads(OUT_PATH.read_text()) if OUT_PATH.exists() else {}
    tab = json.loads((DATA / "s3_decomp_table.json").read_text())
    lin = tab["S3_linear_table"]

    def coeff(key):
        n, dnm = lin[key]
        return Fraction(int(n), int(dnm))

    v_target = {
        (4, 3, 3): coeff("('t3', 4, 3, 3)"),
        (4, 4, 2): coeff("('t3', 4, 4, 2)"),
        (8, 2): coeff("('t2', 8, 2)"),
        (4, 6): coeff("('t2', 4, 6)"),
    }
    assert v_target[(4, 3, 3)] == -2048 and v_target[(4, 4, 2)] == -1024
    assert v_target[(8, 2)] == -512 and v_target[(4, 6)] == -1536
    out["W10_definition"] = {wlabel(k): str(v) for k, v in v_target.items()}

    # cross-check the two known stuffles R7/R8 against qsh
    vec7, _ = build_relation((3,), (4, 3))
    assert vec7 == {(3, 4, 3): 1, (4, 3, 3): 2, (7, 3): 1, (4, 6): 1}, vec7
    vec8, _ = build_relation((4,), (4, 2))
    assert vec8 == {(4, 4, 2): 2, (4, 2, 4): 1, (8, 2): 1, (4, 6): 1}, vec8
    print("R7/R8 qsh cross-check: PASS")

    stage_piv = {}
    for stage in (1, 2, 3):
        t0 = time.time()
        res, pivots, col_pos = run_stage(stage, v_target)
        stage_piv[stage] = (pivots, col_pos)
        res["elapsed_sec"] = round(time.time() - t0, 2)
        out["stage%d" % stage] = res
        print("\n=== Stage %d ===" % stage)
        print("  products=%d unknowns=%d rank=%d quotient=%d (%.1fs)"
              % (res["n_products"], res["n_unknowns"], res["rank"],
                 res["quotient_dim"], res["elapsed_sec"]))
        print("  v in rowspace:", res["v_in_rowspace"])
        if res["residual"]:
            print("  residual: " + " + ".join(
                "(%s)*%s" % (c, l) for l, c in res["residual"].items()))
        if stage <= 2:
            print("  known combination (%d terms):" %
                  len(res["known_combination"]))
            for lbl, c in res["known_combination"].items():
                print("    (%s) * %s" % (c, lbl))
    # ---- canonical-form table for every depth<=3 word (stage-3 pivots,
    #      the strongest stuffle closure) + proportionality scan ----
    pivots3, col_pos3 = stage_piv[3]
    table = {}
    v_d, _ = reduce_vector(v_target, pivots3, col_pos3)
    for wd in sorted(comps(10)):
        if len(wd) not in (2, 3):
            continue
        d_w, acc_w = reduce_vector({wd: Fraction(1)}, pivots3, col_pos3)
        table[wlabel(wd)] = {
            "residual": {wlabel(k): str(c) for k, c in
                         sorted(d_w.items(),
                                key=lambda kv: col_pos3[kv[0]])},
            "knowns": {klabel(k): str(c) for k, c in acc_w.items()},
        }
    out["word_canonical_forms_stage3"] = table
    # is v's class proportional to a single word's class?
    prop_hits = []
    for wd in sorted(comps(10)):
        if len(wd) not in (2, 3):
            continue
        d_w, _ = reduce_vector({wd: Fraction(1)}, pivots3, col_pos3)
        if not d_w:
            continue
        # compare support and ratios with v_d
        if set(d_w) == set(v_d):
            ratios = {v_d[k] / d_w[k] for k in d_w}
            if len(ratios) == 1:
                prop_hits.append((wlabel(wd), str(ratios.pop())))
    out["v_class_proportional_to_single_word"] = prop_hits
    print("\nproportionality scan (v_class == alpha * word_class):",
          prop_hits if prop_hits else "NONE")

    # stage-1 relation list for the numerics verification pass
    out["stage1_relations"] = []
    for u, w in stage_products(1):
        vec, known = build_relation(u, w)
        out["stage1_relations"].append({
            "label": "%s*%s" % (wlabel(u), wlabel(w)),
            "u": list(u), "w": list(w),
            "vec": [[list(k), str(c)] for k, c in sorted(vec.items())],
            "known": [[k[0], list(k[1]), list(k[2]) if k[0] == "prod"
                       else None, str(c)] for k, c in known.items()],
        })
    OUT_PATH.write_text(json.dumps(out, indent=1))
    print("\nSaved:", OUT_PATH)


# ===================================================================
# numerics: value evaluators (Levin-safe; Abel for trailing 1)
# ===================================================================
def numerics_main(dps_low: int = 100, dps_high: int = 200):
    import mpmath

    shared = json.loads(SHARED_CACHE_PATH.read_text())
    vcache = (json.loads(VCACHE_PATH.read_text())
              if VCACHE_PATH.exists() else {})

    def vc_put(key, val, d):
        vcache["%s@%d" % (key, d)] = mpmath.nstr(val, d - 5)
        VCACHE_PATH.write_text(json.dumps(vcache, indent=1))

    def cache_lookup(key, d):
        for store in (vcache, shared):
            for dd in sorted(
                    {int(k.rsplit("@", 1)[1]) for k in store
                     if k.startswith(key + "@")}, reverse=True):
                if dd >= d:
                    return mpmath.mpf(store["%s@%d" % (key, dd)])
        return None

    def t2_raw(b, c):
        def term(j):
            j = int(j)
            return (mpmath.mpf(2 * j + 1) ** -c * mpmath.power(2, -b)
                    * mpmath.zeta(b, mpmath.mpf(j) + mpmath.mpf(3) / 2))
        return mpmath.nsum(term, [0, mpmath.inf], method='levin')

    def t3_raw(b1, b2, b3):
        assert b3 >= 2
        lam3 = (1 - mpmath.power(2, -b3)) * mpmath.zeta(b3)

        def term(j):
            j = int(j)
            tau = mpmath.power(2, -b1) * mpmath.zeta(
                b1, mpmath.mpf(j) + mpmath.mpf(3) / 2)
            pL = lam3 - mpmath.power(2, -b3) * mpmath.zeta(
                b3, mpmath.mpf(j) + mpmath.mpf(1) / 2)
            return mpmath.mpf(2 * j + 1) ** -b2 * tau * pL
        return mpmath.nsum(term, [1, mpmath.inf], method='levin')

    class Trailing1:
        """t3(b1,b2,1) via Abel summation (s3_closure_trailing1.py v3
        design: per-working-precision term arrays, Hurwitz forward
        recurrence, A(k) = A(1) - cumsum)."""
        NTERMS0 = 4096
        KMAX_HARD = 4_000_000

        def __init__(self, b1, b2):
            self.b1, self.b2 = b1, b2
            self._cache = {}

        def _a_direct(self, j):
            j = int(j)
            return (mpmath.mpf(2 * j + 1) ** -self.b2
                    * mpmath.power(2, -self.b1)
                    * mpmath.zeta(self.b1, mpmath.mpf(j) + mpmath.mpf(3) / 2))

        def _new_state(self):
            return {"A1": mpmath.nsum(self._a_direct, [1, mpmath.inf],
                                      method='levin'),
                    "z": mpmath.zeta(self.b1, mpmath.mpf(5) / 2),
                    "j_done": 1, "cum": mpmath.mpf(0), "terms": []}

        def _fill(self, st, k_target):
            if k_target > self.KMAX_HARD:
                raise RuntimeError("k=%d > hard cap" % k_target)
            terms = st["terms"]
            while len(terms) < k_target:
                k = len(terms) + 1
                A_k = st["A1"] - st["cum"]
                terms.append(A_k / (2 * k - 1))
                j = st["j_done"]
                a_j = (mpmath.mpf(2 * j + 1) ** -self.b2
                       * mpmath.power(2, -self.b1) * st["z"])
                st["cum"] += a_j
                st["z"] -= (mpmath.mpf(j) + mpmath.mpf(3) / 2) ** -self.b1
                st["j_done"] += 1

        def term(self, k):
            k = int(k)
            prec = mpmath.mp.prec
            st = self._cache.get(prec)
            if st is None:
                st = self._new_state()
                self._fill(st, self.NTERMS0)
                self._cache[prec] = st
            if k > len(st["terms"]):
                self._fill(st, 2 * k)
            return st["terms"][k - 1]

        def value(self):
            return mpmath.nsum(self.term, [1, mpmath.inf], method='levin')

    def val_word(word, d):
        word = tuple(word)
        if len(word) == 1:
            return (1 - mpmath.power(2, -word[0])) * mpmath.zeta(word[0])
        key = ("t2(%d,%d)" % word if len(word) == 2
               else "t3(%d,%d,%d)" % word)
        v = cache_lookup(key, d)
        if v is not None:
            return v
        t0 = time.time()
        if len(word) == 2:
            v = t2_raw(*word)
        elif len(word) == 3 and word[2] >= 2:
            v = t3_raw(*word)
        elif len(word) == 3:
            v = Trailing1(word[0], word[1]).value()
        else:
            raise ValueError("no evaluator for depth-%d word %s"
                             % (len(word), word))
        vc_put(key, v, d)
        print("    [%s@%d computed, %.0fs]" % (key, d, time.time() - t0),
              flush=True)
        return v

    out = json.loads(OUT_PATH.read_text())

    def known_value(item, d):
        kind, k1, k2, c = item
        c = Fraction(c)
        cv = mpmath.mpf(c.numerator) / c.denominator
        if kind == "prod":
            return cv * val_word(k1, d) * val_word(k2, d)
        return cv * val_word(k1, d)

    # ---- phase A: verify every stage-1 relation at dps_low ----
    dps = dps_low
    mpmath.mp.dps = dps + 20
    gate = mpmath.mpf(10) ** -(dps - 15)
    charter_gate = mpmath.mpf(10) ** -50
    print("=== Phase A: verifying all %d stage-1 stuffle relations at "
          "%d dps ===" % (len(out["stage1_relations"]), dps), flush=True)
    ver = {}
    worst = mpmath.mpf(0)
    for rel in out["stage1_relations"]:
        lhs = mpmath.mpf(0)
        for wd, c in rel["vec"]:
            c = Fraction(c)
            lhs += (mpmath.mpf(c.numerator) / c.denominator
                    * val_word(wd, dps))
        rhs = sum(known_value(item, dps) for item in rel["known"])
        r = abs(lhs - rhs)
        worst = max(worst, r)
        ver[rel["label"]] = mpmath.nstr(r, 4)
        status = "PASS" if r < gate else "FAIL"
        print("  %-22s residual %s  %s"
              % (rel["label"], mpmath.nstr(r, 4), status), flush=True)
    out["stage1_verification"] = {
        "dps": dps, "residuals": ver,
        "worst": mpmath.nstr(worst, 4),
        "all_pass_gate": bool(worst < gate),
        "all_pass_charter_50": bool(worst < charter_gate)}
    OUT_PATH.write_text(json.dumps(out, indent=1))

    # ---- phase B: the derived reduction identity at dps_high ----
    dps = dps_high
    mpmath.mp.dps = dps + 20
    gate = mpmath.mpf(10) ** -(dps - 15)
    Wdef = {(4, 3, 3): -2048, (4, 4, 2): -1024,
            (8, 2): -512, (4, 6): -1536}
    for stage in (1, 2):
        key = "stage%d" % stage
        st = out[key]
        print("=== Phase B: stage-%d derived identity at %d dps ==="
              % (stage, dps), flush=True)
        lhs = mpmath.mpf(0)
        for wd, c in Wdef.items():
            lhs += c * val_word(wd, dps)
        rhs = mpmath.mpf(0)
        for item in st["_known_raw"]:
            rhs += known_value(item, dps)
        for wd, c in st["_residual_raw"]:
            c = Fraction(c)
            rhs += (mpmath.mpf(c.numerator) / c.denominator
                    * val_word(wd, dps))
        r = abs(lhs - rhs)
        status = "PASS" if r < gate else "FAIL"
        print("  W10 minus reduction: residual %s  %s"
              % (mpmath.nstr(r, 4), status), flush=True)
        st["identity_verification"] = {
            "dps": dps, "residual": mpmath.nstr(r, 4),
            "pass_gate": bool(r < gate),
            "pass_charter_50": bool(r < charter_gate)}
        OUT_PATH.write_text(json.dumps(out, indent=1))

    # ---- phase B2: antisymmetric-basis presentation at dps_high ----
    # Exact rewrite of the canonical reduction using the Type-A stuffles
    #   t2(4,4) = (lam4^2-lam8)/2,  t2(2,4) = (lam2 lam4 - lam6 - a42)/2,
    #   t2(6,4) = (lam4 lam6 - lam10 + a64)/2,
    #   t2(2,8) = (lam2 lam8 - lam10 - a82)/2,
    #   t2(7,3) = (lam3 lam7 - lam10 + a73)/2,
    # with a(b,c) := t2(b,c) - t2(c,b):
    #   W10/256 = 4 lam10 - lam2 lam8 - 3 lam4 lam6
    #             - 4 lam3 a(4,3) - 2 lam4 a(4,2)
    #             + a(8,2) + 4 a(7,3) - 3 a(6,4)
    #             - 4 t3(2,4,4) - 8 t3(3,3,4)
    print("=== Phase B2: antisymmetric-basis presentation at %d dps ==="
          % dps, flush=True)

    def lamv(s):
        return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)

    def a(b, c):
        return val_word((b, c), dps) - val_word((c, b), dps)

    lhs = mpmath.mpf(0)
    for wd, c in Wdef.items():
        lhs += c * val_word(wd, dps)
    rhs = 256 * (4 * lamv(10) - lamv(2) * lamv(8) - 3 * lamv(4) * lamv(6)
                 - 4 * lamv(3) * a(4, 3) - 2 * lamv(4) * a(4, 2)
                 + a(8, 2) + 4 * a(7, 3) - 3 * a(6, 4)
                 - 4 * val_word((2, 4, 4), dps)
                 - 8 * val_word((3, 3, 4), dps))
    r = abs(lhs - rhs)
    status = "PASS" if r < gate else "FAIL"
    print("  antisym presentation residual %s  %s"
          % (mpmath.nstr(r, 4), status), flush=True)
    out["antisym_presentation"] = {
        "dps": dps,
        "identity": ("W10/256 = 4*lam(10) - lam(2)*lam(8) - 3*lam(4)*lam(6)"
                     " - 4*lam(3)*a(4,3) - 2*lam(4)*a(4,2)"
                     " + a(8,2) + 4*a(7,3) - 3*a(6,4)"
                     " - 4*t3(2,4,4) - 8*t3(3,3,4);  a(b,c):=t2(b,c)-t2(c,b)"),
        "residual": mpmath.nstr(r, 4),
        "pass_gate": bool(r < gate),
        "pass_charter_50": bool(r < mpmath.mpf(10) ** -50)}
    OUT_PATH.write_text(json.dumps(out, indent=1))
    print("Saved:", OUT_PATH, flush=True)


def assembly_check_main(dps: int = 200):
    """Corollary in the symbolic-assembly track's normalization: its W10
    remainder is W10_assembly = W10 - 1024 t2(7,3).  Substituting
    t2(7,3) = (lam3 lam7 - lam10 + a(7,3))/2 (exact Type-A stuffle) into
    the antisym presentation gives
      W10_assembly/256 = 6 lam10 - lam2 lam8 - 2 lam3 lam7 - 3 lam4 lam6
                         - 4 lam3 a(4,3) - 2 lam4 a(4,2)
                         + a(8,2) + 2 a(7,3) - 3 a(6,4)
                         - 4 t3(2,4,4) - 8 t3(3,3,4).
    Verified here numerically (all values should be cache hits)."""
    import mpmath
    mpmath.mp.dps = dps + 20
    shared = json.loads(SHARED_CACHE_PATH.read_text())
    vcache = json.loads(VCACHE_PATH.read_text())

    def look(key):
        for store in (vcache, shared):
            cands = sorted((int(k.rsplit("@", 1)[1]) for k in store
                            if k.startswith(key + "@")), reverse=True)
            for dd in cands:
                if dd >= dps:
                    return mpmath.mpf(store["%s@%d" % (key, dd)])
        raise KeyError(key)

    def t2(b, c):
        return look("t2(%d,%d)" % (b, c))

    def t3(a_, b, c):
        return look("t3(%d,%d,%d)" % (a_, b, c))

    def lam(s):
        return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)

    def a(b, c):
        return t2(b, c) - t2(c, b)

    lhs = (-2048 * t3(4, 3, 3) - 1024 * t3(4, 4, 2) - 512 * t2(8, 2)
           - 1536 * t2(4, 6) - 1024 * t2(7, 3))
    rhs = 256 * (6 * lam(10) - lam(2) * lam(8) - 2 * lam(3) * lam(7)
                 - 3 * lam(4) * lam(6)
                 - 4 * lam(3) * a(4, 3) - 2 * lam(4) * a(4, 2)
                 + a(8, 2) + 2 * a(7, 3) - 3 * a(6, 4)
                 - 4 * t3(2, 4, 4) - 8 * t3(3, 3, 4))
    r = abs(lhs - rhs)
    ok = r < mpmath.mpf(10) ** -(dps - 15)
    print("assembly-normalized W10 identity residual %s  %s"
          % (mpmath.nstr(r, 4), "PASS" if ok else "FAIL"))
    out = json.loads(OUT_PATH.read_text())
    out["assembly_normalized_identity"] = {
        "dps": dps,
        "identity": ("(W10 - 1024 t2(7,3))/256 = 6 lam(10) - lam(2)lam(8)"
                     " - 2 lam(3)lam(7) - 3 lam(4)lam(6) - 4 lam(3)a(4,3)"
                     " - 2 lam(4)a(4,2) + a(8,2) + 2 a(7,3) - 3 a(6,4)"
                     " - 4 t3(2,4,4) - 8 t3(3,3,4)"),
        "residual": mpmath.nstr(r, 4), "pass": bool(ok)}
    OUT_PATH.write_text(json.dumps(out, indent=1))
    print("Saved:", OUT_PATH)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "symbolic"
    if mode == "symbolic":
        symbolic_main()
    elif mode == "assembly_check":
        assembly_check_main(int(sys.argv[2]) if len(sys.argv) > 2 else 200)
    else:
        numerics_main(int(sys.argv[2]) if len(sys.argv) > 2 else 100,
                      int(sys.argv[3]) if len(sys.argv) > 3 else 200)
