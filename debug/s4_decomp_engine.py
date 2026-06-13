"""S^(4) stage-1: generic exact decomposition engine (k=4 o-space cores).

Decomposes the o-space cores of S^(4) into Hoffman multiple t-values, by a
GENERIC region algebra (replacing the hand-derived 13-region split of the
k=3 engine `s3_decomp_engine.py`, which is re-derived here as a regression).

ORDERING CONVENTION (identical to s3_decomp_engine.py):
    t_m(b1,..,bm) = sum_{o_m > ... > o_1 >= 1, odd} o_m^(-b1) ... o_1^(-bm)
    lam(s)        = sum_{o >= 1 odd} o^(-s)
i.e. the FIRST argument is carried by the LARGEST summation variable.

Objects (all direct sums over INDEPENDENT o_i >= 5 odd, m_ij = min(o_i,o_j),
psi(o) = 8(o^2-1)/o^4):
    C4  = sum m12 m23 m34 psi^4        (4-chain; the only depth-4 source)
    Ce  = sum m12 m23 psi1 psi2 psi3^2
    Cm  = sum m12 m23 psi1 psi2^2 psi3
    H31 = sum m12 psi1 psi2^3
    H22 = sum m12 psi1^2 psi2^2
    C   = sum m12 m23 psi1 psi2 psi3   (k=3 regression)
    G   = sum m12 psi1 psi2^2          (k=3 regression)
    R4o = sum psi^4                     (single Hurwitz row)

Region algebra: for a path of r variables, enumerate all ordered set
partitions (weak orderings; Fubini(3)=13, Fubini(4)=75).  In each region
every min resolves to the smaller class; tied variables share one value.
Each class contributes its Laurent-monomial choices; classes with inverse
exponent 0 (from o^2*psi on a double-min class) are eliminated by exact
odd-counting rules:
    [innermost arg 0]  t(..,b,0)   = (1/2) t(..,b-1) - (1/2) t(..,b)
    [middle arg 0]     t(..a,0,c..) = (1/2) t(..a-1,c..) - (1/2) t(..a,c-1..)
                                      - t(..a,c..)
(count of odds below v is (v-1)/2; strictly between u<v is (v-u)/2 - 1),
applied recursively.  A zero on the largest class never occurs (asserted).

o >= 5 from o >= 1 by inclusion-exclusion in psi(3)-pins (psi(1) = 0): pin
set S at o=3 contributes (-1)^|S| psi(3)^... * 3^{#edges meeting S} times
the product over connected free components of their own o>=1 tables
(min(3,o)=3 on the psi-support).

Verification: every table is checked BIT-EXACTLY in Fraction arithmetic
against the direct o>=5 sum at truncations O = 25 and 41 (symbols truncated
to the same box), plus value-level regression of C and G against the k=3
table `debug/data/s3_decomp_table.json`.

Output: debug/data/s4_decomp_tables.json  (+ census of t4 content)
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from itertools import product as iproduct
from pathlib import Path

PSI = {2: Fraction(8), 4: Fraction(-8)}          # psi(o)  = 8/o^2 - 8/o^4


def lpow(w: dict, k: int) -> dict:
    """k-th power of a Laurent dict {inv_exponent: coeff}."""
    out = {0: Fraction(1)}
    for _ in range(k):
        nxt = {}
        for e1, c1 in out.items():
            for e2, c2 in w.items():
                nxt[e1 + e2] = nxt.get(e1 + e2, Fraction(0)) + c1 * c2
        out = nxt
    return {e: c for e, c in out.items() if c != 0}


def lmul(a: dict, b: dict) -> dict:
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            out[e1 + e2] = out.get(e1 + e2, Fraction(0)) + c1 * c2
    return {e: c for e, c in out.items() if c != 0}


def leval(w: dict, o: int) -> Fraction:
    return sum(c * Fraction(1, o ** e) if e >= 0 else c * o ** (-e)
               for e, c in w.items())


# ---------------------------------------------------------------------------
# symbols: ('lam', s) | ('t', (b1..bm)) | ('one',) | ('prod', sym, sym, ...)
# tables: dict symbol -> Fraction
# ---------------------------------------------------------------------------

def tab_add(dst: dict, sym, c: Fraction) -> None:
    if c == 0:
        return
    dst[sym] = dst.get(sym, Fraction(0)) + c
    if dst[sym] == 0:
        del dst[sym]


def tab_merge(dst: dict, src: dict, scale: Fraction = Fraction(1)) -> None:
    for k, v in src.items():
        tab_add(dst, k, v * scale)


def eliminate_zeros(args: tuple, coeff: Fraction, out: dict) -> None:
    """Register t(args) after recursively eliminating zero exponents."""
    if 0 not in args:
        assert all(b >= 1 for b in args), args
        assert args[0] >= 2, ("divergent leading exponent", args)
        if len(args) == 1:
            tab_add(out, ('lam', args[0]), coeff)
        else:
            tab_add(out, ('t', args), coeff)
        return
    i = args.index(0)          # first zero, position i (0-based)
    m = len(args)
    assert i > 0, ("zero exponent on largest class", args)
    rest = args[:i] + args[i + 1:]
    if i == m - 1:             # innermost class
        up = list(rest)
        up[i - 1] -= 1
        eliminate_zeros(tuple(up), coeff / 2, out)
        eliminate_zeros(tuple(rest), -coeff / 2, out)
    else:                      # middle class: neighbors at i-1 (upper), i (lower in rest)
        up = list(rest)
        up[i - 1] -= 1
        eliminate_zeros(tuple(up), coeff / 2, out)
        lo = list(rest)
        lo[i] -= 1
        eliminate_zeros(tuple(lo), -coeff / 2, out)
        eliminate_zeros(tuple(rest), -coeff, out)


def weak_orderings(r: int):
    """All surjections {0..r-1} -> {0..k-1}, k = 1..r (class 0 = smallest)."""
    for f in iproduct(range(r), repeat=r):
        k = max(f) + 1
        if set(f) == set(range(k)):
            yield f, k


def decompose_o1(weights: list, edges: list) -> dict:
    """Exact table for sum over o_i >= 1 odd of prod_e min * prod_i w_i."""
    r = len(weights)
    out: dict = {}
    for f, k in weak_orderings(r):
        cls_w = [dict(d) for d in [{0: Fraction(1)}] * k]
        cls_w = [{0: Fraction(1)} for _ in range(k)]
        for i, w in enumerate(weights):
            cls_w[f[i]] = lmul(cls_w[f[i]], w)
        nmin = [0] * k
        for (u, v) in edges:
            nmin[min(f[u], f[v])] += 1
        # class j Laurent poly shifted by -nmin[j] (min factor = +1 power)
        polys = []
        for j in range(k):
            polys.append({e - nmin[j]: c for e, c in cls_w[j].items()})
        for combo in iproduct(*[sorted(p.items()) for p in polys]):
            coeff = Fraction(1)
            exps_asc = []
            for (e, c) in combo:
                coeff *= c
                exps_asc.append(e)
            args = tuple(reversed(exps_asc))   # b1 on largest class
            eliminate_zeros(args, coeff, out)
    return out


def decompose_o5(weights: list, edges: list) -> dict:
    """Exact table for the same object with all o_i >= 5, via psi(3)-pins."""
    r = len(weights)
    out: dict = {}
    for pins in iproduct((False, True), repeat=r):
        S = [i for i in range(r) if pins[i]]
        sgn = Fraction(-1) ** len(S)
        scal = Fraction(1)
        for i in S:
            scal *= leval(weights[i], 3)
        if scal == 0:
            continue
        free = [i for i in range(r) if not pins[i]]
        fedges = [(u, v) for (u, v) in edges if u in free and v in free]
        scal *= Fraction(3) ** (len(edges) - len(fedges))
        # connected components of the free subgraph
        comp_of = {i: i for i in free}

        def find(x):
            while comp_of[x] != x:
                comp_of[x] = comp_of[comp_of[x]]
                x = comp_of[x]
            return x

        for (u, v) in fedges:
            ru, rv = find(u), find(v)
            if ru != rv:
                comp_of[max(ru, rv)] = min(ru, rv)
        comps: dict = {}
        for i in free:
            comps.setdefault(find(i), []).append(i)
        factor_tables = []
        for root, mem in sorted(comps.items()):
            sub_w = [weights[i] for i in mem]
            idx = {i: j for j, i in enumerate(mem)}
            sub_e = [(idx[u], idx[v]) for (u, v) in fedges
                     if u in idx and v in idx]
            factor_tables.append(decompose_o1(sub_w, sub_e))
        if not factor_tables:
            tab_add(out, ('one',), sgn * scal)
            continue
        # expand the product of linear tables
        acc = {('one',): sgn * scal}
        for ft in factor_tables:
            nxt = {}
            for sa, ca in acc.items():
                for sb, cb in ft.items():
                    if sa == ('one',):
                        sym = sb
                    elif sa[0] == 'prod':
                        sym = ('prod',) + tuple(sorted(sa[1:] + (sb,),
                                                       key=repr))
                    else:
                        sym = ('prod',) + tuple(sorted((sa, sb), key=repr))
                    tab_add(nxt, sym, ca * cb)
            acc = nxt
        tab_merge(out, acc)
    return out


# ---------------------------------------------------------------------------
# exact truncated evaluation (box o <= O, odd)
# ---------------------------------------------------------------------------

def eval_sym(sym, O: int, cache: dict) -> Fraction:
    if sym in cache:
        return cache[sym]
    odds = range(1, O + 1, 2)
    if sym == ('one',):
        v = Fraction(1)
    elif sym[0] == 'lam':
        v = sum(Fraction(1, o ** sym[1]) for o in odds)
    elif sym[0] == 'prod':
        v = Fraction(1)
        for s in sym[1:]:
            v *= eval_sym(s, O, cache)
    elif sym[0] == 't':
        args = sym[1]
        m = len(args)
        # accumulators A_d(o) = sum over strictly increasing chains below o
        acc = [Fraction(0)] * m   # acc[d] = A_{d+1} at current frontier
        prev = None
        tot = Fraction(0)
        # args reversed: innermost exponent is args[-1]
        rev = tuple(reversed(args))   # rev[0] innermost ... rev[m-1] outermost
        A = [Fraction(0)] * (m - 1)
        for o in odds:
            if prev is not None:
                newA = list(A)
                # extend chains by prev as the d-th smallest
                p = Fraction(prev)
                newA[0] = A[0] + Fraction(1, prev ** rev[0])
                for d in range(1, m - 1):
                    newA[d] = A[d] + Fraction(1, prev ** rev[d]) * A[d - 1]
                A = newA
            tot += Fraction(1, o ** rev[m - 1]) * A[m - 2]
            prev = o
        v = tot
    else:
        raise ValueError(sym)
    cache[sym] = v
    return v


def eval_table(tab: dict, O: int, cache: dict) -> Fraction:
    return sum(c * eval_sym(s, O, cache) for s, c in tab.items())


def direct_o5(weights: list, edges: list, O: int) -> Fraction:
    odds = [o for o in range(5, O + 1, 2)]
    vals = {o: [leval(w, o) for w in weights] for o in odds}
    r = len(weights)
    tot = Fraction(0)
    for combo in iproduct(odds, repeat=r):
        t = Fraction(1)
        for (u, v) in edges:
            t *= min(combo[u], combo[v])
        for i in range(r):
            t *= vals[combo[i]][i]
        tot += t
    return tot


# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.time()
    PSI2, PSI3_, PSI4 = lpow(PSI, 2), lpow(PSI, 3), lpow(PSI, 4)
    CHAIN3 = [(0, 1), (1, 2)]
    CHAIN4 = [(0, 1), (1, 2), (2, 3)]
    PAIR = [(0, 1)]

    OBJECTS = {
        'C':   ([PSI, PSI, PSI], CHAIN3),
        'G':   ([PSI, PSI2], PAIR),
        'C4':  ([PSI, PSI, PSI, PSI], CHAIN4),
        'Ce':  ([PSI, PSI, PSI2], CHAIN3),
        'Cm':  ([PSI, PSI2, PSI], CHAIN3),
        'H31': ([PSI, PSI3_], PAIR),
        'H22': ([PSI2, PSI2], PAIR),
        'R4o': ([PSI4], []),
        'Smin_o': ([PSI, PSI], PAIR),
    }

    OUT = {"convention": "t(b1..bm)=sum_{o_m>..>o_1>=1 odd} o_m^-b1..o_1^-bm;"
                         " first arg = largest variable"}
    tables = {}
    ok_all = True
    for name, (w, e) in OBJECTS.items():
        t1 = time.time()
        tab = decompose_o5(w, e)
        tables[name] = tab
        checks = {}
        for O in (25, 41):
            cache: dict = {}
            d = direct_o5(w, e, O)
            v = eval_table(tab, O, cache)
            checks[O] = (d == v)
            ok_all &= (d == v)
        print("%-4s table %4d symbols   O=25:%s O=41:%s   (%.1fs)"
              % (name, len(tab),
                 "EXACT" if checks[25] else "FAIL",
                 "EXACT" if checks[41] else "FAIL", time.time() - t1))
        OUT["verify_%s" % name] = {str(k): v for k, v in checks.items()}

    # regression: C and G values vs the k=3 table (value-level, O=25/41)
    k3 = json.loads((Path(__file__).parent / "data"
                     / "s3_decomp_table.json").read_text())

    def parse_k3(table_json):
        tab = {}
        for ks, (nu, de) in table_json.items():
            sym = eval(ks)  # noqa: S307 - trusted repo artifact
            if sym[0] == 't3':
                sym2 = ('t', (sym[1], sym[2], sym[3]))
            elif sym[0] == 't2':
                sym2 = ('t', (sym[1], sym[2]))
            elif sym[0] == 'lam':
                sym2 = ('lam', sym[1])
            elif sym[0] == 'lam2':
                sym2 = ('prod', ('lam', sym[1]), ('lam', sym[2]))
            elif sym[0] == 'one':
                sym2 = ('one',)
            else:
                raise ValueError(sym)
            tab[sym2] = tab.get(sym2, Fraction(0)) + Fraction(int(nu), int(de))
        return tab

    reg_ok = True
    for nm, key in (('C', 'C_table'), ('G', 'G_table')):
        old = parse_k3(k3[key])
        for O in (25, 41):
            cache = {}
            a = eval_table(tables[nm], O, cache)
            b = eval_table(old, O, cache)
            reg_ok &= (a == b)
        print("regression %s vs k=3 table: %s" % (nm, "EXACT" if reg_ok else "FAIL"))
    OUT["regression_k3_C_G"] = reg_ok
    ok_all &= reg_ok

    # ---------------- census ----------------
    def census(tab):
        t4 = sorted([s[1] for s in tab if s[0] == 't' and len(s[1]) == 4])
        t3 = sorted([s[1] for s in tab if s[0] == 't' and len(s[1]) == 3])
        t2 = sorted([s[1] for s in tab if s[0] == 't' and len(s[1]) == 2])
        return t4, t3, t2

    all_t4, all_t3, all_t2 = set(), set(), set()
    for name in ('C4', 'Ce', 'Cm', 'H31', 'H22'):
        t4, t3, t2 = census(tables[name])
        all_t4 |= set(t4)
        all_t3 |= set(t3)
        all_t2 |= set(t2)
    tr1 = sorted(a for a in all_t4 if a[3] == 1 and a[2] != 1)
    tr2 = sorted(a for a in all_t4 if a[3] == 1 and a[2] == 1 and a[1] != 1)
    tr3 = sorted(a for a in all_t4 if a[1:] == (1, 1, 1))
    wts = sorted({sum(a) for a in all_t4})
    print("\nCENSUS (k=4 cores):")
    print("  distinct t4: %d   t3: %d   t2: %d"
          % (len(all_t4), len(all_t3), len(all_t2)))
    print("  t4 weights realized: %s  (ceiling 13 predicted)" % wts)
    print("  t4 single-trailing-1  (.,.,b,1): %d   %s" % (len(tr1), tr1))
    print("  t4 double-trailing-1 (.,b,1,1): %d   %s" % (len(tr2), tr2))
    print("  t4 TRIPLE-trailing-1 (b,1,1,1): %d   %s" % (len(tr3), tr3))
    OUT["census"] = {
        "t4_count": len(all_t4), "t3_count": len(all_t3),
        "t2_count": len(all_t2),
        "t4_weights": wts,
        "t4_list": sorted(map(list, all_t4)),
        "t4_single_trailing1": sorted(map(list, tr1)),
        "t4_double_trailing1": sorted(map(list, tr2)),
        "t4_triple_trailing1": sorted(map(list, tr3)),
        "t3_list": sorted(map(list, all_t3)),
        "t2_list": sorted(map(list, all_t2)),
    }

    def serialize(tab):
        return {repr(k): [str(v.numerator), str(v.denominator)]
                for k, v in sorted(tab.items(), key=lambda kv: repr(kv[0]))}

    OUT["tables"] = {nm: serialize(tb) for nm, tb in tables.items()}
    OUT["all_pass"] = ok_all
    outp = Path(__file__).parent / "data" / "s4_decomp_tables.json"
    outp.write_text(json.dumps(OUT, indent=2))
    print("\nALL EXACT: %s   saved %s  (%.1fs)"
          % (ok_all, outp.name, time.time() - t0))


if __name__ == "__main__":
    main()
