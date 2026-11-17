"""Record FX: canonical basis for the depth-2 named bases of debug/fx_reducer.py, and its controls.
(Design: debug/fx_design_memo.md sections 7-8.)  Computes NO scored quantity on M1 / M2r: there only booleans.

WHY.  The reducer names its irreducible residues  B(lam, E, p, [A, B]) = int_A^B e^{lam r} E(r)/(r - p) dr,  E a
product of one or two factors E1(c r + d).  These names are not a free basis.  The relations used here:

IDENTITY (product rule).  Let F be a multiset of factors f = (c_f, d_f), x_f = c_f r + d_f, q_f = -d_f/c_f.  Since
E1'(x) = -e^{-x}/x  and  c_f e^{-x_f}/x_f = e^{-d_f} e^{-c_f r}/(r - q_f),

    d/dr prod_{f in F} E1(x_f) = - sum_{f in F} e^{-d_f} e^{-c_f r}/(r - q_f) prod_{g in F - f} E1(x_g),

and integrating over a panel [A, B] on which no x_f vanishes in the open interval:

    sum_{f distinct} m_f e^{-d_f} B(-c_f, F - f, q_f, [A, B]) = -[prod_F E1(x_f)]_A^B          (I_F)

(m_f the multiplicity of f in F).  Validity: (i) q_f not in the open panel, nor at a finite endpoint other than a
Hadamard finite-part endpoint A = 0 = q_f; (ii) at B = infinity every c_f > 0 (the product vanishes there);
(iii) at A = 0 both sides are Hadamard finite parts (every ln-epsilon power dropped), which is linear and
commutes with the fundamental theorem of calculus: FP int_eps^B g' = g(B) - [constant term of g(eps), ln eps
dropped] -- exactly the reducer's convention for the B atoms (base_value) and for its endpoint terms
(Reducer._boundary_fp).  The identity is checked numerically for every natural set on every panel (control A1).
Every B(lam, E, p) with lam != 0 lies in exactly one identity set: F = E + {f1}, f1 = (-lam, lam p) (then
c_f1 = -lam, q_f1 = p);  the map is an involution (each member determines the same F).

CANONICAL FORM (route (a) of the brief).  A set I_F is APPLIED only if every factor of F is natural for the panel
(occurs among the factors the engine puts on that panel) -- a set with a non-natural factor has a member that never
occurs, so it cannot relate occurring bases.  In an applied set the member with the largest repr is eliminated in
favour of the others plus endpoint products (depth 0, na_exact class).  One pass is enough (the other members are
canonical: their set is the same F).  lam = 0 bases have no identity set (f1 would be E1(0)).
What it does NOT cover: relations outside the product-rule family (dilogarithm-type functional equations of
half-line integrals, unknown period relations).  Those are probed numerically by route (b) below, at the basis
level, on an unscored geometry (controls B1-B4).

ROUTE (b) (numerical, controls only).  Integer relations among a set of constants by LLL on [I | round(W v)]
(flint), W = 10^(dps - 20); a relation is accepted when its coefficients are <= H and it vanishes to 10^-(dps-40)
when re-evaluated at dps + 60.  At the level of F_L itself this route is out of reach: F_L's coefficients reach
~1e130 at L = 30, so a relation would need > n * 130 digits.  It is therefore used on the BASIS (small-height
relations among canonical bases and endpoint products), which is the question (a) answers exactly for the known
family.

CLI
  python debug/fx_indep.py controls        -> debug/data/fx_indep_controls.json (A1-A5, B1-B3, B4, timings)
  python debug/fx_indep.py stage <name>    one stage (A1 A2 A3 A4 A5 B B4)
B4 (stage_B4): LLL on the WHOLE canonical constant set of C1 (L = 2..6, 77 constants), 450 and 550 digits,
|coef| <= 1e4.  Negative: no relation.  Positive: one eliminated identity member is added, together with the
identity's own member and endpoint keys scaled by the e^{e0} of an occurring key (so the identity lies in the span);
the relation found must involve the eliminated member.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import flint
import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import na_exact as NE          # noqa: E402
import fx_reducer as X         # noqa: E402
import fx_exact as FX          # noqa: E402

# unscored control geometries on the matched line (kappa rational), runtime-only additions to FX.CASES
CONTROL_CASES = {
    'C1': dict(a='1', b='9/4', D1='7/5', zX='13/10', zY='1', zZ='3/2', zW='1'),      # a < D1 < b, as M1
    'C2': dict(a='1/4', b='1', D1='7/5', zX='13/10', zY='1', zZ='2', zW='1'),        # a < b < D1, as M2r
}
FX.CASES.update(CONTROL_CASES)
OUTF = HERE / 'data' / 'fx_indep_controls.json'
P_AUX = 1073741441


# ---------------------------------------------------------------------------------------------- natural factors
def natural_factors(case, Ls=range(2, 7)):
    """{panel: set of factors (c, d)} put on the panel by the engine (rho2 E1 factor + psi E1 factors), and whether
    the set is the same for every L in Ls."""
    K = NE.Field(P_AUX)
    eng = FX.Engine(case, K)
    Lm = max(Ls)
    it, kt = NE.bessel_tables(K, Lm + 2)
    Qs = list(eng._Q_steps(Lm)) if eng.matched else None
    per = {}
    for L in Ls:
        for panel in eng.panels:
            fs = set()
            fts = eng.f_terms_matched(L, Qs[L], panel) if eng.matched else eng.f_terms_equal(L, panel, it, kt)
            for t in fts:
                fs.update(t[6])
            pts, _, _ = eng.psi(L, panel, it, kt)
            for t in pts:
                fs.update(t[3])
            per.setdefault(panel, []).append(frozenset(fs))
    stable = all(len(set(v)) == 1 for v in per.values())
    return {panel: set().union(*v) for panel, v in per.items()}, stable


# ---------------------------------------------------------------------------------------------- identity sets
def _valid(F, A, B):
    for c, d in F:
        if c == 0:
            return False
        q = Fr(-d) / c
        if A < q and (B is None or q < B):
            return False
        if q == A and A != 0:
            return False
        if B is not None and q == B:
            return False
        if B is None and c <= 0:
            return False
    return True


def identity_set(lam, E, p, A, B, Fnat=None):
    """Members {f: ((lam_f, E_f, p_f), m_f)} of the product-rule set containing B(lam, E, p, [A, B]), with F, or None
    (lam = 0, a non-natural factor when Fnat is given, or the identity invalid on the panel)."""
    lam, p = Fr(lam), Fr(p)
    if lam == 0:
        return None
    f1 = (-lam, lam * p)
    F = [tuple(Fr(v) for v in f) for f in E] + [f1]
    if Fnat is not None and any(f not in Fnat for f in F):
        return None
    if not _valid(F, A, B):
        return None
    mem = {}
    for f in set(F):
        rest = list(F)
        rest.remove(f)
        mem[f] = ((-f[0], tuple(sorted(rest)), Fr(-f[1]) / f[0]), F.count(f))
    return F, mem


def _bnd(R, K, F, A, B, ck):
    """{CK: coeff} of [prod_F E1(x_f)]_A^B (finite part at A = 0), ck multiplied in."""
    red = X.Reducer(R, A, B)
    one = R.poly([K(1)])
    red._boundary(Fr(0), tuple(sorted(F)), ck, one, None, red.B, 1)
    red._boundary(Fr(0), tuple(sorted(F)), ck, one, None, red.A, -1)
    return red.out


def _srt(atoms):
    return tuple(sorted(atoms, key=repr))


def canonical(res, Fnat, K, stats=None):
    """Rewrite {CK: q} (fx_exact output) to the canonical basis.  Fnat: {(A, B): natural factor set}."""
    R = X.Ring(K)
    out = {}
    for (e0, atoms), v in res.items():
        Bs = [a for a in atoms if a[0] == 'B']
        if not Bs:
            X.add(out, (e0, atoms), v)
            continue
        assert len(Bs) == 1, atoms
        b = Bs[0]
        i = atoms.index(b)
        rest = atoms[:i] + atoms[i + 1:]
        _, lam, E, p, A, Bend = b
        s = identity_set(lam, E, p, A, Bend, Fnat.get((A, Bend)))
        me = (lam, E, p)
        if s is None:
            X.add(out, (e0, atoms), v)
            continue
        F, mem = s
        tgt_f = max(mem, key=lambda f: repr(mem[f][0]))
        if mem[tgt_f][0] != me:
            X.add(out, (e0, atoms), v)
            continue
        if stats is not None:
            stats['eliminated'] = stats.get('eliminated', 0) + 1
        mt, dt = K(mem[tgt_f][1]), tgt_f[1]
        for f, (key, m) in mem.items():
            if f == tgt_f:
                continue
            nb = ('B', key[0], key[1], key[2], A, Bend)
            X.add(out, (e0 + dt - f[1], _srt(rest + (nb,))), -v * K(m) / mt)
        for k, w in _bnd(R, K, F, A, Bend, (e0 + dt, rest)).items():
            X.add(out, k, -v * w / mt)
    return X.canon({k: v for k, v in out.items() if X.nz(v)}, K)


def depth2_keys(res):
    return [k for k in res if any(a[0] == 'B' for a in k[1])]


# ---------------------------------------------------------------------------------------------- route (b): LLL
def int_relation(vals, dps, H, check=None):
    """Smallest-norm integer relation among vals (mp numbers at >= dps digits) with max |coef| <= H, or None.
    check: callable giving the values at higher precision, to confirm a candidate."""
    n = len(vals)
    W = mp.mpf(10) ** (dps - 20)
    rows = []
    for i, v in enumerate(vals):
        r = [0] * (n + 1)
        r[i] = 1
        r[n] = int(mp.nint(W * v))
        rows.append(r)
    M = flint.fmpz_mat(rows).lll()
    best = None
    for i in range(n):
        c = [int(M[i, j]) for j in range(n)]
        if not any(c) or max(abs(x) for x in c) > H:
            continue
        if abs(sum(ci * vi for ci, vi in zip(c, vals))) > mp.mpf(10) ** (-(dps - 40)):
            continue
        if check is not None:
            hv = check()
            if abs(sum(ci * vi for ci, vi in zip(c, hv))) > mp.mpf(10) ** (-(dps - 40)):
                continue
        if best is None or max(map(abs, c)) < max(map(abs, best)):
            best = c
    return best


def key_value(k, cache):
    e0, atoms = k
    t = mp.exp(X._mpq(mp, e0))
    for a in atoms:
        t *= X.atom_value(mp, a, cache)
    return t


# ---------------------------------------------------------------------------------------------- controls
OUT = json.loads(OUTF.read_text()) if OUTF.exists() else {}


def save(k, v):
    OUT[k] = v
    OUTF.write_text(json.dumps(OUT, indent=1, default=str))
    print(k, json.dumps(v, default=str)[:1500], flush=True)


def _natural_sets(Fn, A, B, sizes=(2,)):
    """All valid product-rule sets F (size 2, or 3) from the natural factors of a panel; sets of size 2 include squares."""
    fl = sorted(Fn)
    out = []
    for i, f in enumerate(fl):
        for g in fl[i:]:
            if 2 in sizes and _valid([f, g], A, B):
                out.append([f, g])
            if 3 in sizes:
                for h in fl[fl.index(g):]:
                    if _valid([f, g, h], A, B):
                        out.append([f, g, h])
    return out


def _identity_residual(F, A, B, dps):
    """|sum m e^{-d} B + [prod E1]| / max term, at dps."""
    K = NE.Field(None)
    R = X.Ring(K)
    mp.mp.dps = dps
    cache = {}
    terms = []
    for f in set(F):
        rest = list(F)
        rest.remove(f)
        at = ('B', -f[0], tuple(sorted(rest)), Fr(-f[1]) / f[0], A, B)
        terms.append(F.count(f) * mp.exp(-X._mpq(mp, f[1])) * X.atom_value(mp, at, cache))
    bnd = _bnd(R, K, F, A, B, (Fr(0), ()))
    for k, w in X.canon(bnd, K).items():
        terms.append(X._mpq(mp, Fr(int(w.p), int(w.q))) * key_value(k, cache))
    s = sum(terms)
    return float(abs(s) / max(abs(t) for t in terms))


def stage_A1():
    """Identity check: every valid natural pair set (and the triple sets on one finite panel) on every panel of the
    two control geometries, at 60 digits."""
    res = {}
    for case in CONTROL_CASES:
        Fn, stable = natural_factors(case)
        rows = []
        t = time.time()
        for (A, B), fs in Fn.items():
            sets = _natural_sets(fs, A, B, sizes=(2,))
            if B is not None and A != 0:
                sets += [s for s in _natural_sets(fs, A, B, sizes=(3,))][:3]
            for F in sets:
                rows.append(dict(panel=[str(A), str(B)], F=[[str(c), str(d)] for c, d in F],
                                 rel=_identity_residual(F, A, B, 60)))
        res[case] = dict(factors_stable_L2_6=stable, n_sets=len(rows), max_rel=max(r['rel'] for r in rows),
                         n_triples=sum(len(r['F']) == 3 for r in rows),
                         n_fp_panel=sum(r['panel'][0] == '0' for r in rows),
                         n_inf_panel=sum(r['panel'][1] == 'None' for r in rows),
                         pass_1e_50=max(r['rel'] for r in rows) < 1e-50, seconds=round(time.time() - t, 1))
    save('A1_identity_numeric', res)


def _planted(case, K, seed=11):
    """Planted coefficient dicts on case's panels.  pos: random rational combination of whole identity sets (each
    member with its own e^{-d} prefactor and multiplicity) + depth-0 junk -> must reduce to depth 0.
    neg: the same plus one canonical (non-eliminated) member alone -> a depth-2 key must survive."""
    import random
    rnd = random.Random(seed)
    Fn, _ = natural_factors(case)
    pos, neg_extra = {}, None
    R = X.Ring(K)
    for (A, B), fs in Fn.items():
        for F in _natural_sets(fs, A, B)[:4]:
            w = K(Fr(rnd.randrange(1, 50), rnd.randrange(1, 9)))
            for f in set(F):
                rest = list(F)
                rest.remove(f)
                at = ('B', -f[0], tuple(sorted(rest)), Fr(-f[1]) / f[0], A, B)
                X.add(pos, (-f[1], (at,)), w * K(F.count(f)))
            for k, v in _bnd(R, K, F, A, B, (Fr(0), ())).items():
                X.add(pos, k, w * v)
            if neg_extra is None and len(set(F)) == 2:
                s = identity_set(-F[0][0], (F[1],), Fr(-F[0][1]) / F[0][0], A, B, fs)
                mem = s[1]
                keep = min(mem, key=lambda f: repr(mem[f][0]))
                key = mem[keep][0]
                neg_extra = ((-keep[1], (('B', key[0], key[1], key[2], A, B),)), K(Fr(3, 7)))
    junk = {(Fr(-2), (('E1', Fr(3)),)): K(Fr(5, 3)), (Fr(1, 2), ()): K(-4)}
    for k, v in junk.items():
        X.add(pos, k, v)
    neg = dict(pos)
    X.add(neg, neg_extra[0], neg_extra[1])
    return X.canon(pos, K), X.canon(neg, K), Fn


def stage_A2():
    """Planted controls for route (a): exact (fmpq) and mod p, value preserved numerically at 80 digits."""
    res = {}
    for case in CONTROL_CASES:
        K = NE.Field(None)
        pos, neg, Fn = _planted(case, K)
        st = {}
        cp, cn = canonical(pos, Fn, K, st), canonical(neg, Fn, K)
        mp.mp.dps = 80
        cache = {}

        def val(d):
            return sum(X._mpq(mp, Fr(int(v.p), int(v.q))) * key_value(k, cache) for k, v in d.items())
        vp, vcp, vn, vcn = val(pos), val(cp), val(neg), val(cn)
        Kp = NE.Field(P_AUX)
        posp, negp, _ = _planted(case, Kp)
        modp_ok = (canonical(posp, Fn, Kp) == X.reduce_mod(cp, P_AUX)
                   and canonical(negp, Fn, Kp) == X.reduce_mod(cn, P_AUX))
        res[case] = dict(pos_n_depth2_before=len(depth2_keys(pos)), pos_n_depth2_after=len(depth2_keys(cp)),
                         pos_reduced=len(depth2_keys(cp)) == 0,
                         neg_n_depth2_after=len(depth2_keys(cn)), neg_survives=len(depth2_keys(cn)) >= 1,
                         pos_value_rel=float(abs(vcp / vp - 1)), neg_value_rel=float(abs(vcn / vn - 1)),
                         modp_equals_fmpq_mod_p=modp_ok, eliminated=st.get('eliminated', 0))
    save('A2_planted_route_a', res)


def stage_A3(Lmax=5, dps=350):
    """Real engine output: canonical(F_L) has the same value as F_L (L = 2..Lmax, dps digits).  Control geometries:
    counts recorded.  M1, M2r: BOOLEANS ONLY (value preserved; rewrite exercised)."""
    res = {}
    for case in ('C1', 'C2', 'M1', 'M2r'):
        t = time.time()
        Fn, stable = natural_factors(case)
        K = NE.Field(None)
        Qs = FX.coefficients(case, Lmax, K)
        mp.mp.dps = dps
        cache = {}
        worst, nel = 0.0, 0
        for L in range(2, Lmax + 1):
            st = {}
            c = canonical(Qs[L], Fn, K, st)
            nel += st.get('eliminated', 0)
            a = sum(X._mpq(mp, Fr(int(v.p), int(v.q))) * key_value(k, cache) for k, v in Qs[L].items())
            b = sum(X._mpq(mp, Fr(int(v.p), int(v.q))) * key_value(k, cache) for k, v in c.items())
            worst = max(worst, float(abs(b / a - 1)))
        row = dict(value_preserved=worst < 1e-200, rewrite_exercised=nel > 0, factors_stable=stable,
                   seconds=round(time.time() - t, 1))
        if case in CONTROL_CASES:
            row.update(max_rel=worst, n_eliminated=nel)
        res[case] = row
    save('A3_engine_value_preserved', res)


def stage_A4(Lmax=10):
    """Mod-p canonical == fmpq canonical reduced mod p, control geometries, L <= Lmax (two primes)."""
    res = {}
    for case in CONTROL_CASES:
        Fn, _ = natural_factors(case)
        Kq = NE.Field(None)
        Qs = FX.coefficients(case, Lmax, Kq)
        cq = [canonical(q, Fn, Kq) for q in Qs]
        row = {}
        for p in (1073741441, 1073741381):
            Kp = NE.Field(p)
            Qp = FX.coefficients(case, Lmax, Kp)
            row[str(p)] = all(canonical(Qp[L], Fn, Kp) == X.reduce_mod(cq[L], p) for L in range(Lmax + 1))
        res[case] = row
    save('A4_modp_canonical', res)


def stage_A5():
    """a = b (K4A1, record AB's class): the raw output has no named base, so canonical() is the identity."""
    Fn, _ = natural_factors('K4A1')
    K = NE.Field(None)
    Qs = FX.coefficients('K4A1', 12, K)
    save('A5_a_eq_b_noop', dict(no_named_base=all(not depth2_keys(q) for q in Qs),
                               canonical_is_identity=all(canonical(q, Fn, K) == q for q in Qs)))


def stage_B(dps=300, H=10 ** 6, Lb=4):
    """Route (b) on the BASIS, control geometry C1 only: per panel, the canonical bases occurring in C1's F_L
    (L <= Lb) with their e^{e0} prefactors, plus the panel's endpoint products of every applied identity set and a
    few depth-0 constants.
    B1 (positive): add back an eliminated member -> LLL must return that identity.
    B2 (planted): one number = small-integer combination of a whole identity set + depth-0 junk, with the
       canonical set -> LLL must find a relation involving the planted number.
    B3 (negative): the canonical set alone -> no relation with |coef| <= H (at dps and at dps + 100)."""
    case = 'C1'
    res = dict(case=case, dps=dps, H=H)
    t0 = time.time()
    Fn, _ = natural_factors(case)
    K = NE.Field(None)
    Qs = FX.coefficients(case, Lb, K)
    can = [canonical(q, Fn, K) for q in Qs[2:]]
    keys = set(k for c in can for k in depth2_keys(c))
    by_panel = {}
    for k in keys:
        b = [a for a in k[1] if a[0] == 'B'][0]
        by_panel.setdefault((b[4], b[5]), set()).add(k)
    # the panel with the most canonical bases whose identity sets are applied
    def applied(k):
        b = [a for a in k[1] if a[0] == 'B'][0]
        return identity_set(b[1], b[2], b[3], b[4], b[5], Fn[(b[4], b[5])]) is not None
    panel = max(by_panel, key=lambda pn: (sum(applied(k) for k in by_panel[pn]), len(by_panel[pn])))
    ks = sorted(by_panel[panel], key=repr)
    A, B = panel
    res['panel'] = [str(A), str(B)]
    res['n_canonical_bases'] = len(ks)
    # pick one applied canonical base and its eliminated partner(s)
    pick = next(k for k in ks if applied(k))
    b = [a for a in pick[1] if a[0] == 'B'][0]
    F, mem = identity_set(b[1], b[2], b[3], A, B, Fn[panel])
    R = X.Ring(K)
    bnd = X.canon(_bnd(R, K, F, A, B, (Fr(0), ())), K)
    tgt = max(mem, key=lambda f: repr(mem[f][0]))
    elim_key = (-tgt[1], (('B', mem[tgt][0][0], mem[tgt][0][1], mem[tgt][0][2], A, B),))
    member_keys = [(-f[1], (('B', key[0], key[1], key[2], A, B),)) for f, (key, m) in mem.items()]
    endpoint = sorted(bnd, key=repr)
    junk = [(Fr(0), (('E1', Fr(3)),)), (Fr(-1), ()), (Fr(0), (('g',),))]
    base_list = sorted(set(ks) | set(member_keys) - {elim_key}, key=repr)

    def vals(klist, prec, extra=None):
        mp.mp.dps = prec
        cache = {}
        v = [key_value(k, cache) for k in klist]
        if extra is not None:
            v = [extra(cache)] + v
        return v
    # B1
    kl1 = [elim_key] + base_list + endpoint + junk
    t = time.time()
    rel = int_relation(vals(kl1, dps), dps, H, check=lambda: vals(kl1, dps + 60))
    expect = {elim_key: mem[tgt][1]}
    for f, (key, m) in mem.items():
        if f != tgt:
            expect[(-f[1], (('B', key[0], key[1], key[2], A, B),))] = m
    for k, w in bnd.items():
        expect[k] = expect.get(k, 0) + int(w)
    found = None if rel is None else {kl1[i]: c for i, c in enumerate(rel) if c}
    ok1 = False
    if found is not None:
        s = found[elim_key] // expect[elim_key] if elim_key in found else 0
        ok1 = s != 0 and all(found.get(k, 0) == s * v for k, v in expect.items()) and \
            all(k in expect for k in found)
    res['B1_positive'] = dict(n=len(kl1), relation_found=rel is not None, equals_identity=ok1,
                              seconds=round(time.time() - t, 1))
    # B2 planted
    coefs = {k: 3 * v for k, v in expect.items()}

    def planted(cache):
        return sum(c * key_value(k, cache) for k, c in coefs.items()) + 5 * key_value(junk[0], cache)
    kl2 = base_list + endpoint + junk
    t = time.time()
    rel2 = int_relation(vals(kl2, dps, planted), dps, H, check=lambda: vals(kl2, dps + 60, planted))
    res['B2_planted'] = dict(n=len(kl2) + 1, relation_found=rel2 is not None,
                             involves_planted=rel2 is not None and rel2[0] != 0, seconds=round(time.time() - t, 1))
    # B3 negative
    kl3 = base_list + endpoint + junk
    t = time.time()
    rel3 = int_relation(vals(kl3, dps), dps, H, check=lambda: vals(kl3, dps + 60))
    rel3b = int_relation(vals(kl3, dps + 100), dps + 100, H)
    res['B3_negative'] = dict(n=len(kl3), relation_found=rel3 is not None, relation_found_dps_plus_100=rel3b is not None,
                              seconds=round(time.time() - t, 1))
    res['seconds'] = round(time.time() - t0, 1)
    save('B_route_b_lll', res)


def stage_B4(dps=450, H=10 ** 4, Lb=6):
    """Route (b), whole constant set: every canonical constant of C1's F_L (L = 2..Lb; depth 0 and depth 2 together,
    e^{e0} included).  B4neg: no integer relation with |coef| <= H at dps and dps + 100.  B4pos: the same set plus one
    eliminated member of an applied identity set -> a relation is found and it involves that member."""
    case = 'C1'
    t0 = time.time()
    Fn, _ = natural_factors(case)
    K = NE.Field(None)
    Qs = FX.coefficients(case, Lb, K)
    can = [canonical(q, Fn, K) for q in Qs[2:]]
    ks = sorted(set(k for c in can for k in c), key=repr)
    res = dict(case=case, Lb=Lb, dps=dps, H=H, n=len(ks), n_depth2=sum(1 for k in ks if any(a[0] == 'B' for a in k[1])))
    elim, added = None, []
    for k in ks:                         # an occurring canonical key e^{e0} B_f whose identity set is applied
        e0, atoms = k
        Bs = [a for a in atoms if a[0] == 'B']
        if not Bs or len(atoms) != 1:
            continue
        a = Bs[0]
        s = identity_set(a[1], a[2], a[3], a[4], a[5], Fn[(a[4], a[5])])
        if s is None:
            continue
        F, mem = s
        tgt = max(mem, key=lambda f: repr(mem[f][0]))
        me = next(f for f in mem if mem[f][0] == (a[1], a[2], a[3]))
        sh = e0 + me[1]                  # e^{e0} B_f = e^{sh} (e^{-d_f} B_f): the identity scaled by e^{sh}
        elim = (sh - tgt[1], (('B', mem[tgt][0][0], mem[tgt][0][1], mem[tgt][0][2], a[4], a[5]),))
        R = X.Ring(K)
        for kb in X.canon(_bnd(R, K, F, a[4], a[5], (sh, ())), K):
            if kb not in ks:
                added.append(kb)
        for f, (key, m) in mem.items():
            kf = (sh - f[1], (('B', key[0], key[1], key[2], a[4], a[5]),))
            if f != tgt and kf not in ks:
                added.append(kf)
        break
    res['B4pos_added_keys'] = len(added)

    def vals(kl, prec):
        mp.mp.dps = prec
        cache = {}
        return [key_value(k, cache) for k in kl]
    t = time.time()
    r = int_relation(vals(ks, dps), dps, H, check=lambda: vals(ks, dps + 60))
    r2 = int_relation(vals(ks, dps + 100), dps + 100, H)
    res['B4neg'] = dict(relation_found=r is not None, relation_found_dps_plus_100=r2 is not None,
                        seconds=round(time.time() - t, 1))
    kl = [elim] + ks + added
    t = time.time()
    rp = int_relation(vals(kl, dps), dps, H, check=lambda: vals(kl, dps + 60))
    res['B4pos'] = dict(relation_found=rp is not None, involves_eliminated=rp is not None and rp[0] != 0,
                        seconds=round(time.time() - t, 1))
    res['seconds'] = round(time.time() - t0, 1)
    save('B4_route_b_whole_set', res)


STAGES = dict(A1=stage_A1, A2=stage_A2, A3=stage_A3, A4=stage_A4, A5=stage_A5, B=stage_B, B4=stage_B4)

if __name__ == '__main__':
    if sys.argv[1] == 'controls':
        for s in STAGES:
            t = time.time()
            STAGES[s]()
            print('stage', s, round(time.time() - t, 1), 's', flush=True)
    elif sys.argv[1] == 'stage':
        for s in sys.argv[2:]:
            t = time.time()
            STAGES[s]()
            print('stage', s, round(time.time() - t, 1), 's', flush=True)
