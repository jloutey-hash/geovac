"""Record S1 (to be registered; inventory row 6, singular-points column): controls and the shared pipeline.

S1 asks whether every genuine singular point of the bond-angle operator of (AA|BC) that moves with mu = 2 zeta_A is a
genuine singular point of the operator of one of record AB's 21 constant pieces (F_L = sum_c q_{L,c} c).
This file holds the pipeline the scored driver (debug/s1_run.py) imports, and its positive controls. It computes NO
scored quantity: no piece operator of AB's cases A1/A2 is built here, and nothing is compared with the full operator's
mu-moving points.

Pipeline (all mod one prime p; functions of na_modguess / na_modanalyze / intclass_lib / aabc_exact called unchanged):
  piece_operator(seq)       minimal x-operator of one Legendre-coefficient sequence (intclass_lib.minimal_operator,
                            the gcrd of the independent annihilators, as IC/AK)
  generic_operator(seq)     minimal operator of sum_c w_c seq_c, by the same gcrd routine.  NOT AB's operator:
                            AB's search (generic_operator_ab, first null vector at the smallest findable order) returns
                            at A1 an order-17 left multiple of the order-14 minimal operator (control W0)
  points(c)                 every F_p-irreducible factor of the leading coefficient with its class:
                              M  some exponent not an integer (rational non-integer, not reconstructable, or not
                                 constant over the factor's roots): nontrivial local monodromy
                              L  integer exponents, a logarithm forced (linear factors): unipotent monodromy
                              P  integer exponents, one negative, log-free: a pole, trivial monodromy
                              A  apparent-type, or integer exponents >= 0 with no log (linear factors)
                              I  integer exponents, degree > 1, log test not run
                              X  irregular (non-Fuchsian)
                            genuine = M, X, L, P
                            exponents mod 1 (linear factors; rational ones) for the monodromy comparison.
  weight_fixed(...)         factors of lc common to the generic operators of two independent weight vectors: a factor
                            that moves with the weights is an artifact of the generic combination, not of f.

Modes:
  python debug/s1_controls.py toy <p>        toy LCLM controls (known points; collapse; fire test)
  python debug/s1_controls.py na <p>         NA at a fresh geometry N3 (not AB's): per-piece vs generic points
  python debug/s1_controls.py w0 <p>         A1 full operator under two weight vectors (no mu shift, no pieces),
                                             AB's search and the minimal operator; caches A1's coefficients
  python debug/s1_controls.py summary        -> debug/data/s1_controls.json
"""
from __future__ import annotations

import json
import random
import sys
import time
from math import isqrt
from fractions import Fraction
from pathlib import Path

import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import intclass_lib as I        # noqa: E402
import na_modanalyze as A       # noqa: E402
import na_modguess as G         # noqa: E402

DATA = HERE / 'data'
PRIMES = (1073741441, 1073741381)          # record AK's primes (< 2^31, as intclass_lib requires)
SEEDS = (2026, 7)                          # 2026 = AB's p3a seed


# ------------------------------------------------------------------------------------------------- pipeline
def nm(v, p):
    return flint.nmod(int(v) % p, p)


def piece_operator(seq, p, Rmax=24):
    r = I.minimal_operator([int(v) for v in seq], p, Rmax=Rmax)
    if r is None:
        return None
    return dict(c=r['c'], order=r['order'], degree=r['degree'], found=r['found'])


def generic_seq(seqs, p, seed):
    keys = sorted(seqs, key=str)
    rnd = random.Random(seed)
    w = {k: rnd.randrange(1, p) for k in keys}
    n = min(len(seqs[k]) for k in keys)
    return [sum(w[k] * int(seqs[k][L]) for k in keys) % p for L in range(n)]


def generic_operator_ab(seq, p, Rmax=24, margin=30):
    """AB's search as aabc_p3a_points.py: first null vector at the smallest findable order (NOT minimal: at A1 this
    is an order-17 left multiple of the order-14 minimal operator; kept for the W0 comparison only)."""
    B = G.Basis([nm(v, p) for v in seq], p)
    for R in range(1, Rmax + 1):
        Q0, nul = G.min_degree_bisect(B, R, margin, len(seq))
        if Q0 is not None:
            nv, ks, (N, nrow) = G.nullity(B, R, Q0, margin)
            c = A.standard_form(ks, [N[i, 0] for i in range(N.nrows())], p)
            return dict(c=I.content_reduce(c), order=R, degree=Q0, nullity=nv)
    return None


def generic_operator(seq, p, Rmax=24):
    """Minimal operator of the generic sequence: gcrd of its independent annihilators (intclass_lib.minimal_operator,
    unchanged), the same routine as the pieces."""
    r = I.minimal_operator([int(v) for v in seq], p, Rmax=Rmax)
    if r is None:
        return None
    return dict(c=r['c'], order=r['order'], degree=r['degree'], found=r['found'], nullity=None)


def fkey(f):
    f = f * (1 / f.leading_coefficient())
    return tuple(int(a) for a in f.coeffs())


def points(c, p):
    """Classes: M non-integer exponent; X irregular; L integer exponents with a forced logarithm (linear factors,
    intclass_lib.log_free_count); P integer exponents, a negative one, no logarithm (pole); A apparent-type;
    I integer exponents otherwise (incl. degree > 1 factors, where the log test is not run)."""
    R = len(c) - 1
    out = {}
    for d in I.finite_points(c, p):
        f = d['factor']
        kinds = [k for (k, v), m in d['croots']]
        nonint = [v for (k, v), m in d['croots'] if k == 'rat' and Fraction(v).denominator != 1]
        pos, x0, logfree = None, None, None
        if f.degree() == 1:
            x0 = (-int(f[0]) * pow(int(f[1]), p - 2, p)) % p
            fr = I.ratrec(x0, p, bound=isqrt(p // 2))
            pos = str(fr) if fr is not None else f'mod-p:{x0}'
        if not d['regular']:
            cls = 'X'
        elif d['unaccounted'] > 0 or 'other' in kinds or nonint:
            cls = 'M'
        elif d['apparent_type']:
            cls = 'A'
        elif f.degree() == 1 and len(d['ints']) == R:
            logfree = I.log_free_count(c, x0, p, d['ints'])
            cls = 'L' if logfree < R else ('P' if min(d['ints']) < 0 else 'A')
        else:
            cls = 'I'
        modz = sorted({str(Fraction(v) % 1) for (k, v), m in d['croots'] if k in ('rat', 'int')} |
                      {f'irr:{v}' for (k, v), m in d['croots'] if k == 'other'})
        exps = [str(v) for (k, v), m in d['croots'] for _ in range(m)]
        out[fkey(f)] = dict(deg=f.degree(), cls=cls, regular=d['regular'], modz=modz, log_free=logfree, exps=exps,
                            unaccounted=d['unaccounted'], ints=d['ints'], pos=pos)
    return out


def describe(pts):
    return sorted(((v['deg'], v['cls'], v['pos'], v['modz']) for v in pts.values()), key=str)


def weight_fixed(pts1, pts2):
    return {k: v for k, v in pts1.items() if k in pts2}


GENUINE = ('M', 'X', 'L', 'P')       # non-analytic solutions at the point
MONODROMY = ('M', 'X', 'L')          # nontrivial local monodromy (X: assumed)


def genuine(pts):
    return {k: v for k, v in pts.items() if v['cls'] in GENUINE}


def attribute(full_pts, piece_pts):
    """For each genuine factor of full_pts: the pieces carrying it (a monodromy point M/X/L needs a carrier with
    monodromy there; a pole P needs any genuine carrier), and whether the full operator's exponents mod 1 (linear
    factors) are covered by the union of the carriers' exponents mod 1."""
    out = {}
    for k, v in genuine(full_pts).items():
        ok = MONODROMY if v['cls'] in MONODROMY else GENUINE
        carriers = [name for name, pp in piece_pts.items() if k in pp and pp[k]['cls'] in ok]
        roots_any = [name for name, pp in piece_pts.items() if k in pp]
        cov = None
        if v['deg'] == 1 and carriers:
            fm = {e for e in v['modz'] if e != '0'}
            pm = set().union(*({e for e in piece_pts[n][k]['modz']} for n in carriers))
            cov = fm <= pm
        out[k] = dict(pos=v['pos'], deg=v['deg'], cls=v['cls'], modz=v['modz'], carriers=carriers, lc_roots_in=roots_any,
                      modz_covered=cov)
    return out


def mo2_points(wf, wf_shift, exclude):
    """Genuine weight-fixed points of BOTH sides, minus exclude (fkeys).  A base point moves if its factor is absent
    from the shifted weight-fixed set (it vanished); a shifted point moves if absent from the base set (it appeared,
    at its new position).  Non-moving base points are listed too (moves=False).  integer_only: classes L, P."""
    out = []
    for side, A, B in (('base', wf, wf_shift), ('shift', wf_shift, wf)):
        for k, v in sorted(genuine(A).items(), key=lambda kv: str(kv[1]['pos'])):
            if k in exclude:
                continue
            mv = k not in B
            if side == 'shift' and not mv:
                continue                       # already listed from the base side
            out.append(dict(side=side, pos=v['pos'], deg=v['deg'], cls=v['cls'], exps=v['exps'],
                            unaccounted=v['unaccounted'], moves=mv, integer_only=v['cls'] in ('L', 'P')))
    return out


def mo2_holds(rows):
    return all(r['integer_only'] for r in rows if r['moves'])


def union_genuine(piece_pts):
    u = {}
    for pp in piece_pts.values():
        u.update(genuine(pp))
    return u


# ------------------------------------------------------------------------------------------------- toy controls
def leg_seqs(p, n):
    """Exact Legendre-coefficient sequences with known functions (generating function of P_L):
         a: t^L          <-> (1 - 2xt + t^2)^(-1/2)          point x_t = (1+t^2)/(2t), exponent -1/2
         b: L t^L        <-> t d/dt of a = (xt - t^2)(...)^(-3/2)   x_t exponent -3/2; zero at x = t (apparent)
         c: (2L+1) s^L   <-> (1 - s^2)(1 - 2xs + s^2)^(-3/2) point x_s, exponent -3/2
         d: u^L          <-> point x_u (planted: used only in the fire test, never given as a piece)
         c2: (2L+1) s'^L  the 'mu-shifted' c (s -> s'), to emulate a moving point"""
    t, s, u, s2 = Fraction(1, 3), Fraction(2, 5), Fraction(3, 7), Fraction(3, 8)
    T, S, U, S2 = (I.to_fp(v, p) for v in (t, s, u, s2))
    seqs = dict(a=[pow(T, L, p) for L in range(n)], b=[L * pow(T, L, p) % p for L in range(n)],
                c=[(2 * L + 1) * pow(S, L, p) % p for L in range(n)], d=[pow(U, L, p) for L in range(n)],
                c2=[(2 * L + 1) * pow(S2, L, p) % p for L in range(n)])
    xs = {k: (1 + v * v) / (2 * v) for k, v in (('x_t', t), ('x_s', s), ('x_u', u), ('x_s2', s2))}
    return seqs, xs


def toy(p, n=400):
    t0 = time.time()
    seqs, xs = leg_seqs(p, n)
    res = dict(prime=p, known_points={k: str(v) for k, v in xs.items()})
    pieces = {k: piece_operator(seqs[k], p) for k in ('a', 'b', 'c')}
    ppts = {k: points(v['c'], p) for k, v in pieces.items()}
    res['pieces'] = {k: dict(order=pieces[k]['order'], degree=pieces[k]['degree'], points=describe(ppts[k]))
                     for k in ppts}
    # T1: each piece's genuine points are its known point, with the known exponent mod 1 (1/2)
    exp = dict(a=[str(xs['x_t'])], b=[str(xs['x_t'])], c=[str(xs['x_s'])])
    res['T1_piece_points'] = all(sorted(v['pos'] for v in genuine(ppts[k]).values()) == exp[k] and
                                 all(v['modz'] == ['1/2'] for v in genuine(ppts[k]).values()) for k in exp)
    res['T1b_b_apparent_at_t'] = any(v['pos'] == '1/3' and v['cls'] == 'A' for v in ppts['b'].values())
    # T2: generic of a, b, c under two weights: order 2 (a, b collapse: b = rational multiple of a);
    #     weight-fixed genuine points = {x_t, x_s}; the weight-moving factors are apparent
    sub = {k: seqs[k] for k in ('a', 'b', 'c')}
    gens = [generic_operator(generic_seq(sub, p, sd), p) for sd in SEEDS]
    gpts = [points(g['c'], p) for g in gens]
    wf = weight_fixed(*gpts)
    moving = {k: v for k, v in gpts[0].items() if k not in gpts[1]}
    res['generic'] = [dict(order=g['order'], degree=g['degree'], nullity=g['nullity'], points=describe(pt))
                      for g, pt in zip(gens, gpts)]
    res['T2_order'] = [g['order'] for g in gens]
    res['T2_weight_fixed_genuine'] = sorted(v['pos'] for v in genuine(wf).values())
    res['T2_weight_moving'] = describe(moving)
    res['T2'] = (res['T2_order'] == [2, 2] and res['T2_weight_fixed_genuine'] == sorted([str(xs['x_t']),
                                                                                        str(xs['x_s'])])
                 and all(v['cls'] == 'A' for v in moving.values()) and len(moving) >= 1)
    # T3: attribution and union equality (the S1a / S1c logic) on the toy
    att = attribute(wf, ppts)
    res['T3_attribution'] = {v['pos']: v['carriers'] for v in att.values()}
    res['T3'] = (res['T3_attribution'] == {str(xs['x_t']): ['a', 'b'], str(xs['x_s']): ['c']}
                 and all(v['modz_covered'] for v in att.values())
                 and set(union_genuine(ppts)) == set(genuine(wf)))
    # T4: 'moving' emulation: c -> c2 moves x_s; the moving-and-weight-fixed genuine set is {x_s}
    sub2 = {'a': seqs['a'], 'b': seqs['b'], 'c': seqs['c2']}
    gsh = generic_operator(generic_seq(sub2, p, SEEDS[0]), p)
    psh = points(gsh['c'], p)
    mov = {k: v for k, v in genuine(wf).items() if k not in psh}
    res['T4_moving_genuine'] = sorted(v['pos'] for v in mov.values())
    res['T4'] = res['T4_moving_genuine'] == [str(xs['x_s'])]
    # T5 fire test: the 'full' operator includes d, the pieces do not: S1a must report an uncarried point at x_u
    subf = {'a': seqs['a'], 'b': seqs['b'], 'c': seqs['c'], 'd': seqs['d']}
    gf = [points(generic_operator(generic_seq(subf, p, sd), p)['c'], p) for sd in SEEDS]
    attf = attribute(weight_fixed(*gf), ppts)
    unc = sorted(v['pos'] for v in attf.values() if not v['carriers'])
    res['T5_uncarried'] = unc
    res['T5'] = unc == [str(xs['x_u'])]
    # T6 log class: e_L = t^(L+1)/(L+1) <-> int_0^t (1 - 2xs + s^2)^(-1/2) ds = log((t - x + sqrt(.))/(1 - x)):
    #    x = 1 is a logarithmic point (integer exponents, class L), x_t keeps exponent 1/2 mod 1 (class M), and
    #    x = -1 is logarithmic on the other branch of the square root, log((t - x - sqrt(.))/(1 - x)) (class L)
    T = I.to_fp(Fraction(1, 3), p)
    e = [pow(T, L + 1, p) * pow(L + 1, p - 2, p) % p for L in range(n)]
    pe = points(piece_operator(e, p)['c'], p)
    res['T6_points'] = describe(genuine(pe))
    res['T6'] = sorted((v['pos'], v['cls']) for v in genuine(pe).values()) == [('-1', 'L'), ('1', 'L'), (str(xs['x_t']), 'M')]
    res['time'] = round(time.time() - t0, 1)
    # T7 (MO2 positive control): a mu-moving point WITH a non-integer exponent is reported as such.  The full
    #    combination of a = t^L and c = (2L+1) s^L under two weights; the 'shift' moves t = 1/3 -> 3/10, so x_t moves
    #    and carries exponent 1/2 mod 1: mo2_points must list it with moves=True and integer_only=False.
    T2_ = I.to_fp(Fraction(3, 10), p)
    a2 = [pow(T2_, L, p) for L in range(n)]
    base = {'a': seqs['a'], 'c': seqs['c']}
    shf = {'a': a2, 'c': seqs['c']}
    wf0 = weight_fixed(*[points(generic_operator(generic_seq(base, p, sd), p)['c'], p) for sd in SEEDS])
    wf1 = weight_fixed(*[points(generic_operator(generic_seq(shf, p, sd), p)['c'], p) for sd in SEEDS])
    m2 = mo2_points(wf0, wf1, exclude=())
    res['T7_mo2'] = m2
    xt2 = (1 + Fraction(3, 10) ** 2) / (2 * Fraction(3, 10))
    res['T7'] = (sorted((r['side'], r['pos'], r['integer_only']) for r in m2 if r['moves'])
                 == [('base', str(xs['x_t']), False), ('shift', str(xt2), False)] and not mo2_holds(m2))
    res['all_pass'] = all(res[k] for k in ('T1_piece_points', 'T1b_b_apparent_at_t', 'T2', 'T3', 'T4', 'T5', 'T6',
                                           'T7'))
    print(json.dumps(res, indent=1, default=str))
    (DATA / f's1_ctrl_toy_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    return res


# ------------------------------------------------------------------------------------------------- NA control
N3 = ('7/5', '11/5', '6/5', '1')            # D1, D2, zeta_B (X), zeta_C (H): not N1, N2, A1, A2


def na_control(p, Lmax=1000):
    import na_exact as E
    t0 = time.time()
    K = E.Field(p)
    q = E.coefficients(*N3, Lmax, K)
    keys = sorted(set(k for d in q for k in d), key=str)
    seqs = {k: [int(d.get(k, flint.nmod(0, p))) for d in q] for k in keys}
    tc = round(time.time() - t0, 1)
    D1, D2 = Fraction(N3[0]), Fraction(N3[1])
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    pieces, ppts = {}, {}
    for k in keys:
        t1 = time.time()
        po = piece_operator(seqs[k], p)
        pieces[repr(k)] = dict(order=po['order'], degree=po['degree'], time=round(time.time() - t1, 1))
        ppts[repr(k)] = points(po['c'], p)
        pieces[repr(k)]['points'] = describe(ppts[repr(k)])
    gens = [generic_operator(generic_seq(seqs, p, sd), p) for sd in SEEDS]
    gpts = [points(g['c'], p) for g in gens]
    wf = weight_fixed(*gpts)
    att = attribute(wf, ppts)
    U, F = union_genuine(ppts), genuine(wf)
    res = dict(prime=p, geometry=dict(zip(('D1', 'D2', 'zB', 'zC'), N3)), xstar=str(xs), coeff_time=tc,
               n_pieces=len(keys), pieces=pieces,
               generic=[dict(order=g['order'], degree=g['degree'], found=g['found'], points=describe(pt))
                        for g, pt in zip(gens, gpts)],
               weight_fixed=describe(wf), n_weight_moving=len(set(gpts[0]) - set(gpts[1])),
               attribution={str(v['pos'] or k[:3]): dict(carriers=v['carriers'], modz=v['modz'],
                                                        covered=v['modz_covered']) for k, v in att.items()},
               union_pieces_genuine=describe(U), full_genuine=describe(F),
               union_equals_full=set(U) == set(F))
    # expectations from record NA (measured at N1, N2): order-8 generic operator; x* the only finite point with a
    # non-integer exponent (3/2, 5/2 -> 1/2 mod 1); every genuine point carried by a piece
    res['C1_generic_order8'] = [g['order'] for g in gens] == [8, 8]
    res['C2_full_noninteger_is_xstar'] = sorted(v['pos'] for v in F.values() if v['cls'] == 'M') == [str(xs)]
    res['C2b_full_genuine_classes'] = sorted((v['pos'], v['cls']) for v in F.values())
    res['C3_carried'] = all(v['carriers'] and v['modz_covered'] for v in att.values())
    res['C4_union_equals_full'] = res['union_equals_full']
    res['time'] = round(time.time() - t0, 1)
    print(json.dumps(res, indent=1, default=str))
    (DATA / f's1_ctrl_na_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    return res


# ------------------------------------------------------------------------------------------------- W0 control
def a_coeffs(case, p, Lmax=1000, zA_shift=None):
    """AB's 21 sequences of f (aabc_exact.pieces + combine, unchanged), cached to debug/data."""
    import aabc_exact as X
    import na_exact as E
    tag = '' if zA_shift is None else '_shift'
    fn = DATA / f's1_coeffs_{case}_{p}{tag}.json'
    if fn.exists():
        d = json.loads(fn.read_text())
        return {eval(k): v for k, v in d['seqs'].items()}, d['time']   # noqa: S307 (own file, Fraction keys)
    D1, D2, zB, zC, zA = X.params(case)
    if zA_shift is not None:
        zA = zA + zA_shift
    t0 = time.time()
    K = E.Field(p)
    pc = X.pieces(case, Lmax, K, params_override=(D1, D2, zB, zC, zA))
    f = X.combine(pc, zA, K)
    keys = sorted(set(k for q in f for k in q), key=str)
    seqs = {k: [int(q.get(k, flint.nmod(0, p))) for q in f] for k in keys}
    tt = round(time.time() - t0, 1)
    fn.write_text(json.dumps(dict(case=case, prime=p, zA=str(zA), time=tt, seqs={repr(k): v for k, v in seqs.items()})))
    return seqs, tt


def w0(p, case='A1'):
    """Full operator only, two weight vectors, NO mu shift and NO piece operator.  (i) AB's search (first null vector,
    order 17): which lc factors are weight-fixed; (ii) the minimal operator (gcrd): its order and its non-apparent
    points, and whether any non-apparent factor moves with the weights."""
    t0 = time.time()
    seqs, tc = a_coeffs(case, p)
    res = dict(case=case, prime=p, coeff_time=tc, n_const=len(seqs))
    for tag, fn in (('ab_search', generic_operator_ab), ('minimal', generic_operator)):
        t1 = time.time()
        gens = [fn(generic_seq(seqs, p, sd), p) for sd in SEEDS]
        gpts = [points(g['c'], p) for g in gens]
        wf = weight_fixed(*gpts)
        mov = {k: v for k, v in gpts[0].items() if k not in gpts[1]}
        res[tag] = dict(ops=[dict(order=g['order'], degree=g['degree'], n_factors=len(pt),
                                  factor_degrees=sorted(v['deg'] for v in pt.values())) for g, pt in zip(gens, gpts)],
                        weight_fixed=describe(wf),
                        weight_fixed_classes={c: sum(1 for v in wf.values() if v['cls'] == c) for c in 'MXLPIA'},
                        weight_moving_classes={c: sum(1 for v in mov.values() if v['cls'] == c) for c in 'MXLPIA'},
                        weight_moving_degrees=sorted(v['deg'] for v in mov.values()),
                        weight_fixed_detail=[dict(pos=v['pos'], cls=v['cls'], ints=v['ints'], modz=v['modz'],
                                                  log_free=v['log_free']) for v in wf.values() if v['cls'] != 'A'],
                        seconds_two_ops=round(time.time() - t1, 1))
    res['time'] = round(time.time() - t0, 1)
    print(json.dumps(res, indent=1, default=str))
    (DATA / f's1_ctrl_w0_{case}_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    return res


def summary():
    out = {}
    for f in sorted(DATA.glob('s1_ctrl_*.json')):
        out[f.stem] = json.loads(f.read_text())
    (DATA / 's1_controls.json').write_text(json.dumps(out, indent=1, default=str))
    for k, v in out.items():
        print(k, {kk: vv for kk, vv in v.items() if kk.startswith(('T', 'C', 'all', 'time'))})


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'toy':
        toy(int(sys.argv[2]))
    elif m == 'na':
        na_control(int(sys.argv[2]))
    elif m == 'w0':
        w0(int(sys.argv[2]))
    elif m == 'summary':
        summary()
