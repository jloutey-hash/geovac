"""Record S1 scored driver (PROPOSED; written before registration, NOT run on A1/A2). Inventory row 6.

Object: f(x) = (AA|BC) in x = cos(angle C B A), record AB's cases A1, A2; pieces = AB's 21 sequences q_{L,c}
(L = 0..1000, frozen engine via s1_controls.a_coeffs); full operator = the MINIMAL operator (gcrd of annihilators,
intclass_lib.minimal_operator) of the generic sum_c w_c q_{L,c}; two weight vectors (seeds 2026, 7); the zeta_A
shift is +1/997 (AB's).  Primes: AK's 1073741441, 1073741381.  Pipeline and point classes: debug/s1_controls.py.

Clauses
  MO0 (instrument; VOID on failure): controls T1-T7 + N3 + div rerun on both primes; the shifted coefficient
      sequences differ from the base ones; minimal operator found at order <= 24 (both weights, base and shift);
      minimal (order, degree) equal across weights and primes, base and shift separately; x* a weight-fixed class-M
      point on both sides; AB's search operator right-divisible by the minimal one (pseudo-remainder 0); no genuine
      weight-fixed point (either side) irregular (flag irregular_point) or with unaccounted exponents.
  MO1 (A2, fresh): the only finite point of class M (non-integer exponent) is A2's x* = 89/80; exponents mod 1 there
      {0, 1/2}; its half-integer exponents form a subset of {3/2, 7/2, 9/2}; x*'s exponents equal under both weights.
  MO2 (both cases): a genuine weight-fixed point other than x*, +-1 MOVES if it is in the base set and not the
      shifted one (vanished) or in the shifted set and not the base one (appeared, new position).  Holds iff every
      moving point has integer exponents only (class L or P).  Per point: side, position, class, exponents, moves.
  MO3: S1a carriage: every genuine weight-fixed point of the full operator is carried by a piece of the right class
      (s1_controls.attribute); a miss is VOID (a theorem given correct operators).
  MO4 (reported): attribution table; per-piece carried degree equal across primes.

CLI
  python debug/s1_run.py controls <p>      MO0 controls (T1-T7, N3) on prime p   (~7 min)
  python debug/s1_run.py run <case> <p>    one (case, prime)  -> debug/data/s1_run_<case>_<p>.json
  python debug/s1_run.py summary           -> debug/data/s1_run.json (clauses, gate, timing)
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import intclass_lib as I     # noqa: E402
import s1_controls as C      # noqa: E402
from aabc_exact import params  # noqa: E402

DATA = HERE / 'data'
PRIMES = C.PRIMES
CASES = ('A1', 'A2')
SHIFT = Fraction(1, 997)
HALF_ALLOWED = {'3/2', '7/2', '9/2'}


def fkey_lin(x0, p):
    """monic linear factor x - x0 as a points() key"""
    return ((-I.to_fp(x0, p)) % p, 1)


def rem_zero(L1, L2, p):
    r = I.right_rem(L1, L2, p)
    return len(r) == 1 and r[0].degree() < 0


def full(seqs, p):
    gen = [C.generic_seq(seqs, p, sd) for sd in C.SEEDS]
    ops = [C.generic_operator(g, p) for g in gen]
    if any(o is None for o in ops):
        return None
    pts = [C.points(o['c'], p) for o in ops]
    return dict(gen=gen, ops=ops, pts=pts, wf=C.weight_fixed(*pts))


def run(case, p):
    T = {}
    T0 = time.time()
    seqs, T['coeff'] = C.a_coeffs(case, p)
    seqs_s, T['coeff_shift'] = C.a_coeffs(case, p, zA_shift=SHIFT)
    t = time.time()
    F, Fs = full(seqs, p), full(seqs_s, p)
    T['full_minimal'] = round(time.time() - t, 1)
    res_same = (set(seqs) == set(seqs_s)) and all(seqs[k] == seqs_s[k] for k in seqs)
    D1, D2 = params(case)[:2]
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    res = dict(case=case, prime=p, xstar=str(xs))
    void = []
    if res_same:
        void.append('shifted coefficient sequences equal the base ones')
    if F is None or Fs is None:
        res['void_flags'] = ['minimal operator not found at order <= 24']
        (DATA / f's1_run_{case}_{p}.json').write_text(json.dumps(res, indent=1, default=str))
        return
    res['full'] = [dict(tag=tg, order=o['order'], degree=o['degree'], found=o['found'])
                   for tg, FF in (('base', F), ('shift', Fs)) for o in FF['ops']]
    for tg, FF in (('base', F), ('shift', Fs)):
        if len({(o['order'], o['degree']) for o in FF['ops']}) != 1:
            void.append(f'minimal (order, degree) differs between weights ({tg})')
    # MO0: AB's search operator right-divisible by the minimal one (same generic sequence, seed 2026)
    t = time.time()
    ab = C.generic_operator_ab(F['gen'][0], p)
    res['ab_search'] = None if ab is None else dict(order=ab['order'], degree=ab['degree'], nullity=ab['nullity'])
    res['MO0_ab_right_divisible'] = ab is not None and rem_zero(ab['c'], F['ops'][0]['c'], p)
    T['ab_search'] = round(time.time() - t, 1)
    if not res['MO0_ab_right_divisible']:
        void.append('AB search operator not right-divisible by the minimal operator')
    wf, wfs = F['wf'], Fs['wf']
    kx = fkey_lin(xs, p)
    k1, km1 = fkey_lin(Fraction(1), p), fkey_lin(Fraction(-1), p)
    for tg, W in (('base', wf), ('shift', wfs)):
        if not (kx in W and W[kx]['cls'] == 'M'):
            void.append(f'x* not a weight-fixed class-M point ({tg})')
        for v in C.genuine(W).values():
            if v['cls'] == 'X':
                void.append(f'irregular_point ({tg}, {v["pos"]})')
            if v['unaccounted'] > 0:
                void.append(f'unaccounted exponents ({tg}, {v["pos"]}, deg {v["deg"]})')
    # MO1 (scored for A2; reported for A1)
    Ms = [v for v in wf.values() if v['cls'] == 'M']
    xsd = wf.get(kx, {})
    xs_w = [sorted(P[kx]['exps'], key=Fraction) if kx in P else None for P in F['pts']]
    halves = {e for e in xsd.get('exps', []) if Fraction(e).denominator == 2} if xsd else set()
    res['MO1'] = dict(class_M_points=[(v['pos'], v['exps']) for v in Ms], xstar_exps=xsd.get('exps'),
                      xstar_modz=xsd.get('modz'), xstar_exps_by_weight=xs_w,
                      holds=(len(Ms) == 1 and xs_w[0] is not None and xs_w[0] == xs_w[1] and kx in wf and wf[kx]['cls'] == 'M' and xsd['modz'] == ['0', '1/2']
                             and halves <= HALF_ALLOWED and bool(halves)))
    # MO2
    rows = C.mo2_points(wf, wfs, exclude=(kx, k1, km1))
    res['MO2'] = dict(points=rows, holds=C.mo2_holds(rows), n_moving=sum(r['moves'] for r in rows),
                      n_moving_base_vanished=sum(r['moves'] and r['side'] == 'base' for r in rows),
                      n_moving_shift_appeared=sum(r['moves'] and r['side'] == 'shift' for r in rows))
    # pieces, MO3, MO4
    t = time.time()
    pieces, ppts = {}, {}
    for k, sq in seqs.items():
        po = C.piece_operator(sq, p)
        if po is None:
            pieces[repr(k)] = dict(found=False)
            void.append(f'piece not found: {k!r}')
            continue
        ppts[repr(k)] = C.points(po['c'], p)
        pieces[repr(k)] = dict(found=True, order=po['order'], degree=po['degree'],
                               points=C.describe(ppts[repr(k)]))
    T['pieces'] = round(time.time() - t, 1)
    att = C.attribute(wf, ppts)
    moving = {k for k in C.genuine(wf) if k not in wfs}
    res['MO3'] = dict(holds=all(v['carriers'] and v['modz_covered'] is not False for v in att.values()),
                      uncarried=[(v['pos'], v['cls']) for v in att.values() if not v['carriers']])
    if not res['MO3']['holds']:
        void.append('MO3 carriage miss (instrument)')
    res['MO4'] = dict(attribution=[dict(pos=v['pos'], deg=v['deg'], cls=v['cls'], moves=k in moving,
                                        carriers=v['carriers'], modz_covered=v['modz_covered'])
                                   for k, v in att.items()],
                      carried_degree_by_piece={n: sum(v['deg'] for k, v in att.items() if n in v['carriers'])
                                               for n in ppts},
                      union_pieces_equals_full=set(C.union_genuine(ppts)) == set(C.genuine(wf)))
    res['pieces'] = pieces
    res['void_flags'] = void
    T['total'] = round(time.time() - T0, 1)
    res['timing'] = T
    (DATA / f's1_run_{case}_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: v for k, v in res.items() if k not in ('pieces', 'MO4')}, indent=1, default=str))


def controls(p):
    t = time.time()
    toy = C.toy(p)
    # div: AB-style search operator right-divided by the minimal one; an unrelated minimal operator does not divide
    seqs, xs = C.leg_seqs(p, 400)
    g = C.generic_seq({k: seqs[k] for k in ('a', 'c', 'd')}, p, 2026)
    mn, ab = C.generic_operator(g, p), C.generic_operator_ab(g, p)
    other = C.generic_operator(C.generic_seq({'a': seqs['a'], 'b': seqs['b']}, p, 7), p)
    div = dict(min=(mn['order'], mn['degree']), ab=(ab['order'], ab['degree']),
               divides=rem_zero(ab['c'], mn['c'], p), unrelated_divides=rem_zero(ab['c'], other['c'], p))
    div['pass'] = div['divides'] and not div['unrelated_divides']
    na = C.na_control(p)
    out = dict(prime=p, toy=toy['all_pass'], div=div, n3=all(na[k] for k in ('C1_generic_order8', 'C2_full_noninteger_is_xstar',
                                                                     'C3_carried', 'C4_union_equals_full')),
               time=round(time.time() - t, 1))
    (DATA / f's1_run_controls_{p}.json').write_text(json.dumps(out, indent=1))
    print(out)


def summary():
    def load(f):
        f = DATA / f
        return json.loads(f.read_text()) if f.exists() else None
    ctrl = {p: load(f's1_run_controls_{p}.json') for p in PRIMES}
    runs = {(c, p): load(f's1_run_{c}_{p}.json') for c in CASES for p in PRIMES}
    out = dict(missing=[str(k) for k, v in list(ctrl.items()) + list(runs.items()) if v is None])
    void = list(out['missing'])
    void += [f'controls fail {p}' for p, v in ctrl.items() if v and not (v['toy'] and v['n3'] and v['div']['pass'])]
    for (c, p), r in runs.items():
        if r:
            void += [f'{c} {p}: {f}' for f in r['void_flags']]
    for c in CASES:
        rs = [runs[(c, p)] for p in PRIMES if runs[(c, p)] and 'full' in runs[(c, p)]]
        if len(rs) == 2:
            for tg in ('base', 'shift'):
                if len({(o['order'], o['degree']) for r in rs for o in r['full'] if o['tag'] == tg}) != 1:
                    void.append(f'{c}: minimal (order, degree) differs across primes ({tg})')
            d = [r['MO4']['carried_degree_by_piece'] for r in rs]
            out[f'MO4_{c}_degrees_agree'] = d[0] == d[1]
    ok = [r for r in runs.values() if r and 'MO1' in r]
    out['MO0'] = 'HELD' if not void else 'MISSED'
    out['MO1'] = 'HELD' if all(runs[('A2', p)] and runs[('A2', p)].get('MO1', {}).get('holds') for p in PRIMES) \
        else 'MISSED'
    out['MO2'] = 'HELD' if len(ok) == 4 and all(r['MO2']['holds'] for r in ok) else 'MISSED'
    out['MO2_points'] = {f'{r["case"]} {r["prime"]}': [(x['side'], x['pos'], x['cls'], x['exps'],
                                                        'y' if x['moves'] else 'n') for x in r['MO2']['points']]
                         for r in ok}
    out['n_moving'] = {f'{r["case"]} {r["prime"]}': r['MO2']['n_moving'] for r in ok}
    nmov_case = {c: min((runs[(c, p)] or {}).get('MO2', {}).get('n_moving', 0) for p in PRIMES) for c in CASES}
    out['MO3'] = 'HELD' if len(ok) == 4 and all(r['MO3']['holds'] for r in ok) else 'MISSED'
    if out['MO3'] != 'HELD':
        void.append('MO3 missed (instrument)')
    out['void_flags'] = void
    out['timing'] = {f'{r["case"]} {r["prime"]}': r['timing'] for r in ok}
    out['timing_controls'] = {p: v['time'] for p, v in ctrl.items() if v}
    if void:
        gate = 'VOID'
    elif out['MO2'] == 'MISSED':
        gate = 'STOP'          # a genuine mu-moving point with a non-integer exponent exists
    elif out['MO1'] == 'HELD' and any(n >= 1 for n in nmov_case.values()):
        gate = 'GO'
    else:
        gate = 'BORDERLINE'    # incl. n_moving = 0 in both cases (MO2 vacuous)
    out['gate'] = gate
    (DATA / 's1_run.json').write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'run':
        run(sys.argv[2], int(sys.argv[3]))
    elif m == 'controls':
        controls(int(sys.argv[2]))
    elif m == 'summary':
        summary()
