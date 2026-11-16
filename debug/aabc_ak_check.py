"""Record AK (docs/predictions/2026-11-16_aa_bc_class_decided.md): record AB's P4 re-run on primes that split the
mu-carrying components' rates at infinity.  Verdict logic is IC's (intclass_score.operator_sets / component_verdict,
intclass_lib.peel) and the erratum check is IN's (intclass_in_check.analyse_full), all called unchanged; the
coefficients are AB's frozen engine (aabc_exact.pieces / combine), exactly as aabc_run.p4.

  python debug/aabc_ak_check.py edge <p>    -> the slope-one edge of each component's w-operator (no verdicts)
  python debug/aabc_ak_check.py run <p>     -> debug/data/aabc_ak_<p>.json
  python debug/aabc_ak_check.py score       -> debug/data/aabc_ak_score.json

Controls (AK-0) are IC's and IN's own drivers, unchanged, on the new primes:
  python debug/intclass_controls.py <p>; python debug/intclass_in_check.py control <p>
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X                  # noqa: E402
import intclass_in_check as INC         # noqa: E402
import intclass_lib as I                # noqa: E402
import intclass_score as S              # noqa: E402
import na_exact as E                    # noqa: E402

DATA = HERE / 'data'
AK_PRIMES = (1073741441, 1073741381)    # 615, -105, -69, 11 all squares; = 1 mod 4
AB_PRIMES = (1073741789, 1073741833)
CASE = 'A1'


def components(p, L=1200):
    K = E.Field(p)
    pc = X.pieces(CASE, L, K)
    f = X.combine(pc, X.params(CASE)[4], K)
    na_keys = set(k for q in pc['NA'] for k in q)
    keys = sorted(set(k for q in f for k in q) - na_keys, key=repr)
    return {k: [int(q.get(k, flint.nmod(0, p))) for q in f] for k in keys}


def edge(p):
    D1, D2 = X.params(CASE)[:2]
    a0 = I.to_fp((D1 ** 2 + D2 ** 2) / (2 * D1 * D2), p)
    a2 = I.to_fp(-1 / (2 * D1 * D2), p)
    out = {}
    for k, seq in components(p).items():
        ops = S.operator_sets(seq, p, 24)
        cw = I.pullback_quadratic(ops['A'], a0, a2, p)
        inf = I.infinity_data(cw, p)
        out[repr(k)] = dict(edge_degree=inf['edge_degree'], unsplit=inf['lam_unsplit'],
                            roots=[(I.signed(r, p), m) for r, m in inf['lam_roots']])
        print(repr(k), out[repr(k)], flush=True)
    (DATA / f'aabc_ak_edge_{p}.json').write_text(json.dumps(out, indent=1, default=str))


def run(p):
    D1, D2 = X.params(CASE)[:2]
    T0 = time.time()
    comps = components(p)
    res = dict(case=CASE, prime=p, L=1200, coeff_time=round(time.time() - T0, 1), components={})
    print(CASE, p, len(comps), 'mu-carrying components', res['coeff_time'], 's', flush=True)
    for k, seq in comps.items():
        t1 = time.time()
        ops = S.operator_sets(seq, p, 24)
        if ops is None:
            rec = dict(verdict='NOT-FOUND')
        else:
            rA = INC.analyse_full(ops['A'], D1, D2, p)
            rB = INC.analyse_full(ops['B'], D1, D2, p) if ops['B'] is not None else None
            rec = dict(R0=ops['R0'], found={str(a): b for a, b in ops['found'].items()}, A=rA, B=rB,
                       verdict=S.component_verdict(rA, rB))
            rec['erratum_pass'] = (rec['verdict'] != 'STUCK') or (rA['erratum']['passed'] and rB['erratum']['passed'])
        rec['time'] = round(time.time() - t1, 1)
        res['components'][repr(k)] = rec
        print(p, repr(k), rec['verdict'], rec.get('erratum_pass'), (rec.get('A') or {}).get('residual'),
              rec['time'], 's', flush=True)
        (DATA / f'aabc_ak_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    res['time'] = round(time.time() - T0, 1)
    (DATA / f'aabc_ak_{p}.json').write_text(json.dumps(res, indent=1, default=str))


def score():
    out = {}
    ctrl = all((DATA / f).exists() and json.loads((DATA / f).read_text())['all_pass'] for p in AK_PRIMES
               for f in (f'intclass_controls_{p}.json', f'intclass_in_controls_{p}.json'))
    out['AK-0'] = 'HELD' if ctrl else 'MISSED'
    files = [DATA / f'aabc_ak_{p}.json' for p in AK_PRIMES] + [DATA / f'aabc_p4_{p}.json' for p in AB_PRIMES]
    if not all(f.exists() for f in files):
        out.update({'AK-1': 'UNDECIDED', 'AK-2': 'MISSING', 'AK-3': 'MISSING', 'gate': 'BORDERLINE' if ctrl else 'VOID',
                    'missing': [str(f) for f in files if not f.exists()]})
        (DATA / 'aabc_ak_score.json').write_text(json.dumps(out, indent=1, default=str))
        print(out)
        return
    per = [json.loads(f.read_text())['components'] for f in files[:2]]
    ab = [json.loads(f.read_text())['components'] for f in files[2:]]
    keys = sorted(set(per[0]) | set(per[1]))
    v = {k: [(d.get(k, {}).get('verdict', 'NOT-FOUND'), d.get(k, {}).get('erratum_pass')) for d in per] for k in keys}
    out['components'] = v
    out['n_components'] = len(keys)
    clean = [k for k, ((a, ea), (b, eb)) in v.items() if a == b == 'STUCK' and ea and eb]
    out['clean_stuck_both'] = len(clean)
    if clean:
        out['AK-1'] = 'HELD'
    elif all(a == b == 'TRIANGULAR' for (a, _), (b, _) in v.values()):
        out['AK-1'] = 'MISSED'
    else:
        out['AK-1'] = 'UNDECIDED'
    flags = {k: sorted({f for d in per for s in ('A', 'B') for f in (((d.get(k) or {}).get(s) or {}).get('residual')
                                                                      or {}).get('flags', []) if 'lambda edge' in f})
             for k in keys}
    out['edge_flags'] = {k: f for k, f in flags.items() if f}
    found_all = all(a != 'NOT-FOUND' and b != 'NOT-FOUND' for (a, _), (b, _) in v.values())
    out['AK-2'] = 'HELD' if (len(keys) == 15 and found_all and not out['edge_flags']) else 'MISSED'

    def order(d, k):
        r = (d.get(k, {}).get('A') or {}).get('residual')
        return r['order'] if r else None
    mism = {k: ([order(d, k) for d in per], [order(d, k) for d in ab]) for k in keys
            if len({order(d, k) for d in per + ab}) != 1}
    out['order_mismatch'] = mism
    stuckish = ('STUCK', 'STUCK-UNCONFIRMED')     # "STUCK (any flag)", as the record states
    out['AK-3'] = 'HELD' if (not mism and all(a in stuckish and b in stuckish for (a, _), (b, _) in v.values()))         else 'MISSED'
    if out['AK-0'] != 'HELD':
        gate = 'VOID'
    elif out['AK-1'] == 'HELD' and out['AK-2'] == 'HELD' and out['AK-3'] == 'HELD':
        gate = 'GO'
    elif out['AK-1'] == 'MISSED':
        gate = 'STOP'
    else:
        gate = 'BORDERLINE'
    out['gate'] = gate
    (DATA / 'aabc_ak_score.json').write_text(json.dumps(out, indent=1, default=str))
    for k in ('AK-0', 'n_components', 'clean_stuck_both', 'AK-1', 'edge_flags', 'AK-2', 'order_mismatch', 'AK-3',
              'gate'):
        print(k, out[k])


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'edge':
        edge(int(sys.argv[2]))
    elif m == 'run':
        run(int(sys.argv[2]))
    elif m == 'score':
        score()
