"""Peeling diagnostic (UNSCORED; C1 calibration geometry and AB controls only).  fxk_design_memo.md section 7.

For a hard key k (no operator within FXK's budget) and a found operator A: g = (w^Ra A) f_k, w = 1 - x^2, computed
on the Legendre-coefficient sequence through A's T-form (x and T = w d/dx act locally on the coefficients).  If
gcrd(A, L_k) has order s, g's minimal operator has order ord L_k - s, so it may fit the window where L_k does not.
The search on g is FX's unchanged (fx_run.minimal_windowed, R <= 32, window L <= 1300 where the rows allow, the
rest held out).  The composite M o (w^Ra A) annihilates f_k (checked on the rows, by code separate from the search).

Refuses M1 / M2r: only fxk_seq_calib_C1_<p>.json and s1_coeffs_A1_<p>.json are read.  No frozen file is written.

CLI   python debug/peel_diag.py ctrl <p>     AB controls   -> data/peel_diag_ctrl_ab_<p>.json
      python debug/peel_diag.py run <p>      C1            -> data/peel_diag_C1_<p>.json (+ .jsonl progress, resumable)
      python debug/peel_diag.py summary      both primes   -> data/peel_diag_C1_summary.json
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from fractions import Fraction
from math import comb
from pathlib import Path

import flint
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import fxk_run as K            # noqa: E402
R, I, A = K.R, K.I, K.A
DATA = HERE / 'data'
PRIMES = R.PRIMES
WINDOW, RMAX, MIN_HOLD = 1300, 32, 60
_NS = {'Fraction': Fraction, '__builtins__': {}}


# ------------------------------------------------------------------------------- sequence-side operator application
class Leg:
    def __init__(self, p: int, n: int):
        self.p = p
        inv = lambda v: pow(int(v) % p, p - 2, p)                     # noqa: E731
        ls = range(n)
        self.xa = np.array([l * inv(2 * l - 1) % p for l in ls], dtype=np.int64)
        self.xb = np.array([(l + 1) * inv(2 * l + 3) % p for l in ls], dtype=np.int64)
        self.ta = np.array([(l - 1) * l * inv(2 * l - 1) % p for l in ls], dtype=np.int64)
        self.tb = np.array([(l + 1) * (l + 2) * inv(2 * l + 3) % p for l in ls], dtype=np.int64)

    def X(self, s):
        n, p = len(s) - 1, self.p
        out = self.xb[:n] * s[1:n + 1] % p
        out[1:] = (out[1:] + self.xa[1:n] * s[:n - 1]) % p
        return out

    def T(self, s):
        n, p = len(s) - 1, self.p
        out = self.tb[:n] * s[1:n + 1] % p
        out[1:] = (out[1:] - self.ta[1:n] * s[:n - 1]) % p
        return out


def t_form(c, p: int):
    """w^R sum_k c_k D^k = sum_k q_k(x) T^k; returns the q_k (nmod_poly)."""
    Rr = len(c) - 1
    w = flint.nmod_poly([1, 0, -1], p)
    Tp = A.T_powers(Rr, p)
    Lp = [w ** Rr * ck for ck in c]
    q = [None] * (Rr + 1)
    for k in range(Rr, -1, -1):
        qk, r = divmod(Lp[k], w ** k)
        assert r.degree() < 0
        q[k] = qk
        for m, b in Tp[k].items():
            Lp[m] = Lp[m] - qk * b
    assert all(x.degree() < 0 for x in Lp)
    return q


def apply_seq(c, seq, p: int):
    """Legendre coefficients of (w^R L) f, on every row the T-form reaches."""
    q = [[int(a) for a in qk.coeffs()] for qk in t_form(c, p)]
    lg = Leg(p, len(seq) + 2)
    v = np.array([int(x) % p for x in seq], dtype=np.int64)
    Tv = [v]
    for _ in range(len(q) - 1):
        Tv.append(lg.T(Tv[-1]))
    parts = []
    for i, qi in enumerate(q):
        if not qi:
            continue
        r = (qi[-1] * Tv[i]) % p
        for a in reversed(qi[:-1]):
            r = lg.X(r)
            r = (r + a * Tv[i][:len(r)]) % p
        parts.append(r)
    n = min(len(r) for r in parts)
    tot = np.zeros(n, dtype=np.int64)
    for r in parts:
        tot = (tot + r[:n]) % p
    return tot


def wpow_op(c, p):
    w = flint.nmod_poly([1, 0, -1], p)
    return [w ** (len(c) - 1) * ck for ck in c]


def compose(a, b, p):
    """a o b in D-form."""
    out = [I.zp(p) for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        if ai.degree() < 0:
            continue
        for k, t in enumerate(I.d_power_compose(b, i, p)):
            out[k] += ai * t
    return out


def shape(c):
    return [len(c) - 1, max(x.degree() for x in c)]


def peel_once(Aop, seq, p, window=WINDOW):
    """Apply A, search.  Returns record (+ quotient operator M under 'c' when found)."""
    g = [int(v) for v in apply_seq(Aop, seq, p)]
    n = len(g)
    Lg = min(window, n - 1 - MIN_HOLD)
    rec = dict(rows=n, window=Lg, zero=not any(g))
    if rec['zero']:
        return rec
    t = time.time()
    r = R.minimal_windowed(g, p, Lg, RMAX)
    rec.update(found=r['found'], seconds=round(time.time() - t, 1))
    if r['found']:
        rec.update(q_order=r['order'], q_degree=r['degree'], search=list(r['search'][:2]), verified=r['verified'],
                   held_out=r['min_held_out'])
        rec['c'] = r['c']
    return rec


def composite_check(M, Aop, seq, p):
    """N = M o (w^Ra A); N f = 0 on the rows (separate code path from the search)."""
    N = I.content_reduce(compose(M, wpow_op(Aop, p), p))
    res = apply_seq(N, seq, p)
    return N, dict(shape=shape(N), rows=len(res), annihilates=not res.any())


# ------------------------------------------------------------------------------- controls (AB, A1)
def ctrl(p):
    T0 = time.time()
    seqs = K._load_seqs(DATA / f's1_coeffs_A1_{p}.json')
    seqs = {k: s[:1001] for k, s in seqs.items()}
    pos = json.loads((DATA / f'fxk_ctrl_pos_ab_{p}.json').read_text())['per_key']
    SH = dict(big=(10, 76), f4=(4, 26), s6=(6, 36), n5=(5, 16), e5=(5, 28))     # fingerprints are prime-dependent
    ops = {}
    for name, sh in SH.items():
        k = sorted(kk for kk, v in pos.items() if (v['order'], v['degree']) == sh)[0]
        r = R.minimal_windowed(seqs[k], p, 1000, RMAX)
        c = I.content_reduce(r['c'])
        assert K.fingerprint(c) == pos[k]['fp'] and tuple(shape(c)) == sh, (name, shape(c))
        ops[name] = (k, c)
    one = flint.nmod_poly([1], p)
    toy = [flint.nmod_poly([3, 1], p), one]                     # D + x + 3: divides nothing here
    out = dict(prime=p, tests={})

    def run(tag, Aop, target, expect_drop, window=None):
        k, L = ops[target]
        s = len(I.gcrd(Aop, L, p)) - 1 if len(Aop) > 1 else 0
        g_ord = len(L) - 1 - s
        rec = peel_once(Aop, seqs[k], p, window if window else 10 ** 9)
        M = rec.pop('c', None)
        rec.update(target_shape=shape(L), peel_shape=shape(Aop), gcrd_order=s, expected_q_order=g_ord,
                   expect_drop=expect_drop)
        if M is not None:
            N, cc = composite_check(M, Aop, seqs[k], p)
            rec['composite'] = cc
            rec['composite_right_divisible_by_target'] = K.right_divides(N, L, p)
            rec['ok'] = bool(rec['q_order'] == g_ord and rec['verified'] and cc['annihilates']
                             and rec['composite_right_divisible_by_target']
                             and ((rec['q_order'] < len(L) - 1) == expect_drop))
        else:
            rec['ok'] = False
        out['tests'][tag] = rec
        print(tag, rec, flush=True)

    run('P1_peel_4_26_from_10_76', ops['f4'][1], 'big', True)          # right factor: quotient order 6
    run('P2_peel_4_26_from_6_36', ops['f4'][1], 's6', True)            # right factor: quotient order 2
    run('P3_peel_5_28_from_10_76', ops['e5'][1], 'big', True)          # shares order 4, does not divide: 6
    run('N1_peel_4_26_from_5_16', ops['f4'][1], 'n5', False)           # gcrd 0: order stays 5
    run('N2_peel_toy_from_6_36', toy, 's6', False)                     # gcrd 0: order stays 6
    run('N3_peel_toy_from_10_76', toy, 'big', False)                   # gcrd 0: order stays 10
    # fire: a composite built with the WRONG quotient must fail the row check
    k, L = ops['big']
    bad = I.content_reduce(compose(toy, wpow_op(ops['f4'][1], p), p))
    out['fire_wrong_quotient_caught'] = bool(apply_seq(bad, seqs[k], p).any())
    out['pass'] = all(v['ok'] for v in out['tests'].values()) and out['fire_wrong_quotient_caught']
    out['seconds'] = round(time.time() - T0, 1)
    (DATA / f'peel_diag_ctrl_ab_{p}.json').write_text(json.dumps(out, indent=1, default=str))
    print('PASS' if out['pass'] else 'FAIL', out['seconds'])


# ------------------------------------------------------------------------------- C1
def key_type(k: str) -> str:
    key = eval(k, _NS)                                                 # noqa: S307  (repo data)
    lets = [a[0] for a in key[1]]
    if 'B' in lets:
        return 'depth2_B'
    if len(lets) >= 2:
        return 'E1xE1'
    if lets == ['ln']:
        return 'ln'
    if lets == ['E1']:
        return 'single_E1'
    if not lets:
        return 'pure_exp'
    return 'other:' + ','.join(lets)


def classes(seqs, keys, p):
    cls = {}
    for k in keys:
        s = [int(v) % p for v in seqs[k]]
        f = next(v for v in s if v)
        inv = pow(f, p - 2, p)
        sig = hashlib.sha1(str([v * inv % p for v in s]).encode()).hexdigest()[:16]
        cls.setdefault(sig, []).append(k)
    return cls


def run(p):
    T0 = time.time()
    seqs = K._load_seqs(DATA / f'fxk_seq_calib_C1_{p}.json')
    cal = json.loads((DATA / f'fxk_ctrl_calib_C1_{p}.json').read_text())['per_key']
    found = {k: v for k, v in cal.items() if v.get('found')}
    unfound = sorted(k for k, v in cal.items() if v.get('in_scope') and not v.get('found'))
    byfp = {}
    for k, v in found.items():
        byfp.setdefault(v['fp'], []).append(k)
    ops = {}
    for fp, ks in sorted(byfp.items()):
        k = sorted(ks)[0]
        r = R.minimal_windowed(seqs[k], p, K.LMAX - K.HOLD, RMAX)
        c = I.content_reduce(r['c'])
        assert K.fingerprint(c) == fp, (fp, K.fingerprint(c))
        ops[fp] = c
        print('op', fp, shape(c), flush=True)
    cls = classes(seqs, unfound, p)
    prog = DATA / f'peel_diag_C1_{p}.jsonl'
    done = {}
    if prog.exists():
        for line in prog.read_text().splitlines():
            d = json.loads(line)
            done[(d['class'], d['op'])] = d
    order = sorted(ops, key=lambda fp: (len(ops[fp]), max(x.degree() for x in ops[fp])))
    quot = {}
    for sig, ks in sorted(cls.items(), key=lambda kv: kv[1][0]):
        rep = ks[0]
        for fp in order:
            if (sig, fp) in done:
                continue
            rec = peel_once(ops[fp], seqs[rep], p)
            M = rec.pop('c', None)
            rec.update({'class': sig, 'op': fp, 'rep': rep, 'peel_shape': shape(ops[fp])})
            if M is not None:
                N, cc = composite_check(M, ops[fp], seqs[rep], p)
                rec['composite'] = cc
                quot[(sig, fp)] = N
            with prog.open('a') as f:
                f.write(json.dumps(rec, default=str) + '\n')
            done[(sig, fp)] = rec
            print(p, key_type(rep), rep[:60], fp, {x: rec.get(x) for x in
                  ('found', 'q_order', 'q_degree', 'search', 'verified', 'held_out', 'seconds')}, flush=True)
    res = dict(prime=p, case='C1', window=WINDOW, RMAX=RMAX, n_unfound=len(unfound), n_classes=len(cls),
               ops={fp: shape(c) for fp, c in ops.items()}, classes={})
    for sig, ks in cls.items():
        hits = {fp: done[(sig, fp)] for fp in order if done[(sig, fp)].get('found')}
        ent = dict(keys=ks, types=sorted({key_type(k) for k in ks}), n_keys=len(ks),
                   findable=bool(hits), zero_under=[fp for fp in order if done[(sig, fp)].get('zero')],
                   hits={fp: {x: h.get(x) for x in ('q_order', 'q_degree', 'search', 'verified', 'held_out',
                                                    'peel_shape', 'composite')} for fp, h in hits.items()})
        # the key's own operator: gcrd of the composites available in this process (>= 2 gives a real reduction)
        Ns = [quot[(sig, fp)] for fp in hits if (sig, fp) in quot]
        if Ns:
            g = Ns[0]
            for o in Ns[1:]:
                g = I.gcrd(g, o, p)
            g = I.content_reduce(g)
            ent['key_operator_candidate'] = dict(shape=shape(g), n_composites=len(Ns), fp=K.fingerprint(g),
                                                 annihilates=not apply_seq(g, seqs[ks[0]], p).any())
        res['classes'][sig] = ent
    res['seconds'] = round(time.time() - T0, 1)
    (DATA / f'peel_diag_C1_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    print('done', p, res['seconds'])


def summary():
    rs = [json.loads((DATA / f'peel_diag_C1_{p}.json').read_text()) for p in PRIMES]
    out = dict(case='C1', primes=list(PRIMES), by_type={}, classes={})
    keyres = {}
    for r in rs:
        for sig, e in r['classes'].items():
            for k in e['keys']:
                keyres.setdefault(k, []).append(e)
    for k, es in sorted(keyres.items()):
        t = key_type(k)
        b = out['by_type'].setdefault(t, dict(n_unfound=0, findable_both=0, findable_one=0))
        b['n_unfound'] += 1
        nf = sum(e['findable'] for e in es)
        b['findable_both'] += nf == 2
        b['findable_one'] += nf == 1
        out['classes'][k] = dict(type=t, findable=[e['findable'] for e in es],
                                 hits=[{fp: (h['peel_shape'], [h['q_order'], h['q_degree']], h['verified'],
                                             h['composite']['annihilates']) for fp, h in e['hits'].items()}
                                       for e in es],
                                 key_op=[e.get('key_operator_candidate') for e in es])
    out['n_unfound'] = len(keyres)
    out['n_findable_both'] = sum(v['findable_both'] for v in out['by_type'].values())
    (DATA / 'peel_diag_C1_summary.json').write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: v for k, v in out.items() if k != 'classes'}, indent=1))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'ctrl':
        ctrl(int(sys.argv[2]))
    elif cmd == 'run':
        run(int(sys.argv[2]))
    elif cmd == 'summary':
        summary()
    else:
        raise SystemExit(__doc__)
