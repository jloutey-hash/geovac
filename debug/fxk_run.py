"""Record FXK driver (PROPOSED; written before registration, NOT run on M1 / M2r beyond the L <= 300 pilot).

FXK is the follow-up to record FX (docs/predictions/2026-11-17_four_centre_recurrence.md, STOP on FX3).  FX searched
for one recurrence of a GENERIC random combination sum_k w_k c_k(L) of the canonical coefficient sequences of
G_L = F_L/pi^2 = sum_k c_k(L) key_k and found none at R <= 32 on L <= 1300.  FXK asks instead, per canonical key k,
for the minimal bond-angle operator L_k of f_k(x) = sum_L c_k(L) P_L(x), and reads F's structure through them:

  per key      (order, degree) of the minimal operator (fx_run.minimal_windowed: NA's bisection on the window
               L <= LMAX - HOLD, every annihilator at R0..R0+EXTRA verified on ALL rows, gcrd = MO's construction);
  structure    distinct operators (after content_reduce, i.e. monic, content-free), the F_p-rank of the key x L
               coefficient matrix, proportional key classes, pairwise gcrd orders of the distinct operators, and the
               ORDER of their LCLM (= dim of the sum of solution spaces, read as the rank of truncated power-series
               solution bases at two random regular points; the lclm annihilates F(x) = sum_k key_k f_k(x) and the
               generic combination, so it bounds FX's generic order from above; max_k ord L_k bounds it from below);
  contacts     at x_Z*, x_W*, +1, -1, for every distinct operator: present as a factor of the leading coefficient,
               class (s1_controls.points), and the candidate indicial roots (na_modanalyze.local), exactly as FX4 does
               (fx_run.analyse).

Frozen FX files are imported, never edited: fx_run (minimal_windowed, analyse, coeffs, landau, generic), and through
it na_modguess, na_modanalyze, intclass_lib, s1_controls, fx_exact, fx_indep.

CLI
  python debug/fxk_run.py ctrl_pos ab <p>         AB A1 frozen pieces (s1_coeffs_A1), window = all L <= 1000 (as FX's AB
                                                  posctrl); vs MO's per-piece shapes; lclm order vs generic (14, 111)
  python debug/fxk_run.py ctrl_pos k4a1 <p>       K4A1 (a = b) canonical keys through the FX engine, L <= 1200 (FX's
                                                  saved e2e control sequences), window 1000 + 200 held out
  python debug/fxk_run.py ctrl_series <p>         lclm-order instrument on toy operators with known answers
  python debug/fxk_run.py pilot M1 <p> [Lmax]     PILOT (not scored): recompute c_k(L) for L <= Lmax (default 300),
                                                  save, check equal to the prefix of FX's saved sequences, per-key
                                                  search on window Lmax - 50 (+ 50 held out).  NO contact analysis
                                                  (keeps FXK4 unseen)
  python debug/fxk_run.py score <case> <p> --registered     scored (M1, M2r): FX's saved fx_seq_<case>_<p>.json
                                                  (sha256 checked against SEQ_SHA, frozen at registration) or, if
                                                  absent, recomputed (~7000 s); refuses without --registered
  python debug/fxk_run.py summary                 -> data/fxk_run.json (clauses, gate)
"""
from __future__ import annotations

import hashlib
import json
import random
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import flint
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import fx_run as R             # noqa: E402  (frozen; brings G, A, I, S, FX, C, NE)
I, S, A, G = R.I, R.S, R.A, R.G

DATA = HERE / 'data'
PRIMES = R.PRIMES
LMAX, HOLD, RMAX, MARGIN, EXTRA = 1500, 200, 32, 30, 2
SCORED = ('M1', 'M2r')
SERIES_SEEDS = (11, 12)                 # two random regular points for the lclm-order rank
# frozen FX dependencies whose code the per-key pipeline runs (sha256 as in FX's registered record)
DEP_SHA = {'fx_run.py': '63ba74df1653c7824f644aca20e8ebf1ff03c225c9f1e5506c283b5ed06fc5c6',
           'na_modanalyze.py': 'de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005',
           's1_controls.py': '1c50016620c5f6291ac735028d9c87a658d8a9c226f28f08a1c667eba9436687'}
GCRD_CAP = 12                           # pairwise gcrd table only if at most this many distinct operators
# sha256 of FX's saved scored sequences (written by FX's scored run 2026-11-17, untracked at design time; the PM
# commits them with the registration).  score() refuses a file whose hash differs.
SEQ_SHA = {('M1', 1073741441): 'edcf845b3186660750e59e595724967ddb97e6903906ecbedf8d9d2ec4cbc3e1',
           ('M1', 1073741381): 'faf2f135acc4f6b7166af44a24734360d36f39dd0448cb0fa58743b8acc0526e',
           ('M2r', 1073741441): '55b8dfdf882c2a06ff765eb395d95244491f5129063172fa56fe31caa00e8576',
           ('M2r', 1073741381): 'e6374838388e38404cdcf9af7a90a19237a4d50de39d418b06e7f557239d489c'}


# ---------------------------------------------------------------------------------------------- utilities
def sha(fn):
    return hashlib.sha256(Path(fn).read_bytes()).hexdigest()


def fingerprint(c):
    c = I.content_reduce(c)
    h = hashlib.sha1()
    for ck in c:
        h.update(str([int(a) for a in ck.coeffs()]).encode())
        h.update(b'|')
    return h.hexdigest()[:16]


def span_structure(seqs, p, Lmax=None):
    """F_p-rank of the key x L matrix (L <= Lmax) and the classes of keys with proportional sequences."""
    keys = sorted(seqs)
    n = min(len(seqs[k]) for k in keys) if Lmax is None else Lmax + 1
    M = flint.nmod_mat(len(keys), n, [int(seqs[k][L]) % p for k in keys for L in range(n)], p)
    rank = M.rank()
    cls = {}
    for k in keys:
        s = [int(v) % p for v in seqs[k][:n]]
        f = next((v for v in s if v), None)
        if f is None:
            continue
        inv = pow(f, p - 2, p)
        sig = hashlib.sha1(str([v * inv % p for v in s]).encode()).hexdigest()[:16]
        cls.setdefault(sig, []).append(k)
    return dict(n_keys=len(keys), L_used=n - 1, rank=rank, n_prop_classes=len(cls),
                prop_class_sizes=sorted((len(v) for v in cls.values()), reverse=True))


# ---------------------------------------------------------------------------------------------- lclm order
def _shift_ints(ck, x0, p):
    if ck.degree() < 0:
        return []
    q = ck(flint.nmod_poly([x0 % p, 1], p))
    return [int(a) for a in q.coeffs()]


def series_basis(c, p, x0, N):
    """R power-series solutions y = sum_n y_n t^n, t = x - x0, of sum_i c_i(x) d^i y = 0 at a regular point x0
    (c_R(x0) != 0), initial data the unit vectors; returns an (R, N) int64 array mod p."""
    Rr = len(c) - 1
    a = [_shift_ints(ck, x0, p) for ck in c]
    a0 = a[Rr][0] if a[Rr] else 0
    if a0 % p == 0:
        raise ValueError('x0 singular')
    Y = np.zeros((N, Rr), dtype=np.int64)
    for j in range(Rr):
        Y[j, j] = 1

    def ff(k, i):
        r = 1
        for m in range(i):
            r = r * (k - m) % p
        return r
    for n in range(0, N - Rr):
        s = np.zeros(Rr, dtype=np.int64)
        for i in range(Rr + 1):
            for m, am in enumerate(a[i]):
                if am == 0 or (i == Rr and m == 0):
                    continue
                k = n + i - m
                if k < 0 or k < i:
                    continue
                co = am * ff(k, i) % p
                if co:
                    s = (s + co * Y[k]) % p
        den = a0 * ff(n + Rr, Rr) % p
        Y[n + Rr] = (-s % p) * pow(den, p - 2, p) % p
    return Y.T


def lclm_order(ops, p, seeds=SERIES_SEEDS, extra=40):
    """dim of the sum of the operators' solution spaces = order of their LCLM; rank of the stacked truncated series
    bases at random regular points (two seeds; must agree).  N = sum of orders + extra terms."""
    ops = list(ops)
    if not ops:
        return dict(order=0, ranks=[])
    tot = sum(len(c) - 1 for c in ops)
    N = tot + extra
    ranks = []
    for sd in seeds:
        rnd = random.Random(sd)
        while True:
            x0 = rnd.randrange(2, p - 2)
            if all(int(c[-1](flint.nmod(x0, p))) for c in ops):
                break
        rows = np.vstack([series_basis(c, p, x0, N) for c in ops])
        M = flint.nmod_mat(rows.shape[0], N, [int(v) for v in rows.flatten()], p)
        ranks.append(M.rank())
    return dict(order=ranks[0], ranks=ranks, agree=len(set(ranks)) == 1, sum_orders=tot, N=N)


def gcrd_table(ops, p):
    names = sorted(ops)
    if len(names) > GCRD_CAP:
        return dict(skipped=f'{len(names)} distinct operators > cap {GCRD_CAP}')
    tab = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            g = I.gcrd(ops[a], ops[b], p)
            tab[f'{a}|{b}'] = len(g) - 1
    return tab


def right_divides(L1, L2, p):
    r = I.right_rem(L1, L2, p)
    return len(r) == 1 and r[0].degree() < 0


# ---------------------------------------------------------------------------------------------- per-key search
def per_key(seqs, p, Lg, Rmax=RMAX, log=print):
    """Minimal operator of every in-scope key on the window L <= Lg, verified on all rows of the sequence."""
    recs, ops = {}, {}
    for k in sorted(seqs):
        s = seqs[k]
        if not any(int(v) % p for v in s[2:]):
            recs[k] = dict(in_scope=False)
            continue
        t = time.time()
        r = R.minimal_windowed([int(v) % p for v in s], p, Lg, Rmax)
        rec = dict(in_scope=True, found=r['found'])
        if r['found']:
            rec.update(search=list(r['search'][:2]), order=r['order'], degree=r['degree'], verified=r['verified'],
                       min_held_out=r['min_held_out'], fp=fingerprint(r['c']))
            ops[k] = I.content_reduce(r['c'])
        else:
            rec['Qtop'] = {Rr: v[1] for Rr, v in r['trade'].items()}     # no operator of order R and degree <= Qtop
        rec['seconds'] = round(time.time() - t, 1)
        recs[k] = rec
        log(p, k[:70], {x: rec.get(x) for x in ('found', 'order', 'degree', 'search', 'verified', 'seconds')})
    return recs, ops


def structure(recs, ops, p):
    distinct = {}
    for k, c in ops.items():
        distinct.setdefault(recs[k]['fp'], c)
    shapes = {}
    for k, r in recs.items():
        if r.get('found'):
            shapes.setdefault(str((r['order'], r['degree'])), []).append(r['fp'])
    lc = lclm_order(distinct.values(), p)
    return dict(n_in_scope=sum(r['in_scope'] for r in recs.values()), n_found=len(ops), n_distinct=len(distinct),
                shapes={s: dict(n_keys=len(v), n_distinct=len(set(v))) for s, v in sorted(shapes.items())},
                max_order=max((len(c) - 1 for c in distinct.values()), default=None),
                lclm=lc, gcrd=gcrd_table(distinct, p)), distinct


def contacts(distinct, p, marks):
    out = {}
    for fp, c in distinct.items():
        _, at = R.analyse(c, p, marks)
        out[fp] = at
    return out


def contact_reading(cont, names=('xZ', 'xW')):
    """Per contact: carriers (distinct operators with the point present, class M), the union of their half-integer
    candidate roots, and whether 7/2 is a root of every carrier / of some carrier."""
    res = {}
    for m in names:
        car = {fp: at[m] for fp, at in cont.items() if at[m]['present'] and at[m]['cls'] == 'M'}
        pres = {fp: at[m]['cls'] for fp, at in cont.items() if at[m]['present']}
        half = sorted({r for v in car.values() for r in v['roots'] if '/' in r}, key=Fr)
        res[m] = dict(n_present=len(pres), classes=sorted(pres.values()), n_carriers=len(car), half_roots=half,
                      seven_half_every=bool(car) and all('7/2' in v['roots'] for v in car.values()),
                      seven_half_some=any('7/2' in v['roots'] for v in car.values()),
                      lowest_half_per_carrier=sorted({min((r for r in v['roots'] if '/' in r), key=Fr, default='none')
                                                      for v in car.values()}))
    return res


# ---------------------------------------------------------------------------------------------- runs
def _save(name, rec):
    (DATA / name).write_text(json.dumps(rec, indent=1, default=str))
    print(name, json.dumps(rec, default=str)[:3000], flush=True)


def _load_seqs(fn):
    d = json.loads(Path(fn).read_text())
    if 'seqs' in d:
        d = d['seqs']
    return {k: [int(v) for v in s] for k, s in d.items()}


def ctrl_pos(which, p):
    """Positive control: AB's class (a = b), where the minimal generic operator is known: (14, 111), 7/2 at x*."""
    T0 = time.time()
    if which == 'ab':
        seqs = _load_seqs(DATA / f's1_coeffs_A1_{p}.json')
        import aabc_exact as AB
        D1, D2 = AB.params('A1')[:2]
        Lm, Lg = 1000, 1000
        mo = json.loads((DATA / f's1_run_A1_{p}.json').read_text())['pieces']
        mo_shapes = sorted(str((v['order'], v['degree'])) for v in mo.values())
    else:
        seqs = _load_seqs(DATA / f'fx_seqctl_e2e_K4A1_{p}.json')
        pr = R.FX.params('K4A1')
        D1, D2 = pr['D1'], pr['a']
        Lm, Lg = 1200, 1000
        mo_shapes = None
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    seqs = {k: s[:Lm + 1] for k, s in seqs.items()}
    span = span_structure(seqs, p)
    recs, ops = per_key(seqs, p, Lg)
    st, distinct = structure(recs, ops, p)
    cont = contacts(distinct, p, dict(xZ=xs, xW=xs, p1=Fr(1), m1=Fr(-1)))
    cr = contact_reading(cont, ('xZ',))
    fire = contact_reading(contacts(distinct, p, dict(xZ=xs + Fr(1, 7), xW=xs + Fr(1, 7), p1=Fr(1), m1=Fr(-1))), ('xZ',))
    gen = R.minimal_windowed(R.generic(seqs, p, 2026), p, Lg)
    gshape = (gen['order'], gen['degree']) if gen['found'] else None
    divs = {fp: right_divides(gen['c'], c, p) for fp, c in distinct.items()} if gen['found'] else None
    key_shapes = sorted(str((r['order'], r['degree'])) for r in recs.values() if r.get('found'))
    rec = dict(which=which, prime=p, xstar=str(xs), Lmax=Lm, window=Lg, span=span, structure=st,
               per_key={k: v for k, v in recs.items()}, contacts=cont, contact_reading=cr,
               generic_minimal=gshape, generic_right_divisible_by_each=divs,
               mo_piece_shapes=mo_shapes, key_shapes=key_shapes)
    rec['all_found'] = st['n_found'] == st['n_in_scope']
    rec['all_verified'] = all(r.get('verified') for r in recs.values() if r.get('found'))
    # the generic minimal operator annihilates g = sum w_k f_k, whose solutions lie in the sum of the keys' solution
    # spaces: its order is <= the lclm order, and adding it to the keys' operators must not raise the lclm order
    # (a key operator that missed part of f_k's solution space would make this rise)
    lc_plus = lclm_order(list(distinct.values()) + [I.content_reduce(gen['c'])], p) if gen['found'] else None
    rec['lclm_with_generic'] = lc_plus
    rec['generic_le_lclm'] = gshape is not None and gshape[0] <= st['lclm']['order'] and st['lclm']['agree']
    rec['lclm_closed_under_generic'] = lc_plus is not None and lc_plus['order'] == st['lclm']['order'] and lc_plus['agree']
    rec['fire_shifted_contact'] = fire
    shape_of = {r['fp']: (r['order'], r['degree']) for r in recs.values() if r.get('found')}
    rec['carrier_shapes'] = {fp: dict(shape=shape_of[fp], half_roots=[x for x in at['xZ']['roots'] if '/' in x])
                             for fp, at in cont.items() if at['xZ']['present']}
    if gen['found']:            # fire: drop the largest key operator; the generic must then RAISE the lclm order
        big = max(distinct, key=lambda fp: len(distinct[fp]))
        rest = [c for fp, c in distinct.items() if fp != big]
        a, b = lclm_order(rest, p), lclm_order(rest + [I.content_reduce(gen['c'])], p)
        rec['fire_drop_largest'] = dict(dropped=shape_of[big], lclm_rest=a['order'], lclm_rest_plus_generic=b['order'],
                                        caught=b['order'] > a['order'])
    rec['fire_caught'] = fire['xZ']['n_carriers'] == 0 and not fire['xZ']['seven_half_some']
    rec['generic_is_AB'] = gshape == (14, 111)
    rec['xstar_7_2'] = cr['xZ']['seven_half_some']
    if which == 'ab':
        rec['shapes_equal_MO'] = key_shapes == mo_shapes
    rec['pass'] = bool(rec['all_found'] and rec['all_verified'] and rec['generic_le_lclm'] and rec['generic_is_AB']
                       and rec['lclm_closed_under_generic'] and rec['xstar_7_2'] and rec['fire_caught']
                       and rec.get('fire_drop_largest', {}).get('caught')
                       and (which != 'ab' or rec['shapes_equal_MO']))
    rec['seconds'] = round(time.time() - T0, 1)
    _save(f'fxk_ctrl_pos_{which}_{p}.json', rec)


def ctrl_series(p):
    """lclm-order instrument on toy operators: (D - 1), (D - 2) -> 2; (D - 1) twice -> 1; D^2 - 1 with D - 1 -> 2;
    a Legendre-type order 2 (1 - x^2) D^2 - 2x D + 6 with (D - 1) -> 3; the K4A1/AB checks are the real control."""
    z = lambda cs: flint.nmod_poly(cs, p)                          # noqa: E731
    one, D1_, D2_ = [z([-1]), z([1])], [z([-2]), z([1])], [z([-1]), z([0]), z([1])]
    leg = [z([6]), z([0, -2]), z([1, 0, -1])]
    cases = dict(two_exps=([one, D1_], 2), same=([one, one], 1), contained=([D2_, one], 2), legendre_exp=([leg, one], 3),
                 single_leg=([leg], 2))
    res = {k: dict(expected=e, got=lclm_order(ops, p)) for k, (ops, e) in cases.items()}
    for v in res.values():
        v['ok'] = v['got']['order'] == v['expected'] and v['got']['agree']
    res['pass'] = all(v['ok'] for v in res.values() if isinstance(v, dict))
    _save(f'fxk_ctrl_series_{p}.json', res)


def pilot(case, p, Lmax=300, hold=50):
    """NOT the scored run.  Recompute c_k(L), L <= Lmax; save; reproducibility vs the prefix of FX's saved
    sequences; per-key search on the window Lmax - hold.  No contact / singular-point analysis."""
    T0 = time.time()
    cf = R.coeffs(case, p, Lmax, log=lambda *a: None)
    seqs = {repr(k): [int(v) for v in s] for k, s in cf['seqs'].items()}
    (DATA / f'fxk_seq_pilot_{case}_{p}_L{Lmax}.json').write_text(json.dumps(seqs))
    t_coeff = round(time.time() - T0, 1)
    saved = DATA / f'fx_seq_{case}_{p}.json'
    repro = None
    if saved.exists():
        full = json.loads(saved.read_text())
        repro = set(full) == set(seqs) and all(full[k][:Lmax + 1] == seqs[k] for k in seqs)
        del full
    span = span_structure(seqs, p)
    recs, ops = per_key(seqs, p, Lmax - hold)
    found = {k: r for k, r in recs.items() if r.get('found')}
    distinct = {}
    for k, c in ops.items():
        distinct.setdefault(recs[k]['fp'], c)
    lc = lclm_order(distinct.values(), p) if distinct else None
    qtop = next((r['Qtop'] for r in recs.values() if r.get('in_scope') and not r.get('found')), None)
    rec = dict(case=case, prime=p, Lmax=Lmax, window=Lmax - hold, hold=hold, t_coeff=t_coeff,
               reproduces_FX_prefix=repro, span=span, n_in_scope=sum(r['in_scope'] for r in recs.values()),
               n_found=len(found), found_shapes=sorted(str((r['order'], r['degree'])) for r in found.values()),
               found_keys={k: dict(order=r['order'], degree=r['degree'], search=r['search'], verified=r['verified'],
                                   held_out=r['min_held_out']) for k, r in found.items()},
               n_distinct_found=len(distinct), lclm_of_found=lc, qtop_ceiling_example=qtop,
               per_key_seconds=round(sum(r['seconds'] for r in recs.values()), 1),
               note='pilot: no singular-point or indicial analysis performed on M1 (FXK4 kept unseen)')
    rec['seconds'] = round(time.time() - T0, 1)
    _save(f'fxk_ctrl_pilot_{case}_{p}_L{Lmax}.json', rec)


def score(case, p):
    T0 = time.time()
    res = dict(case=case, prime=p, LMAX=LMAX, HOLD=HOLD, RMAX=RMAX, instrument_error=None)
    try:
        for name, want in DEP_SHA.items():          # a changed dependency -> exception -> VOID
            got = sha(HERE / name)
            if got != want:
                raise RuntimeError(f'{name} sha256 {got} != pinned {want}')
        res['dep_sha_ok'] = True
        fn = DATA / f'fx_seq_{case}_{p}.json'
        if fn.exists():
            h = sha(fn)
            res['seq_source'], res['seq_sha'] = str(fn.name), h
            if SEQ_SHA.get((case, p)) != h:
                raise RuntimeError(f'sequence file hash {h} != registered {SEQ_SHA[(case, p)]}')
            seqs = _load_seqs(fn)
        else:
            cf = R.coeffs(case, p, LMAX)
            seqs = {repr(k): [int(v) for v in s] for k, s in cf['seqs'].items()}
            (DATA / f'fxk_seq_{case}_{p}.json').write_text(json.dumps(seqs))
            res['seq_source'] = 'recomputed'
        if any(len(s) != LMAX + 1 for s in seqs.values()):
            raise RuntimeError('sequence length != LMAX + 1')
        res['span'] = span_structure(seqs, p)
        recs, ops = per_key(seqs, p, LMAX - HOLD)
        st, distinct = structure(recs, ops, p)
        lm = R.landau(case)
        cont = contacts(distinct, p, dict(xZ=lm['xZ'], xW=lm['xW'], p1=Fr(1), m1=Fr(-1)))
        res.update(landau={k: str(v) for k, v in lm.items()}, per_key=recs, structure=st, contacts=cont,
                   contact_reading=contact_reading(cont))
    except Exception:                                       # noqa: BLE001  (recorded -> VOID)
        import traceback
        res['instrument_error'] = traceback.format_exc()
        print(res['instrument_error'])
    res['t_total'] = round(time.time() - T0, 1)
    _save(f'fxk_run_{case}_{p}.json', res)


# ---------------------------------------------------------------------------------------------- summary / gate
def _j(name):
    fn = DATA / name
    return json.loads(fn.read_text()) if fn.exists() else None


def fxk0():
    c = {}
    for p in PRIMES:
        for w in ('ab', 'k4a1'):
            c[f'pos_{w}_{p}'] = bool((_j(f'fxk_ctrl_pos_{w}_{p}.json') or {}).get('pass'))
        c[f'series_{p}'] = bool((_j(f'fxk_ctrl_series_{p}.json') or {}).get('pass'))
    c['firetest_gate'] = bool((_j('fxk_ctrl_firetest_gate.json') or {}).get('all_pass'))
    return all(c.values()), c


def _key_at(r, k, m):
    """(shape, present, class, half-integer candidate roots) of key k's operator at mark m in run r; None if the
    key's operator was not found on that prime."""
    pk = r['per_key'].get(k, {})
    if not pk.get('found'):
        return None
    at = r['contacts'][pk['fp']][m]
    return dict(shape=(pk['order'], pk['degree']), present=bool(at['present']), cls=at['cls'],
                half=sorted((x for x in at['roots'] if '/' in x), key=Fr) if at['present'] else [])


def case_summary(case, r1, r2):
    """Clauses and gate of the draft record (docs/predictions/2026-11-17_fxk_per_key_recurrence.md; revised after C1 and the cold review).
    A key's operator is PRESENT at a contact if the contact is a factor of its leading coefficient, in any class
    (M, L, A, ...): FX4's condition.  Class does not enter the gate; the half-integer candidate roots do.
    FXK1 (reported) coverage: keys found on both primes with the same (order, degree); complete = every key.
    FXK2 each contact is present in some found key operator (consistently on both primes).
    FXK3 at each contact, some operator in which it is present has 7/2 among its half-integer candidate roots.
    FXK4 (reported) lclm of the found operators (and whether the two series points agree), max order, distinct
         operators, gcrd table, classes at +-1.
    Cross-prime consistency (VOID on failure): for each contact, every key whose operator is present there on either
    prime is found on both primes with the same (order, degree), present on both, with the same half-integer root
    set.  Shape differences of keys present at no contact are reported (FXK1), not gated.
    Gate: VOID (exception, held-out failure, cross-prime inconsistency) > STOP (a: neither contact present in any
    found operator; b: a contact present and no operator present there has 7/2) > GO (FXK2 and FXK3) > BORDERLINE."""
    void = [f"instrument exception on prime {r['prime']}" for r in (r1, r2) if r.get('instrument_error')]
    if void:
        return dict(void_flags=void, gate='VOID')
    runs = (r1, r2)
    keys = set(r1['per_key']) | set(r2['per_key'])
    scope = [k for k in keys if any(r['per_key'].get(k, {}).get('in_scope') for r in runs)]
    found_both = [k for k in scope if all(r['per_key'].get(k, {}).get('found') for r in runs)]
    same_shape = [k for k in found_both if (r1['per_key'][k]['order'], r1['per_key'][k]['degree']) ==
                  (r2['per_key'][k]['order'], r2['per_key'][k]['degree'])]
    unverified = sorted({k for k in scope for r in runs if r['per_key'].get(k, {}).get('found')
                         and not r['per_key'][k]['verified']})
    short_hold = sorted({k for k in scope for r in runs if r['per_key'].get(k, {}).get('found')
                         and r['per_key'][k]['min_held_out'] < HOLD})
    k1 = dict(n_scope=len(scope), n_found=[r['structure']['n_found'] for r in runs], n_found_both=len(found_both),
              n_same_shape=len(same_shape), shapes=[r['structure']['shapes'] for r in runs])
    k1['complete'] = len(found_both) == len(scope) and len(same_shape) == len(found_both)
    per_contact, inconsistent = {}, []
    for m in ('xZ', 'xW'):
        rows = {k: [_key_at(r, k, m) for r in runs] for k in scope}
        present = [k for k, v in rows.items() if any(x and x['present'] for x in v)]
        for k in present:
            a, b = rows[k]
            if a is None or b is None or a['shape'] != b['shape'] or not (a['present'] and b['present']) \
                    or a['half'] != b['half']:
                inconsistent.append(dict(contact=m, key=k, p1=a, p2=b))
        per_contact[m] = dict(present_keys={k: rows[k][0] for k in present},
                              n_present=len(present),
                              has_7_2=any('7/2' in (rows[k][0] or {}).get('half', []) for k in present),
                              half_roots=sorted({h for k in present for x in rows[k] if x for h in x['half']}, key=Fr))
    k2 = dict(n_present={m: v['n_present'] for m, v in per_contact.items()})
    k2['holds'] = all(v['n_present'] >= 1 for v in per_contact.values())
    k3 = dict(half_roots={m: v['half_roots'] for m, v in per_contact.items()},
              has_7_2={m: v['has_7_2'] for m, v in per_contact.items()})
    k3['holds'] = k2['holds'] and all(v['has_7_2'] for v in per_contact.values())
    lc = [r['structure']['lclm'] for r in runs]
    k4 = dict(lclm_orders=[x['order'] for x in lc], lclm_series_points_agree=[x.get('agree') for x in lc],
              max_orders=[r['structure']['max_order'] for r in runs],
              n_distinct=[r['structure']['n_distinct'] for r in runs], span=[r.get('span') for r in runs],
              pm1=[{fp: {m: at[m]['cls'] for m in ('p1', 'm1')} for fp, at in r['contacts'].items()} for r in runs],
              contacts=per_contact)
    if unverified:
        void.append(f'annihilator failed held-out rows: {unverified[:5]}')
    if short_hold:
        void.append(f'fewer than HOLD held-out rows: {short_hold[:5]}')
    if inconsistent:
        void.append(f'cross-prime inconsistency at a contact: {inconsistent[:3]}')
    rec = dict(FXK1=k1, FXK2=k2, FXK3=k3, FXK4=k4, void_flags=void)
    none_reached = all(v['n_present'] == 0 for v in per_contact.values())
    present_no72 = any(v['n_present'] >= 1 and not v['has_7_2'] for v in per_contact.values())
    if void:
        gate = 'VOID'
    elif none_reached or present_no72:
        gate = 'STOP'
    elif k2['holds'] and k3['holds']:
        gate = 'GO'
    else:
        gate = 'BORDERLINE'
    rec['stop_reason'] = ('no contact reached' if none_reached else
                          'a contact present, no operator present there has 7/2' if present_no72 else None)
    rec['gate'] = gate
    return rec


PRECEDENCE = ('INCOMPLETE', 'VOID', 'STOP', 'BORDERLINE', 'GO')


def summary():
    ok0, checks = fxk0()
    out = dict(FXK0=dict(holds=ok0, failed=[k for k, v in checks.items() if not v]))
    for case in SCORED:
        runs = {p: _j(f'fxk_run_{case}_{p}.json') for p in PRIMES}
        if any(v is None for v in runs.values()):
            out[case] = dict(gate='INCOMPLETE')
            continue
        out[case] = case_summary(case, runs[PRIMES[0]], runs[PRIMES[1]])
    gates = [out[c]['gate'] for c in SCORED]
    if not ok0 and 'INCOMPLETE' not in gates:
        gates.append('VOID')
    out['overall'] = next(g for g in PRECEDENCE if g in gates)
    _save('fxk_run.json', out)


def calib(case, p, Lmax=LMAX, hold=HOLD):
    """Calibration control on an UNSCORED matched geometry (C1/C2 of fx_indep): sequences through FX's engine to
    Lmax (saved as fxk_seq_calib_<case>_<p>.json; never an FX file), then the full scored pipeline, contacts included."""
    if case in SCORED:
        raise SystemExit('calib is for unscored geometries only')
    T0 = time.time()
    res = dict(case=case, prime=p, LMAX=Lmax, HOLD=hold, RMAX=RMAX, instrument_error=None)
    try:
        fn = DATA / f'fxk_seq_calib_{case}_{p}.json'
        if fn.exists():
            seqs = _load_seqs(fn)
            res['t_coeff'] = 'cached'
        else:
            cf = R.coeffs(case, p, Lmax, log=lambda *a: print(*a, flush=True))
            seqs = {repr(k): [int(v) for v in s] for k, s in cf['seqs'].items()}
            fn.write_text(json.dumps(seqs))
            res['t_coeff'] = cf['t_cum']
        res['seq_sha'] = sha(fn)
        res['span'] = span_structure(seqs, p)
        recs, ops = per_key(seqs, p, Lmax - hold, log=lambda *a: print(*a, flush=True))
        st, distinct = structure(recs, ops, p)
        lm = R.landau(case)
        cont = contacts(distinct, p, dict(xZ=lm['xZ'], xW=lm['xW'], p1=Fr(1), m1=Fr(-1)))
        shape_of = {r['fp']: (r['order'], r['degree']) for r in recs.values() if r.get('found')}
        res.update(landau={k: str(v) for k, v in lm.items()}, per_key=recs, structure=st, contacts=cont,
                   contact_reading=contact_reading(cont),
                   carrier_shapes={m: {fp: dict(shape=shape_of[fp], cls=at[m]['cls'],
                                                half_roots=[x for x in at[m]['roots'] if '/' in x])
                                       for fp, at in cont.items() if at[m]['present']} for m in ('xZ', 'xW')},
                   t_search=round(sum(r['seconds'] for r in recs.values()), 1))
    except Exception:                                       # noqa: BLE001
        import traceback
        res['instrument_error'] = traceback.format_exc()
        print(res['instrument_error'])
    res['t_total'] = round(time.time() - T0, 1)
    _save(f'fxk_ctrl_calib_{case}_{p}.json', res)


def firetest_gate():
    """Synthetic run pairs through case_summary.  Each key is given per prime as (found, shape, {mark: at}); a mark
    missing from the dict is absent.  Each case names the wrong answer it must catch."""
    def at(cls='M', roots=('7/2',)):
        return dict(present=True, cls=cls, roots=list(roots))
    absent = dict(present=False, cls=None, roots=[])

    def run(p, keys, verified=True, err=None):
        pk, cont = {}, {}
        for name, (found, shape, marks) in keys.items():
            if not found:
                pk[name] = dict(in_scope=True, found=False)
                continue
            fp = f'{name}_{p}'
            pk[name] = dict(in_scope=True, found=True, order=shape[0], degree=shape[1], verified=verified,
                            min_held_out=HOLD, fp=fp)
            cont[fp] = {m: marks.get(m, absent) for m in ('xZ', 'xW', 'p1', 'm1')}
        st = dict(n_found=sum(1 for v in pk.values() if v['found']), shapes={}, lclm=dict(order=4, agree=True),
                  max_order=4, n_distinct=len(cont))
        return dict(prime=p, instrument_error=err, per_key=pk, contacts=cont, structure=st, span={})
    P1, P2 = PRIMES
    S0 = (15, 77)
    both = dict(xZ=at(), xW=at())

    def pair(k1, k2=None, **kw):
        return run(P1, k1, **kw), run(P2, k2 if k2 is not None else k1)
    cases = dict(
        go=(pair({'c': (True, S0, both)}), 'GO'),
        go_incomplete=(pair({'c': (True, S0, both), 'u': (False, None, {})}), 'GO'),
        go_two_carriers=(pair({'z': (True, S0, dict(xZ=at())), 'w': (True, (10, 44), dict(xW=at()))}), 'GO'),
        complete_no_7_2=(pair({'c': (True, S0, dict(xZ=at(roots=('3/2',)), xW=at(roots=('3/2',))))}), 'STOP'),
        one_contact_no_7_2=(pair({'c': (True, S0, dict(xZ=at(), xW=at(roots=('3/2',))))}), 'STOP'),
        no_contact=(pair({'c': (True, S0, {})}), 'STOP'),
        one_contact_only=(pair({'c': (True, S0, dict(xZ=at()))}), 'BORDERLINE'),
        present_L_no_7_2=(pair({'c': (True, S0, dict(xZ=at('L', ()), xW=at()))}), 'STOP'),
        present_A_no_7_2=(pair({'c': (True, S0, dict(xZ=at(), xW=at('A', ())))}), 'STOP'),
        present_L_with_7_2=(pair({'c': (True, S0, dict(xZ=at('L'), xW=at()))}), 'GO'),
        xZ_p1_only_xW_p2_only=(pair({'c': (True, S0, dict(xZ=at()))}, {'c': (True, S0, dict(xW=at()))}), 'VOID'),
        carrier_shape_differs=(pair({'c': (True, S0, both)}, {'c': (True, (15, 78), both)}), 'VOID'),
        carrier_roots_differ=(pair({'c': (True, S0, both)},
                                   {'c': (True, S0, dict(xZ=at(roots=('3/2',)), xW=at()))}), 'VOID'),
        carrier_found_one_prime=(pair({'c': (True, S0, both)}, {'c': (False, None, {})}), 'VOID'),
        noncarrier_shape_differs=(pair({'c': (True, S0, both), 'n': (True, (6, 21), {})},
                                       {'c': (True, S0, both), 'n': (True, (6, 22), {})}), 'GO'),
        unverified=((run(P1, {'c': (True, S0, both)}, verified=False), run(P2, {'c': (True, S0, both)})), 'VOID'),
        exception=((run(P1, {'c': (True, S0, both)}, err='boom'), run(P2, {'c': (True, S0, both)})), 'VOID'))
    res = {}
    for name, ((a, b), want) in cases.items():
        got = case_summary('T', a, b)['gate']
        res[name] = dict(expected=want, got=got, ok=got == want)
    res['notes'] = ('present_L_with_7_2 is a consistency case only: s1_controls.points classes a point with a '
                    'half-integer root as M, so L/A with 7/2 does not arise; the gate reads roots, not classes. '
                    'Series-point (seed) disagreement of the LCLM is reported in FXK4, not gated.')
    res['all_pass'] = all(v['ok'] for v in res.values() if isinstance(v, dict))
    _save('fxk_ctrl_firetest_gate.json', res)

if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'ctrl_pos':
        ctrl_pos(sys.argv[2], int(sys.argv[3]))
    elif m == 'ctrl_series':
        ctrl_series(int(sys.argv[2]))
    elif m == 'pilot':
        if sys.argv[2] in SCORED and len(sys.argv) > 4 and int(sys.argv[4]) > 300:
            raise SystemExit('pilot is capped at L <= 300 on a scored case')
        pilot(sys.argv[2], int(sys.argv[3]), *(int(a) for a in sys.argv[4:5]))
    elif m == 'score':
        if '--registered' not in sys.argv:
            raise SystemExit('scored run: register record FXK first, then pass --registered')
        score(sys.argv[2], int(sys.argv[3]))
    elif m == 'summary':
        summary()
    elif m == 'calib':
        calib(sys.argv[2], int(sys.argv[3]))
    elif m == 'firetest_gate':
        firetest_gate()
    else:
        raise SystemExit(m)
