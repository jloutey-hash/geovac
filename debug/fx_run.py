"""Record FX scoring driver (PROPOSED; written before registration, NOT run on M1 / M2r).
Design: debug/fx_design_memo.md sections 5, 7, 8.  Object: G_L = F_L/pi^2 of the four-centre 1s integral (XY|ZW) on
FC's matched slice (debug/fx_exact.py), rewritten to the canonical basis (debug/fx_indep.py, route (a)).

Pipeline per (case, prime), all mod p (record AK's primes, < 2^31 as intclass_lib requires):
  1. coefficients L = 0..LMAX (Engine mod p) -> canonical() per L -> q-sequences {key: [q_L]}.
  2. FX1/FX2 bookkeeping (no reading by the operator): per key the L >= 2 where q_L = 0 mod p; per L whether some
     depth-2 key is nonzero; the atom alphabet test (alphabet(case), the section-3 rule mechanised).
  3. FX3: for seeds 2026 and 7, the generic sequence sum_c w_c q_{L,c}; NA's bisection search (na_modguess, unchanged)
     on the guess window L <= LMAX - HOLD for the smallest order R <= RMAX; every annihilator at R0..R0+2 verified on
     ALL rows to LMAX (>= HOLD rows beyond the window); minimal operator = gcrd of them (intclass_lib.gcrd, as MO).
  4. FX4: at x_Z*, x_W* (fc_engine.landau): class (s1_controls.points) and the indicial roots among the candidates
     (na_modanalyze.local, AB's routine): '7/2' a root at both points.
  5. FX5: genuine (class M/X/L/P) weight-fixed finite points of the minimal operator; positions of linear ones
     reconstructed by CRT over the two primes (summary).
Clauses and gate: summary() (memo section 8).

CLI
  python debug/fx_run.py score <case> <p> --registered   scored run (M1, M2r) -> data/fx_run_<case>_<p>.json
                                                        (+ data/fx_seq_<case>_<p>.json); refuses without the flag
  python debug/fx_run.py summary                        -> data/fx_run.json (clauses, gate)
  python debug/fx_run.py posctrl na <p>                 NA at S1's N3 geometry vs the frozen s1_controls.json
  python debug/fx_run.py posctrl ab <p>                 AB A1 frozen coefficients (s1_coeffs_A1) vs frozen W0
  python debug/fx_run.py posctrl k4a1 <p>               K4A1 (a = b) through the FX engine + canonical, L <= 1200
  python debug/fx_run.py ctrl_alphabet                  alphabet + L-independence on C1, C2 (unscored), L <= 60
  python debug/fx_run.py timing                         canonical-rewrite cost on C1 (unscored), L <= 300
  python debug/fx_run.py e2e_run <p>                    K4A1 through score() (L <= 1200) -> data/fx_runctl_e2e_K4A1_<p>.json
  python debug/fx_run.py e2e_summary                    -> data/fx_runctl_e2e_k4a1.json (as_expected)
  python debug/fx_run.py firetest_fx1                   planted late/sporadic/terminal keys -> data/fx_runctl_firetest_fx1.json
FX0 = the pass fields of fx_gate_controls.json, fx_indep_controls.json and the fx_runctl_*.json controls (fx0_check).
"""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction as Fr
from math import isqrt
from pathlib import Path

import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import na_exact as NE           # noqa: E402
import na_modguess as G         # noqa: E402
import na_modanalyze as A       # noqa: E402
import intclass_lib as I        # noqa: E402
import s1_controls as S         # noqa: E402
import fx_reducer as X          # noqa: E402
import fx_exact as FX           # noqa: E402
import fx_indep as C            # noqa: E402  (adds the unscored control geometries C1, C2 to FX.CASES)

DATA = HERE / 'data'
PRIMES = (1073741441, 1073741381)
SEEDS = (2026, 7)
LMAX, HOLD, RMAX, MARGIN, EXTRA = 1500, 200, 32, 30, 2
SCORED = ('M1', 'M2r')
CAND = {str(e): e for e in range(-3, 14)}
CAND.update({f'{2 * h + 1}/2': Fr(2 * h + 1, 2) for h in range(-3, 12)})


def landau(case):
    p = FX.params(case)
    D1, a, b = p['D1'], p['a'], p['b']
    xZ, xW = (D1 ** 2 + a ** 2) / (2 * D1 * a), (D1 ** 2 + b ** 2) / (2 * D1 * b)
    rZ, rW = min(D1, a) / max(D1, a), min(D1, b) / max(D1, b)
    return dict(xZ=xZ, xW=xW, dom='Z' if rZ > rW else ('W' if rW > rZ else 'tie'))


# ---------------------------------------------------------------------------------------------- FX1 alphabet
def structure(case, Ls=range(2, 7)):
    """Per panel: the (lam, E) pairs the engine puts on the reducer, the poles, the natural factors; psi rates."""
    K = NE.Field(C.P_AUX)
    eng = FX.Engine(case, K)
    Lm = max(Ls)
    it, kt = NE.bessel_tables(K, Lm + 2)
    Qs = list(eng._Q_steps(Lm)) if eng.matched else None
    out, psi_rates = {}, set()
    for L in Ls:
        for inner in (True, False):
            hi, lo = eng._rho1_parts(L, inner, it, kt)
            psi_rates |= {b for (_, b) in list(hi) + list(lo)}
        for panel in eng.panels:
            d = out.setdefault(panel, dict(pairs=set(), poles={Fr(0)}, factors=set()))
            fts = eng.f_terms_matched(L, Qs[L], panel) if eng.matched else eng.f_terms_equal(L, panel, it, kt)
            pts, _, _ = eng.psi(L, panel, it, kt)
            for (_, lam, N, Ii, pole, Kp, E) in fts:
                d['factors'].update(E)
                if pole is not None and Kp:
                    d['poles'].add(Fr(pole))
                for (_, b, _, E2) in pts:
                    d['factors'].update(E2)
                    d['pairs'].add((Fr(lam + b), tuple(sorted(E + E2))))
    return out, psi_rates


def _primes_of(v):
    v = abs(Fr(v))
    s = set()
    for n in (v.numerator, v.denominator):
        if n > 1:
            s |= {int(q) for q, _ in flint.fmpz(n).factor()}
    return s


def alphabet(case):
    st, psi_rates = structure(case)
    D1 = FX.params(case)['D1']
    Bs, E1s, lnp = set(), set(), set()
    for (Aa, Bb), d in st.items():
        fac = d['factors']
        poles = set(d['poles']) | {Fr(-dd) / c for c, dd in fac}
        pairs, todo = set(d['pairs']), list(d['pairs'])
        while todo:
            lam, E = todo.pop()
            for i, (c, dd) in enumerate(E):
                nw = (lam - c, E[:i] + E[i + 1:])
                if nw not in pairs:
                    pairs.add(nw)
                    todo.append(nw)
        ends = [x for x in (Aa, Bb) if x is not None and x != 0]
        for lam, E in pairs:
            if lam != 0:
                lnp |= _primes_of(lam)
            for p in poles:
                if E:
                    Bs.add(('B', lam, E, p, Aa, Bb))
                    s = C.identity_set(lam, E, p, Aa, Bb, fac)
                    if s is not None:
                        for key, m in s[1].values():
                            Bs.add(('B', key[0], key[1], key[2], Aa, Bb))
                elif lam != 0:                       # depth-0 base: E1(-lam (x - p)) at both finite ends
                    for x in (v for v in (Aa, Bb) if v is not None):
                        if x != p:
                            E1s.add(-lam * (x - p))
            for x in ends:
                for p in poles:
                    if x != p:
                        lnp |= _primes_of(x - p)
        for c, dd in fac:
            for x in ends:
                E1s.add(c * x + dd)
            if Aa == 0:
                if dd != 0:
                    E1s.add(dd)
                else:
                    lnp |= _primes_of(c)
    for b in psi_rates:
        E1s.add(-b * D1)
        lnp |= _primes_of(b)
    return dict(B=Bs, E1=E1s, lnp=lnp)


def outside(keys, alph):
    bad = set()
    for _, atoms in keys:
        for a in atoms:
            if a[0] == 'B' and a in alph['B']:
                continue
            if a[0] == 'E1' and a[1] in alph['E1']:
                continue
            if a[0] == 'ln' and int(a[1]) in alph['lnp']:
                continue
            if a[0] == 'g':
                continue
            bad.add(a)
    return bad


# ---------------------------------------------------------------------------------------------- coefficients
def accumulate(stream, Lmax):
    """(L, canonical dict) stream -> padded q-sequences.  Zeros are derived from the PADDED sequences at the end, so a key
    first appearing at L0 > 2 is recorded as zero on 2..L0-1 (the incremental version missed that)."""
    seqs, no_d2 = {}, []
    for L, c in stream:
        for k in c:
            if k not in seqs:
                seqs[k] = [0] * L
        for k, s in seqs.items():
            s.append(int(c[k]) if k in c else 0)
        if L >= 2 and not C.depth2_keys(c):
            no_d2.append(L)
    zeros = {k: [L for L in range(2, len(s)) if s[L] == 0] for k, s in seqs.items()}
    zeros = {k: v for k, v in zeros.items() if v}
    first = {k: next((L for L, v in enumerate(s) if v), None) for k, s in seqs.items()}
    return seqs, zeros, no_d2, first


def coeffs(case, p, Lmax, log=print):
    """Canonical q-sequences mod p; per-key zero lists (L >= 2); per-L depth-2 presence; timing."""
    K = NE.Field(p)
    Fn, stable = C.natural_factors(case)
    t0 = time.time()
    tc = [0.0]
    marks = {}

    def stream():
        for L, q in FX.Engine(case, K).run(Lmax):
            t = time.time()
            c = C.canonical(q, Fn, K)
            tc[0] += time.time() - t
            if L in (100, 300, 600, 1000, 1500) or L == Lmax:
                marks[L] = round(time.time() - t0, 1)
                log(case, p, 'L', L, 'cum s', marks[L], 'canonical s', round(tc[0], 1))
            yield L, c
    seqs, zeros, no_d2, first = accumulate(stream(), Lmax)
    return dict(seqs=seqs, zeros=zeros, no_depth2_L=no_d2, first_nonzero=first, natural_stable=stable,
                t_cum=marks, t_canonical=round(tc[0], 1))


def _runs(Ls):
    """sorted ints -> list of maximal runs [(a, b)]"""
    out = []
    for L in Ls:
        if out and L == out[-1][1] + 1:
            out[-1] = (out[-1][0], L)
        else:
            out.append((L, L))
    return out


def fx1_classify(per_prime, Lmax):
    """per_prime: [(zeros {key: [L >= 2]}, set of keys of that run)].  A key is present at L iff nonzero on at least
    one prime.  Keys in scope: nonzero at some L >= 2 (on some prime).  For each, the L in 2..Lmax absent on both
    primes, split into an INITIAL run (starting at 2), a TERMINAL run (ending at Lmax) and INTERIOR runs.
    FX1's L-independence fails iff some key in scope has an initial or terminal run; interior zeros are reported."""
    full = set(range(2, Lmax + 1))
    keys = set().union(*(ks for _, ks in per_prime))
    absent_by = {}
    for k in keys:
        a = None
        for z, ks in per_prime:
            zk = set(z.get(k, [])) if k in ks else full
            a = zk if a is None else a & zk
        absent_by[k] = sorted(a)
    scope = [k for k in keys if len(absent_by[k]) < len(full)]
    only01 = sorted(map(str, (k for k in keys if k not in scope)))
    initial, terminal, interior = {}, {}, {}
    for k in scope:
        for a, b in _runs(absent_by[k]):
            if a == 2:
                initial[k] = (a, b)
            elif b == Lmax:
                terminal[k] = (a, b)
            else:
                interior.setdefault(k, []).append((a, b))
    return dict(L_independent=not initial and not terminal, n_scope=len(scope), n_only_L01=len(only01),
                initial={str(k): v for k, v in initial.items()}, terminal={str(k): v for k, v in terminal.items()},
                n_keys_interior=len(interior), n_interior_zero_L=sum(b - a + 1 for v in interior.values() for a, b in v),
                max_interior_run=max((b - a + 1 for v in interior.values() for a, b in v), default=0),
                interior_sample={str(k): v[:5] for k, v in list(interior.items())[:10]})


def fx2_classify(per_prime_nod2, Lmax):
    """L >= 2 with no nonzero depth-2 key on either prime: uncapped count, and whether depth 2 is absent at ALL L."""
    nd = sorted(set.intersection(*(set(x) for x in per_prime_nod2)))
    return dict(holds=not nd, n_L_without_depth2=len(nd), absent_at_all_L=len(nd) == Lmax - 1,
                runs=_runs(nd)[:20], n_runs=len(_runs(nd)))


# ---------------------------------------------------------------------------------------------- operator search
def generic(seqs, p, seed):
    keys = sorted(seqs, key=repr)
    rnd = random.Random(seed)
    w = {k: rnd.randrange(1, p) for k in keys}
    n = min(len(seqs[k]) for k in keys)
    return [sum(w[k] * seqs[k][L] for k in keys) % p for L in range(n)]


def _verify(Bf, keys, vec, p):
    cols = [Bf.col(i, j) for i, j in keys]
    nrow = min(len(c) for c in cols)
    zero = flint.nmod(0, p)
    for l in range(nrow):
        s = zero
        for a, c in zip(vec, cols):
            s += a * c[l]
        if int(s):
            return False, nrow
    return True, nrow


def minimal_windowed(seq, p, Lg, Rmax=RMAX, margin=MARGIN, extra=EXTRA):
    """NA's bisection on the window L <= Lg; annihilators at R0..R0+extra verified on every row of the full sequence;
    minimal operator = gcrd (intclass_lib.minimal_operator's construction, with a guess window and verification)."""
    nm = [flint.nmod(int(v) % p, p) for v in seq]
    Bg, Bf = G.Basis(nm[:Lg + 1], p), G.Basis(nm, p)
    trade, R0 = {}, None
    for R in range(1, Rmax + 1):
        Q0, info = G.min_degree_bisect(Bg, R, margin, Lg + 1)
        trade[R] = (Q0, info)
        if Q0 is not None:
            R0 = R
            break
    if R0 is None:
        return dict(found=False, trade=trade)
    ops, found, ver = [], [], []
    for R in range(R0, R0 + extra + 1):
        Q0, nul = G.min_degree_bisect(Bg, R, margin, Lg + 1)
        if Q0 is None:
            continue
        nv, keys, (N, nrow_g) = G.nullity(Bg, R, Q0, margin)
        found.append((R, Q0, nv))
        for col in range(nv):
            vec = [N[i, col] for i in range(len(keys))]
            ok, nrow_f = _verify(Bf, keys, vec, p)
            ver.append(dict(R=R, ok=ok, rows_guess=nrow_g, rows_checked=nrow_f, held_out=nrow_f - nrow_g))
            ops.append(A.standard_form(keys, vec, p))
    g = ops[0]
    for o in ops[1:]:
        g = I.gcrd(g, o, p)
    g = I.content_reduce(g)
    return dict(found=True, search=(R0, trade[R0][0], trade[R0][1]), annihilators=found, verify=ver,
                verified=all(v['ok'] for v in ver), min_held_out=min(v['held_out'] for v in ver),
                c=g, order=len(g) - 1, degree=max(x.degree() for x in g), trade=trade)


def fkey_lin(x0, p):
    return ((-I.to_fp(x0, p)) % p, 1)


def analyse(c, p, marks):
    """points() of the operator; at each marked rational point: present?, class, exponents, local() candidate roots."""
    pts = S.points(c, p)
    at = {}
    for name, x0 in marks.items():
        k = fkey_lin(x0, p)
        loc = A.local(c, flint.nmod(I.to_fp(x0, p), p), p,
                      {n: flint.nmod(I.to_fp(v, p), p) for n, v in CAND.items()})
        at[name] = dict(x=str(x0), present=k in pts, cls=pts[k]['cls'] if k in pts else None,
                        exps=pts[k]['exps'] if k in pts else None, regular=loc['regular'],
                        roots=sorted([n for n, v in loc.get('candidate_roots', {}).items() if v], key=Fr))
    return pts, at


def op_block(seqs, p, Lg, marks, seeds=SEEDS, Rmax=RMAX):
    out, pts_by = [], []
    for sd in seeds:
        t = time.time()
        r = minimal_windowed(generic(seqs, p, sd), p, Lg, Rmax)
        rec = dict(seed=sd, found=r['found'], trade={k: list(v) for k, v in r['trade'].items()})
        if r['found']:
            pts, at = analyse(r['c'], p, marks)
            pts_by.append(pts)
            rec.update(search=r['search'], annihilators=r['annihilators'], verified=r['verified'],
                       min_held_out=r['min_held_out'], order=r['order'], degree=r['degree'], marks=at)
        rec['seconds'] = round(time.time() - t, 1)
        out.append(rec)
    wf = S.weight_fixed(*pts_by) if len(pts_by) == 2 else {}
    gen = S.genuine(wf)
    mk = {fkey_lin(x, p) for x in marks.values()} | {fkey_lin(Fr(1), p), fkey_lin(Fr(-1), p)}
    extra = [dict(deg=v['deg'], cls=v['cls'], pos=v['pos'], root=(None if v['deg'] != 1 else
                  (-k[0]) % p), exps=v['exps']) for k, v in gen.items() if k not in mk]
    return dict(ops=out, weight_fixed=S.describe(wf), genuine_extra=extra,
                n_weight_moving=len(set(pts_by[0]) - set(pts_by[1])) if len(pts_by) == 2 else None)


# ---------------------------------------------------------------------------------------------- scored run
def score(case, p, Lmax=LMAX, hold=HOLD, rmax=RMAX, tag='fx_run'):
    """Scored run (tag 'fx_run', M1/M2r) or a control run under another tag.  Any exception is recorded in the run
    file (instrument_error), which summary() turns into VOID."""
    T0 = time.time()
    lg = open(DATA / f'{tag}_{case}_{p}.log', 'w')

    def log(*a):
        s = ' '.join(map(str, a))
        print(s, flush=True)
        lg.write(s + '\n')
        lg.flush()
    res = dict(case=case, prime=p, LMAX=Lmax, HOLD=hold, RMAX=rmax, instrument_error=None)
    try:
        cf = coeffs(case, p, Lmax, log)
        seqs = cf['seqs']
        (DATA / f'{tag.replace("run", "seq")}_{case}_{p}.json').write_text(json.dumps({repr(k): v for k, v in seqs.items()}))
        alph = alphabet(case)
        bad = outside(seqs, alph)
        lm = landau(case)
        res.update(landau={k: str(v) for k, v in lm.items()},
                   n_keys=len(seqs), zeros={repr(k): v for k, v in cf['zeros'].items()},
                   first_nonzero={repr(k): v for k, v in cf['first_nonzero'].items()},
                   no_depth2_L=cf['no_depth2_L'], natural_stable=cf['natural_stable'],
                   outside_alphabet=sorted(map(repr, bad)), t_coeffs=cf['t_cum'], t_canonical=cf['t_canonical'])
        log('coefficients done', round(time.time() - T0), 's;', len(seqs), 'keys')
        t = time.time()
        res['op'] = op_block(seqs, p, Lmax - hold, dict(xZ=lm['xZ'], xW=lm['xW']), Rmax=rmax)
        res['t_op'] = round(time.time() - t, 1)
    except Exception:                                       # noqa: BLE001  (recorded -> VOID)
        import traceback
        res['instrument_error'] = traceback.format_exc()
        log(res['instrument_error'])
    res['t_total'] = round(time.time() - T0, 1)
    (DATA / f'{tag}_{case}_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    log('done', res['t_total'], 's')


def _crt_pos(roots):
    """{prime: [root mod p]} -> rationals reconstructed from pairs (both primes), height <= 1e6."""
    (p1, r1), (p2, r2) = sorted(roots.items())
    M = p1 * p2
    out = set()
    for a in r1:
        for b in r2:
            x = (a * p2 * pow(p2, -1, p1) + b * p1 * pow(p1, -1, p2)) % M
            q = I.ratrec(x, M, bound=10 ** 6)        # height <= 1e6: spurious pairs pass w.p. ~1e-6
            if q is not None and I.to_fp(Fr(q), p1) == a and I.to_fp(Fr(q), p2) == b:
                out.add(Fr(q))
    return out


PRECEDENCE = ('INCOMPLETE', 'VOID', 'STOP', 'BORDERLINE', 'GO')     # overall = first of these held by any case


def _j(name):
    fn = DATA / name
    return json.loads(fn.read_text()) if fn.exists() else None


def fx0_check():
    """FX0: every instrument control's pass field, read from its frozen output file.  Returns (ok, {check: bool})."""
    c = {}
    g = _j('fx_gate_controls.json') or {}
    v1 = g.get('V1_exact_vs_fc_numeric', {})
    for case in ('M1', 'M2r'):
        for lev in ('A', 'B'):
            c[f'V1_{case}_{lev}_pass_1e_25'] = bool(v1.get(case, {}).get(lev, {}).get('pass_1e_25'))
    c['V2_all_equal_L_le_40'] = bool(g.get('V2_a_eq_b_vs_aabc_exact_A1', {}).get('all_equal_L_le_40'))
    v3 = g.get('V3_modp_vs_fmpq', {})
    for case in ('M1', 'M2r', 'K4A1'):
        c[f'V3_{case}'] = bool(v3.get(case)) and all(r.get('all_equal_L_le_40') for r in v3[case].values()) \
            and len(v3[case]) == 2
    v4 = g.get('V4_NA_8_31', {})
    c['V4_reproduced'] = len(v4) == 2 and all(r.get('reproduced') for r in v4.values())
    d = _j('fx_indep_controls.json') or {}
    for case in ('C1', 'C2'):
        c[f'A1_{case}'] = bool(d.get('A1_identity_numeric', {}).get(case, {}).get('pass_1e_50'))
        a2 = d.get('A2_planted_route_a', {}).get(case, {})
        c[f'A2_{case}'] = bool(a2.get('pos_reduced') and a2.get('neg_survives') and a2.get('modp_equals_fmpq_mod_p'))
        c[f'A4_{case}'] = bool(d.get('A4_modp_canonical', {}).get(case)) and all(d['A4_modp_canonical'][case].values())
    for case in ('C1', 'C2', 'M1', 'M2r'):
        a3 = d.get('A3_engine_value_preserved', {}).get(case, {})
        c[f'A3_{case}'] = bool(a3.get('value_preserved') and a3.get('rewrite_exercised'))
    a5 = d.get('A5_a_eq_b_noop', {})
    c['A5'] = bool(a5.get('no_named_base') and a5.get('canonical_is_identity'))
    bb = d.get('B_route_b_lll', {})
    c['B1_equals_identity'] = bool(bb.get('B1_positive', {}).get('equals_identity'))
    c['B2_involves_planted'] = bool(bb.get('B2_planted', {}).get('involves_planted'))
    c['B3_no_relation'] = bb.get('B3_negative', {}).get('relation_found') is False and \
        bb.get('B3_negative', {}).get('relation_found_dps_plus_100') is False
    b4 = d.get('B4_route_b_whole_set', {})
    c['B4neg_no_relation'] = b4.get('B4neg', {}).get('relation_found') is False and \
        b4.get('B4neg', {}).get('relation_found_dps_plus_100') is False
    c['B4pos_involves_eliminated'] = bool(b4.get('B4pos', {}).get('involves_eliminated'))
    al = _j('fx_runctl_ctrl_alphabet.json') or {}
    for case in ('C1', 'C2'):
        r = al.get(case, {})
        c[f'alphabet_{case}'] = r.get('n_outside') == 0 and bool(r.get('L_independent'))
    for p in PRIMES:
        na = _j(f'fx_runctl_posctrl_na_{p}.json') or {}
        c[f'posctrl_na_{p}'] = bool(na.get('reproduces_frozen') and na.get('xstar_has_3_2')
                                   and na.get('verified') and all(na['verified']))
        ab = _j(f'fx_runctl_posctrl_ab_{p}.json') or {}
        c[f'posctrl_ab_{p}'] = bool(ab.get('reproduces_frozen') and ab.get('xstar_has_7_2'))
        k4 = _j(f'fx_runctl_posctrl_k4a1_{p}.json') or {}
        c[f'posctrl_k4a1_{p}'] = bool(k4.get('minimal_equals_AB_A1') and k4.get('verified') and all(k4['verified'])
                                     and all('7/2' in x['roots'] for x in k4.get('xstar', [])) and k4.get('xstar'))
    e2e = _j('fx_runctl_e2e_k4a1.json') or {}
    c['e2e_k4a1_as_expected'] = bool(e2e.get('as_expected'))
    return all(c.values()), c


def case_summary(case, r1, r2, Lmax=LMAX, hold=HOLD):
    lm = landau(case)
    void = []
    for r in (r1, r2):
        if r.get('instrument_error'):
            void.append(f"instrument exception on prime {r['prime']}")
    if void:
        return dict(void_flags=void, gate='VOID')
    fx1 = fx1_classify([(r['zeros'], set(r['first_nonzero'])) for r in (r1, r2)], Lmax)
    fx1.update(inside_alphabet=not r1['outside_alphabet'] and not r2['outside_alphabet'],
               n_outside=len(set(r1['outside_alphabet']) | set(r2['outside_alphabet'])))
    fx1['holds'] = fx1['L_independent'] and fx1['inside_alphabet']
    fx2 = fx2_classify([r1['no_depth2_L'], r2['no_depth2_L']], Lmax)
    ops = [o for r in (r1, r2) for o in r['op']['ops']]
    shapes = {(o.get('order'), o.get('degree')) for o in ops}
    fx3 = dict(found_all=all(o['found'] for o in ops), verified_all=all(o.get('verified') for o in ops),
               min_held_out=min((o.get('min_held_out') or 0) for o in ops),
               minimal_shapes=sorted(map(str, shapes)), search=sorted({str(o.get('search')) for o in ops}))
    fx3['holds'] = fx3['found_all'] and fx3['verified_all'] and fx3['min_held_out'] >= hold and len(shapes) == 1
    fx4 = dict(per_run=[o.get('marks') for o in ops])
    fx4['holds'] = fx3['found_all'] and all(
        o['marks'][m]['present'] and '7/2' in o['marks'][m]['roots'] for o in ops for m in ('xZ', 'xW'))
    ext = {r['prime']: r['op']['genuine_extra'] for r in (r1, r2)}
    nonlin = any(e['deg'] != 1 for v in ext.values() for e in v)
    pos = _crt_pos({q: [e['root'] for e in v if e['deg'] == 1] for q, v in ext.items()})
    n_lin = max(sum(e['deg'] == 1 for e in v) for v in ext.values())
    xd = lm['xZ'] if lm['dom'] == 'Z' else lm['xW']
    inside = sorted(str(x) for x in pos if 1 < x < xd)
    fx5 = dict(extra=ext, positions=sorted(map(str, pos)), unresolved=nonlin or len(pos) < n_lin, in_1_xdom=inside,
               in_mirror=sorted(str(x) for x in pos if -xd < x < -1))
    # tri-state: False if a resolved genuine point lies in (1, x_dom*); None if none does but some genuine
    # factor is non-linear or not reconstructable at height <= 1e6 (two primes cannot place it); else True
    fx5['holds'] = (None if not fx3['found_all'] else False if inside else None if fx5['unresolved'] else True)
    fx5['strict_holds'] = fx3['found_all'] and not any(ext.values())
    rec = dict(FX1=fx1, FX2=fx2, FX3=fx3, FX4=fx4, FX5=fx5)
    for f, bad in (('natural factors L-dependent', not (r1['natural_stable'] and r2['natural_stable'])),
                   ('an annihilator failed held-out rows', any(o['found'] and not o['verified'] for o in ops))):
        if bad:
            void.append(f)
    rec['void_flags'] = void
    if void:
        gate = 'VOID'
    elif fx3['holds'] and fx1['holds'] and fx2['holds'] and fx4['holds']:
        gate = 'GO'
    elif (not any(o['found'] for o in r1['op']['ops']) and not any(o['found'] for o in r2['op']['ops'])) or \
            any(o['found'] and o['marks'][m]['present'] and '7/2' not in o['marks'][m]['roots']
                for o in ops for m in ('xZ', 'xW')):
        gate = 'STOP'               # no operator on either prime, or a contact point present without 7/2
    else:
        gate = 'BORDERLINE'
    rec['gate'] = gate
    return rec


def summary():
    ok0, checks = fx0_check()
    out = dict(FX0=dict(holds=ok0, failed=[k for k, v in checks.items() if not v], checks=checks))
    for case in SCORED:
        runs = {p: _j(f'fx_run_{case}_{p}.json') for p in PRIMES}
        if any(v is None for v in runs.values()):
            out[case] = dict(gate='INCOMPLETE', primes_present=[p for p, v in runs.items() if v])
            continue
        out[case] = case_summary(case, runs[PRIMES[0]], runs[PRIMES[1]])
    gates = [out[c]['gate'] for c in SCORED]
    if not ok0 and 'INCOMPLETE' not in gates:
        gates.append('VOID')
    out['overall'] = next(g for g in PRECEDENCE if g in gates)
    out['overall_rule'] = ('overall = first of ' + ' > '.join(PRECEDENCE) + ' held by a case; FX0 failure forces VOID; '
                           'GO needs every case GO')
    (DATA / 'fx_run.json').write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: (v.get('gate') if isinstance(v, dict) and 'gate' in v else v) for k, v in out.items()
                      if k != 'FX0'}, default=str), 'FX0', ok0, out['FX0']['failed'])


def e2e_k4a1():
    """End-to-end control: score() + case_summary() on K4A1 (a = b, record AB's class), L <= 1200, window 1000.
    Expected: FX2 False at every L (no depth-2 in AB's class), FX3 minimal (14, 111), FX4 holds, gate not GO."""
    runs = {p: _j(f'fx_runctl_e2e_K4A1_{p}.json') for p in PRIMES}
    rec = case_summary('K4A1', runs[PRIMES[0]], runs[PRIMES[1]], Lmax=1200, hold=HOLD)
    rec['as_expected'] = (rec.get('gate') == 'BORDERLINE' and rec['FX2']['absent_at_all_L'] and rec['FX1']['holds']
                          and rec['FX3']['holds'] and rec['FX3']['minimal_shapes'] == ['(14, 111)'] and rec['FX4']['holds'])
    _ctrl_save('e2e_k4a1', rec)


def firetest_fx1():
    """Planted keys through the real accumulate() + fx1_classify() + fx2_classify():
         A always nonzero; B first appears at L = 10 (late); C sporadic zero at L = 50 on both primes;
         D nonzero only at L in {0, 1}; E zero at L = 30 on prime 1 only; F vanishes from L = 90 on (terminal);
         depth-2 key G absent at L = 40 only (FX2 isolated; also an interior zero of G for FX1)."""
    Lmax = 100
    G = (Fr(0), (('B', Fr(1), ((Fr(1), Fr(0)),), Fr(0), Fr(0), Fr(1)),))
    A_, B_, C_, D_, E_, F_ = ((Fr(i), ()) for i in range(6))

    def stream(prime):
        for L in range(Lmax + 1):
            c = {A_: 1, G: 1}
            if L >= 10:
                c[B_] = 1
            if L != 50:
                c[C_] = 1
            if L <= 1:
                c[D_] = 1
            if not (prime == 0 and L == 30):
                c[E_] = 1
            if L < 90:
                c[F_] = 1
            if L == 40:
                del c[G]
            yield L, c
    per, nod2, per_old = [], [], []
    for prime in (0, 1):
        seqs, zeros, nd, first = accumulate(stream(prime), Lmax)
        per.append((zeros, set(first)))
        nod2.append(nd)
        # the pre-fix incremental rule: a zero recorded only from the key's creation L on
        created = {k: next(L for L, v in enumerate(sq) if v) for k, sq in seqs.items()}
        per_old.append(({k: [L for L in v if L >= created[k]] for k, v in zeros.items()}, set(first)))
    r = fx1_classify(per, Lmax)
    r2 = fx2_classify(nod2, Lmax)
    r_old = fx1_classify(per_old, Lmax)
    res = dict(fx1=r, fx2=r2, old_rule_initial=r_old['initial'])
    res['catches_late_key'] = str(B_) in r['initial'] and not r['L_independent']
    res['catches_terminal_key'] = str(F_) in r['terminal']
    res['sporadic_classified_interior'] = r['n_keys_interior'] == 2 and str(C_) in r['interior_sample'] \
        and str(C_) not in r['initial'] and str(C_) not in r['terminal']
    res['only_L01_excluded'] = r['n_only_L01'] == 1 and str(D_) not in r['initial'] and str(D_) not in r['terminal']
    res['one_prime_zero_not_absent'] = str(E_) not in r['interior_sample']
    res['fx2_isolated'] = r2['n_L_without_depth2'] == 1 and not r2['absent_at_all_L'] and not r2['holds']
    res['old_rule_would_miss_late_key'] = str(B_) not in r_old['initial']
    res['all_pass'] = all(res[k] for k in ('catches_late_key', 'catches_terminal_key', 'sporadic_classified_interior',
                                            'only_L01_excluded', 'one_prime_zero_not_absent', 'fx2_isolated',
                                            'old_rule_would_miss_late_key'))
    _ctrl_save('firetest_fx1', res)


# ---------------------------------------------------------------------------------------------- controls
def _ctrl_save(name, rec):
    (DATA / f'fx_runctl_{name}.json').write_text(json.dumps(rec, indent=1, default=str))
    print(name, json.dumps(rec, default=str)[:2500], flush=True)


def _norm(x):
    return sorted(json.dumps(v) for v in json.loads(json.dumps(x, default=str)))


def posctrl(which, p):
    T0 = time.time()
    frozen = json.loads((DATA / 's1_controls.json').read_text())
    if which == 'na':
        D1, D2 = Fr(7, 5), Fr(11, 5)
        Lm = 1000
        q = NE.coefficients(*S.N3, Lm, NE.Field(p))
        keys = sorted(set(k for d in q for k in d), key=str)
        seqs = {k: [int(d.get(k, flint.nmod(0, p))) for d in q] for k in keys}
        tc = round(time.time() - T0, 1)
        xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
        blk = op_block(seqs, p, Lm - HOLD, dict(xstar=xs))
        fz = frozen[f's1_ctrl_na_{p}']
        ops = blk['ops']
        rec = dict(coeff_s=tc, search=[o.get('search') for o in ops], minimal=[(o.get('order'), o.get('degree')) for o in ops],
                   verified=[o.get('verified') for o in ops], held_out=[o.get('min_held_out') for o in ops],
                   xstar=[o['marks']['xstar'] for o in ops], genuine_extra=blk['genuine_extra'],
                   weight_fixed=blk['weight_fixed'])
        rec['reproduces_frozen'] = (all(tuple(o['search'][:2]) == (8, 31) for o in ops)
                                    and all((o['order'], o['degree']) == (fz['generic'][0]['order'], fz['generic'][0]['degree'])
                                            for o in ops)
                                    and _norm(blk['weight_fixed']) == _norm(fz['weight_fixed']))
        rec['xstar_has_3_2'] = all('3/2' in o['marks']['xstar']['roots'] for o in ops)
        rec['strict_fx5_analog_fires'] = bool(blk['genuine_extra'])
    elif which == 'ab':
        seqs, tc = S.a_coeffs('A1', p)
        import aabc_exact as AB
        D1, D2 = AB.params('A1')[:2]
        xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
        Lm = len(next(iter(seqs.values()))) - 1
        blk = op_block(seqs, p, Lm, dict(xstar=xs))             # AB's setting: the whole L <= 1000 as the window
        fz = frozen[f's1_ctrl_w0_A1_{p}']
        ops = blk['ops']
        rec = dict(search=[o.get('search') for o in ops], minimal=[(o.get('order'), o.get('degree')) for o in ops],
                   xstar=[o['marks']['xstar'] for o in ops], weight_fixed=blk['weight_fixed'])
        rec['reproduces_frozen'] = (all(tuple(o['search'][:2]) == (17, 48) for o in ops)
                                    and all((o['order'], o['degree']) == (14, 111) for o in ops)
                                    and _norm(blk['weight_fixed']) == _norm(fz['minimal']['weight_fixed']))
        rec['xstar_has_7_2'] = all('7/2' in o['marks']['xstar']['roots'] for o in ops)
    elif which == 'k4a1':
        cf = coeffs('K4A1', p, 1200)
        pr = FX.params('K4A1')
        # AB A1: X = B at 0, Y = C at D1, Z = W = A at D2 = a
        xs = (pr['D1'] ** 2 + pr['a'] ** 2) / (2 * pr['D1'] * pr['a'])
        blk = op_block(cf['seqs'], p, 1000, dict(xstar=xs))
        ops = blk['ops']
        rec = dict(coeff_s=cf['t_cum'], search=[o.get('search') for o in ops],
                   minimal=[(o.get('order'), o.get('degree')) for o in ops], verified=[o.get('verified') for o in ops],
                   held_out=[o.get('min_held_out') for o in ops], xstar=[o['marks']['xstar'] for o in ops],
                   weight_fixed=blk['weight_fixed'], n_keys=len(cf['seqs']))
        rec['minimal_equals_AB_A1'] = all((o['order'], o['degree']) == (14, 111) for o in ops)
    else:
        raise SystemExit(which)
    rec['seconds'] = round(time.time() - T0, 1)
    _ctrl_save(f'posctrl_{which}_{p}', rec)


def ctrl_alphabet(Lmax=60):
    """Alphabet rule and L-independence on the unscored control geometries (instrument calibration)."""
    res = {}
    for case in ('C1', 'C2'):
        t = time.time()
        alph = alphabet(case)
        per = {}
        for p in PRIMES:
            cf = coeffs(case, p, Lmax, log=lambda *a: None)
            per[p] = cf
        keys = set(per[PRIMES[0]]['seqs']) | set(per[PRIMES[1]]['seqs'])
        bad = outside(keys, alph)
        f1 = fx1_classify([(per[p]['zeros'], set(per[p]['seqs'])) for p in PRIMES], Lmax)
        f2 = fx2_classify([per[p]['no_depth2_L'] for p in PRIMES], Lmax)
        res[case] = dict(n_keys=len(keys), n_outside=len(bad), outside=sorted(map(repr, bad))[:30],
                         L_independent=f1['L_independent'], fx1=f1, fx2=f2,
                         alphabet_sizes={k: len(v) for k, v in alph.items()},
                         seconds=round(time.time() - t, 1))
    _ctrl_save('ctrl_alphabet', res)


def timing(Lmax=300):
    p = PRIMES[0]
    t = time.time()
    cf = coeffs('C1', p, Lmax, log=lambda *a: None)
    _ctrl_save('timing_C1', dict(prime=p, Lmax=Lmax, t_cum=cf['t_cum'], t_canonical=cf['t_canonical'],
                                 wall=round(time.time() - t, 1), n_keys=len(cf['seqs'])))


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'score':
        if sys.argv[2] in SCORED and '--registered' not in sys.argv:
            raise SystemExit('scored run: register record FX first, then pass --registered')
        score(sys.argv[2], int(sys.argv[3]))
    elif m == 'summary':
        summary()
    elif m == 'posctrl':
        posctrl(sys.argv[2], int(sys.argv[3]))
    elif m == 'ctrl_alphabet':
        ctrl_alphabet()
    elif m == 'timing':
        timing()
    elif m == 'e2e_run':                    # e2e_run <p>: K4A1 control run, L <= 1200, window 1000, R <= 20
        score('K4A1', int(sys.argv[2]), Lmax=1200, hold=HOLD, rmax=20, tag='fx_runctl_e2e')
    elif m == 'e2e_summary':
        e2e_k4a1()
    elif m == 'firetest_fx1':
        firetest_fx1()
