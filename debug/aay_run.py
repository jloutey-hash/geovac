"""Scoring driver for record AY (docs/predictions/2026-11-16_aa_bc_a_side_d_mneq0.md): m != 0 3d pairs on A,
DNX, DNY, DXY, DXX (engine debug/aay_exact.py; DA recomputed as a control).  Written before the record is committed;
the scored modes are run only after the registration commit has landed.

    python debug/aay_run.py toys
        amplitude toys at L = 218..222 (no engine): ta = P_L coefficients of (x* - x)^{7/2}, tb2 = P_L^2 coefficients
        of (1 - x^2)(x* - x)^{7/2}, tb4 = P_L^4 coefficients of (1 - x^2)^2 (x* - x)^{7/2}, tc = P_L coefficients of
        (x* - x)^{9/2} -> data/aay_toys.json   (control, run before registering)
    python debug/aay_run.py p1p3 [Lmax=400] [dps1=3000] [dps2=3600] [case]
        exact coefficients over Q to Lmax, built once, evaluated at both precisions; summed values at the referee
        points (R2 n = 34, data/aay_controls.json); the re-expansion check; beta_eff on coefficient / rho^L for every
        part (f0, f2, f4) and the total c; amplitude ratios; signs; tails for P0 -> data/aay_p1p3[_<case>].json
    python debug/aay_run.py p2 <prime> [Lmax=1000]
        constant sets mod <prime>, per integral and per part (f0, f2 or f4, union) -> data/aay_p2_<prime>.json
No pass/fail bands are coded here: the record sets them.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aay_exact as Y  # noqa: E402
import aabc_exact as X  # noqa: E402
import aabc_run as ABR  # noqa: E402
import na_exact as E  # noqa: E402

DATA = HERE / 'data'
POINTS = {'A1': ('1/2', '-1/3'), 'A2': ('0', '4/5')}
NAMES = ('DNX', 'DNY', 'DXY', 'DXX', 'DA')
PART_M = {'DNX': (0, 2), 'DNY': (0, 2), 'DXY': (0, 4), 'DXX': (0, 4), 'DA': (0,)}


def mpq(q) -> mp.mpf:
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def xstar(case: str):
    D1, D2 = Y.params(case)[:2]
    return mpq(Fraction(D1 ** 2 + D2 ** 2) / (2 * D1 * D2))


# ---------------------------------------------------------------- toys
def toys() -> None:
    mp.mp.dps = 120
    nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(9, mp.mp.prec)
    Lmax = 228

    def anti(cv):
        gv = [mp.mpf(0)] * len(cv)
        for L in range(1, len(cv) - 1):
            gv[L] = cv[L - 1] / (2 * L - 1) - cv[L + 1] / (2 * L + 3)
        return gv
    out = {}
    for case in Y.CASES:
        xs = xstar(case)
        rec = dict(x_star=mp.nstr(xs, 30), x_star_sq_minus_1=mp.nstr(xs * xs - 1, 30),
                   x_star_sq_minus_1_squared=mp.nstr((xs * xs - 1) ** 2, 30))
        for al_s in ('7/2', '9/2'):
            al = mpq(Fraction(al_s))
            cv = [mp.mpf(0)] * (Lmax + 2)
            for x, w in nodes:
                hv = w * (xs - x) ** al
                p0, p1 = mp.mpf(1), x
                cv[0] += hv
                cv[1] += hv * x
                for n in range(1, Lmax + 1):
                    p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
                    cv[n + 1] += hv * p1
            cv = [(2 * n + 1) / mp.mpf(2) * v for n, v in enumerate(cv)]
            b2 = anti(anti(cv))
            b4 = anti(anti(b2))
            tag = {'7/2': '', '9/2': '_92'}[al_s]
            rec['ta' + tag] = {L: mp.nstr(cv[L], 40) for L in range(218, 223)}
            rec['tb2' + tag] = {L: mp.nstr(b2[L], 40) for L in range(218, 223)}
            rec['tb4' + tag] = {L: mp.nstr(b4[L], 40) for L in range(218, 223)}
        # cross-check against record AS's amp_toys (same geometry)
        amp = json.loads((DATA / 'aaa_controls.json').read_text())['amp_toys'][case]
        rec['vs_AS_amp_toys'] = {k: mp.nstr(max(abs(mp.mpf(rec[k2][L]) / mp.mpf(amp[k][str(L)]) - 1)
                                                for L in range(218, 223)), 3)
                                 for k, k2 in (('ta', 'ta'), ('tb', 'tb2'), ('tc', 'ta_92'))}
        out[case] = rec
        print(case, rec['vs_AS_amp_toys'], flush=True)
    out['definition'] = ('P_L^m = (1 - x^2)^{m/2} d^m P_L/dx^m, no CS phase; P_L^m coefficients of (1 - x^2)^{m/2} h = '
                         'Legendre coefficients of h integrated m times; 768-node Gauss-Legendre, 120 digits')
    (DATA / 'aay_toys.json').write_text(json.dumps(out, indent=1))


# ---------------------------------------------------------------- associated Legendre, no CS phase
def plm(Lmax: int, m: int, x) -> list:
    """P_L^m(x) = (1 - x^2)^{m/2} d^m P_L/dx^m, L = 0..Lmax (zero for L < m)."""
    out = [mp.mpf(0)] * (Lmax + 1)
    if m > Lmax:
        return out
    pmm = mp.mpf(1)
    for k in range(1, m + 1):
        pmm *= (2 * k - 1)
    pmm *= (1 - x * x) ** (mp.mpf(m) / 2)
    out[m] = pmm
    if m + 1 <= Lmax:
        out[m + 1] = x * (2 * m + 1) * pmm
    for L in range(m + 1, Lmax):
        out[L + 1] = ((2 * L + 1) * x * out[L] - (L + m) * out[L - 1]) / (L - m + 1)
    return out


def referee() -> dict:
    d = json.loads((DATA / 'aay_controls.json').read_text())['referee']
    out = {}
    for k, v in d.items():
        if '_x=' not in k:
            continue
        case, integral, _ = k.split('_', 2)
        out[(case, integral, v['x'])] = dict(R2_34=v['R2']['34'], R2_30=v['R2']['30'], R1=v['R1']['120,120,144'])
    return out


# ---------------------------------------------------------------- P1, P3, P4, P0
def p1p3(Lmax: int = 400, dps1: int = 3000, dps2: int = 3600, cases=None) -> None:
    K = E.Field(None)
    ref = referee()
    amp = json.loads((DATA / 'aay_toys.json').read_text())
    dpss = (dps1, dps2)
    out = {}
    outf = DATA / ('aay_p1p3.json' if not cases else f'aay_p1p3_{"_".join(cases)}.json')
    for case in (cases or Y.CASES):
        t0 = time.time()
        eng = Y.Engine(case, K, Lmax + 2)
        D2 = eng.D2
        tw = {}
        cval = {d: {} for d in dpss}
        V = {d: {nm: {m: [mp.mpf(0)] * (Lmax + 1) for m in PART_M[nm]} for nm in NAMES} for d in dpss}
        for L in range(Lmax + 1):
            Gs = {kern: eng.G_kernel(kern, L, Y.NMAX[kern]) for kern in ('U4', 'U2', 'U0')}
            for nm in NAMES:
                for m, rec in Y.channel_jets(nm, L).items():
                    acc = {}
                    for kern, jet in rec.items():
                        for (j, p), cc in jet.items():
                            a = K(cc) * K.pw(D2, -p)
                            for k, v in Gs[kern][j].items():
                                E.add(acc, k, a * v)
                    q = {k: v for k, v in acc.items() if v != 0}
                    for d in dpss:
                        with mp.workdps(d):
                            if d not in tw:
                                tw[d] = Y.twopiN(case)
                            cs = cval[d]
                            tot = mp.mpf(0)
                            for k, v in q.items():
                                if k not in cs:
                                    cs[k] = E.value(k, mp)
                                tot += mp.mpf(int(v.p)) / int(v.q) * cs[k]
                            V[d][nm][m][L] = tw[d] * tot
            for k in [k for k in eng.cache if k[2] == L]:
                del eng.cache[k]
            eng.dens_cache.pop(L, None)
            if L % 50 == 0:
                print(case, 'L', L, round(time.time() - t0), 's', flush=True)
        eng_s = round(time.time() - t0)
        res = {}
        for d in dpss:
            mp.mp.dps = d
            xs = xstar(case)
            rho = 1 / (xs + mp.sqrt(xs * xs - 1))
            A_ = amp[case]
            for nm in NAMES:
                f = V[d][nm]
                rec = dict(dps=d, Lmax=Lmax, engine_s=eng_s, values={}, diffs={}, beta={}, signs={}, tails={})
                c = Y.total_coeffs(f) if nm != 'DA' else f[0]
                for xs_ in POINTS[case]:
                    xr = mpq(Fraction(xs_))
                    val = mp.mpf(0)
                    tail = mp.mpf(0)
                    for m, fm in f.items():
                        P = plm(Lmax, m, xr)
                        val += sum(fm[L] * P[L] for L in range(Lmax + 1))
                        tail = max(tail, abs(fm[Lmax] * P[Lmax]))
                    rec['values'][xs_] = mp.nstr(val, 30)
                    if nm != 'DA':
                        P0 = plm(len(c) - 1, 0, xr)
                        rec.setdefault('total_reexpansion_check', {})[xs_] = mp.nstr(
                            abs(sum(c[j] * P0[j] for j in range(len(c))) - val), 3)
                    r = ref.get((case, nm, xs_))
                    rec['diffs'][xs_] = dict(tail=mp.nstr(tail, 5))
                    if r:
                        rec['diffs'][xs_].update(R2_34=float(val - mp.mpf(r['R2_34'])),
                                                 R2_30=float(val - mp.mpf(r['R2_30'])),
                                                 R1=float(val - mp.mpf(r['R1'])))
                series = {f'f{m}': fm for m, fm in f.items()}
                if nm != 'DA':
                    series['c'] = c[:Lmax + 1]
                for n, F in series.items():
                    G = [F[L] / rho ** L for L in range(Lmax + 1)]
                    lo = {'f0': 1, 'f2': 3, 'f4': 5, 'c': 1}[n]
                    top = Lmax if n != 'c' else Lmax - 160      # c_j needs the f tail beyond j
                    rec['beta'][n] = {L: float(mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1)))
                                      for L in sorted(set(list(range(20, top, 20)) + [220] + [Lmax - 20]))
                                      if L - 1 >= lo and L + 1 <= top and G[L - 1] != 0}
                    rec['signs'][n] = dict(L_ge_20=sorted(set(int(mp.sign(v)) for v in F[20:top + 1])),
                                           L_ge_5=sorted(set(int(mp.sign(v)) for v in F[5:top + 1])))
                    rec['tails'][n] = {L: mp.nstr(F[L], 25) for L in (219, 220, 221, Lmax - 2, Lmax - 1, Lmax)
                                       if 0 <= L <= (Lmax if n != 'c' else min(221, top))}
                # amplitude ratios: (f0 / ta) / (f_m / tb_m), m = 2 or 4
                if nm != 'DA' and Lmax >= 400:
                    m = PART_M[nm][1]
                    tb = A_['tb2'] if m == 2 else A_['tb4']
                    rec['amp_ratio'] = {L: float((f[0][L] / mp.mpf(A_['ta'][str(L)])) / (f[m][L] / mp.mpf(tb[str(L)])))
                                        for L in (219, 220, 221)}
                    rec['c_over_ta'] = {L: float(c[L] / mp.mpf(A_['ta'][str(L)])) for L in (219, 220, 221)}
                    ct = Y.total_coeffs({mm: list(fm[:381]) + [mp.mpf(0)] * (Lmax - 380) for mm, fm in f.items()})
                    rec['c_truncation_380_vs_400'] = {L: mp.nstr(abs(ct[L] / c[L] - 1), 3) for L in (219, 220, 221)}
                    rec['c_over_tc'] = {L: float(c[L] / mp.mpf(A_['ta_92'][str(L)])) for L in (219, 220, 221)}
                res[f'{case}_{nm}_{d}'] = rec
                print(case, nm, d, json.dumps({k: rec[k] for k in ('values', 'diffs', 'beta')}, indent=0)[:1500],
                      flush=True)
        # P0: self-differences between the two precisions (tails and values)
        p0 = {}
        with mp.workdps(max(dpss)):
            for nm in NAMES:
                a, b = res[f'{case}_{nm}_{dps1}'], res[f'{case}_{nm}_{dps2}']
                worst = mp.mpf(0)
                for n in a['tails']:
                    for L, v in a['tails'][n].items():
                        u, w = mp.mpf(v), mp.mpf(b['tails'][n][L])
                        worst = max(worst, mp.mpf(0) if u == w else (abs(w / u - 1) if u != 0 else mp.inf))
                for xs_ in a['values']:
                    u, w = mp.mpf(a['values'][xs_]), mp.mpf(b['values'][xs_])
                    worst = max(worst, abs(w / u - 1))
                p0[nm] = mp.nstr(worst, 3)
        out.update(res)
        out[f'{case}_P0_worst_relative'] = p0
        print(case, 'P0', p0, flush=True)
        outf.write_text(json.dumps(out, indent=1, default=str))


# ---------------------------------------------------------------- P2
def p2(p: int, Lmax: int = 1000) -> None:
    out = {}
    K = E.Field(p)
    for case in Y.CASES:
        na_keys = set(k for q in X.pieces(case, 40, K, which=('NA',))['NA'] for k in q)
        C = ABR.candidate_set(case, na_keys)
        ab21 = set(json.loads((DATA / f'aabc_p2_{p}.json').read_text())[case]['constants'])
        eng = Y.Engine(case, K, Lmax + 2)
        for nm in NAMES[:4]:
            t0 = time.time()
            ch = Y.coefficients(case, Lmax, K, nm, eng=eng)
            sec = round(time.time() - t0)
            m2 = PART_M[nm][1]
            series = {'f0': ch[0], f'f{m2}': ch[m2]}
            series['union'] = [set(x) | set(y) for x, y in zip(ch[0], ch[m2])]
            recs = {}
            for name, f in series.items():
                sets = [frozenset(q) for q in f]
                L0 = 2 if name != 'f4' else 4
                S0 = sets[L0]
                changes = [L for L in range(L0 + 1, len(sets)) if sets[L] != S0]
                union = set().union(*sets[L0:])
                first_fixed = next((L for L in range(len(sets)) if all(sets[M] == S0 for M in range(L, len(sets)))),
                                   None)
                u = set(map(repr, union))
                recs[name] = dict(L0=L0, n_constants_L0=len(S0), L_where_set_changes=changes[:20], fixed=not changes,
                                  first_L_fixed=first_fixed, sizes_L0_to_5=[len(sets[L]) for L in range(6)],
                                  n_union=len(union), inside_C=not (union - C),
                                  outside_C=sorted(map(repr, union - C)),
                                  gamma_present=[repr(k) for k in union if k[1] == 'g'],
                                  equals_AB21=u == ab21, vs_AB21=dict(missing=sorted(ab21 - u), extra=sorted(u - ab21)))
            out[f'{case}_{nm}'] = dict(seconds=sec, Lmax=Lmax, series=recs)
            print(case, nm, p, sec, 's', json.dumps(recs, indent=1), flush=True)
            (DATA / f'aay_p2_{p}.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    cmd = sys.argv[1]
    nums = [int(a) for a in sys.argv[2:] if a.isdigit()]
    cs = [a for a in sys.argv[2:] if a in Y.CASES]
    if cmd == 'toys':
        toys()
    elif cmd == 'p1p3':           # p1p3 [Lmax dps1 dps2] [A1|A2]  (one case -> data/aay_p1p3_<case>.json)
        p1p3(*nums, cases=cs or None)
    elif cmd == 'p2':
        p2(nums[0], *(nums[1:2] or [1000]))
