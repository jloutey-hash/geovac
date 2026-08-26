"""Guarded, decoy-calibrated, weight-graded PSLQ on the high-precision T2 value.

Protocol (corpus standard: sprint_routeC_corner_sigma2_pslq_memo.md + routeC_T2_pslq_decisive.py):
  * every constant built as mpf (never int**-int);
  * a same-magnitude STRUCTURELESS DECOY tested against every basis at every precision;
  * >= 2 working precisions; a hit must be stable across them;
  * maxcoeff CALIBRATED per basis: 10^(dps/n) is the height at which a random real acquires a
    spurious relation against an n-element ring (PSLQ vector length n+1), so we search only to 10^(dps/n - margin).  A
    "no relation" is then a genuine HEIGHT-BOUNDED negative; if the budget collapses below 10 the
    leg is reported UNDERPOWERED (the precision cannot resolve a ring that large) rather than NEG.
Usage: python debug/beta2_t2_pslq.py <T2 decimal string> <ndig_trusted>
"""
from __future__ import annotations
import sys, itertools
import mpmath as mp
sys.path.insert(0, 'debug')
from routeC_T2_pslq_decisive import guarded_pslq


def consts(dps):
    mp.mp.dps = dps + 40
    pi = +mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))       # K(1/2)
    P8 = (mp.sqrt(1 + mp.sqrt(2)) * mp.gamma(mp.mpf(1) / 8) * mp.gamma(mp.mpf(3) / 8)
          / (mp.mpf(2) ** (mp.mpf(13) / 4) * mp.sqrt(pi)))          # disc-8 period
    d = dict(pi=pi, vp=varpi, P8=P8, G=+mp.catalan, ln2=mp.log(2), z3=mp.zeta(3))
    mp.mp.dps = dps
    return {k: +v for k, v in d.items()}


WT = dict(pi=1, vp=1, P8=1, G=2, ln2=1, z3=3)
SIGNED = dict(pi=False, vp=True, P8=True, G=False, ln2=False, z3=False)


def graded(dps, names, wmax, signed_override=None):
    C = consts(dps)
    mp.mp.dps = dps
    ring = {'1': mp.mpf(1)}
    ranges = []
    for n in names:
        sg = SIGNED[n] if signed_override is None else signed_override.get(n, SIGNED[n])
        hi = wmax // WT[n]
        ranges.append(range(-hi, hi + 1) if sg else range(0, hi + 1))
    for exps in itertools.product(*ranges):
        wt = sum(abs(e) * WT[names[i]] for i, e in enumerate(exps))
        if 1 <= wt <= wmax:
            key = '*'.join(f'{names[i]}^{e}' for i, e in enumerate(exps) if e)
            val = mp.mpf(1)
            for i, e in enumerate(exps):
                val *= C[names[i]] ** e
            ring[key] = val
    return ring


def custom_leg(label, target_str, decoy_str, ring_fn, dps_list, margin=2):
    print(chr(10) + "-" * 92); print(f"LEG {label}"); print("-" * 92)
    verdicts = []
    for dps in dps_list:
        ring = ring_fn(dps); n = len(ring)
        mp.mp.dps = dps + 20
        tgt = mp.mpf(target_str); dec = mp.mpf(decoy_str)
        mp.mp.dps = dps
        expo = mp.mpf(dps) / n   # PSLQ vector length is n+1 (target + ring)
        mc_exp = min(int(mp.floor(expo)) - margin, 12)
        maxcoeff = 10 ** mc_exp if mc_exp >= 1 else 0
        print(f"  dps={dps} dim={n} fp-height=10^{mp.nstr(expo,3)}  search maxcoeff=10^{mc_exp}")
        if maxcoeff == 0:
            print("    -> UNDERPOWERED"); verdicts.append('UNDERPOWERED'); continue
        real = guarded_pslq(+tgt, ring, dps, maxcoeff, "REAL")
        deco = guarded_pslq(+dec, ring, dps, maxcoeff, "DECOY")
        hr = None if real is None else max(abs(x) for x in real)
        hd = None if deco is None else max(abs(x) for x in deco)
        v = (f'DECISIVE-NEG(h<=1e{mc_exp})' if real is None
             else ('CANDIDATE' if (hd is None or hd > 8 * hr) else 'UNDERPOWERED(decoy matched)'))
        print(f"    -> {v}   (REAL h={hr}, DECOY h={hd})")
        verdicts.append(v)
    print(f"  ==> {label}: {verdicts}")
    return verdicts


def leg(label, target_str, decoy_str, names, wmax, dps_list, signed_override=None, margin=2):
    print("\n" + "-" * 92)
    print(f"LEG {label}: ring {{{','.join(names)}}} wt<={wmax}")
    print("-" * 92)
    verdicts = []
    for dps in dps_list:
        ring = graded(dps, names, wmax, signed_override)
        n = len(ring)
        mp.mp.dps = dps + 20
        tgt = mp.mpf(target_str); dec = mp.mpf(decoy_str)
        mp.mp.dps = dps
        expo = mp.mpf(dps) / n   # PSLQ vector length is n+1 (target + ring)
        mc_exp = min(int(mp.floor(expo)) - margin, 12)   # cap: a relation of height >1e12 is not
        maxcoeff = 10 ** mc_exp if mc_exp >= 1 else 0    # a 'closed form' in any useful sense
        print(f"  dps={dps} dim={n} fp-height=10^{mp.nstr(expo,3)}  search maxcoeff=10^{mc_exp}")
        if maxcoeff == 0:
            print("    -> UNDERPOWERED (no height budget: ring dim too large for this precision)")
            verdicts.append('UNDERPOWERED')
            continue
        real = guarded_pslq(+tgt, ring, dps, maxcoeff, "REAL")
        deco = guarded_pslq(+dec, ring, dps, maxcoeff, "DECOY")
        hr = None if real is None else max(abs(x) for x in real)
        hd = None if deco is None else max(abs(x) for x in deco)
        if real is None:
            v = f'DECISIVE-NEG(h<=1e{mc_exp})'
        elif hd is None or hd > 8 * hr:
            v = 'CANDIDATE'
        else:
            v = 'UNDERPOWERED(decoy matched)'
        print(f"    -> {v}   (REAL h={hr}, DECOY h={hd})")
        verdicts.append(v)
    agree = 'STABLE' if len(set(verdicts)) == 1 else 'UNSTABLE-ACROSS-PRECISION'
    print(f"  ==> {label}: {verdicts}  [{agree}]")
    return verdicts


def main():
    T2s = sys.argv[1]
    nd = int(sys.argv[2])
    mp.mp.dps = nd + 25
    T2 = mp.mpf(T2s)
    W = T2 * mp.pi / 8
    decT = T2 * (1 + mp.mpf(10) ** (-11)) + mp.euler / mp.mpf(10) ** 5
    decW = W * (1 + mp.mpf(10) ** (-11)) + mp.euler / mp.mpf(10) ** 5
    Ts, Ws = mp.nstr(T2, nd + 15), mp.nstr(W, nd + 15)
    dTs, dWs = mp.nstr(decT, nd + 15), mp.nstr(decW, nd + 15)
    d1, d2, d3 = nd - 8, nd - 4, nd - 1
    print(f"T2 = {Ts}")
    print(f"W = T2*pi/8 = {Ws}")
    print(f"trusted digits nd={nd}; dps grid {d1},{d2},{d3}")
    U = {'vp': False, 'P8': False}
    for name, tgt, dec in (('W', Ws, dWs), ('T2', Ts, dTs)):
        print("\n" + "=" * 92)
        print(f"TARGET = {name}")
        print("=" * 92)
        leg(f'{name}/A disc-4 wt<=1 {{1,pi,vp}}',   tgt, dec, ['pi', 'vp'],       1, [d1, d3], U)
        leg(f'{name}/B {{1,pi,vp,G}}',              tgt, dec, ['pi', 'vp', 'G'],  2, [d1, d3], U)
        leg(f'{name}/C disc-8 wt<=1',               tgt, dec, ['pi', 'vp', 'P8'], 1, [d1, d3], U)
        leg(f'{name}/D disc-4 wt<=2',               tgt, dec, ['pi', 'vp'],       2, [d1, d3], U)
        leg(f'{name}/E disc-4 wt<=3',               tgt, dec, ['pi', 'vp'],       3, [d1, d3], U)
        leg(f'{name}/F disc-8 wt<=2',               tgt, dec, ['pi', 'vp', 'P8'], 2, [d1, d3], U)
        leg(f'{name}/G corrected wt<=2 NO-G',       tgt, dec, ['pi', 'vp'],       2, [d1, d3])
        leg(f'{name}/H corrected wt<=2 +G',         tgt, dec, ['pi', 'vp', 'G'],  2, [d1, d3])
        leg(f'{name}/I corrected wt<=3 NO-G',       tgt, dec, ['pi', 'vp'],       3, [d1, d3])
        leg(f'{name}/J corrected wt<=3 +G  [pre-registered T-2 ring]',
            tgt, dec, ['pi', 'vp', 'G'], 3, [d2, d3])
        leg(f'{name}/K corrected wt<=3 +G +ln2',    tgt, dec, ['pi', 'vp', 'G', 'ln2'], 3, [d2, d3])
        # targeted, LOW-DIMENSION beta(2) probes (maximum height budget)
        def mk(keys):
            def f(dps):
                C = consts(dps); mp.mp.dps = dps
                r = {'1': mp.mpf(1)}
                for k_ in keys:
                    v = mp.mpf(1)
                    for tok in k_.split('*'):
                        b_, e_ = tok.split('^'); v *= C[b_] ** int(e_)
                    r[k_] = v
                return r
            return f
        for keys in (['G^1'], ['pi^1', 'G^1'], ['vp^1', 'G^1'], ['vp^2', 'G^1'],
                     ['pi^1*vp^1', 'G^1'], ['pi^1*vp^-1', 'G^1'], ['pi^2', 'G^1'],
                     ['pi^1', 'vp^1', 'vp^-1', 'G^1']):
            custom_leg(f'{name}/G-probe {{1,' + ','.join(keys) + '}}', tgt, dec, mk(keys), [d1, d3])


if __name__ == '__main__':
    main()
