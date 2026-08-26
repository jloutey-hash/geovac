"""
A-vs-C verdict from the balanced-LiH max_n = 2, 3, 4 series.

Two read-outs, both protocol-identical across max_n:

  (i) 3-point central difference on the symmetric h=0.1 stencil
      {R_true-h, R_true, R_true+h}: tilt and curvature AT R_true.
      This is the ABC memo's registered quantity ("balanced LiH curvature"),
      and it needs only 3 points, so it is available at every max_n.

  (ii) quartic fit over the full 6-point decider grid: the curve's own minimum
      (R_eq, omega_e) -- the chemistry-error memo's quantity.

The verdict test is NOT "did the number move" but "is the remaining movement
budget under the observed geometric decay large enough to reach the true value".
"""
from __future__ import annotations
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from davidson_pes import fit_shape, K_TRUE, R_TRUE, W_E_TRUE

H = 0.1


def series(tag='faithful'):
    out = {}
    for max_n in (2, 3, 4):
        fn = f'debug/data/davidson_pes_n{max_n}_decider.json'
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        pts = {round(r['R'], 4): r[tag]['E'] for r in d['rows'] if tag in r}
        rec = {'n_pts': len(pts), 'pts': pts}
        need = (round(R_TRUE - H, 4), round(R_TRUE, 4), round(R_TRUE + H, 4))
        if all(k in pts for k in need):
            Em, E0, Ep = pts[need[0]], pts[need[1]], pts[need[2]]
            rec['E_Rtrue'] = E0
            rec['tilt3'] = (Ep - Em) / (2 * H)
            rec['curv3'] = (Ep - 2 * E0 + Em) / H ** 2
            rec['curv3_over_k'] = rec['curv3'] / K_TRUE
        if len(pts) >= 5:
            Rs = sorted(pts); Es = [pts[r] for r in Rs]
            f = fit_shape(Rs, Es, order=4)
            rec['fit4'] = f
        out[max_n] = rec
    return out


def model_robustness(name, vals, target, unit=''):
    """The verdict must not rest on one extrapolation model. Fit three laws to the
    same 3 points and report the max_n->inf limit each predicts, plus the
    model-free 'how many more shells at the CURRENT step size' figure."""
    v2, v3, v4 = vals
    d1, d2 = v3 - v2, v4 - v3
    print(f"  extrapolation-model robustness for {name}:")
    # (a) geometric
    r = d2 / d1
    print(f"    geometric (ratio r={r:.3f})        -> limit {v4 + d2*r/(1-r):+.5f}{unit}")
    # (b) power law v = a + b n^-p, p solved from the observed step ratio
    from scipy.optimize import brentq

    def f(p):
        return ((4.0 ** -p - 3.0 ** -p) / (3.0 ** -p - 2.0 ** -p)) - r
    try:
        p = brentq(f, 0.05, 12.0)
        b = d1 / (3.0 ** -p - 2.0 ** -p)
        print(f"    power law  n^-p  (p={p:.2f})       -> limit {v4 - b*4.0**-p:+.5f}{unit}")
    except Exception as e:                                     # pragma: no cover
        print(f"    power law: not solvable ({e})")
    # (c) model-free floor: 1/n_max decay is SLOWER than observed; use it as a
    #     conservative upper bound on the remaining tail
    tail_1overn = d2 * 3.0      # sum_{n>=4} of a c/(n(n+1)) tail, relative to d2
    print(f"    conservative 1/n_max tail bound     -> limit {v4 + tail_1overn:+.5f}{unit}")
    need = target - v4
    if d2 != 0 and not np.isnan(need):
        print(f"    MODEL-FREE: reaching the target needs {need:+.4f}; the current "
              f"step is {d2:+.5f} => {abs(need/d2):.0f} MORE SHELLS at an "
              f"undecaying step size (the step is in fact decaying by {1-r:.0%}/shell"
              + (", and in the WRONG direction" if np.sign(d2) != np.sign(need) else "")
              + ")")


def geometric_report(name, vals, target, unit=''):
    """vals = [v(n=2), v(n=3), v(n=4)]; target = the physically correct value."""
    d1, d2 = vals[1] - vals[0], vals[2] - vals[1]
    r = d2 / d1 if d1 != 0 else float('nan')
    remaining = d2 * r / (1 - r) if abs(r) < 1 else float('nan')
    limit = vals[2] + remaining
    gap_now = vals[2] - target
    print(f"  {name}")
    print(f"    n_max=2,3,4 : {vals[0]:+.6f}{unit}  {vals[1]:+.6f}{unit}  {vals[2]:+.6f}{unit}")
    print(f"    steps       : 2->3 {d1:+.6f}   3->4 {d2:+.6f}   ratio r = {r:+.3f}")
    if abs(r) < 1:
        print(f"    ALL remaining movement under geometric decay = {remaining:+.6f}")
        print(f"    extrapolated max_n->inf limit = {limit:+.6f}{unit}   "
              f"(target {target:+.6f}{unit})")
        need = target - vals[2]
        print(f"    distance still to travel to reach target = {need:+.6f}  "
              f"=> budget covers {abs(remaining/need)*100:5.1f}% of it"
              if need != 0 else "")
        print(f"    residual gap at the extrapolated limit = {limit - target:+.6f}{unit} "
              f"({abs(limit-target)/abs(gap_now)*100:.1f}% of the n_max=4 gap survives)")
    else:
        print("    |r| >= 1: not geometrically convergent -- no extrapolation")


if __name__ == '__main__':
    print("=" * 100)
    print("A-vs-C VERDICT -- balanced LiH, max_n = 2, 3, 4")
    print("=" * 100)
    print("Registered predictions (debug/sprint_abc_connections_test_memo.md +")
    print("debug/aha_track4a_findings.md Probe 2):")
    print("  A (irreducible free-side wall): curvature/omega_e stays FROZEN near")
    print("     the n_max=2,3 value (omega_e ~ +45%, curv/k_true ~ 2.1-2.2x).")
    print("  C (slow basis-response error) : it RELAXES back toward the true")
    print("     stiffness (omega_e -> +0%, curv/k_true -> 1.0x).")
    print(f"\nk_true = {K_TRUE:.6f} Ha/bohr^2 ; omega_e_true = {W_E_TRUE} cm^-1\n")

    for tag in ('faithful', 'corrected'):
        S = series(tag)
        if not all(k in S and 'curv3' in S[k] for k in (2, 3, 4)):
            print(f"[{tag}] incomplete series -- have "
                  f"{[k for k in S if 'curv3' in S[k]]}")
            continue
        print("=" * 100)
        print(f"[{tag}]  " + ("shipped library sign convention"
                              if tag == 'faithful' else "sign-corrected physics"))
        print("=" * 100)
        print(f"  {'max_n':>5s} {'E(R_true)':>17s} {'tilt(R_true)':>13s} "
              f"{'curv(R_true)':>13s} {'curv/k_true':>12s}")
        for n in (2, 3, 4):
            s = S[n]
            print(f"  {n:5d} {s['E_Rtrue']:+17.12f} {s['tilt3']:+13.6f} "
                  f"{s['curv3']:+13.6f} {s['curv3_over_k']:12.4f}")
        print()
        geometric_report("E(R_true)  [Ha]  -- the energy leg, for reference",
                         [S[n]['E_Rtrue'] for n in (2, 3, 4)],
                         target=float('nan'))
        print()
        geometric_report("curv(R_true)/k_true  -- THE REGISTERED DECIDER",
                         [S[n]['curv3_over_k'] for n in (2, 3, 4)], target=1.0, unit='x')
        model_robustness("curv(R_true)/k_true",
                         [S[n]['curv3_over_k'] for n in (2, 3, 4)], target=1.0, unit='x')
        print()
        geometric_report("tilt(R_true)  [Ha/bohr]",
                         [S[n]['tilt3'] for n in (2, 3, 4)], target=0.0)
        model_robustness("tilt(R_true)",
                         [S[n]['tilt3'] for n in (2, 3, 4)], target=0.0)
        have4 = all('fit4' in S[n] and S[n]['fit4']['own_min'] for n in (2, 3, 4))
        if have4:
            print()
            geometric_report("omega_e at own minimum  [cm^-1]",
                             [S[n]['fit4']['own_min']['omega_e_cm1'] for n in (2, 3, 4)],
                             target=W_E_TRUE)
            print()
            geometric_report("R_eq  [bohr]",
                             [S[n]['fit4']['own_min']['R_eq'] for n in (2, 3, 4)],
                             target=R_TRUE)
        else:
            print("\n  [6-point grid not complete at every max_n -- own-minimum "
                  "(omega_e, R_eq) legs deferred]")
            for n in (2, 3, 4):
                print(f"    max_n={n}: {S[n]['n_pts']} points on disk")
        print()
