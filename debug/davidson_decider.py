"""
The banked A-vs-C decider: balanced-LiH curvature / omega_e vs max_n.

Reads the davidson_pes_n{2,3,4}_decider.json sweeps and prints the registered
comparison. Registered predictions, quoted from
debug/sprint_abc_connections_test_memo.md ("The residual A-vs-C question"):

  "The distinguishing question -- irreducible free-side wall (A) vs
   slow-but-eventual basis-response error (C) -- turns on whether the curvature
   converges as max_n -> infinity. Over n_max 2->3 it does not converge (frozen),
   consistent with both. Deciding needs n_max >= 4 ... Named decider (banked,
   not run): balanced LiH curvature at n_max=4."

and the operational reading from debug/aha_track4a_findings.md (Probe 2):

  "if curvature/omega_e at n_max=4 stays frozen near the n_max=2,3 value
   (omega_e ~ +45%, curv/k_true ~ 2.1-2.2x) -> supports A (irreducible wall).
   If it relaxes back toward the true stiffness (omega_e -> +0%, curv/k -> 1.0x)
   -> supports C (slow basis-response error, healing)."
"""
from __future__ import annotations
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from davidson_pes import fit_shape, K_TRUE, R_TRUE, W_E_TRUE

E_EXACT = -8.0706        # LiH BO total, experiment-derived (pk_amplification_lih.json)


def load(max_n, grid='decider'):
    fn = f'debug/data/davidson_pes_n{max_n}_{grid}.json'
    if not os.path.exists(fn):
        return None
    return json.load(open(fn))


if __name__ == '__main__':
    print("=" * 104)
    print("BANKED A-vs-C DECIDER -- balanced LiH curvature / omega_e vs max_n")
    print("=" * 104)
    print(f"k_true = {K_TRUE:.6f} Ha/bohr^2   omega_e_true = {W_E_TRUE} cm^-1   "
          f"R_true = {R_TRUE} bohr\n")

    rows = []
    for max_n in (2, 3, 4):
        d = load(max_n)
        if d is None:
            print(f"  [max_n={max_n}: no decider sweep on disk]")
            continue
        for tag in ('faithful', 'corrected'):
            pts = [(r['R'], r[tag]['E'], r[tag]['converged'], r[tag]['residual'])
                   for r in d['rows'] if tag in r]
            if len(pts) < 4:
                print(f"  [max_n={max_n} {tag}: only {len(pts)} points -- skipped]")
                continue
            Rs = [p[0] for p in pts]; Es = [p[1] for p in pts]
            f3 = fit_shape(Rs, Es, order=3)
            f4 = fit_shape(Rs, Es, order=4) if len(pts) >= 5 else None
            rows.append(dict(max_n=max_n, tag=tag, n_pts=len(pts),
                             all_conv=all(p[2] for p in pts),
                             max_res=max(p[3] for p in pts),
                             f3=f3, f4=f4, Rs=Rs, Es=Es))

    hdr = (f"{'max_n':>5s} {'mode':>10s} {'n':>2s} {'tilt(Rt)':>10s} {'curv(Rt)':>10s} "
           f"{'/k_true':>8s} {'R_eq':>7s} {'R_eq%':>7s} {'omega_e':>8s} {'w_e%':>7s} "
           f"{'curv@min/k':>10s}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        f = r['f3']; om = f['own_min']
        print(f"{r['max_n']:5d} {r['tag']:>10s} {r['n_pts']:2d} "
              f"{f['tilt_at_Rtrue']:+10.6f} {f['curv_at_Rtrue']:+10.6f} "
              f"{f['curv_over_ktrue_at_Rtrue']:8.3f} "
              + (f"{om['R_eq']:7.4f} {om['R_eq_err_pct']:+7.2f} "
                 f"{om['omega_e_cm1']:8.0f} {om['omega_e_err_pct']:+7.1f} "
                 f"{om['curv_over_ktrue']:10.3f}" if om else " " * 44))
    print("\n(cubic fit over the uniform h=0.1 decider grid; "
          "'Rt' = R_true = 3.015 bohr)")

    print("\nquartic-fit cross-check (same points, order=4):")
    for r in rows:
        if r['f4'] is None:
            continue
        f = r['f4']; om = f['own_min']
        print(f"  max_n={r['max_n']} {r['tag']:>10s}: tilt={f['tilt_at_Rtrue']:+.6f} "
              f"curv={f['curv_at_Rtrue']:+.6f} ({f['curv_over_ktrue_at_Rtrue']:.3f}x)"
              + (f"  R_eq={om['R_eq']:.4f} omega_e={om['omega_e_cm1']:.0f} "
                 f"({om['omega_e_err_pct']:+.1f}%)" if om else ""))

    print("\nconvergence hygiene:")
    for r in rows:
        print(f"  max_n={r['max_n']} {r['tag']:>10s}: all Davidson converged="
              f"{r['all_conv']}, worst |r| = {r['max_res']:.2e}")

    print("\nenergies at R_true (E_coupled convention -- carries the R-independent "
          "core double-count):")
    for r in rows:
        i = r['Rs'].index(R_TRUE) if R_TRUE in r['Rs'] else None
        if i is not None:
            print(f"  max_n={r['max_n']} {r['tag']:>10s}: {r['Es'][i]:+.12f}")

    # ---------------- verdict arithmetic -------------------------------
    print("\n" + "=" * 104)
    print("VERDICT ARITHMETIC")
    print("=" * 104)
    for tag in ('faithful', 'corrected'):
        seq = [r for r in rows if r['tag'] == tag]
        if len(seq) < 3:
            continue
        print(f"\n[{tag}]")
        c = [r['f3']['curv_over_ktrue_at_Rtrue'] for r in seq]
        w = [r['f3']['own_min'].get('omega_e_err_pct') for r in seq]
        k = [r['f3']['own_min'].get('curv_over_ktrue') for r in seq]
        n = [r['max_n'] for r in seq]
        print(f"  curv(R_true)/k_true : " + "  ".join(f"n{a}={b:.3f}x" for a, b in zip(n, c)))
        print(f"  curv(own min)/k_true: " + "  ".join(f"n{a}={b:.3f}x" for a, b in zip(n, k)))
        print(f"  omega_e error       : " + "  ".join(f"n{a}={b:+.1f}%" for a, b in zip(n, w)))
        if len(c) >= 3:
            d23 = c[1] - c[0]; d34 = c[2] - c[1]
            print(f"  step 2->3 = {d23:+.4f}x ; step 3->4 = {d34:+.4f}x  "
                  f"(healing would need a LARGE negative step toward 1.000x)")
            print(f"  omega_e steps: 2->3 {w[1]-w[0]:+.1f} pp ; 3->4 {w[2]-w[1]:+.1f} pp "
                  f"(healing target: -45 pp)")
            frac = abs(k[2] - 1.0) / abs(k[0] - 1.0) if k[0] is not None else None
            if frac is not None:
                print(f"  residual stiffness excess at n_max=4 relative to n_max=2: "
                      f"{frac*100:.1f}% (100% = fully frozen, 0% = fully healed)")
