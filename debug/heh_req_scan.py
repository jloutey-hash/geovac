"""HeH+ R_eq PoC: does the variational-prolate (H2-style) recipe give the correct
bond length for a 2-electron HETERONUCLEAR system, i.e. NO balanced-style drift?

If yes, it is direct evidence that the prolate route fixes geometry for exactly
the class of frozen-core LiH (2 valence e, heteronuclear) -- before we build the
Li core-shield. Reference: HeH+ X^1Sigma+ R_e = 1.4632 bohr (0.7743 Angstrom),
E_tot(R_e) ~ -2.97869 Ha.

Reuses today's engine (debug/prolate_r12_mpf.assemble_hetero + solve_canonical_mpf,
r12 p={0,1}, full parity). We scan R, optionally optimize the basis exponent alpha
at each R (the variational "cloud tightening"), fit E_tot(R), and read off R_eq.

Run from root:  python debug/heh_req_scan.py [jmax] [lmax] [l_neumann] [dps]
"""
import os
import sys
import time
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac import prolate_recondition as pr        # noqa: E402
import prolate_r12_mpf as m                          # noqa: E402
from heh_probe import build_basis_full, ZA, ZB       # noqa: E402
from r12ci_first_energy import solve_canonical       # noqa: E402  (float64 canonical)

R_REF = 1.4632          # HeH+ equilibrium bond length (bohr) -- the target
E_REF = -2.97869        # HeH+ total energy near R_e (Ha)

JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 1
L_NEUMANN = int(sys.argv[3]) if len(sys.argv) > 3 else 16
DPS = int(sys.argv[4]) if len(sys.argv) > 4 else 30

ALPHAS = [1.3, 1.6, 1.9]   # optimum is low/shallow (verified: wider grid picks ~1.4-1.9)
# 5th arg: comma-separated R grid override (for a focused convergence check)
if len(sys.argv) > 5:
    R_GRID = [float(x) for x in sys.argv[5].split(',')]
else:
    R_GRID = [1.10, 1.25, 1.4632, 1.65, 1.85, 2.10]


def E_tot(R, alpha):
    # float64 canonical solve == mpf solve here (verified: -2.967037 both at (2,2));
    # assembly stays mpf-accurate, only the eigensolve is float64 (the mpf eig was
    # the n^3 bottleneck).
    basis = build_basis_full(JMAX, LMAX, alpha, p_set=(0, 1))
    S, H = m.assemble_hetero(basis, R, alpha, ZA, ZB,
                             l_neumann=L_NEUMANN, dps=DPS)   # float64 out
    e = solve_canonical(S, H)[0]
    return e + ZA * ZB / R


def main():
    t0 = time.time()
    print(f"HeH+ R_eq scan  (j,l)=({JMAX},{LMAX})  l_neumann={L_NEUMANN} dps={DPS}  "
          f"ref R_e={R_REF} bohr", flush=True)
    Emin, best_a = [], []
    for R in R_GRID:
        vals = []
        for a in ALPHAS:
            e = E_tot(R, a)
            vals.append(e)
        i = int(np.argmin(vals))
        Emin.append(vals[i]); best_a.append(ALPHAS[i])
        print(f"  R={R:.4f}  E_min={vals[i]:.6f}  (a*={ALPHAS[i]})  "
              f"all={[round(v,5) for v in vals]}  [{time.time()-t0:.0f}s]", flush=True)

    R = np.array(R_GRID); E = np.array(Emin)
    p = np.poly1d(np.polyfit(R - R_REF, E, min(4, len(R) - 1)))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_REF for r in roots if ddp(r) > 0 and R.min() < r + R_REF < R.max()]
    Req = float(min(cand, key=lambda rr: p(rr - R_REF))) if cand else float('nan')
    err = (Req - R_REF) / R_REF * 100 if np.isfinite(Req) else float('nan')

    print("\n" + "=" * 66)
    print(f"R_eq(computed) = {Req:.4f} bohr    ref = {R_REF} bohr    "
          f"drift = {err:+.1f}%")
    print(f"E_min(grid)    = {E.min():.6f} Ha   ref = {E_REF} Ha")
    print("=" * 66)
    print("VERDICT:", "NO DRIFT -- prolate fixes geometry" if abs(err) < 5
          else ("MILD DRIFT" if abs(err) < 12 else "DRIFTS"))

    out = {'R_grid': R_GRID, 'E_min': [float(x) for x in Emin], 'best_alpha': best_a,
           'R_eq': Req, 'R_eq_err_pct': err, 'R_ref': R_REF, 'E_ref': E_REF,
           'basis': [JMAX, LMAX], 'l_neumann': L_NEUMANN, 'dps': DPS}
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/heh_req_scan.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"[saved] debug/data/heh_req_scan.json  ({time.time()-t0:.0f}s)")


if __name__ == '__main__':
    main()
