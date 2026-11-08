"""UNREGISTERED follow-up to WA1/WA2 (does not change any score): finer grid (level 2) for higher l, and a fit with a
1/l correction, log|F_l| = a - p log l + l log r + c/l. Also the local ratio r_l = (F_{l+1}/F_l) (l/(l+1))^{-p}."""
import sys
import numpy as np
sys.path.insert(0, 'debug')
from wa_channels import F_coeffs

for (a, b) in [(1.809, 2.5), (1.809, 1.0), (1.809, 1.809)]:
    F1 = F_coeffs(a, b, 1.0, 1.0, 60, level=1)
    F2 = F_coeffs(a, b, 1.0, 1.0, 60, level=2)
    floor = np.abs(F2 - F1)
    blk = []
    for l in range(15, 61):
        if abs(F2[l]) >= 100 * floor[l]: blk.append(l)
        else: break
    ls = np.array(blk, dtype=float)
    y = np.log(np.abs(F2[blk]))
    A = np.vstack([np.ones_like(ls), -np.log(ls), ls, 1 / ls]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    lo, hi = min(a, b), max(a, b)
    print(f"D=({a},{b}) level2 usable l {blk[0]}..{blk[-1]}: with 1/l term r={np.exp(c[2]):.5f} (D</D> {lo/hi:.5f},"
          f" rel {np.exp(c[2])/(lo/hi)-1:+.2e}) p={c[1]:.3f}")
    # tail-only fit (last 12 usable), no 1/l
    tl = ls[-12:]; ty = y[-12:]
    A2 = np.vstack([np.ones_like(tl), -np.log(tl), tl]).T
    c2, *_ = np.linalg.lstsq(A2, ty, rcond=None)
    print(f"    tail fit l {int(tl[0])}..{int(tl[-1])}: r={np.exp(c2[2]):.5f} (rel {np.exp(c2[2])/(lo/hi)-1:+.2e}) p={c2[1]:.3f}")
